use lib_genome_kit::amino::AminoAcid;
use lib_genome_kit::blosum::Blosum;

/// We invert the original recursive problem, starting with the base cases (gaps)
/// in the lowest levels (first row/col) and work our way up from smaller solutions
/// to a bigger ones
///
/// - https://www.cs.otago.ac.nz/cosc348/alignments/Lecture05_GlobalAlignment.pdf
/// - https://www.ncbi.nlm.nih.gov/nuccore/NC_050012.1?report=fasta
fn align_dp_table(
  seq_a: &[AminoAcid],
  seq_b: &[AminoAcid],
  gap_penalty: i32,
  table: &impl Blosum,
) -> i32 {
  let m = seq_a.len();
  let n = seq_b.len();

  // m+1 rows and n+1 cols to represent 0 elements in either sequence...
  let mut dp = vec![0; (m + 1) * (n + 1)];
  // lil helper to compute index in 1d dp vector, take care to use real width
  let idx = |i: usize, j: usize| i * (n + 1) + j;

  // aligning an empty sequence with a not empty sequence incurs a cumulative
  // gap penalty, this is our lowest level (base case)
  for i in 0..=m {
    // rows corresponds to seq_a
    dp[idx(i, 0)] = (i as i32) * gap_penalty;
  }
  for j in 0..=n {
    // cols correspond to seq_b
    dp[idx(0, j)] = (j as i32) * gap_penalty;
  }

  // start above gaps (1..), each cell dp[i][j] represents best alignment score for
  // prefixes seq_a[0..i] and seq_b[0..j]
  for i in 1..=m {
    for j in 1..=n {
      // identical to recursive solution, just upside down
      let score = table.score(seq_a[i - 1], seq_b[j - 1]);
      let diag = dp[idx(i - 1, j - 1)] + score;
      let up = dp[idx(i - 1, j)] + gap_penalty;
      let left = dp[idx(i, j - 1)] + gap_penalty;
      // we depend on diagonal, cell above, cell to left
      dp[idx(i, j)] = diag.max(up).max(left);
    }
  }

  // bottom right holds best alignment score...
  dp[idx(m, n)]
}

/// We want to be able to run the needleman-wunsch algorithm on any of the tables,
/// and blosum tables are all kinda the same, so we just default implement it...
pub trait Needleman {
  fn needleman_wunsch(seq_a: &[AminoAcid], seq_b: &[AminoAcid]) -> i32;
}

// blanket impl that actually gives tables access to the needleman algorithm
impl<B> Needleman for B 
where 
  B: Blosum + Default
{
  #[inline]
  fn needleman_wunsch(seq_a: &[AminoAcid], seq_b: &[AminoAcid]) -> i32 {
    align_dp_table(seq_a, seq_b, -5, &Self::default())
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  use lib_genome_kit::genome::Genome;
  // we only test with `Blosum62` as of now...
  use lib_genome_kit::blosum::Blosum62;

  mod needleman {
    use super::*;

    #[test]
    fn needleman_wunsch_empty_sequences() {
      let seq1 = Genome::from("".chars());
      let seq2 = Genome::from("".chars());
      // both empty, 0 * gap = 0
      let score = Blosum62::needleman_wunsch(&seq1, &seq2);
      assert_eq!(score, 0);
    }

    #[test]
    fn needleman_wunsch_one_empty_sequence() {
      let seq1 = Genome::from("".chars());
      let seq2 = Genome::from("ARN".chars());
      // if a sequence is empty, we run the gap for the rest of the opposing
      // sequence (3 chars * -5 gap)
      let score = Blosum62::needleman_wunsch(&seq1, &seq2);
      assert_eq!(score, -15);
    }

    #[test]
    fn needleman_wunsch_single_char_match() {
      let seq1 = Genome::from("A".chars());
      let seq2 = Genome::from("A".chars());
      // score of 4 for A-A match
      let score = Blosum62::needleman_wunsch(&seq1, &seq2);
      assert_eq!(score, 4);
    }

    #[test]
    fn needleman_wunsch_single_char_mismatch() {
      let seq1 = Genome::from("A".chars());
      let seq2 = Genome::from("R".chars());
      // A-R alignment score (-1) better than gap (-5)
      let score = Blosum62::needleman_wunsch(&seq1, &seq2);
      assert_eq!(score, -1);
    }

    #[test]
    fn needleman_wunsch_longer_sequence() {
      let seq1 = Genome::from("PLEASANTLY".chars());
      let seq2 = Genome::from("MEANLY".chars());
      let score = Blosum62::needleman_wunsch(&seq1, &seq2);
      assert_eq!(score, 8);
    }
  }
}
