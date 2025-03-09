use std::collections::HashMap;

use lib_genome_kit::amino::AminoAcid;
use lib_genome_kit::blosum::Blosum;

/// Eay optimization opportunities here, we use a recursive take/skip approach,
/// but a bottom-up DP table could speedup the entire search drastically...
///
/// - https://www.cs.otago.ac.nz/cosc348/alignments/Lecture05_GlobalAlignment.pdf
/// - https://www.ncbi.nlm.nih.gov/nuccore/NC_050012.1?report=fasta
fn align_recurse(
  seq_a: &[AminoAcid],
  seq_b: &[AminoAcid],
  pos @ (i, j): (usize, usize),
  gap_penalty: i32,
  table: &impl Blosum,
  cache: &mut HashMap<(usize, usize), i32>,
) -> i32 {
  // take/skip rationale:
  // at any given moment, we want to handle 3 cases;
  // - we want to take both bases
  // - we want to skip one of seq_a
  // - we want to skip one of seq_b

  // base cases;
  // when we exceeded the length, we pad rest of remaining sequence with the gap
  if i >= seq_a.len() {
    return (seq_b.len().saturating_sub(j) as i32) * gap_penalty;
  }
  if j >= seq_b.len() {
    return (seq_a.len().saturating_sub(i) as i32) * gap_penalty;
  }

  // fetch from cache
  if let Some(score) = cache.get(&pos).copied() {
    return score;
  }

  // now we take 3 paths, and choose the best...
  let score = table.score(seq_a[i], seq_b[j]);
  let both = align_recurse(seq_a, seq_b, (i + 1, j + 1), gap_penalty, table, cache) + score;
  let skip_b = align_recurse(seq_a, seq_b, (i + 1, j), gap_penalty, table, cache) + gap_penalty;
  let skip_a = align_recurse(seq_a, seq_b, (i, j + 1), gap_penalty, table, cache) + gap_penalty;
  let best_so_far = both.max(skip_b).max(skip_a);
  // insert our result so far into the cache, we will probably be here again...
  cache.insert(pos, best_so_far);
  best_so_far
}

/// We want to be able to run the needleman-wunsch algorithm on any of the tables,
/// and blosum tables are all kinda the same, so we just default implement it...
pub trait Needleman: Blosum {
  fn needleman_wunsch(seq_a: &[AminoAcid], seq_b: &[AminoAcid]) -> i32
  where
    Self: Default,
  {
    let mut cache = HashMap::new();
    align_recurse(seq_a, seq_b, (0, 0), -5, &Self::default(), &mut cache)
  }
}

// blanket impl that actually gives tables access to the needleman algorithm
impl<B> Needleman for B where B: Blosum + Default {}

#[cfg(test)]
mod tests {
  use super::*;

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
