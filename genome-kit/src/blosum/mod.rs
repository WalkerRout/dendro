use crate::amino::AminoAcid;

pub mod blosum62;
pub use blosum62::*;

pub mod blosum45;
pub use blosum45::*;

pub(crate) const ROWS: usize = 24;
pub(crate) const COLS: usize = 24;

/// A trait to abstract over kinds of BLOSUM tables...
///
/// - https://ftp.ncbi.nlm.nih.gov/blast/matrices/
pub trait Blosum {
  /// Returns the substitution score for two amino acids
  fn score(&self, a: AminoAcid, b: AminoAcid) -> i32;
}

// all tables are pretty much the same in terms of their size and positions,
// so we provide a generic way to access a table index given some acid...
#[inline]
pub(crate) const fn acid_to_index(acid: AminoAcid) -> usize {
  match acid {
    AminoAcid::Alanine => 0,
    AminoAcid::Arginine => 1,
    AminoAcid::Asparagine => 2,
    AminoAcid::AsparticAcid => 3,
    AminoAcid::Cysteine => 4,
    AminoAcid::Glutamine => 5,
    AminoAcid::GlutamicAcid => 6,
    AminoAcid::Glycine => 7,
    AminoAcid::Histidine => 8,
    AminoAcid::Isoleucine => 9,
    AminoAcid::Leucine => 10,
    AminoAcid::Lysine => 11,
    AminoAcid::Methionine => 12,
    AminoAcid::Phenylalanine => 13,
    AminoAcid::Proline => 14,
    AminoAcid::Serine => 15,
    AminoAcid::Threonine => 16,
    AminoAcid::Tryptophan => 17,
    AminoAcid::Tyrosine => 18,
    AminoAcid::Valine => 19,
    AminoAcid::Asx => 20,
    AminoAcid::Glx => 21,
    AminoAcid::Unknown => 22,
    AminoAcid::Stop => 23,
  }
}

#[inline]
pub(crate) const fn score_for(matrix: &[i32; ROWS * COLS], a: AminoAcid, b: AminoAcid) -> i32 {
  let (i, j) = (acid_to_index(a), acid_to_index(b));
  matrix[i * COLS + j]
}
