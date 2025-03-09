use std::ops::{Deref, DerefMut};

use crate::amino::AminoAcid;

/// A genome is a list of amino acids for something like COX3
#[derive(Debug, Clone, PartialEq)]
pub struct Genome(Vec<AminoAcid>);

impl<T, I> From<I> for Genome
where
  T: Into<AminoAcid>,
  I: IntoIterator<Item = T>,
{
  /// We should be able to get a genome from any list of things treatable as
  /// amino acids...
  fn from(iter: I) -> Self {
    Self(iter.into_iter().map(Into::into).collect())
  }
}

// our `Genome` struct is a sort of smart-pointer wrapper for the actual amino
// acids it contains...

impl Deref for Genome {
  type Target = Vec<AminoAcid>;

  fn deref(&self) -> &Self::Target {
    &self.0
  }
}

impl DerefMut for Genome {
  fn deref_mut(&mut self) -> &mut Self::Target {
    &mut self.0
  }
}
