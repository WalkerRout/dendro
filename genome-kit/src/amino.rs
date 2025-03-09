/// All relevant amino acids...
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AminoAcid {
  Alanine,       // A
  Arginine,      // R
  Asparagine,    // N
  AsparticAcid,  // D
  Cysteine,      // C
  Glutamine,     // Q
  GlutamicAcid,  // E
  Glycine,       // G
  Histidine,     // H
  Isoleucine,    // I
  Leucine,       // L
  Lysine,        // K
  Methionine,    // M
  Phenylalanine, // F
  Proline,       // P
  Serine,        // S
  Threonine,     // T
  Tryptophan,    // W
  Tyrosine,      // Y
  Valine,        // V
  Asx,           // B (Aspartic Acid or Asparagine)
  Glx,           // Z (Glutamic Acid or Glutamine)
  Unknown,       // X (Unknown or unimportant)
  Stop,          // * (Stop codon)
}

impl From<char> for AminoAcid {
  fn from(c: char) -> Self {
    match c {
      'A' => AminoAcid::Alanine,
      'R' => AminoAcid::Arginine,
      'N' => AminoAcid::Asparagine,
      'D' => AminoAcid::AsparticAcid,
      'C' => AminoAcid::Cysteine,
      'Q' => AminoAcid::Glutamine,
      'E' => AminoAcid::GlutamicAcid,
      'G' => AminoAcid::Glycine,
      'H' => AminoAcid::Histidine,
      'I' => AminoAcid::Isoleucine,
      'L' => AminoAcid::Leucine,
      'K' => AminoAcid::Lysine,
      'M' => AminoAcid::Methionine,
      'F' => AminoAcid::Phenylalanine,
      'P' => AminoAcid::Proline,
      'S' => AminoAcid::Serine,
      'T' => AminoAcid::Threonine,
      'W' => AminoAcid::Tryptophan,
      'Y' => AminoAcid::Tyrosine,
      'V' => AminoAcid::Valine,
      'B' => AminoAcid::Asx,
      'Z' => AminoAcid::Glx,
      'X' => AminoAcid::Unknown,
      '*' => AminoAcid::Stop,
      _ => AminoAcid::Unknown, // Default case for invalid characters
    }
  }
}
