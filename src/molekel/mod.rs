use qc_file_parsers::xyz::Xyz;
use super::atom::Atom;
/// Represents a molecule in the _atoms in molecules_ sense.
pub struct Molecule {
    /// Molecules have a total molar mass.
    pub molar_mass: f32,
    /// Molecules are commonly interpreted as a collection of building atoms.
    pub building_atoms: Vec<Atom>,
}

impl From<Xyz> for Molecule {
    /// Cast Xyz struct into a Molecule.
    fn from(value: Xyz) -> Self {
        let mut atoms = Vec::with_capacity(value.number_of_atoms);
        value.lines.iter().for_each(|xyz_line| atoms.push(Atom::from(xyz_line.clone())));
        let molar_mass = atoms.iter().fold(0.0_f32, |at, at2| at + at2.atomic_mass);
        Self { molar_mass, building_atoms: atoms }
    }
}


