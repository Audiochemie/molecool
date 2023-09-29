use super::atom::Atom;
use coordinate_systems::DistanceTo;
use qc_file_parsers::xyz::Xyz;
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
        value
            .lines
            .iter()
            .for_each(|xyz_line| atoms.push(Atom::from(xyz_line.clone())));
        let molar_mass = atoms.iter().fold(0.0_f32, |at, at2| at + at2.atomic_mass);
        Self {
            molar_mass,
            building_atoms: atoms,
        }
    }
}

pub fn get_distances(mol: &Molecule) -> Vec<f32> {
    let num_atms = mol.building_atoms.len();
    let mut distance_vector: Vec<f32> = Vec::with_capacity((num_atms * num_atms - num_atms) / 2);
    for i_atm in 0..mol.building_atoms.len() - 1 {
        let ref_atm_coor = &mol.building_atoms[i_atm].coordinates;
        for i_nxt_atm in i_atm+1..mol.building_atoms.len() {
            let nxt_atm_coor = &mol.building_atoms[i_nxt_atm].coordinates;
            distance_vector.push(ref_atm_coor.distance_to(nxt_atm_coor));
        }
    }
    distance_vector
}
