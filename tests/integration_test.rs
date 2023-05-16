mod file_setup;
#[cfg(test)]
mod tests {
    use crate::file_setup;
    use approx::assert_relative_eq;
    use molecool::{Atom, Coordinates, Molecule};
    use molecool::pse_data::PSE_MASSES;
    use nalgebra::Point3;
    use qc_file_parsers::Xyz;
    #[test]
    fn test_atom_from_symbolic() {
        let mut test_file = file_setup::setup_allene_symbolic().unwrap();
        let test_parsed = Xyz::new(&mut test_file, "Ang").unwrap();
        let test_atom = Atom::from(test_parsed.lines[0].clone());
        assert_eq!(test_atom.z_value, 6);
        assert_eq!(
            test_atom.coordinates,
            Coordinates::Cartesian(Point3::new(0.0_f32, 0.0_f32, 1.889_726_f32))
        );
        assert_eq!(test_atom.atomic_mass, PSE_MASSES[6]);
    }
    #[test]
    fn test_molecule_from_numeric() {
        let mut test_file = file_setup::setup_acetaldehyde_numeric().unwrap();
        let test_parsed = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let expected_mol = PSE_MASSES[8] + 2.0_f32 * PSE_MASSES[6] + 4.0_f32 * PSE_MASSES[1];
        assert_relative_eq!(mol.molar_mass, expected_mol)
    }
}
