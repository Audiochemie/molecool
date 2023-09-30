mod file_setup;
#[cfg(test)]
mod tests {
    use std::f32::consts::FRAC_1_PI;
    use std::iter::zip;

    use crate::file_setup;
    use approx::assert_relative_eq;
    use coordinate_systems::cartesian::Cartesian;
    use molecool::pse_data::PSE_MASSES;
    use molecool::{
        atom::Atom,
        molekel::{get_distances, get_oop, Molecule},
    };
    use nalgebra::Point3;
    use qc_file_parsers::xyz::Xyz;
    #[test]
    fn test_atom_from_symbolic() {
        let mut test_file = file_setup::setup_allene_symbolic().unwrap();
        let test_parsed = Xyz::new(&mut test_file, "Ang").unwrap();
        let test_atom = Atom::from(test_parsed.lines[0].clone());
        assert_eq!(test_atom.z_value, 6);
        assert_eq!(
            test_atom.coordinates,
            Cartesian(Point3::new(0.0_f32, 0.0_f32, 1.889_726_f32))
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
    #[test]
    fn test_get_distances() {
        let mut test_file = file_setup::setup_acetaldehyde_numeric().unwrap();
        let test_parsed = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let distance_vec = get_distances(&mol);
        let expected = Vec::from([
            2.84511_f32,
            4.55395_f32,
            4.19912_f32,
            2.06517_f32,
            2.07407_f32,
            2.07407_f32,
            2.29803_f32,
            2.09811_f32,
            4.04342_f32,
            4.05133_f32,
            4.05133_f32,
            3.81330_f32,
            4.84040_f32,
            5.89151_f32,
            5.89151_f32,
            5.87463_f32,
            4.83836_f32,
            4.83836_f32,
            3.38971_f32,
            3.38971_f32,
            3.33994_f32,
        ]);
        let natms = mol.building_atoms.len();
        let iter_space = (natms * natms - natms) / 2;
        assert_eq!(iter_space, expected.len());
        for i in 0..iter_space {
            assert_relative_eq!(distance_vec[i].norm(), expected[i], epsilon = 0.00001);
        }
    }
    #[test]
    fn test_get_oop() {
        let mut test_file = file_setup::setup_acetaldehyde_numeric().unwrap();
        let test_parsed = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let quadruples: [(usize, usize, usize, usize); 4] =
            [(0, 5, 4, 6), (0, 4, 5, 6), (6, 0, 5, 4), (4, 1, 0, 5)];
        let expected: [f32; 4] = [19.939_726_f32, -19.850_523_f32, -31.064_344_f32, -53.651_534_f32];
        for (angle, quadruple) in zip(expected, quadruples) {
            let h = &mol.building_atoms[quadruple.0];
            let l = &mol.building_atoms[quadruple.1];
            let c = &mol.building_atoms[quadruple.2];
            let r = &mol.building_atoms[quadruple.3];
            assert_relative_eq!(
                angle,
                get_oop(l, c, r, h) * 180.0_f32 * FRAC_1_PI,
                epsilon = 0.00001
            );
        }
    }
}
