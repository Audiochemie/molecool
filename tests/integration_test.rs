mod file_setup;
#[cfg(test)]
mod tests {
    use std::iter::zip;

    use crate::file_setup;
    use approx::assert_relative_eq;
    use coordinate_systems::cartesian::Cartesian;
    use molecool::pse_data::PSE_MASSES;
    use molecool::{
        atom::Atom,
        molekel::{get_centre_of_mass, get_distances, get_oop, Molecule},
    };
    use nalgebra::Point3;
    use qc_file_parsers::xyz::Xyz;

    #[test]
    fn test_atom_from_symbolic() {
        let mut test_file = file_setup::setup_allene_symbolic().unwrap();
        let test_parsed = Xyz::new(&mut test_file, "Ang").unwrap();
        test_parsed.lines.iter().for_each(|l| println!("{:?}", l));
        let test_atom = Atom::from(test_parsed.lines[0].clone());
        assert_eq!(test_atom.z_value, 6);
        assert_eq!(
            test_atom.coordinates,
            Cartesian(Point3::new(0.0_f32, 0.0_f32, 1.889_726_f32))
        );
        assert_eq!(test_atom.atomic_mass, PSE_MASSES[6]);
        let test_atom = Atom::from(test_parsed.lines[6].clone());
        assert_eq!(test_atom.z_value, 1);
        assert_eq!(
            test_atom.coordinates,
            Cartesian(Point3::<f32>::new(-3.495_993_1, 1.157_216_1, 3.046_942))
        );
        assert_eq!(test_atom.atomic_mass, PSE_MASSES[1]);
    }

    #[test]
    fn test_molecule_from_numeric() {
        let mut test_file = file_setup::setup_acetaldehyde_numeric().unwrap();
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let expected_mol = (PSE_MASSES[8] + 2.0 * PSE_MASSES[6] + 4.0 * PSE_MASSES[1]) as f32;
        assert_relative_eq!(mol.molar_mass, expected_mol)
    }

    #[test]
    fn test_molecule_from_symbolic() {
        // Test acetaldehyde
        let mut test_file = file_setup::setup_acetaldehyde_symbolic().unwrap();
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let expected_mol = (PSE_MASSES[8] + 2.0 * PSE_MASSES[6] + 4.0 * PSE_MASSES[1]) as f32;
        assert_relative_eq!(mol.molar_mass, expected_mol);
        // Test allene
        let mut test_file = file_setup::setup_allene_symbolic().unwrap();
        let test_parsed = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let test_atom = &mol.building_atoms[0];
        assert_eq!(test_atom.z_value, 6);
        assert_eq!(
            test_atom.coordinates,
            Cartesian(Point3::new(0.0_f32, 0.0_f32, 1.889_726_f32))
        );
        assert_eq!(test_atom.atomic_mass, PSE_MASSES[6]);
    }

    #[test]
    fn test_get_distances() {
        let mut test_file = file_setup::setup_acetaldehyde_numeric().unwrap();
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "Ang").unwrap();
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
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let quadruples: [(usize, usize, usize, usize); 4] =
            [(0, 5, 4, 6), (0, 4, 5, 6), (6, 0, 5, 4), (4, 1, 0, 5)];
        let expected: [f32; 4] = [19.939_726, -19.850_523, -31.064_344, -53.651_534];
        for (angle, quadruple) in zip(expected, quadruples) {
            let h = &mol.building_atoms[quadruple.0];
            let l = &mol.building_atoms[quadruple.1];
            let c = &mol.building_atoms[quadruple.2];
            let r = &mol.building_atoms[quadruple.3];
            assert_relative_eq!(angle, get_oop(l, c, r, h).to_degrees(),);
        }
    }

    #[test]
    fn test_get_com() {
        let mut test_file = file_setup::setup_allene_symbolic().unwrap();
        let test_parsed: Xyz<f64> = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let com: Point3<f64> = get_centre_of_mass(&mol);
        let expected: Point3<f64> = Point3::new(0.0, 0.0, 1.889_725_99);
        assert_relative_eq!(com, expected, epsilon = 1e-8);

        let mut test_file = file_setup::setup_acetaldehyde_numeric().unwrap();
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "Ang").unwrap();
        let mol = Molecule::from(test_parsed);
        let com: Point3<f32> = get_centre_of_mass(&mol);
        let expected: Point3<f32> = Point3::new(0.644_949_26, 0.0, 2.316_637_92);
        assert_relative_eq!(com, expected, epsilon = 1e-3);
    }
}
