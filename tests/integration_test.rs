mod file_setup;
#[cfg(test)]
mod water_integrals;
#[cfg(test)]
mod tests {
    use std::iter::zip;

    use crate::file_setup;
    use approx::assert_relative_eq;
    use coordinate_systems::cartesian::Cartesian;

    use molecool::molekel::principle_moments_of_inertia;
    use molecool::molekel::properties::{mass_weight_hessian, Hessian};
    use molecool::pse_data::PSE_MASSES;
    use molecool::{
        atom::Atom,
        molekel::{get_centre_of_mass, get_distances, get_i_tensor, get_oop, Molecule},
    };
    use nalgebra::{DMatrix, DVector, Matrix3, Point3, Vector3};
    use qc_file_parsers::array_text::parse_text_into_matrix;
    use qc_file_parsers::xyz::Xyz;

    #[test]
    fn test_atom_from_symbolic() {
        let mut test_file = file_setup::setup_allene_symbolic().unwrap();
        let test_parsed = Xyz::new(&mut test_file, "bohr").unwrap();
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
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "bohr").unwrap();
        let mol = Molecule::from(test_parsed);
        let expected_mol = (PSE_MASSES[8] + 2.0 * PSE_MASSES[6] + 4.0 * PSE_MASSES[1]) as f32;
        assert_relative_eq!(mol.molar_mass, expected_mol);
        let mut test_file = file_setup::setup_chlorbenzene_geo().unwrap();
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "bohr").unwrap();
        let mol = Molecule::from(test_parsed);
        let expected_mol = (PSE_MASSES[6] * 4.0 + 7.0 * PSE_MASSES[1] + PSE_MASSES[17]) as f32;
        assert_relative_eq!(mol.molar_mass, expected_mol);
        assert_eq!(
            mol.building_atoms[11].coordinates,
            Cartesian(Point3::<f32>::new(
                -3.121_360_3000,
                -0.895_719_0610,
                -0.095_875_3352
            ))
        );
    }

    #[test]
    fn test_molecule_from_symbolic() {
        // Test acetaldehyde
        let mut test_file = file_setup::setup_acetaldehyde_symbolic().unwrap();
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "bohr").unwrap();
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
            assert_relative_eq!(distance_vec[i].norm(), expected[i], epsilon = 1e-5);
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
        let test_parsed: Xyz<f64> = Xyz::new(&mut test_file, "bohr").unwrap();
        let mol = Molecule::from(test_parsed);
        let com: Point3<f64> = get_centre_of_mass(&mol);
        let expected: Point3<f64> = Point3::new(0.0, 0.0, 1.889_725_99);
        assert_relative_eq!(com, expected, epsilon = 1e-8);

        let mut test_file = file_setup::setup_acetaldehyde_numeric().unwrap();
        let test_parsed: Xyz<f64> = Xyz::new(&mut test_file, "bohr").unwrap();
        let mol = Molecule::from(test_parsed);
        let com: Point3<f64> = get_centre_of_mass(&mol);
        let expected: Point3<f64> = Point3::new(0.644_949_26, 0.0, 2.316_638);
        assert_relative_eq!(com, expected, epsilon = 1e-3);
    }

    #[test]
    fn test_i_tensor() {
        let mut test_file = file_setup::setup_allene_symbolic().unwrap();
        let test_parsed: Xyz<f32> = Xyz::new(&mut test_file, "bohr").unwrap();
        let mol = Molecule::from(test_parsed);
        let i_tensor = get_i_tensor(&mol);
        let expected = Matrix3::<f32>::new(
            10.797_024, 0.0000000, 0.0000000, 0.0000000, 210.867_28, 0.0000000, 0.0000000,
            0.0000000, 210.867_28,
        );
        assert_relative_eq!(i_tensor, expected, epsilon = 1e0);
        let comput = principle_moments_of_inertia(i_tensor);
        assert_relative_eq!(
            comput,
            Vector3::<f32>::new(10.797_024, 210.867_28, 210.867_28),
            epsilon = 1e0
        );
    }

    #[test]
    fn test_mass_weight_hessian() {
        let mut test_file = file_setup::setup_water_geo().unwrap();
        let test_parsed: Xyz<f64> = Xyz::new(&mut test_file, "bohr").unwrap();
        let test_mol = Molecule::from(test_parsed);
        let mut hess_file = file_setup::setup_water_hessian().unwrap();
        let parse_matrix: DMatrix<f64> = parse_text_into_matrix(&mut hess_file, " ", 9, 9);
        let mut hessian = Hessian(parse_matrix);
        mass_weight_hessian(&test_mol, &mut hessian);
        let expected: DMatrix<f64> = DMatrix::from_iterator(
            9,
            9,
            [
                0.0057996, 0.0000000, 0.0000000, -0.0115523, 0.0000000, 0.0000000, -0.0115523,
                0.0000000, 0.0000000, 0.0000000, 0.0198271, 0.0000000, 0.0000000, -0.0394937,
                0.0199304, 0.0000000, -0.0394937, -0.0199304, 0.0000000, 0.0000000, 0.0175112,
                0.0000000, 0.0086617, -0.0348807, 0.0000000, -0.0086617, -0.0348807, -0.0115523,
                0.0000000, 0.0000000, 0.0510672, 0.0000000, 0.0000000, -0.0050452, 0.0000000,
                0.0000000, 0.0000000, -0.0394937, 0.0086617, 0.0000000, 0.1716643, -0.0569527,
                0.0000000, -0.0143291, 0.0224462, 0.0000000, 0.0199304, -0.0348807, 0.0000000,
                -0.0569527, 0.1258526, 0.0000000, -0.0224462, 0.0131055, -0.0115523, 0.0000000,
                0.0000000, -0.0050452, 0.0000000, 0.0000000, 0.0510672, 0.0000000, 0.0000000,
                0.0000000, -0.0394937, -0.0086617, 0.0000000, -0.0143291, -0.0224462, 0.0000000,
                0.1716643, 0.0569527, 0.0000000, -0.0199304, -0.0348807, 0.0000000, 0.0224462,
                0.0131055, 0.0000000, 0.0569527, 0.1258526,
            ],
        );
        let mut computed_eigvals = hessian.0.symmetric_eigenvalues();
        computed_eigvals
            .as_mut_slice()
            .sort_by(|a, b| b.partial_cmp(a).unwrap());
        for (h, e) in hessian.0.cast::<f64>().iter().zip(expected.iter()) {
            assert_relative_eq!(h, e, epsilon = 1e-5);
        }
        let expected_eigvals: DVector<f64> = DVector::from_iterator(
            9,
            [
                0.2351542439,
                0.2107113210,
                0.1317512832,
                0.0561123974,
                0.0547551476,
                0.0518216614,
                0.0000000000,
                0.0000000000,
                0.0000000000,
            ],
        );
        assert_relative_eq!(expected_eigvals, computed_eigvals, epsilon = 1e-5);

        let mut test_file = file_setup::setup_chlorbenzene_geo().unwrap();
        let test_parsed: Xyz<f64> = Xyz::new(&mut test_file, "bohr").unwrap();
        let test_mol = Molecule::from(test_parsed);
        let mut hess_file = file_setup::setup_chlorbenzene_hessian().unwrap();
        let parse_matrix: DMatrix<f64> = parse_text_into_matrix(&mut hess_file, " ", 36, 36);
        let mut hessian = Hessian(parse_matrix);
        mass_weight_hessian(&test_mol, &mut hessian);
        let expected_eigvals: DVector<f64> = DVector::from_iterator(
            36,
            [
                0.5504194388,
                0.5372530704,
                0.5351184811,
                0.5255730312,
                0.5062966607,
                0.4992575932,
                0.4845191221,
                0.1594678267,
                0.1246207278,
                0.1240628138,
                0.1124389759,
                0.1102183592,
                0.0940432225,
                0.0897847299,
                0.0836000497,
                0.0690595364,
                0.0640720123,
                0.0595822760,
                0.0569922250,
                0.0512084824,
                0.0498908695,
                0.0410337099,
                0.0283830220,
                0.0213722484,
                0.0088837125,
                0.0065288923,
                0.0046443780,
                0.0028070052,
                0.0022837821,
                0.0000001155,
                0.0000000871,
                0.0000000051,
                0.0000000000,
                0.0000000000,
                0.0000000000,
                -0.0004144295,
            ],
        );
        let mut computed_eigvals = hessian.0.eigenvalues().unwrap();
        computed_eigvals
            .as_mut_slice()
            .sort_by(|a, b| b.partial_cmp(a).unwrap());
        for (c, ex) in computed_eigvals.iter().zip(expected_eigvals.iter()) {
            assert_relative_eq!(c, ex, epsilon = 1e-3);
        }
    }
}
