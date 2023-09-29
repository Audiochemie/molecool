use super::atom::Atom;
use coordinate_systems::DistanceTo;
use nalgebra::Vector3;
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

/// Function to retrive the distance Vector between  two atoms a and b
///
/// # Arguments
///
/// * `a` - First atom
/// * `b` - Second atom

pub fn get_distance_vec_between(a: &Atom, b: &Atom) -> Vector3<f32> {
    a.coordinates.distance_to(&b.coordinates)
}

/// Function to retrive angle between three build atoms.
///
/// # Arguments
///
/// * `l` - atom to the left of center atom.
/// * `c` - center atom.
/// * `r` - atom to the right of center atom.
///
pub fn get_angle_between(l: &Atom, c: &Atom, r: &Atom) -> f32 {
    let e_lc = c.coordinates.distance_to(&l.coordinates).normalize();
    let e_rc = c.coordinates.distance_to(&r.coordinates).normalize();

    e_lc.dot(&e_rc).acos()
}

/// Function to retrive all distance vectors for all atom pairs in the molecule.
///
/// # Arguments
///
/// * `mol` - molecule containing atoms.
pub fn get_distances(mol: &Molecule) -> Vec<Vector3<f32>> {
    let num_atms = mol.building_atoms.len();
    let mut distance_vector: Vec<Vector3<f32>> =
        Vec::with_capacity((num_atms * num_atms - num_atms) / 2);
    for i_atm in 0..mol.building_atoms.len() - 1 {
        let ref_atm_coor = &mol.building_atoms[i_atm].coordinates;
        for i_nxt_atm in i_atm + 1..mol.building_atoms.len() {
            let nxt_atm_coor = &mol.building_atoms[i_nxt_atm].coordinates;
            distance_vector.push(ref_atm_coor.distance_to(nxt_atm_coor));
        }
    }
    distance_vector
}

#[cfg(test)]
mod unit_tests {
    use std::f32::{consts::FRAC_1_PI, EPSILON};

    use super::{get_angle_between, get_distance_vec_between};
    use crate::atom::Atom;
    use approx::assert_relative_eq;
    use nalgebra::Point3;

    #[test]
    fn test_get_diff_vector() {
        let a = Atom {
            z_value: 6,
            atomic_mass: 12.0,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                0.0_f32, 0.0_f32, 0.0_f32,
            )),
        };
        let b = Atom {
            z_value: 6,
            atomic_mass: 12.0,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                0.0_f32,
                0.0_f32,
                2.845_112_f32,
            )),
        };
        let v = get_distance_vec_between(&a, &b);
        assert_eq!((0.0_f32, 0.0_f32, -2.845_112_f32), (v.x, v.y, v.z))
    }
    #[test]
    fn test_get_angle_between_atoms() {
        let l = Atom {
            z_value: 6,
            atomic_mass: 12.0,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                0.0_f32, 0.0_f32, 0.0_f32,
            )),
        };
        let c = Atom {
            z_value: 6,
            atomic_mass: 12.0,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                0.0_f32,
                0.0_f32,
                2.845_112_f32,
            )),
        };
        let r = Atom {
            z_value: 8,
            atomic_mass: 15.99,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                1.899_115_f32,
                0.0_f32,
                4.139_062_f32,
            )),
        };

        let d12  = get_distance_vec_between(&l, &c).norm();
        let e12  = get_distance_vec_between(&l, &c).normalize();
        let d13 = get_distance_vec_between(&l, &r).norm();
        assert_relative_eq!(d13, 4.553_95_f32, epsilon=0.00001);
        assert_relative_eq!(d12, 2.845_112_f32, epsilon=0.00001);
        assert_relative_eq!(e12.norm(), 1.0_f32, epsilon=0.00001);
        let a = get_angle_between(&l, &c, &r);
        assert_relative_eq!(a * 180.0_f32 * FRAC_1_PI, 124.268_31_f32, epsilon=0.000001 )
    }
}
