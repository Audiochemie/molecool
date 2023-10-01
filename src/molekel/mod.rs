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

/// Function to retrive out-of-plane angle between four building atoms.
///
/// # Arguments
///
/// * `l` - atom to the left of center atom.
/// * `c` - center atom.
/// * `r` - atom to the right of center atom.
/// * `h` - _hovering_ atom.
///
pub fn get_oop(l: &Atom, c: &Atom, r: &Atom, h: &Atom) -> f32 {
    let e_lc = c.coordinates.distance_to(&l.coordinates).normalize();
    let e_rc = c.coordinates.distance_to(&r.coordinates).normalize();
    let e_hc = c.coordinates.distance_to(&h.coordinates).normalize();
    let denom = get_angle_between(l, c, r).sin();
    let n_vec_lcr = e_lc.cross(&e_rc);
    if n_vec_lcr.norm() == 0.0_f32 {
        // TODO Change Return to Result and handle properly.
        panic!("Normal vector is zero vector => e_lc || e_rc: unit vector are parallel.");
    }
    (n_vec_lcr * 1.0_f32 / denom).dot(&e_hc).asin()
}

/// Function to retrive torsional angle between four building atoms.
///
/// # Arguments
///
/// * `c1` - atom connected to b1
/// * `b1` - bonding atom1
/// * `b2` - bonding atom2
/// * `c2` - atom connected to b2
///
pub fn get_tors(c1: &Atom, b1: &Atom, b2: &Atom, c2: &Atom) -> f32 {
    let u1 = get_distance_vec_between(c1, b1);
    let u2 = get_distance_vec_between(b1, b2);
    let u3 = get_distance_vec_between(b2, c2);
    let cross_1 = u1.cross(&u2);
    let cross_2 = u2.cross(&u3);
    let denom = 1.0_f32 / (cross_1.norm() * cross_2.norm());
    let dot_denom = (cross_1).dot(&cross_2) * denom;
    dot_denom.acos()
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
    use std::f32::consts::FRAC_1_PI;

    use super::{get_angle_between, get_distance_vec_between, get_oop};
    use crate::{atom::Atom, molekel::get_tors};
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

        let d12 = get_distance_vec_between(&l, &c).norm();
        let e12 = get_distance_vec_between(&l, &c).normalize();
        let d13 = get_distance_vec_between(&l, &r).norm();
        assert_relative_eq!(d13, 4.553_95_f32, epsilon = 0.00001);
        assert_relative_eq!(d12, 2.845_112_f32, epsilon = 0.00001);
        assert_relative_eq!(e12.norm(), 1.0_f32, epsilon = 0.00001);
        let a = get_angle_between(&l, &c, &r);
        assert_relative_eq!(a.to_degrees(), 124.268_31_f32, epsilon = 0.000001)
    }

    #[test]
    fn test_oop() {
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
        let h = Atom {
            z_value: 1,
            atomic_mass: 1.09,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                -1.894_048_f32,
                0.0_f32,
                3.747_688_7_f32,
            )),
        };
        assert_relative_eq!(0.0_f32, get_oop(&l, &c, &r, &h))
    }
    #[test]
    fn test_torsional_angle() {
        let c2 = Atom {
            z_value: 6,
            atomic_mass: 12.0,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                0.0_f32, 0.0_f32, 0.0_f32,
            )),
        };
        let b2 = Atom {
            z_value: 6,
            atomic_mass: 12.0,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                0.0_f32,
                0.0_f32,
                2.845_112_f32,
            )),
        };
        let b1 = Atom {
            z_value: 8,
            atomic_mass: 15.99,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                1.899_115_f32,
                0.0_f32,
                4.139_062_f32,
            )),
        };
        let c1 = Atom {
            z_value: 1,
            atomic_mass: 1.09,
            coordinates: coordinate_systems::cartesian::Cartesian(Point3::new(
                -1.894_048_f32,
                0.0_f32,
                3.747_688_7_f32,
            )),
        };
        assert_eq!(180.0_f32, get_tors(&c1, &b1, &b2, &c2).to_degrees());
        assert_eq!(180.0_f32, get_tors(&c2, &b2, &b1, &c1).to_degrees())
    }
}
