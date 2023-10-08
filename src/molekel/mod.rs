use super::atom::Atom;
use nalgebra::{Point3, Vector3};
use num::{traits::AsPrimitive, zero, Float};

use coordinate_systems::DistanceTo;
use qc_file_parsers::xyz::Xyz;
/// Represents a molecule in the _atoms in molecules_ sense.
pub struct Molecule<T>
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + nalgebra::ClosedDiv
        + num::cast::AsPrimitive<T>
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    /// Molecules have a total molar mass.
    pub molar_mass: T,
    /// Molecules are commonly interpreted as a collection of building atoms.
    pub building_atoms: Vec<Atom<T>>,
}

impl<T> From<Xyz<T>> for Molecule<T>
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + num::cast::AsPrimitive<T>
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + nalgebra::ClosedDiv
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
    f32: num::traits::AsPrimitive<T>,
    f64: num::traits::AsPrimitive<T>,
{
    /// Cast Xyz struct into a Molecule.
    fn from(value: Xyz<T>) -> Self {
        let mut atoms = Vec::with_capacity(value.number_of_atoms);
        value
            .lines
            .iter()
            .for_each(|xyz_line| atoms.push(Atom::from(xyz_line.clone())));
        let molar_mass: T = atoms
            .iter()
            .fold(0.0.as_(), |at, at2| at + at2.atomic_mass.as_());
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
pub fn get_distance_vec_between<T>(a: &Atom<T>, b: &Atom<T>) -> Vector3<T>
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + std::ops::SubAssign
        + nalgebra::ComplexField<RealField = T>
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
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
pub fn get_angle_between<T>(l: &Atom<T>, c: &Atom<T>, r: &Atom<T>) -> T
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + std::ops::SubAssign
        + nalgebra::ComplexField<RealField = T>
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    let e_lc = c.coordinates.distance_to(&l.coordinates).normalize();
    let e_rc = c.coordinates.distance_to(&r.coordinates).normalize();
    num::Float::acos(e_lc.dot(&e_rc))
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
pub fn get_oop<T>(l: &Atom<T>, c: &Atom<T>, r: &Atom<T>, h: &Atom<T>) -> T
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + std::ops::SubAssign
        + nalgebra::ComplexField<RealField = T>
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    let e_lc = c.coordinates.distance_to(&l.coordinates).normalize();
    let e_rc = c.coordinates.distance_to(&r.coordinates).normalize();
    let e_hc = c.coordinates.distance_to(&h.coordinates).normalize();
    let denom = num::Float::sin(get_angle_between(l, c, r));
    let n_vec_lcr = e_lc.cross(&e_rc);
    if n_vec_lcr.norm() == zero() {
        // TODO Change Return to Result and handle properly.
        panic!("Normal vector is zero vector => e_lc || e_rc: unit vector are parallel.");
    }
    num::Float::asin((n_vec_lcr).dot(&e_hc).div(denom))
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
pub fn get_tors<T>(c1: &Atom<T>, b1: &Atom<T>, b2: &Atom<T>, c2: &Atom<T>) -> T
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + std::ops::SubAssign
        + nalgebra::ComplexField<RealField = T>
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    let u1 = get_distance_vec_between(c1, b1);
    let u2 = get_distance_vec_between(b1, b2);
    let u3 = get_distance_vec_between(b2, c2);
    let cross_1 = u1.cross(&u2);
    let cross_2 = u2.cross(&u3);
    let denom = cross_1.norm() * cross_2.norm();
    let dot_denom = (cross_1).dot(&cross_2).div(denom);
    num::Float::acos(dot_denom)
}

/// Function to retrive all distance vectors for all atom pairs in the molecule.
///
/// # Arguments
///
/// * `mol` - molecule containing atoms.
pub fn get_distances<T>(mol: &Molecule<T>) -> Vec<Vector3<T>>
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + std::ops::SubAssign
        + nalgebra::ComplexField<RealField = T>
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + num::cast::AsPrimitive<T>
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
{
    let num_atms = mol.building_atoms.len();
    let mut distance_vector: Vec<Vector3<T>> =
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

/// Computes the centre of mass (com) of the Molecule.
/// com is computed as the sum of all mass weighted coordinates.
/// # Arguments
///
///  *  mol - Reference to Molecule struct.
///
pub fn get_centre_of_mass<T>(mol: &Molecule<T>) -> Point3<T>
where
    T: Float
        + std::fmt::Debug
        + std::str::FromStr
        + num::cast::AsPrimitive<T>
        + nalgebra::ClosedMul
        + nalgebra::ClosedAdd
        + nalgebra::ClosedDiv
        + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug,
    f32: num::traits::AsPrimitive<T>,
    f64: num::traits::AsPrimitive<T>,
{
    mol.building_atoms.iter().fold(Point3::<T>::origin(), |mut o, at| {
        let mw = at.coordinates.0.map(|e| e*at.atomic_mass.as_()/mol.molar_mass);
        o.x +=  mw.x;
        o.y +=  mw.y;
        o.z +=  mw.z;
        o
    })
}

#[cfg(test)]
mod unit_tests {
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
