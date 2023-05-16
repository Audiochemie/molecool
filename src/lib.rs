use crate::pse_data::PSE_MASSES;
use nalgebra::Point3;
use qc_file_parsers::{xyzline::symbol::PSE_SYMBOLS, XyzLine};

pub mod pse_data;

/// Collects different coordinate types commonly used in quantum chemistry.
#[derive(Debug, PartialEq)]
pub enum Coordinates {
    /// Cartesian coordinates are a collection of points in Euclidean space.
    Cartesian(Point3<f32>),
    /// Internal coordinates are displacement coordinates, such as distances, angles, ...
    Internal,
}

/// Represents an Atom as found in the periodic system
pub struct Atom {
    /// The core charge
    pub z_value: usize,
    /// The atomic mass.
    pub atomic_mass: f32,
    /// The coordinates in space.
    pub coordinates: Coordinates,
}

/// Represents a molecule in the _atoms in molecules_ sense.
pub struct Molecule {
    /// Molecules have a total molar mass.
    pub molar_mass: f32,
    /// Molecules are commonly interpreted as a collection of building atoms.
    pub building_atoms: Vec<Atom>,
}

impl From<XyzLine> for Atom {
    /// Conversion function of a line in an Xyz file (see there) to an Atom struct.
    fn from(value: XyzLine) -> Self {
        match value {
            XyzLine::Numeric(n) => Self {
                coordinates: Coordinates::Cartesian(n.xyz),
                z_value: n.z_value,
                atomic_mass: PSE_MASSES[n.z_value],
            },
            XyzLine::Symbolic(s) => {
                let z_value = PSE_SYMBOLS.iter().position(|&sym| sym == s.symbol).unwrap();
                Self {
                    coordinates: Coordinates::Cartesian(s.xyz),
                    z_value,
                    atomic_mass: PSE_MASSES[z_value],
                }
            }
        }
    }
}
