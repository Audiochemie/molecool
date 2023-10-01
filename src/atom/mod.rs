use num::Float;
use super::pse_data::PSE_MASSES;
use coordinate_systems::cartesian::Cartesian;
use qc_file_parsers::xyz::{xyzline::symbol::PSE_SYMBOLS, XyzLine};

/// Represents an Atom as found in the periodic system
pub struct Atom<T>
where
    T: Float + std::fmt::Debug + std::str::FromStr + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug
{
    /// The core charge
    pub z_value: usize,
    /// The atomic mass.
    pub atomic_mass: f32,
    /// The coordinates in space.
    pub coordinates: Cartesian<T>,
}

impl<T> From<XyzLine<T>> for Atom<T> 
where
    T: Float + std::fmt::Debug + std::str::FromStr + 'static,
    <T as std::str::FromStr>::Err: std::fmt::Debug
{
    /// Conversion function of a line in an Xyz file (see there) to an Atom struct.
    fn from(value: XyzLine<T>) -> Self {
        match value {
            XyzLine::Numeric(n) => Self {
                coordinates: Cartesian(n.xyz),
                z_value: n.z_value,
                atomic_mass: PSE_MASSES[n.z_value],
            },
            XyzLine::Symbolic(s) => {
                let z_value = PSE_SYMBOLS.iter().position(|&sym| sym == s.symbol).unwrap();
                Self {
                    coordinates: Cartesian(s.xyz),
                    z_value,
                    atomic_mass: PSE_MASSES[z_value],
                }
            }
        }
    }
}
