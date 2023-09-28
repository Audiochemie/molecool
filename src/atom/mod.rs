use super::pse_data::PSE_MASSES;
use coordinate_systems::Coordinates;
use qc_file_parsers::xyz::{xyzline::symbol::PSE_SYMBOLS, XyzLine};
/// Represents an Atom as found in the periodic system
pub struct Atom {
    /// The core charge
    pub z_value: usize,
    /// The atomic mass.
    pub atomic_mass: f32,
    /// The coordinates in space.
    pub coordinates: Coordinates,
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
