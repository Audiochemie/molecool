use nalgebra::Point3;
use qc_file_parsers::{xyzline::symbol::PSE_SYMBOLS, XyzLine};

pub enum Coordinates {
    Cartesian(Point3<f32>),
    Internal,
}


pub struct Atom {
    pub z_value: usize,
    pub atomic_mass: f32,
    pub coordinates: Coordinates,
}

pub struct Molecule {
    pub molar_mass: f32,
    pub atoms: Atom,
}

impl From<XyzLine> for Atom {
    fn from(value: XyzLine) -> Self {
        match value {
            XyzLine::Numeric(n) => Self {
                coordinates: Coordinates::Cartesian(n.xyz),
                z_value: n.z_value,
                atomic_mass: 0.0_f32,
            },
            XyzLine::Symbolic(s) => Self {
                coordinates: Coordinates::Cartesian(s.xyz),
                z_value: PSE_SYMBOLS.iter().position(|&sym| sym == s.symbol).unwrap(),
                atomic_mass: 0.0_f32,
            },
        }
    }
}
