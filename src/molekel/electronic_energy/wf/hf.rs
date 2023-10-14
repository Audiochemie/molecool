//! This module implements the Hartree-Fock Hamiltonian (a.k.a) Fockian and it's eigenvalues
//! (orbital energies) and eigenfunctions (Molecular orbitals)

use std::io::BufRead;

use nalgebra::DMatrix;
use qc_file_parsers::format_string;

pub struct Fockian {
    /// The nuclear-nuclear potential energy in E\(_h\) -> Constant shift to energy
    pub v_nuc_nuc: f64,
}


impl Fockian {
    fn new(v_nuc_nuc: f64) -> Self {
        Self { v_nuc_nuc }
    }
}
