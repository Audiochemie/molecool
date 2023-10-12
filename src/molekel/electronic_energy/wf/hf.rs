//! This module implements the Hartree-Fock Hamiltonian (a.k.a) Fockian and it's eigenvalues
//! (orbital energies) and eigenfunctions (Molecular orbitals)

use std::io::BufRead;

struct Fockian {
    /// The nuclear-nuclear potential energy in E\(_h\) -> Constant shift to energy
    v_nuc_nuc: f64,
}

/// Function to get the repulsion potential of the nuclei V\(_{NM}\) from a file.
/// The file shall only contain the numeric value of V\(_{NM}\) as a floating point number.
/// # Arguments
///  * `file` - The file to read V\(_{NM}\) from.
///
/// # Panics
///  If the file contains anything else than a single floating point number.
pub fn get_v_nuc_nuc_from_file<I: BufRead>(file: &mut I) -> f64 {
    let buf = &mut String::new();
    file.read_to_string(buf).unwrap();
    buf.trim().parse::<f64>().unwrap()
}
