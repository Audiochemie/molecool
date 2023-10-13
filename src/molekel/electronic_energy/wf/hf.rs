//! This module implements the Hartree-Fock Hamiltonian (a.k.a) Fockian and it's eigenvalues
//! (orbital energies) and eigenfunctions (Molecular orbitals)

use std::io::BufRead;

use nalgebra::DMatrix;

pub struct Fockian {
    /// The nuclear-nuclear potential energy in E\(_h\) -> Constant shift to energy
    pub v_nuc_nuc: f64,
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

/// Function to read the overlap  matrix \(\mathbf{S}\) from a file.
/// The file shall contain the indices in S, \(\mu\) and \(\nu\),
/// as well as the float value representing the matrix elements \(S_{\mu\nu}\) of the overlap matrix.
/// # Arguments:
///  * `file` - The file to read the \(S_{\mu\nu}\) from.
///  * `n_mo` - The number of molecular orbitals
pub fn get_s_mn_from_file<I: BufRead>(file: &mut I, n_mo: usize) -> DMatrix<f64> {
    let mut s_mn = DMatrix::<f64>::zeros(n_mo, n_mo);
    let mut line_iter = file.lines().map(|l| l.unwrap());
    let first_line = line_iter.next().unwrap();
    let mut first_split = first_line.split_whitespace().take(3);
    let m: usize = first_split.next().unwrap().parse::<usize>().unwrap();
    let n: usize = first_split.next().unwrap().parse::<usize>().unwrap();
    let val: f64 = first_split.next().unwrap().parse::<f64>().unwrap();
    let mut index_shift: usize = 0;
    // If the first read index is one, we need to shift by -1 as Rust indexing starts at 0
    if m == 1 {
        index_shift = 1
    }
    s_mn[(m - index_shift, n - index_shift)] = val;
    for l in line_iter {
        let mut split = l.split_whitespace().take(3);
        let m: usize = split.next().unwrap().parse::<usize>().unwrap() - index_shift;
        let n: usize = split.next().unwrap().parse::<usize>().unwrap() - index_shift;
        let val: f64 = split.next().unwrap().parse::<f64>().unwrap();
        s_mn[(m, n)] = val;
    }
    s_mn.lower_triangle()
}

impl Fockian {
    fn new(v_nuc_nuc: f64) -> Self {
        Self { v_nuc_nuc }
    }
}
