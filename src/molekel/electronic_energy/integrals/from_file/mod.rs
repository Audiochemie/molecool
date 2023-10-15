use nalgebra::{DMatrix, DVector};
use qc_file_parsers::format_string::{parse_fortran_formatted_buf, ParsedValue};
use std::{io::BufRead, usize};

const PRECOMPUTE: usize = 1024 * 8;

/// Helper function to check if the thing parsed from a file can be used to index an Array.
/// # Arguments
///  * `parsed_thing` - Iterator over ParsedValue enums.
fn _check_valid_index(parsed_thing: Option<&ParsedValue>) -> usize
where
{
    match parsed_thing {
        Some(ParsedValue::In(i)) => *i as usize,
        _ => panic!("Given index is not an integer!"),
    }
}

/// Helper function to  build a composite index in the form of (m * m + m) / 2 + v
/// For m > v. First term is start of block, v is shift inside block
/// Uses precomputed index map for values <= 8192.
/// Example:
///
/// # Arguments
///  * `ind1` - first index
///  * `ind2` - Second index
///
fn _build_composite_index(ind1: usize, ind2: usize) -> usize
where
{
    let mut precomputed_offsets: [usize; PRECOMPUTE] = [0; PRECOMPUTE];
    for (i, e) in precomputed_offsets
        .iter_mut()
        .enumerate()
        .take(PRECOMPUTE)
        .skip(1)
    {
        *e = (i * i + i) / 2
    }
    if ind1 > 8192 || ind2 > 8192 {
        unimplemented!("Computation still needs to be implemented. Sorry.")
    }
    match ind1 > ind2 {
        true => precomputed_offsets[ind1] + ind2,
        false => precomputed_offsets[ind2] + ind1,
    }
}

/// Function to get the repulsion potential of the nuclei V\(_{NM}\) from a file.
/// The file shall only contain the numeric value of V\(_{NM}\) as a floating point number.
///
/// # Arguments:
///
///  * `file` - The file to read V\(_{NM}\) from.
///
/// # Panics:
///  
/// If nothing was parsed successful.
///
///  If the file contains anything else than a single floating point number.
pub fn get_v_nuc_nuc_from_fortran_format_file<I: BufRead>(file: &mut I) -> f64 {
    let parsed = parse_fortran_formatted_buf(file).unwrap();
    match parsed.first() {
        Some(p) => match p.iter().next() {
            Some(ParsedValue::Fl(f)) => *f,
            _ => {
                panic!("Parsed value is not a float. Assuming this is wrong for V_{{nuc,nuc}}!")
            }
        },
        None => {
            panic!("Nothing could be parsed from file")
        }
    }
}

/// Function to read the overlap  matrix \(\mathbf{S}\) from a file.
/// The file shall contain the indices in S, \(\mu\) and \(\nu\),
/// as well as the float value representing the matrix elements \(S_{\mu\nu}\) of the overlap matrix.
/// Therefore, the Fortran format shall be something like, `(2i4,f21.21)`
/// As the matrix is symmetric only a view of the lower diagonal is returned.
/// # Arguments:
///
///  * `file` - The file to read the \(S_{\mu\nu}\) from. MUST start with a fortran format string.
///  * `n_mo` - The number of molecular orbitals
///
pub fn get_s_mn_from_fortran_format_file<I: BufRead>(file: &mut I, n_mo: usize) -> DMatrix<f64> {
    let mut s_mn = DMatrix::<f64>::zeros(n_mo, n_mo);
    let parsed_data: Vec<Vec<ParsedValue>> = parse_fortran_formatted_buf(file).unwrap();
    let mut parsed_lines_iter = parsed_data.iter();
    let mut indices = parsed_lines_iter.next().unwrap().iter();
    let m = _check_valid_index(indices.next());
    let n = _check_valid_index(indices.next());
    let index_shift = match m {
        1 => 1,
        0 => 0,
        _ => panic!("First index should be one or zero! Is {}", m),
    };
    let val = match parsed_lines_iter.next().unwrap().iter().next().unwrap() {
        ParsedValue::Fl(f) => *f,
        _ => panic!("First value is not a float!"),
    };
    s_mn[(m - index_shift, n - index_shift)] = val;
    let mut mn: [usize; 2] = [usize::MAX, usize::MAX];
    let mut val: f64 = std::f64::NEG_INFINITY;
    for what_remains in parsed_lines_iter {
        for (index, iterable) in what_remains.iter().enumerate() {
            if let ParsedValue::In(i) = iterable {
                // Here it is CRUCIAL to we only have three entries in the order index, index, value
                mn[index] = *i as usize
            } else if let ParsedValue::Fl(f) = iterable {
                val = *f
            }
        }
        s_mn[(mn[0] - index_shift, mn[1] - index_shift)] = val;
    }
    s_mn.lower_triangle()
}

/// Function to read the one-particle kinetic energy matrix \(\mathbf{T}\) from a file.
/// The file shall contain the indices in T, \(\mu\) and \(\nu\),
/// as well as the float value representing the matrix elements \(T_{\mu\nu}\) of the kinetic
/// energy matrix.
/// Therefore, the Fortran format shall be something like, `(2i4,f21.21)`
/// As the matrix is symmetric only a view of the lower diagonal is returned.
/// # Arguments:
///
///  * `file` - The file to read the \(S_{\mu\nu}\) from. MUST start with a fortran format string.
///  * `n_mo` - The number of molecular orbitals
///
pub fn get_t_mn_from_fortran_format_file<I: BufRead>(file: &mut I, n_mo: usize) -> DMatrix<f64> {
    let mut t_mn = DMatrix::<f64>::zeros(n_mo, n_mo);
    let parsed_data: Vec<Vec<ParsedValue>> = parse_fortran_formatted_buf(file).unwrap();
    let mut parsed_lines_iter = parsed_data.iter();
    let mut indices = parsed_lines_iter.next().unwrap().iter();
    let m = _check_valid_index(indices.next());
    let n = _check_valid_index(indices.next());
    let index_shift = match m {
        1 => 1,
        0 => 0,
        _ => panic!("First index should be one or zero! Is {}", m),
    };
    let val = match parsed_lines_iter.next().unwrap().iter().next().unwrap() {
        ParsedValue::Fl(f) => *f,
        _ => panic!("First value is not a float!"),
    };
    t_mn[(m - index_shift, n - index_shift)] = val;
    let mut mn: [usize; 2] = [usize::MAX, usize::MAX];
    let mut val: f64 = std::f64::NEG_INFINITY;
    for what_remains in parsed_lines_iter {
        for (index, iterable) in what_remains.iter().enumerate() {
            if let ParsedValue::In(i) = iterable {
                // Here it is CRUCIAL to we only have three entries in the order index, index, value
                mn[index] = *i as usize
            } else if let ParsedValue::Fl(f) = iterable {
                val = *f
            }
        }
        t_mn[(mn[0] - index_shift, mn[1] - index_shift)] = val;
    }
    t_mn.lower_triangle()
}

/// Function to read the one-particle nucleus-particle interaction matrix \(\mathbf{V_{Ne}}\) from a file.
/// The file shall contain the indices in V, \(\mu\) and \(\nu\),
/// as well as the float value representing the matrix elements \(V_{\mu\nu}\) of the matrix.
/// Therefore, the Fortran format shall be something like, `(2i4,f21.21)`
/// As the matrix is symmetric only a view of the lower diagonal is returned.
/// # Arguments:
///
///  * `file` - The file to read the \(V_{\mu\nu}\) from. MUST start with a fortran format string.
///  * `n_mo` - The number of molecular orbitals
///
pub fn get_v_mn_from_fortran_format_file<I: BufRead>(file: &mut I, n_mo: usize) -> DMatrix<f64> {
    let mut v_mn = DMatrix::<f64>::zeros(n_mo, n_mo);
    let parsed_data: Vec<Vec<ParsedValue>> = parse_fortran_formatted_buf(file).unwrap();
    let mut parsed_lines_iter = parsed_data.iter();
    let mut indices = parsed_lines_iter.next().unwrap().iter();
    let m = _check_valid_index(indices.next());
    let n = _check_valid_index(indices.next());
    let index_shift = match m {
        1 => 1,
        0 => 0,
        _ => panic!("First index should be one or zero! Is {}", m),
    };
    let val = match parsed_lines_iter.next().unwrap().iter().next().unwrap() {
        ParsedValue::Fl(f) => *f,
        _ => panic!("First value is not a float!"),
    };
    v_mn[(m - index_shift, n - index_shift)] = val;
    let mut mn: [usize; 2] = [usize::MAX, usize::MAX];
    let mut val: f64 = std::f64::NEG_INFINITY;
    for what_remains in parsed_lines_iter {
        for (index, iterable) in what_remains.iter().enumerate() {
            if let ParsedValue::In(i) = iterable {
                // Here it is CRUCIAL to we only have three entries in the order index, index, value
                mn[index] = *i as usize
            } else if let ParsedValue::Fl(f) = iterable {
                val = *f
            }
        }
        v_mn[(mn[0] - index_shift, mn[1] - index_shift)] = val;
    }
    v_mn.lower_triangle()
}

/// Function to read the two-particle coulomb interaction matrix \(\mathbf{V_{ee}}\) from a file.
/// The file shall contain the indices in V, \(\mu\), \(\nu\), \(\sigma\) and \(\tau\)
/// as well as the float value representing the matrix elements \(V_{\mu\nu\sigma\tau}\) of the matrix.
/// Therefore, the Fortran format shall be something like, `(4i4,f21.21)`
/// At the moment it is assumed that only the permutationally unique indices are given,
/// as the Coulomb-Integrals in the LCAO ansatz have an eight-fold permutational symmetry for real
/// orbitals. This is rooted in the symmetry of the Columb-potential, as
/// \(r_{12}^{-1}= r_{21}^{-1}\), and the fact that the particles are indistinguishable.
/// Neglecting symmstry there are at most M^2 index pairs \(\left(\mu\nu\right)\) per particle,
/// where M is the number of available one-particle basis states, resulting in M^4 indices in total.
/// Considering the operator symmetry brings the sets of index pairs \(A, B\) per particle down to
/// a size of unique pairs including repetitions \(N = \frac{M * (M + 1)}{2}\) each.
/// \(\forall a_i \in A\) exist N \( b_i \in B\)
///
/// # Arguments:
///
///  * `file` - The file to read the \(V_{\mu\nu}\) from. MUST start with a fortran format string.
///  * `n_mo` - The number of molecular orbitals
///
pub fn get_v_ee_from_fortran_format_file<I: BufRead>(file: &mut I, n_mo: usize) -> DVector<f64> {
    let pairs = (n_mo * n_mo + n_mo) / 2;
    let max_index = (pairs * pairs + pairs) / 2;
    let mut vee = DVector::zeros(max_index);
    let parsed_data: Vec<Vec<ParsedValue>> = parse_fortran_formatted_buf(file).unwrap();
    let mut parsed_lines_iter = parsed_data.iter();
    let mut indices = parsed_lines_iter.next().unwrap().iter();
    let mut m: usize = _check_valid_index(indices.next());
    let index_shift = match m {
        1 => 1,
        0 => 0,
        _ => panic!("First index should be one or zero! Is {}", m),
    };
    m -= index_shift;
    let mut n: usize = _check_valid_index(indices.next()) - index_shift;
    let mut s: usize = _check_valid_index(indices.next()) - index_shift;
    let mut t: usize = _check_valid_index(indices.next()) - index_shift;
    let val = match parsed_lines_iter.next().unwrap().iter().next().unwrap() {
        ParsedValue::Fl(f) => *f,
        _ => panic!("First value is not a float!"),
    };
    let mut mn: usize = _build_composite_index(m, n);
    let mut st: usize = _build_composite_index(s, t);
    let mut mnst: usize = _build_composite_index(mn, st);
    vee[mnst] = val;
    let mut val: f64 = std::f64::NEG_INFINITY;
    for what_remains in parsed_lines_iter {
        let mut this_line = what_remains.iter();
        match this_line.next() {
            Some(ParsedValue::In(i)) => {
                m = *i as usize - index_shift;
                n = _check_valid_index(this_line.next()) - index_shift;
                s = _check_valid_index(this_line.next()) - index_shift;
                t = _check_valid_index(this_line.next()) - index_shift;
                mn = _build_composite_index(m, n);
                st = _build_composite_index(s, t);
                mnst = _build_composite_index(mn, st);
            }
            Some(ParsedValue::Fl(f)) => val = *f,
            _ => {
                panic!("Data not valid!")
            }
        };
        vee[mnst] = val;
    }
    vee
}

#[cfg(test)]
mod unit_tests {
    use super::_build_composite_index;
    #[test]
    fn test_composite_index() {
        let m = 9;
        let n = 10;
        assert_eq!(_build_composite_index(m, n), 64);
        let m = 14;
        let n = 13;
        assert_eq!(_build_composite_index(m, n), 118);
        let mn = 105;
        let st = 32;
        assert_eq!(_build_composite_index(mn, st), 5597);
    }
}
