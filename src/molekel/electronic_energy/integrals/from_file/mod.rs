use nalgebra::DMatrix;
use qc_file_parsers::format_string::{self, parse_fortran_formatted_buf};
use std::{io::BufRead, usize};

/// Function to get the repulsion potential of the nuclei V\(_{NM}\) from a file.
/// The file shall only contain the numeric value of V\(_{NM}\) as a floating point number.
///
/// # Arguments:
///
///  * `file` - The file to read V\(_{NM}\) from.
///
/// # Panics:
///
///  If the file contains anything else than a single floating point number.
pub fn get_v_nuc_nuc_from_fortran_format_file<I: BufRead>(file: &mut I) -> f64 {
    let parsed = parse_fortran_formatted_buf(file).unwrap();
    match parsed.first() {
        Some(p) => match p.iter().next() {
            Some(format_string::ParsedValue::Fl(f)) => *f,
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
/// # Arguments:
///
///  * `file` - The file to read the \(S_{\mu\nu}\) from. MUST start with a fortran format string.
///  * `n_mo` - The number of molecular orbitals
///
pub fn get_s_mn_from_fortran_format_file<I: BufRead>(file: &mut I, n_mo: usize) -> DMatrix<f64> {
    let mut s_mn = DMatrix::<f64>::zeros(n_mo, n_mo);
    let parsed_data: Vec<Vec<format_string::ParsedValue>> =
        format_string::parse_fortran_formatted_buf(file).unwrap();
    let mut parsed_lines_iter = parsed_data.iter();
    let mut indices = parsed_lines_iter.next().unwrap().iter();
    let m: usize = match indices.next() {
        Some(format_string::ParsedValue::In(i)) => *i as usize,
        None => {
            panic!("There is no ParsedValue here, i.e. second entry is not an integer! No Index => No matrix. :-(")
        }
        _ => {
            panic!("Given index is not an integer!")
        }
    };
    let n: usize = match indices.next() {
        Some(format_string::ParsedValue::In(i)) => *i as usize,
        None => {
            panic!(
                "There is no ParsedValue here, i.e. second entry is not an integer! No Index => No matrix. :-("
            )
        }
        _ => {
            panic!("Given index is not an integer!")
        }
    };
    let index_shift = match m {
        1 => 1,
        0 => 0,
        _ => panic!("First index should be one or zero! Is {}", m),
    };
    let val = match parsed_lines_iter.next().unwrap().iter().next().unwrap() {
        format_string::ParsedValue::Fl(f) => *f,
        _ => panic!("First value is not a float!"),
    };
    s_mn[(m - index_shift, n - index_shift)] = val;
    let mut mn: [usize; 2] = [usize::MAX, usize::MAX];
    let mut val: f64 = std::f64::NEG_INFINITY;
    for what_remains in parsed_lines_iter {
        for (index, iterable) in what_remains.iter().enumerate() {
            if let format_string::ParsedValue::In(i) = iterable {
                // Here it is CRUCIAL to we only have three entries in the order index, index, value
                mn[index] = *i as usize
            } else if let format_string::ParsedValue::Fl(f) = iterable {
                val = *f
            }
        }
        s_mn[(mn[0] - index_shift, mn[1] - index_shift)] = val;
    }
    s_mn.lower_triangle()
}

/// Function to read the one-particle kinetic energy matrix \(\mathbf{T}\) from a file.
/// The file shall contain the indices in T, \(\mu\) and \(\nu\),
/// as well as the float value representing the matrix elements \(T_{\mu\nu}\) of the overlap matrix.
/// Therefore, the Fortran format shall be something like, `(2i4,f21.21)`
/// # Arguments:
///
///  * `file` - The file to read the \(S_{\mu\nu}\) from. MUST start with a fortran format string.
///  * `n_mo` - The number of molecular orbitals
///
pub fn get_t_mn_from_fortran_format_file<I: BufRead>(file: &mut I, n_mo: usize) -> DMatrix<f64> {
    let mut t_mn = DMatrix::<f64>::zeros(n_mo, n_mo);
    let parsed_data: Vec<Vec<format_string::ParsedValue>> =
        format_string::parse_fortran_formatted_buf(file).unwrap();
    let mut parsed_lines_iter = parsed_data.iter();
    let mut indices = parsed_lines_iter.next().unwrap().iter();
    let m: usize = match indices.next() {
        Some(format_string::ParsedValue::In(i)) => *i as usize,
        None => {
            panic!("There is no ParsedValue here, i.e. second entry is not an integer! No Index => No matrix. :-(")
        }
        _ => {
            panic!("Given index is not an integer!")
        }
    };
    let n: usize = match indices.next() {
        Some(format_string::ParsedValue::In(i)) => *i as usize,
        None => {
            panic!(
                "There is no ParsedValue here, i.e. second entry is not an integer! No Index => No matrix. :-("
            )
        }
        _ => {
            panic!("Given index is not an integer!")
        }
    };
    let index_shift = match m {
        1 => 1,
        0 => 0,
        _ => panic!("First index should be one or zero! Is {}", m),
    };
    let val = match parsed_lines_iter.next().unwrap().iter().next().unwrap() {
        format_string::ParsedValue::Fl(f) => *f,
        _ => panic!("First value is not a float!"),
    };
    t_mn[(m - index_shift, n - index_shift)] = val;
    let mut mn: [usize; 2] = [usize::MAX, usize::MAX];
    let mut val: f64 = std::f64::NEG_INFINITY;
    for what_remains in parsed_lines_iter {
        for (index, iterable) in what_remains.iter().enumerate() {
            if let format_string::ParsedValue::In(i) = iterable {
                // Here it is CRUCIAL to we only have three entries in the order index, index, value
                mn[index] = *i as usize
            } else if let format_string::ParsedValue::Fl(f) = iterable {
                val = *f
            }
        }
        t_mn[(mn[0] - index_shift, mn[1] - index_shift)] = val;
    }
    t_mn.lower_triangle()
}

/// Function to read the one-particle nucleus-particle interaction matrix \(\mathbf{V_{Ne}}\) from a file.
/// The file shall contain the indices in V, \(\mu\) and \(\nu\),
/// as well as the float value representing the matrix elements \(V_{\mu\nu}\) of the overlap matrix.
/// Therefore, the Fortran format shall be something like, `(2i4,f21.21)`
/// # Arguments:
///
///  * `file` - The file to read the \(V_{\mu\nu}\) from. MUST start with a fortran format string.
///  * `n_mo` - The number of molecular orbitals
///
pub fn get_v_mn_from_fortran_format_file<I: BufRead>(file: &mut I, n_mo: usize) -> DMatrix<f64> {
    let mut v_mn = DMatrix::<f64>::zeros(n_mo, n_mo);
    let parsed_data: Vec<Vec<format_string::ParsedValue>> =
        format_string::parse_fortran_formatted_buf(file).unwrap();
    let mut parsed_lines_iter = parsed_data.iter();
    let mut indices = parsed_lines_iter.next().unwrap().iter();
    let m: usize = match indices.next() {
        Some(format_string::ParsedValue::In(i)) => *i as usize,
        None => {
            panic!("There is no ParsedValue here, i.e. second entry is not an integer! No Index => No matrix. :-(")
        }
        _ => {
            panic!("Given index is not an integer!")
        }
    };
    let n: usize = match indices.next() {
        Some(format_string::ParsedValue::In(i)) => *i as usize,
        None => {
            panic!(
                "There is no ParsedValue here, i.e. second entry is not an integer! No Index => No matrix. :-("
            )
        }
        _ => {
            panic!("Given index is not an integer!")
        }
    };
    let index_shift = match m {
        1 => 1,
        0 => 0,
        _ => panic!("First index should be one or zero! Is {}", m),
    };
    let val = match parsed_lines_iter.next().unwrap().iter().next().unwrap() {
        format_string::ParsedValue::Fl(f) => *f,
        _ => panic!("First value is not a float!"),
    };
    v_mn[(m - index_shift, n - index_shift)] = val;
    let mut mn: [usize; 2] = [usize::MAX, usize::MAX];
    let mut val: f64 = std::f64::NEG_INFINITY;
    for what_remains in parsed_lines_iter {
        for (index, iterable) in what_remains.iter().enumerate() {
            if let format_string::ParsedValue::In(i) = iterable {
                // Here it is CRUCIAL to we only have three entries in the order index, index, value
                mn[index] = *i as usize
            } else if let format_string::ParsedValue::Fl(f) = iterable {
                val = *f
            }
        }
        v_mn[(mn[0] - index_shift, mn[1] - index_shift)] = val;
    }
    v_mn.lower_triangle()
}

