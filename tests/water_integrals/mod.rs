use crate::file_setup;
use approx::assert_relative_eq;
use molecool::molekel::electronic_energy::integrals::from_file;
use nalgebra::DMatrix;

#[test] fn test_get_v_nuc_nuc_from_file() {
    let mut test_file = file_setup::setup_v_nuc_nuc_water().unwrap();
    let test_v_nuc_nuc = from_file::get_v_nuc_nuc_from_fortran_format_file(&mut test_file);
    assert_eq!(test_v_nuc_nuc, 8.002_367_061_810_450)
}

#[test]
#[should_panic(expected = "Given index is not an integer!")]
fn test_get_smn_from_faulty_file() {
    let mut test_file = file_setup::setup_s_mn_water_faulty().unwrap();
    from_file::get_s_mn_from_fortran_format_file(&mut test_file, 14);
}

#[test]
#[should_panic(expected = "Given index is not an integer!")]
fn test_get_smn_from_faulty_file2() {
    let mut test_file = file_setup::setup_s_mn_water_faulty2().unwrap();
    from_file::get_s_mn_from_fortran_format_file(&mut test_file, 14);
}

#[test]
fn test_get_smn_from_file() {
    let mut test_file = file_setup::setup_s_mn_water().unwrap();
    let test_s_mn = from_file::get_s_mn_from_fortran_format_file(&mut test_file, 14);
    let expected_s_mn = [
        (1, 1, 1.000000000000000),
        (2, 1, 0.838055657090899),
        (2, 2, 1.000000000000000),
        (3, 1, 0.339351358965910),
        (3, 2, 0.432187282413285),
        (3, 3, 1.000000000000000),
        (4, 1, 0.160966484410815),
        (4, 2, 0.194384266485776),
        (4, 3, 0.776448695132864),
        (4, 4, 1.000000000000000),
        (5, 1, 0.000000000000000),
        (5, 2, -0.000000000000000),
        (5, 3, 0.000000000000000),
        (5, 4, 0.000000000000000),
        (5, 5, 1.000000000000000),
        (6, 1, -0.000000000000000),
        (6, 2, -0.000000000000000),
        (6, 3, -0.000000000000000),
        (6, 4, 0.000000000000000),
        (6, 5, 0.000000000000000),
        (6, 6, 1.000000000000000),
        (7, 1, 0.000000000000000),
        (7, 2, 0.000000000000000),
        (7, 3, 0.000000000000000),
        (7, 4, -0.000000000000000),
        (7, 5, -0.000000000000000),
        (7, 6, 0.000000000000000),
        (7, 7, 1.000000000000000),
        (8, 1, 0.000000000000000),
        (8, 2, 0.000000000000000),
        (8, 3, 0.000000000000000),
        (8, 4, 0.000000000000000),
        (8, 5, 0.505998479748195),
        (8, 6, 0.000000000000000),
        (8, 7, 0.000000000000000),
        (8, 8, 1.000000000000000),
        (9, 1, 0.000000000000000),
        (9, 2, 0.000000000000000),
        (9, 3, 0.000000000000000),
        (9, 4, 0.000000000000000),
        (9, 5, 0.000000000000000),
        (9, 6, 0.505998479748195),
        (9, 7, 0.000000000000000),
        (9, 8, 0.000000000000000),
        (9, 9, 1.000000000000000),
        (10, 1, 0.000000000000000),
        (10, 2, 0.000000000000000),
        (10, 3, 0.000000000000000),
        (10, 4, 0.000000000000000),
        (10, 5, 0.000000000000000),
        (10, 6, 0.000000000000000),
        (10, 7, 0.505998479748195),
        (10, 8, 0.000000000000000),
        (10, 9, 0.000000000000000),
        (10, 10, 1.000000000000000),
        (11, 1, 0.017995102963445),
        (11, 2, 0.019973601579310),
        (11, 3, 0.159559998572033),
        (11, 4, 0.339597749855903),
        (11, 5, 0.207078353555747),
        (11, 6, 0.161787341193781),
        (11, 7, -0.000000000000000),
        (11, 8, 0.423644139341113),
        (11, 9, 0.330987076821021),
        (11, 10, -0.000000000000000),
        (11, 11, 1.000000000000000),
        (12, 1, 0.055288573218445),
        (12, 2, 0.065316334525039),
        (12, 3, 0.327930464915055),
        (12, 4, 0.598191944670414),
        (12, 5, 0.160497791954347),
        (12, 6, 0.125394617939996),
        (12, 7, -0.000000000000000),
        (12, 8, 0.449149362937940),
        (12, 9, 0.350913941418061),
        (12, 10, -0.000000000000000),
        (12, 11, 0.683060190047392),
        (12, 12, 1.000000000000000),
        (13, 1, 0.017995102963445),
        (13, 2, 0.019973601579310),
        (13, 3, 0.159559998572033),
        (13, 4, 0.339597749855903),
        (13, 5, -0.207078353555747),
        (13, 6, 0.161787341193781),
        (13, 7, -0.000000000000000),
        (13, 8, -0.423644139341114),
        (13, 9, 0.330987076821021),
        (13, 10, -0.000000000000000),
        (13, 11, 0.020896239214577),
        (13, 12, 0.148173235303287),
        (13, 13, 1.000000000000000),
        (14, 1, 0.055288573218445),
        (14, 2, 0.065316334525039),
        (14, 3, 0.327930464915055),
        (14, 4, 0.598191944670414),
        (14, 5, -0.160497791954347),
        (14, 6, 0.125394617939996),
        (14, 7, -0.000000000000000),
        (14, 8, -0.449149362937940),
        (14, 9, 0.350913941418061),
        (14, 10, -0.000000000000000),
        (14, 11, 0.148173235303287),
        (14, 12, 0.385559366906972),
        (14, 13, 0.683060190047392),
        (14, 14, 1.000000000000000),
    ];
    let mut expctd_smn: DMatrix<f64> = DMatrix::zeros(14, 14);
    for (c, r, v) in expected_s_mn {
        expctd_smn[(c - 1, r - 1)] = v;
    }
    for (_i, (ex, comp)) in expctd_smn
        .lower_triangle()
        .iter()
        .zip(test_s_mn.iter())
        .enumerate()
    {
        assert_relative_eq!(ex, comp)
    }
}

#[test]
fn test_get_tmn_from_file() {
    let mut test_file = file_setup::setup_t_mn_water().unwrap();
    let test_t_mn = from_file::get_t_mn_from_fortran_format_file(&mut test_file, 14);
    assert_relative_eq!(52.614894976915927_f64, test_t_mn[(0, 0)]);
    assert_relative_eq!(0.266400000000000_f64, test_t_mn[(13, 13)]);
    assert_relative_eq!(0.085062970728081_f64, test_t_mn[(10, 5)]);
}

#[test]
fn test_get_vmn_from_file() {
    let mut test_file = file_setup::setup_v_nuc_el_water().unwrap();
    let test_t_mn = from_file::get_v_mn_from_fortran_format_file(&mut test_file, 14);
    assert_relative_eq!(-82.499042646076362, test_t_mn[(0, 0)]);
    assert_relative_eq!(-4.517555775212115, test_t_mn[(13, 13)]);
    assert_relative_eq!(-1.238061146392203, test_t_mn[(10, 5)]);
}

#[test]
fn test_get_eri_from_file() {
    let mut test_file = file_setup::setup_eri_water().unwrap();
    let test_eri = from_file::get_v_ee_from_fortran_format_file(&mut test_file, 14);
    assert_eq!(test_eri.len(), 5565);
    assert_eq!(test_eri[0], 6.354189973248753);
    assert_eq!(test_eri[18],0.472289079534761);
    assert_eq!(test_eri[5515], 0.007712593683891);
    assert_eq!(test_eri[1887], 0.025526245601437);
    assert_eq!(test_eri[5564], 0.475528488258028);
}
