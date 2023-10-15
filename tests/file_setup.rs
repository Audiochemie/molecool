//! Wrapper functions to setup required data for integration tests of qc_file_parsers
use std::{
    fs::File,
    io::{BufReader, Result},
};

pub fn setup_allene_symbolic() -> Result<BufReader<File>> {
    let test_file = File::open("tests/test_allene_symbolic.xyz")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_acetaldehyde_numeric() -> Result<BufReader<File>> {
    let test_file = File::open("tests/test_acetaldehyde_numeric.xyz")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_acetaldehyde_symbolic() -> Result<BufReader<File>> {
    let test_file = File::open("tests/test_acetaldehyde_numeric.xyz")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_water_geo() -> Result<BufReader<File>> {
    let test_file = File::open("tests/h2o_geom.xyz")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_water_hessian() -> Result<BufReader<File>> {
    let test_file = File::open("tests/h2o_hessian.txt")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_chlorbenzene_geo() -> Result<BufReader<File>> {
    let test_file = File::open("tests/cl_b_geom.xyz")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_chlorbenzene_hessian() -> Result<BufReader<File>> {
    let test_file = File::open("tests/cl_b_hessian.txt")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_v_nuc_nuc_water() -> Result<BufReader<File>> {
    let test_file = File::open("tests/water_integrals/h2o_vnn.dat")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_v_nuc_el_water() -> Result<BufReader<File>> {
    let test_file = File::open("tests/water_integrals/h2o_ve.dat")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_s_mn_water() -> Result<BufReader<File>> {
    let test_file = File::open("tests/water_integrals/h2o_smn.dat")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_t_mn_water() -> Result<BufReader<File>> {
    let test_file = File::open("tests/water_integrals/h2o_t.dat")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_s_mn_water_faulty() -> Result<BufReader<File>> {
    let test_file = File::open("tests/water_integrals/h2o_smn_faulty.dat")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_s_mn_water_faulty2() -> Result<BufReader<File>> {
    let test_file = File::open("tests/water_integrals/h2o_smn_faulty2.dat")?;
    Ok(BufReader::new(test_file))
}

pub fn setup_eri_water() -> Result<BufReader<File>> {
    let test_file = File::open("tests/water_integrals/h2o_eri.dat")?;
    Ok(BufReader::new(test_file))
}
