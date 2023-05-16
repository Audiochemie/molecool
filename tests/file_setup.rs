use std::{fs::File, io::{BufReader,Result}}; 
/// Wrapper function to setup required data for integration tests of data module
///
pub fn setup_allene_symbolic() -> Result<BufReader<File>> {
    let test_file = File::open("tests/test_allene_symbolic.xyz")?;
    Ok(BufReader::new(test_file))
}
pub fn setup_acetaldehyde_numeric() -> Result<BufReader<File>> {
    let test_file = File::open("tests/test_acetaldehyde_numeric.xyz")?;
    Ok(BufReader::new(test_file))
}

