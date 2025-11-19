use pyo3::prelude::*;

mod bam_counter;
use bam_counter::BamCounter;

/// Simple test function to verify PyO3 is working
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// WASP2 Rust acceleration module
#[pymodule]
fn wasp2_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<BamCounter>()?;
    Ok(())
}
