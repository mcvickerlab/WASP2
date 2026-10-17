fn main() {
    pyo3_build_config::add_extension_module_link_args();
    pyo3_build_config::use_pyo3_cfgs();
    if matches!(std::env::var("CARGO_CFG_TARGET_OS").as_deref(), Ok("macos")) {
        // Reserve space for library-path relocation when wheels are repackaged by Conda.
        println!("cargo:rustc-link-arg-cdylib=-Wl,-headerpad_max_install_names");
    }
}
