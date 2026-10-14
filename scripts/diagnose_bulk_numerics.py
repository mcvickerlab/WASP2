"""Compare fused/unfused Rust Cephes arithmetic with the installed SciPy.

Builds two isolated probe modules from the production source in a temporary
crate. The repository's source and numerical policy are never modified.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import platform
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import scipy
from scipy.special import gamma, gammaln
from scipy.stats import betabinom

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    source = (ROOT / "rust/src/scipy_compat.rs").read_text()
    start = source.index("#[inline]\nfn cephes_mul_add")
    end = source.index("\n}\n", start) + 3
    implementation = source[start:end]
    signature = "#[inline]\nfn cephes_mul_add(a: f64, b: f64, c: f64) -> f64"
    fused = source.replace(implementation, signature + " { a.mul_add(b, c) }\n")
    plain = source.replace(implementation, signature + " { a * b + c }\n")
    wrapper = "\npub(crate) fn gamma_probe(x: f64) -> (f64, f64) { (gamma_positive(x), lgamma_positive(x)) }\n"
    rust_main = r"""
mod fused;
mod plain;
fn main() {
    for x in [0.1, 9.548148325980558, 16.97709256606237, 1000.0, 1e8, 5e11] {
        let f = fused::gamma_probe(x);
        let p = plain::gamma_probe(x);
        println!("gamma,{x:?},{:?},{:?},{:?},{:?}", f.0, p.0, f.1, p.1);
    }
    for rho in [1e-12, 1e-9, 0.01, 0.1, 0.99] {
        for k in [0, 10, 20] {
            let a = 0.5 * (1.0-rho)/rho;
            let f = fused::betabinom_logpmf(k, 20, a, a).unwrap();
            let p = plain::betabinom_logpmf(k, 20, a, a).unwrap();
            println!("pmf,{k},20,{a:?},{a:?},{f:?},{p:?}");
        }
    }
    let (a,b) = (8.977092566062371,5.5481483259805575);
    println!("pmf,8,12,{a:?},{b:?},{:?},{:?}",
        fused::betabinom_logpmf(8,12,a,b).unwrap(),
        plain::betabinom_logpmf(8,12,a,b).unwrap());
    // Keep both reduction functions live when compiling with warning checks.
    assert_eq!(fused::numpy_sum(&[1.0,2.0]),plain::numpy_sum(&[1.0,2.0]));
}
"""
    with tempfile.TemporaryDirectory(prefix="wasp2-numerical-probe-") as directory:
        temp = Path(directory)
        (temp / "src").mkdir()
        (temp / "Cargo.toml").write_text(
            '[package]\nname="wasp2-numerical-probe"\nversion="0.0.0"\nedition="2021"\n'
            '[dependencies]\nanyhow="1"\n'
        )
        (temp / "src/fused.rs").write_text(fused + wrapper)
        (temp / "src/plain.rs").write_text(plain + wrapper)
        (temp / "src/main.rs").write_text(rust_main)
        env = dict(os.environ)
        env["CARGO_TARGET_DIR"] = str(ROOT / "rust/target/numerical-probe")
        completed = subprocess.run(
            ["cargo", "run", "--release", "--manifest-path", str(temp / "Cargo.toml")],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
    rows = []
    logs = []
    for row in csv.reader(io.StringIO(completed.stdout)):
        if row[0] == "gamma":
            x, fg, pg, fl, pl = map(float, row[1:])
            rows.append(
                {
                    "kind": "gamma",
                    "x": x,
                    "python_gamma": float(gamma(x)),
                    "fused_gamma": fg,
                    "plain_gamma": pg,
                    "python_lgamma": float(gammaln(x)),
                    "fused_lgamma": fl,
                    "plain_lgamma": pl,
                }
            )
        else:
            k, n, a, b, fused, plain_value = map(float, row[1:])
            reference = float(betabinom.logpmf(k, n, a, b))
            logs.append(reference)
            rows.append(
                {
                    "kind": "logpmf",
                    "k": k,
                    "n": n,
                    "alpha": a,
                    "beta": b,
                    "python": reference,
                    "fused": fused,
                    "plain": plain_value,
                }
            )
    values = np.array(logs * 32)
    scalar_exp = np.array([math.exp(x) for x in values])
    vector_exp = np.exp(values)
    numpy_so = Path(np.core._multiarray_umath.__file__)
    report = {
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "numpy_binary_sha256": hashlib.sha256(numpy_so.read_bytes()).hexdigest(),
        "numpy_cpu_features": np.core._multiarray_umath.__cpu_features__,
        "scalar_vector_exp_differences": int(np.count_nonzero(scalar_exp != vector_exp)),
        "scalar_vector_exp_max_delta": float(np.max(np.abs(scalar_exp - vector_exp))),
        "rows": rows,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
