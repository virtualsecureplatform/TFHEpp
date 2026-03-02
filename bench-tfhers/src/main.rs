// pbs_tfhers — tfhe-rs boolean NAND gate benchmark for comparison with TFHEpp.
//
// Uses tfhe-rs DEFAULT_PARAMETERS which matches the updated TFHEpp tfhe-rs.hpp:
//   lwe_dimension  : 805   (n)
//   glwe_dimension : 3     (k)
//   polynomial_size: 512   (N)
//   pbs_base_log   : 10    (Bgbit)
//   pbs_level      : 2     (l)
//   ks_base_log    : 3     (basebit)
//   ks_level       : 5     (t)
//   lwe_std_dev    : 5.86e-6
//   glwe_std_dev   : 9.32e-10
//   security       : 132-bit, p-fail = 2^-64.344
//
// Outputs one line: "<ms-per-gate>ms"
// Followed by "Passed" if correctness checks pass.
//
// Matching C++ benchmark: test/pbs_tfhers.cpp

use std::time::Instant;
use tfhe::boolean::prelude::*;

fn main() {
    const NUM_TESTS: usize = 1000;

    // Key generation uses DEFAULT_PARAMETERS internally.
    // This is the slow step (bootstrapping key generation); not timed.
    let (client_key, server_key) = gen_keys();

    // Prepare random inputs
    let inputs_a: Vec<bool> = (0..NUM_TESTS).map(|i| i % 2 == 0).collect();
    let inputs_b: Vec<bool> = (0..NUM_TESTS).map(|i| (i / 2) % 2 == 0).collect();

    let ct_a: Vec<Ciphertext> = inputs_a.iter().map(|&b| client_key.encrypt(b)).collect();
    let ct_b: Vec<Ciphertext> = inputs_b.iter().map(|&b| client_key.encrypt(b)).collect();

    // Warmup: one gate to ensure JIT / cache effects are settled
    let _ = server_key.nand(&ct_a[0], &ct_b[0]);

    // Timed benchmark
    let start = Instant::now();
    let results: Vec<Ciphertext> = ct_a
        .iter()
        .zip(ct_b.iter())
        .map(|(a, b)| server_key.nand(a, b))
        .collect();
    let elapsed = start.elapsed();

    let ms_per_gate = elapsed.as_secs_f64() * 1000.0 / NUM_TESTS as f64;
    println!("{:.3}ms", ms_per_gate);

    // Correctness check
    for (i, ct_res) in results.iter().enumerate() {
        let expected = !(inputs_a[i] && inputs_b[i]);
        let got = client_key.decrypt(ct_res);
        assert_eq!(
            got, expected,
            "NAND({}, {}) = {} but expected {}",
            inputs_a[i], inputs_b[i], got, expected
        );
    }
    println!("Passed");
}
