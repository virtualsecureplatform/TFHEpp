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
// Outputs BR, IKS, and full NAND timings, followed by "Passed" if correct.
//
// Matching C++ benchmark: test/pbs_tfhers.cpp

use std::time::Instant;
use tfhe::boolean::prelude::*;
use tfhe::core_crypto::prelude::*;

fn main() {
    const NUM_TESTS: usize = 1000;
    // 1/8 in torus32 encoding (same as tfhe-rs internal PLAINTEXT_TRUE)
    const PLAINTEXT_TRUE: u32 = 1 << 29;

    // Force single-threaded execution so Rayon's parallel key switch does not
    // give tfhe-rs an unfair advantage over TFHEpp's sequential implementation.
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .unwrap();

    // Key generation uses DEFAULT_PARAMETERS internally.
    // This is the slow step (bootstrapping key generation); not timed.
    let (client_key, server_key) = gen_keys();

    // Prepare random inputs
    let inputs_a: Vec<bool> = (0..NUM_TESTS).map(|i| i % 2 == 0).collect();
    let inputs_b: Vec<bool> = (0..NUM_TESTS).map(|i| (i / 2) % 2 == 0).collect();

    let ct_a: Vec<Ciphertext> = inputs_a.iter().map(|&b| client_key.encrypt(b)).collect();
    let ct_b: Vec<Ciphertext> = inputs_b.iter().map(|&b| client_key.encrypt(b)).collect();

    // --- Full NAND benchmark ---
    let _ = server_key.nand(&ct_a[0], &ct_b[0]); // warmup
    let nand_start = Instant::now();
    let nand_results: Vec<Ciphertext> = ct_a
        .iter()
        .zip(ct_b.iter())
        .map(|(a, b)| server_key.nand(a, b))
        .collect();
    let nand_elapsed = nand_start.elapsed();

    // Extract raw LWE ciphertexts from Ciphertext::Encrypted
    let raw_a: Vec<&LweCiphertextOwned<u32>> = ct_a
        .iter()
        .map(|ct| match ct {
            Ciphertext::Encrypted(lwe) => lwe,
            _ => panic!("expected encrypted ciphertext"),
        })
        .collect();
    let raw_b: Vec<&LweCiphertextOwned<u32>> = ct_b
        .iter()
        .map(|ct| match ct {
            Ciphertext::Encrypted(lwe) => lwe,
            _ => panic!("expected encrypted ciphertext"),
        })
        .collect();

    // Decompose server_key into raw BSK and KSK for split benchmarks.
    // DEFAULT_PARAMETERS uses EncryptionKeyChoice::Small → BootstrapKeyswitch
    // order: PBS first (small→big LWE), then KS (big→small LWE).
    let (fourier_bsk, ksk, _pbs_order) = server_key.into_raw_parts();

    // Build accumulator (trivial GLWE: mask=0, body=PLAINTEXT_TRUE).
    // This replicates the boolean engine's LUT construction.
    let accumulator = {
        let mut acc = GlweCiphertextOwned::new(
            0u32,
            fourier_bsk.glwe_size(),
            fourier_bsk.polynomial_size(),
            CiphertextModulus::new_native(),
        );
        acc.get_mut_body().as_mut().fill(PLAINTEXT_TRUE);
        acc
    };

    // Prepare NAND inputs: -(ct_a + ct_b) + 1/8
    let small_lwe_size = raw_a[0].lwe_size();
    let prepared: Vec<LweCiphertextOwned<u32>> = (0..NUM_TESTS)
        .map(|i| {
            let mut buf =
                LweCiphertext::new(0u32, small_lwe_size, CiphertextModulus::new_native());
            lwe_ciphertext_add(&mut buf, raw_a[i], raw_b[i]);
            lwe_ciphertext_opposite_assign(&mut buf);
            lwe_ciphertext_plaintext_add_assign(&mut buf, Plaintext(PLAINTEXT_TRUE));
            buf
        })
        .collect();

    // Warmup for split benchmarks
    {
        let big_lwe_size = fourier_bsk.output_lwe_dimension().to_lwe_size();
        let mut tmp_big =
            LweCiphertext::new(0u32, big_lwe_size, CiphertextModulus::new_native());
        programmable_bootstrap_lwe_ciphertext(
            &prepared[0],
            &mut tmp_big,
            &accumulator,
            &fourier_bsk,
        );
        let mut tmp_small =
            LweCiphertext::new(0u32, small_lwe_size, CiphertextModulus::new_native());
        keyswitch_lwe_ciphertext(&ksk, &tmp_big, &mut tmp_small);
    }

    // --- BR benchmark ---
    let big_lwe_size = fourier_bsk.output_lwe_dimension().to_lwe_size();
    let mut big_lwe: Vec<LweCiphertextOwned<u32>> = (0..NUM_TESTS)
        .map(|_| LweCiphertext::new(0u32, big_lwe_size, CiphertextModulus::new_native()))
        .collect();

    let br_start = Instant::now();
    for i in 0..NUM_TESTS {
        programmable_bootstrap_lwe_ciphertext(
            &prepared[i],
            &mut big_lwe[i],
            &accumulator,
            &fourier_bsk,
        );
    }
    let br_elapsed = br_start.elapsed();

    // --- IKS benchmark ---
    // KS output has the same dimension as the original small LWE
    let mut small_lwe: Vec<LweCiphertextOwned<u32>> = (0..NUM_TESTS)
        .map(|_| LweCiphertext::new(0u32, small_lwe_size, CiphertextModulus::new_native()))
        .collect();

    let iks_start = Instant::now();
    for i in 0..NUM_TESTS {
        keyswitch_lwe_ciphertext(&ksk, &big_lwe[i], &mut small_lwe[i]);
    }
    let iks_elapsed = iks_start.elapsed();

    // Output
    let br_ms = br_elapsed.as_secs_f64() * 1000.0 / NUM_TESTS as f64;
    let iks_ms = iks_elapsed.as_secs_f64() * 1000.0 / NUM_TESTS as f64;
    let nand_ms = nand_elapsed.as_secs_f64() * 1000.0 / NUM_TESTS as f64;
    println!("BR: {:.3}ms", br_ms);
    println!("IKS: {:.3}ms", iks_ms);
    println!("NAND: {:.3}ms", nand_ms);

    // Correctness check (split BR+IKS)
    for i in 0..NUM_TESTS {
        let expected = !(inputs_a[i] && inputs_b[i]);
        let ct_res = Ciphertext::Encrypted(small_lwe[i].clone());
        let got: bool = client_key.decrypt(&ct_res);
        assert_eq!(
            got, expected,
            "BR+IKS NAND({}, {}) = {} but expected {}",
            inputs_a[i], inputs_b[i], got, expected
        );
    }
    // Correctness check (full NAND)
    for (i, ct_res) in nand_results.iter().enumerate() {
        let expected = !(inputs_a[i] && inputs_b[i]);
        let got: bool = client_key.decrypt(ct_res);
        assert_eq!(
            got, expected,
            "NAND({}, {}) = {} but expected {}",
            inputs_a[i], inputs_b[i], got, expected
        );
    }
    println!("Passed");
}
