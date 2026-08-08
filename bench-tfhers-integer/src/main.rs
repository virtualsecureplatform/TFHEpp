use std::hint::black_box;
use std::time::{Duration, Instant};

use tfhe::prelude::*;
use tfhe::shortint::parameters::PARAM_MESSAGE_2_CARRY_2_KS_PBS_GAUSSIAN_2M64;
use tfhe::{generate_keys, set_server_key, ConfigBuilder, FheUint128, FheUint64};

const RUNS: usize = 5;
const LHS: u64 = 0xfedc_ba98_7654_3210;
const RHS: u64 = 0x1234_5678_9abc_def0;

fn milliseconds(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1_000.0
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn sample_stddev(values: &[f64]) -> f64 {
    let average = mean(values);
    let squared_error = values
        .iter()
        .map(|value| (value - average).powi(2))
        .sum::<f64>();
    (squared_error / (values.len() - 1) as f64).sqrt()
}

fn main() {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .expect("the global Rayon pool must not already be initialized");

    let config =
        ConfigBuilder::with_custom_parameters(PARAM_MESSAGE_2_CARRY_2_KS_PBS_GAUSSIAN_2M64).build();
    let (client_key, server_key) = generate_keys(config);
    set_server_key(server_key);

    let lhs = FheUint64::encrypt(LHS, &client_key);
    let rhs = FheUint64::encrypt(RHS, &client_key);
    let expected = u128::from(LHS) * u128::from(RHS);

    // Widening is necessary because FheUint64 multiplication is reduced
    // modulo 2^64, whereas the TFHE2CLPX pipeline returns all 128 product bits.
    let warm_lhs = FheUint128::cast_from(lhs.clone());
    let warm_rhs = FheUint128::cast_from(rhs.clone());
    let warm_product = black_box(&warm_lhs * &warm_rhs);
    let warm_decrypted: u128 = warm_product.decrypt(&client_key);
    assert_eq!(warm_decrypted, expected);

    let mut widening_ms = Vec::with_capacity(RUNS);
    let mut multiplication_ms = Vec::with_capacity(RUNS);
    let mut total_ms = Vec::with_capacity(RUNS);

    for run in 0..RUNS {
        let widening_start = Instant::now();
        let wide_lhs = FheUint128::cast_from(lhs.clone());
        let wide_rhs = FheUint128::cast_from(rhs.clone());
        let widening = widening_start.elapsed();

        let multiplication_start = Instant::now();
        let product = black_box(&wide_lhs * &wide_rhs);
        let multiplication = multiplication_start.elapsed();

        let decrypted: u128 = product.decrypt(&client_key);
        assert_eq!(decrypted, expected, "incorrect product in run {run}");

        let widening = milliseconds(widening);
        let multiplication = milliseconds(multiplication);
        widening_ms.push(widening);
        multiplication_ms.push(multiplication);
        total_ms.push(widening + multiplication);
        println!(
            "run={} widening_ms={:.3} multiplication_ms={:.3} total_ms={:.3}",
            run + 1,
            widening,
            multiplication,
            widening + multiplication
        );
    }

    println!("tfhe_rs_version=1.6.1");
    println!("tfhe_rs_git=99ded5fc5222b0ed24f4fc5420e0a04c56b5c88a");
    println!("rayon_threads=1");
    println!("operand_bits=64 output_bits=128 runs={RUNS}");
    println!("parameter=PARAM_MESSAGE_2_CARRY_2_KS_PBS_GAUSSIAN_2M64");
    println!(
        "widening_mean_ms={:.3} widening_stddev_ms={:.3}",
        mean(&widening_ms),
        sample_stddev(&widening_ms)
    );
    println!(
        "multiplication_mean_ms={:.3} multiplication_stddev_ms={:.3}",
        mean(&multiplication_ms),
        sample_stddev(&multiplication_ms)
    );
    println!(
        "total_mean_ms={:.3} total_stddev_ms={:.3}",
        mean(&total_ms),
        sample_stddev(&total_ms)
    );
    println!("Passed");
}
