// fft-tfhers-bench — tfhe-rs FFT benchmark for comparison with TFHEpp bench_fft.
//
// Benchmarks individual FFT operations at N=512 (matching tfhe-rs DEFAULT_PARAMETERS):
//   - Forward FFT  (coefficient -> Fourier domain)
//   - Backward FFT (Fourier -> coefficient domain)
//   - Pointwise complex multiply in Fourier domain
//   - Full polynomial multiplication (forward + mul + backward)
//   - IFFT + Mul + FFT pipeline (ExternalProduct core)
//
// Output format matches Google Benchmark style for easy comparison.

use std::hint::black_box;
use std::time::Instant;

use dyn_stack::{PodBuffer, PodStack};
use tfhe::core_crypto::fft_impl::fft64::math::fft::Fft;
use tfhe::core_crypto::prelude::polynomial_algorithms::polynomial_wrapping_mul;
use tfhe::core_crypto::prelude::*;
use tfhe_fft::c64;

const N: usize = 512;
const WARMUP: usize = 1000;
const ITERS: usize = 10000;
const REPS: usize = 5;

fn pin_to_core(core_id: usize) {
    #[cfg(target_os = "linux")]
    {
        unsafe {
            let mut cpuset = std::mem::zeroed::<libc::cpu_set_t>();
            libc::CPU_SET(core_id, &mut cpuset);
            libc::sched_setaffinity(
                0,
                std::mem::size_of::<libc::cpu_set_t>(),
                &cpuset,
            );
        }
    }
}

/// Run a benchmark: warmup, then REPS repetitions of ITERS iterations.
fn bench<F: FnMut()>(name: &str, mut f: F) {
    for _ in 0..WARMUP {
        f();
    }

    let mut times_ns = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let start = Instant::now();
        for _ in 0..ITERS {
            f();
        }
        let elapsed = start.elapsed();
        times_ns.push(elapsed.as_nanos() as f64 / ITERS as f64);
    }

    let mean = times_ns.iter().sum::<f64>() / REPS as f64;
    let variance =
        times_ns.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / REPS as f64;
    let stddev = variance.sqrt();
    let cv = if mean > 0.0 {
        stddev / mean * 100.0
    } else {
        0.0
    };

    println!(
        "{:<30} {:>8.1} ns  +/- {:>5.1} ns  (cv={:.1}%)",
        name, mean, stddev, cv
    );
}

fn main() {
    // Force single-threaded
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .unwrap();
    pin_to_core(0);

    println!("=== tfhe-rs FFT benchmark (N={}) ===", N);
    println!(
        "  {} warmup, {} iterations x {} repetitions",
        WARMUP, ITERS, REPS
    );
    println!();

    let poly_size = PolynomialSize(N);
    let fft = Fft::new(poly_size);

    // Allocate scratch buffer (enough for both forward and backward)
    let fwd_req = fft.as_view().forward_scratch();
    let bwd_req = fft.as_view().backward_scratch();
    let max_req = fwd_req.or(bwd_req);
    let mut scratch_buf = PodBuffer::new(max_req);

    // Prepare input data
    let mut rng = 42u64;
    let mut next_u64 = || -> u64 {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };

    let mut poly_data: Vec<u64> = (0..N).map(|_| next_u64()).collect();
    let poly = Polynomial::from_container(poly_data.as_mut_slice());

    let mut fourier = FourierPolynomial::new(poly_size);

    // Pre-compute a fourier polynomial for mul benchmarks
    {
        let stack = PodStack::new(&mut scratch_buf);
        fft.as_view()
            .forward_as_torus(fourier.as_mut_view(), poly.as_view(), stack);
    }

    let mut fourier_a = FourierPolynomial::new(poly_size);
    let mut fourier_b = FourierPolynomial::new(poly_size);
    fourier_a.data.clone_from(&fourier.data);
    fourier_b.data.clone_from(&fourier.data);

    // --- Benchmark: Forward FFT (coefficient -> Fourier domain) ---
    bench("Forward FFT (IFFT equiv)", || {
        let stack = PodStack::new(&mut scratch_buf);
        fft.as_view()
            .forward_as_torus(fourier.as_mut_view(), poly.as_view(), stack);
        black_box(&fourier);
    });

    // --- Benchmark: Backward FFT (Fourier -> coefficient domain) ---
    let mut out_data: Vec<u64> = vec![0u64; N];
    let mut out_poly = Polynomial::from_container(out_data.as_mut_slice());

    bench("Backward FFT (FFT equiv)", || {
        let stack = PodStack::new(&mut scratch_buf);
        fft.as_view()
            .backward_as_torus(out_poly.as_mut_view(), fourier_a.as_view(), stack);
        black_box(&out_poly);
    });

    // --- Benchmark: Pointwise complex multiply in Fourier domain ---
    let mut fourier_res = FourierPolynomial::new(poly_size);
    bench("MulInFD equiv", || {
        let r = &mut fourier_res.data;
        let a = &fourier_a.data;
        let b = &fourier_b.data;
        for i in 0..r.len() {
            let (ar, ai) = (a[i].re, a[i].im);
            let (br, bi) = (b[i].re, b[i].im);
            r[i] = c64 {
                re: ar * br - ai * bi,
                im: ar * bi + ai * br,
            };
        }
        black_box(&*r);
    });

    // --- Benchmark: FMA in Fourier domain ---
    let mut fourier_acc = FourierPolynomial::new(poly_size);
    bench("FMAInFD equiv", || {
        for v in fourier_acc.data.iter_mut() {
            *v = c64 { re: 0.0, im: 0.0 };
        }
        let acc = &mut fourier_acc.data;
        let a = &fourier_a.data;
        let b = &fourier_b.data;
        for i in 0..acc.len() {
            let (ar, ai) = (a[i].re, a[i].im);
            let (br, bi) = (b[i].re, b[i].im);
            acc[i].re += ar * br - ai * bi;
            acc[i].im += ar * bi + ai * br;
        }
        black_box(&*acc);
    });

    // --- Benchmark: Full polynomial multiplication ---
    let mut poly_a_data: Vec<u64> = (0..N).map(|_| next_u64()).collect();
    let mut poly_b_data: Vec<u64> = (0..N).map(|_| next_u64()).collect();
    let mut poly_res_data: Vec<u64> = vec![0u64; N];
    let poly_a = Polynomial::from_container(poly_a_data.as_mut_slice());
    let poly_b = Polynomial::from_container(poly_b_data.as_mut_slice());
    let mut poly_res = Polynomial::from_container(poly_res_data.as_mut_slice());

    bench("PolyMul (wrapping_mul)", || {
        polynomial_wrapping_mul(&mut poly_res, &poly_a, &poly_b);
        black_box(&poly_res);
    });

    // --- Benchmark: IFFT + Mul + FFT pipeline (ExternalProduct core) ---
    let mut pipeline_poly_data: Vec<u64> = (0..N).map(|_| next_u64()).collect();
    let pipeline_poly =
        Polynomial::from_container(pipeline_poly_data.as_mut_slice());
    let mut pipeline_fourier = FourierPolynomial::new(poly_size);
    let mut pipeline_out_data: Vec<u64> = vec![0u64; N];
    let mut pipeline_out =
        Polynomial::from_container(pipeline_out_data.as_mut_slice());

    bench("IFFT+Mul+FFT pipeline", || {
        // Forward: coefficient -> Fourier
        let stack = PodStack::new(&mut scratch_buf);
        fft.as_view().forward_as_torus(
            pipeline_fourier.as_mut_view(),
            pipeline_poly.as_view(),
            stack,
        );
        // Pointwise multiply
        let f = &mut pipeline_fourier.data;
        let b = &fourier_b.data;
        for i in 0..f.len() {
            let (ar, ai) = (f[i].re, f[i].im);
            let (br, bi) = (b[i].re, b[i].im);
            f[i] = c64 {
                re: ar * br - ai * bi,
                im: ar * bi + ai * br,
            };
        }
        // Backward: Fourier -> coefficient
        let stack = PodStack::new(&mut scratch_buf);
        fft.as_view().backward_as_torus(
            pipeline_out.as_mut_view(),
            pipeline_fourier.as_view(),
            stack,
        );
        black_box(&pipeline_out);
    });

    println!();
    println!("Done.");
}
