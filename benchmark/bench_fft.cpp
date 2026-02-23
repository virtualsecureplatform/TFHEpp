// Pin benchmark to a specific CPU core for stable measurements
#include <sched.h>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>

#include "../include/tfhe++.hpp"
#include "google-benchmark/include/benchmark/benchmark.h"

// Pin the process to CPU core 0 (3D V-Cache CCD on Ryzen 9950X3D)
static void pin_to_core(int core_id)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core_id, &cpuset);
    if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuset) != 0)
        fprintf(stderr, "Warning: sched_setaffinity failed\n");
}

// Global initialization: pin to core 0 before any benchmarks run
struct PinInit {
    PinInit() { pin_to_core(0); }
} pin_init;

// ============================================================
// Benchmark: IFFT (execute_reverse_torus32) for lvl1
// ============================================================
static void BM_IFFT_lvl1(benchmark::State& state)
{
    constexpr int N = TFHEpp::lvl1param::n;
    alignas(64) std::array<uint32_t, N> input;
    alignas(64) TFHEpp::PolynomialInFD<TFHEpp::lvl1param> output;

    std::mt19937 rng(42);
    for (auto& v : input) v = rng();

    for (auto _ : state) {
        fftplvl1.execute_reverse_torus32(output.data(), input.data());
        benchmark::DoNotOptimize(output.data());
    }
}

// ============================================================
// Benchmark: FFT (execute_direct_torus32) for lvl1
// ============================================================
static void BM_FFT_lvl1(benchmark::State& state)
{
    constexpr int N = TFHEpp::lvl1param::n;
    alignas(64) TFHEpp::PolynomialInFD<TFHEpp::lvl1param> input;
    alignas(64) std::array<uint32_t, N> output;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1e9, 1e9);
    for (auto& v : input) v = dist(rng);

    for (auto _ : state) {
        fftplvl1.execute_direct_torus32(output.data(), input.data());
        benchmark::DoNotOptimize(output.data());
    }
}

// ============================================================
// Benchmark: MulInFD (pointwise complex multiplication)
// ============================================================
static void BM_MulInFD_lvl1(benchmark::State& state)
{
    constexpr int N = TFHEpp::lvl1param::n;
    alignas(64) TFHEpp::PolynomialInFD<TFHEpp::lvl1param> a, b, res;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& v : a) v = dist(rng);
    for (auto& v : b) v = dist(rng);

    for (auto _ : state) {
        TFHEpp::MulInFD<N>(res, a, b);
        benchmark::DoNotOptimize(res.data());
    }
}

// ============================================================
// Benchmark: FMAInFD (pointwise complex FMA - hot path)
// ============================================================
static void BM_FMAInFD_lvl1(benchmark::State& state)
{
    constexpr int N = TFHEpp::lvl1param::n;
    alignas(64) TFHEpp::PolynomialInFD<TFHEpp::lvl1param> a, b, res;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& v : a) v = dist(rng);
    for (auto& v : b) v = dist(rng);

    for (auto _ : state) {
        // Reset accumulator each iteration to avoid FP overflow
        memset(res.data(), 0, sizeof(res));
        TFHEpp::FMAInFD<N>(res, a, b);
        benchmark::DoNotOptimize(res.data());
    }
}

// ============================================================
// Benchmark: Full PolyMul pipeline for lvl1
// ============================================================
static void BM_PolyMul_lvl1(benchmark::State& state)
{
    constexpr int N = TFHEpp::lvl1param::n;
    alignas(64) TFHEpp::Polynomial<TFHEpp::lvl1param> a, b, res;

    std::mt19937 rng(42);
    for (auto& v : a) v = rng();
    for (auto& v : b) v = rng();

    for (auto _ : state) {
        TFHEpp::PolyMul<TFHEpp::lvl1param>(res, a, b);
        benchmark::DoNotOptimize(res.data());
    }
}

// ============================================================
// Benchmark: TwistIFFT + MulInFD + TwistFFT (simulates ExternalProduct core)
// ============================================================
static void BM_IFFTMulFFT_lvl1(benchmark::State& state)
{
    constexpr int N = TFHEpp::lvl1param::n;
    alignas(64) TFHEpp::Polynomial<TFHEpp::lvl1param> a, res;
    alignas(64) TFHEpp::PolynomialInFD<TFHEpp::lvl1param> fftb;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& v : a) v = rng();
    for (auto& v : fftb) v = dist(rng);

    for (auto _ : state) {
        alignas(64) TFHEpp::PolynomialInFD<TFHEpp::lvl1param> ffta;
        TFHEpp::TwistIFFT<TFHEpp::lvl1param>(ffta, a);
        TFHEpp::MulInFD<N>(ffta, fftb);
        TFHEpp::TwistFFT<TFHEpp::lvl1param>(res, ffta);
        benchmark::DoNotOptimize(res.data());
    }
}

BENCHMARK(BM_IFFT_lvl1)
    ->Iterations(10000)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_FFT_lvl1)
    ->Iterations(10000)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_MulInFD_lvl1)
    ->Iterations(10000)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_FMAInFD_lvl1)
    ->Iterations(10000)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_PolyMul_lvl1)
    ->Iterations(10000)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_IFFTMulFFT_lvl1)
    ->Iterations(10000)
    ->Repetitions(5)
    ->DisplayAggregatesOnly(true);
BENCHMARK_MAIN();
