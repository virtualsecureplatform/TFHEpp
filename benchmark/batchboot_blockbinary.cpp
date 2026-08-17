// End-to-end operational-throughput benchmark for the two implementations.
//
// Build this target twice: once with the default parameter set (BatchBoot) and
// once with -DUSE_BLOCK_BINARY=ON (block-binary NAND).
// The configurations deliberately stay separate: the current block-binary
// lvl1 key has k=2, while generic BatchRingSwitch requires a k=1 source.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

#include "../include/tfhe++.hpp"
#include "google-benchmark/include/benchmark/benchmark.h"

#ifndef BATCHBOOT_BENCHMARK_REPETITIONS
#define BATCHBOOT_BENCHMARK_REPETITIONS 7
#endif

#ifdef USE_BLOCK_BINARY

static void BM_BlockBinaryNAND(benchmark::State& state)
{
    using bkP = TFHEpp::lvl01param;
    using iksP = TFHEpp::lvl10param;

    TFHEpp::SecretKey secret_key;
    TFHEpp::EvalKey eval_key;
    TFHEpp::TLWE<TFHEpp::lvl0param> left;
    TFHEpp::TLWE<TFHEpp::lvl0param> right;
    TFHEpp::TLWE<TFHEpp::lvl0param> result;
    bool initialized = false;

    for (auto _ : state) {
        if (!initialized) {
            // Google Benchmark starts its timer on entry to the loop, so pause
            // it here (rather than before the loop) for setup and warm-up.
            state.PauseTiming();
            eval_key.emplacebkfft<bkP>(secret_key);
            eval_key.emplacesubiksk<iksP>(secret_key);
            TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
                left, TFHEpp::lvl0param::μ, TFHEpp::lvl0param::α,
                secret_key.key.getSubset<TFHEpp::lvl0param>());
            TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
                right, -TFHEpp::lvl0param::μ, TFHEpp::lvl0param::α,
                secret_key.key.getSubset<TFHEpp::lvl0param>());
            TFHEpp::HomNAND<bkP, TFHEpp::lvl1param::μ, iksP>(result, left,
                                                             right, eval_key);
            initialized = true;
            state.ResumeTiming();
        }
        TFHEpp::HomNAND<bkP, TFHEpp::lvl1param::μ, iksP>(result, left, right,
                                                         eval_key);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations());
    state.SetLabel("one Boolean NAND / PBS");
}

BENCHMARK(BM_BlockBinaryNAND)
    ->Iterations(1)
    ->Repetitions(BATCHBOOT_BENCHMARK_REPETITIONS)
    ->DisplayAggregatesOnly(true);

#else

namespace {
// Parameter-Selection's 128-bit candidate source: public fixed-weight 64
// binary RLWE over the 64-bit level-2 ring.  The accumulator remains the
// default dense ternary level-2 key.  Keeping it a distinct type prevents the
// sparse source distribution from being accidentally used for the target.
struct SecureBatchSourceP : TFHEpp::lvl2param {
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_max = 1;
    static const inline double α = std::pow(2.0, -25);
};

using SourceP = SecureBatchSourceP;
using TargetP = TFHEpp::lvl2param;
#ifndef BATCHBOOT_BENCHMARK_SLOTS
#define BATCHBOOT_BENCHMARK_SLOTS 8
#endif
constexpr std::uint32_t slots = BATCHBOOT_BENCHMARK_SLOTS;
static_assert(std::has_single_bit(slots) && slots <= SourceP::n);
constexpr std::uint32_t source_weight = 64;
using ModuleP = TFHEpp::BatchRingSwitchP<SourceP, slots>;

void BM_BatchBoot(benchmark::State& state)
{
    constexpr std::uint32_t input_bits = 2;
    constexpr std::uint32_t stride = SourceP::n / slots;

    TFHEpp::Key<SourceP> sparse_key{};
    TFHEpp::Key<TargetP> accumulator_key{};
    TFHEpp::Polynomial<SourceP> plaintext{};
    TFHEpp::TRLWE<SourceP> packed_input;
    TFHEpp::Key<ModuleP> module_key{};
    TFHEpp::BatchBootKey<ModuleP, TargetP> boot_key;
    TFHEpp::AnnihilateKey<TargetP> automorphism_keys;
    const auto identity = TFHEpp::MakeBatchBootTestVector<TargetP, input_bits>(
        [](const std::uint32_t value) { return value; });
    TFHEpp::TRLWE<ModuleP> module_input;
    TFHEpp::TRLWE<TargetP> result;
    bool initialized = false;

    for (auto _ : state) {
        if (!initialized) {
            state.PauseTiming();
            // The support is part of the secret. Sample it uniformly without
            // replacement, matching the fixed-weight sparse-binary screen.
            std::array<std::uint32_t, SourceP::n> positions{};
            std::iota(positions.begin(), positions.end(), std::uint32_t{0});
            for (std::uint32_t selected = 0; selected < source_weight;
                 selected++) {
                std::uniform_int_distribution<std::uint32_t> supportgen(
                    selected, SourceP::n - 1);
                const std::uint32_t replacement = supportgen(TFHEpp::generator);
                std::swap(positions[selected], positions[replacement]);
                sparse_key[positions[selected]] = 1;
            }
            TFHEpp::keyGen<TargetP>(accumulator_key);
            for (std::uint32_t i = 0; i < slots; i++)
                plaintext[i * stride] =
                    static_cast<SourceP::T>(i & 1)
                    << (std::numeric_limits<SourceP::T>::digits - input_bits);
            TFHEpp::trlweSymEncrypt<SourceP>(packed_input, plaintext, 0.0,
                                             sparse_key);
            module_key =
                TFHEpp::BatchRingSwitchSecret<SourceP, slots>(sparse_key);
            TFHEpp::BatchBootKeyGen<ModuleP, TargetP>(boot_key, module_key,
                                                      accumulator_key);
            TFHEpp::annihilatekeygen<TargetP>(automorphism_keys,
                                              accumulator_key);
            TFHEpp::BatchRingSwitch<SourceP, slots>(module_input, packed_input);
            TFHEpp::BatchBootstrap<ModuleP, TargetP>(
                result, module_input, identity.polynomial, boot_key,
                automorphism_keys, identity.exponent_bias);
            const auto phase =
                TFHEpp::trlwePhase<TargetP>(result, accumulator_key);
            for (std::uint32_t i = 0; i < slots; i++)
                if (TFHEpp::BatchTorusDecode<input_bits>(phase[i * stride]) !=
                    (i & 1))
                    throw std::runtime_error(
                        "secure BatchBoot benchmark correctness check failed");
            initialized = true;
            state.ResumeTiming();
        }
        TFHEpp::BatchRingSwitch<SourceP, slots>(module_input, packed_input);
        TFHEpp::BatchBootstrap<ModuleP, TargetP>(
            result, module_input, identity.polynomial, boot_key,
            automorphism_keys, identity.exponent_bias);
        benchmark::DoNotOptimize(result);
    }
    state.SetItemsProcessed(state.iterations() * slots);
    state.SetLabel("128-bit screen: h=64 source, " + std::to_string(slots) +
                   " packed 1-bit outputs / job");
}
}  // namespace

BENCHMARK(BM_BatchBoot)
    ->Iterations(1)
    ->Repetitions(BATCHBOOT_BENCHMARK_REPETITIONS)
    ->DisplayAggregatesOnly(true);

#endif

BENCHMARK_MAIN();
