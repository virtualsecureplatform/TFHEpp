#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <type_traits>

// Concrete parameter records from Jain--Lin--Liu--Saha, ePrint 2026/1730.
//
// These are intentionally *not* aliases for TFHEpp's active lvl0/lvl1
// parameters.  The paper uses explicit Z_q/Z_Q moduli, while TFHEpp's current
// gate-bootstrap path uses a power-of-two torus and derives its blind-rotation
// modulus from the target ring dimension.  Activating the records before that
// representation gap is implemented would silently select different crypto
// parameters from the ones analysed in the paper.
namespace shallowboot {

enum class SecretDistribution {
    GeneralSparseBinary,
    StructuredSparseBinary,
    GeneralSparseTernary,
    StructuredSparseTernary,
    BinaryNTT,
};

struct Parameters {
    SecretDistribution secret_distribution;
    std::uint32_t lwe_dimension;
    std::uint32_t lwe_modulus_log2;
    std::uint32_t lwe_hamming_weight;
    std::uint32_t rlwe_dimension;
    std::uint32_t rlwe_modulus_log2;
    std::uint32_t intermediate_modulus_log2;
};

// Table 3: general sparse binary keys.  These require the PBC compiler from
// Algorithm 2; a plain one-hot block schedule is not sufficient.
inline constexpr Parameters general_binary_gate_std128 = {
    SecretDistribution::GeneralSparseBinary, 1024, 9, 43, 1024, 28, 12};
inline constexpr Parameters general_binary_gate_medium = {
    SecretDistribution::GeneralSparseBinary, 1024, 9, 31, 1024, 28, 12};
inline constexpr Parameters general_binary_function_std128 = {
    SecretDistribution::GeneralSparseBinary, 2048, 12, 70, 2048, 54, 35};

// Table 4: structured binary keys.  Their one-hot secret distribution matches
// TFHEpp's current StructuredSparseBlindRotate frontend once a compatible
// explicit Z_q/Z_Q parameter family is available.
inline constexpr Parameters structured_binary_gate_std128 = {
    SecretDistribution::StructuredSparseBinary, 1024, 9, 64, 1024, 28, 12};
inline constexpr Parameters structured_binary_gate_medium = {
    SecretDistribution::StructuredSparseBinary, 990, 9, 45, 1024, 28, 12};
inline constexpr Parameters structured_binary_function_std128 = {
    SecretDistribution::StructuredSparseBinary, 2040, 12, 102, 2048, 54,
    35};

// Parameter-Selection candidate for a post-key-switch blind-rotation key.
// Its source proxy screens at 130.9 bits and the existing accumulator proxy
// at 129.6 bits; see estimates/shallowboot_structured_std128_results.md.
// It is not executable until a q-aware dimension key switch is implemented.
inline constexpr Parameters structured_binary_short_gate_std128 = {
    SecretDistribution::StructuredSparseBinary, 576, 9, 64, 1024, 32, 12};

// Section 7.2 / Algorithm 3.  The larger modulus is represented by two
// stages: 105 bits before the first modulus switch, then approximately 50
// bits.  It requires the new Binary-NTT RLWE assumption and is not yet an
// executable TFHEpp parameter family.
struct LowDepthParameters : Parameters {
    double lwe_noise_stddev;
    double rlwe_noise_stddev;
    std::uint32_t multiplication_window;
    std::uint32_t key_switch_digits;
};

// Paper Section 7.2 row.  The local full sparse-secret screen currently
// reports 112.9 bits for its LWE component; retain it as a reproduction target
// rather than enabling it as a validated TFHEpp parameter set.
inline constexpr LowDepthParameters binary_ntt_paper_std128 = {
    {SecretDistribution::BinaryNTT, 1450, 9, 29, 4096, 105, 12},
    3.2,
    0.75,
    3,
    4,
};

// Local source-security candidate: changing the sparse source to h=37 gives
// a 133.3-bit uniform-sparse proxy under the vendored estimator.  PBC uses
// c=3 and k=h+3=40 leaves; the 16-by-4 schedule keeps its extra layer above
// the modulus boundary and retains two post-boundary product layers.
inline constexpr LowDepthParameters binary_ntt_source_screened = {
    {SecretDistribution::BinaryNTT, 1450, 9, 37, 4096, 105, 12},
    3.2,
    0.75,
    3,
    4,
};

// Concrete Section 7.2 execution schedule: three Binary-NTT multiplication
// layers at the initial ~105-bit modulus, then two post-switch layers near 50
// bits.  The RNS implementation uses two 62-bit primes to represent the
// initial level; exact modulus-down/key-switch arithmetic remains experimental.
struct LowDepthExecutionSchedule {
    std::uint32_t initial_modulus_bits;
    std::uint32_t post_switch_modulus_bits;
    std::array<std::uint32_t, 2> multiplication_windows;
    std::uint32_t initial_rns_primes;
};

inline constexpr LowDepthExecutionSchedule binary_ntt_std128_schedule = {
    105, 50, {8, 4}, 2};

inline constexpr LowDepthExecutionSchedule
    binary_ntt_source_screened_schedule = {105, 50, {16, 4}, 2};

// Security- and noise-screened QH-SS candidate for the direct-RLWE PBC tree.
// Input, refreshed-output, and ring LWE proxies exceed 133, 138, and 188
// classical bits. The R_4 selector tree uses modular-inverse LSB-to-MSB
// conversion before its dense sign LUT. The refreshed key uses h=50 at the
// 15-bit ring-to-LWE key-switch boundary; the paper's 12-bit h=37 choice does
// not retain enough key-switch noise margin in this implementation.
struct LowDepthSecureRNSParameters {
    std::uint32_t ring_dimension;
    std::array<std::uint32_t, 4> modulus_bits;
    std::array<std::uint32_t, 4> multiplication_windows;
    std::uint32_t pbc_copies;
    std::uint32_t pbc_slack;
    std::uint32_t plaintext_modulus;
    std::uint32_t key_switch_modulus_bits;
    std::uint32_t refreshed_lwe_dimension;
    std::uint32_t refreshed_lwe_hamming_weight;
};

inline constexpr LowDepthSecureRNSParameters qh_ss_source_screened_rns = {
    8192, {151, 69, 52, 36}, {8, 2, 2, 2}, 3, 3, 4, 15, 1450, 60};

struct LowDepthStructuredSecureRNSParameters {
    LowDepthSecureRNSParameters execution;
    std::uint32_t source_lwe_dimension;
    std::uint32_t source_lwe_hamming_weight;
    std::uint32_t structured_block_width;
    bool seeded_pbc_key;
};

inline constexpr LowDepthStructuredSecureRNSParameters
    qh_ss_structured_source_rns = {
        {8192, {151, 69, 52, 36}, {8, 2, 2, 2},
         1, 0, 4, 15, 1450, 60},
        1024, 64, 16, true};


// Prospective mapping of the source-screened Algorithm-3 row onto TFHEpp's
// 128-bit torus / Double-Decomposition representation.  The 105-bit native
// high level is embedded with a 23-bit scale, so its sigma'=0.75 noise maps to
// an absolute torus standard deviation of 0.75*2^23.  This preserves the
// paper-level relative noise and a coherent wide secret through a future
// coefficient-domain boundary.  The quadratic-hint two-component FullDD
// product is implemented and unit-tested; PBC/tree and boundary integration
// are still required before this parameter can be activated.
struct LowDepthDoubleDecompositionParameters : LowDepthParameters {
    std::uint32_t torus_modulus_bits;
    std::uint32_t native_high_modulus_bits;
    std::uint32_t embedding_scale_bits;
    std::uint32_t decomposition_limbs;
    std::uint32_t decomposition_limb_bits;
};

inline constexpr LowDepthDoubleDecompositionParameters
    binary_ntt_source_screened_dd = {
        binary_ntt_source_screened,
        128,
        105,
        23,
        8,
        16};

constexpr bool isStructuredOneHotCompatible(const Parameters &parameters)
{
    return parameters.secret_distribution ==
               SecretDistribution::StructuredSparseBinary &&
           parameters.lwe_hamming_weight != 0 &&
           parameters.lwe_dimension % parameters.lwe_hamming_weight == 0;
}

static_assert(isStructuredOneHotCompatible(structured_binary_gate_std128));
static_assert(isStructuredOneHotCompatible(structured_binary_gate_medium));
static_assert(isStructuredOneHotCompatible(
    structured_binary_function_std128));
static_assert(isStructuredOneHotCompatible(
    structured_binary_short_gate_std128));

// Executable structured-sparse gate-bootstrap candidate.  Its n, q, and h
// are the paper's Table 4 STD-128 setting, represented exactly in a 32-bit
// torus as Z_512 values shifted left by 23 bits.  The sigma = 3.2 choice is
// the paper's Section 7.2 LWE noise setting; Table 4 leaves it to OpenFHE.
// The current TFHEpp RLWE backend stores the N = 1024 accumulator in a 32-bit
// torus (Q = 2^32), rather than the paper's Q = 2^28.  This is an executable
// compatibility candidate, not a claim that the paper's RLWE security
// analysis transfers unchanged.
struct structured_binary_std128_lweparam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = 1;
    static constexpr std::uint32_t n = 1024;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t q = 1U << 9;
    static constexpr std::uint32_t qbit = 9;
    static constexpr std::uint32_t sparse_hamming_weight = 64;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    // Section 7.2 gives sigma = 3.2 for the q = 512 LWE layer.  This is
    // normalized to the enclosing 32-bit torus.
    static constexpr double alpha = 3.2 / q;
    static constexpr double α = alpha;
    using T = std::uint32_t;
    static constexpr std::make_signed_t<T> μ =
        static_cast<std::make_signed_t<T>>(q / 8)
        << (std::numeric_limits<T>::digits - qbit);
    static constexpr std::uint32_t plain_modulus = 8;
    static constexpr double Δ = static_cast<double>(q) / plain_modulus;
};

struct structured_binary_std128_gateparam {
    using domainP = structured_binary_std128_lweparam;
    using targetP = lvl1param;
    static constexpr std::uint32_t blind_rotation_modulus = domainP::q;
    static constexpr std::uint32_t Addends = 1;
};

}  // namespace shallowboot
