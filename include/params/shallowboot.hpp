#pragma once

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

inline constexpr LowDepthParameters binary_ntt_std128 = {
    {SecretDistribution::BinaryNTT, 1450, 9, 29, 4096, 105, 12},
    3.2,
    0.75,
    3,
    4,
};

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
