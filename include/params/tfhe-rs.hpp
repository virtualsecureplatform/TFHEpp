#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

// Parameters matched to tfhe-rs DEFAULT_PARAMETERS (boolean) and
// PARAMETERS_ERROR_PROB_2_POW_MINUS_165 (lvl2).
//
// References:
//   https://github.com/zama-ai/tfhe-rs/blob/main/tfhe/src/boolean/parameters/params.rs
//
// DEFAULT_PARAMETERS (used for lvl0 / lvl1 / lvl10param):
//   lwe_dimension:   805   (lvl0::n)
//   glwe_dimension:    3   (lvl1::k)
//   polynomial_size: 512   (lvl1::n)
//   pbs_base_log:     10   (lvl1::Bgbit / l-params)
//   pbs_level:          2   (lvl1::l)
//   ks_base_log:        3   (lvl10param::basebit)
//   ks_level:           5   (lvl10param::t)
//   lwe_std_dev:  5.8615896642671336e-06
//   glwe_std_dev: 9.315272083503367e-10
//   security: 132-bit, p-fail = 2^-64.344
//   encryption_key_choice: Small (ciphertexts live under the n=805 LWE key)
//
// PARAMETERS_ERROR_PROB_2_POW_MINUS_165 (used for lvl2 / lvl20param):
//   glwe_dimension:   2   (lvl2::k)
//   polynomial_size: 1024 (lvl2::n)
//   pbs_base_log:    10   (lvl2::Bgbit)
//   pbs_level:        2   (lvl2::l)
//   ks_base_log:      3   (lvl20param::basebit)
//   ks_level:         5   (lvl20param::t)
//   glwe_std_dev: 9.313225746198247e-10
//   p-fail = 2^-165.434

constexpr bool isternary = false;

struct lvl0param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 805;  // LWE dimension
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    // StandardDev from tfhe-rs DEFAULT_PARAMETERS lwe_noise_distribution
    static constexpr double α = 5.8615896642671336e-06;
    using T = uint32_t;  // Torus representation
    static constexpr std::make_signed_t<T> μ =
        1U << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

// Dummy
struct lvlhalfparam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 760;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -17);  // fresh noise
    using T = uint32_t;                                  // Torus representation
    static constexpr std::make_signed_t<T> μ =
        1U << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

// Matched to tfhe-rs DEFAULT_PARAMETERS GLWE:
//   glwe_dimension=3, polynomial_size=512, pbs_base_log=10, pbs_level=2
//   glwe_std_dev = 9.315272083503367e-10
struct lvl1param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;  // binary keys (matching tfhe-rs)
    static constexpr std::uint32_t nbit =
        9;  // dimension must be a power of 2 for ease of polynomial
            // multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // 512
    static constexpr std::uint32_t k = 3;
    static constexpr std::uint32_t lₐ = 2;   // pbs_level = 2
    static constexpr std::uint32_t l = 2;
    static constexpr std::uint32_t Bgₐbit = 10;  // pbs_base_log = 10
    static constexpr std::uint32_t Bgbit = 10;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    // StandardDev from tfhe-rs DEFAULT_PARAMETERS glwe_noise_distribution
    static const inline double α = 9.315272083503367e-10;
    using T = uint64_t;  // Torus representation
    static constexpr std::make_signed_t<T> μ = 1ULL << 61;
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        2 * static_cast<double>(1ULL << (std::numeric_limits<T>::digits - 1)) /
        plain_modulus;
    // Double Decomposition (bivariate representation) parameters
    // For now, set to trivial values (no actual second decomposition)
    static constexpr std::uint32_t l̅ = 1;  // auxiliary decomposition levels
    static constexpr std::uint32_t l̅ₐ = l̅;
    static constexpr std::uint32_t B̅gbit =
        std::numeric_limits<T>::digits;  // full coefficient width
    static constexpr std::uint32_t B̅gₐbit = B̅gbit;
};

// Annihilate-Homomorphism variant of lvl1param.
// Uses the same ring (k, N) but different decomposition parameters, as in
// the 128bit.hpp reference.
struct AHlvl1param {
    using baseP = lvl1param;
    static constexpr int32_t key_value_max = baseP::key_value_max;
    static constexpr int32_t key_value_min = baseP::key_value_min;
    static constexpr std::uint32_t nbit = baseP::nbit;
    static constexpr std::uint32_t n = baseP::n;
    static constexpr std::uint32_t k = baseP::k;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 5;
    static constexpr std::uint32_t Bgₐbit = 5;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist = baseP::errordist;
    static const inline double α = baseP::α;
    using T = typename baseP::T;
    static constexpr std::make_signed_t<T> μ = baseP::μ;
    static constexpr uint32_t plain_modulus = baseP::plain_modulus;
    static constexpr double Δ = baseP::Δ;
    static constexpr std::uint32_t l̅ = baseP::l̅;
    static constexpr std::uint32_t l̅ₐ = baseP::l̅ₐ;
    static constexpr std::uint32_t B̅gbit = baseP::B̅gbit;
    static constexpr std::uint32_t B̅gₐbit = baseP::B̅gₐbit;
};

// Matched to tfhe-rs PARAMETERS_ERROR_PROB_2_POW_MINUS_165 GLWE:
//   glwe_dimension=2, polynomial_size=1024, pbs_base_log=10, pbs_level=2
//   glwe_std_dev = 9.313225746198247e-10
struct lvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static const std::uint32_t nbit = 10;  // dimension must be a power of 2 for
                                           // ease of polynomial multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // 1024
    static constexpr std::uint32_t k = 2;
    static constexpr std::uint32_t lₐ = 2;   // pbs_level = 2
    static constexpr std::uint32_t l = 2;
    static constexpr std::uint32_t Bgₐbit = 10;  // pbs_base_log = 10
    static constexpr std::uint32_t Bgbit = 10;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    // StandardDev from tfhe-rs PARAMETERS_ERROR_PROB_2_POW_MINUS_165
    // glwe_noise_distribution
    static const inline double α = 9.313225746198247e-10;
    using T = uint64_t;  // Torus representation
    static constexpr std::make_signed_t<T> μ = 1ULL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ = μ;
    // Double Decomposition (bivariate representation) parameters
    // For now, set to trivial values (no actual second decomposition)
    static constexpr std::uint32_t l̅ = 1;  // auxiliary decomposition levels
    static constexpr std::uint32_t l̅ₐ = l̅;
    static constexpr std::uint32_t B̅gbit =
        std::numeric_limits<T>::digits;  // full coefficient width
    static constexpr std::uint32_t B̅gₐbit = B̅gbit;
};

// Annihilate-Homomorphism variant of lvl2param.
struct AHlvl2param {
    using baseP = lvl2param;
    static constexpr int32_t key_value_max = baseP::key_value_max;
    static constexpr int32_t key_value_min = baseP::key_value_min;
    static const std::uint32_t nbit = baseP::nbit;
    static constexpr std::uint32_t n = baseP::n;
    static constexpr std::uint32_t k = baseP::k;
    static constexpr std::uint32_t lₐ = 5;
    static constexpr std::uint32_t l = 5;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist = baseP::errordist;
    static const inline double α = baseP::α;
    using T = typename baseP::T;
    static constexpr std::make_signed_t<T> μ = baseP::μ;
    static constexpr uint32_t plain_modulus = baseP::plain_modulus;
    static constexpr double Δ = baseP::Δ;
    static constexpr std::uint32_t l̅ = baseP::l̅;
    static constexpr std::uint32_t l̅ₐ = baseP::l̅ₐ;
    static constexpr std::uint32_t B̅gbit = baseP::B̅gbit;
    static constexpr std::uint32_t B̅gₐbit = baseP::B̅gₐbit;
};

// Dummy
struct lvl3param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static const std::uint32_t nbit = 13;  // dimension must be a power of 2 for
    // ease of polynomial multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // 8192
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -47);  // fresh noise
    using T = uint64_t;                                  // Torus representation
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulusbit = 31;
    static constexpr uint64_t plain_modulus = 1ULL << plain_modulusbit;
    static constexpr double Δ = 1ULL << (64 - plain_modulusbit - 1);
    // Double Decomposition (bivariate representation) parameters
    // For now, set to trivial values (no actual second decomposition)
    static constexpr std::uint32_t l̅ = 1;  // auxiliary decomposition levels
    static constexpr std::uint32_t l̅ₐ = l̅;
    static constexpr std::uint32_t B̅gbit =
        std::numeric_limits<T>::digits;  // full coefficient width
    static constexpr std::uint32_t B̅gₐbit = B̅gbit;
};

// Dummy
using lvl4param = lvl3param;

// Key Switching parameters

// Matched to tfhe-rs DEFAULT_PARAMETERS:
//   ks_base_log=3, ks_level=5
struct lvl10param {
    static constexpr std::uint32_t t = 5;  // ks_level = 5
    static constexpr std::uint32_t basebit =
        3;  // ks_base_log = 3
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;  // key noise
    using domainP = lvl1param;
    using targetP = lvl0param;
};

// Dummy
struct lvl1hparam {
    static constexpr std::uint32_t t =
        10;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        3;  // how many bit should be encrypted in keyswitching key
    static const inline double α = lvlhalfparam::α;  // key noise
    using domainP = lvl1param;
    using targetP = lvlhalfparam;
};

struct lvl11param {
    static constexpr std::uint32_t t = 6;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        4;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl1param;
    using targetP = lvl1param;
};

// Matched to tfhe-rs PARAMETERS_ERROR_PROB_2_POW_MINUS_165:
//   ks_base_log=3, ks_level=5
struct lvl20param {
    static constexpr std::uint32_t t = 5;  // ks_level = 5
    static constexpr std::uint32_t basebit =
        3;  // ks_base_log = 3
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvl0param;
};

// Dummy
struct lvl2hparam {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvlhalfparam::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvlhalfparam;
};

struct lvl21param {
    static constexpr std::uint32_t t = 24;  // number of addition in
                                            // keyswitching
    static constexpr std::uint32_t basebit =
        1;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvl1param;
};

struct lvl22param {
    static constexpr std::uint32_t t = 8;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        4;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl2param::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvl2param;
};

// Dummy
struct lvl31param {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl3param;
    using targetP = lvl1param;
};

// Dummy
using lvl41param = lvl31param;
