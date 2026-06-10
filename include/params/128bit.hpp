#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

struct lvl0param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 636;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = 0.000'092'511'997'467'675'6;  // fresh noise
    using T = uint16_t;  // Torus representation
    static constexpr std::make_signed_t<T> μ =
        1LL << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvlhalfparam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 760;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -17);  // fresh noise
    using T = uint32_t;                                 // Torus representation
    static constexpr T μ = 1U << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 32;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvl1param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit =
        10;  // dimension must be a power of 2 for ease of polynomial
             // multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 3;
    static constexpr std::uint32_t lₐ = l;
    static constexpr std::uint32_t Bgbit = 6;
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -25);  // fresh noise
    using T = uint32_t;                                 // Torus representation
    static constexpr std::make_signed_t<T> μ = 1 << 29;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
    // Double Decomposition (bivariate representation) parameters
    // For now, set to trivial values (no actual second decomposition)
    static constexpr std::uint32_t l̅ = 1;  // auxiliary decomposition levels
    static constexpr std::uint32_t l̅ₐ = l̅;
    static constexpr std::uint32_t B̅gbit =
        std::numeric_limits<T>::digits;  // full coefficient width
    static constexpr std::uint32_t B̅gₐbit = B̅gbit;
};

struct AHlvl1param {
    using baseP = lvl1param;
    static constexpr int32_t key_value_max = baseP::key_value_max;
    static constexpr int32_t key_value_min = baseP::key_value_min;
    static constexpr std::uint32_t nbit = baseP::nbit;
    static constexpr std::uint32_t n = baseP::n;  // dimension
    static constexpr std::uint32_t k = baseP::k;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 5;
    static constexpr std::uint32_t Bgₐbit = 5;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist = baseP::errordist;
    static const inline double α = baseP::α;  // fresh noise
    using T = typename baseP::T;              // Torus representation
    static constexpr std::make_signed_t<T> μ = baseP::μ;
    static constexpr uint32_t plain_modulus = baseP::plain_modulus;
    static constexpr double Δ = baseP::Δ;
    // Double Decomposition parameters inherited from baseP
    static constexpr std::uint32_t l̅ = baseP::l̅;
    static constexpr std::uint32_t l̅ₐ = baseP::l̅ₐ;
    static constexpr std::uint32_t B̅gbit = baseP::B̅gbit;
    static constexpr std::uint32_t B̅gₐbit = baseP::B̅gₐbit;
};

struct lvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static const std::uint32_t nbit = 11;  // dimension must be a power of 2 for
                                           // ease of polynomial multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t lₐ = l;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -51);  // fresh noise
    using T = uint64_t;                                 // Torus representation
    static constexpr std::make_signed_t<T> μ = 1LL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << (std::numeric_limits<T>::digits - 4));
    // Double Decomposition (bivariate representation) parameters
    // For now, set to trivial values (no actual second decomposition)
    static constexpr std::uint32_t l̅ = 1;  // auxiliary decomposition levels
    static constexpr std::uint32_t l̅ₐ = l̅;
    static constexpr std::uint32_t B̅gbit =
        std::numeric_limits<T>::digits;  // full coefficient width
    static constexpr std::uint32_t B̅gₐbit = B̅gbit;
};

struct AHlvl2param {
    using baseP = lvl2param;
    static constexpr int32_t key_value_max = baseP::key_value_max;
    static constexpr int32_t key_value_min = baseP::key_value_min;
    static constexpr std::uint32_t nbit = baseP::nbit;
    static constexpr std::uint32_t n = baseP::n;  // dimension
    static constexpr std::uint32_t k = baseP::k;
    static constexpr std::uint32_t lₐ = 5;
    static constexpr std::uint32_t l = 5;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist = baseP::errordist;
    static const inline double α = baseP::α;  // fresh noise
    using T = typename baseP::T;              // Torus representation
    static constexpr std::make_signed_t<T> μ = baseP::μ;
    static constexpr uint32_t plain_modulus = baseP::plain_modulus;
    static constexpr double Δ = baseP::Δ;
    // Double Decomposition parameters inherited from baseP
    static constexpr std::uint32_t l̅ = baseP::l̅;
    static constexpr std::uint32_t l̅ₐ = baseP::l̅ₐ;
    static constexpr std::uint32_t B̅gbit = baseP::B̅gbit;
    static constexpr std::uint32_t B̅gₐbit = baseP::B̅gₐbit;
};

// lvl3param with 128-bit Torus and non-trivial Double Decomposition
// Double decomposition structure:
//   - Primary decomposition (l, Bgbit): Decomposes plaintext by μ in TRGSW gadget
//   - Auxiliary decomposition (l̅, B̅gbit): Decomposes TRLWE ciphertext coefficients
//     in the external product. Must cover full 128-bit coefficient.
// Numerical safety for FFT (double mantissa): Bgbit + B̅gbit + nbit + 3 < 53.
// Full limb cover: l̅ * B̅gbit = 128 (same for nonce part).
struct lvl3param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static const std::uint32_t nbit = 12;  // dimension must be a power of 2 for
    // ease of polynomial multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // dimension = 4096
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 21;
    static constexpr std::uint32_t Bgₐbit = 21;
    static constexpr uint32_t Bg = 1U << Bgbit;
    static constexpr uint32_t Bgₐ = 1U << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -105);  // fresh noise
    using T = __uint128_t;                               // Torus representation
    static constexpr T μ = static_cast<T>(1) << 125;
    static constexpr uint32_t plain_modulusbit = 31;
    static constexpr T plain_modulus = static_cast<T>(1) << plain_modulusbit;
    static constexpr double Δ =
        static_cast<double>(static_cast<T>(1) << (128 - plain_modulusbit - 1));
    // Double Decomposition (bivariate representation) parameters
    // Auxiliary decomposition must cover full 128-bit ciphertext coefficients
    // l̅ * B̅gbit = 8 * 16 = 128 bits
    static constexpr std::uint32_t l̅ = 8;     // auxiliary decomposition levels
    static constexpr std::uint32_t l̅ₐ = 8;
    static constexpr std::uint32_t B̅gbit = 16;    // 2^16 limb base (covers 128-bit T with l̅=8)
    static constexpr std::uint32_t B̅gₐbit = 16;
};

// lvl3simdparam: like lvl3param but with prime plain_modulus enabling SIMD slots.
// t = 114689 = 7*2^14+1, prime, t ≡ 1 (mod 8192 = 2n), supports n=4096 SIMD slots.
// Primitive 2n-th root of unity: ψ = 80720, ψ_inv = 7887, n_inv = 114661 (all mod t).
// Security: identical ring/noise/DD to lvl3param.
// FFT safety: Bgbit + B̅gbit + nbit + 3 = 21 + 16 + 12 + 3 = 52 < 53 ✓
//
// BFV scaling: Δ = floor(Q/t) ≈ 2^111.19 (NOT a power of 2).
// Encrypt uses floor(m·Q/t) per coefficient (≤ 1 unit rounding error).
// Decrypt uses round(phase·t/Q) (exact BFV decoding).
struct lvl3simdparam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 12;
    static constexpr std::uint32_t n = 1 << nbit;  // 4096
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 21;
    static constexpr std::uint32_t Bgₐbit = 21;
    static constexpr uint32_t Bg = 1U << Bgbit;
    static constexpr uint32_t Bgₐ = 1U << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -105);  // same as lvl3param
    using T = __uint128_t;
    static constexpr T μ = static_cast<T>(1) << 125;
    // plain_modulusbit kept for DD decomposition compatibility (used by existing
    // code paths like ExternalProduct).  For BFV scaling, use delta_int below.
    static constexpr uint32_t plain_modulusbit = 18;
    static constexpr T plain_modulus = static_cast<T>(114689); // prime, 7*2^14+1
    static constexpr double Δ =
        static_cast<double>(static_cast<T>(1) << (128 - plain_modulusbit));
    // Δ_int = floor(Q/t) = floor(2^128 / 114689).
    // Computed as: (2^128-1)/t + correction.  Since t is odd and does not
    // divide 2^128, floor(2^128/t) = floor((2^128-1)/t).
    static constexpr T delta_int =
        static_cast<T>(-1) / plain_modulus;  // UINT128_MAX / t
    // Remainder: Q mod t.  Used for exact floor(m·Q/t) encoding.
    // Q_mod_t = (UINT128_MAX % t) + 1, handling the +1 from 2^128 = UINT128_MAX+1.
    // If UINT128_MAX % t == t-1, then Q mod t = 0 (t divides 2^128), but t is odd so this won't happen.
    static constexpr uint64_t Q_mod_t =
        static_cast<uint64_t>(static_cast<T>(-1) % plain_modulus) + 1;
    static constexpr uint64_t bfv_bootstrap_digit_error_bound = 15;
    static constexpr int bfv_bootstrap_linear_bsgs_step = 64;
    // DD parameters — identical to lvl3param
    static constexpr std::uint32_t l̅ = 8;
    static constexpr std::uint32_t l̅ₐ = 8;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;
};

struct lvl4param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static const std::uint32_t nbit = 13;  // dimension must be a power of 2 for
    // ease of polynomial multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -51);  // fresh noise
    using T = __uint128_t;                                 // Torus representation
    static constexpr T μ = static_cast<T>(1) << 125;
    static constexpr uint32_t plain_modulusbit = 31;
    static constexpr T plain_modulus = static_cast<T>(1) << plain_modulusbit;
    static constexpr double Δ =
        static_cast<double>(static_cast<T>(1) << (128 - plain_modulusbit - 1));
    // Double Decomposition (bivariate representation) parameters
    // Trivial values (no actual second decomposition)
    static constexpr std::uint32_t l̅ = 8;  // auxiliary decomposition levels
    static constexpr std::uint32_t l̅ₐ = 8;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;
};

// lvl5param: BFV-oriented DD multi-limb torus scaffold.
//
// Q = 2^448 is represented as 7 little-endian 64-bit limbs.  The SIMD
// plaintext modulus is p = 786433 = 3*2^18+1, prime and congruent to 1 mod
// 2n for n = 2^14.  This makes the parameter usable for BFV SIMD slot
// experiments once the remaining multi-limb FFT/key paths are wired.
struct lvl5param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 14;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 5;
    static constexpr std::uint32_t l = 5;
    static constexpr std::uint32_t Bgbit = 19;
    static constexpr std::uint32_t Bgₐbit = 19;
    using T = MultiLimbUInt<7>;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bgₐ = T{1} << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -425);
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint64_t plain_modulus_u64 = 786433;
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{plain_modulus_u64};
    static constexpr double Δ = 0.0;  // BFV code should use delta_int.
    static constexpr T delta_int =
        std::numeric_limits<T>::max() / plain_modulus_u64;
    static constexpr uint64_t Q_mod_t =
        (std::numeric_limits<T>::max() % plain_modulus_u64) + 1;
    static constexpr uint64_t bfv_bootstrap_digit_error_bound = 15;
    static constexpr int bfv_bootstrap_linear_bsgs_step = 128;
    static constexpr std::uint32_t l̅ = 28;
    static constexpr std::uint32_t l̅ₐ = 28;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;

    static constexpr uint64_t simd_modulus = plain_modulus_u64;
    static constexpr uint64_t simd_psi = 585160;
    static constexpr uint64_t simd_psi_inv = 253771;
    static constexpr uint64_t simd_n_inv = 786385;
};

// lvl6param: n = 2^15 multi-limb torus scaffold for dense CKKS bootstrap work.
//
// This intentionally keeps lvl5param at n = 2^14 and provides the 2^15 ring as
// the next level.  Q = 2^896 is represented as 14 little-endian 64-bit limbs,
// which is enough for the tuned 896-bit CKKS bootstrap schedule while avoiding
// the impractical absolute EvalMod noise created by the earlier 1108-bit
// schedule.  The CKKS encryption noise is tuned close to the small integer
// Gaussian noise used by RNS CKKS libraries: at logQ = 896, α = 2^-886 gives
// σ = 2^10, with a 3.2 floor at lower active levels.
struct lvl6param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 15;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 5;
    static constexpr std::uint32_t l = 5;
    // CKKS bootstrapping EvalMod is sensitive to pre-EvalMod automorphism
    // noise.  A 7-bit auxiliary base gives the lvl6 linear key switches more
    // precision for the practical dense CKKS schedule while still covering the
    // 896-bit torus exactly at compile time.
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t Bgₐbit = 16;
    using T = MultiLimbUInt<14>;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bgₐ = T{1} << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -886);
    // Low active CKKS levels would otherwise round α*2^LogQ to zero.  This
    // floor keeps direct lvl6 CKKS encryption/keygen non-degenerate while the
    // high-level bootstrap keys still use the security-tuned relative α above.
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint64_t plain_modulus_u64 = 786433;
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{plain_modulus_u64};
    static constexpr double Δ = 0.0;
    static constexpr T delta_int =
        std::numeric_limits<T>::max() / plain_modulus_u64;
    static constexpr uint64_t Q_mod_t =
        (std::numeric_limits<T>::max() % plain_modulus_u64) + 1;
    static constexpr uint64_t bfv_bootstrap_digit_error_bound = 15;
    static constexpr int bfv_bootstrap_linear_bsgs_step = 128;
    static constexpr std::uint32_t l̅ = 128;
    static constexpr std::uint32_t l̅ₐ = 128;
    static constexpr std::uint32_t B̅gbit = 7;
    static constexpr std::uint32_t B̅gₐbit = 7;

    static constexpr uint64_t simd_modulus = plain_modulus_u64;
    static constexpr uint64_t simd_psi = 108788;
    static constexpr uint64_t simd_psi_inv = 295516;
    static constexpr uint64_t simd_n_inv = 786409;
};

// Key Switching parameters
struct lvl10param {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;  // key noise
    using domainP = lvl1param;
    using targetP = lvl0param;
};

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

struct lvl20param {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvl0param;
};

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
    static constexpr std::uint32_t t = 8;  // number of addition in
                                           // keyswitching
    static constexpr std::uint32_t basebit =
        3;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvl1param;
};

struct lvl22param {
    static constexpr std::uint32_t t =
        38;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        1;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl2param::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvl2param;
};

struct lvl31param {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl3param;
    using targetP = lvl1param;
};

struct lvl41param {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl4param;
    using targetP = lvl1param;
};
