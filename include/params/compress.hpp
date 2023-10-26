#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

constexpr bool isternary = false;

struct lvl0param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 636;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static constexpr inline double α =
        0.000'092'511'997'467'675'6;  // fresh noise, 2^{-13.4}
    using T = uint32_t;               // Torus representation
    static constexpr T μ = 1U << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvl1param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit =
        9;  // dimension must be a power of 2 for ease of polynomial
            // multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // dimension
    static constexpr std::uint32_t k = 2;
    static constexpr std::uint32_t l = 2;
    static constexpr std::uint32_t Bgbit = 8;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::CenteredBinomial;
    static constexpr uint η = 3;
    using T = uint32_t;  // Torus representation
    static constexpr T q = 40960001;
    static constexpr uint qbit = 27;
    static constexpr T μ = 1U << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static const std::uint32_t nbit = 9;  // dimension must be a power of 2 for
                                          // ease of polynomial multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // dimension
    static constexpr std::uint32_t k = 3;
    static constexpr std::uint32_t l = 3;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::CenteredBinomial;
    static constexpr uint η = 3;
    using T = uint64_t;  // Torus representation
    static constexpr T q = 1ULL << 48;
    static constexpr uint qbit = 48;
    static constexpr T μ = q / 8;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ = μ;
};

// Key Switching parameters
struct lvl10param {
    static constexpr std::uint32_t t = 5;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    using domainP = lvl1param;
    using targetP = lvl0param;
};

struct lvl11param {
    static constexpr std::uint32_t t = 6;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        4;  // how many bit should be encrypted in keyswitching key
    using domainP = lvl1param;
    using targetP = lvl1param;
};

struct lvl20param {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    using domainP = lvl2param;
    using targetP = lvl0param;
};

struct lvl21param {
    static constexpr std::uint32_t t = 8;  // number of addition in
                                           // keyswitching
    static constexpr std::uint32_t basebit =
        3;  // how many bit should be encrypted in keyswitching key
    using domainP = lvl2param;
    using targetP = lvl1param;
};

struct lvl22param {
    static constexpr std::uint32_t t =
        38;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        1;  // how many bit should be encrypted in keyswitching key
    using domainP = lvl2param;
    using targetP = lvl2param;
};
