#pragma once

#include <cmath>
#include <cstdint>

constexpr bool isternary = false;

struct lvl0param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 630;
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -15);
    using T = uint32_t;
    static constexpr T μ = 1U << 29;
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
    using T = uint32_t;                                 // Torus representation
    static constexpr T μ = 1U << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvl1param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr std::uint32_t nbit = 10;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 3;
    static constexpr std::uint32_t Bgbit = 7;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -25);
    using T = uint32_t;
    static constexpr T μ = 1U << 29;
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static const std::uint32_t nbit = 11;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -44);
    using T = uint64_t;
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ = μ;
};

// Dummy
struct lvl3param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static const std::uint32_t nbit = 13;  // dimension must be a power of 2 for
    // ease of polynomial multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -47);  // fresh noise
    using T = uint64_t;                                 // Torus representation
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulusbit = 31;
    static constexpr uint64_t plain_modulus = 1ULL << plain_modulusbit;
    static constexpr double Δ = 1ULL << (64 - plain_modulusbit - 1);
};

// Dummy
struct lvl11param {
    static constexpr std::uint32_t t = 0;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        0;  // how many bit should be encrypted in keyswitching key
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;  // key noise
    using domainP = lvl1param;
    using targetP = lvl0param;
};

struct lvl10param {
    static constexpr std::uint32_t t = 8;
    static constexpr std::uint32_t basebit = 2;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;
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

struct lvl21param {
    static constexpr std::uint32_t t = 10;
    static constexpr std::uint32_t basebit = 3;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2, -29);
    using domainP = lvl2param;
    using targetP = lvl1param;
};

// Dummy
struct lvl20param {
    static constexpr std::uint32_t t = 0;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        0;  // how many bit should be encrypted in keyswitching key
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

// Dummy
struct lvl22param {
    static constexpr std::uint32_t t = 0;
    static constexpr std::uint32_t basebit = 0;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = lvl2param::α;
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