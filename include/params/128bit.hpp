#pragma once

#import <cmath>
#import <cstdint>
#import <limits>

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
    static constexpr std::uint32_t lₐ = 3;
    static constexpr std::uint32_t l = 3;
    static constexpr std::uint32_t Bgbit = 6;
    static constexpr std::uint32_t Bgₐbit = 6;
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
};

struct Annihilatelvl1param {
    static constexpr int32_t key_value_max = lvl1param::key_value_max;
    static constexpr int32_t key_value_min = lvl1param::key_value_min;
    static constexpr std::uint32_t nbit = lvl1param::nbit;
    static constexpr std::uint32_t n = lvl1param::n;  // dimension
    static constexpr std::uint32_t k = lvl1param::k;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 5;
    static constexpr std::uint32_t Bgₐbit = 5;
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
};

struct lvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static const std::uint32_t nbit = 11;  // dimension must be a power of 2 for
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
    using T = uint64_t;                                 // Torus representation
    static constexpr std::make_signed_t<T> μ = 1LL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << (std::numeric_limits<T>::digits - 4));
};

struct lvl3param {
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
    using T = uint64_t;                                 // Torus representation
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulusbit = 31;
    static constexpr uint64_t plain_modulus = 1ULL << plain_modulusbit;
    static constexpr double Δ = 1ULL << (64 - plain_modulusbit - 1);
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