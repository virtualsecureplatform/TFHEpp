#pragma once

#include <cmath>
#include <cstdint>

struct lvl0param {
    static constexpr std::uint32_t n = 500;
    static constexpr std::uint32_t k = 1;
    static const inline double α = 2.44e-5;
    using T = uint32_t;
    static constexpr T μ = 1U << 29;
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvl1param {
    static constexpr std::uint32_t nbit = 10;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 2;
    static constexpr std::uint32_t Bgbit = 10;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static const inline double α = 3.73e-9;
    using T = uint32_t;
    static constexpr T μ = 1U << 29;
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

struct lvl2param {
    static const std::uint32_t nbit = 11;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static const inline double α = std::pow(2.0, -44);
    using T = uint64_t;
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ = μ;
};

struct lvl10param {
    static constexpr std::uint32_t t = 8;
    static constexpr std::uint32_t basebit = 2;
    static const inline double α = lvl0param::α;
    using domainP = lvl1param;
    using targetP = lvl0param;
};

// dummy
struct lvl11param {
    static constexpr std::uint32_t t = 0;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        0;  // how many bit should be encrypted in keyswitching key
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl1param;
    using targetP = lvl1param;
};

struct lvl21param {
    static constexpr std::uint32_t t = 10;
    static constexpr std::uint32_t basebit = 3;
    static const inline double α = std::pow(2, -31);
    using domainP = lvl2param;
    using targetP = lvl1param;
};

// Dummy
struct lvl20param {
    static constexpr std::uint32_t t = 0;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        0;  // how many bit should be encrypted in keyswitching key
    static const inline double α = lvl0param::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvl0param;
};

// Dummy
struct lvl22param {
    static constexpr std::uint32_t t = 0;
    static constexpr std::uint32_t basebit = 0;
    static const inline double α = lvl2param::α;
    using domainP = lvl2param;
    using targetP = lvl2param;
};