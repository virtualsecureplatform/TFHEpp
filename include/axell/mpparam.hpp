#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

#ifdef USE_80BIT_SECURITY
// IGG16.hpp
struct lvlMparam {
    static constexpr std::uint32_t nbit = 10;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 2;
    static constexpr std::uint32_t Bgbit = 10;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static const inline double α = 3.73e-9;
    using T = uint32_t;
    // static constexpr T μ = 0x15555555;
    static constexpr T μ =
        (1ULL << std::numeric_limits<typename lvlMparam::T>::digits) / 12;
};
#elif defined(USE_CGGI19)
// IGG19.hpp
struct lvlMparam {
    static constexpr std::uint32_t nbit = 10;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 3;
    static constexpr std::uint32_t Bgbit = 7;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static const inline double α = std::pow(2.0, -25);
    using T = uint32_t;
    // static constexpr T μ = 0x15555555;
    static constexpr T μ =
        (1ULL << std::numeric_limits<typename lvlMparam::T>::digits) / 12;
    // dummy
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};
#else
// 128bit.hpp
struct lvlMparam {
    static constexpr std::uint32_t nbit = 10;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 3;
    static constexpr std::uint32_t Bgbit = 6;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static const inline double α = std::pow(2.0, -25);
    using T = uint32_t;
    // static constexpr T μ = 0x15555555;
    static constexpr T μ =
        (1ULL << std::numeric_limits<typename lvlMparam::T>::digits) / 12;
    // dummy
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};
#endif

struct lvlM0param {
    static constexpr std::uint32_t t = 7;
    static constexpr std::uint32_t basebit = 2;
    static const inline double α = lvl0param::α;
    using domainP = lvlMparam;
    using targetP = lvl0param;
};

struct lvl0Mparam {
    using domainP = lvl0param;
    using targetP = lvlMparam;
};
