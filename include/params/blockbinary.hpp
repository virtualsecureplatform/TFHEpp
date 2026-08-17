#pragma once

#include <cstdint>
#include <limits>

// Block-binary changes only the lvl0 secret distribution.  Keeping the
// 128-bit parameter topology from lvlhalf onward lets the extended schemes
// (CLPX, BFV, and CKKS) share their tested parameter assumptions.
constexpr bool isternary = false;

struct lvl0param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 630;
    static constexpr std::uint32_t ell = 2;
    static_assert(n % ell == 0,
                  "block-binary dimension must be divisible by ell");
    static constexpr std::uint32_t k = 1;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static constexpr inline double α =
        0.000'092'511'997'467'675'6;  // fresh noise, 2^{-13.4}
    using T = uint16_t;
    static constexpr std::make_signed_t<T> μ =
        1 << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

#define TFHEPP_BLOCK_BINARY_LVL0
#include "128bit.hpp"
#undef TFHEPP_BLOCK_BINARY_LVL0
