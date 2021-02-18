#pragma once

#include <randen.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <params.hpp>
#include <random>

namespace TFHEpp {
static randen::Randen<uint64_t> generator;

// Double to Torus(32bit fixed-point number)
inline uint32_t dtot32(double d)
{
    return int32_t(int64_t((d - int64_t(d)) * (1LL << 32)));
}

// Modular Gaussian Distribution over Torus(32bit fixed-point number)
inline uint32_t gaussian32(uint32_t message, double α)
{
    std::normal_distribution<double> distribution(
        0., α);  // TODO: can we create a global distrib of param 1 and multiply
                 // by sigma?
    double err = distribution(generator);
    return message + dtot32(err);
}

// Modular Gaussian Distribution over Torus(64bit fixed-point number)
inline uint64_t gaussian64(uint64_t center, double stdev)
{
    static const double _2p64 = pow(2., 64);
    std::normal_distribution<double> distribution(0., 1.0);
    const double val = stdev * distribution(generator) * _2p64;
    const uint64_t ival = static_cast<uint64_t>(val);
    return ival + center;
}

template <uint32_t Mbit = DEF_Nbit + 1>
inline uint32_t modSwitchFromTorus32(uint32_t phase)
{
    return (phase + (1U << (31 - Mbit))) >> (32 - Mbit);
}

template <uint32_t Msize = 2 * DEF_nbar>
inline uint64_t modSwitchFromTorus64(uint32_t phase)
{
    uint64_t interv = ((1ULL << 63) / Msize) * 2;  // width of each intervall
    uint64_t half_interval = interv / 2;  // begin of the first intervall

    // Mod Switching (as in modSwitchFromTorus32)
    uint64_t temp =
        (static_cast<uint64_t>(phase) << 32) + half_interval;  // RIVEDI
    return temp / interv;
}

template <uint32_t N>
inline void MulInFD(array<double, N> &res, const array<double, N> &a,
                    const array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        double aimbim = a[i + N / 2] * b[i + N / 2];
        double arebim = a[i] * b[i + N / 2];
        res[i] = std::fma(a[i], b[i], -aimbim);
        res[i + N / 2] = std::fma(a[i + N / 2], b[i], arebim);
    }
}

template <uint32_t N>
inline void FMAInFD(array<double, N> &res, const array<double, N> &a,
                    const array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        res[i] = std::fma(a[i + N / 2], b[i + N / 2], -res[i]);
        res[i] = std::fma(a[i], b[i], -res[i]);
        res[i + N / 2] = std::fma(a[i], b[i + N / 2], res[i + N / 2]);
        res[i + N / 2] = std::fma(a[i + N / 2], b[i], res[i + N / 2]);
    }
}
}  // namespace TFHEpp
