#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <params.hpp>
#include <random>

namespace TFHEpp {
using namespace std;
static random_device generator;

inline uint32_t dtot32(double d)
{
    return int32_t(int64_t((d - int64_t(d)) * (1L << 32)));
}

inline uint32_t gaussian32(uint32_t message, double α)
{
    normal_distribution<double> distribution(
        0., α);  // TODO: can we create a global distrib of param 1 and multiply
                 // by sigma?
    double err = distribution(generator);
    return message + dtot32(err);
}

inline int64_t gaussian64(int64_t center, double stdev)
{
    static const double _2p64 = pow(2., 64);
    normal_distribution<double> distribution(0., 1.0);
    double val = stdev * distribution(generator) * _2p64;
    int64_t ival = (int64_t)val;
    return ival + center;
}

template <uint32_t Msize = 2 * DEF_N>
inline uint32_t modSwitchFromTorus32(uint32_t phase)
{
    uint64_t interv = ((1UL << 63) / Msize) * 2;  // width of each intervall
    uint64_t half_interval = interv / 2;  // begin of the first intervall
    uint64_t phase64 = (uint64_t(phase) << 32) + half_interval;
    // floor to the nearest multiples of interv
    return static_cast<uint32_t>(phase64 / interv);
}

template <uint32_t N>
inline void MulInFD(array<double, N> &res, const array<double, N> &a,
                    const array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        double aimbim = a[i + N / 2] * b[i + N / 2];
        double arebim = a[i] * b[i + N / 2];
        res[i] = a[i] * b[i] - aimbim;
        res[i + N / 2] = a[i + N / 2] * b[i] + arebim;
    }
}

template <uint32_t N>
inline void FMAInFD(array<double, N> &res, const array<double, N> &a,
                    const array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        res[i] = a[i + N / 2] * b[i + N / 2] - res[i];
        res[i] = a[i] * b[i] - res[i];
        res[i + N / 2] += a[i] * b[i + N / 2];
        res[i + N / 2] += a[i + N / 2] * b[i];
    }
}
}  // namespace TFHEpp