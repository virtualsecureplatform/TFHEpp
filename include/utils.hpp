#pragma once

#include <randen.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <params.hpp>
#include <random>
#include <immintrin.h>

namespace TFHEpp {
using namespace std;
static randen::Randen<uint64_t> generator;

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

inline uint64_t gaussian64(uint64_t center, double stdev)
{
    static const double _2p64 = pow(2., 64);
    normal_distribution<double> distribution(0., 1.0);
    const double val = stdev * distribution(generator) * _2p64;
    const uint64_t ival = static_cast<uint64_t>(val);
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

template <uint32_t Msize = 2 * DEF_nbar>
inline uint64_t modSwitchFromTorus64(uint32_t phase)
{
    uint64_t interv = ((1UL << 63) / Msize) * 2;  // width of each intervall
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
    for (int i = 0; i < N / 2; i+=4) {
        __m256d are,aim,bre,bim,resre,resim;
        are = _mm256_load_pd(a.data()+i);
        aim = _mm256_load_pd(a.data()+i+N/2);
        bre = _mm256_load_pd(b.data()+i);
        bim = _mm256_load_pd(b.data()+i+N/2);
        resre = _mm256_fmsub_pd(are,bre,_mm256_mul_pd(aim,bim));
        resim = _mm256_fmadd_pd(aim,bre,_mm256_mul_pd(are,bim));

        _mm256_store_pd(res.data()+i,resre);
        _mm256_store_pd(res.data()+i+N/2,resim);
    }
}

template <uint32_t N>
inline void FMAInFD(array<double, N> &res, const array<double, N> &a,
                    const array<double, N> &b)
{
    for (int i = 0; i < N / 2; i+=4) {
        __m256d are,aim,bre,bim,resre,resim;
        are = _mm256_load_pd(a.data()+i);
        aim = _mm256_load_pd(a.data()+i+N/2);
        bre = _mm256_load_pd(b.data()+i);
        bim = _mm256_load_pd(b.data()+i+N/2);
        resre = _mm256_load_pd(res.data()+i);
        resim = _mm256_load_pd(res.data()+i+N/2);
        resre = _mm256_add_pd(_mm256_fmsub_pd(are,bre,_mm256_mul_pd(aim,bim)),resre);
        resim = _mm256_add_pd(_mm256_fmadd_pd(aim,bre,_mm256_mul_pd(are,bim)),resim);

        _mm256_store_pd(res.data()+i,resre);
        _mm256_store_pd(res.data()+i+N/2,resim);
    }
}
}  // namespace TFHEpp