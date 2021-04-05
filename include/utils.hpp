#pragma once

#include <randen.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <params.hpp>
#include <random>

namespace TFHEpp {
static thread_local std::random_device trng;
static thread_local randen::Randen<uint64_t> generator(trng());

// https://qiita.com/saka1_p/items/e8c4dfdbfa88449190c5
template <typename T>
constexpr bool false_v = false;

// Double to Torus(32bit fixed-point number)
inline uint32_t dtot32(double d)
{
    return int32_t(int64_t((d - int64_t(d)) * (1LL << 32)));
}

// Modular Gaussian Distribution over Torus
template <class P>
typename P::T ModularGaussian(typename P::T center, double stdev);

template <class P>
typename P::T modSwitchFromTorus(uint32_t phase);

template <uint32_t N>
void MulInFD(array<double, N> &res, const array<double, N> &a,
             const array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        double aimbim = a[i + N / 2] * b[i + N / 2];
        double arebim = a[i] * b[i + N / 2];
        res[i] = std::fma(a[i], b[i], -aimbim);
        res[i + N / 2] = std::fma(a[i + N / 2], b[i], arebim);
    }
}

// removing inline seems to be faster in my environment.
template <uint32_t N>
void FMAInFD(array<double, N> &res, const array<double, N> &a,
             const array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        res[i] = std::fma(a[i + N / 2], b[i + N / 2], -res[i]);
        res[i] = std::fma(a[i], b[i], -res[i]);
        res[i + N / 2] = std::fma(a[i], b[i + N / 2], res[i + N / 2]);
        res[i + N / 2] = std::fma(a[i + N / 2], b[i], res[i + N / 2]);
    }
}

template <class P>
void PolynomialMulByXai(Polynomial<P> &res, const Polynomial<P> &poly,
                        const typename P::T a);

template <class P>
void PolynomialMulByXaiMinusOne(Polynomial<P> &res, const Polynomial<P> &poly,
                                const typename P::T a);

inline void PolynomialMulByXaiMinusOnelvl1(Polynomial<lvl1param> &res,
                                           const Polynomial<lvl1param> &poly,
                                           const typename lvl1param::T a)
{
    PolynomialMulByXaiMinusOne<lvl1param>(res, poly, a);
}

#define TFHEPP_EXPLICIT_INST_WRT_LVL0_1_2(fun) \
    fun(lvl0param);                            \
    fun(lvl1param);                            \
    fun(lvl2param);
#define TFHEPP_EXPLICIT_INST_WRT_LVL1_2(fun) \
    fun(lvl1param);                          \
    fun(lvl2param);
#define TFHEPP_EXPLICIT_INST_WRT_LVL01_02(fun) \
    fun(lvl01param);                           \
    fun(lvl02param);
#define TFHEPP_EXPLICIT_INST_WRT_LVL10_21_22(fun) \
    fun(lvl10param);                              \
    fun(lvl21param);                              \
    fun(lvl22param);
#define TFHEPP_EXPLICIT_INST_WRT_LVL0221_0222(fun) \
    fun(lvl02param, lvl21param);                   \
    fun(lvl02param, lvl22param);
}  // namespace TFHEpp