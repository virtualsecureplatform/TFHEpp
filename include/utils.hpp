#pragma once

#include <randen.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <params.hpp>
#include <random>
#include <functional>

namespace TFHEpp {
static thread_local randen::Randen<uint64_t> generator;
// static thread_local std::random_device generator;

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
inline typename P::T ModularGaussian(typename P::T center, double stdev)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>) {
        // 32bit fixed-point number version
        std::normal_distribution<double> distribution(
            0., stdev);  // TODO: can we create a global distrib of param 1 and
                         // multiply by sigma?
        double err = distribution(generator);
        return center + dtot32(err);
    }
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        // 64bit fixed-point number version
        static const double _2p64 = std::pow(2., 64);
        std::normal_distribution<double> distribution(0., 1.0);
        const double val = stdev * distribution(generator) * _2p64;
        const uint64_t ival = static_cast<typename P::T>(val);
        return ival + center;
    }
    else
        static_assert(false_v<typename P::T>, "Undefined Modular Gaussian!");
}

template <class P>
inline typename P::T modSwitchFromTorus(uint32_t phase)
{
    constexpr uint32_t Mbit = P::nbit + 1;
    constexpr uint32_t Msize = 1U << Mbit;
    if constexpr (std::is_same_v<typename P::T, uint32_t>)
        return (phase + (1U << (31 - Mbit))) >> (32 - Mbit);
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        typename P::T interv =
            ((1ULL << 63) / Msize) * 2;  // width of each intervall
        typename P::T half_interval =
            interv / 2;  // begin of the first intervall

        // Mod Switching (as in modSwitchFromTorus32)
        typename P::T temp = (static_cast<typename P::T>(phase) << 32) +
                             half_interval;  // RIVEDI
        return temp / interv;
    }
    else
        static_assert(false_v<typename P::T>, "Undefined modSwitchFromTorus!");
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
inline void PolynomialMulByXai(Polynomial<P> &res,
                               const Polynomial<P> &poly,
                               const typename P::T a)
{
    if (a == 0)
        return;
    else if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a];
    }
    else {
        const typename P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa];
    }
}

template <class P>
inline void PolynomialMulByXaiMinusOne(Polynomial<P> &res,
                                       const Polynomial<P> &poly,
                                       const typename P::T a)
{
    if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n] - poly[i];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a] - poly[i];
    }
    else {
        const typename P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n] - poly[i];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa] - poly[i];
    }
}

static void PolynomialMulByXaiMinusOnelvl1(Polynomial<lvl1param> &res,
                                       const Polynomial<lvl1param> &poly,
                                       const typename lvl1param::T a){
                                           PolynomialMulByXaiMinusOne<lvl1param>(res,poly,a);
                                       }
}  // namespace TFHEpp