#pragma once

#ifdef USE_RANDEN
#include <randen.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <random>

#include "params.hpp"

namespace TFHEpp {
#ifdef USE_RANDEN
static thread_local std::random_device trng;
static thread_local randen::Randen<uint64_t> generator(trng());
#else
static thread_local std::random_device generator;
#endif

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
        std::normal_distribution<double> distribution(0., stdev);
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

// https://stackoverflow.com/questions/21191307/minimum-number-of-bits-to-represent-a-given-int
template <uint32_t data>
constexpr int bits_needed()
{
    uint32_t value = data;
    int bits = 0;
    for (int bit_test = 16; bit_test > 0; bit_test >>= 1) {
        if (value >> bit_test != 0) {
            bits += bit_test;
            value >>= bit_test;
        }
    }
    return bits + value;
}

template <uint32_t N>
inline void MulInFD(std::array<double, N> &res, const std::array<double, N> &a,
                    const std::array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        double aimbim = a[i + N / 2] * b[i + N / 2];
        double arebim = a[i] * b[i + N / 2];
        res[i] = std::fma(a[i], b[i], -aimbim);
        res[i + N / 2] = std::fma(a[i + N / 2], b[i], arebim);
    }
}

// Be careful about memory accesss (We assume b has relatively high memory access cost)
template <uint32_t N>
inline void FMAInFD(std::array<double, N> &res, const std::array<double, N> &a,
             const std::array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        res[i] = std::fma(a[i], b[i], res[i]);
        res[i + N / 2] = std::fma(a[i + N / 2], b[i], res[i + N / 2]);
    }
    for (int i = 0; i < N / 2; i++) {
        res[i + N / 2] = std::fma(a[i], b[i + N / 2], res[i + N / 2]);
        res[i] -= a[i + N / 2] * b[i + N / 2];
    }
}

template <class P>
inline void PolynomialMulByXai(Polynomial<P> &res, const Polynomial<P> &poly,
                               const typename P::T a)
{
    if (a == 0)
        res = poly;
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

// calcurate τ_d
template <class P>
inline void Automorphism(Polynomial<P> &res, const Polynomial<P> &poly,
                         const uint d)
{
    res = {};
    constexpr uint Nmask = (1ULL << (P::nbit)) - 1;
    constexpr uint signmask = 1ULL << (P::nbit);
    for (uint i = 0; i < P::n; i++) {
        const uint index = i * d;
        if (index & signmask)
            res[index & Nmask] -= poly[i];
        else
            res[index & Nmask] += poly[i];
    }
}
}  // namespace TFHEpp
