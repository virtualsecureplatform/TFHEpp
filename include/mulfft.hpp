#pragma once

#include <spqlios-fft.h>

#include <params.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;

template<class P>
inline void TwistFFT(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr(P::T==uint32_t) fftplvl1.execute_direct_torus32(res.data(), a.data());
    else if constexpr(P::T==uint64_t) fftplvl2.execute_direct_torus64(res.data(), a.data());
    else static_assert(false_v<T>, "Undefined TwistFFT!");
}

template<class P>
inline void TwistFFT(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr(P::T==uint32_t) fftplvl1.execute_reverse_torus32(res.data(), a.data());
    else if constexpr(P::T==uint64_t) fftplvl2.execute_direct_torus64(res.data(), a.data());
    else static_assert(false_v<T>, "Undefined TwistIFFT!");
}

template<class P>
inline void PolyMul(Polynomial<P> &res, const Polynomial<P> &a,
                        const Polynomial<P> &b)
{
    if constexpr(P::T==uint32_t){
        PolynomialInFD<P> ffta;
        TwistIFFT<P>(ffta, a);
        PolynomialInFD<P> fftb;
        TwistIFFT<P>(fftb, b);
        MulInFD<P::n>(ffta, ffta, fftb);
        TwistFFT<P>(res, ffta);
    }else if constexpr(P::T==uint64_t){
        for (int i = 0; i < P::n; i++) {
            P::T ri = 0;
            for (int j = 0; j <= i; j++)
                ri += static_cast<make_signed<P::T>>(a[j]) * b[i - j];
            for (int j = i + 1; j < DEF_nbar; j++)
                ri -= static_cast<make_signed<P::T>>(a[j]) * b[P::n + i - j];
            res[i] = ri;
        }
    }
    else static_assert(false_v<T>, "Undefined PolyMul!");
}
}  // namespace TFHEpp