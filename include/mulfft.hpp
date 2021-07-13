#pragma once

#include <spqlios-fft.h>

#include "cuhe++.hpp"
#include "params.hpp"
#include "utils.hpp"

namespace TFHEpp {

inline const std::array<std::array<cuHEpp::INTorus, TFHEpp::lvl1param::n>, 2>
    ntttwistlvl1 = cuHEpp::TwistGen<TFHEpp::lvl1param::nbit>();
inline const std::array<std::array<cuHEpp::INTorus, TFHEpp::lvl1param::n>, 2>
    ntttablelvl1 = cuHEpp::TableGen<TFHEpp::lvl1param::nbit>();

template <class P>
inline void TwistNTT(Polynomial<P> &res, const PolynomialNTT<P> &a)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>)
        cuHEpp::TwistNTTlvl1<typename TFHEpp::lvl1param::T,
                             TFHEpp::lvl1param::nbit>(res, a, ntttablelvl1[0],
                                                      ntttwistlvl1[0]);
    // else if constexpr (std::is_same_v<typename P::T, uint64_t>)
    //     fftplvl2.execute_direct_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistNTT!");
}

template <class P>
inline void TwistFFT(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>)
        fftplvl1.execute_direct_torus32(res.data(), a.data());
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        fftplvl2.execute_direct_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistFFT!");
}

template <class P>
inline void TwistINTT(PolynomialNTT<P> &res, const Polynomial<P> &a)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>)
        cuHEpp::TwistINTTlvl1<typename TFHEpp::lvl1param::T,
                              TFHEpp::lvl1param::nbit>(res, a, ntttablelvl1[1],
                                                       ntttwistlvl1[1]);
    // else if constexpr (std::is_same_v<typename P::T, uint64_t>)
    //     fftplvl2.execute_reverse_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistINTT!");
}

template <class P>
inline void TwistIFFT(PolynomialInFD<P> &res, const Polynomial<P> &a)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>)
        fftplvl1.execute_reverse_torus32(res.data(), a.data());
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        fftplvl2.execute_reverse_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistIFFT!");
}

template <class P>
inline void PolyMul(Polynomial<P> &res, const Polynomial<P> &a,
                    const Polynomial<P> &b)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>) {
        PolynomialInFD<P> ffta;
        TwistIFFT<P>(ffta, a);
        PolynomialInFD<P> fftb;
        TwistIFFT<P>(fftb, b);
        MulInFD<P::n>(ffta, ffta, fftb);
        TwistFFT<P>(res, ffta);
    }
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        for (int i = 0; i < P::n; i++) {
            typename P::T ri = 0;
            for (int j = 0; j <= i; j++)
                ri +=
                    static_cast<typename std::make_signed<typename P::T>::type>(
                        a[j]) *
                    b[i - j];
            for (int j = i + 1; j < P::n; j++)
                ri -=
                    static_cast<typename std::make_signed<typename P::T>::type>(
                        a[j]) *
                    b[P::n + i - j];
            res[i] = ri;
        }
    }
    else
        static_assert(false_v<typename P::T>, "Undefined PolyMul!");
}
}  // namespace TFHEpp