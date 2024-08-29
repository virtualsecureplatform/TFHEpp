#pragma once

#include <memory>

#include "INTorus.hpp"
#ifdef USE_FFTW3
#include <fft_processor_fftw.h>
#elif USE_MKL
#include <fft_processor_mkl.hpp>
#elif USE_SPQLIOX_AARCH64
#include <fft_processor_spqliox_aarch64.h>
#elif USE_CONCRETE_FFT
#include <fft_processor_concrete.hpp>
#else
#include <fft_processor_spqlios.h>
#endif
#ifdef USE_HEXL
#include "hexl/hexl.hpp"
#endif

#ifdef USE_INTERLEAVED_FORMAT
#include <complex>
#endif

#include "cuhe++.hpp"
#include "params.hpp"
#include "utils.hpp"

namespace TFHEpp {

inline const std::unique_ptr<
    const std::array<std::array<cuHEpp::INTorus, lvl1param::n>, 2>>
    ntttwistlvl1 = cuHEpp::TwistGen<lvl1param::nbit>();
inline const std::unique_ptr<
    const std::array<std::array<cuHEpp::INTorus, lvl1param::n>, 2>>
    ntttablelvl1 = cuHEpp::TableGen<lvl1param::nbit>();
inline const std::unique_ptr<
    const std::array<std::array<cuHEpp::INTorus, lvl2param::n>, 2>>
    ntttwistlvl2 = cuHEpp::TwistGen<lvl2param::nbit>();
inline const std::unique_ptr<
    const std::array<std::array<cuHEpp::INTorus, lvl2param::n>, 2>>
    ntttablelvl2 = cuHEpp::TableGen<lvl2param::nbit>();
inline const std::unique_ptr<
    std::array<std::array<raintt::SWord, lvl1param::n>, 2>>
    raintttwist = raintt::TwistGen<lvl1param::nbit, 3>();
inline const std::unique_ptr<
    std::array<std::array<std::array<raintt::SWord, lvl1param::n>, 2>, 2>>
    raintttable = raintt::TableGen<lvl1param::nbit>();
#ifdef USE_HEXL
// Biggest prime number less than 2^30 and staisfies 1 mod 2N.
constexpr uint64_t lvl1P = 1073707009;
#endif

template <class P>
inline void TwistNTT(Polynomial<P> &res, PolynomialNTT<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>)
#ifdef USE_HEXL
    {
        std::array<uint64_t, lvl1param::n> temp;
        static intel::hexl::NTT nttlvl1(lvl1param::n, lvl1P);
        nttlvl1.ComputeInverse(temp.data(), &(a[0].value), 1, 1);
        for (int i = 0; i < lvl1param::n; i++) res[i] = (temp[i] << 32) / lvl1P;
    }
#else
        cuHEpp::TwistNTT<typename lvl1param::T, lvl1param::nbit>(
            res, a, (*ntttablelvl1)[0], (*ntttwistlvl1)[0]);
#endif
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        cuHEpp::TwistNTT<typename lvl2param::T, lvl2param::nbit>(
            res, a, (*ntttablelvl2)[0], (*ntttwistlvl2)[0]);
    }
    else
        static_assert(false_v<typename P::T>, "Undefined TwistNTT!");
}

template <class P>
inline void TwistFFT(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>) {
        if constexpr (std::is_same_v<typename P::T, uint32_t>)
            // if constexpr(hasq<P>)
            // fftplvl1.execute_direct_torus32_q(res.data(), a.data(), P::q);
            // else
            fftplvl1.execute_direct_torus32(res.data(), a.data());
        if constexpr (std::is_same_v<typename P::T, uint64_t>)
            fftplvl1.execute_direct_torus64(res.data(), a.data());
    }
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        fftplvl2.execute_direct_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistFFT!");
}

template <class P>
inline void TwistFFTrescale(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>) {
        if constexpr (std::is_same_v<typename P::T, uint32_t>)
            fftplvl1.execute_direct_torus32_rescale(res.data(), a.data(), P::Δ);
        else if constexpr (std::is_same_v<typename P::T, uint64_t>)
            fftplvl1.execute_direct_torus64_rescale(res.data(), a.data(), P::Δ);
    }
    else if constexpr (std::is_same_v<P, lvl2param>)
        fftplvl2.execute_direct_torus64_rescale(res.data(), a.data(), P::Δ);
    else
        static_assert(false_v<typename P::T>, "Undefined TwistFFT!");
}

template <class P>
inline void TwistINTT(PolynomialNTT<P> &res, const Polynomial<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>)
#ifdef USE_HEXL
    {
        std::array<uint64_t, lvl1param::n> temp;
        for (int i = 0; i < lvl1param::n; i++)
            temp[i] = (lvl1P * static_cast<uint64_t>(a[i])) >> 32;
        static intel::hexl::NTT nttlvl1(lvl1param::n, lvl1P);
        nttlvl1.ComputeForward(&(res[0].value), temp.data(), 1, 1);
    }
#else
        cuHEpp::TwistINTT<typename P::T, P::nbit>(res, a, (*ntttablelvl1)[1],
                                                  (*ntttwistlvl1)[1]);
#endif
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        cuHEpp::TwistINTT<typename P::T, P::nbit>(res, a, (*ntttablelvl2)[1],
                                                  (*ntttwistlvl2)[1]);
    else
        static_assert(false_v<typename P::T>, "Undefined TwistINTT!");
}

template <class P>
inline void TwistIFFT(PolynomialInFD<P> &res, const Polynomial<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>) {
        if constexpr (std::is_same_v<typename P::T, uint32_t>)
            fftplvl1.execute_reverse_torus32(res.data(), a.data());
        if constexpr (std::is_same_v<typename P::T, uint64_t>)
            fftplvl1.execute_reverse_torus64(res.data(), a.data());
    }
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        fftplvl2.execute_reverse_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistIFFT!");
}

template <class P>
inline void TwistIFFTUInt(PolynomialInFD<P> &res, const Polynomial<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>) {
        if constexpr (std::is_same_v<typename P::T, uint32_t>)
            fftplvl1.execute_reverse_uint(res.data(), a.data());
        // if constexpr (std::is_same_v<typename P::T, uint64_t>)
        // fftplvl1.execute_reverse_torus64(res.data(), a.data());
    }
    // else if constexpr (std::is_same_v<typename P::T, uint64_t>)
    // fftplvl2.execute_reverse_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistIFFT!");
}

template <uint32_t N>
inline void MulInFD(std::array<double, N> &res, const std::array<double, N> &b)
{
#ifdef USE_INTERLEAVED_FORMAT
    for (int i = 0; i < N / 2; i++) {
        const std::complex tmp = std::complex(res[2 * i], res[2 * i + 1]) *
                                 std::complex(b[2 * i], b[2 * i + 1]);
        res[2 * i] = tmp.real();
        res[2 * i + 1] = tmp.imag();
    }
#else
    for (int i = 0; i < N / 2; i++) {
        double aimbim = res[i + N / 2] * b[i + N / 2];
        double arebim = res[i] * b[i + N / 2];
        res[i] = std::fma(res[i], b[i], -aimbim);
        res[i + N / 2] = std::fma(res[i + N / 2], b[i], arebim);
    }
#endif
}

template <uint32_t N>
inline void MulInFD(std::array<double, N> &res, const std::array<double, N> &a,
                    const std::array<double, N> &b)
{
#ifdef USE_INTERLEAVED_FORMAT
    for (int i = 0; i < N / 2; i++) {
        const std::complex tmp = std::complex(a[2 * i], a[2 * i + 1]) *
                                 std::complex(b[2 * i], b[2 * i + 1]);
        res[2 * i] = tmp.real();
        res[2 * i + 1] = tmp.imag();
    }
#else
    // for (int i = 0; i < N / 2; i++) {
    //     double aimbim = a[i + N / 2] * b[i + N / 2];
    //     double arebim = a[i] * b[i + N / 2];
    //     res[i] = std::fma(a[i], b[i], -aimbim);
    //     res[i + N / 2] = std::fma(a[i + N / 2], b[i], arebim);
    // }

    // for (int i = 0; i < N / 2; i++) {
    //     res[i] = a[i + N / 2] * b[i + N / 2];
    //     res[i + N / 2] = a[i] * b[i + N / 2];
    // }
    // for (int i = 0; i < N / 2; i++) {
    //     res[i] = std::fma(a[i], b[i], -res[i]);
    //     res[i + N / 2] = std::fma(a[i + N / 2], b[i], res[i + N / 2]);
    // }

    // for (int i = 0; i < N / 2; i++) {
    //     double arebre = a[i] * b[i];
    //     double aimbre = a[i + N/2] * b[i];
    //     res[i] = std::fma(- a[i + N / 2] , b[i + N / 2],arebre);
    //     res[i + N / 2] = std::fma(a[i], b[i  + N / 2], aimbre);
    // }

    for (int i = 0; i < N / 2; i++) {
        res[i] = a[i] * b[i];
        res[i + N / 2] = a[i + N / 2] * b[i];
    }
    for (int i = 0; i < N / 2; i++) {
        res[i + N / 2] += a[i] * b[i + N / 2];
        res[i] -= a[i + N / 2] * b[i + N / 2];
    }
#endif
}

// Be careful about memory accesss (We assume b has relatively high memory
// access cost)
template <uint32_t N>
inline void FMAInFD(std::array<double, N> &res, const std::array<double, N> &a,
                    const std::array<double, N> &b)
{
#ifdef USE_INTERLEAVED_FORMAT
    for (int i = 0; i < N / 2; i++) {
        std::complex tmp = std::complex(a[2 * i], a[2 * i + 1]) *
                           std::complex(b[2 * i], b[2 * i + 1]);
        res[2 * i] += tmp.real();
        res[2 * i + 1] += tmp.imag();
    }
#else
    for (int i = 0; i < N / 2; i++) {
        res[i] = std::fma(a[i], b[i], res[i]);
        res[i + N / 2] = std::fma(a[i + N / 2], b[i], res[i + N / 2]);
    }
    for (int i = 0; i < N / 2; i++) {
        res[i + N / 2] = std::fma(a[i], b[i + N / 2], res[i + N / 2]);
        res[i] -= a[i + N / 2] * b[i + N / 2];
    }
// for (int i = 0; i < N / 2; i++) {
//     res[i] = std::fma(a[i + N / 2], b[i + N / 2], -res[i]);
//     res[i] = std::fma(a[i], b[i], -res[i]);
//     res[i + N / 2] = std::fma(a[i], b[i + N / 2], res[i + N / 2]);
//     res[i + N / 2] = std::fma(a[i + N / 2], b[i], res[i + N / 2]);
// }
#endif
}

template <class P>
inline void PolyMul(Polynomial<P> &res, const Polynomial<P> &a,
                    const Polynomial<P> &b)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>) {
        alignas(64) PolynomialInFD<P> ffta;
        TwistIFFT<P>(ffta, a);
        alignas(64) PolynomialInFD<P> fftb;
        TwistIFFT<P>(fftb, b);
        MulInFD<P::n>(ffta, fftb);
        TwistFFT<P>(res, ffta);
        // alignas(64) PolynomialInFD<P> fftres;
        // MulInFD<P::n>(fftres, ffta, fftb);
        // TwistFFT<P>(res, fftres);
    }
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        // Naieve
        // for (int i = 0; i < P::n; i++) {
        //     typename P::T ri = 0;
        //     for (int j = 0; j <= i; j++)
        //         ri +=
        //             static_cast<typename std::make_signed<typename
        //             P::T>::type>(
        //                 a[j]) *
        //             b[i - j];
        //     for (int j = i + 1; j < P::n; j++)
        //         ri -=
        //             static_cast<typename std::make_signed<typename
        //             P::T>::type>(
        //                 a[j]) *
        //             b[P::n + i - j];
        //     res[i] = ri;
        // }
        PolynomialNTT<P> ntta, nttb;
        TwistINTT<P>(ntta, a);
        TwistINTT<P>(nttb, b);
        for (int i = 0; i < P::n; i++) ntta[i] *= nttb[i];
        TwistNTT<P>(res, ntta);
    }
    else {
        // Naieve
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
}

template <class P>
inline void PolyMulRescaleUnsigned(Polynomial<P> &res, const Polynomial<P> &a,
                                   const Polynomial<P> &b)
{
    // if constexpr (std::is_same_v<typename P::T, uint32_t>) {
    PolynomialInFD<P> ffta, fftb;
    TwistIFFTUInt<P>(ffta, a);
    TwistIFFTUInt<P>(fftb, b);
    MulInFD<P::n>(ffta, fftb);
    TwistFFTrescale<P>(res, ffta);
    // }
    // else
    //     static_assert(false_v<typename P::T>, "Undefined PolyMul!");
}

template <class P>
inline void PolyMulNaive(Polynomial<P> &res, const Polynomial<P> &a,
                         const Polynomial<P> &b)
{
    for (int i = 0; i < P::n; i++) {
        typename P::T ri = 0;
        for (int j = 0; j <= i; j++)
            ri += static_cast<typename std::make_signed<typename P::T>::type>(
                      a[j]) *
                  b[i - j];
        for (int j = i + 1; j < P::n; j++)
            ri -= static_cast<typename std::make_signed<typename P::T>::type>(
                      a[j]) *
                  b[P::n + i - j];
        res[i] = ri;
    }
}

template <class P>
inline void PolyMulNaieveRescaleUnsigned(Polynomial<P> &res,
                                         const Polynomial<P> &a,
                                         const Polynomial<P> &b)
{
    for (int i = 0; i < P::n; i++) {
        __int128_t ri = 0;
        for (int j = 0; j <= i; j++)
            ri += static_cast<__int128_t>(a[j]) *
                  static_cast<__int128_t>(b[i - j]);
        for (int j = i + 1; j < P::n; j++)
            ri -= static_cast<__int128_t>(a[j]) *
                  static_cast<__int128_t>(b[P::n + i - j]);
        // res[i] = static_cast<typename P::T>((ri) >>
        // (std::numeric_limits<typename P::T>::digits - 3));
        res[i] = static_cast<typename P::T>((ri) >> 29);
    }
}

template <class P>
std::unique_ptr<std::array<PolynomialInFD<P>, 2 * P::n>> XaittGen()
{
    std::unique_ptr<std::array<PolynomialInFD<P>, 2 * P::n>> xaitt =
        std::make_unique<std::array<PolynomialInFD<P>, 2 * P::n>>();
    for (int i = 0; i < 2 * P::n; i++) {
        std::array<typename P::T, P::n> xai = {};
        xai[0] = -1;
        if (i < P::n)
            xai[i] += 1;
        else
            xai[i - P::n] -= 1;
        TwistIFFT<P>((*xaitt)[i], xai);
    }
    return xaitt;
}

template <class P>
std::unique_ptr<std::array<PolynomialNTT<P>, 2 * P::n>> XaittGenNTT()
{
    std::unique_ptr<std::array<PolynomialNTT<P>, 2 * P::n>> xaitt =
        std::make_unique<std::array<PolynomialNTT<P>, 2 * P::n>>();
    for (int i = 0; i < 2 * P::n; i++) {
        std::array<typename P::T, P::n> xai = {};
        xai[0] = -1;
        if (i < P::n)
            xai[i] += 1;
        else
            xai[i - P::n] -= 1;
        TwistINTT<P>((*xaitt[i]), xai);
    }
    return xaitt;
}

#if defined(USE_TERNARY) || defined(USE_KEY_BUNDLE)
alignas(64) static const std::unique_ptr<
    const std::array<PolynomialInFD<lvl1param>, 2 * lvl1param::n>> xaittlvl1 =
    XaittGen<lvl1param>();
alignas(64) static const std::unique_ptr<
    const std::array<PolynomialInFD<lvl2param>, 2 * lvl2param::n>> xaittlvl2 =
    XaittGen<lvl2param>();
#endif
#ifdef USE_TERNARY
template <class P>
inline void PolynomialMulByXaiMinusOneInFD(PolynomialInFD<P> &res,
                                           const PolynomialInFD<P> &poly,
                                           const int a)
{
    const int mod = a % (2 * P::n);
    const int index = mod > 0 ? mod : mod + (2 * P::n);
    if constexpr (std::is_same_v<P, lvl1param>) {
        MulInFD<P::n>(res, poly, (*xaittlvl1)[index]);
    }
}
#endif
}  // namespace TFHEpp