#pragma once

#include <array>
#include <cstdint>

#include "./mulfft.hpp"
#include "./params.hpp"
#include "./trlwe.hpp"

namespace TFHEpp {
using namespace std;

template <class P>
void DecompositionPolynomial(DecomposedPolynomial<P> &decpoly,
                             const Polynomial<P> &poly, const int digit);

template <class P>
void DecompositionPolynomialFFT(DecomposedPolynomialInFD<P> &decpolyfft,
                                const Polynomial<P> &poly, const int digit);

template <class P>
void trgswfftExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                             const TRGSWFFT<P> &trgswfft);

template <class P>
constexpr array<typename P::T, P::l> hgen()
{
    array<typename P::T, P::l> h{};
    for (int i = 0; i < P::l; i++)
        h[i] = 1ULL << (numeric_limits<typename P::T>::digits -
                        (i + 1) * P::Bgbit);
    return h;
}

template<class P>
inline TRGSWFFT<P> ApplyFFT2trgsw(const TRGSW<P> &trgsw){
    TRGSWFFT<P> trgswfft;
    for (int i = 0; i < 2 * P::l; i++)
        for (int j = 0; j < 2; j++) TwistIFFT<P>(trgswfft[i][j], trgsw[i][j]);
    return trgswfft;
}

template <class P>
TRGSWFFT<P> ApplyFFT2trgsw(const TRGSW<P> &trgsw);

template <class P>
TRGSW<P> trgswSymEncrypt(const typename make_signed<typename P::T>::type p,
                         const double α, const Key<P> &key);

template <class P>
TRGSWFFT<P> trgswfftSymEncrypt(
    const typename make_signed<typename P::T>::type p, const double α,
    const Key<P> &key);
}  // namespace TFHEpp