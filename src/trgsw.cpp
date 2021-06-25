#include "trgsw.hpp"

#include <limits>
#include <type_traits>

namespace TFHEpp {
using namespace std;

template <class P>
constexpr typename P::T offsetgen()
{
    typename P::T offset = 0;
    for (int i = 1; i <= P::l; i++)
        offset +=
            P::Bg / 2 *
            (1ULL << (numeric_limits<typename P::T>::digits - i * P::Bgbit));
    return offset;
}

template <class P>
void DecompositionPolynomial(DecomposedPolynomial<P> &decpoly,
                             const Polynomial<P> &poly, const int digit)
{
    constexpr typename P::T offset = offsetgen<P>();
    constexpr typename P::T roundoffset =
        1ULL << (std::numeric_limits<typename P::T>::digits - P::l * P::Bgbit -
                 1);
    constexpr typename P::T mask =
        static_cast<typename P::T>((1ULL << P::Bgbit) - 1);
    constexpr typename P::T halfBg = (1ULL << (P::Bgbit - 1));

    for (int i = 0; i < P::n; i++) {
        decpoly[i] = (((poly[i] + offset + roundoffset) >>
                       (numeric_limits<typename P::T>::digits -
                        (digit + 1) * P::Bgbit)) &
                      mask) -
                     halfBg;
    }
}
#define INST(P)                                                       \
    template void DecompositionPolynomial<P>(                         \
        DecomposedPolynomial<P> & decpoly, const Polynomial<P> &poly, \
        const int digit)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
void DecompositionPolynomialFFT(DecomposedPolynomialInFD<P> &decpolyfft,
                                const Polynomial<P> &poly, const int digit)
{
    DecomposedPolynomial<P> decpoly;
    DecompositionPolynomial<P>(decpoly, poly, digit);
    TwistIFFT<P>(decpolyfft, decpoly);
}
#define INST(P)                                                              \
    template void DecompositionPolynomialFFT<P>(                             \
        DecomposedPolynomialInFD<P> & decpolyfft, const Polynomial<P> &poly, \
        const int digit)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
void trgswfftExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                             const TRGSWFFT<P> &trgswfft)
{
    DecomposedPolynomialInFD<P> decpolyfft;
    DecompositionPolynomialFFT<P>(decpolyfft, trlwe[0], 0);
    TRLWEInFD<P> restrlwefft;
    MulInFD<P::n>(restrlwefft[0], decpolyfft, trgswfft[0][0]);
    MulInFD<P::n>(restrlwefft[1], decpolyfft, trgswfft[0][1]);
    for (int i = 1; i < P::l; i++) {
        DecompositionPolynomialFFT<P>(decpolyfft, trlwe[0], i);
        FMAInFD<P::n>(restrlwefft[0], decpolyfft, trgswfft[i][0]);
        FMAInFD<P::n>(restrlwefft[1], decpolyfft, trgswfft[i][1]);
    }
    for (int i = 0; i < P::l; i++) {
        DecompositionPolynomialFFT<P>(decpolyfft, trlwe[1], i);
        FMAInFD<P::n>(restrlwefft[0], decpolyfft, trgswfft[i + P::l][0]);
        FMAInFD<P::n>(restrlwefft[1], decpolyfft, trgswfft[i + P::l][1]);
    }
    TwistFFT<P>(res[0], restrlwefft[0]);
    TwistFFT<P>(res[1], restrlwefft[1]);
}
#define INST(P)                               \
    template void trgswfftExternalProduct<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWFFT<P> &trgswfft)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
TRGSWFFT<P> ApplyFFT2trgsw(const TRGSW<P> &trgsw)
{
    TRGSWFFT<P> trgswfft;
    for (int i = 0; i < 2 * P::l; i++)
        for (int j = 0; j < 2; j++) TwistIFFT<P>(trgswfft[i][j], trgsw[i][j]);
    return trgswfft;
}
#define INST(P) template TRGSWFFT<P> ApplyFFT2trgsw<P>(const TRGSW<P> &trgsw)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
TRGSWNTT<P> ApplyNTT2trgsw(const TRGSW<P> &trgsw)
{
    TRGSWNTT<P> trgswntt;
    for (int i = 0; i < 2 * P::l; i++)
        for (int j = 0; j < 2; j++) TwistINTT<P>(trgswntt[i][j], trgsw[i][j]);
    return trgswntt;
}
#define INST(P) template TRGSWNTT<P> ApplyNTT2trgsw<P>(const TRGSW<P> &trgsw)
// TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
INST(lvl1param);
#undef INST

template <class P>
TRGSWNTT<P> TRGSW2NTT(const TRGSW<P> &trgsw)
{
    TRGSWNTT<P> trgswntt;
    for (int i = 0; i < 2 * P::l; i++)
        for (int j = 0; j < 2; j++) {
            PolynomialNTT<P> temp;
            TwistINTT<P>(temp, trgsw[i][j]);
            for (uint32_t k = 0; k < P::n; k++)
                trgswntt[i][j][k] = temp[cuHEpp::BitReverse<P::nbit>(k)];
        }
    return trgswntt;
}
#define INST(P) template TRGSWNTT<P> TRGSW2NTT<P>(const TRGSW<P> &trgsw)
// TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
INST(lvl1param);
#undef INST

template <class P>
TRGSW<P> trgswSymEncrypt(const typename make_signed<typename P::T>::type p,
                         const double α, const Key<P> &key)
{
    constexpr array<typename P::T, P::l> h = hgen<P>();

    TRGSW<P> trgsw;
    for (TRLWE<P> &trlwe : trgsw) trlwe = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::l; i++) {
        trgsw[i][0][0] += static_cast<typename P::T>(p) * h[i];
        trgsw[i + P::l][1][0] += static_cast<typename P::T>(p) * h[i];
    }
    return trgsw;
}
#define INST(P)                                                            \
    template TRGSW<P> trgswSymEncrypt<P>(                                  \
        const typename make_signed<typename P::T>::type p, const double α, \
        const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
TRGSWFFT<P> trgswfftSymEncrypt(
    const typename make_signed<typename P::T>::type p, const double α,
    const Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, α, key);
    return ApplyFFT2trgsw<P>(trgsw);
}
#define INST(P)                                                            \
    template TRGSWFFT<P> trgswfftSymEncrypt<P>(                            \
        const typename make_signed<typename P::T>::type p, const double α, \
        const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST
}  // namespace TFHEpp
