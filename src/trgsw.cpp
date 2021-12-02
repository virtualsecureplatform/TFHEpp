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
void DecompositionPolynomialNTT(DecomposedPolynomialNTT<P> &decpolyntt,
                                const Polynomial<P> &poly, const int digit)
{
    DecomposedPolynomial<P> decpoly;
    DecompositionPolynomial<P>(decpoly, poly, digit);
    TwistINTT<P>(decpolyntt, decpoly);
}
#define INST(P)                                                             \
    template void DecompositionPolynomialNTT<P>(                            \
        DecomposedPolynomialNTT<P> & decpolyntt, const Polynomial<P> &poly, \
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
    for (int m = 0; m < P::k+1; m++)
        MulInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[0][m]);
    for (int i = 1; i < P::l; i++) {
        DecompositionPolynomialFFT<P>(decpolyfft, trlwe[0], i);
        for (int m = 0; m < P::k+1; m++)
            FMAInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[i][m]);
    }
    for (int k = 1; k < P::k+1; k++){
        for (int i = 0; i < P::l; i++) {
            DecompositionPolynomialFFT<P>(decpolyfft, trlwe[k], i);
            for (int m = 0; m < P::k+1; m++)
                FMAInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[i + k*P::l][m]);
        }
    }
    for (int k = 0; k < P::k+1; k++)
        TwistFFT<P>(res[k], restrlwefft[k]);
}
#define INST(P)                               \
    template void trgswfftExternalProduct<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWFFT<P> &trgswfft)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
void trgswnttExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                             const TRGSWNTT<P> &trgswntt)
{
    DecomposedPolynomialNTT<P> decpolyntt;
    DecompositionPolynomialNTT<P>(decpolyntt, trlwe[0], 0);
    TRLWENTT<P> restrlwentt;
    for (int i = 0; i < P::n; i++)
        restrlwentt[0][i] = decpolyntt[i] * trgswntt[0][0][i];
    for (int i = 0; i < P::n; i++)
        restrlwentt[1][i] = decpolyntt[i] * trgswntt[0][1][i];
    for (int i = 1; i < P::l; i++) {
        DecompositionPolynomialNTT<P>(decpolyntt, trlwe[0], i);
        for (int m = 0; m < P::k+1; m++)
                for (int j = 0; j < P::n; j++)
                    restrlwentt[m][j] += decpolyntt[j] * trgswntt[i][m][j];
    }
    for (int k = 1; k < P::k+1; k++)
        for (int i = 0; i < P::l; i++) {
            DecompositionPolynomialNTT<P>(decpolyntt, trlwe[k], i);
            for (int m = 0; m < P::k+1; m++)
                for (int j = 0; j < P::n; j++)
                    restrlwentt[m][j] += decpolyntt[j] * trgswntt[i + k*P::l][m][j];
        }
    for (int k = 0; k < P::k+1; k++)
        TwistNTT<P>(res[k], restrlwentt[k]);
}
#define INST(P)                               \
    template void trgswnttExternalProduct<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWNTT<P> &trgswntt)
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
    for (int i = 0; i < (P::k+1) * P::l; i++)
        for (int j = 0; j < P::k+1; j++) TwistINTT<P>(trgswntt[i][j], trgsw[i][j]);
    return trgswntt;
}
#define INST(P) template TRGSWNTT<P> ApplyNTT2trgsw<P>(const TRGSW<P> &trgsw)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
TRGSWNTT<P> TRGSW2NTT(const TRGSW<P> &trgsw)
{
    TRGSWNTT<P> trgswntt;
    for (int i = 0; i < (P::k+1) * P::l; i++)
        for (int j = 0; j < P::k+1; j++) {
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
TRGSW<P> trgswSymEncrypt(const Polynomial<P> &p, const double α,
                         const Key<P> &key)
{
    constexpr std::array<typename P::T, P::l> h = hgen<P>();

    TRGSW<P> trgsw;
    for (TRLWE<P> &trlwe : trgsw) trlwe = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::l; i++) {
        for (int k = 0; k < P::k+1; k++){
            for (int j = 0; j < P::n; j++) {
                trgsw[i + k*P::l][k][j] += static_cast<typename P::T>(p[j]) * h[i];
            }
        }
    }
    return trgsw;
}
#define INST(P)                                                  \
    template TRGSW<P> trgswSymEncrypt<P>(const Polynomial<P> &p, \
                                         const double α, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
TRGSWFFT<P> trgswfftSymEncrypt(const Polynomial<P> &p, const double α,
                               const Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, α, key);
    return ApplyFFT2trgsw<P>(trgsw);
}
#define INST(P)                                 \
    template TRGSWFFT<P> trgswfftSymEncrypt<P>( \
        const Polynomial<P> &p, const double α, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
TRGSWNTT<P> trgswnttSymEncrypt(const Polynomial<P> &p, const double α,
                               const Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, α, key);
    return ApplyNTT2trgsw<P>(trgsw);
}
#define INST(P)                                 \
    template TRGSWNTT<P> trgswnttSymEncrypt<P>( \
        const Polynomial<P> &p, const double α, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST
}  // namespace TFHEpp
