#include "detwfa.hpp"

namespace TFHEpp {
template <class P>
void CMUXFFT(TRLWE<P> &res, const TRGSWFFT<P> &cs, const TRLWE<P> &c1,
             const TRLWE<P> &c0)
{
    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) res[k][i] = c1[k][i] - c0[k][i];
    trgswfftExternalProduct<P>(res, res, cs);
    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) res[k][i] += c0[k][i];
}
#define INST(P)                                                     \
    template void CMUXFFT<P>(TRLWE<P> & res, const TRGSWFFT<P> &cs, \
                             const TRLWE<P> &c1, const TRLWE<P> &c0)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
void CMUXFFTwithPolynomialMulByXaiMinusOne(TRLWE<P> &acc, const TRGSWFFT<P> &cs,
                                           const typename P::T a)
{
    TRLWE<P> temp;
    for (int k = 0; k < P::k + 1; k++)
        PolynomialMulByXaiMinusOne<P>(temp[k], acc[k], a);
    trgswfftExternalProduct<P>(temp, temp, cs);
    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) acc[k][i] += temp[k][i];
}
#define INST(P)                                             \
    template void CMUXFFTwithPolynomialMulByXaiMinusOne<P>( \
        TRLWE<P> & acc, const TRGSWFFT<P> &cs, const typename P::T a)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
void CMUXNTTwithPolynomialMulByXaiMinusOne(TRLWE<P> &acc, const TRGSWNTT<P> &cs,
                                           const typename P::T a)
{
    TRLWE<P> temp;
    for (int k = 0; k < P::k + 1; k++)
        PolynomialMulByXaiMinusOne<P>(temp[k], acc[k], a);
    trgswnttExternalProduct<P>(temp, temp, cs);
    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) acc[k][i] += temp[k][i];
}
#define INST(P)                                             \
    template void CMUXNTTwithPolynomialMulByXaiMinusOne<P>( \
        TRLWE<P> & acc, const TRGSWNTT<P> &cs, const typename P::T a)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

}  // namespace TFHEpp
