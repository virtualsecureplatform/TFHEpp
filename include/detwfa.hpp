#pragma once

#include "trgsw.hpp"

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

template <class P>
TRGSWFFT<P> TRGSWFFTOneGen()
{
    constexpr std::array<typename P::T, P::l> h = hgen<P>();

    TRGSW<P> trgsw = {};
    for (int i = 0; i < P::l; i++) {
        for (int k = 0; k < P::k + 1; k++) {
            trgsw[i + k * P::l][k][0] += h[i];
        }
    }
    return ApplyFFT2trgsw<P>(trgsw);
}

alignas(64) const TRGSWFFT<lvl1param> trgswonelvl1 =
    TRGSWFFTOneGen<lvl1param>();
alignas(64) const TRGSWFFT<lvl2param> trgswonelvl2 =
    TRGSWFFTOneGen<lvl2param>();

template <class bkP>
void CMUXFFTwithPolynomialMulByXaiMinusOne(
    TRLWE<typename bkP::targetP> &acc,
    const BootstrappingKeyElementFFT<bkP> &cs, const int a)
{
    if constexpr (bkP::domainP::key_value_diff == 1) {
        alignas(64) TRLWE<typename bkP::targetP> temp;
        for (int k = 0; k < bkP::targetP::k + 1; k++)
            PolynomialMulByXaiMinusOne<typename bkP::targetP>(temp[k], acc[k],
                                                              a);
        trgswfftExternalProduct<typename bkP::targetP>(temp, temp, cs[0]);
        for (int k = 0; k < bkP::targetP::k + 1; k++)
            for (int i = 0; i < bkP::targetP::n; i++) acc[k][i] += temp[k][i];
    }
    else {
#ifdef USE_TERNARY_CMUX
        alignas(32) TRGSWFFT<typename bkP::targetP> trgsw;
        if constexpr (std::is_same_v<typename bkP::targetP, lvl1param>) {
            trgsw = trgswonelvl1;
        }
        else if constexpr (std::is_same_v<typename bkP::targetP, lvl2param>) {
            trgsw = trgswonelvl2;
        }
        int count = 0;
        alignas(32) PolynomialInFD<typename bkP::targetP> poly;
        for (int i = bkP::domainP::key_value_min;
             i <= bkP::domainP::key_value_max; i++) {
            if (i != 0) {
                for (int j = 0; j < (bkP::targetP::k + 1) * bkP::targetP::l;
                     j++) {
                    for (int k = 0; k < bkP::targetP::k + 1; k++) {
                        PolynomialMulByXaiMinusOneInFD<typename bkP::targetP>(
                            poly, cs[count][j][k], a * i);
                        for (int l = 0; l < bkP::targetP::n; l++)
                            trgsw[j][k][l] += poly[l];
                    }
                }
                count++;
            }
        }
        trgswfftExternalProduct<typename bkP::targetP>(acc, acc, trgsw);
#else
        alignas(32) TRLWE<typename bkP::targetP> temp;
        int count = 0;
        for (int i = bkP::domainP::key_value_min;
             i <= bkP::domainP::key_value_max; i++) {
            if (i != 0) {
                const int mod = (a * i) % (2 * bkP::targetP::n);
                const int index = mod > 0 ? mod : mod + (2 * bkP::targetP::n);
                for (int k = 0; k < bkP::targetP::k + 1; k++)
                    PolynomialMulByXaiMinusOne<typename bkP::targetP>(
                        temp[k], acc[k], index);
                trgswfftExternalProduct<typename bkP::targetP>(temp, temp,
                                                               cs[count]);
                for (int k = 0; k < bkP::targetP::k + 1; k++)
                    for (int i = 0; i < bkP::targetP::n; i++)
                        acc[k][i] += temp[k][i];
                count++;
            }
        }
#endif
    }
}

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

template <class P>
void CMUXRAINTTwithPolynomialMulByXaiMinusOne(TRLWE<P> &acc,
                                              const TRGSWRAINTT<P> &cs,
                                              const typename P::T a)
{
    TRLWE<P> temp;
    for (int k = 0; k < P::k + 1; k++)
        PolynomialMulByXaiMinusOne<P>(temp[k], acc[k], a);
    trgswrainttExternalProduct<P>(temp, temp, cs);
    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) acc[k][i] += temp[k][i];
}

}  // namespace TFHEpp