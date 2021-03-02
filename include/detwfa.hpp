#pragma once

#include <trgsw.hpp>
#include <utils.hpp>

namespace TFHEpp {
void CMUXFFTlvl1(TRLWE<lvl1param> &res, const TRGSWFFT<lvl1param> &cs,
                 const TRLWE<lvl1param> &c1, const TRLWE<lvl1param> &c0);

template<class P>
inline void CMUXFFTwithPolynomialMulByXaiMinusOne(TRLWE<P> &acc, const TRGSWFFT<P> &cs, const typename P::T a)
{
    TRLWE<P> temp;
    PolynomialMulByXaiMinusOne<P>(temp[0],acc[0],a);
    PolynomialMulByXaiMinusOne<P>(temp[1],acc[1],a);
    trgswfftExternalProduct<P>(temp, temp, cs);
    for (int i = 0; i < P::n; i++) {
        acc[0][i] += temp[0][i];
        acc[1][i] += temp[1][i];
    }
}

void CMUXFFTwithPolynomialMulByXaiMinusOnelvl1(TRLWE<lvl1param> &acc, const TRGSWFFT<lvl1param> &cs, const typename lvl1param::T a);

}