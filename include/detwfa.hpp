#pragma once

#include "trgsw.hpp"

namespace TFHEpp {
template <class P>
void CMUXFFT(TRLWE<P> &res, const TRGSWFFT<P> &cs, const TRLWE<P> &c1,
             const TRLWE<P> &c0);

template <class bkP>
void CMUXFFTwithPolynomialMulByXaiMinusOne(
    TRLWE<typename bkP::targetP> &acc,
    const BootstrappingKeyElementFFT<bkP> &cs, const int a);
template <class P>
void CMUXNTTwithPolynomialMulByXaiMinusOne(TRLWE<P> &acc, const TRGSWNTT<P> &cs,
                                           const typename P::T a);

}  // namespace TFHEpp