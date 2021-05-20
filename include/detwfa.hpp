#pragma once

#include "./trgsw.hpp"
#include "./utils.hpp"

namespace TFHEpp {
template <class P>
void CMUXFFT(TRLWE<P> &res, const TRGSWFFT<P> &cs, const TRLWE<P> &c1,
             const TRLWE<P> &c0);

template <class P>
void CMUXFFTwithPolynomialMulByXaiMinusOne(TRLWE<P> &acc, const TRGSWFFT<P> &cs,
                                           const typename P::T a);

}  // namespace TFHEpp