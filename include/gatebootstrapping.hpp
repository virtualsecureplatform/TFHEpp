#pragma once

#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;

void PolynomialMulByXaiMinusOnelvl1(Polynomial<lvl1param> &res,
                                    const Polynomial<lvl1param> &poly,
                                    const lvl1param::T a);
void GateBootstrappingTLWE2TRLWEFFTlvl01(TRLWE<lvl1param> &acc, const TLWE<lvl0param> &tlwe,
                                         const GateKey &gk);
void GateBootstrappingTLWE2TLWEFFTlvl01(TLWE<lvl1param> &res, const TLWE<lvl0param> &tlwe,
                                        const GateKey &gk);
void GateBootstrappingTLWE2TLWEFFTlvl02(
    TLWE<lvl2param> &res, const TLWE<lvl0param> &tlwe,
    const BootStrappingKeyFFT<lvl02param> &bkfftlvl02, const lvl2param::T μs2);
void GateBootstrapping(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe, const GateKey &gk);
}  // namespace TFHEpp