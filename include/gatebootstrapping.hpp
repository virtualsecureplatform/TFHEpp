#pragma once

#include <array>
#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;
  
void PolynomialMulByXaiMinusOnelvl1(Polynomiallvl1 &res,
                                           const Polynomiallvl1 &poly,
                                           const uint32_t a);
void GateBootstrappingTLWE2TLWEFFTlvl01(TLWElvl1 &res, const TLWElvl0 &tlwe,
                                        const GateKey &gk);
void GateBootstrapping(TLWElvl0 &res, const TLWElvl0 &tlwe, const GateKey &gk);
void CircuitBootstrappingFFT(TRGSWFFTlvl1 &trgswfft, const TLWElvl0 &tlwe,
                             const CircuitKey &ck);
}  // namespace TFHEpp