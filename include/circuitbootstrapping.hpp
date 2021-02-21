#pragma once

#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;

void CircuitBootstrappingFFT(TRGSWFFTlvl1 &trgswfft, const TLWElvl0 &tlwe,
                             const CircuitKey &ck);
void CircuitBootstrappingFFTInv(TRGSWFFTlvl1 &invtrgswfft, const TLWElvl0 &tlwe,
                                const CircuitKey &ck);
void CircuitBootstrappingFFTwithInv(TRGSWFFTlvl1 &trgswfft,
                                    TRGSWFFTlvl1 &invtrgswfft,
                                    const TLWElvl0 &tlwe, const CircuitKey &ck);
}  // namespace TFHEpp