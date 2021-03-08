#pragma once

#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;

void CircuitBootstrappingFFTlvl01(TRGSWFFT<lvl1param> &trgswfft,
                                  const TLWE<lvl0param> &tlwe,
                                  const CircuitKey<lvl02param, lvl21param> &ck);
void CircuitBootstrappingFFTlvl02(TRGSWFFT<lvl2param> &trgswfft,
                                  const TLWE<lvl0param> &tlwe,
                                  const CircuitKey<lvl02param, lvl22param> &ck);
void CircuitBootstrappingFFTInvlvl01(
    TRGSWFFT<lvl1param> &invtrgswfft, const TLWE<lvl0param> &tlwe,
    const CircuitKey<lvl02param, lvl21param> &ck);
void CircuitBootstrappingFFTwithInvlvl01(
    TRGSWFFT<lvl1param> &trgswfft, TRGSWFFT<lvl1param> &invtrgswfft,
    const TLWE<lvl0param> &tlwe, const CircuitKey<lvl02param, lvl21param> &ck);
}  // namespace TFHEpp