#pragma once

#include <array>
#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;

void GateBootstrappingTLWE2TLWEFFTlvl01(TLWElvl1 &res,
                                               const TLWElvl0 &tlwe,
                                               const CloudKey &ck);
void GateBootstrapping(TLWElvl0 &res, const TLWElvl0 &tlwe, const CloudKey &ck);
void CircuitBootstrappingFFT(TRGSWFFTlvl1 &trgswfft, const TLWElvl0 &tlwe, const CloudKey &ck);
}  // namespace TFHEpp