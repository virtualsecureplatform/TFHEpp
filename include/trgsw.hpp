#pragma once

#include <array>
#include <cstdint>
#include <params.hpp>

namespace TFHEpp {
using namespace std;

TRGSWFFTlvl1 trgswfftSymEncryptlvl1(int32_t p, double α, Keylvl1 &key);
void trgswfftExternalProductlvl1(TRLWElvl1 &trlwe,
                                 const TRGSWFFTlvl1 &trgswfft);
TRGSWFFTlvl2 trgswfftSymEncryptlvl2(int64_t p, double α, Keylvl2 &key);
void trgswfftExternalProductlvl2(TRLWElvl2 &trlwe,
                                 const TRGSWFFTlvl2 &trgswfft);
}  // namespace TFHEpp