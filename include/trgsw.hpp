#pragma once

#include <array>
#include <cstdint>
#include <params.hpp>

namespace TFHEpp {
using namespace std;

TRGSWlvl1 trgswSymEncryptlvl1(int32_t p, double α, Keylvl1 &key);
TRGSWFFTlvl1 trgswfftSymEncryptlvl1(int32_t p, double α, Keylvl1 &key);
void trgswfftExternalProductlvl1(TRLWElvl1 &res, const TRLWElvl1 &trlwe,
                                 const TRGSWFFTlvl1 &trgswfft);

TRGSWlvl2 trgswSymEncryptlvl2(int64_t p, double α, Keylvl2 &key);
TRGSWFFTlvl2 trgswfftSymEncryptlvl2(int64_t p, double α, Keylvl2 &key);
void trgswfftExternalProductlvl2(TRLWElvl2 &res, const TRLWElvl2 &trlwe,
                                 const TRGSWFFTlvl2 &trgswfft);
}  // namespace TFHEpp