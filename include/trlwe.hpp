#pragma once

#include <array>

namespace TFHEpp {
using namespace std;
TRLWElvl1 trlweSymEncryptZerolvl1(const double α, const Keylvl1 &key);
TRLWElvl1 trlweSymEncryptlvl1(const Polynomiallvl1 &p, const double α,
                              const Keylvl1 &key);
array<bool, DEF_N> trlweSymDecryptlvl1(const TRLWElvl1 &c, const Keylvl1 &key);
void SampleExtractIndexlvl1(TLWElvl1 &tlwe, const TRLWElvl1 &trlwe,
                            const int index);
}  // namespace TFHEpp