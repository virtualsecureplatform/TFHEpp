#pragma once

#include <array>

namespace TFHEpp {
using namespace std;
TRLWElvl1 trlweSymEncryptZerolvl1(const double α, const Keylvl1 &key);
TRLWElvl1 trlweSymEncryptlvl1(const Polynomiallvl1 &p, const double α,
                              const Keylvl1 &key);
TRLWElvl2 trlweSymEncryptZerolvl2(const double α, const Keylvl2 &key);
TRLWElvl2 trlweSymEncryptlvl2(const Polynomiallvl2 &p, const double α,
                              const Keylvl2 &key);
array<bool, DEF_N> trlweSymDecryptlvl1(const TRLWElvl1 &c, const Keylvl1 &key);
array<bool, DEF_nbar> trlweSymDecryptlvl2(const TRLWElvl2 &c, const Keylvl2 &key);
void SampleExtractIndexlvl1(TLWElvl1 &tlwe, const TRLWElvl1 &trlwe,
                            const int index);
void SampleExtractIndexlvl2(TLWElvl2 &tlwe, const TRLWElvl2 &trlwe,
                            const int index);
}  // namespace TFHEpp