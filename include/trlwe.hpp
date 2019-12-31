#pragma once

#include <array>

namespace TFHEpp {
using namespace std;
array<array<uint32_t, DEF_N>, 2> trlweSymEncryptZerolvl1(
    const double α, const array<uint32_t, DEF_N> &key);
array<array<uint32_t, DEF_N>, 2> trlweSymEncryptlvl1(
    const array<uint32_t, DEF_N> &p, const double α,
    const array<uint32_t, DEF_N> &key);
array<bool, DEF_N> trlweSymDecryptlvl1(
    const array<array<uint32_t, DEF_N>, 2> &c,
    const array<uint32_t, DEF_N> &key);
void SampleExtractIndexlvl1(array<uint32_t, DEF_N + 1> &tlwe,
                            const array<array<uint32_t, DEF_N>, 2> &trlwe,
                            const int index);
}  // namespace TFHEpp