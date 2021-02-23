#pragma once

#include <array>

namespace TFHEpp {
using namespace std;
TRLWE<lvl1param> trlweSymEncryptZerolvl1(const double α, const Key<lvl1param> &key);
TRLWE<lvl1param> trlweSymEncryptlvl1(const Polynomial<lvl1param> &p, const double α,
                              const Key<lvl1param> &key);
TRLWE<lvl2param> trlweSymEncryptZerolvl2(const double α, const Key<lvl2param> &key);
TRLWE<lvl2param> trlweSymEncryptlvl2(const Polynomial<lvl2param> &p, const double α,
                              const Key<lvl2param> &key);
array<bool, lvl1param::n> trlweSymDecryptlvl1(const TRLWE<lvl1param> &c, const Key<lvl1param> &key);
array<bool, lvl2param::n> trlweSymDecryptlvl2(const TRLWE<lvl2param> &c,
                                          const Key<lvl2param> &key);
void SampleExtractIndexlvl1(TLWE<lvl1param> &tlwe, const TRLWE<lvl1param> &trlwe,
                            const int index);
void SampleExtractIndexlvl2(TLWE<lvl2param> &tlwe, const TRLWE<lvl2param> &trlwe,
                            const int index);
}  // namespace TFHEpp