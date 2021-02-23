#pragma once

#include <array>
#include <cstdint>
#include <params.hpp>

namespace TFHEpp {
using namespace std;

TRGSWFFT<lvl1param> trgswfftSymEncryptlvl1(make_signed<lvl1param::T> p, double α, Key<lvl1param> &key);
void trgswfftExternalProductlvl1(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &trlwe,
                                 const TRGSWFFT<lvl1param> &trgswfft);

TRGSWFFT<lvl2param> trgswfftSymEncryptlvl2(make_signed<lvl2param::T> p, double α, Key<lvl2param> &key);
void trgswfftExternalProductlvl2(TRLWE<lvl2param> &res, const TRLWE<lvl2param> &trlwe,
                                 const TRGSWFFT<lvl2param> &trgswfft);
}  // namespace TFHEpp