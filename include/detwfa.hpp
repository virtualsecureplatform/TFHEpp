#pragma once

#include <trgsw.hpp>

namespace TFHEpp {
void CMUXFFTlvl1(TRLWE<lvl1param> &res, const TRGSWFFT<lvl1param> &cs, const TRLWE<lvl1param> &c1,
                 const TRLWE<lvl1param> &c0);
}