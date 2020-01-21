#pragma once

#include <trgsw.hpp>

namespace TFHEpp {
void CMUXFFTlvl1(TRLWElvl1 &res, const TRGSWFFTlvl1 &cs, const TRLWElvl1 &c1,
                 const TRLWElvl1 &c0);
}