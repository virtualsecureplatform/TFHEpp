#pragma once

#include <array>
#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;

void GateBootstrapping(TLWElvl0 &res, const TLWElvl0 &tlwe, CloudKey &ck);
}  // namespace TFHEpp