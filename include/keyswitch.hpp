#pragma once

#include <array>
#include <params.hpp>

namespace TFHEpp {
using namespace std;

void IdentityKeySwitchlvl10(TLWElvl0 &res, TLWElvl1 &tlwe,
                            KeySwitchingKey &ksk);
}  // namespace TFHEpp