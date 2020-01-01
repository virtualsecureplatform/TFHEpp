#pragma once

#include <array>
#include <params.hpp>

namespace TFHEpp {
using namespace std;

void IdentityKeySwitchlvl10(TLWElvl0 &res, TLWElvl1 &tlwe,
                            const KeySwitchingKey &ksk);
void PrivKeySwitchlvl21(TRLWElvl1 &res, const TLWElvl2 &tlwe, const int u,
                        const PrivKeySwitchKey &privksk);
}  // namespace TFHEpp