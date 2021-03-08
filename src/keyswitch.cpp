#include <cloudkey.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {

void IdentityKeySwitchlvl10(TLWE<lvl0param> &res, const TLWE<lvl1param> &tlwe,
                            const KeySwitchingKey<lvl10param> &ksk)
{
    IdentityKeySwitch<lvl10param>(res, tlwe, ksk);
}

void PrivKeySwitchlvl21(TRLWE<lvl1param> &res, const TLWE<lvl2param> &tlwe,
                        const int u,
                        const PrivKeySwitchKey<lvl21param> &privksk)
{
    PrivKeySwitch<lvl21param>(res, tlwe, u, privksk);
}

void PrivKeySwitchlvl22(TRLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe,
                        const int u,
                        const PrivKeySwitchKey<lvl22param> &privksk)
{
    PrivKeySwitch<lvl22param>(res, tlwe, u, privksk);
}
}  // namespace TFHEpp
