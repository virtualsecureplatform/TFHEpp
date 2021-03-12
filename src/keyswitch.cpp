#include <cloudkey.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {

void IdentityKeySwitchlvl10(TLWE<lvl0param> &res, const TLWE<lvl1param> &tlwe,
                            const KeySwitchingKey<lvl10param> &ksk)
{
    IdentityKeySwitch<lvl10param>(res, tlwe, ksk);
}
}  // namespace TFHEpp
