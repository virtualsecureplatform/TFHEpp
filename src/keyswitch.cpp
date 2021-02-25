#include <cloudkey.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {

void IdentityKeySwitchlvl10(TLWE<lvl0param> &res, const TLWE<lvl1param> &tlwe,
                            const KeySwitchingKey<lvl10param> &ksk)
{
    IdentityKeySwitch<lvl10param>(res,tlwe,ksk);
}

template <class P>
void PrivKeySwitch(TRLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe, const int u,
                        const PrivKeySwitchKey<P> &privksk)
{
    constexpr uint32_t mask = (1 << P::basebit) - 1;
    constexpr uint64_t prec_offset =
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits - (1 + P::basebit * P::t));

    res = {};
    for (int i = 0; i <= P::domainP::n; i++) {
        typename P::domainP::T aibar = tlwe[i] + prec_offset;

        for (int j = 0; j < P::t; j++) {
            const typename P::domainP::T aij =
                (aibar >> (std::numeric_limits<typename P::domainP::T>::digits - (j + 1) * P::basebit)) & mask;

            if (aij != 0) {
                for (int p = 0; p < P::targetP::n; p++) {
                    res[0][p] -= privksk[u][i][j][aij - 1][0][p];
                    res[1][p] -= privksk[u][i][j][aij - 1][1][p];
                }
            }
        }
    }
}

void PrivKeySwitchlvl21(TRLWE<lvl1param> &res, const TLWE<lvl2param> &tlwe, const int u,
                        const PrivKeySwitchKey<lvl21param> &privksk)
{
    PrivKeySwitch<lvl21param>(res,tlwe,u,privksk);
}

void PrivKeySwitchlvl22(TRLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe, const int u,
                        const PrivKeySwitchKey<lvl22param> &privksk)
{
    PrivKeySwitch<lvl22param>(res,tlwe,u,privksk);
}
}  // namespace TFHEpp
