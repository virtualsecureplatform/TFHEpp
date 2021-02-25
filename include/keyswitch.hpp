#pragma once

#include <array>
#include <params.hpp>

namespace TFHEpp {
using namespace std;

template<class P>
inline void IdentityKeySwitch(TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
                            const KeySwitchingKey<P> &ksk)
{
    const typename P::domainP::T prec_offset = 1ULL << (numeric_limits<typename P::domainP::T>::digits - (1 + P::basebit * P::t));
    const uint32_t mask = (1U << P::basebit) - 1;
    res = {};
    res[P::targetP::n] = tlwe[P::domainP::n];
    for (int i = 0; i < P::domainP::n; i++) {
        const uint32_t aibar = tlwe[i] + prec_offset;
        for (int j = 0; j <P::t; j++) {
            const uint32_t aij = (aibar >> (numeric_limits<typename P::domainP::T>::digits - (j + 1) * P::basebit)) & mask;
            if (aij != 0)
                for (int k = 0; k <= P::targetP::n; k++)
                    res[k] -= ksk[i][j][aij - 1][k];
        }
    }
}
void IdentityKeySwitchlvl10(TLWE<lvl0param> &res, const TLWE<lvl1param> &tlwe,
                            const KeySwitchingKey<lvl10param> &ksk);
void PrivKeySwitchlvl21(TRLWE<lvl1param> &res, const TLWE<lvl2param> &tlwe, const int u,
                        const PrivKeySwitchKey<lvl21param> &privksk);
void PrivKeySwitchlvl22(TRLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe, const int u,
                        const PrivKeySwitchKey<lvl22param> &privksk);
}  // namespace TFHEpp