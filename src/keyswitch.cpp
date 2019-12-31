#include <cloudkey.hpp>

namespace TFHEpp {
void IdentityKeySwitchlvl10(TLWElvl0 &res, TLWElvl1 &tlwe,
                            const KeySwitchingKey &ksk)
{
    const uint32_t prec_offset = 1 << (32 - (1 + DEF_basebit * DEF_t));
    const uint32_t mask = (1U << DEF_basebit) - 1;
    res = {};
    res[DEF_n] = tlwe[DEF_N];
    for (int i = 0; i < DEF_N; i++) {
        const uint32_t aibar = tlwe[i] + prec_offset;
        for (int j = 0; j < DEF_t; j++) {
            const uint32_t aij = (aibar >> (32 - (j + 1) * DEF_basebit)) & mask;
            if (aij != 0)
                for (int k = 0; k <= DEF_n; k++) res[k] -= ksk[i][j][aij-1][k];
        }
    }
}
}  // namespace TFHEpp
