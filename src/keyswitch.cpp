#include <cloudkey.hpp>

namespace TFHEpp {
void IdentityKeySwitchlvl10(TLWElvl0 &res, const TLWElvl1 &tlwe,
                            const KeySwitchingKey &ksk)
{
    constexpr uint32_t prec_offset = 1U << (32 - (1 + DEF_basebit * DEF_t));
    constexpr uint32_t mask = (1U << DEF_basebit) - 1;
    res = {};
    res[DEF_n] = tlwe[DEF_N];
    for (int i = 0; i < DEF_N; i++) {
        const uint32_t aibar = tlwe[i] + prec_offset;
        for (int j = 0; j < DEF_t; j++) {
            const uint32_t aij = (aibar >> (32 - (j + 1) * DEF_basebit)) & mask;
            if (aij != 0)
                for (int k = 0; k <= DEF_n; k++)
                    res[k] -= ksk[i][j][aij - 1][k];
        }
    }
}

void PrivKeySwitchlvl21(TRLWElvl1 &res, const TLWElvl2 &tlwe, const int u,
                        const PrivKeySwitchKey &privksk)
{
    const uint32_t mask = (1 << DEF_basebitlvl21) - 1;
    const uint64_t prec_offset = 1UL
                                 << (64 - (1 + DEF_basebitlvl21 * DEF_tbar));

    res = {};
    for (int i = 0; i <= DEF_nbar; i++) {
        uint64_t aibar = tlwe[i] + prec_offset;

        for (int j = 0; j < DEF_tbar; j++) {
            const uint64_t aij =
                (aibar >> (64 - (j + 1) * DEF_basebitlvl21)) & mask;

            if (aij != 0) {
                for (int p = 0; p < DEF_N; p++) {
                    res[0][p] -= privksk[u][i][j][aij - 1][0][p];
                    res[1][p] -= privksk[u][i][j][aij - 1][1][p];
                }
            }
        }
    }
}
}  // namespace TFHEpp
