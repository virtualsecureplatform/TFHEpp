#include <cloudkey.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {

void IdentityKeySwitchlvl10(TLWE<lvl0param> &res, const TLWE<lvl1param> &tlwe,
                            const KeySwitchingKey<lvl10param> &ksk)
{
    IdentityKeySwitch<lvl10param>(res,tlwe,ksk);
}

void PrivKeySwitchlvl21(TRLWElvl1 &res, const TLWElvl2 &tlwe, const int u,
                        const PrivKeySwitchKey &privksk)
{
    constexpr uint32_t mask = (1 << DEF_basebitlvl21) - 1;
    constexpr uint64_t prec_offset =
        1ULL << (64 - (1 + DEF_basebitlvl21 * DEF_tbar));

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

void PrivKeySwitchlvl22(TRLWElvl2 &res, const TLWElvl2 &tlwe, const int u,
                        const PrivKeySwitchlvl22Key &privksk)
{
    constexpr uint32_t mask = (1 << DEF_basebitlvl22) - 1;
    constexpr uint64_t prec_offset =
        1ULL << (64 - (1 + DEF_basebitlvl22 * DEF_tlvl22));

    res = {};
    for (int i = 0; i <= DEF_nbar; i++) {
        uint64_t aibar = tlwe[i] + prec_offset;

        for (int j = 0; j < DEF_tlvl22; j++) {
            const uint64_t aij =
                (aibar >> (64 - (j + 1) * DEF_basebitlvl22)) & mask;

            if (aij != 0) {
                for (int p = 0; p < DEF_nbar; p++) {
                    res[0][p] -= privksk[u][i][j][aij - 1][0][p];
                    res[1][p] -= privksk[u][i][j][aij - 1][1][p];
                }
            }
        }
    }
}
}  // namespace TFHEpp
