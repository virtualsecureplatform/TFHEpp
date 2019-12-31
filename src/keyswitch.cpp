#include <cloudkey.hpp>

namespace TFHEpp {
void IdentityKeySwitchlvl10(
    array<uint32_t, DEF_n + 1> &res, array<uint32_t, DEF_N + 1> &tlwe,
    array<
        array<array<array<uint32_t, DEF_n + 1>, (1 << DEF_basebit) - 1>, DEF_t>,
        DEF_N> &ksk)
{
    const uint32_t prec_offset = 1 << (32 - (1 + DEF_basebit * DEF_t));
    const uint32_t mask = (1U << DEF_basebit) - 1;
    res = {};
    res[DEF_n] = tlwe[DEF_N];
    for (int i = 0; i < DEF_n; i++) {
        const uint32_t aibar = tlwe[i] + prec_offset;
        for (int j = 0; j < DEF_t; j++) {
            const uint32_t aij = (aibar >> (32 - (j + 1) * DEF_basebit)) & mask;
            if (aibar != 0)
                for (int k = 0; k <= DEF_n; k++) res[k] -= ksk[i][j][aij][k];
        }
    }
}
}  // namespace TFHEpp
