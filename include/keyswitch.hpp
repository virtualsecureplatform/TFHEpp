#pragma once

#include <array>
#include <params.hpp>

namespace TFHEpp {
using namespace std;

void IdentityKeySwitchlvl10(
    array<uint32_t, DEF_n + 1> &res, array<uint32_t, DEF_N + 1> &tlwe,
    array<
        array<array<array<uint32_t, DEF_n + 1>, (1 << DEF_basebit) - 1>, DEF_t>,
        DEF_N> &ksk);
}  // namespace TFHEpp