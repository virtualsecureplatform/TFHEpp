#pragma once

#include <array>
#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;

void GateBootstrapping(array<uint32_t, DEF_n + 1> &res,
                       const array<uint32_t, DEF_n + 1> &tlwe, CloudKey &ck);
}  // namespace TFHEpp