#pragma once

#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;
void HomNAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
}  // namespace TFHEpp