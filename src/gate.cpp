#include <gatebootstrapping.hpp>

namespace TFHEpp {
void HomNAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] - cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, ck);
}
}  // namespace TFHEpp