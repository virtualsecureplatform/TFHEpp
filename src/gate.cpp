#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {
void HomNAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] - cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomNOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] - cb[i];
    res[DEF_n] -= 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomXNOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -2 * (ca[i] + cb[i]);
    res[DEF_n] -= 1U << 30;
    GateBootstrapping(res, res, ck);
}

void HomAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i] + cb[i];
    res[DEF_n] -= 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
           const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i] + cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomXOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = 2 * (ca[i] + cb[i]);
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomNOT(TLWElvl0 &res, const TLWElvl0 &ca)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i];
}

void HomANDNY(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
              const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] + cb[i];
    res[DEF_n] -= 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomANDYN(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
              const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i] - cb[i];
    res[DEF_n] -= 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomORNY(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] + cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomORYN(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i] - cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, ck);
}

void HomMUX(TLWElvl0 &res, const TLWElvl0 &cs, const TLWElvl0 &c1,
            const TLWElvl0 &c0, const CloudKey &ck)
{
    TLWElvl0 temp;
    for (int i = 0; i <= DEF_n; i++) temp[i] = cs[i] + c1[i];
    for (int i = 0; i <= DEF_n; i++) res[i] = -cs[i] + c0[i];
    temp[DEF_n] -= 1U << 29;
    res[DEF_n] += 1U << 29;
    TLWElvl1 and1;
    TLWElvl1 and0;
    GateBootstrappingTLWE2TLWEFFTlvl01(and1, temp, ck);
    GateBootstrappingTLWE2TLWEFFTlvl01(and0, res, ck);

    for (int i = 0; i <= DEF_n; i++) and1[i] += and0[i];
    and1[DEF_N] += 1U << 29;
    IdentityKeySwitchlvl10(res, and1, ck.ksk);
}

}  // namespace TFHEpp