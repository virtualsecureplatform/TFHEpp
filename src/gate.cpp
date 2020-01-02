#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {
void HomNAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] - cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomNOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] - cb[i];
    res[DEF_n] -= 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomXNOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -2 * (ca[i] + cb[i]);
    res[DEF_n] -= 1U << 30;
    GateBootstrapping(res, res, gk);
}

void HomAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i] + cb[i];
    res[DEF_n] -= 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
           const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i] + cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomXOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = 2 * (ca[i] + cb[i]);
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomNOT(TLWElvl0 &res, const TLWElvl0 &ca)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i];
}

void HomCOPY(TLWElvl0 &res, const TLWElvl0 &ca)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i];
}

void HomCONSTANTONE(TLWElvl0 &res)
{
    res = {};
    res[DEF_n] = DEF_μ;
}

void HomCONSTANTZERO(TLWElvl0 &res)
{
    res = {};
    res[DEF_n] = -DEF_μ;
}

void HomANDNY(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
              const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] + cb[i];
    res[DEF_n] -= 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomANDYN(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
              const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i] - cb[i];
    res[DEF_n] -= 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomORNY(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = -ca[i] + cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomORYN(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const GateKey &gk)
{
    for (int i = 0; i <= DEF_n; i++) res[i] = ca[i] - cb[i];
    res[DEF_n] += 1U << 29;
    GateBootstrapping(res, res, gk);
}

void HomMUX(TLWElvl0 &res, const TLWElvl0 &cs, const TLWElvl0 &c1,
            const TLWElvl0 &c0, const GateKey &gk)
{
    TLWElvl0 temp;
    for (int i = 0; i <= DEF_n; i++) temp[i] = cs[i] + c1[i];
    for (int i = 0; i <= DEF_n; i++) res[i] = -cs[i] + c0[i];
    temp[DEF_n] -= 1U << 29;
    res[DEF_n] -= 1U << 29;
    TLWElvl1 and1;
    TLWElvl1 and0;
    GateBootstrappingTLWE2TLWEFFTlvl01(and1, temp, gk);
    GateBootstrappingTLWE2TLWEFFTlvl01(and0, res, gk);

    for (int i = 0; i <= DEF_N; i++) and1[i] += and0[i];
    and1[DEF_N] += 1U << 29;
    IdentityKeySwitchlvl10(res, and1, gk.ksk);
}

}  // namespace TFHEpp