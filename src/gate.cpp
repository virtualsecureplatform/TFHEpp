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

// cs?c1:c0
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

void HomMUXwoSE(TRLWElvl1 &res, const TLWElvl0 &cs, const TLWElvl0 &c1,
                const TLWElvl0 &c0, const GateKey &gk)
{
    TLWElvl0 temp1;
    TLWElvl0 temp0;
    for (int i = 0; i <= DEF_n; i++) temp1[i] = cs[i] + c1[i];
    for (int i = 0; i <= DEF_n; i++) temp0[i] = -cs[i] + c0[i];
    temp1[DEF_n] -= 1U << 29;
    temp0[DEF_n] -= 1U << 29;
    TRLWElvl1 and0;
    GateBootstrappingTLWE2TRLWEFFTlvl01(res, temp1, gk);
    GateBootstrappingTLWE2TRLWEFFTlvl01(and0, temp0, gk);

    for (int i = 0; i < DEF_N; i++) {
        res[0][i] += and0[0][i];
        res[1][i] += and0[1][i];
    };
    res[1][0] += 1U << 29;
}

void ExtractSwitchAndHomMUX(TRLWElvl1 &res, const TRLWElvl1 &csr,
                            const TRLWElvl1 &c1r, const TRLWElvl1 &c0r,
                            const GateKey &gk)
{
    TLWElvl1 templvl1;
    TLWElvl0 cs, c1, c0;
    SampleExtractIndexlvl1(templvl1, csr, 0);
    IdentityKeySwitchlvl10(cs, templvl1, gk.ksk);
    SampleExtractIndexlvl1(templvl1, c1r, 0);
    IdentityKeySwitchlvl10(c1, templvl1, gk.ksk);
    SampleExtractIndexlvl1(templvl1, c0r, 0);
    IdentityKeySwitchlvl10(c0, templvl1, gk.ksk);

    for (int i = 0; i <= DEF_n; i++) c1[i] += cs[i];
    for (int i = 0; i <= DEF_n; i++) c0[i] -= cs[i];
    c1[DEF_n] -= 1U << 29;
    c0[DEF_n] -= 1U << 29;
    TRLWElvl1 and0;
    GateBootstrappingTLWE2TRLWEFFTlvl01(res, c1, gk);
    GateBootstrappingTLWE2TRLWEFFTlvl01(and0, c0, gk);

    for (int i = 0; i < DEF_N; i++) {
        res[0][i] += and0[0][i];
        res[1][i] += and0[1][i];
    };
    res[1][0] += 1U << 29;
}

}  // namespace TFHEpp