#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {

// No input
void HomCONSTANTONE(TLWE<lvl1param> &res)
{
    res = {};
    res[lvl1param::n] = lvl1param::μ;
}

void HomCONSTANTZERO(TLWE<lvl1param> &res)
{
    res = {};
    res[lvl1param::n] = -lvl1param::μ;
}

// 1 input
void HomNOT(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca)
{
    for (int i = 0; i <= lvl1param::n; i++) res[i] = -ca[i];
}

void HomCOPY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca)
{
    for (int i = 0; i <= lvl1param::n; i++) res[i] = ca[i];
}

template <int casign, int cbsign, typename lvl1param::T offset>
inline void HomGate(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                    const TLWE<lvl1param> &cb, const GateKey &gk)
{
    for (int i = 0; i <= lvl1param::n; i++)
        res[i] = casign * ca[i] + cbsign * cb[i];
    res[lvl1param::n] += offset;
    GateBootstrapping(res, res, gk);
}

void HomNAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<-1, -1, lvl1param::μ>(res, ca, cb, gk);
}

void HomNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<-1, -1, -lvl1param::μ>(res, ca, cb, gk);
}

void HomXNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<-2, -2, -2 * lvl1param::μ>(res, ca, cb, gk);
}

void HomAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<1, 1, -lvl1param::μ>(res, ca, cb, gk);
}

void HomOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
           const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<1, 1, lvl1param::μ>(res, ca, cb, gk);
}

void HomXOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<2, 2, 2 * lvl1param::μ>(res, ca, cb, gk);
}

void HomANDNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<-1, 1, -lvl1param::μ>(res, ca, cb, gk);
}

void HomANDYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<1, -1, -lvl1param::μ>(res, ca, cb, gk);
}

void HomORNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<-1, 1, lvl1param::μ>(res, ca, cb, gk);
}

void HomORYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const GateKey &gk)
{
    HomGate<1, -1, lvl1param::μ>(res, ca, cb, gk);
}

// 3input
// cs?c1:c0
void HomMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
            const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
            const GateKey &gk)
{
    TLWE<lvl1param> temp;
    for (int i = 0; i <= lvl1param::n; i++) temp[i] = cs[i] + c1[i];
    for (int i = 0; i <= lvl1param::n; i++) res[i] = -cs[i] + c0[i];
    temp[lvl1param::n] -= lvl1param::μ;
    res[lvl1param::n] -= lvl1param::μ;
    TLWE<lvl0param> and1, and0;
    IdentityKeySwitch<lvl10param>(and1, temp, gk.ksk);
    IdentityKeySwitch<lvl10param>(and0, res, gk.ksk);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        temp, and1, gk.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, and0, gk.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());

    for (int i = 0; i <= lvl1param::n; i++) res[i] += temp[i];
    res[lvl1param::n] += lvl1param::μ;
}

void HomNMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
             const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
             const GateKey &gk)
{
    TLWE<lvl1param> temp;
    for (int i = 0; i <= lvl1param::n; i++) temp[i] = cs[i] + c1[i];
    for (int i = 0; i <= lvl1param::n; i++) res[i] = -cs[i] + c0[i];
    temp[lvl1param::n] -= lvl1param::μ;
    res[lvl1param::n] -= lvl1param::μ;
    TLWE<lvl0param> and1, and0;
    IdentityKeySwitch<lvl10param>(and1, temp, gk.ksk);
    IdentityKeySwitch<lvl10param>(and0, res, gk.ksk);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        temp, and1, gk.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, and0, gk.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());

    for (int i = 0; i <= lvl1param::n; i++) res[i] = -res[i] - temp[i];
    res[lvl1param::n] -= lvl1param::μ;
}

template <class iksP, class bkP>
void HomMUXwoSE(TRLWE<typename bkP::targetP> &res,
                const TLWE<typename iksP::domainP> &cs,
                const TLWE<typename iksP::domainP> &c1,
                const TLWE<typename iksP::domainP> &c0,
                const KeySwitchingKey<iksP> &iksk,
                const BootstrappingKeyFFT<bkP> &bkfft)
{
    TLWE<typename iksP::domainP> temp1;
    TLWE<typename iksP::domainP> temp0;
    for (int i = 0; i <= iksP::domainP::n; i++) temp1[i] = cs[i] + c1[i];
    for (int i = 0; i <= iksP::domainP::n; i++) temp0[i] = -cs[i] + c0[i];
    temp1[iksP::domainP::n] -= iksP::domainP::μ;
    temp0[iksP::domainP::n] -= iksP::domainP::μ;
    TLWE<lvl0param> and1, and0;
    IdentityKeySwitch<iksP>(and1, temp1, iksk);
    IdentityKeySwitch<iksP>(and0, temp0, iksk);
    TRLWE<typename bkP::targetP> and0trlwe;
    BlindRotate<bkP>(res, and1, bkfft,
                     μpolygen<typename bkP::targetP, bkP::targetP::μ>());
    BlindRotate<bkP>(and0trlwe, and0, bkfft,
                     μpolygen<typename bkP::targetP, bkP::targetP::μ>());

    for (int i = 0; i < bkP::targetP::n; i++) {
        res[0][i] += and0trlwe[0][i];
        res[1][i] += and0trlwe[1][i];
    };
    res[1][0] += bkP::targetP::μ;
}
#define INST(iksP, bkP)                         \
    template void HomMUXwoSE<iksP, bkP>(        \
        TRLWE<typename bkP::targetP> & res,     \
        const TLWE<typename iksP::domainP> &cs, \
        const TLWE<typename iksP::domainP> &c1, \
        const TLWE<typename iksP::domainP> &c0, \
        const KeySwitchingKey<iksP> &iksk,      \
        const BootstrappingKeyFFT<bkP> &bkfft)
TFHEPP_EXPLICIT_INSTANTIATION_GATE(INST)
#undef INST

void ExtractSwitchAndHomMUX(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &csr,
                            const TRLWE<lvl1param> &c1r,
                            const TRLWE<lvl1param> &c0r, const GateKey &gk)
{
    TLWE<lvl1param> templvl1;
    TLWE<lvl0param> cs, c1, c0;
    SampleExtractIndex<lvl1param>(templvl1, csr, 0);
    IdentityKeySwitch<lvl10param>(cs, templvl1, gk.ksk);
    SampleExtractIndex<lvl1param>(templvl1, c1r, 0);
    IdentityKeySwitch<lvl10param>(c1, templvl1, gk.ksk);
    SampleExtractIndex<lvl1param>(templvl1, c0r, 0);
    IdentityKeySwitch<lvl10param>(c0, templvl1, gk.ksk);

    for (int i = 0; i <= lvl0param::n; i++) c1[i] += cs[i];
    for (int i = 0; i <= lvl0param::n; i++) c0[i] -= cs[i];
    c1[lvl0param::n] -= lvl0param::μ;
    c0[lvl0param::n] -= lvl0param::μ;
    TRLWE<lvl1param> and0;
    BlindRotate<lvl01param>(res, c1, gk.bkfftlvl01,
                            μpolygen<lvl1param, lvl1param::μ>());
    BlindRotate<lvl01param>(and0, c0, gk.bkfftlvl01,
                            μpolygen<lvl1param, lvl1param::μ>());

    for (int i = 0; i < lvl1param::n; i++) {
        res[0][i] += and0[0][i];
        res[1][i] += and0[1][i];
    };
    res[1][0] += lvl1param::n;
}

}  // namespace TFHEpp