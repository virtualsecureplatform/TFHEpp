#pragma once

#include "cloudkey.hpp"
#include "keyswitch.hpp"
#include "gatebootstrapping.hpp"

namespace TFHEpp {
template <class P, int casign, int cbsign, uint64_t offset>
inline void HomGate(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
                    const EvalKey &ek)
{
    for (int i = 0; i <= P::k * P::n; i++)
        res[i] = casign * ca[i] + cbsign * cb[i];
    res[P::k * P::n] += offset;
    if constexpr (std::is_same_v<P, lvl1param>) GateBootstrapping<lvl10param,lvl01param,lvl1param::μ>(res, res, ek);
    else if constexpr (std::is_same_v<P, lvl0param>) GateBootstrapping<lvl01param,lvl1param::μ,lvl10param>(res, res, ek);
}

// No input
template <class P = lvl1param>
void HomCONSTANTONE(TLWE<P> &res)
{
    res = {};
    res[P::k * P::n] = P::μ;
}

template <class P = lvl1param>
void HomCONSTANTZERO(TLWE<P> &res)
{
    res = {};
    res[P::k * P::n] = -P::μ;
}

// 1 input
template <class P = lvl1param>
void HomNOT(TLWE<P> &res, const TLWE<P> &ca)
{
    for (int i = 0; i <= P::k * P::n; i++) res[i] = -ca[i];
}
template <class P = lvl1param>
void HomCOPY(TLWE<P> &res, const TLWE<P> &ca)
{
    for (int i = 0; i <= P::k * P::n; i++) res[i] = ca[i];
}
template <class P = lvl1param>
void HomNAND(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek)
{
    HomGate<P, -1, -1, lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomNOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek)
{
    HomGate<P, -1, -1, -lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomXNOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek)
{
    HomGate<P, -2, -2, -2 * lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomAND(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek)
{
    HomGate<P, 1, 1, -lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
           const EvalKey &ek)
{
    HomGate<P, 1, 1, lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomXOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek)
{
    HomGate<P, 2, 2, 2 * lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomANDNY(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
              const EvalKey &ek)
{
    HomGate<P, -1, 1, -lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomANDYN(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
              const EvalKey &ek)
{
    HomGate<P, 1, -1, -lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomORNY(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek)
{
    HomGate<P, -1, 1, lvl1param::μ>(res, ca, cb, ek);
}
template <class P = lvl1param>
void HomORYN(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek)
{
    HomGate<P, 1, -1, lvl1param::μ>(res, ca, cb, ek);
}

// 3input
// cs?c1:c0
template <class P = lvl1param>
void HomMUX(TLWE<P> &res, const TLWE<P> &cs, const TLWE<P> &c1,
            const TLWE<P> &c0, const EvalKey &ek)
{
    TLWE<P> temp;
    for (int i = 0; i <= P::k * P::n; i++) temp[i] = cs[i] + c1[i];
    for (int i = 0; i <= P::k * P::n; i++) res[i] = -cs[i] + c0[i];
    temp[P::k * P::n] -= P::μ;
    res[P::k * P::n] -= P::μ;
    if constexpr (std::is_same_v<P, lvl1param>) {
        TLWE<lvl0param> and1, and0;
        IdentityKeySwitch<lvl10param>(and1, temp, *ek.iksklvl10);
        IdentityKeySwitch<lvl10param>(and0, res, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(
            temp, and1, *ek.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(
            res, and0, *ek.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
        for (int i = 0; i <= P::k * lvl1param::n; i++) res[i] += temp[i];
        res[P::k * P::n] += P::μ;
    }
    else if constexpr (std::is_same_v<P, lvl0param>) {
        TLWE<lvl1param> and1, and0;
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(
            and1, temp, *ek.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(
            and0, res, *ek.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
        for (int i = 0; i <= lvl1param::k * lvl1param::n; i++)
            and0[i] += and1[i];
        IdentityKeySwitch<lvl10param>(res, and0, *ek.iksklvl10);
        res[P::k * P::n] += P::μ;
    }
}
template <class P = lvl1param>
void HomNMUX(TLWE<P> &res, const TLWE<P> &cs, const TLWE<P> &c1,
             const TLWE<P> &c0, const EvalKey &ek)
{
    HomMUX<P>(res, cs, c1, c0, ek);
    for (int i = 0; i <= P::k * P::n; i++) res[i] = -res[i];
}
template <class bkP>
void HomMUXwoIKSandSE(TRLWE<typename bkP::targetP> &res,
                      const TLWE<typename bkP::domainP> &cs,
                      const TLWE<typename bkP::domainP> &c1,
                      const TLWE<typename bkP::domainP> &c0, const EvalKey &ek)
{
    TLWE<typename bkP::domainP> temp1;
    TLWE<typename bkP::domainP> temp0;
    for (int i = 0; i <= bkP::domainP::n; i++) temp1[i] = cs[i] + c1[i];
    for (int i = 0; i <= bkP::domainP::n; i++) temp0[i] = -cs[i] + c0[i];
    temp1[lvl0param::n] -= bkP::domainP::μ;
    temp0[lvl0param::n] -= bkP::domainP::μ;
    TRLWE<typename bkP::targetP> and0;
    BlindRotate<bkP>(res, temp1, ek.getbkfft<bkP>(),
                     μpolygen<typename bkP::targetP, bkP::targetP::μ>());
    BlindRotate<bkP>(and0, temp0, ek.getbkfft<bkP>(),
                     μpolygen<typename bkP::targetP, bkP::targetP::μ>());

    for (int i = 0; i < bkP::targetP::n; i++) {
        res[0][i] += and0[0][i];
        res[1][i] += and0[1][i];
    };
    res[1][0] += bkP::targetP::μ;
}
template <class iksP, class bkP>
void HomMUXwoSE(TRLWE<typename bkP::targetP> &res,
                const TLWE<typename iksP::domainP> &cs,
                const TLWE<typename iksP::domainP> &c1,
                const TLWE<typename iksP::domainP> &c0, const EvalKey &ek)
{
    TLWE<typename iksP::domainP> temp1;
    TLWE<typename iksP::domainP> temp0;
    for (int i = 0; i <= iksP::domainP::n; i++) temp1[i] = cs[i] + c1[i];
    for (int i = 0; i <= iksP::domainP::n; i++) temp0[i] = -cs[i] + c0[i];
    temp1[iksP::domainP::n] -= iksP::domainP::μ;
    temp0[iksP::domainP::n] -= iksP::domainP::μ;
    TLWE<lvl0param> and1, and0;
    IdentityKeySwitch<iksP>(and1, temp1, ek.getiksk<iksP>());
    IdentityKeySwitch<iksP>(and0, temp0, ek.getiksk<iksP>());
    TRLWE<typename bkP::targetP> and0trlwe;
    BlindRotate<bkP>(res, and1, ek.getbkfft<bkP>(),
                     μpolygen<typename bkP::targetP, bkP::targetP::μ>());
    BlindRotate<bkP>(and0trlwe, and0, ek.getbkfft<bkP>(),
                     μpolygen<typename bkP::targetP, bkP::targetP::μ>());

    for (int i = 0; i < bkP::targetP::n; i++) {
        res[0][i] += and0trlwe[0][i];
        res[1][i] += and0trlwe[1][i];
    };
    res[1][0] += bkP::targetP::μ;
}

void ExtractSwitchAndHomMUX(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &csr,
                            const TRLWE<lvl1param> &c1r,
                            const TRLWE<lvl1param> &c0r, const EvalKey &ek);

#include "externs/gate.hpp"
}  // namespace TFHEpp