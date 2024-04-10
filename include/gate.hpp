#pragma once

#include "cloudkey.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"

namespace TFHEpp {
template <class brP, typename brP::targetP::T μ, class iksP, int casign,
          int cbsign, std::make_signed_t<typename brP::domainP::T> offset>
inline void HomGate(TLWE<typename iksP::targetP> &res,
                    const TLWE<typename brP::domainP> &ca,
                    const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    for (int i = 0; i <= brP::domainP::k * brP::domainP::n; i++)
        res[i] = casign * ca[i] + cbsign * cb[i];
    res[brP::domainP::k * brP::domainP::n] += offset;
    GateBootstrapping<brP, μ, iksP>(res, res, ek);
}

template <class iksP, class brP, typename brP::targetP::T μ, int casign,
          int cbsign, std::make_signed_t<typename iksP::domainP::T> offset>
inline void HomGate(TLWE<typename brP::targetP> &res,
                    const TLWE<typename iksP::domainP> &ca,
                    const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    for (int i = 0; i <= iksP::domainP::k * iksP::domainP::n; i++)
        res[i] = casign * ca[i] + cbsign * cb[i];
    res[iksP::domainP::k * iksP::domainP::n] += offset;
    GateBootstrapping<iksP, brP, μ>(res, res, ek);
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

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomNAND(TLWE<typename iksP::targetP> &res,
             const TLWE<typename brP::domainP> &ca,
             const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, -1, -1, brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomNAND(TLWE<typename brP::targetP> &res,
             const TLWE<typename iksP::domainP> &ca,
             const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, -1, -1, iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomNOR(TLWE<typename iksP::targetP> &res,
            const TLWE<typename brP::domainP> &ca,
            const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, -1, -1, -brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomNOR(TLWE<typename brP::targetP> &res,
            const TLWE<typename iksP::domainP> &ca,
            const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, -1, -1, -iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomXNOR(TLWE<typename iksP::targetP> &res,
             const TLWE<typename brP::domainP> &ca,
             const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, -2, -2, -2 * brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomXNOR(TLWE<typename brP::targetP> &res,
             const TLWE<typename iksP::domainP> &ca,
             const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, -2, -2, -2 * iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomAND(TLWE<typename iksP::targetP> &res,
            const TLWE<typename brP::domainP> &ca,
            const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, 1, 1, -brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomAND(TLWE<typename brP::targetP> &res,
            const TLWE<typename iksP::domainP> &ca,
            const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, 1, 1, -iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomOR(TLWE<typename iksP::targetP> &res,
           const TLWE<typename brP::domainP> &ca,
           const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, 1, 1, brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomOR(TLWE<typename brP::targetP> &res,
           const TLWE<typename iksP::domainP> &ca,
           const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, 1, 1, iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomXOR(TLWE<typename iksP::targetP> &res,
            const TLWE<typename brP::domainP> &ca,
            const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, 2, 2, 2 * brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomXOR(TLWE<typename brP::targetP> &res,
            const TLWE<typename iksP::domainP> &ca,
            const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, 2, 2, 2 * iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomANDNY(TLWE<typename iksP::targetP> &res,
              const TLWE<typename brP::domainP> &ca,
              const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, -1, 1, -brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomANDNY(TLWE<typename brP::targetP> &res,
              const TLWE<typename iksP::domainP> &ca,
              const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, -1, 1, -iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomANDYN(TLWE<typename iksP::targetP> &res,
              const TLWE<typename brP::domainP> &ca,
              const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, 1, -1, -brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomANDYN(TLWE<typename brP::targetP> &res,
              const TLWE<typename iksP::domainP> &ca,
              const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, 1, -1, -iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomORNY(TLWE<typename iksP::targetP> &res,
             const TLWE<typename brP::domainP> &ca,
             const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, -1, 1, brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomORNY(TLWE<typename brP::targetP> &res,
             const TLWE<typename iksP::domainP> &ca,
             const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, -1, 1, iksP::domainP::μ>(res, ca, cb, ek);
}

template <class brP = lvl01param, typename brP::targetP::T μ = lvl1param::μ,
          class iksP = lvl10param>
void HomORYN(TLWE<typename iksP::targetP> &res,
             const TLWE<typename brP::domainP> &ca,
             const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
{
    HomGate<brP, μ, iksP, 1, -1, brP::domainP::μ>(res, ca, cb, ek);
}
template <class iksP = lvl10param, class brP = lvl01param,
          typename brP::targetP::T μ = lvl1param::μ>
void HomORYN(TLWE<typename brP::targetP> &res,
             const TLWE<typename iksP::domainP> &ca,
             const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
{
    HomGate<iksP, brP, μ, 1, -1, iksP::domainP::μ>(res, ca, cb, ek);
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

template <class brP, typename brP::targetP::T μ = brP::targetP::μ>
void HomMUXwoSE(TRLWE<typename brP::targetP> &res,
                const TLWE<typename brP::domainP> &cs,
                const TLWE<typename brP::domainP> &c1,
                const TLWE<typename brP::domainP> &c0, const EvalKey &ek)
{
    TLWE<typename brP::domainP> and1, and0;
    for (int i = 0; i <= brP::domainP::k * brP::domainP::n; i++)
        and1[i] = cs[i] + c1[i];
    for (int i = 0; i <= brP::domainP::k * brP::domainP::n; i++)
        and0[i] = -cs[i] + c0[i];
    and1[brP::domainP::k * brP::domainP::n] -= brP::domainP::μ;
    and0[brP::domainP::k * brP::domainP::n] -= brP::domainP::μ;
    TRLWE<typename brP::targetP> and0trlwe;
    BlindRotate<brP>(res, and1, ek.getbkfft<brP>(),
                     μpolygen<typename brP::targetP, brP::targetP::μ>());
    BlindRotate<brP>(and0trlwe, and0, ek.getbkfft<brP>(),
                     μpolygen<typename brP::targetP, brP::targetP::μ>());

    for (int i = 0; i < brP::targetP::k * brP::targetP::n; i++) {
        res[0][i] += and0trlwe[0][i];
        res[1][i] += and0trlwe[1][i];
    };
    res[1][0] += brP::targetP::μ;
}

template <class iksP, class brP, typename brP::targetP::T μ = brP::targetP::μ>
void HomMUXwoSE(TRLWE<typename brP::targetP> &res,
                const TLWE<typename iksP::domainP> &cs,
                const TLWE<typename iksP::domainP> &c1,
                const TLWE<typename iksP::domainP> &c0, const EvalKey &ek)
{
    TLWE<typename iksP::domainP> temp1, temp0;
    for (int i = 0; i <= iksP::domainP::k * iksP::domainP::n; i++)
        temp1[i] = cs[i] + c1[i];
    for (int i = 0; i <= iksP::domainP::k * iksP::domainP::n; i++)
        temp0[i] = -cs[i] + c0[i];
    temp1[iksP::domainP::k * iksP::domainP::n] -= iksP::domainP::μ;
    temp0[iksP::domainP::k * iksP::domainP::n] -= iksP::domainP::μ;
    TLWE<typename iksP::targetP> and1, and0;
    IdentityKeySwitch<iksP>(and1, temp1, ek.getiksk<iksP>());
    IdentityKeySwitch<iksP>(and0, temp0, ek.getiksk<iksP>());
    TRLWE<typename brP::targetP> and0trlwe;
    BlindRotate<brP>(res, and1, ek.getbkfft<brP>(),
                     μpolygen<typename brP::targetP, brP::targetP::μ>());
    BlindRotate<brP>(and0trlwe, and0, ek.getbkfft<brP>(),
                     μpolygen<typename brP::targetP, brP::targetP::μ>());

    for (int i = 0; i < brP::targetP::k * brP::targetP::n; i++) {
        res[0][i] += and0trlwe[0][i];
        res[1][i] += and0trlwe[1][i];
    };
    res[1][0] += brP::targetP::μ;
}

void ExtractSwitchAndHomMUX(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &csr,
                            const TRLWE<lvl1param> &c1r,
                            const TRLWE<lvl1param> &c0r, const EvalKey &ek);

}  // namespace TFHEpp