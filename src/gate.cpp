#include <gate.hpp>
#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {

// No input

template <class P>
void HomCONSTANTONE(TLWE<lvl1param> &res)
{
    res = {};
    res[lvl1param::n] = P::μ;
}
#define INST(P)                                               \
    template void HomCONSTANTONE<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomCONSTANTZERO(TLWE<lvl1param> &res)
{
    res = {};
    res[lvl1param::n] = -lvl1param::μ;
}
#define INST(P)                                               \
    template void HomCONSTANTZERO<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvlMparam);
#undef INST

// 1 input
void HomNOT(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca)
{
    for (int i = 0; i <= lvl1param::n; i++) res[i] = -ca[i];
}

void HomCOPY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca)
{
    for (int i = 0; i <= lvl1param::n; i++) res[i] = ca[i];
}

template <class P, int casign, int cbsign, typename P::T offset>
inline void HomGate(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                    const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    for (int i = 0; i <= lvl1param::n; i++)
        res[i] = casign * ca[i] + cbsign * cb[i];
    res[lvl1param::n] += offset;
    GateBootstrapping<P::μ>(res, res, ek);
}

template <class P>
void HomNAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    HomGate<P, -1, -1, P::μ>(res, ca, cb, ek);
}
#define INST(P)                                                \
    template void HomNAND<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    HomGate<P, -1, -1, -P::μ>(res, ca, cb, ek);
}
#define INST(P)                                                \
    template void HomNOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomXNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    if constexpr (std::is_same_v<P, lvl1param>) {
        HomGate<P,-2, -2, -2 * lvl1param::μ>(res, ca, cb, ek);
    }
    else if constexpr (std::is_same_v<P, lvlMparam>) {
        HomGate<P,-3, -3, -3 * lvlMparam::μ>(res, ca, cb, ek);
    }
}
#define INST(P)                                                \
    template void HomXNOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    HomGate<P, 1, 1, -P::μ>(res, ca, cb, ek);
}
#define INST(P)                                                \
    template void HomAND<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
           const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    HomGate<P,1, 1, P::μ>(res, ca, cb, ek);
}
#define INST(P)                                                \
    template void HomOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomXOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    if constexpr (std::is_same_v<P, lvl1param>) {
        HomGate<P, 2, 2, 2 * lvl1param::μ>(res, ca, cb, ek);
    }
    else if constexpr (std::is_same_v<P, lvlMparam>) {
        HomGate<P, 3, 3, 3 * lvlMparam::μ>(res, ca, cb, ek);
    }
}
#define INST(P)                                               \
    template void HomXOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomANDNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    HomGate<P, -1, 1, -P::μ>(res, ca, cb, ek);
}
#define INST(P)                                               \
    template void HomANDNY<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomANDYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    HomGate<P, 1, -1, -P::μ>(res, ca, cb, ek);
}
#define INST(P)                                               \
    template void HomANDYN<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST


template <class P>
void HomORNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    HomGate<P, -1, 1, P::μ>(res, ca, cb, ek);
}
#define INST(P)                                               \
    template void HomORNY<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST


template <class P>
void HomORYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    HomGate<P, 1, -1, P::μ>(res, ca, cb, ek);
}
#define INST(P)                                               \
    template void HomORYN<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST


// 3input
// cs?c1:c0
template <class P>
void HomMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
            const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
            const EvalKey &ek)
{
    if constexpr (std::is_same_v<P, lvl1param>) {
        TLWE<lvl1param> temp;
        for (int i = 0; i <= lvl1param::n; i++) temp[i] = cs[i] + c1[i];
        for (int i = 0; i <= lvl1param::n; i++) res[i] = -cs[i] + c0[i];
        temp[lvl1param::n] -= P::μ;
        res[lvl1param::n] -= P::μ;
        TLWE<lvl0param> and1, and0;
        IdentityKeySwitch<lvl10param>(and1, temp, *ek.iksklvl10);
        IdentityKeySwitch<lvl10param>(and0, res, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(
            temp, and1, *ek.bkfftlvl01, μpolygen<lvl1param, P::μ>());
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(
            res, and0, *ek.bkfftlvl01, μpolygen<lvl1param, P::μ>());

        for (int i = 0; i <= lvl1param::n; i++) res[i] += temp[i];
        res[lvl1param::n] += P::μ;
    }
    else if constexpr (std::is_same_v<P, lvlMparam>) {
        TLWE<lvlMparam> temp;
        HomANDNY<lvlMparam>(temp, cs, c0, ek);
        HomAO3<lvlMparam>(res, c1, cs, temp, ek);
    }
}
#define INST(P)                                                               \
    template void HomMUX<P>(TLWE<lvl1param> & res, const TLWE<lvl1param> &cs, \
                            const TLWE<lvl1param> &c1,                        \
                            const TLWE<lvl1param> &c0, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

template <class P>
void HomNMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
             const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
             const EvalKey &ek)
{
    HomMUX<P>(res, cs, c1, c0, ek);
    HomNOT(res, res);
}
#define INST(P)                                                                \
    template void HomNMUX<P>(TLWE<lvl1param> & res, const TLWE<lvl1param> &cs, \
                             const TLWE<lvl1param> &c1,                        \
                             const TLWE<lvl1param> &c0, const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

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
#define INST(iksP, bkP)                         \
    template void HomMUXwoSE<iksP, bkP>(        \
        TRLWE<typename bkP::targetP> & res,     \
        const TLWE<typename iksP::domainP> &cs, \
        const TLWE<typename iksP::domainP> &c1, \
        const TLWE<typename iksP::domainP> &c0, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE(INST)
#undef INST

void ExtractSwitchAndHomMUX(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &csr,
                            const TRLWE<lvl1param> &c1r,
                            const TRLWE<lvl1param> &c0r, const EvalKey &ek)
{
    TLWE<lvl1param> templvl1;
    TLWE<lvl0param> cs, c1, c0;
    SampleExtractIndex<lvl1param>(templvl1, csr, 0);
    IdentityKeySwitch<lvl10param>(cs, templvl1, *ek.iksklvl10);
    SampleExtractIndex<lvl1param>(templvl1, c1r, 0);
    IdentityKeySwitch<lvl10param>(c1, templvl1, *ek.iksklvl10);
    SampleExtractIndex<lvl1param>(templvl1, c0r, 0);
    IdentityKeySwitch<lvl10param>(c0, templvl1, *ek.iksklvl10);

    for (int i = 0; i <= lvl0param::n; i++) c1[i] += cs[i];
    for (int i = 0; i <= lvl0param::n; i++) c0[i] -= cs[i];
    c1[lvl0param::n] -= lvl0param::μ;
    c0[lvl0param::n] -= lvl0param::μ;
    TRLWE<lvl1param> and0;
    BlindRotate<lvl01param>(res, c1, *ek.bkfftlvl01,
                            μpolygen<lvl1param, lvl1param::μ>());
    BlindRotate<lvl01param>(and0, c0, *ek.bkfftlvl01,
                            μpolygen<lvl1param, lvl1param::μ>());

    for (int i = 0; i < lvl1param::n; i++) {
        res[0][i] += and0[0][i];
        res[1][i] += and0[1][i];
    };
    res[1][0] += lvl1param::n;
}

}  // namespace TFHEpp
