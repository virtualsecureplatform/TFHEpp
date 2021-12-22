#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>

namespace TFHEpp {

// No input
template <class P>
void HomCONSTANTONE(TLWE<P> &res)
{
    res = {};
    res[P::k * P::n] = P::μ;
}
#define INST(P) template void HomCONSTANTONE<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomCONSTANTZERO(TLWE<P> &res)
{
    res = {};
    res[P::k * P::n] = -P::μ;
}
#define INST(P) template void HomCONSTANTZERO<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvl0param);
#undef INST

// 1 input
template <class P>
void HomNOT(TLWE<P> &res, const TLWE<P> &ca)
{
    for (int i = 0; i <= P::k * P::n; i++) res[i] = -ca[i];
}
#define INST(P) template void HomNOT<P>(TLWE<P> & res, const TLWE<P> &ca)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomCOPY(TLWE<P> &res, const TLWE<P> &ca)
{
    for (int i = 0; i <= P::k * P::n; i++) res[i] = ca[i];
}
#define INST(P) template void HomCOPY<P>(TLWE<P> & res, const TLWE<P> &ca)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P, int casign, int cbsign, typename P::T offset>
inline void HomGate(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
                    const EvalKey &ek)
{
    for (int i = 0; i <= P::k * P::n; i++)
        res[i] = casign * ca[i] + cbsign * cb[i];
    res[P::k * P::n] += offset;
    GateBootstrapping(res, res, ek);
}

template <class P>
void HomNAND(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek)
{
    HomGate<P, -1, -1, lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                                \
    template void HomNAND<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomNOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek)
{
    HomGate<P, -1, -1, -lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                               \
    template void HomNOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomXNOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek)
{
    HomGate<P, -2, -2, -2 * lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                                \
    template void HomXNOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomAND(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek)
{
    HomGate<P, 1, 1, -lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                               \
    template void HomAND<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
           const EvalKey &ek)
{
    HomGate<P, 1, 1, lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                              \
    template void HomOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                           const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomXOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek)
{
    HomGate<P, 2, 2, 2 * lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                               \
    template void HomXOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomANDNY(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
              const EvalKey &ek)
{
    HomGate<P, -1, 1, -lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                                 \
    template void HomANDNY<P>(TLWE<P> & res, const TLWE<P> &ca, \
                              const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomANDYN(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
              const EvalKey &ek)
{
    HomGate<P, 1, -1, -lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                                 \
    template void HomANDYN<P>(TLWE<P> & res, const TLWE<P> &ca, \
                              const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomORNY(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek)
{
    HomGate<P, -1, 1, lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                                \
    template void HomORNY<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomORYN(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek)
{
    HomGate<P, 1, -1, lvl1param::μ>(res, ca, cb, ek);
}
#define INST(P)                                                \
    template void HomORYN<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

// 3input
// cs?c1:c0
template <class P>
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
#define INST(P)                                                   \
    template void HomMUX<P>(TLWE<P> & res, const TLWE<P> &cs,     \
                            const TLWE<P> &c1, const TLWE<P> &c0, \
                            const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

template <class P>
void HomNMUX(TLWE<P> &res, const TLWE<P> &cs, const TLWE<P> &c1,
             const TLWE<P> &c0, const EvalKey &ek)
{
    HomMUX<P>(res, cs, c1, c0, ek);
    for (int i = 0; i <= P::k * P::n; i++) res[i] = -res[i];
}
#define INST(P)                                                    \
    template void HomNMUX<P>(TLWE<P> & res, const TLWE<P> &cs,     \
                             const TLWE<P> &c1, const TLWE<P> &c0, \
                             const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

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
#define INST(bkP)                                                              \
    template void HomMUXwoIKSandSE<bkP>(TRLWE<typename bkP::targetP> & res,    \
                                        const TLWE<typename bkP::domainP> &cs, \
                                        const TLWE<typename bkP::domainP> &c1, \
                                        const TLWE<typename bkP::domainP> &c0, \
                                        const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
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