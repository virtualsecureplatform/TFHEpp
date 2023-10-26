#include <gate.hpp>

namespace TFHEpp {

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

#define INST(P) template void HomCONSTANTONE<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) template void HomCONSTANTZERO<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) template void HomNOT<P>(TLWE<P> & res, const TLWE<P> &ca)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) template void HomCOPY<P>(TLWE<P> & res, const TLWE<P> &ca)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(iksP, brP, μ)                      \
    template void HomNAND<iksP, brP, μ>(        \
        TLWE<typename brP::targetP> & res,      \
        const TLWE<typename iksP::domainP> &ca, \
        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                     \
    template void HomNAND<brP, μ, iksP>(TLWE<typename iksP::targetP> & res,    \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, \
                                        const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                     \
    template void HomNOR<iksP, brP, μ>(TLWE<typename brP::targetP> & res,      \
                                       const TLWE<typename iksP::domainP> &ca, \
                                       const TLWE<typename iksP::domainP> &cb, \
                                       const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                    \
    template void HomNOR<brP, μ, iksP>(TLWE<typename iksP::targetP> & res,    \
                                       const TLWE<typename brP::domainP> &ca, \
                                       const TLWE<typename brP::domainP> &cb, \
                                       const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                      \
    template void HomXNOR<iksP, brP, μ>(        \
        TLWE<typename brP::targetP> & res,      \
        const TLWE<typename iksP::domainP> &ca, \
        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                     \
    template void HomXNOR<brP, μ, iksP>(TLWE<typename iksP::targetP> & res,    \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, \
                                        const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                     \
    template void HomAND<iksP, brP, μ>(TLWE<typename brP::targetP> & res,      \
                                       const TLWE<typename iksP::domainP> &ca, \
                                       const TLWE<typename iksP::domainP> &cb, \
                                       const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                    \
    template void HomAND<brP, μ, iksP>(TLWE<typename iksP::targetP> & res,    \
                                       const TLWE<typename brP::domainP> &ca, \
                                       const TLWE<typename brP::domainP> &cb, \
                                       const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                    \
    template void HomOR<iksP, brP, μ>(TLWE<typename brP::targetP> & res,      \
                                      const TLWE<typename iksP::domainP> &ca, \
                                      const TLWE<typename iksP::domainP> &cb, \
                                      const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                   \
    template void HomOR<brP, μ, iksP>(TLWE<typename iksP::targetP> & res,    \
                                      const TLWE<typename brP::domainP> &ca, \
                                      const TLWE<typename brP::domainP> &cb, \
                                      const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                     \
    template void HomXOR<iksP, brP, μ>(TLWE<typename brP::targetP> & res,      \
                                       const TLWE<typename iksP::domainP> &ca, \
                                       const TLWE<typename iksP::domainP> &cb, \
                                       const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                    \
    template void HomXOR<brP, μ, iksP>(TLWE<typename iksP::targetP> & res,    \
                                       const TLWE<typename brP::domainP> &ca, \
                                       const TLWE<typename brP::domainP> &cb, \
                                       const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                      \
    template void HomANDNY<iksP, brP, μ>(       \
        TLWE<typename brP::targetP> & res,      \
        const TLWE<typename iksP::domainP> &ca, \
        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                     \
    template void HomANDNY<brP, μ, iksP>(      \
        TLWE<typename iksP::targetP> & res,    \
        const TLWE<typename brP::domainP> &ca, \
        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                      \
    template void HomANDYN<iksP, brP, μ>(       \
        TLWE<typename brP::targetP> & res,      \
        const TLWE<typename iksP::domainP> &ca, \
        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                     \
    template void HomANDYN<brP, μ, iksP>(      \
        TLWE<typename iksP::targetP> & res,    \
        const TLWE<typename brP::domainP> &ca, \
        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                      \
    template void HomORNY<iksP, brP, μ>(        \
        TLWE<typename brP::targetP> & res,      \
        const TLWE<typename iksP::domainP> &ca, \
        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                     \
    template void HomORNY<brP, μ, iksP>(TLWE<typename iksP::targetP> & res,    \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, \
                                        const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                      \
    template void HomORYN<iksP, brP, μ>(        \
        TLWE<typename brP::targetP> & res,      \
        const TLWE<typename iksP::domainP> &ca, \
        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                     \
    template void HomORYN<brP, μ, iksP>(TLWE<typename iksP::targetP> & res,    \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, \
                                        const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(P)                                                   \
    template void HomMUX<P>(TLWE<P> & res, const TLWE<P> &cs,     \
                            const TLWE<P> &c1, const TLWE<P> &c0, \
                            const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                    \
    template void HomNMUX<P>(TLWE<P> & res, const TLWE<P> &cs,     \
                             const TLWE<P> &c1, const TLWE<P> &c0, \
                             const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(bkP)                                                              \
    template void HomMUXwoIKSandSE<bkP>(TRLWE<typename bkP::targetP> & res,    \
                                        const TLWE<typename bkP::domainP> &cs, \
                                        const TLWE<typename bkP::domainP> &c1, \
                                        const TLWE<typename bkP::domainP> &c0, \
                                        const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(iksP, bkP)                         \
    template void HomMUXwoSE<iksP, bkP>(        \
        TRLWE<typename bkP::targetP> & res,     \
        const TLWE<typename iksP::domainP> &cs, \
        const TLWE<typename iksP::domainP> &c1, \
        const TLWE<typename iksP::domainP> &c0, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE(INST)
#undef INST

}  // namespace TFHEpp