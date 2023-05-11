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

#define INST(P)                                                \
    template void HomNAND<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                               \
    template void HomNOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                \
    template void HomXNOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                               \
    template void HomAND<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                              \
    template void HomOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                           const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                               \
    template void HomXOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                 \
    template void HomANDNY<P>(TLWE<P> & res, const TLWE<P> &ca, \
                              const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                 \
    template void HomANDYN<P>(TLWE<P> & res, const TLWE<P> &ca, \
                              const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                \
    template void HomORNY<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                \
    template void HomORYN<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
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