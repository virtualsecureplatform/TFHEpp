#pragma once

#define INST(P) extern template void HomCONSTANTONE<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) extern template void HomCONSTANTZERO<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) extern template void HomNOT<P>(TLWE<P> & res, const TLWE<P> &ca)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) extern template void HomCOPY<P>(TLWE<P> & res, const TLWE<P> &ca)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                \
    extern template void HomNAND<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                               \
    extern template void HomNOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                \
    extern template void HomXNOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                               \
    extern template void HomAND<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                              \
    extern template void HomOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                           const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                               \
    extern template void HomXOR<P>(TLWE<P> & res, const TLWE<P> &ca, \
                            const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                 \
    extern template void HomANDNY<P>(TLWE<P> & res, const TLWE<P> &ca, \
                              const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                 \
    extern template void HomANDYN<P>(TLWE<P> & res, const TLWE<P> &ca, \
                              const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                \
    extern template void HomORNY<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                \
    extern template void HomORYN<P>(TLWE<P> & res, const TLWE<P> &ca, \
                             const TLWE<P> &cb, const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                   \
    extern template void HomMUX<P>(TLWE<P> & res, const TLWE<P> &cs,     \
                            const TLWE<P> &c1, const TLWE<P> &c0, \
                            const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                    \
    extern template void HomNMUX<P>(TLWE<P> & res, const TLWE<P> &cs,     \
                             const TLWE<P> &c1, const TLWE<P> &c0, \
                             const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST


#define INST(bkP)                                                              \
    extern template void HomMUXwoIKSandSE<bkP>(TRLWE<typename bkP::targetP> & res,    \
                                        const TLWE<typename bkP::domainP> &cs, \
                                        const TLWE<typename bkP::domainP> &c1, \
                                        const TLWE<typename bkP::domainP> &c0, \
                                        const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST


#define INST(iksP, bkP)                         \
    extern template void HomMUXwoSE<iksP, bkP>(        \
        TRLWE<typename bkP::targetP> & res,     \
        const TLWE<typename iksP::domainP> &cs, \
        const TLWE<typename iksP::domainP> &c1, \
        const TLWE<typename iksP::domainP> &c0, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE(INST)
#undef INST