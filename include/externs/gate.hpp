#pragma once
#include"../gate.hpp"

namespace TFHEpp{
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

#define INST(iksP, brP, μ)                                                \
    extern template void HomNAND<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomNAND<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomNOR<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomNOR<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomXNOR<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomXNOR<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomAND<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomAND<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomOR<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomOR<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomXOR<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomXOR<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomANDNY<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomANDNY<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomANDYN<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomANDYN<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomORNY<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomORNY<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
#undef INST

#define INST(iksP, brP, μ)                                                \
    extern template void HomORYN<iksP, brP, μ>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST
#define INST(brP, μ, iksP)                                                \
    extern template void HomORYN<brP, μ, iksP>(TLWE<typename iksP::targetP> &res, \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(INST)
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
}