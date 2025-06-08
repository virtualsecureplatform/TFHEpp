#pragma once
#import "../gate.hpp"

namespace TFHEpp{
#undef INST

#undef INST

#undef INST

#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(iksP, brP, μ)                                                \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
#undef INST
#define INST(brP, μ, iksP)                                                \
                                        const TLWE<typename brP::domainP> &ca, \
                                        const TLWE<typename brP::domainP> &cb, const EvalKey &ek)
#undef INST

#define INST(P)                                                   \
                            const TLWE<P> &c1, const TLWE<P> &c0, \
                            const EvalKey &ek)
#undef INST

#define INST(P)                                                    \
                             const TLWE<P> &c1, const TLWE<P> &c0, \
                             const EvalKey &ek)
#undef INST


#define INST(bkP)                                                              \
                                        const TLWE<typename bkP::domainP> &cs, \
                                        const TLWE<typename bkP::domainP> &c1, \
                                        const TLWE<typename bkP::domainP> &c0, \
                                        const EvalKey &ek)
#undef INST


#define INST(iksP, bkP)                         \
        TRLWE<typename bkP::targetP> & res,     \
        const TLWE<typename iksP::domainP> &cs, \
        const TLWE<typename iksP::domainP> &c1, \
        const TLWE<typename iksP::domainP> &c0, const EvalKey &ek)
#undef INST
}