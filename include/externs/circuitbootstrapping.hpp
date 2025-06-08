#pragma once
#import "../circuitbootstrapping.hpp"

namespace TFHEpp{
#define INST(iksP, bkP, privksP)                            \
        TRGSW<typename privksP::targetP> & trgsw,           \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
#undef INST

#define INST(iksP, bkP, privksP)                               \
        TRGSWFFT<typename privksP::targetP> & trgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
#undef INST

#define INST(iksP, bkP, privksP)                            \
        TRGSW<typename privksP::targetP> & trgsw,           \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
#undef INST

#define INST(iksP, bkP, privksP)                               \
        TRGSWFFT<typename privksP::targetP> & trgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
#undef INST

#define INST(iksP, bkP, privksP)                                  \
        TRGSWFFT<typename privksP::targetP> & invtrgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
#undef INST

#define INST(iksP, bkP, privksP)                                      \
        TRGSWFFT<typename privksP::targetP> & trgswfft,               \
        TRGSWFFT<typename privksP::targetP> & invtrgswfft,            \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
#undef INST
}