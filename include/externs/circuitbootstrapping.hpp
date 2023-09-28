#pragma once
#include"../circuitbootstrapping.hpp"

namespace TFHEpp{
#define INST(iksP, bkP, privksP)                            \
    extern template void CircuitBootstrapping<iksP, bkP, privksP>( \
        TRGSW<typename privksP::targetP> & trgsw,           \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

#define INST(iksP, bkP, privksP)                               \
    extern template void CircuitBootstrappingFFT<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & trgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

#define INST(iksP, bkP, privksP)                            \
    extern template void CircuitBootstrappingSub<iksP, bkP, privksP>( \
        TRGSW<typename privksP::targetP> & trgsw,           \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING_SUBIKS(INST)
#undef INST

#define INST(iksP, bkP, privksP)                               \
    extern template void CircuitBootstrappingSubFFT<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & trgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING_SUBIKS(INST)
#undef INST

#define INST(iksP, bkP, privksP)                                  \
    extern template void CircuitBootstrappingFFTInv<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & invtrgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

#define INST(iksP, bkP, privksP)                                      \
    extern template void CircuitBootstrappingFFTwithInv<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & trgswfft,               \
        TRGSWFFT<typename privksP::targetP> & invtrgswfft,            \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST
}