#include <circuitbootstrapping.hpp>
namespace TFHEpp {
#define INST(iksP, bkP, privksP)                            \
    template void CircuitBootstrapping<iksP, bkP, privksP>( \
        TRGSW<typename privksP::targetP> & trgsw,           \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

#define INST(iksP, bkP, privksP)                               \
    template void CircuitBootstrappingFFT<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & trgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

#define INST(iksP, bkP, privksP)                               \
    template void CircuitBootstrappingSub<iksP, bkP, privksP>( \
        TRGSW<typename privksP::targetP> & trgsw,              \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING_SUBIKS(INST)
#undef INST

#define INST(iksP, bkP, privksP)                                  \
    template void CircuitBootstrappingSubFFT<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & trgswfft,           \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING_SUBIKS(INST)
#undef INST

#define INST(iksP, bkP, privksP)                                  \
    template void CircuitBootstrappingFFTInv<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & invtrgswfft,        \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

#define INST(iksP, bkP, privksP)                                      \
    template void CircuitBootstrappingFFTwithInv<iksP, bkP, privksP>( \
        TRGSWFFT<typename privksP::targetP> & trgswfft,               \
        TRGSWFFT<typename privksP::targetP> & invtrgswfft,            \
        const TLWE<typename iksP::domainP> &tlwe, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST
}  // namespace TFHEpp