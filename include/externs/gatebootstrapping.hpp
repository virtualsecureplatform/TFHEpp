#pragma once
#include"../gatebootstrapping.hpp"

namespace TFHEpp{
#define INST(P)                                     \
    extern template void GateBootstrappingTLWE2TLWEFFT<P>( \
        TLWE<typename P::targetP> & res,            \
        const TLWE<typename P::domainP> &tlwe,      \
        const BootstrappingKeyFFT<P> &bkfft,        \
        const Polynomial<typename P::targetP> &testvector)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST
}