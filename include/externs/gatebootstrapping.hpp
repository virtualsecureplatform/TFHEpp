#pragma once
#import "../gatebootstrapping.hpp"

namespace TFHEpp{
#define INST(P)                                     \
        TLWE<typename P::targetP> & res,            \
        const TLWE<typename P::domainP> &tlwe,      \
        const BootstrappingKeyFFT<P> &bkfft,        \
        const Polynomial<typename P::targetP> &testvector)
#undef INST
}