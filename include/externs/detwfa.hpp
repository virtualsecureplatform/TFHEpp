#pragma once
#import "../detwfa.hpp"

namespace TFHEpp{
#define INST(P)                                                     \
                             const TRLWE<P> &c1, const TRLWE<P> &c0)
#undef INST

#define INST(bkP)                                             \
        TRLWE<typename bkP::targetP> & acc,                   \
        const BootstrappingKeyElementFFT<bkP> &cs, const int a)
#undef INST

#define INST(P)                                             \
        TRLWE<P> & acc, const TRGSWNTT<P> &cs, const typename P::T a)
#undef INST
}