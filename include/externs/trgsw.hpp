#pragma once
#import "../trgsw.hpp"

namespace TFHEpp{
#define INST(P)                                                       \
        DecomposedPolynomial<P> & decpoly, const Polynomial<P> &poly)
#undef INST

#define INST(P)                                                       \
        DecomposedNoncePolynomial<P> & decpoly, const Polynomial<P> &poly)
#undef INST

#define INST(P)                               \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWFFT<P> &trgswfft)
#undef INST

#define INST(P)                               \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWNTT<P> &trgswntt)
#undef INST

#define INST(P) template TRGSWFFT<P> ApplyFFT2trgsw<P>(const TRGSW<P> &trgsw)
#undef INST

#undef INST

// TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P)                                                  \
                                         const Key<P> &key)
#undef INST

#define INST(P)                                 \
        const Polynomial<P> &p, const Key<P> &key)
#undef INST

#define INST(P)                                 \
        const Polynomial<P> &p, const Key<P> &key)
#undef INST
}