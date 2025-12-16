#pragma once
#include"../trgsw.hpp"

namespace TFHEpp{
#define INST(P)                                                       \
    extern template void Decomposition<P>(                         \
        DecomposedPolynomial<P> & decpoly, const Polynomial<P> &poly)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P)                                                       \
    extern template void NonceDecomposition<P>(                         \
        DecomposedNoncePolynomial<P> & decpoly, const Polynomial<P> &poly)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P)                          \
    extern template void ExternalProduct<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWFFT<P> &trgswfft)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P)                          \
    extern template void ExternalProduct<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWNTT<P> &trgswntt)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P) template TRGSWFFT<P> ApplyFFT2trgsw<P>(const TRGSW<P> &trgsw)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P) extern template TRGSWNTT<P> ApplyNTT2trgsw<P>(const TRGSW<P> &trgsw)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P) extern template TRGSWNTT<P> TRGSW2NTT<P>(const TRGSW<P> &trgsw)
// TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
INST(lvl1param);
#undef INST

#define INST(P)                                                  \
    extern template TRGSW<P> trgswSymEncrypt<P>(const Polynomial<P> &p, \
                                         const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P)                                                            \
    extern template void trgswSymEncrypt<P>(TRGSWFFT<P> &trgswfft,          \
                                            const Polynomial<P> &p,        \
                                            const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P)                                                            \
    extern template void trgswSymEncrypt<P>(TRGSWNTT<P> &trgswntt,         \
                                            const Polynomial<P> &p,        \
                                            const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST
}