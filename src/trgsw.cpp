#import "trgsw.hpp"

namespace TFHEpp {
#define INST(P)                                                       \
    template void Decomposition<P>(DecomposedPolynomial<P> & decpoly, \
                                   const Polynomial<P> &poly)
#undef INST

#define INST(P)                          \
    template void NonceDecomposition<P>( \
        DecomposedNoncePolynomial<P> & decpoly, const Polynomial<P> &poly)
#undef INST

#define INST(P)                               \
    template void trgswfftExternalProduct<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWFFT<P> &trgswfft)
#undef INST

#define INST(P)                               \
    template void trgswnttExternalProduct<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWNTT<P> &trgswntt)
#undef INST

#define INST(P) template TRGSWNTT<P> ApplyNTT2trgsw<P>(const TRGSW<P> &trgsw)
#undef INST

#define INST(P) template TRGSWNTT<P> TRGSW2NTT<P>(const TRGSW<P> &trgsw)
// TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
INST(lvl1param);
#undef INST

#define INST(P)                                                  \
    template TRGSW<P> trgswSymEncrypt<P>(const Polynomial<P> &p, \
                                         const Key<P> &key)
#undef INST

#define INST(P)                                                        \
    template TRGSWFFT<P> trgswfftSymEncrypt<P>(const Polynomial<P> &p, \
                                               const Key<P> &key)
#undef INST

#define INST(P)                                                        \
    template TRGSWNTT<P> trgswnttSymEncrypt<P>(const Polynomial<P> &p, \
                                               const Key<P> &key)
#undef INST
}  // namespace TFHEpp
