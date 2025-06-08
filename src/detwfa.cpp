#import <detwfa.hpp>

namespace TFHEpp {
#define INST(P)                                                     \
    template void CMUXFFT<P>(TRLWE<P> & res, const TRGSWFFT<P> &cs, \
                             const TRLWE<P> &c1, const TRLWE<P> &c0)
#undef INST

#define INST(bkP)                                             \
    template void CMUXFFTwithPolynomialMulByXaiMinusOne<bkP>( \
        TRLWE<typename bkP::targetP> & acc,                   \
        const BootstrappingKeyElementFFT<bkP> &cs, const int a)
#undef INST

#define INST(P)                                             \
    template void CMUXNTTwithPolynomialMulByXaiMinusOne<P>( \
        TRLWE<P> & acc, const TRGSWNTT<P> &cs, const typename P::T a)
#undef INST
}  // namespace TFHEpp