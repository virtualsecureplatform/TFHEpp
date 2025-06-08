#import <keyswitch.hpp>

namespace TFHEpp {

#define INST(P)                                                               \
    template void IdentityKeySwitch<P>(TLWE<typename P::targetP> & res,       \
                                       const TLWE<typename P::domainP> &tlwe, \
                                       const KeySwitchingKey<P> &ksk)
#undef INST

#define INST(P)                                \
    template void SubsetIdentityKeySwitch<P>(  \
        TLWE<typename P::targetP> & res,       \
        const TLWE<typename P::domainP> &tlwe, \
        const SubsetKeySwitchingKey<P> &ksk)
#undef INST

#define INST(P)                                                           \
    template void TLWE2TRLWEIKS<P>(TRLWE<typename P::targetP> & res,      \
                                   const TLWE<typename P::domainP> &tlwe, \
                                   const TLWE2TRLWEIKSKey<P> &iksk)
#undef INST

#define INST(P)                                                      \
    template void EvalAuto<P>(TRLWE<P> & res, const TRLWE<P> &trlwe, \
                              const int d, const EvalAutoKey<P> &autokey)
#undef INST

#define INST(P)                              \
    template void AnnihilateKeySwitching<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const AnnihilateKey<P> &ahk)
#undef INST

#define INST(P)                                                           \
    template void PrivKeySwitch<P>(TRLWE<typename P::targetP> & res,      \
                                   const TLWE<typename P::domainP> &tlwe, \
                                   const PrivateKeySwitchingKey<P> &privksk)
#undef INST

#define INST(P)                                \
    template void SubsetPrivKeySwitch<P>(      \
        TRLWE<typename P::targetP> & res,      \
        const TLWE<typename P::targetP> &tlwe, \
        const SubsetPrivateKeySwitchingKey<P> &privksk)
#undef INST

}  // namespace TFHEpp
