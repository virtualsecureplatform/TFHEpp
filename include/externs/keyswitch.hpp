#pragma once
#import "../keyswitch.hpp"

namespace TFHEpp{
#define INST(P)                                                               \
                                       const TLWE<typename P::domainP> &tlwe, \
                                       const KeySwitchingKey<P> &ksk)
#undef INST

#define INST(P)                                                               \
                                       const TLWE<typename P::domainP> &tlwe, \
                                       const SubsetKeySwitchingKey<P> &ksk)
#undef INST

#define INST(P)                                                           \
                                   const TLWE<typename P::domainP> &tlwe, \
                                   const TLWE2TRLWEIKSKey<P> &iksk)
#undef INST

#define INST(P)                                                      \
                              const int d, const EvalAutoKey<P> &autokey)
#undef INST

#define INST(P)                              \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const AnnihilateKey<P> &ahk)
#undef INST

#define INST(P)                                                           \
                                   const TLWE<typename P::domainP> &tlwe, \
                                   const PrivateKeySwitchingKey<P> &privksk)
#undef INST

#define INST(P)                                                           \
                                   const TLWE<typename P::targetP> &tlwe, \
                                   const SubsetPrivateKeySwitchingKey<P> &privksk)
#undef INST
}