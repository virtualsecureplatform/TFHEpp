#import <cloudkey.hpp>

namespace TFHEpp {
#define INST(P) \
    template void bkgen<P>(BootstrappingKey<P> & bk, const SecretKey& sk)
#undef INST

#define INST(P)                                               \
    template void bkfftgen<P>(BootstrappingKeyFFT<P> & bkfft, \
                              const SecretKey& sk)
#undef INST

#define INST(P)                                               \
    template void bknttgen<P>(BootstrappingKeyNTT<P> & bkntt, \
                              const SecretKey& sk)
#undef INST

#define INST(P) \
    template void ikskgen<P>(KeySwitchingKey<P> & ksk, const SecretKey& sk)
#undef INST

#define INST(P)                                                 \
    template void subikskgen<P>(SubsetKeySwitchingKey<P> & ksk, \
                                const SecretKey& sk)
#undef INST

#define INST(P)                                                              \
    template void privkskgen<P>(PrivateKeySwitchingKey<P> & ksk,             \
                                const Polynomial<typename P::targetP>& func, \
                                const SecretKey& sk)
#undef INST

#define INST(P)                                \
    template void subprivkskgen<P>(            \
        SubsetPrivateKeySwitchingKey<P> & ksk, \
        const Polynomial<typename P::targetP>& func, const SecretKey& sk)
#undef INST

#define INST(P) template void EvalKey::emplacebk<P>(const SecretKey& sk)
#undef INST

#define INST(P) template void EvalKey::emplacebkfft<P>(const SecretKey& sk)
#undef INST

#define INST(P) template void EvalKey::emplacebk2bkfft<P>()
#undef INST

#define INST(P) template void EvalKey::emplacebk2bkntt<P>()
#undef INST

#define INST(P) template void EvalKey::emplacebkntt<P>(const SecretKey& sk)
#undef INST

#define INST(P) template void EvalKey::emplaceiksk<P>(const SecretKey& sk)
#undef INST

#define INST(P) template void EvalKey::emplacesubiksk<P>(const SecretKey& sk)
#undef INST

#define INST(P)                                                              \
    template void EvalKey::emplaceprivksk<P>(                                \
        const std::string& key, const Polynomial<typename P::targetP>& func, \
        const SecretKey& sk)
#undef INST

#define INST(P) template void EvalKey::emplaceprivksk4cb<P>(const SecretKey& sk)
#undef INST

#define INST(P)                                                              \
    template void EvalKey::emplacesubprivksk<P>(                             \
        const std::string& key, const Polynomial<typename P::targetP>& func, \
        const SecretKey& sk)
#undef INST

#define INST(P) \
    template void EvalKey::emplacesubprivksk4cb<P>(const SecretKey& sk)
#undef INST

#define INST(P) template BootstrappingKey<P>& EvalKey::getbk<P>() const
#undef INST

#define INST(P) template BootstrappingKeyFFT<P>& EvalKey::getbkfft<P>() const
#undef INST

#define INST(P) template BootstrappingKeyNTT<P>& EvalKey::getbkntt<P>() const
#undef INST

#define INST(P) template KeySwitchingKey<P>& EvalKey::getiksk<P>() const
#undef INST

#define INST(P) \
    template SubsetKeySwitchingKey<P>& EvalKey::getsubiksk<P>() const
#undef INST

#define INST(P)                                                 \
    template PrivateKeySwitchingKey<P>& EvalKey::getprivksk<P>( \
        const std::string& key) const
#undef INST

#define INST(P)                                                          \
    template SubsetPrivateKeySwitchingKey<P>& EvalKey::getsubprivksk<P>( \
        const std::string& key) const
#undef INST
}  // namespace TFHEpp