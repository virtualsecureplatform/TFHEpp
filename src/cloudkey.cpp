#include <cloudkey.hpp>

namespace TFHEpp {
#define INST(P) \
    template void bkgen<P>(BootstrappingKey<P> & bk, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P)                                               \
    template void bkfftgen<P>(BootstrappingKeyFFT<P> & bkfft, \
                              const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P)                                               \
    template void bknttgen<P>(BootstrappingKeyNTT<P> & bkntt, \
                              const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) \
    template void ikskgen<P>(KeySwitchingKey<P> & ksk, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P)                                                 \
    template void subikskgen<P>(SubsetKeySwitchingKey<P> & ksk, \
                                const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P)                                                              \
    template void privkskgen<P>(PrivateKeySwitchingKey<P> & ksk,             \
                                const Polynomial<typename P::targetP>& func, \
                                const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P)                                \
    template void subprivkskgen<P>(            \
        SubsetPrivateKeySwitchingKey<P> & ksk, \
        const Polynomial<typename P::targetP>& func, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) template void EvalKey::emplacebk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) template void EvalKey::emplacebkfft<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) template void EvalKey::emplacebk2bkfft<P>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) template void EvalKey::emplacebk2bkntt<P>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) template void EvalKey::emplacebkntt<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) template void EvalKey::emplaceiksk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P) template void EvalKey::emplacesubiksk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P)                                                              \
    template void EvalKey::emplaceprivksk<P>(                                \
        const std::string& key, const Polynomial<typename P::targetP>& func, \
        const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) template void EvalKey::emplaceprivksk4cb<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P)                                                              \
    template void EvalKey::emplacesubprivksk<P>(                             \
        const std::string& key, const Polynomial<typename P::targetP>& func, \
        const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) \
    template void EvalKey::emplacesubprivksk4cb<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) template BootstrappingKey<P>& EvalKey::getbk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) template BootstrappingKeyFFT<P>& EvalKey::getbkfft<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) template BootstrappingKeyNTT<P>& EvalKey::getbkntt<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) template KeySwitchingKey<P>& EvalKey::getiksk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P) \
    template SubsetKeySwitchingKey<P>& EvalKey::getsubiksk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P)                                                 \
    template PrivateKeySwitchingKey<P>& EvalKey::getprivksk<P>( \
        const std::string& key) const
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P)                                                          \
    template SubsetPrivateKeySwitchingKey<P>& EvalKey::getsubprivksk<P>( \
        const std::string& key) const
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST
}  // namespace TFHEpp