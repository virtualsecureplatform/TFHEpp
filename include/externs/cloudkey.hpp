#pragma once
#include"../cloudkey.hpp"

namespace TFHEpp{
#define INST(P) \
    extern template void bkgen<P>(BootstrappingKey<P> & bk, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P)                                               \
    extern template void bkfftgen<P>(BootstrappingKeyFFT<P> & bkfft, \
                              const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P)                                               \
    extern template void bknttgen<P>(BootstrappingKeyNTT<P> & bkntt, \
                              const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) \
    extern template void ikskgen<P>(KeySwitchingKey<P> & ksk, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P) \
    extern template void subikskgen<P>(SubsetKeySwitchingKey<P> & ksk, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P)                                                              \
    extern template void privkskgen<P>(PrivateKeySwitchingKey<P> & ksk,             \
                                const Polynomial<typename P::targetP>& func, \
                                const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P)                                                              \
    extern template void subprivkskgen<P>(SubsetPrivateKeySwitchingKey<P> & ksk,             \
                                const Polynomial<typename P::targetP>& func, \
                                const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplacebk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplacebkfft<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplacebk2bkfft<P>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplacebk2bkntt<P>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplacebkntt<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplaceiksk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplacesubiksk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P)                                                              \
    extern template void EvalKey::emplaceprivksk<P>(                                \
        const std::string& key, const Polynomial<typename P::targetP>& func, \
        const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplaceprivksk4cb<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplaceahk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_ANNIHILATE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplacecbsk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_ANNIHILATE(INST)
#undef INST

#define INST(P)                                                              \
    extern template void EvalKey::emplacesubprivksk<P>(                                \
        const std::string& key, const Polynomial<typename P::targetP>& func, \
        const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) extern template void EvalKey::emplacesubprivksk4cb<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) extern template BootstrappingKey<P>& EvalKey::getbk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template BootstrappingKeyFFT<P>& EvalKey::getbkfft<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template BootstrappingKeyNTT<P>& EvalKey::getbkntt<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template KeySwitchingKey<P>& EvalKey::getiksk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P) extern template SubsetKeySwitchingKey<P>& EvalKey::getsubiksk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P)                                                 \
    extern template PrivateKeySwitchingKey<P>& EvalKey::getprivksk<P>( \
        const std::string& key) const
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P)                                                 \
    extern template SubsetPrivateKeySwitchingKey<P>& EvalKey::getsubprivksk<P>( \
        const std::string& key) const
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) extern template AnnihilateKey<P>& EvalKey::getahk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_ANNIHILATE(INST)
#undef INST

#define INST(P) extern template CBswitchingKey<P>& EvalKey::getcbsk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_ANNIHILATE(INST)
#undef INST

// Extern template definitions for get<T>() helper functions
#define INST(P) extern template auto& EvalKey::get<BootstrappingKey<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template auto& EvalKey::get<BootstrappingKeyFFT<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template auto& EvalKey::get<BootstrappingKeyNTT<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(P) extern template auto& EvalKey::get<KeySwitchingKey<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P) extern template auto& EvalKey::get<SubsetKeySwitchingKey<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TLWE(INST)
#undef INST

#define INST(P) extern template auto& EvalKey::get<AnnihilateKey<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_ANNIHILATE(INST)
#undef INST

#define INST(P) extern template auto& EvalKey::get<CBswitchingKey<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_ANNIHILATE(INST)
#undef INST

// Extern template definitions for get_map<T>() helper functions
#define INST(P) extern template auto& EvalKey::get_map<PrivateKeySwitchingKey<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

#define INST(P) extern template auto& EvalKey::get_map<SubsetPrivateKeySwitchingKey<P>>()
TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

}