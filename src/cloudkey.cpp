#include "cloudkey.hpp"

namespace TFHEpp {

template <class P>
void bkgen(BootstrappingKey<P> &bk, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bk[i] = trgswSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P) \
    template void bkgen<P>(BootstrappingKey<P> & bk, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void bkfftgen(BootstrappingKeyFFT<P> &bkfft, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bkfft[i] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P)                                               \
    template void bkfftgen<P>(BootstrappingKeyFFT<P> & bkfft, \
                              const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void bknttgen(BootstrappingKeyNTT<P> &bkntt, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bkntt[i] = trgswnttSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P)                                               \
    template void bknttgen<P>(BootstrappingKeyNTT<P> & bkntt, \
                              const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void ikskgen(KeySwitchingKey<P> &ksk, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (uint32_t k = 0; k < (1 << P::basebit) - 1; k++)
                ksk[i][j][k] = tlweSymEncrypt<typename P::targetP>(
                    sk.key.get<typename P::domainP>()[i] * (k + 1) *
                        (1ULL
                         << (numeric_limits<typename P::targetP::T>::digits -
                             (j + 1) * P::basebit)),
                    P::α, sk.key.get<typename P::targetP>());
}
#define INST(P) \
    template void ikskgen<P>(KeySwitchingKey<P> & ksk, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

template <class P>
void privkskgen(PrivateKeySwitchingKey<P> &privksk,
                       const Polynomial<typename P::targetP>& func,
                       const SecretKey &sk)
{
    std::array<typename P::domainP::T, P::domainP::n + 1> key;
    for (int i = 0; i < P::domainP::n; i++) key[i] = sk.key.lvl2[i];
    key[P::domainP::n] = -1;
#pragma omp parallel for collapse(3)
    for (int i = 0; i <= P::domainP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (typename P::targetP::T u = 0; u < (1 << P::basebit) - 1; u++) {
                TRLWE<typename P::targetP> c =
                    trlweSymEncryptZero<typename P::targetP>(
                        P::α, sk.key.get<typename P::targetP>());
                for (int k = 0; k < P::targetP::n; k++)
                    c[1][k] +=
                        (u + 1) * func[k] * key[i]
                        << (numeric_limits<typename P::targetP::T>::digits -
                            (j + 1) * P::basebit);
                privksk[i][j][u] = c;
            }
}
#define INST(P) \
    template void privkskgen<P>(PrivateKeySwitchingKey<P> & ksk, \
        const Polynomial<typename P::targetP>& func,             \
        const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

GateKeywoFFT::GateKeywoFFT(const SecretKey &sk)
{
    // Generete bkfft
    bkgen<lvl01param>(bklvl01, sk);

    // Generete ksk
    ikskgen<lvl10param>(ksk, sk);
}

GateKey::GateKey(const SecretKey &sk)
{
    // Generete bkfft
    bkfftgen<lvl01param>(bkfftlvl01, sk);

    // Generete ksk
    ikskgen<lvl10param>(ksk, sk);
}

GateKey::GateKey(const GateKeywoFFT &gkwofft)
{
    for (int i = 0; i < lvl01param::domainP::n; i++)
        bkfftlvl01[i] = ApplyFFT2trgsw<lvl1param>(gkwofft.bklvl01[i]);
    ksk = gkwofft.ksk;
}

GateKeyNTT::GateKeyNTT(const GateKeywoFFT &gkwofft)
{
    for (int i = 0; i < lvl01param::domainP::n; i++)
        bknttlvl01[i] = ApplyNTT2trgsw<lvl1param>(gkwofft.bklvl01[i]);
    ksk = gkwofft.ksk;
}

template<class P>
void EvalKey::emplacebkfft(const SecretKey& sk){
    if constexpr (std::is_same_v<P, lvl01param>){
        bkfftlvl01 = std::make_unique<BootstrappingKeyFFT<lvl01param>>();
        bkfftgen<lvl01param>(*bkfftlvl01, sk);
    }
    else if constexpr (std::is_same_v<P, lvl02param>){
        bkfftlvl02 = std::make_unique<BootstrappingKeyFFT<lvl02param>>();
        bkfftgen<lvl02param>(*bkfftlvl02, sk);
    }
}
#define INST(P)                                     \
    template void EvalKey::emplacebkfft<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template<class P>
void EvalKey::emplaceiksk(const SecretKey& sk){
    if constexpr (std::is_same_v<P, lvl10param>){
        iksklvl10 = std::make_unique<KeySwitchingKey<lvl10param>>();
        ikskgen<lvl10param>(*iksklvl10, sk);
    }
    else if constexpr (std::is_same_v<P, lvl20param>){
        iksklvl20 = std::make_unique<KeySwitchingKey<lvl20param>>();
        ikskgen<lvl20param>(*iksklvl20, sk);
    }
    else if constexpr (std::is_same_v<P, lvl21param>){
        iksklvl21 = std::make_unique<KeySwitchingKey<lvl21param>>();
        ikskgen<lvl21param>(*iksklvl21, sk);
    }
}
#define INST(P)                                     \
    template void EvalKey::emplaceiksk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

template<class P>
void EvalKey::emplaceprivksk(const std::string &key, const Polynomial<typename P::targetP>& func, const SecretKey& sk){
    if constexpr (std::is_same_v<P, lvl21param>){
        privksklvl21[key] = std::make_unique<PrivateKeySwitchingKey<lvl21param>>();
        privkskgen<lvl21param>(*privksklvl21[key], func, sk);
    }
    else if constexpr (std::is_same_v<P, lvl22param>){
        privksklvl22[key] = std::make_unique<PrivateKeySwitchingKey<lvl22param>>();
        privkskgen<lvl22param>(*privksklvl22[key], func, sk);
    }
}
#define INST(P)                                     \
    template void EvalKey::emplaceprivksk<P>(const std::string &key, const Polynomial<typename P::targetP>& func, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

template<class P>
BootstrappingKeyFFT<P>& EvalKey::getbkfft() const{
    if constexpr (std::is_same_v<P, lvl01param>){
        return *bkfftlvl01;
    }
    else if constexpr (std::is_same_v<P, lvl02param>){
        return *bkfftlvl02;
    }
}
#define INST(P)                                     \
    template BootstrappingKeyFFT<P>& EvalKey::getbkfft<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template<class P>
KeySwitchingKey<P>& EvalKey::getiksk() const{
    if constexpr (std::is_same_v<P, lvl10param>){
        return *iksklvl10;
    }
    else if constexpr (std::is_same_v<P, lvl11param>){
        return *iksklvl11;
    }
    else if constexpr (std::is_same_v<P, lvl20param>){
        return *iksklvl20;
    }
    else if constexpr (std::is_same_v<P, lvl21param>){
        return *iksklvl21;
    }
    else if constexpr (std::is_same_v<P, lvl22param>){
        return *iksklvl22;
    }
}
#define INST(P)                                     \
    template KeySwitchingKey<P>& EvalKey::getiksk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

template<class P>
PrivateKeySwitchingKey<P>& EvalKey::getprivksk(const std::string &key) const{
    if constexpr (std::is_same_v<P, lvl11param>){
        return *(privksklvl11.at(key));
    }
    else if constexpr (std::is_same_v<P, lvl21param>){
        return *(privksklvl21.at(key));
    }
    else if constexpr (std::is_same_v<P, lvl22param>){
        return *(privksklvl22.at(key));
    }
}
#define INST(P)                                     \
    template PrivateKeySwitchingKey<P>& EvalKey::getprivksk<P>(const std::string &key) const
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

}  // namespace TFHEpp
