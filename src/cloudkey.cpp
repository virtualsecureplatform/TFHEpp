#include "cloudkey.hpp"

#include <string>

namespace TFHEpp {

template <class P>
void bkgen(BootstrappingKey<P>& bk, const SecretKey& sk)
{
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bk[i] = trgswSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P) \
    template void bkgen<P>(BootstrappingKey<P> & bk, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void bkfftgen(BootstrappingKeyFFT<P>& bkfft, const SecretKey& sk)
{
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bkfft[i] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P)                                               \
    template void bkfftgen<P>(BootstrappingKeyFFT<P> & bkfft, \
                              const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void bknttgen(BootstrappingKeyNTT<P>& bkntt, const SecretKey& sk)
{
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bkntt[i] = trgswnttSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P)                                               \
    template void bknttgen<P>(BootstrappingKeyNTT<P> & bkntt, \
                              const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void ikskgen(KeySwitchingKey<P>& ksk, const SecretKey& sk)
{
    for (int l = 0; l < P::domainP::k; l++)
        for (int i = 0; i < P::domainP::n; i++)
            for (int j = 0; j < P::t; j++)
                for (uint32_t k = 0; k < (1 << P::basebit) - 1; k++)
                    ksk[l * P::domainP::n + i][j][k] =
                        tlweSymEncrypt<typename P::targetP>(
                            sk.key.get<
                                typename P::domainP>()[l * P::domainP::n + i] *
                                (k + 1) *
                                (1ULL << (numeric_limits<
                                              typename P::targetP::T>::digits -
                                          (j + 1) * P::basebit)),
                            P::α, sk.key.get<typename P::targetP>());
}
#define INST(P) \
    template void ikskgen<P>(KeySwitchingKey<P> & ksk, const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

template <class P>
void privkskgen(PrivateKeySwitchingKey<P>& privksk,
                const Polynomial<typename P::targetP>& func,
                const SecretKey& sk)
{
    std::array<typename P::domainP::T, P::domainP::n + 1> key;
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++)
        key[i] = sk.key.get<typename P::domainP>()[i];
    key[P::domainP::k * P::domainP::n] = -1;
#pragma omp parallel for collapse(3)
    for (int i = 0; i <= P::domainP::k * P::domainP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (typename P::targetP::T u = 0; u < (1 << P::basebit) - 1; u++) {
                TRLWE<typename P::targetP> c =
                    trlweSymEncryptZero<typename P::targetP>(
                        P::α, sk.key.get<typename P::targetP>());
                for (int k = 0; k < P::targetP::n; k++)
                    c[P::targetP::k][k] +=
                        (u + 1) * func[k] * key[i]
                        << (numeric_limits<typename P::targetP::T>::digits -
                            (j + 1) * P::basebit);
                privksk[i][j][u] = c;
            }
}
#define INST(P)                                                              \
    template void privkskgen<P>(PrivateKeySwitchingKey<P> & ksk,             \
                                const Polynomial<typename P::targetP>& func, \
                                const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

template <class P>
void EvalKey::emplacebk(const SecretKey& sk)
{
    if constexpr (std::is_same_v<P, lvl01param>) {
        bklvl01 = std::make_unique<BootstrappingKey<lvl01param>>();
        bkgen<lvl01param>(*bklvl01, sk);
    }
    else if constexpr (std::is_same_v<P, lvl02param>) {
        bklvl02 = std::make_unique<BootstrappingKey<lvl02param>>();
        bkgen<lvl02param>(*bklvl02, sk);
    }
}
#define INST(P) template void EvalKey::emplacebk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void EvalKey::emplacebkfft(const SecretKey& sk)
{
    if constexpr (std::is_same_v<P, lvl01param>) {
        bkfftlvl01 = std::make_unique<BootstrappingKeyFFT<lvl01param>>();
        bkfftgen<lvl01param>(*bkfftlvl01, sk);
    }
    else if constexpr (std::is_same_v<P, lvl02param>) {
        bkfftlvl02 = std::make_unique<BootstrappingKeyFFT<lvl02param>>();
        bkfftgen<lvl02param>(*bkfftlvl02, sk);
    }
}
#define INST(P) template void EvalKey::emplacebkfft<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void EvalKey::emplacebk2bkfft()
{
    if constexpr (std::is_same_v<P, lvl01param>) {
        bkfftlvl01 = std::make_unique<BootstrappingKeyFFT<lvl01param>>();
        for (int i = 0; i < lvl01param::domainP::n; i++)
            (*bkfftlvl01)[i] = ApplyFFT2trgsw<lvl1param>((*bklvl01)[i]);
    }
    else if constexpr (std::is_same_v<P, lvl02param>) {
        bkfftlvl02 = std::make_unique<BootstrappingKeyFFT<lvl02param>>();
        for (int i = 0; i < lvl02param::domainP::n; i++)
            (*bkfftlvl02)[i] = ApplyFFT2trgsw<lvl2param>((*bklvl02)[i]);
    }
}
#define INST(P) template void EvalKey::emplacebk2bkfft<P>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void EvalKey::emplacebk2bkntt()
{
    if constexpr (std::is_same_v<P, lvl01param>) {
        bknttlvl01 = std::make_unique<BootstrappingKeyNTT<lvl01param>>();
        for (int i = 0; i < lvl01param::domainP::n; i++)
            (*bknttlvl01)[i] = ApplyNTT2trgsw<lvl1param>((*bklvl01)[i]);
    }
    else if constexpr (std::is_same_v<P, lvl02param>) {
        bknttlvl02 = std::make_unique<BootstrappingKeyNTT<lvl02param>>();
        for (int i = 0; i < lvl02param::domainP::n; i++)
            (*bknttlvl02)[i] = ApplyNTT2trgsw<lvl2param>((*bklvl02)[i]);
    }
}
#define INST(P) template void EvalKey::emplacebk2bkntt<P>()
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void EvalKey::emplacebkntt(const SecretKey& sk)
{
    if constexpr (std::is_same_v<P, lvl01param>) {
        bknttlvl01 = std::make_unique<BootstrappingKeyNTT<lvl01param>>();
        bknttgen<lvl01param>(*bknttlvl01, sk);
    }
    else if constexpr (std::is_same_v<P, lvl02param>) {
        bknttlvl02 = std::make_unique<BootstrappingKeyNTT<lvl02param>>();
        bknttgen<lvl02param>(*bknttlvl02, sk);
    }
}
#define INST(P) template void EvalKey::emplacebkntt<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void EvalKey::emplaceiksk(const SecretKey& sk)
{
    if constexpr (std::is_same_v<P, lvl10param>) {
        iksklvl10 = std::make_unique<KeySwitchingKey<lvl10param>>();
        ikskgen<lvl10param>(*iksklvl10, sk);
    }
    else if constexpr (std::is_same_v<P, lvl20param>) {
        iksklvl20 = std::make_unique<KeySwitchingKey<lvl20param>>();
        ikskgen<lvl20param>(*iksklvl20, sk);
    }
    else if constexpr (std::is_same_v<P, lvl21param>) {
        iksklvl21 = std::make_unique<KeySwitchingKey<lvl21param>>();
        ikskgen<lvl21param>(*iksklvl21, sk);
    }
}
#define INST(P) template void EvalKey::emplaceiksk<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

template <class P>
void EvalKey::emplaceprivksk(const std::string& key,
                             const Polynomial<typename P::targetP>& func,
                             const SecretKey& sk)
{
    if constexpr (std::is_same_v<P, lvl21param>) {
        privksklvl21[key] =
            std::make_unique<PrivateKeySwitchingKey<lvl21param>>();
        privkskgen<lvl21param>(*privksklvl21[key], func, sk);
    }
    else if constexpr (std::is_same_v<P, lvl22param>) {
        privksklvl22[key] =
            std::make_unique<PrivateKeySwitchingKey<lvl22param>>();
        privkskgen<lvl22param>(*privksklvl22[key], func, sk);
    }
}
#define INST(P)                                                              \
    template void EvalKey::emplaceprivksk<P>(                                \
        const std::string& key, const Polynomial<typename P::targetP>& func, \
        const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

template <class P>
void EvalKey::emplaceprivksk4cb(const SecretKey& sk)
{
    for (int k = 0; k < P::targetP::k; k++) {
        Polynomial<typename P::targetP> partkey;
        for (int i = 0; i < P::targetP::n; i++)
            partkey[i] =
                -sk.key.get<typename P::targetP>()[k * P::targetP::n + i];
        emplaceprivksk<P>("privksk4cb_" + std::to_string(k), partkey, sk);
    }
    emplaceprivksk<P>("privksk4cb_" + std::to_string(P::targetP::k), {1}, sk);
}
#define INST(P) template void EvalKey::emplaceprivksk4cb<P>(const SecretKey& sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

template <class P>
BootstrappingKey<P>& EvalKey::getbk() const
{
    if constexpr (std::is_same_v<P, lvl01param>) {
        return *bklvl01;
    }
    else if constexpr (std::is_same_v<P, lvl02param>) {
        return *bklvl02;
    }
}
#define INST(P) template BootstrappingKey<P>& EvalKey::getbk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
BootstrappingKeyFFT<P>& EvalKey::getbkfft() const
{
    if constexpr (std::is_same_v<P, lvl01param>) {
        return *bkfftlvl01;
    }
    else if constexpr (std::is_same_v<P, lvl02param>) {
        return *bkfftlvl02;
    }
}
#define INST(P) template BootstrappingKeyFFT<P>& EvalKey::getbkfft<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
BootstrappingKeyNTT<P>& EvalKey::getbkntt() const
{
    if constexpr (std::is_same_v<P, lvl01param>) {
        return *bknttlvl01;
    }
    else if constexpr (std::is_same_v<P, lvl02param>) {
        return *bknttlvl02;
    }
}
#define INST(P) template BootstrappingKeyNTT<P>& EvalKey::getbkntt<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
KeySwitchingKey<P>& EvalKey::getiksk() const
{
    if constexpr (std::is_same_v<P, lvl10param>) {
        return *iksklvl10;
    }
    else if constexpr (std::is_same_v<P, lvl11param>) {
        return *iksklvl11;
    }
    else if constexpr (std::is_same_v<P, lvl20param>) {
        return *iksklvl20;
    }
    else if constexpr (std::is_same_v<P, lvl21param>) {
        return *iksklvl21;
    }
    else if constexpr (std::is_same_v<P, lvl22param>) {
        return *iksklvl22;
    }
}
#define INST(P) template KeySwitchingKey<P>& EvalKey::getiksk<P>() const
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

template <class P>
PrivateKeySwitchingKey<P>& EvalKey::getprivksk(const std::string& key) const
{
    if constexpr (std::is_same_v<P, lvl11param>) {
        return *(privksklvl11.at(key));
    }
    else if constexpr (std::is_same_v<P, lvl21param>) {
        return *(privksklvl21.at(key));
    }
    else if constexpr (std::is_same_v<P, lvl22param>) {
        return *(privksklvl22.at(key));
    }
}
#define INST(P)                                                 \
    template PrivateKeySwitchingKey<P>& EvalKey::getprivksk<P>( \
        const std::string& key) const
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

}  // namespace TFHEpp
