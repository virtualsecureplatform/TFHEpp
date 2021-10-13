#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>

#include "key.hpp"
#include "params.hpp"
#include "tlwe.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class P>
void bkgen(BootstrappingKey<P> &bk, const SecretKey &sk);

template <class P>
void bkfftgen(BootstrappingKeyFFT<P> &bkfft, const SecretKey &sk);

template <class P>
void bknttgen(BootstrappingKeyNTT<P> &bkntt, const SecretKey &sk);

template <class P>
void tlwe2trlweikskkgen(TLWE2TRLWEIKSKey<P> &iksk, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (uint32_t k = 0; k < (1 << P::basebit) - 1; k++) {
                Polynomial<typename P::targetP> p = {};
                p[0] =
                    sk.key.get<typename P::domainP>()[i] * (k + 1) *
                    (1ULL << (numeric_limits<typename P::targetP::T>::digits -
                              (j + 1) * P::basebit));
                iksk[i][j][k] = trlweSymEncrypt<typename P::targetP>(
                    p, P::α, sk.key.get<typename P::targetP>());
            }
}

template <class P>
inline void annihilatekeyegen(AnnihilateKey<P> &ahk, const SecretKey &sk)
{
    for (int i = 0; i < P::nbit; i++) {
        Polynomial<P> autokey;
        Automorphism<P>(autokey, sk.key.get<P>(), (1 << (P::nbit - i)) + 1);
        ahk[i] = trgswfftSymEncrypt<P>(autokey, P::α, sk.key.get<P>());
    }
}

template <class P>
void ikskgen(KeySwitchingKey<P> &ksk, const SecretKey &sk);

template <class P>
inline void privkskgen(PrivateKeySwitchingKey<P> &privksk,
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

template <class P>
inline relinKey<P> relinKeygen(const Key<P> &key)
{
    constexpr std::array<typename P::T, P::l> h = hgen<P>();

    Polynomial<P> keysquare;
    PolyMulNaieve<P>(keysquare, key, key);
    relinKey<P> relinkey;
    for (TRLWE<P> &ctxt : relinkey) ctxt = trlweSymEncryptZero<P>(P::α, key);
    for (int i = 0; i < P::l; i++)
        for (int j = 0; j < P::n; j++)
            relinkey[i][1][j] +=
                static_cast<typename P::T>(keysquare[j]) * h[i];
    return relinkey;
}

template <class P>
inline relinKeyFFT<P> relinKeyFFTgen(const Key<P> &key)
{
    relinKey<P> relinkey = relinKeygen<P>(key);
    relinKeyFFT<P> relinkeyfft;
    for (int i = 0; i < P::l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFT<P>(relinkeyfft[i][j], relinkey[i][j]);
    return relinkeyfft;
}

struct GateKeywoFFT {
    BootstrappingKey<lvl01param> bklvl01;
    KeySwitchingKey<lvl10param> ksk;
    GateKeywoFFT(const SecretKey &sk);
    GateKeywoFFT() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(ksk, bklvl01);
    }
};

struct GateKey {
    BootstrappingKeyFFT<lvl01param> bkfftlvl01;
    KeySwitchingKey<lvl10param> ksk;
    GateKey(const SecretKey &sk);
    GateKey(const GateKeywoFFT &gkwofft);
    GateKey() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(ksk, bkfftlvl01);
    }
};

struct GateKeyNTT {
    BootstrappingKeyNTT<lvl01param> bknttlvl01;
    KeySwitchingKey<lvl10param> ksk;
    // GateKey(const SecretKey &sk);
    GateKeyNTT(const GateKeywoFFT &gkwofft);
    GateKeyNTT() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(ksk, bknttlvl01);
    }
};

template <class bsP, class privksP>
struct CircuitKey {
    BootstrappingKeyFFT<bsP> bkfft;
    std::array<PrivateKeySwitchingKey<privksP>, 2> privksk;
    CircuitKey(const SecretKey &sk);
    CircuitKey() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(privksk, bkfft);
    }
};

template <class CBbsP, class CBprivksP, class CMksP>
struct CloudKey {
    GateKey gk;
    CircuitKey<CBbsP, CBprivksP> ck;
    KeySwitchingKey<CMksP> ksk;
    lweParams params;
    CloudKey(SecretKey sk) : gk(sk), ck(sk) { ikskgen<CMksP>(ksk, sk); params = sk.params;}
    CloudKey() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(gk.ksk, gk.bkfftlvl01, ck.privksk, ck.bkfft, params);
    }
};

struct EvalKey {
    lweParams params;
    std::unique_ptr<BootstrappingKeyFFT<lvl01param>> bkfftlvl01;
    std::unique_ptr<BootstrappingKeyFFT<lvl02param>> bkfftlvl02;
    std::unique_ptr<KeySwitchingKey<lvl10param>> iksklvl10;
    std::unique_ptr<KeySwitchingKey<lvl11param>> iksklvl11;
    std::unique_ptr<KeySwitchingKey<lvl20param>> iksklvl20;
    std::unique_ptr<KeySwitchingKey<lvl21param>> iksklvl21;
    std::unique_ptr<KeySwitchingKey<lvl22param>> iksklvl22;
    std::unordered_map<std::string,std::unique_ptr<PrivateKeySwitchingKey<lvl10param>>> privksklvl10;
    std::unordered_map<std::string,std::unique_ptr<PrivateKeySwitchingKey<lvl11param>>> privksklvl11;
    std::unordered_map<std::string,std::unique_ptr<PrivateKeySwitchingKey<lvl20param>>> privksklvl20;
    std::unordered_map<std::string,std::unique_ptr<PrivateKeySwitchingKey<lvl21param>>> privksklvl21;
    std::unordered_map<std::string,std::unique_ptr<PrivateKeySwitchingKey<lvl22param>>> privksklvl22;

    EvalKey(SecretKey sk) {params = sk.params;}
    EvalKey() {}

    template<class P> void emplacebkfft(const SecretKey &sk);
    template<class P> void emplaceiksk(const SecretKey &sk);
    template<class P> void emplaceprivksk(const std::string &key, const Polynomial<typename P::targetP>& func, const SecretKey &sk);
    template<class P, uint index> void emplaceprivksk(const SecretKey &sk){
        if constexpr(index == 0){
            emplaceprivksk<P>("identity",{1},sk);
        }else if constexpr(index == 1){
            TFHEpp::Polynomial<typename P::targetP> poly;
            for (int i = 0; i < P::targetP::n; i++)
                poly[i] = -sk.key.get<typename P::targetP>()[i];
            emplaceprivksk<P>("secret key",poly,sk);
        }else{
            static_assert(false_v<P>, "Not a predefined function for Private Key Switching!");
        }
    }

    template<class P> BootstrappingKeyFFT<P>& getbkfft() const;
    template<class P> KeySwitchingKey<P>& getiksk() const;
    template<class P> PrivateKeySwitchingKey<P>& getprivksk(const std::string &key) const;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(params,bkfftlvl01,bkfftlvl02,iksklvl10,iksklvl20,iksklvl21,privksklvl21,privksklvl22);
    }
};
}  // namespace TFHEpp