#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>

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
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++)
            partkey[i] = sk.key.get<P>()[0 * P::n + i];
        Automorphism<P>(autokey, partkey, (1 << (P::nbit - i)) + 1);
        ahk[i] = trgswfftSymEncrypt<P>(autokey, P::α, sk.key.get<P>());
    }
}

template <class P>
void ikskgen(KeySwitchingKey<P> &ksk, const SecretKey &sk);

template <class P>
void privkskgen(PrivateKeySwitchingKey<P> &privksk,
                const Polynomial<typename P::targetP> &func,
                const SecretKey &sk);

template <class P>
inline relinKey<P> relinKeygen(const Key<P> &key)
{
    constexpr std::array<typename P::T, P::l> h = hgen<P>();

    Polynomial<P> keysquare;
    std::array<typename P::T, P::n> partkey;
    for (int i = 0; i < P::n; i++) partkey[i] = key[0 * P::n + i];
    PolyMulNaive<P>(keysquare, partkey, partkey);
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

struct EvalKey {
    lweParams params;
    std::unique_ptr<BootstrappingKey<lvl01param>> bklvl01;
    std::unique_ptr<BootstrappingKey<lvl02param>> bklvl02;
    std::unique_ptr<BootstrappingKeyFFT<lvl01param>> bkfftlvl01;
    std::unique_ptr<BootstrappingKeyFFT<lvl02param>> bkfftlvl02;
    std::unique_ptr<BootstrappingKeyNTT<lvl01param>> bknttlvl01;
    std::unique_ptr<BootstrappingKeyNTT<lvl02param>> bknttlvl02;
    std::unique_ptr<KeySwitchingKey<lvl10param>> iksklvl10;
    std::unique_ptr<KeySwitchingKey<lvl11param>> iksklvl11;
    std::unique_ptr<KeySwitchingKey<lvl20param>> iksklvl20;
    std::unique_ptr<KeySwitchingKey<lvl21param>> iksklvl21;
    std::unique_ptr<KeySwitchingKey<lvl22param>> iksklvl22;
    std::unique_ptr<SubsetKeySwitchingKey<lvl21param>> subiksklvl21;
    std::unordered_map<std::string,
                       std::unique_ptr<PrivateKeySwitchingKey<lvl11param>>>
        privksklvl11;
    std::unordered_map<std::string,
                       std::unique_ptr<PrivateKeySwitchingKey<lvl21param>>>
        privksklvl21;
    std::unordered_map<std::string,
                       std::unique_ptr<PrivateKeySwitchingKey<lvl22param>>>
        privksklvl22;
    std::unordered_map<std::string,
                       std::unique_ptr<SubsetPrivateKeySwitchingKey<lvl21param>>>
        subprivksklvl21;

    EvalKey(SecretKey sk) { params = sk.params; }
    EvalKey() {}

    template <class P>
    void emplacebk(const SecretKey &sk);
    template <class P>
    void emplacebkfft(const SecretKey &sk);
    template <class P>
    void emplacebkntt(const SecretKey &sk);
    template <class P>
    void emplacebk2bkfft();
    template <class P>
    void emplacebk2bkntt();
    template <class P>
    void emplaceiksk(const SecretKey &sk);
    template <class P>
    void emplacesubiksk(const SecretKey &sk);
    template <class P>
    void emplaceprivksk(const std::string &key,
                        const Polynomial<typename P::targetP> &func,
                        const SecretKey &sk);
    template <class P>
    void emplacesubprivksk(const std::string &key,
                        const Polynomial<typename P::targetP> &func,
                        const SecretKey &sk);
    template <class P>
    void emplaceprivksk4cb(const SecretKey &sk);
    template <class P>
    void emplacesubprivksk4cb(const SecretKey &sk);

    template <class P>
    BootstrappingKey<P> &getbk() const;
    template <class P>
    BootstrappingKeyFFT<P> &getbkfft() const;
    template <class P>
    BootstrappingKeyNTT<P> &getbkntt() const;
    template <class P>
    KeySwitchingKey<P> &getiksk() const;
    template <class P>
    SubsetKeySwitchingKey<P> &getsubiksk() const;
    template <class P>
    PrivateKeySwitchingKey<P> &getprivksk(const std::string &key) const;
    template <class P>
    SubsetPrivateKeySwitchingKey<P> &getsubprivksk(const std::string &key) const;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(params, bklvl01, bklvl02, bkfftlvl01, bkfftlvl02, bknttlvl01,
                bknttlvl02, iksklvl10, iksklvl11, iksklvl20, iksklvl21,
                iksklvl22, privksklvl11, privksklvl21, privksklvl22);
    }
};
}  // namespace TFHEpp