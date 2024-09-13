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
void bkgen(BootstrappingKey<P>& bk, const Key<typename P::domainP>& domainkey,
           const Key<typename P::targetP>& targetkey)
{
    Polynomial<typename P::targetP> plainpoly = {};
#ifdef USE_KEY_BUNDLE
    for (int i = 0; i < P::domainP::k * P::domainP::n / P::Addends; i++) {
        plainpoly[0] =
            static_cast<int32_t>(domainkey[2 * i] * domainkey[2 * i + 1]);
        bk[i][0] = trgswSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, targetkey);
        plainpoly[0] =
            static_cast<int32_t>(domainkey[2 * i] * (1 - domainkey[2 * i + 1]));
        bk[i][1] = trgswSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, targetkey);
        plainpoly[0] =
            static_cast<int32_t>((1 - domainkey[2 * i]) * domainkey[2 * i + 1]);
        bk[i][2] = trgswSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, targetkey);
    }
#else
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        int count = 0;
        for (int j = P::domainP::key_value_min; j <= P::domainP::key_value_max;
             j++) {
            if (j != 0) {
                plainpoly[0] = domainkey[i] == j;
                bk[i][count] =
                    trgswSymEncrypt<typename P::targetP>(plainpoly, targetkey);
                count++;
            }
        }
    }
#endif
}

template <class P>
void bkgen(BootstrappingKey<P>& bk, const SecretKey& sk)
{
    bkgen<P>(bk, sk.key.get<typename P::domainP>(),
             sk.key.get<typename P::targetP>());
}

template <class P>
void bkfftgen(BootstrappingKeyFFT<P>& bkfft,
              const Key<typename P::domainP>& domainkey,
              const Key<typename P::targetP>& targetkey)
{
    Polynomial<typename P::targetP> plainpoly = {};
#ifdef USE_KEY_BUNDLE
    for (int i = 0; i < P::domainP::k * P::domainP::n / P::Addends; i++) {
        plainpoly[0] =
            static_cast<int32_t>(domainkey[2 * i] * domainkey[2 * i + 1]);
        bkfft[i][0] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, targetkey);
        plainpoly[0] =
            static_cast<int32_t>(domainkey[2 * i] * (1 - domainkey[2 * i + 1]));
        bkfft[i][1] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, targetkey);
        plainpoly[0] =
            static_cast<int32_t>((1 - domainkey[2 * i]) * domainkey[2 * i + 1]);
        bkfft[i][2] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, targetkey);
    }
#else
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        int count = 0;
        for (int j = P::domainP::key_value_min; j <= P::domainP::key_value_max;
             j++) {
            if (j != 0) {
                plainpoly[0] = domainkey[i] == j;
                bkfft[i][count] = trgswfftSymEncrypt<typename P::targetP>(
                    plainpoly, targetkey);
                count++;
            }
        }
    }
#endif
}

template <class P>
void bkfftgen(BootstrappingKeyFFT<P>& bkfft, const SecretKey& sk)
{
    bkfftgen<P>(bkfft, sk.key.get<typename P::domainP>(),
                sk.key.get<typename P::targetP>());
}

template <class P>
void bknttgen(BootstrappingKeyNTT<P>& bkntt,
              const Key<typename P::domainP>& domainkey,
              const Key<typename P::targetP>& targetkey)
{
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = domainkey[i];
        bkntt[i] =
            trgswnttSymEncrypt<typename P::targetP>(plainpoly, targetkey);
    }
}

template <class P>
void bknttgen(BootstrappingKeyNTT<P>& bkntt, const SecretKey& sk)
{
    bknttgen<P>(bkntt, sk.key.get<typename P::domainP>(),
                sk.key.get<typename P::targetP>());
}

template <class P>
void bkrainttgen(BootstrappingKeyRAINTT<P>& bkraintt,
                 const Key<typename P::domainP>& domainkey,
                 const Key<typename P::targetP>& targetkey)
{
#pragma omp parallel for
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = domainkey[i];
        bkraintt[i] =
            trgswrainttSymEncrypt<typename P::targetP>(plainpoly, targetkey);
    }
}

template <class P>
void bkrainttgen(BootstrappingKeyRAINTT<P>& bkraintt, const SecretKey& sk)
{
    bkrainttgen<P>(bkraintt, sk.key.get<typename P::domainP>(),
                   sk.key.get<typename P::targetP>());
}

template <class P>
void tlwe2trlweikskgen(TLWE2TRLWEIKSKey<P>& iksk,
                       const Key<typename P::domainP>& domainkey,
                       const Key<typename P::targetP>& targetkey)
{
    for (int i = 0; i < P::domainP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (uint32_t k = 0; k < (1 << P::basebit) - 1; k++) {
                Polynomial<typename P::targetP> p = {};
                p[0] =
                    domainkey[i] * (k + 1) *
                    (1ULL << (numeric_limits<typename P::targetP::T>::digits -
                              (j + 1) * P::basebit));
                iksk[i][j][k] =
                    trlweSymEncrypt<typename P::targetP>(p, targetkey);
            }
}

template <class P>
void tlwe2trlweikskgen(TLWE2TRLWEIKSKey<P>& iksk, const SecretKey& sk)
{
    tlwe2trlweikskgen<P>(iksk, sk.key.get<typename P::domainP>(),
                         sk.key.get<typename P::targetP>());
}

template <class P>
void evalautokeygen(EvalAutoKey<P>& eak, const uint d, const Key<P>& key)
{
    for(int j = 0; j < P::k; j++) {
        Polynomial<P> autokey;
        std::array<typename P::T, P::n> partkey;
        for (int k = 0; k < P::n; k++) partkey[k] = key[j * P::n + k];
        Automorphism<P>(autokey, partkey, d);
        eak[j] = halftrgswfftSymEncrypt<P>(autokey, key);
    }
}

template <class P>
void annihilatekeygen(AnnihilateKey<P>& ahk, const Key<P>& key)
{
    for (int i = 0; i < P::nbit; i++)
        evalautokeygen<P>(ahk[i], (1 << (P::nbit - i)) + 1, key);
}

template <class P>
void annihilatekeygen(AnnihilateKey<P>& ahk, const SecretKey& sk)
{
    annihilatekeygen<P>(ahk, sk.key.get<P>());
}

template <class P>
void ikskgen(KeySwitchingKey<P>& ksk, const Key<typename P::domainP>& domainkey,
             const Key<typename P::targetP>& targetkey)
{
    for (int l = 0; l < P::domainP::k; l++)
        for (int i = 0; i < P::domainP::n; i++)
            for (int j = 0; j < P::t; j++)
                for (uint32_t k = 0; k < 1U << (P::basebit - 1); k++)
                    ksk[l * P::domainP::n + i][j][k] =
                        tlweSymEncrypt<typename P::targetP>(
                            domainkey[l * P::domainP::n + i] * (k + 1) *
                                (1ULL << (numeric_limits<
                                              typename P::targetP::T>::digits -
                                          (j + 1) * P::basebit)),
                            targetkey);
}

template <class P>
void ikskgen(KeySwitchingKey<P>& ksk, const SecretKey& sk)
{
    ikskgen<P>(ksk, sk.key.get<typename P::domainP>(),
               sk.key.get<typename P::targetP>());
}

template <class P>
void privkskgen(PrivateKeySwitchingKey<P>& privksk,
                const Polynomial<typename P::targetP>& func,
                const Key<typename P::domainP>& domainkey,
                const Key<typename P::targetP>& targetkey)
{
    std::array<typename P::domainP::T, P::domainP::k * P::domainP::n + 1> key;
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++)
        key[i] = domainkey[i];
    key[P::domainP::k * P::domainP::n] = -1;
#pragma omp parallel for collapse(3)
    for (int i = 0; i <= P::domainP::k * P::domainP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (typename P::targetP::T u = 0; u < (1 << (P::basebit - 1));
                 u++) {
                TRLWE<typename P::targetP> c =
                    trlweSymEncryptZero<typename P::targetP>(targetkey);
                for (int k = 0; k < P::targetP::n; k++)
                    c[P::targetP::k][k] +=
                        (u + 1) * func[k] * key[i]
                        << (numeric_limits<typename P::targetP::T>::digits -
                            (j + 1) * P::basebit);
                privksk[i][j][u] = c;
            }
}

template <class P>
void privkskgen(PrivateKeySwitchingKey<P>& privksk,
                const Polynomial<typename P::targetP>& func,
                const SecretKey& sk)
{
    privkskgen<P>(privksk, func, sk.key.get<typename P::domainP>(),
                  sk.key.get<typename P::targetP>());
}

template <class P>
void subikskgen(SubsetKeySwitchingKey<P>& ksk,
                const Key<typename P::domainP>& domainkey)
{
    Key<typename P::targetP> subkey;
    for (int i = 0; i < P::targetP::n; i++) subkey[i] = domainkey[i];
    for (int i = 0;
         i < P::domainP::k * P::domainP::n - P::targetP::k * P::targetP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (uint32_t k = 0; k < (1 << P::basebit) - 1; k++)
                ksk[i][j][k] = tlweSymEncrypt<typename P::targetP>(
                    domainkey[P::targetP::k * P::targetP::n + i] * (k + 1) *
                        (1ULL
                         << (numeric_limits<typename P::targetP::T>::digits -
                             (j + 1) * P::basebit)),
                    subkey);
}

template <class P>
void subikskgen(SubsetKeySwitchingKey<P>& ksk, const SecretKey& sk)
{
    subikskgen<P>(ksk, sk.key.get<typename P::domainP>());
}

template <class P>
relinKey<P> relinKeygen(const Key<P>& key)
{
    constexpr std::array<typename P::T, P::l> h = hgen<P>();

    Polynomial<P> keysquare;
    std::array<typename P::T, P::n> partkey;
    for (int i = 0; i < P::n; i++) partkey[i] = key[0 * P::n + i];
    PolyMulNaive<P>(keysquare, partkey, partkey);
    relinKey<P> relinkey;
    for (TRLWE<P>& ctxt : relinkey) ctxt = trlweSymEncryptZero<P>(key);
    for (int i = 0; i < P::l; i++)
        for (int j = 0; j < P::n; j++)
            relinkey[i][1][j] +=
                static_cast<typename P::T>(keysquare[j]) * h[i];
    return relinkey;
}

template <class P>
void subprivkskgen(SubsetPrivateKeySwitchingKey<P>& privksk,
                   const Polynomial<typename P::targetP>& func,
                   const Key<typename P::domainP>& domainkey,
                   const Key<typename P::targetP>& targetkey)
{
    std::array<typename P::targetP::T, P::targetP::k * P::targetP::n + 1> key;
    for (int i = 0; i < P::targetP::k * P::targetP::n; i++)
        key[i] = domainkey[i];
    key[P::targetP::k * P::targetP::n] = -1;
#pragma omp parallel for collapse(3)
    for (int i = 0; i <= P::targetP::k * P::targetP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (typename P::targetP::T u = 0; u < (1 << P::basebit) - 1; u++) {
                TRLWE<typename P::targetP> c =
                    trlweSymEncryptZero<typename P::targetP>(targetkey);
                for (int k = 0; k < P::targetP::n; k++)
                    c[P::targetP::k][k] +=
                        (u + 1) * func[k] * key[i]
                        << (numeric_limits<typename P::targetP::T>::digits -
                            (j + 1) * P::basebit);
                privksk[i][j][u] = c;
            }
}

template <class P>
void subprivkskgen(SubsetPrivateKeySwitchingKey<P>& privksk,
                   const Polynomial<typename P::targetP>& func,
                   const SecretKey& sk)
{
    subprivkskgen<P>(privksk, func, sk.key.get<typename P::domainP>(),
                     sk.key.get<typename P::targetP>());
}

template <class P>
relinKeyFFT<P> relinKeyFFTgen(const Key<P>& key)
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
    // BootstrapingKey
    std::shared_ptr<BootstrappingKey<lvl01param>> bklvl01;
    std::shared_ptr<BootstrappingKey<lvlh1param>> bklvlh1;
    std::shared_ptr<BootstrappingKey<lvl02param>> bklvl02;
    std::shared_ptr<BootstrappingKey<lvlh2param>> bklvlh2;
    // BoostrappingKeyFFT
    std::shared_ptr<BootstrappingKeyFFT<lvl01param>> bkfftlvl01;
    std::shared_ptr<BootstrappingKeyFFT<lvlh1param>> bkfftlvlh1;
    std::shared_ptr<BootstrappingKeyFFT<lvl02param>> bkfftlvl02;
    std::shared_ptr<BootstrappingKeyFFT<lvlh2param>> bkfftlvlh2;
    // BootstrappingKeyNTT
    std::shared_ptr<BootstrappingKeyNTT<lvl01param>> bknttlvl01;
    std::shared_ptr<BootstrappingKeyNTT<lvlh1param>> bknttlvlh1;
    std::shared_ptr<BootstrappingKeyNTT<lvl02param>> bknttlvl02;
    std::shared_ptr<BootstrappingKeyNTT<lvlh2param>> bknttlvlh2;
    // KeySwitchingKey
    std::shared_ptr<KeySwitchingKey<lvl10param>> iksklvl10;
    std::shared_ptr<KeySwitchingKey<lvl1hparam>> iksklvl1h;
    std::shared_ptr<KeySwitchingKey<lvl20param>> iksklvl20;
    std::shared_ptr<KeySwitchingKey<lvl2hparam>> iksklvl2h;
    std::shared_ptr<KeySwitchingKey<lvl21param>> iksklvl21;
    std::shared_ptr<KeySwitchingKey<lvl22param>> iksklvl22;
    std::shared_ptr<KeySwitchingKey<lvl31param>> iksklvl31;
    // SubsetKeySwitchingKey
    std::shared_ptr<SubsetKeySwitchingKey<lvl21param>> subiksklvl21;
    // PrivateKeySwitchingKey
    std::unordered_map<std::string,
                       std::shared_ptr<PrivateKeySwitchingKey<lvl11param>>>
        privksklvl11;
    std::unordered_map<std::string,
                       std::shared_ptr<PrivateKeySwitchingKey<lvl21param>>>
        privksklvl21;
    std::unordered_map<std::string,
                       std::shared_ptr<PrivateKeySwitchingKey<lvl22param>>>
        privksklvl22;
    // SubsetPrivateKeySwitchingKey
    std::unordered_map<
        std::string, std::shared_ptr<SubsetPrivateKeySwitchingKey<lvl21param>>>
        subprivksklvl21;

    EvalKey(SecretKey sk) { params = sk.params; }
    EvalKey() {}

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(params, bklvl01, bklvlh1, bklvl02, bklvlh2, bkfftlvl01,
                bkfftlvlh1, bkfftlvl02, bkfftlvlh2, bknttlvl01, bknttlvlh1,
                bknttlvl02, bknttlvlh2, iksklvl10, iksklvl1h, iksklvl20,
                iksklvl21, iksklvl22, iksklvl31, privksklvl11, privksklvl21,
                privksklvl22);
    }

    // emplace keys
    template <class P>
    void emplacebk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            bklvl01 =
                std::make_unique_for_overwrite<BootstrappingKey<lvl01param>>();
            bkgen<lvl01param>(*bklvl01, sk);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            bklvlh1 =
                std::make_unique_for_overwrite<BootstrappingKey<lvlh1param>>();
            bkgen<lvlh1param>(*bklvlh1, sk);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            bklvl02 =
                std::make_unique_for_overwrite<BootstrappingKey<lvl02param>>();
            bkgen<lvl02param>(*bklvl02, sk);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            bklvlh2 =
                std::make_unique_for_overwrite<BootstrappingKey<lvlh2param>>();
            bkgen<lvlh2param>(*bklvlh2, sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacebkfft(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            bkfftlvl01 = std::unique_ptr<BootstrappingKeyFFT<lvl01param>>(
                new (std::align_val_t(64)) BootstrappingKeyFFT<lvl01param>());
            bkfftgen<lvl01param>(*bkfftlvl01, sk);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            bkfftlvlh1 = std::unique_ptr<BootstrappingKeyFFT<lvlh1param>>(
                new (std::align_val_t(64)) BootstrappingKeyFFT<lvlh1param>());
            bkfftgen<lvlh1param>(*bkfftlvlh1, sk);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            bkfftlvl02 = std::unique_ptr<BootstrappingKeyFFT<lvl02param>>(
                new (std::align_val_t(64)) BootstrappingKeyFFT<lvl02param>());
            bkfftgen<lvl02param>(*bkfftlvl02, sk);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            bkfftlvlh2 = std::unique_ptr<BootstrappingKeyFFT<lvlh2param>>(
                new (std::align_val_t(64)) BootstrappingKeyFFT<lvlh2param>());
            bkfftgen<lvlh2param>(*bkfftlvlh2, sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacebkntt(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            bknttlvl01 = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvl01param>>();
            bknttgen<lvl01param>(*bknttlvl01, sk);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            bknttlvlh1 = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvlh1param>>();
            bknttgen<lvlh1param>(*bknttlvlh1, sk);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            bknttlvl02 = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvl02param>>();
            bknttgen<lvl02param>(*bknttlvl02, sk);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            bknttlvlh2 = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvlh2param>>();
            bknttgen<lvlh2param>(*bknttlvlh2, sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacebk2bkfft()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            bkfftlvl01 = std::make_unique_for_overwrite<
                BootstrappingKeyFFT<lvl01param>>();
            for (int i = 0; i < lvl01param::domainP::n; i++)
                (*bkfftlvl01)[i][0] =
                    ApplyFFT2trgsw<lvl1param>((*bklvl01)[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            bkfftlvlh1 = std::make_unique_for_overwrite<
                BootstrappingKeyFFT<lvlh1param>>();
            for (int i = 0; i < lvlh1param::domainP::n; i++)
                (*bkfftlvlh1)[i][0] =
                    ApplyFFT2trgsw<lvl1param>((*bklvlh1)[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            bkfftlvl02 = std::make_unique_for_overwrite<
                BootstrappingKeyFFT<lvl02param>>();
            for (int i = 0; i < lvl02param::domainP::n; i++)
                (*bkfftlvl02)[i][0] =
                    ApplyFFT2trgsw<lvl2param>((*bklvl02)[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            bkfftlvlh2 = std::make_unique_for_overwrite<
                BootstrappingKeyFFT<lvlh2param>>();
            for (int i = 0; i < lvlh2param::domainP::n; i++)
                (*bkfftlvlh2)[i][0] =
                    ApplyFFT2trgsw<lvl2param>((*bklvlh2)[i][0]);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacebk2bkntt()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            bknttlvl01 = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvl01param>>();
            for (int i = 0; i < lvl01param::domainP::n; i++)
                (*bknttlvl01)[i] = ApplyNTT2trgsw<lvl1param>((*bklvl01)[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            bknttlvlh1 = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvlh1param>>();
            for (int i = 0; i < lvlh1param::domainP::n; i++)
                (*bknttlvlh1)[i] = ApplyNTT2trgsw<lvl1param>((*bklvlh1)[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            bknttlvl02 = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvl02param>>();
            for (int i = 0; i < lvl02param::domainP::n; i++)
                (*bknttlvl02)[i] = ApplyNTT2trgsw<lvl2param>((*bklvl02)[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            bknttlvlh2 = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvlh2param>>();
            for (int i = 0; i < lvlh2param::domainP::n; i++)
                (*bknttlvlh2)[i] = ApplyNTT2trgsw<lvl2param>((*bklvlh2)[i][0]);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplaceiksk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl10param>) {
            iksklvl10 = std::unique_ptr<KeySwitchingKey<lvl10param>>(
                new (std::align_val_t(64)) KeySwitchingKey<lvl10param>());
            ikskgen<lvl10param>(*iksklvl10, sk);
        }
        else if constexpr (std::is_same_v<P, lvl1hparam>) {
            iksklvl1h = std::unique_ptr<KeySwitchingKey<lvl1hparam>>(
                new (std::align_val_t(64)) KeySwitchingKey<lvl1hparam>());
            ikskgen<lvl1hparam>(*iksklvl1h, sk);
        }
        else if constexpr (std::is_same_v<P, lvl20param>) {
            iksklvl20 = std::unique_ptr<KeySwitchingKey<lvl20param>>(
                new (std::align_val_t(64)) KeySwitchingKey<lvl20param>());
            ikskgen<lvl20param>(*iksklvl20, sk);
        }
        else if constexpr (std::is_same_v<P, lvl2hparam>) {
            iksklvl2h =
                std::make_unique_for_overwrite<KeySwitchingKey<lvl2hparam>>();
            ikskgen<lvl2hparam>(*iksklvl2h, sk);
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            iksklvl21 = std::unique_ptr<KeySwitchingKey<lvl21param>>(
                new (std::align_val_t(64)) KeySwitchingKey<lvl21param>());
            ikskgen<lvl21param>(*iksklvl21, sk);
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            iksklvl22 = std::unique_ptr<KeySwitchingKey<lvl22param>>(
                new (std::align_val_t(64)) KeySwitchingKey<lvl22param>());
            ikskgen<lvl22param>(*iksklvl22, sk);
        }
        else if constexpr (std::is_same_v<P, lvl31param>) {
            iksklvl31 = std::unique_ptr<KeySwitchingKey<lvl31param>>(
                new (std::align_val_t(64)) KeySwitchingKey<lvl31param>());
            ikskgen<lvl31param>(*iksklvl31, sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacesubiksk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            subiksklvl21 = std::make_unique_for_overwrite<
                SubsetKeySwitchingKey<lvl21param>>();
            subikskgen<lvl21param>(*subiksklvl21, sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplaceprivksk(const std::string& key,
                        const Polynomial<typename P::targetP>& func,
                        const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl11param>) {
            privksklvl11[key] =
                std::unique_ptr<PrivateKeySwitchingKey<lvl11param>>(new (
                    std::align_val_t(64)) PrivateKeySwitchingKey<lvl11param>());
            privkskgen<lvl11param>(*privksklvl11[key], func, sk);
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            privksklvl21[key] =
                std::unique_ptr<PrivateKeySwitchingKey<lvl21param>>(new (
                    std::align_val_t(64)) PrivateKeySwitchingKey<lvl21param>());
            privkskgen<lvl21param>(*privksklvl21[key], func, sk);
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            privksklvl22[key] =
                std::unique_ptr<PrivateKeySwitchingKey<lvl22param>>(new (
                    std::align_val_t(64)) PrivateKeySwitchingKey<lvl22param>());
            privkskgen<lvl22param>(*privksklvl22[key], func, sk);
        }
        else
            static_assert(false_v<typename P::targetP::T>,
                          "Not predefined parameter!");
    }
    template <class P>
    void emplacesubprivksk(const std::string& key,
                           const Polynomial<typename P::targetP>& func,
                           const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            subprivksklvl21[key] = std::make_unique_for_overwrite<
                SubsetPrivateKeySwitchingKey<lvl21param>>();
            subprivkskgen<lvl21param>(*subprivksklvl21[key], func, sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplaceprivksk4cb(const SecretKey& sk)
    {
        for (int k = 0; k < P::targetP::k; k++) {
            Polynomial<typename P::targetP> partkey;
            for (int i = 0; i < P::targetP::n; i++)
                partkey[i] =
                    -sk.key.get<typename P::targetP>()[k * P::targetP::n + i];
            emplaceprivksk<P>("privksk4cb_" + std::to_string(k), partkey, sk);
        }
        emplaceprivksk<P>("privksk4cb_" + std::to_string(P::targetP::k), {1},
                          sk);
    }
    template <class P>
    void emplacesubprivksk4cb(const SecretKey& sk)
    {
        for (int k = 0; k < P::targetP::k; k++) {
            Polynomial<typename P::targetP> partkey;
            for (int i = 0; i < P::targetP::n; i++)
                partkey[i] =
                    -sk.key.get<typename P::targetP>()[k * P::targetP::n + i];
            emplacesubprivksk<P>("subprivksk4cb_" + std::to_string(k), partkey,
                                 sk);
        }
        emplacesubprivksk<P>("subprivksk4cb_" + std::to_string(P::targetP::k),
                             {1}, sk);
    }

    // get keys
    template <class P>
    BootstrappingKey<P>& getbk() const
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return *bklvl01;
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return *bklvlh1;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *bklvl02;
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return *bklvlh2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    BootstrappingKeyFFT<P>& getbkfft() const
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return *bkfftlvl01;
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return *bkfftlvlh1;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *bkfftlvl02;
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return *bkfftlvlh2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    BootstrappingKeyNTT<P>& getbkntt() const
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return *bknttlvl01;
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return *bknttlvlh1;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *bknttlvl02;
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return *bknttlvlh2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    KeySwitchingKey<P>& getiksk() const
    {
        if constexpr (std::is_same_v<P, lvl10param>) {
            return *iksklvl10;
        }
        else if constexpr (std::is_same_v<P, lvl1hparam>) {
            return *iksklvl1h;
        }
        else if constexpr (std::is_same_v<P, lvl20param>) {
            return *iksklvl20;
        }
        else if constexpr (std::is_same_v<P, lvl2hparam>) {
            return *iksklvl2h;
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            return *iksklvl21;
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            return *iksklvl22;
        }
        else if constexpr (std::is_same_v<P, lvl31param>) {
            return *iksklvl31;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    SubsetKeySwitchingKey<P>& getsubiksk() const
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return *subiksklvl21;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    PrivateKeySwitchingKey<P>& getprivksk(const std::string& key) const
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
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    SubsetPrivateKeySwitchingKey<P>& getsubprivksk(const std::string& key) const
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return *(subprivksklvl21.at(key));
        }
        else
            static_assert(false_v<typename P::targetP::T>,
                          "Not predefined parameter!");
    }
};

}  // namespace TFHEpp