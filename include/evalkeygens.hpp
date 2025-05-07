#pragma once

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
    for (int j = 0; j < P::k; j++) {
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
        evalautokeygen<P>(ahk[i], (1 << (i + 1)) + 1, key);
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
    for (int i = 0; i < P::targetP::k * P::targetP::n; i++)
        subkey[i] = domainkey[i];
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
    Polynomial<P> keysquare;
    std::array<typename P::T, P::n> partkey;
    for (int i = 0; i < P::n; i++) partkey[i] = key[0 * P::n + i];
    PolyMulNaive<P>(keysquare, partkey, partkey);
    relinKey<P> relinkey;
    for (TRLWE<P>& ctxt : relinkey) ctxt = trlweSymEncryptZero<P>(key);
    halftrgswhadd<P>(relinkey, keysquare);
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
}  // namespace TFHEpp