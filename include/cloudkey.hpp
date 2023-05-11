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
void bkgen(BootstrappingKey<P>& bk, const SecretKey& sk)
{
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        int count = 0;
        for (int j = P::domainP::key_value_min; j <= P::domainP::key_value_max;
             j++) {
            plainpoly[0] = sk.key.get<typename P::domainP>()[i] == j;
            if (j != 0) {
                plainpoly[0] = sk.key.get<typename P::domainP>()[i];
                bk[i][count] = trgswSymEncrypt<typename P::targetP>(
                    plainpoly, P::targetP::α,
                    sk.key.get<typename P::targetP>());
                count++;
            }
        }
    }
}

template <class P>
void bkfftgen(BootstrappingKeyFFT<P>& bkfft, const SecretKey& sk)
{
    Polynomial<typename P::targetP> plainpoly = {};
#ifdef USE_KEY_BUNDLE
    for (int i = 0; i < P::domainP::k * P::domainP::n / P::Addends; i++) {
        plainpoly[0] =
            static_cast<int32_t>(sk.key.get<typename P::domainP>()[2 * i] *
                                 sk.key.get<typename P::domainP>()[2 * i + 1]);
        bkfft[i][0] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
        plainpoly[0] = static_cast<int32_t>(
            sk.key.get<typename P::domainP>()[2 * i] *
            (1 - sk.key.get<typename P::domainP>()[2 * i + 1]));
        bkfft[i][1] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
        plainpoly[0] = static_cast<int32_t>(
            (1 - sk.key.get<typename P::domainP>()[2 * i]) *
            sk.key.get<typename P::domainP>()[2 * i + 1]);
        bkfft[i][2] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
#else
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        int count = 0;
        for (int j = P::domainP::key_value_min; j <= P::domainP::key_value_max;
             j++) {
            if (j != 0) {
                plainpoly[0] = sk.key.get<typename P::domainP>()[i] == j;
                bkfft[i][count] = trgswfftSymEncrypt<typename P::targetP>(
                    plainpoly, P::targetP::α,
                    sk.key.get<typename P::targetP>());
                count++;
            }
        }
    }
#endif
}

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

template <class P>
void privkskgen(PrivateKeySwitchingKey<P>& privksk,
                const Polynomial<typename P::targetP>& func,
                const SecretKey& sk)
{
    std::array<typename P::domainP::T, P::domainP::k * P::domainP::n + 1> key;
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

template <class P>
void subikskgen(SubsetKeySwitchingKey<P>& ksk, const SecretKey& sk)
{
    Key<typename P::targetP> subkey;
    for(int i = 0; i < P::targetP::n; i++) subkey[i] = sk.key.get<typename P::domainP>()[i];
    for (int i = 0; i < P::domainP::k * P::domainP::n - P::targetP::k * P::targetP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (uint32_t k = 0; k < (1 << P::basebit) - 1; k++)
                ksk[i][j][k] =
                    tlweSymEncrypt<typename P::targetP>(
                        sk.key.get<
                            typename P::domainP>()[P::targetP::k * P::targetP::n + i] *
                            (k + 1) *
                            (1ULL << (numeric_limits<
                                            typename P::targetP::T>::digits -
                                        (j + 1) * P::basebit)),
                        P::α, subkey);
}

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
void subprivkskgen(SubsetPrivateKeySwitchingKey<P>& privksk,
                const Polynomial<typename P::targetP>& func,
                const SecretKey& sk)
{
    std::array<typename P::targetP::T, P::targetP::k * P::targetP::n + 1> key;
    for (int i = 0; i < P::targetP::k * P::targetP::n; i++)
        key[i] = sk.key.get<typename P::domainP>()[i];
    key[P::targetP::k * P::targetP::n] = -1;
#pragma omp parallel for collapse(3)
    for (int i = 0; i <= P::targetP::k * P::targetP::n; i++)
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

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(params, bklvl01, bklvl02, bkfftlvl01, bkfftlvl02, bknttlvl01,
                bknttlvl02, iksklvl10, iksklvl11, iksklvl20, iksklvl21,
                iksklvl22, privksklvl11, privksklvl21, privksklvl22);
    }

    //emplace keys
    template <class P>
    void emplacebk(const SecretKey& sk)
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
    template <class P>
    void emplacebkfft(const SecretKey& sk)
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
    template <class P>
    void emplacebkntt(const SecretKey& sk)
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
    template <class P>
    void emplacebk2bkfft()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            bkfftlvl01 = std::make_unique<BootstrappingKeyFFT<lvl01param>>();
            for (int i = 0; i < lvl01param::domainP::n; i++)
                (*bkfftlvl01)[i][0] = ApplyFFT2trgsw<lvl1param>((*bklvl01)[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            bkfftlvl02 = std::make_unique<BootstrappingKeyFFT<lvl02param>>();
            for (int i = 0; i < lvl02param::domainP::n; i++)
                (*bkfftlvl02)[i][0] = ApplyFFT2trgsw<lvl2param>((*bklvl02)[i][0]);
        }
    }
    template <class P>
    void emplacebk2bkntt()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            bknttlvl01 = std::make_unique<BootstrappingKeyNTT<lvl01param>>();
            for (int i = 0; i < lvl01param::domainP::n; i++)
                (*bknttlvl01)[i] = ApplyNTT2trgsw<lvl1param>((*bklvl01)[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            bknttlvl02 = std::make_unique<BootstrappingKeyNTT<lvl02param>>();
            for (int i = 0; i < lvl02param::domainP::n; i++)
                (*bknttlvl02)[i] = ApplyNTT2trgsw<lvl2param>((*bklvl02)[i][0]);
        }
    }
    template <class P>
    void emplaceiksk(const SecretKey& sk)
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
    template <class P>
    void emplacesubiksk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            subiksklvl21 = std::make_unique<SubsetKeySwitchingKey<lvl21param>>();
            subikskgen<lvl21param>(*subiksklvl21, sk);
        }
    }
    template <class P>
    void emplaceprivksk(const std::string& key,
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
    template <class P>
    void emplacesubprivksk(const std::string& key,
                                const Polynomial<typename P::targetP>& func,
                                const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            subprivksklvl21[key] =
                std::make_unique<SubsetPrivateKeySwitchingKey<lvl21param>>();
            subprivkskgen<lvl21param>(*subprivksklvl21[key], func, sk);
        }
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
        emplaceprivksk<P>("privksk4cb_" + std::to_string(P::targetP::k), {1}, sk);
    }
    template <class P>
    void emplacesubprivksk4cb(const SecretKey& sk)
    {
        for (int k = 0; k < P::targetP::k; k++) {
            Polynomial<typename P::targetP> partkey;
            for (int i = 0; i < P::targetP::n; i++)
                partkey[i] =
                    -sk.key.get<typename P::targetP>()[k * P::targetP::n + i];
            emplacesubprivksk<P>("subprivksk4cb_" + std::to_string(k), partkey, sk);
        }
        emplacesubprivksk<P>("subprivksk4cb_" + std::to_string(P::targetP::k), {1}, sk);
    }

    //get keys
    template <class P>
    BootstrappingKey<P>& getbk() const
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return *bklvl01;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *bklvl02;
        }
    }
    template <class P>
    BootstrappingKeyFFT<P>& getbkfft() const
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return *bkfftlvl01;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *bkfftlvl02;
        }
    }
    template <class P>
    BootstrappingKeyNTT<P>& getbkntt() const
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return *bknttlvl01;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *bknttlvl02;
        }
    }
    template <class P>
    KeySwitchingKey<P>& getiksk() const
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
    template <class P>
    SubsetKeySwitchingKey<P>& getsubiksk() const
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return *subiksklvl21;
        }
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
    }
    template <class P>
    SubsetPrivateKeySwitchingKey<P>& getsubprivksk(const std::string& key) const
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return *(subprivksklvl21.at(key));
        }
    }
};

#include "externs/cloudkey.hpp"
}  // namespace TFHEpp