#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>

#include "evalkeygens.hpp"

namespace TFHEpp {

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
    // AnnihilateKey
    std::shared_ptr<AnnihilateKey<lvl1param>> ahklvl1;
    std::shared_ptr<AnnihilateKey<lvl2param>> ahklvl2;
    // CBswitchKey
    std::shared_ptr<CBswitchingKey<lvl1param>> cbsklvl1;
    std::shared_ptr<CBswitchingKey<lvl2param>> cbsklvl2;

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
            subprivksklvl21[key] =
                std::unique_ptr<SubsetPrivateKeySwitchingKey<lvl21param>>(new (
                    std::align_val_t(64)) SubsetPrivateKeySwitchingKey<lvl21param>());
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

    template <class P>
    void emplaceahk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl1param>) {
            ahklvl1 = std::make_unique_for_overwrite<AnnihilateKey<lvl1param>>();
            annihilatekeygen<lvl1param>(*ahklvl1, sk);
        }
        else if constexpr (std::is_same_v<P, lvl2param>) {
            ahklvl2 = std::make_unique_for_overwrite<AnnihilateKey<lvl2param>>();
            annihilatekeygen<lvl2param>(*ahklvl2, sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    void emplacecbsk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl1param>) {
            cbsklvl1 = std::make_unique_for_overwrite<CBswitchingKey<lvl1param>>();
            for (int i = 0; i < lvl1param::k; i++){
                Polynomial<P> partkey;
                for (int j = 0; j < P::n; j++)
                    partkey[j] =
                        -sk.key.get<P>()[i * P::n + j];
                (*cbsklvl1)[i] = trgswfftSymEncrypt<P>(partkey,sk.key.get<P>());
            }
        }
        else if constexpr (std::is_same_v<P, lvl2param>) {
            cbsklvl2 = std::make_unique_for_overwrite<CBswitchingKey<lvl2param>>();
            for (int i = 0; i < lvl2param::k; i++){
                Polynomial<P> partkey;
                for (int j = 0; j < P::n; j++)
                    partkey[j] =
                        -sk.key.get<P>()[i * P::n + j];
                (*cbsklvl2)[i] = trgswfftSymEncrypt<P>(partkey,sk.key.get<P>());
            }
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
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
    template <class P>
    AnnihilateKey<P>& getahk() const
    {
        if constexpr (std::is_same_v<P, lvl1param>) {
            return *ahklvl1;
        }
        else if constexpr (std::is_same_v<P, lvl2param>) {
            return *ahklvl2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    CBswitchingKey<P>& getcbsk() const
    {
        if constexpr (std::is_same_v<P, lvl1param>) {
            return *cbsklvl1;
        }
        else if constexpr (std::is_same_v<P, lvl2param>) {
            return *cbsklvl2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
};

}  // namespace TFHEpp