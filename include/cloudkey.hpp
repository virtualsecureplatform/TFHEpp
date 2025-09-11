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
    std::shared_ptr<AnnihilateKey<AHlvl1param>> ahklvl1;
    std::shared_ptr<AnnihilateKey<AHlvl2param>> ahklvl2;
    // CBswitchKey
    std::shared_ptr<CBswitchingKey<AHlvl1param>> cbsklvl1;
    std::shared_ptr<CBswitchingKey<AHlvl2param>> cbsklvl2;

    EvalKey(SecretKey sk) { params = sk.params; }
    EvalKey() {}

    // ptr* functions for accessing shared_ptr members
    template <class P>
    std::shared_ptr<BootstrappingKey<P>>& ptrbk()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return bklvl01;
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return bklvlh1;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return bklvl02;
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return bklvlh2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    std::shared_ptr<BootstrappingKeyFFT<P>>& ptrbkfft()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return bkfftlvl01;
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return bkfftlvlh1;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return bkfftlvl02;
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return bkfftlvlh2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    std::shared_ptr<BootstrappingKeyNTT<P>>& ptrbkntt()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return bknttlvl01;
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return bknttlvlh1;
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return bknttlvl02;
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return bknttlvlh2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    std::shared_ptr<KeySwitchingKey<P>>& ptriksk()
    {
        if constexpr (std::is_same_v<P, lvl10param>) {
            return iksklvl10;
        }
        else if constexpr (std::is_same_v<P, lvl1hparam>) {
            return iksklvl1h;
        }
        else if constexpr (std::is_same_v<P, lvl20param>) {
            return iksklvl20;
        }
        else if constexpr (std::is_same_v<P, lvl2hparam>) {
            return iksklvl2h;
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            return iksklvl21;
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            return iksklvl22;
        }
        else if constexpr (std::is_same_v<P, lvl31param>) {
            return iksklvl31;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    std::shared_ptr<SubsetKeySwitchingKey<P>>& ptrsubiksk()
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return subiksklvl21;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    std::unordered_map<std::string, std::shared_ptr<PrivateKeySwitchingKey<P>>>& ptrprivksk()
    {
        if constexpr (std::is_same_v<P, lvl11param>) {
            return privksklvl11;
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            return privksklvl21;
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            return privksklvl22;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    std::unordered_map<std::string, std::shared_ptr<SubsetPrivateKeySwitchingKey<P>>>& ptrsubprivksk()
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return subprivksklvl21;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    std::shared_ptr<AnnihilateKey<P>>& ptrahk()
    {
        if constexpr (std::is_same_v<P, AHlvl1param>) {
            return ahklvl1;
        }
        else if constexpr (std::is_same_v<P, AHlvl2param>) {
            return ahklvl2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    std::shared_ptr<CBswitchingKey<P>>& ptrcbsk()
    {
        if constexpr (std::is_same_v<P, AHlvl1param>) {
            return cbsklvl1;
        }
        else if constexpr (std::is_same_v<P, AHlvl2param>) {
            return cbsklvl2;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

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
        ptrbk<P>() = std::make_unique_for_overwrite<BootstrappingKey<P>>();
        bkgen<P>(*ptrbk<P>(), sk);
    }
    template <class P>
    void emplacebkfft(const SecretKey& sk)
    {
        ptrbkfft<P>() = std::make_unique_for_overwrite<BootstrappingKeyFFT<P>>();
        bkfftgen<P>(*ptrbkfft<P>(), sk);
    }
    template <class P>
    void emplacebkntt(const SecretKey& sk)
    {
        ptrbkntt<P>() = std::make_unique_for_overwrite<BootstrappingKeyNTT<P>>();
        bknttgen<P>(*ptrbkntt<P>(), sk);
    }
    template <class P>
    void emplacebk2bkfft()
    {
        ptrbkfft<P>() = std::make_unique_for_overwrite<BootstrappingKeyFFT<P>>();
        for (int i = 0; i < P::domainP::n; i++) {
            if constexpr (std::is_same_v<P, lvl01param> || std::is_same_v<P, lvlh1param>) {
                (*ptrbkfft<P>())[i][0] = ApplyFFT2trgsw<lvl1param>((*ptrbk<P>())[i][0]);
            }
            else if constexpr (std::is_same_v<P, lvl02param> || std::is_same_v<P, lvlh2param>) {
                (*ptrbkfft<P>())[i][0] = ApplyFFT2trgsw<lvl2param>((*ptrbk<P>())[i][0]);
            }
            else {
                static_assert(false_v<typename P::T>, "Not predefined parameter!");
            }
        }
    }
    template <class P>
    void emplacebk2bkntt()
    {
        ptrbkntt<P>() = std::make_unique_for_overwrite<BootstrappingKeyNTT<P>>();
        for (int i = 0; i < P::domainP::n; i++) {
            if constexpr (std::is_same_v<P, lvl01param> || std::is_same_v<P, lvlh1param>) {
                (*ptrbkntt<P>())[i] = ApplyNTT2trgsw<lvl1param>((*ptrbk<P>())[i][0]);
            }
            else if constexpr (std::is_same_v<P, lvl02param> || std::is_same_v<P, lvlh2param>) {
                (*ptrbkntt<P>())[i] = ApplyNTT2trgsw<lvl2param>((*ptrbk<P>())[i][0]);
            }
            else {
                static_assert(false_v<typename P::T>, "Not predefined parameter!");
            }
        }
    }
    template <class P>
    void emplaceiksk(const SecretKey& sk)
    {
        ptriksk<P>() = std::make_unique_for_overwrite<KeySwitchingKey<P>>();
        ikskgen<P>(*ptriksk<P>(), sk);
    }
    template <class P>
    void emplacesubiksk(const SecretKey& sk)
    {
        ptrsubiksk<P>() = std::make_unique_for_overwrite<SubsetKeySwitchingKey<P>>();
        subikskgen<P>(*ptrsubiksk<P>(), sk);
    }
    template <class P>
    void emplaceprivksk(const std::string& key,
                        const Polynomial<typename P::targetP>& func,
                        const SecretKey& sk)
    {
        ptrprivksk<P>()[key] =
            std::unique_ptr<PrivateKeySwitchingKey<P>>(new (
                std::align_val_t(64)) PrivateKeySwitchingKey<P>());
        privkskgen<P>(*ptrprivksk<P>()[key], func, sk);
    }
    template <class P>
    void emplacesubprivksk(const std::string& key,
                           const Polynomial<typename P::targetP>& func,
                           const SecretKey& sk)
    {
        ptrsubprivksk<P>()[key] =
            std::make_unique_for_overwrite<SubsetPrivateKeySwitchingKey<P>>();
        subprivkskgen<P>(*ptrsubprivksk<P>()[key], func, sk);
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
        ptrahk<P>() = std::make_unique_for_overwrite<AnnihilateKey<P>>();
        annihilatekeygen<P>(*ptrahk<P>(), sk);
    }

    template <class P>
    void emplacecbsk(const SecretKey& sk)
    {
        ptrcbsk<P>() = std::make_unique_for_overwrite<CBswitchingKey<P>>();
        for (int i = 0; i < P::k; i++) {
            Polynomial<P> partkey;
            for (int j = 0; j < P::n; j++)
                partkey[j] = -sk.key.get<P>()[i * P::n + j];
            (*ptrcbsk<P>())[i] = trgswfftSymEncrypt<P>(partkey, sk.key.get<P>());
        }
    }

    // get keys
    template <class P>
    BootstrappingKey<P>& getbk() const
    {
        return *const_cast<EvalKey*>(this)->ptrbk<P>();
    }
    template <class P>
    BootstrappingKeyFFT<P>& getbkfft() const
    {
        return *const_cast<EvalKey*>(this)->ptrbkfft<P>();
    }
    template <class P>
    BootstrappingKeyNTT<P>& getbkntt() const
    {
        return *const_cast<EvalKey*>(this)->ptrbkntt<P>();
    }
    template <class P>
    KeySwitchingKey<P>& getiksk() const
    {
        return *const_cast<EvalKey*>(this)->ptriksk<P>();
    }
    template <class P>
    SubsetKeySwitchingKey<P>& getsubiksk() const
    {
        return *const_cast<EvalKey*>(this)->ptrsubiksk<P>();
    }
    template <class P>
    PrivateKeySwitchingKey<P>& getprivksk(const std::string& key) const
    {
        return *(const_cast<EvalKey*>(this)->ptrprivksk<P>().at(key));
    }
    template <class P>
    SubsetPrivateKeySwitchingKey<P>& getsubprivksk(const std::string& key) const
    {
        return *(const_cast<EvalKey*>(this)->ptrsubprivksk<P>().at(key));
    }
    template <class P>
    AnnihilateKey<P>& getahk() const
    {
        return *const_cast<EvalKey*>(this)->ptrahk<P>();
    }
    template <class P>
    CBswitchingKey<P>& getcbsk() const
    {
        return *const_cast<EvalKey*>(this)->ptrcbsk<P>();
    }
};

}  // namespace TFHEpp