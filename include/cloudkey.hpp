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

    // const versions of ptr* functions
    template <class P>
    const std::shared_ptr<BootstrappingKey<P>>& ptrbk() const
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
    const std::shared_ptr<BootstrappingKeyFFT<P>>& ptrbkfft() const
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
    const std::shared_ptr<BootstrappingKeyNTT<P>>& ptrbkntt() const
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
    const std::shared_ptr<KeySwitchingKey<P>>& ptriksk() const
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
    const std::shared_ptr<SubsetKeySwitchingKey<P>>& ptrsubiksk() const
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return subiksklvl21;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    const std::unordered_map<std::string, std::shared_ptr<PrivateKeySwitchingKey<P>>>& ptrprivksk() const
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
    const std::unordered_map<std::string, std::shared_ptr<SubsetPrivateKeySwitchingKey<P>>>& ptrsubprivksk() const
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return subprivksklvl21;
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    const std::shared_ptr<AnnihilateKey<P>>& ptrahk() const
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
    const std::shared_ptr<CBswitchingKey<P>>& ptrcbsk() const
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
        if constexpr (std::is_same_v<P, lvl01param>) {
            ptrbk<lvl01param>() =
                std::make_unique_for_overwrite<BootstrappingKey<lvl01param>>();
            bkgen<lvl01param>(*ptrbk<lvl01param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            ptrbk<lvlh1param>() =
                std::make_unique_for_overwrite<BootstrappingKey<lvlh1param>>();
            bkgen<lvlh1param>(*ptrbk<lvlh1param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            ptrbk<lvl02param>() =
                std::make_unique_for_overwrite<BootstrappingKey<lvl02param>>();
            bkgen<lvl02param>(*ptrbk<lvl02param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            ptrbk<lvlh2param>() =
                std::make_unique_for_overwrite<BootstrappingKey<lvlh2param>>();
            bkgen<lvlh2param>(*ptrbk<lvlh2param>(), sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacebkfft(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            ptrbkfft<lvl01param>() = std::make_unique_for_overwrite<BootstrappingKeyFFT<lvl01param>>();
            bkfftgen<lvl01param>(*ptrbkfft<lvl01param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            ptrbkfft<lvlh1param>() = std::make_unique_for_overwrite<BootstrappingKeyFFT<lvlh1param>>();
            bkfftgen<lvlh1param>(*ptrbkfft<lvlh1param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            ptrbkfft<lvl02param>() = std::make_unique_for_overwrite<BootstrappingKeyFFT<lvl02param>>();
            bkfftgen<lvl02param>(*ptrbkfft<lvl02param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            ptrbkfft<lvlh2param>() = std::make_unique_for_overwrite<BootstrappingKeyFFT<lvlh2param>>();
            bkfftgen<lvlh2param>(*ptrbkfft<lvlh2param>(), sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacebkntt(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            ptrbkntt<lvl01param>() = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvl01param>>();
            bknttgen<lvl01param>(*ptrbkntt<lvl01param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            ptrbkntt<lvlh1param>() = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvlh1param>>();
            bknttgen<lvlh1param>(*ptrbkntt<lvlh1param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            ptrbkntt<lvl02param>() = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvl02param>>();
            bknttgen<lvl02param>(*ptrbkntt<lvl02param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            ptrbkntt<lvlh2param>() = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvlh2param>>();
            bknttgen<lvlh2param>(*ptrbkntt<lvlh2param>(), sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacebk2bkfft()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            ptrbkfft<lvl01param>() = std::make_unique_for_overwrite<
                BootstrappingKeyFFT<lvl01param>>();
            for (int i = 0; i < lvl01param::domainP::n; i++)
                (*ptrbkfft<lvl01param>())[i][0] =
                    ApplyFFT2trgsw<lvl1param>((*ptrbk<lvl01param>())[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            ptrbkfft<lvlh1param>() = std::make_unique_for_overwrite<
                BootstrappingKeyFFT<lvlh1param>>();
            for (int i = 0; i < lvlh1param::domainP::n; i++)
                (*ptrbkfft<lvlh1param>())[i][0] =
                    ApplyFFT2trgsw<lvl1param>((*ptrbk<lvlh1param>())[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            ptrbkfft<lvl02param>() = std::make_unique_for_overwrite<
                BootstrappingKeyFFT<lvl02param>>();
            for (int i = 0; i < lvl02param::domainP::n; i++)
                (*ptrbkfft<lvl02param>())[i][0] =
                    ApplyFFT2trgsw<lvl2param>((*ptrbk<lvl02param>())[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            ptrbkfft<lvlh2param>() = std::make_unique_for_overwrite<
                BootstrappingKeyFFT<lvlh2param>>();
            for (int i = 0; i < lvlh2param::domainP::n; i++)
                (*ptrbkfft<lvlh2param>())[i][0] =
                    ApplyFFT2trgsw<lvl2param>((*ptrbk<lvlh2param>())[i][0]);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacebk2bkntt()
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            ptrbkntt<lvl01param>() = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvl01param>>();
            for (int i = 0; i < lvl01param::domainP::n; i++)
                (*ptrbkntt<lvl01param>())[i] = ApplyNTT2trgsw<lvl1param>((*ptrbk<lvl01param>())[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            ptrbkntt<lvlh1param>() = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvlh1param>>();
            for (int i = 0; i < lvlh1param::domainP::n; i++)
                (*ptrbkntt<lvlh1param>())[i] = ApplyNTT2trgsw<lvl1param>((*ptrbk<lvlh1param>())[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            ptrbkntt<lvl02param>() = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvl02param>>();
            for (int i = 0; i < lvl02param::domainP::n; i++)
                (*ptrbkntt<lvl02param>())[i] = ApplyNTT2trgsw<lvl2param>((*ptrbk<lvl02param>())[i][0]);
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            ptrbkntt<lvlh2param>() = std::make_unique_for_overwrite<
                BootstrappingKeyNTT<lvlh2param>>();
            for (int i = 0; i < lvlh2param::domainP::n; i++)
                (*ptrbkntt<lvlh2param>())[i] = ApplyNTT2trgsw<lvl2param>((*ptrbk<lvlh2param>())[i][0]);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplaceiksk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl10param>) {
            ptriksk<lvl10param>() = std::make_unique_for_overwrite<KeySwitchingKey<lvl10param>>();
            ikskgen<lvl10param>(*ptriksk<lvl10param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl1hparam>) {
            ptriksk<lvl1hparam>() = std::make_unique_for_overwrite<KeySwitchingKey<lvl1hparam>>();
            ikskgen<lvl1hparam>(*ptriksk<lvl1hparam>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl20param>) {
            ptriksk<lvl20param>() = std::make_unique_for_overwrite<KeySwitchingKey<lvl20param>>();
            ikskgen<lvl20param>(*ptriksk<lvl20param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl2hparam>) {
            ptriksk<lvl2hparam>() =
                std::make_unique_for_overwrite<KeySwitchingKey<lvl2hparam>>();
            ikskgen<lvl2hparam>(*ptriksk<lvl2hparam>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            ptriksk<lvl21param>() = std::make_unique_for_overwrite<KeySwitchingKey<lvl21param>>();
            ikskgen<lvl21param>(*ptriksk<lvl21param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            ptriksk<lvl22param>() = std::make_unique_for_overwrite<KeySwitchingKey<lvl22param>>();
            ikskgen<lvl22param>(*ptriksk<lvl22param>(), sk);
        }
        else if constexpr (std::is_same_v<P, lvl31param>) {
            ptriksk<lvl31param>() = std::make_unique_for_overwrite<KeySwitchingKey<lvl31param>>();
            ikskgen<lvl31param>(*ptriksk<lvl31param>(), sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    void emplacesubiksk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            ptrsubiksk<lvl21param>() = std::make_unique_for_overwrite<
                SubsetKeySwitchingKey<lvl21param>>();
            subikskgen<lvl21param>(*ptrsubiksk<lvl21param>(), sk);
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
            ptrprivksk<lvl11param>()[key] =
                std::unique_ptr<PrivateKeySwitchingKey<lvl11param>>(new (
                    std::align_val_t(64)) PrivateKeySwitchingKey<lvl11param>());
            privkskgen<lvl11param>(*ptrprivksk<lvl11param>()[key], func, sk);
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            ptrprivksk<lvl21param>()[key] =
                std::unique_ptr<PrivateKeySwitchingKey<lvl21param>>(new (
                    std::align_val_t(64)) PrivateKeySwitchingKey<lvl21param>());
            privkskgen<lvl21param>(*ptrprivksk<lvl21param>()[key], func, sk);
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            ptrprivksk<lvl22param>()[key] =
                std::unique_ptr<PrivateKeySwitchingKey<lvl22param>>(new (
                    std::align_val_t(64)) PrivateKeySwitchingKey<lvl22param>());
            privkskgen<lvl22param>(*ptrprivksk<lvl22param>()[key], func, sk);
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
            ptrsubprivksk<lvl21param>()[key] =
                std::make_unique_for_overwrite<SubsetPrivateKeySwitchingKey<lvl21param>>();
            subprivkskgen<lvl21param>(*ptrsubprivksk<lvl21param>()[key], func, sk);
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
        if constexpr (std::is_same_v<P, AHlvl1param>) {
            ptrahk<AHlvl1param>() =
                std::make_unique_for_overwrite<AnnihilateKey<AHlvl1param>>();
            annihilatekeygen<AHlvl1param>(*ptrahk<AHlvl1param>(), sk);
        }
        else if constexpr (std::is_same_v<P, AHlvl2param>) {
            ptrahk<AHlvl2param>() =
                std::make_unique_for_overwrite<AnnihilateKey<AHlvl2param>>();
            annihilatekeygen<AHlvl2param>(*ptrahk<AHlvl2param>(), sk);
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }

    template <class P>
    void emplacecbsk(const SecretKey& sk)
    {
        if constexpr (std::is_same_v<P, AHlvl1param>) {
            ptrcbsk<AHlvl1param>() =
                std::make_unique_for_overwrite<CBswitchingKey<AHlvl1param>>();
            for (int i = 0; i < P::k; i++) {
                Polynomial<P> partkey;
                for (int j = 0; j < P::n; j++)
                    partkey[j] = -sk.key.get<P>()[i * P::n + j];
                (*ptrcbsk<AHlvl1param>())[i] =
                    trgswfftSymEncrypt<P>(partkey, sk.key.get<P>());
            }
        }
        else if constexpr (std::is_same_v<P, AHlvl2param>) {
            ptrcbsk<AHlvl2param>() =
                std::make_unique_for_overwrite<CBswitchingKey<AHlvl2param>>();
            for (int i = 0; i < P::k; i++) {
                Polynomial<P> partkey;
                for (int j = 0; j < P::n; j++)
                    partkey[j] = -sk.key.get<P>()[i * P::n + j];
                (*ptrcbsk<AHlvl2param>())[i] =
                    trgswfftSymEncrypt<P>(partkey, sk.key.get<P>());
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
            return *ptrbk<lvl01param>();
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return *ptrbk<lvlh1param>();
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *ptrbk<lvl02param>();
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return *ptrbk<lvlh2param>();
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    BootstrappingKeyFFT<P>& getbkfft() const
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return *ptrbkfft<lvl01param>();
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return *ptrbkfft<lvlh1param>();
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *ptrbkfft<lvl02param>();
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return *ptrbkfft<lvlh2param>();
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    BootstrappingKeyNTT<P>& getbkntt() const
    {
        if constexpr (std::is_same_v<P, lvl01param>) {
            return *ptrbkntt<lvl01param>();
        }
        else if constexpr (std::is_same_v<P, lvlh1param>) {
            return *ptrbkntt<lvlh1param>();
        }
        else if constexpr (std::is_same_v<P, lvl02param>) {
            return *ptrbkntt<lvl02param>();
        }
        else if constexpr (std::is_same_v<P, lvlh2param>) {
            return *ptrbkntt<lvlh2param>();
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    KeySwitchingKey<P>& getiksk() const
    {
        if constexpr (std::is_same_v<P, lvl10param>) {
            return *ptriksk<lvl10param>();
        }
        else if constexpr (std::is_same_v<P, lvl1hparam>) {
            return *ptriksk<lvl1hparam>();
        }
        else if constexpr (std::is_same_v<P, lvl20param>) {
            return *ptriksk<lvl20param>();
        }
        else if constexpr (std::is_same_v<P, lvl2hparam>) {
            return *ptriksk<lvl2hparam>();
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            return *ptriksk<lvl21param>();
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            return *ptriksk<lvl22param>();
        }
        else if constexpr (std::is_same_v<P, lvl31param>) {
            return *ptriksk<lvl31param>();
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    SubsetKeySwitchingKey<P>& getsubiksk() const
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return *ptrsubiksk<lvl21param>();
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    PrivateKeySwitchingKey<P>& getprivksk(const std::string& key) const
    {
        if constexpr (std::is_same_v<P, lvl11param>) {
            return *(ptrprivksk<lvl11param>().at(key));
        }
        else if constexpr (std::is_same_v<P, lvl21param>) {
            return *(ptrprivksk<lvl21param>().at(key));
        }
        else if constexpr (std::is_same_v<P, lvl22param>) {
            return *(ptrprivksk<lvl22param>().at(key));
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    SubsetPrivateKeySwitchingKey<P>& getsubprivksk(const std::string& key) const
    {
        if constexpr (std::is_same_v<P, lvl21param>) {
            return *(ptrsubprivksk<lvl21param>().at(key));
        }
        else
            static_assert(false_v<typename P::targetP::T>,
                          "Not predefined parameter!");
    }
    template <class P>
    AnnihilateKey<P>& getahk() const
    {
        if constexpr (std::is_same_v<P, AHlvl1param>) {
            return *ptrahk<AHlvl1param>();
        }
        else if constexpr (std::is_same_v<P, AHlvl2param>) {
            return *ptrahk<AHlvl2param>();
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
    template <class P>
    CBswitchingKey<P>& getcbsk() const
    {
        if constexpr (std::is_same_v<P, AHlvl1param>) {
            return *ptrcbsk<AHlvl1param>();
        }
        else if constexpr (std::is_same_v<P, AHlvl2param>) {
            return *ptrcbsk<AHlvl2param>();
        }
        else
            static_assert(false_v<typename P::T>, "Not predefined parameter!");
    }
};

}  // namespace TFHEpp