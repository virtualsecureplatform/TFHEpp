#pragma once

#include <cstdint>

#include "cloudkey.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"

namespace TFHEpp {

template <class midP, std::uint32_t l, std::uint32_t Bgbit>
constexpr Polynomial<midP> CBtestvector()
{
    Polynomial<midP> poly;
    constexpr uint32_t bitwidth = bits_needed<l - 1>();
    for (int i = 0; i < (midP::n >> bitwidth); i++)
        for (int j = 0; j < (1 << bitwidth); j++)
            poly[(i << bitwidth) + j] =
                1ULL << (std::numeric_limits<typename midP::T>::digits -
                         (j + 1) * Bgbit - 1);
    return poly;
}

template <class midP, class targetP>
constexpr Polynomial<midP> CBtestvector()
{
    return CBtestvector<midP, targetP::l, targetP::Bgbit>();
}

template <class P>
constexpr Polynomial<typename P::domainP> CBtestvector()
{
    return CBtestvector<typename P::domainP, typename P::targetP>();
}

template <class bkP, class privksP>
void CircuitBootstrapping(TRGSW<typename privksP::targetP>& trgsw,
                          const TLWE<typename bkP::domainP>& tlwe,
                          const EvalKey& ek)
{
    using targetP = typename privksP::targetP;
    static_assert(targetP::l̅ == 1 && targetP::l̅ₐ == 1,
                  "Circuit bootstrapping does not yet produce DD TRGSW rows");

    if constexpr (targetP::l == targetP::lₐ &&
                  targetP::Bgbit == targetP::Bgₐbit) {
        alignas(64) std::array<TLWE<typename bkP::targetP>, targetP::l> temp;
        GateBootstrappingManyLUT<bkP, targetP::l>(
            temp, tlwe, ek.getbkfft<bkP>(), CBtestvector<privksP>());
        for (int i = 0; i < targetP::l; i++) {
            temp[i][privksP::domainP::k * privksP::domainP::n] +=
                1ULL << (numeric_limits<typename privksP::domainP::T>::digits -
                         (i + 1) * targetP::Bgbit - 1);
            for (int k = 0; k < targetP::k + 1; k++)
                PrivKeySwitch<privksP>(
                    trgsw[i + k * targetP::l], temp[i],
                    ek.getprivksk<privksP>("privksk4cb_" +
                                           std::to_string(k)));
        }
    }
    else {
        alignas(64) std::array<TLWE<typename bkP::targetP>, targetP::lₐ>
            nonce_temp;
        GateBootstrappingManyLUT<bkP, targetP::lₐ>(
            nonce_temp, tlwe, ek.getbkfft<bkP>(),
            CBtestvector<typename privksP::domainP, targetP::lₐ,
                         targetP::Bgₐbit>());
        for (int i = 0; i < targetP::lₐ; i++) {
            nonce_temp[i][privksP::domainP::k * privksP::domainP::n] +=
                1ULL << (numeric_limits<typename privksP::domainP::T>::digits -
                         (i + 1) * targetP::Bgₐbit - 1);
            for (int k = 0; k < targetP::k; k++)
                PrivKeySwitch<privksP>(
                    trgsw[i + k * targetP::lₐ], nonce_temp[i],
                    ek.getprivksk<privksP>("privksk4cb_" +
                                           std::to_string(k)));
        }

        alignas(64) std::array<TLWE<typename bkP::targetP>, targetP::l>
            body_temp;
        GateBootstrappingManyLUT<bkP, targetP::l>(
            body_temp, tlwe, ek.getbkfft<bkP>(),
            CBtestvector<typename privksP::domainP, targetP::l,
                         targetP::Bgbit>());
        for (int i = 0; i < targetP::l; i++) {
            body_temp[i][privksP::domainP::k * privksP::domainP::n] +=
                1ULL << (numeric_limits<typename privksP::domainP::T>::digits -
                         (i + 1) * targetP::Bgbit - 1);
            PrivKeySwitch<privksP>(
                trgsw[i + targetP::k * targetP::lₐ], body_temp[i],
                ek.getprivksk<privksP>("privksk4cb_" +
                                       std::to_string(targetP::k)));
        }
    }
}

// Circuit bootstrap into the subset secret-key chain.  The blind rotation
// produces level-2 TLWEs under the subset key.  Reduce each one to the
// level-1 subset key once, then use subset private key switching to construct
// the TRGSW rows under that same key.
template <class bkP, class privksP>
void CircuitBootstrappingSubset(
    TRGSW<typename privksP::targetP>& trgsw,
    const TLWE<typename bkP::domainP>& tlwe, const EvalKey& ek)
{
    using targetP = typename privksP::targetP;
    using midP = typename bkP::targetP;
    static_assert(std::is_same_v<midP, typename privksP::domainP>);
    static_assert(targetP::l̅ == 1 && targetP::l̅ₐ == 1,
                  "Circuit bootstrapping does not yet produce DD TRGSW rows");

    auto subset_private_switch = [&](auto& out, const auto& in,
                                     const std::string& key) {
        TLWE<targetP> reduced;
        SubsetIdentityKeySwitch<privksP>(
            reduced, in, ek.getsubiksk<privksP>());
        SubsetPrivKeySwitch<privksP>(
            out, reduced, ek.getsubprivksk<privksP>(key));
    };

    if constexpr (targetP::l == targetP::lₐ &&
                  targetP::Bgbit == targetP::Bgₐbit) {
        alignas(64) std::array<TLWE<midP>, targetP::l> temp;
        GateBootstrappingManyLUT<bkP, targetP::l>(
            temp, tlwe, ek.getbkfft<bkP>(), CBtestvector<privksP>());
        for (int i = 0; i < targetP::l; i++) {
            temp[i][midP::k * midP::n] +=
                1ULL << (numeric_limits<typename midP::T>::digits -
                         (i + 1) * targetP::Bgbit - 1);
            TLWE<targetP> reduced;
            SubsetIdentityKeySwitch<privksP>(
                reduced, temp[i], ek.getsubiksk<privksP>());
            for (int k = 0; k < targetP::k + 1; k++)
                SubsetPrivKeySwitch<privksP>(
                    trgsw[i + k * targetP::l], reduced,
                    ek.getsubprivksk<privksP>(
                        "subprivksk4cb_" + std::to_string(k)));
        }
    }
    else {
        alignas(64) std::array<TLWE<midP>, targetP::lₐ> nonce_temp;
        GateBootstrappingManyLUT<bkP, targetP::lₐ>(
            nonce_temp, tlwe, ek.getbkfft<bkP>(),
            CBtestvector<midP, targetP::lₐ, targetP::Bgₐbit>());
        for (int i = 0; i < targetP::lₐ; i++) {
            nonce_temp[i][midP::k * midP::n] +=
                1ULL << (numeric_limits<typename midP::T>::digits -
                         (i + 1) * targetP::Bgₐbit - 1);
            TLWE<targetP> reduced;
            SubsetIdentityKeySwitch<privksP>(
                reduced, nonce_temp[i], ek.getsubiksk<privksP>());
            for (int k = 0; k < targetP::k; k++)
                SubsetPrivKeySwitch<privksP>(
                    trgsw[i + k * targetP::lₐ], reduced,
                    ek.getsubprivksk<privksP>(
                        "subprivksk4cb_" + std::to_string(k)));
        }

        alignas(64) std::array<TLWE<midP>, targetP::l> body_temp;
        GateBootstrappingManyLUT<bkP, targetP::l>(
            body_temp, tlwe, ek.getbkfft<bkP>(),
            CBtestvector<midP, targetP::l, targetP::Bgbit>());
        for (int i = 0; i < targetP::l; i++) {
            body_temp[i][midP::k * midP::n] +=
                1ULL << (numeric_limits<typename midP::T>::digits -
                         (i + 1) * targetP::Bgbit - 1);
            subset_private_switch(
                trgsw[i + targetP::k * targetP::lₐ], body_temp[i],
                "subprivksk4cb_" + std::to_string(targetP::k));
        }
    }
}

// https://eprint.iacr.org/2024/1318
template <class brP, class ahP>
void AnnihilateCircuitBootstrapping(TRGSW<typename brP::targetP>& trgsw,
                                    const TLWE<typename brP::domainP>& tlwe,
                                    const EvalKey& ek)
{
    static_assert(brP::targetP::k == ahP::k,
                  "brP::targetP::k must be equal to ahP::k");
    alignas(64) std::array<TLWE<typename brP::targetP>, brP::targetP::l> temp;
    GateBootstrappingManyLUT<brP, brP::targetP::l>(
        temp, tlwe, ek.getbkfft<brP>(),
        CBtestvector<typename brP::targetP, typename brP::targetP>());
    for (int i = 0; i < brP::targetP::l; i++) {
        temp[i][brP::targetP::k * brP::targetP::n] +=
            1ULL << (numeric_limits<typename brP::targetP::T>::digits -
                     (i + 1) * brP::targetP::Bgbit - 1);
        TRLWE<typename brP::targetP> temptrlwe;
        InvSampleExtractIndex<typename brP::targetP>(temptrlwe, temp[i], 0);
        AnnihilateKeySwitching<ahP>(
            trgsw[i + brP::targetP::k * brP::targetP::l], temptrlwe,
            ek.getahk<ahP>());
        // Scheme Switching
        for (int k = 0; k < brP::targetP::k; k++)
            ExternalProduct<ahP>(trgsw[i + k * brP::targetP::l],
                                 trgsw[i + brP::targetP::k * brP::targetP::l],
                                 ek.getcbsk<ahP>()[k]);
    }
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrapping(TRGSW<typename privksP::targetP>& trgsw,
                          const TLWE<typename iksP::domainP>& tlwe,
                          const EvalKey& ek)
{
    TLWE<typename bkP::domainP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek.getiksk<iksP>());
    CircuitBootstrapping<bkP, privksP>(trgsw, tlwelvl0, ek);
}

template <class iksP, class brP, class ahP>
void AnnihilateCircuitBootstrapping(TRGSW<typename brP::targetP>& trgsw,
                                    const TLWE<typename iksP::domainP>& tlwe,
                                    const EvalKey& ek)
{
    TLWE<typename brP::domainP> tlwelvl0;
    EvalIdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek);
    AnnihilateCircuitBootstrapping<brP, ahP>(trgsw, tlwelvl0, ek);
}

template <class brP, class privksP>
void CircuitBootstrapping(TRGSWFFT<typename privksP::targetP>& trgswfft,
                          const TLWE<typename brP::domainP>& tlwe,
                          const EvalKey& ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<brP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
}

template <class brP, class privksP>
void CircuitBootstrappingSubset(
    TRGSWFFT<typename privksP::targetP>& trgswfft,
    const TLWE<typename brP::domainP>& tlwe, const EvalKey& ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrappingSubset<brP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrapping(TRGSWFFT<typename privksP::targetP>& trgswfft,
                          const TLWE<typename iksP::domainP>& tlwe,
                          const EvalKey& ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<iksP, bkP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
}

template <class brP, class ahP>
void AnnihilateCircuitBootstrapping(TRGSWFFT<typename brP::targetP>& trgswfft,
                                    const TLWE<typename brP::domainP>& tlwe,
                                    const EvalKey& ek)
{
    alignas(64) TRGSW<typename brP::targetP> trgsw;
    AnnihilateCircuitBootstrapping<brP, ahP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename brP::targetP>(trgswfft, trgsw);
}

template <class iksP, class brP, class ahP>
void AnnihilateCircuitBootstrapping(TRGSWFFT<typename brP::targetP>& trgswfft,
                                    const TLWE<typename iksP::domainP>& tlwe,
                                    const EvalKey& ek)
{
    alignas(64) TRGSW<typename brP::targetP> trgsw;
    AnnihilateCircuitBootstrapping<iksP, brP, ahP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename brP::targetP>(trgswfft, trgsw);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingSub(TRGSW<typename privksP::targetP>& trgsw,
                             const TLWE<typename iksP::domainP>& tlwe,
                             const EvalKey& ek)
{
    alignas(64) TLWE<typename bkP::domainP> tlwelvl0;
    SubsetIdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek.getsubiksk<iksP>());
    CircuitBootstrapping<bkP, privksP>(trgsw, tlwelvl0, ek);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingSub(TRGSWFFT<typename privksP::targetP>& trgswfft,
                             const TLWE<typename iksP::domainP>& tlwe,
                             const EvalKey& ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrappingSub<iksP, bkP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
}

template <class brP, class privksP>
void CircuitBootstrappingInv(TRGSWFFT<typename privksP::targetP>& invtrgswfft,
                             const TLWE<typename brP::domainP>& tlwe,
                             const EvalKey& ek)
{
    alignas(64) TLWE<typename brP::domainP> invtlwe;
    // HomNot
    for (int i = 0; i <= brP::domainP::k * brP::domainP::n; i++)
        invtlwe[i] = -tlwe[i];
    CircuitBootstrapping<brP, privksP>(invtrgswfft, invtlwe, ek);
}

template <class brP, class privksP>
void CircuitBootstrappingSubsetInv(
    TRGSWFFT<typename privksP::targetP>& invtrgswfft,
    const TLWE<typename brP::domainP>& tlwe, const EvalKey& ek)
{
    alignas(64) TLWE<typename brP::domainP> invtlwe;
    // HomNot
    for (int i = 0; i <= brP::domainP::k * brP::domainP::n; i++)
        invtlwe[i] = -tlwe[i];
    CircuitBootstrappingSubset<brP, privksP>(invtrgswfft, invtlwe, ek);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingInv(TRGSWFFT<typename privksP::targetP>& invtrgswfft,
                             const TLWE<typename iksP::domainP>& tlwe,
                             const EvalKey& ek)
{
    alignas(64) TLWE<typename iksP::domainP> invtlwe;
    // HomNot
    for (int i = 0; i <= iksP::domainP::k * iksP::domainP::n; i++)
        invtlwe[i] = -tlwe[i];
    CircuitBootstrapping<iksP, bkP, privksP>(invtrgswfft, invtlwe, ek);
}

template <class brP, class privksP>
void CircuitBootstrappingWithInv(
    TRGSWFFT<typename privksP::targetP>& trgswfft,
    TRGSWFFT<typename privksP::targetP>& invtrgswfft,
    const TLWE<typename brP::domainP>& tlwe, const EvalKey& ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<brP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
    for (auto& row : trgsw)
        for (auto& polynomial : row)
            for (auto& coefficient : polynomial) coefficient *= -1;
    trgswhoneadd<typename privksP::targetP>(trgsw);
    ApplyFFT2trgsw<typename privksP::targetP>(invtrgswfft, trgsw);
}

template <class brP, class privksP>
void CircuitBootstrappingSubsetWithInv(
    TRGSWFFT<typename privksP::targetP>& trgswfft,
    TRGSWFFT<typename privksP::targetP>& invtrgswfft,
    const TLWE<typename brP::domainP>& tlwe, const EvalKey& ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrappingSubset<brP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
    for (auto& row : trgsw)
        for (auto& polynomial : row)
            for (auto& coefficient : polynomial) coefficient *= -1;
    trgswhoneadd<typename privksP::targetP>(trgsw);
    ApplyFFT2trgsw<typename privksP::targetP>(invtrgswfft, trgsw);
}

template <class iksP, class bkP, class privksP>
void CircuitBootstrappingWithInv(
    TRGSWFFT<typename privksP::targetP>& trgswfft,
    TRGSWFFT<typename privksP::targetP>& invtrgswfft,
    const TLWE<typename iksP::domainP>& tlwe, const EvalKey& ek)
{
    alignas(64) TRGSW<typename privksP::targetP> trgsw;
    CircuitBootstrapping<iksP, bkP, privksP>(trgsw, tlwe, ek);
    ApplyFFT2trgsw<typename privksP::targetP>(trgswfft, trgsw);
    for (auto& row : trgsw)
        for (auto& polynomial : row)
            for (auto& coefficient : polynomial) coefficient *= -1;
    trgswhoneadd<typename privksP::targetP>(trgsw);
    ApplyFFT2trgsw<typename privksP::targetP>(invtrgswfft, trgsw);
}

}  // namespace TFHEpp
