#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <tfhe++.hpp>

namespace {

// A compact, self-contained parameter pair keeps this regression test small;
// production CKKS->TFHE keys are intentionally much larger.
struct ToyCKKSParam {
    static constexpr int32_t key_value_min = -1;
    static constexpr int32_t key_value_max = 1;
    static constexpr std::uint32_t nbit = 4;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    using T = __uint128_t;
    static constexpr double α = 0.0;
    static constexpr std::uint32_t l = 1;
    static constexpr std::uint32_t lₐ = l;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = Bg;
    static constexpr T μ = T{1} << 125;
    static constexpr std::uint32_t l̅ = 8;
    static constexpr std::uint32_t l̅ₐ = l̅;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = B̅gbit;
};

struct ToyTFHEParam {
    static constexpr int32_t key_value_min = -1;
    static constexpr int32_t key_value_max = 1;
    static constexpr std::uint32_t n = 4;
    static constexpr std::uint32_t k = 1;
    using T = std::uint64_t;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static constexpr double α = 0.0;
};

// This deliberately does not divide the eight real slots in EvalCKKSParam.
// It exercises the rectangular FHEW-to-CKKS layout that has historically been
// easy to get wrong when the LWE dimension and packed batch do not align.
struct NonDivTFHEParam {
    static constexpr int32_t key_value_min = -1;
    static constexpr int32_t key_value_max = 1;
    static constexpr std::uint32_t n = 3;
    static constexpr std::uint32_t k = 1;
    using T = std::uint64_t;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static constexpr double α = 0.0;
};

// EvalMod uses CKKS ciphertext multiplication, which needs an FFT-capable
// torus rather than the __uint128_t-only compact forward-switch fixture.
struct EvalCKKSParam {
    static constexpr int32_t key_value_min = -1;
    static constexpr int32_t key_value_max = 1;
    static constexpr std::uint32_t nbit = 4;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    using T = TFHEpp::MultiLimbUInt<5>;
    static constexpr double α = 0.0;
    static constexpr std::uint32_t l = 1;
    static constexpr std::uint32_t lₐ = l;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bgₐ = Bg;
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr std::uint32_t l̅ = 20;
    static constexpr std::uint32_t l̅ₐ = l̅;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = B̅gbit;
};

using CKKS = ToyCKKSParam;
using TFHE = ToyTFHEParam;
using PackedCKKS = TFHEpp::ckkslvl3param;
constexpr std::uint32_t log_q = 60;
constexpr std::uint32_t log_delta = 20;
constexpr std::uint32_t basebit = 2;
constexpr std::uint32_t levels = 30;
using Ct = TFHEpp::CKKSCiphertext<CKKS, log_q, log_delta>;
using SwitchKey =
    TFHEpp::CKKSToTFHEKeySwitchKey<CKKS, TFHE, basebit, levels>;

void fill_key(TFHEpp::Key<CKKS> &key)
{
    for (std::uint32_t i = 0; i < CKKS::n; i++)
        key[i] = static_cast<CKKS::T>(static_cast<int>(i % 3) - 1);
}

void fill_key(TFHEpp::Key<TFHE> &key)
{
    for (std::uint32_t i = 0; i < TFHE::n; i++)
        key[i] = static_cast<TFHE::T>(static_cast<int>(i % 3) - 1);
}

template <class P>
void fill_ckks_key(TFHEpp::Key<P> &key)
{
    for (std::uint32_t i = 0; i < P::n; i++)
        key[i] = static_cast<typename P::T>(static_cast<int>(i % 3) - 1);
}

std::int64_t centred(std::uint64_t value)
{
    return static_cast<std::int64_t>(value);
}

}  // namespace

int main()
{
    TFHEpp::Key<CKKS> ckks_key{};
    TFHEpp::Key<TFHE> tfhe_key{};
    fill_key(ckks_key);
    fill_key(tfhe_key);

    std::array<double, CKKS::n> values{};
    values.fill(0.0);
    values[0] = 0.25;
    values[3] = -0.5;

    Ct ct{};
    TFHEpp::ckksEncrypt<CKKS, log_q, log_delta>(ct, values, ckks_key,
                                                  {0.0, 0});

    SwitchKey switch_key{};
    TFHEpp::CKKSToTFHEKeyGen<CKKS, TFHE, basebit, levels>(
        switch_key, ckks_key, tfhe_key);

    constexpr auto output_scale =
        TFHEpp::CKKSToTFHEOutputScaleBits<CKKS, TFHE, log_q, log_delta>;
    static_assert(output_scale == 24);
    for (const std::uint32_t index : {0U, 3U}) {
        TFHEpp::TLWE<TFHE> output{};
        TFHEpp::CKKSCoeffToTFHE<CKKS, TFHE, log_q, log_delta, basebit,
                                 levels>(output, ct, index, switch_key);
        const auto phase =
            centred(TFHEpp::tlweSymPhase<TFHE>(output, tfhe_key));
        const auto expected = static_cast<std::int64_t>(
            std::llround(values[index] * std::ldexp(1.0, output_scale)));
        if (std::llabs(phase - expected) > 2) {
            std::cerr << "coefficient " << index << " got=" << phase
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    // The compact forward path follows OpenFHE's intermediate-RLWE switch:
    // switch the CKKS ring secret once to a padded TFHE secret, then extract
    // TLWEs without a per-coefficient TLWE switching key.
    TFHEpp::Key<EvalCKKSParam> ring_ckks_key{};
    fill_ckks_key<EvalCKKSParam>(ring_ckks_key);
    std::array<double, EvalCKKSParam::n> ring_values{};
    ring_values.fill(0.0);
    ring_values[0] = values[0];
    ring_values[3] = values[3];
    TFHEpp::CKKSCiphertext<EvalCKKSParam, log_q, log_delta> ring_ct{};
    TFHEpp::ckksEncrypt<EvalCKKSParam, log_q, log_delta>(
        ring_ct, ring_values, ring_ckks_key, {0.0, 0});
    TFHEpp::CKKSToTFHERingSwitchKey<EvalCKKSParam, TFHE, log_q>
        ring_switch_key{};
    TFHEpp::CKKSToTFHERingSwitchKeyGen<EvalCKKSParam, TFHE, log_q>(
        ring_switch_key, ring_ckks_key, tfhe_key, {0.0, 0});
    constexpr auto ring_output_scale =
        TFHEpp::CKKSToTFHEOutputScaleBits<EvalCKKSParam, TFHE, log_q,
                                           log_delta>;
    for (const std::uint32_t index : {0U, 3U}) {
        TFHEpp::TLWE<TFHE> output{};
        TFHEpp::CKKSCoeffToTFHEViaRingSwitch<EvalCKKSParam, TFHE, log_q,
                                              log_delta>(
            output, ring_ct, index, ring_switch_key);
        const auto phase =
            centred(TFHEpp::tlweSymPhase<TFHE>(output, tfhe_key));
        const auto expected = static_cast<std::int64_t>(
            std::llround(ring_values[index] *
                         std::ldexp(1.0, ring_output_scale)));
        if (std::llabs(phase - expected) > 2) {
            std::cerr << "compact coefficient " << index << " got=" << phase
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }
    std::vector<TFHEpp::TLWE<TFHE>> compact_outputs;
    TFHEpp::CKKSCoeffsToTFHEViaRingSwitch<EvalCKKSParam, TFHE, log_q,
                                           log_delta>(compact_outputs, ring_ct,
                                                      4, ring_switch_key);
    for (const std::uint32_t index : {0U, 3U}) {
        const auto phase = centred(
            TFHEpp::tlweSymPhase<TFHE>(compact_outputs[index], tfhe_key));
        const auto expected = static_cast<std::int64_t>(
            std::llround(ring_values[index] *
                         std::ldexp(1.0, ring_output_scale)));
        if (std::llabs(phase - expected) > 2) {
            std::cerr << "batched compact coefficient " << index
                      << " got=" << phase << " expected=" << expected
                      << std::endl;
            return 1;
        }
    }
    const std::vector<std::uint32_t> compact_indices{3, 0};
    std::vector<TFHEpp::TLWE<TFHE>> selected_compact_outputs;
    TFHEpp::CKKSCoeffIndicesToTFHEViaRingSwitch<
        EvalCKKSParam, TFHE, log_q, log_delta>(
        selected_compact_outputs, ring_ct, compact_indices, ring_switch_key);
    for (std::size_t i = 0; i < compact_indices.size(); i++) {
        const auto phase = centred(TFHEpp::tlweSymPhase<TFHE>(
            selected_compact_outputs[i], tfhe_key));
        const auto expected = static_cast<std::int64_t>(std::llround(
            ring_values[compact_indices[i]] *
            std::ldexp(1.0, ring_output_scale)));
        if (std::llabs(phase - expected) > 2) {
            std::cerr << "selected compact coefficient " << compact_indices[i]
                      << " got=" << phase << " expected=" << expected
                      << std::endl;
            return 1;
        }
    }

    // Full forward switch: packed CKKS slots are decoded to coefficients,
    // sample-extracted, modulus-lifted, and key-switched into TFHE TLWEs.
    constexpr std::uint32_t forward_log_q = 70;
    constexpr std::uint32_t forward_log_delta = 12;
    constexpr std::uint32_t forward_plain_log_delta = 10;
    constexpr std::uint32_t forward_decoded_log_q =
        forward_log_q - forward_plain_log_delta;
    constexpr std::uint32_t forward_basebit = 2;
    constexpr std::uint32_t forward_levels = 30;
    using ForwardCt = TFHEpp::CKKSCiphertext<EvalCKKSParam, forward_log_q,
                                             forward_log_delta>;
    using ForwardDecoded = TFHEpp::CKKSPlainMulResult<
        EvalCKKSParam, forward_log_q, forward_log_delta,
        forward_plain_log_delta>;
    TFHEpp::Key<EvalCKKSParam> forward_ckks_key{};
    fill_ckks_key<EvalCKKSParam>(forward_ckks_key);
    TFHEpp::CKKSSlotVector<EvalCKKSParam> forward_slots{};
    forward_slots.fill({0.0, 0.0});
    forward_slots[0] = {0.125, 0.0};
    forward_slots[1] = {-0.0625, 0.0};
    ForwardCt forward_ct{};
    TFHEpp::ckksSlotEncrypt<EvalCKKSParam, forward_log_q, forward_log_delta>(
        forward_ct, forward_slots, forward_ckks_key, {0.0, 0});
    auto forward_input_gk = std::make_unique<
        TFHEpp::CKKSHybridSparseGaloisKey<EvalCKKSParam, forward_log_q>>();
    auto forward_output_gk = std::make_unique<
        TFHEpp::CKKSHybridSparseGaloisKey<EvalCKKSParam,
                                           forward_decoded_log_q>>();
    TFHEpp::CKKSSlotsToCoeffsHybridSparseGaloisKeyGen<
        EvalCKKSParam, forward_log_q, forward_log_delta,
        forward_plain_log_delta>(*forward_input_gk, *forward_output_gk,
                                  forward_ckks_key, 8, {0.0, 0});
    using ForwardSwitchKey = TFHEpp::CKKSToTFHEKeySwitchKey<
        EvalCKKSParam, TFHE, forward_basebit, forward_levels>;
    ForwardSwitchKey forward_switch_key{};
    TFHEpp::CKKSToTFHEKeyGen<EvalCKKSParam, TFHE, forward_basebit,
                              forward_levels>(forward_switch_key,
                                              forward_ckks_key, tfhe_key);
    ForwardDecoded forward_decoded{};
    TFHEpp::CKKSSlotsToCoeffs<EvalCKKSParam, forward_log_q,
                              forward_log_delta, forward_plain_log_delta>(
        forward_decoded, forward_ct, 8, *forward_input_gk,
        *forward_output_gk);
    std::array<double, EvalCKKSParam::n> forward_plain{};
    TFHEpp::ckksDecrypt<EvalCKKSParam>(forward_plain, forward_decoded,
                                        forward_ckks_key);
    std::vector<TFHEpp::TLWE<TFHE>> forward_lwes;
    TFHEpp::CKKSSlotsToTFHE<
        EvalCKKSParam, TFHE, forward_log_q, forward_log_delta,
        forward_plain_log_delta, forward_basebit, forward_levels>(
        forward_lwes, forward_ct, 4, 8, *forward_input_gk, *forward_output_gk,
        forward_switch_key);
    constexpr auto forward_output_scale =
        TFHEpp::CKKSToTFHEOutputScaleBits<
            EvalCKKSParam, TFHE, forward_decoded_log_q, forward_log_delta>;
    static_assert(forward_output_scale == 16);
    for (std::uint32_t i = 0; i < forward_lwes.size(); i++) {
        const auto got = centred(TFHEpp::tlweSymPhase<TFHE>(forward_lwes[i],
                                                            tfhe_key));
        const auto expected = static_cast<std::int64_t>(std::llround(
            forward_plain[i] * std::ldexp(1.0, forward_output_scale)));
        if (std::llabs(got - expected) > 8) {
            std::cerr << "packed CKKS-to-TFHE mismatch coefficient=" << i
                      << " got=" << got << " expected=" << expected
                      << std::endl;
            return 1;
        }
    }
    TFHEpp::CKKSToTFHERingSwitchKey<EvalCKKSParam, TFHE,
                                     forward_decoded_log_q>
        forward_ring_switch_key{};
    TFHEpp::CKKSToTFHERingSwitchKeyGen<EvalCKKSParam, TFHE,
                                        forward_decoded_log_q>(
        forward_ring_switch_key, forward_ckks_key, tfhe_key, {0.0, 0});
    std::vector<TFHEpp::TLWE<TFHE>> forward_compact_lwes;
    TFHEpp::CKKSSlotsToTFHEViaRingSwitch<
        EvalCKKSParam, TFHE, forward_log_q, forward_log_delta,
        forward_plain_log_delta>(
        forward_compact_lwes, forward_ct, 4, 8, *forward_input_gk,
        *forward_output_gk, forward_ring_switch_key);
    for (std::uint32_t i = 0; i < forward_compact_lwes.size(); i++) {
        const auto got = centred(TFHEpp::tlweSymPhase<TFHE>(
            forward_compact_lwes[i], tfhe_key));
        const auto expected = static_cast<std::int64_t>(std::llround(
            forward_plain[i] * std::ldexp(1.0, forward_output_scale)));
        if (std::llabs(got - expected) > 8) {
            std::cerr << "compact packed CKKS-to-TFHE mismatch coefficient="
                      << i << " got=" << got << " expected=" << expected
                      << std::endl;
            return 1;
        }
    }

    // Match OpenFHE's sparse packed-slot extraction gap.  With four logical
    // slots in this N=16 ring, the decoded coefficient indices are 0, 2, 4,
    // and 6 rather than consecutive coefficients.
    std::vector<TFHEpp::TLWE<TFHE>> forward_spaced_direct;
    std::vector<TFHEpp::TLWE<TFHE>> forward_spaced_compact;
    TFHEpp::CKKSSlotsToTFHEOpenFHEPacked<
        EvalCKKSParam, TFHE, forward_log_q, forward_log_delta,
        forward_plain_log_delta, forward_basebit, forward_levels>(
        forward_spaced_direct, forward_ct, 2, 4, 8, *forward_input_gk,
        *forward_output_gk, forward_switch_key);
    TFHEpp::CKKSSlotsToTFHEViaRingSwitchOpenFHEPacked<
        EvalCKKSParam, TFHE, forward_log_q, forward_log_delta,
        forward_plain_log_delta>(
        forward_spaced_compact, forward_ct, 2, 4, 8, *forward_input_gk,
        *forward_output_gk, forward_ring_switch_key);
    for (std::uint32_t i = 0; i < 2; i++) {
        const std::uint32_t coefficient = 2 * i;
        const auto expected = static_cast<std::int64_t>(std::llround(
            forward_plain[coefficient] *
            std::ldexp(1.0, forward_output_scale)));
        const auto direct = centred(
            TFHEpp::tlweSymPhase<TFHE>(forward_spaced_direct[i], tfhe_key));
        const auto compact = centred(TFHEpp::tlweSymPhase<TFHE>(
            forward_spaced_compact[i], tfhe_key));
        if (std::llabs(direct - expected) > 8 ||
            std::llabs(compact - expected) > 8) {
            std::cerr << "OpenFHE-gap CKKS-to-TFHE mismatch coefficient="
                      << coefficient << " direct=" << direct
                      << " compact=" << compact << " expected=" << expected
                      << std::endl;
            return 1;
        }
    }

    // OpenFHE performs a second modulus switch after SlotToCoeff.  Exercise
    // the matching TFHEpp extraction-level API with a distinct Q'.
    constexpr std::uint32_t forward_extraction_log_q = 56;
    constexpr auto forward_reduced_output_scale =
        TFHEpp::CKKSToTFHEOutputScaleBits<
            EvalCKKSParam, TFHE, forward_extraction_log_q,
            forward_log_delta>;
    static_assert(forward_reduced_output_scale == 20);
    std::vector<TFHEpp::TLWE<TFHE>> forward_reduced_lwes;
    TFHEpp::CKKSSlotsToTFHEAtLevel<
        EvalCKKSParam, TFHE, forward_log_q, forward_log_delta,
        forward_plain_log_delta, forward_extraction_log_q, forward_basebit,
        forward_levels>(
        forward_reduced_lwes, forward_ct, 4, 8, *forward_input_gk,
        *forward_output_gk, forward_switch_key);
    TFHEpp::CKKSToTFHERingSwitchKey<EvalCKKSParam, TFHE,
                                     forward_extraction_log_q>
        forward_reduced_ring_switch_key{};
    TFHEpp::CKKSToTFHERingSwitchKeyGen<EvalCKKSParam, TFHE,
                                        forward_extraction_log_q>(
        forward_reduced_ring_switch_key, forward_ckks_key, tfhe_key,
        {0.0, 0});
    std::vector<TFHEpp::TLWE<TFHE>> forward_reduced_compact_lwes;
    TFHEpp::CKKSSlotsToTFHEViaRingSwitchAtLevel<
        EvalCKKSParam, TFHE, forward_log_q, forward_log_delta,
        forward_plain_log_delta, forward_extraction_log_q>(
        forward_reduced_compact_lwes, forward_ct, 4, 8, *forward_input_gk,
        *forward_output_gk, forward_reduced_ring_switch_key);
    for (std::uint32_t i = 0; i < forward_reduced_lwes.size(); i++) {
        const auto expected = static_cast<std::int64_t>(std::llround(
            forward_plain[i] *
            std::ldexp(1.0, forward_reduced_output_scale)));
        const auto direct = centred(TFHEpp::tlweSymPhase<TFHE>(
            forward_reduced_lwes[i], tfhe_key));
        const auto compact = centred(TFHEpp::tlweSymPhase<TFHE>(
            forward_reduced_compact_lwes[i], tfhe_key));
        if (std::llabs(direct - expected) > 8 ||
            std::llabs(compact - expected) > 8) {
            std::cerr << "reduced-level CKKS-to-TFHE mismatch coefficient="
                      << i << " direct=" << direct << " compact=" << compact
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    // The actual EvalCKKStoFHEW composition uses Q' and the sparse extraction
    // stride together.  Verify both direct and compact combined wrappers.
    std::vector<TFHEpp::TLWE<TFHE>> forward_reduced_spaced_direct;
    std::vector<TFHEpp::TLWE<TFHE>> forward_reduced_spaced_compact;
    TFHEpp::CKKSSlotsToTFHEAtLevelOpenFHEPacked<
        EvalCKKSParam, TFHE, forward_log_q, forward_log_delta,
        forward_plain_log_delta, forward_extraction_log_q, forward_basebit,
        forward_levels>(
        forward_reduced_spaced_direct, forward_ct, 2, 4, 8,
        *forward_input_gk, *forward_output_gk, forward_switch_key);
    TFHEpp::CKKSSlotsToTFHEViaRingSwitchAtLevelOpenFHEPacked<
        EvalCKKSParam, TFHE, forward_log_q, forward_log_delta,
        forward_plain_log_delta, forward_extraction_log_q>(
        forward_reduced_spaced_compact, forward_ct, 2, 4, 8,
        *forward_input_gk, *forward_output_gk,
        forward_reduced_ring_switch_key);
    for (std::uint32_t i = 0; i < 2; i++) {
        const std::uint32_t coefficient = 2 * i;
        const auto expected = static_cast<std::int64_t>(std::llround(
            forward_plain[coefficient] *
            std::ldexp(1.0, forward_reduced_output_scale)));
        const auto direct = centred(TFHEpp::tlweSymPhase<TFHE>(
            forward_reduced_spaced_direct[i], tfhe_key));
        const auto compact = centred(TFHEpp::tlweSymPhase<TFHE>(
            forward_reduced_spaced_compact[i], tfhe_key));
        if (std::llabs(direct - expected) > 8 ||
            std::llabs(compact - expected) > 8) {
            std::cerr << "Q'-level OpenFHE-gap CKKS-to-TFHE mismatch "
                      << "coefficient=" << coefficient << " direct="
                      << direct << " compact=" << compact
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    TFHEpp::Key<PackedCKKS> packed_ckks_key{};
    fill_ckks_key<PackedCKKS>(packed_ckks_key);

    // Reverse-direction core: encrypt the TFHE secret under CKKS and form the
    // unwrapped b-a*s phase.  Its value is checked modulo one, because the
    // periodic EvalMod stage is deliberately separate from this primitive.
    using ReverseCt =
        TFHEpp::CKKSPlainMulResult<PackedCKKS, log_q, log_delta, log_delta>;
    TFHEpp::TFHEToCKKSSwitchKey<PackedCKKS, TFHE, log_q, log_delta>
        reverse_switch_key{};
    TFHEpp::TFHEToCKKSKeyGen<PackedCKKS, TFHE, log_q, log_delta>(
        reverse_switch_key, packed_ckks_key, tfhe_key, {0.0, 0});
    TFHEpp::TLWE<TFHE> lwe_input{};
    TFHEpp::tlweSymEncrypt<TFHE>(lwe_input, std::uint64_t{1} << 60, 0.0,
                                  tfhe_key);
    ReverseCt reverse_phase{};
    TFHEpp::TFHEToCKKSPhase<PackedCKKS, TFHE, log_q, log_delta, log_delta>(
        reverse_phase, lwe_input, reverse_switch_key);
    std::array<double, PackedCKKS::n> reverse_plain{};
    TFHEpp::ckksDecrypt<PackedCKKS>(reverse_plain, reverse_phase,
                                    packed_ckks_key);
    const double expected_phase = std::ldexp(
        static_cast<double>(TFHEpp::tlweSymPhase<TFHE>(lwe_input, tfhe_key)),
        -64);
    if (std::abs(std::remainder(reverse_plain[0] - expected_phase, 1.0)) >
        2e-4) {
        std::cerr << "TFHE-to-CKKS phase mismatch got=" << reverse_plain[0]
                  << " expected(mod 1)=" << expected_phase << std::endl;
        return 1;
    }

    // Complete reverse switch: normalize the encrypted phase and apply the
    // bounded-cosine periodic reduction.  The toy ring keeps this regression
    // small while still exercising the same EvalMod composition used by a
    // production CKKS parameter set.
    constexpr std::uint32_t eval_log_q = 150;
    constexpr std::uint32_t eval_log_delta = 15;
    constexpr std::uint32_t eval_phase_plain_log_delta = 10;
    constexpr std::uint32_t eval_coeff_log_delta = 10;
    constexpr std::size_t eval_degree = 15;
    constexpr std::uint32_t eval_double_angle = 2;
    constexpr std::uint32_t eval_phase_bound = TFHE::n + 1;
    constexpr std::uint32_t eval_log_message_ratio = 2;
    constexpr std::uint32_t eval_phase_log_q =
        eval_log_q - eval_phase_plain_log_delta;
    using EvalKey = TFHEpp::TFHEToCKKSPackedEvalModKey<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta, eval_phase_plain_log_delta,
        eval_coeff_log_delta, eval_degree, eval_double_angle>;
    using EvalOut = TFHEpp::CKKSEvalModBoundedCosResult<
        EvalCKKSParam, eval_phase_log_q, eval_log_delta, eval_coeff_log_delta,
        eval_degree, eval_double_angle>;
    TFHEpp::Key<EvalCKKSParam> eval_ckks_key{};
    fill_ckks_key<EvalCKKSParam>(eval_ckks_key);
    auto eval_key = std::make_unique<EvalKey>();
    TFHEpp::TFHEToCKKSPackedEvalModKeyGen<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta, eval_phase_plain_log_delta,
        eval_coeff_log_delta, eval_degree, eval_double_angle>(
        *eval_key, eval_ckks_key, tfhe_key, eval_phase_bound,
        eval_log_message_ratio, {0.0, 0});
    auto eval_input_gk =
        std::make_unique<
            TFHEpp::CKKSSparseGaloisKey<EvalCKKSParam, eval_log_q>>();
    auto eval_output_gk =
        std::make_unique<
            TFHEpp::CKKSSparseGaloisKey<EvalCKKSParam, eval_phase_log_q>>();
    TFHEpp::TFHEToCKKSPackedPhaseSparseGaloisKeyGen<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta>(*eval_input_gk, *eval_output_gk,
                                    eval_ckks_key, {0.0, 0});
    if (TFHEpp::CKKSRotationKeyIndexSetCount<EvalCKKSParam>(
            eval_input_gk->available) != 2 ||
        TFHEpp::CKKSRotationKeyIndexSetCount<EvalCKKSParam>(
            eval_output_gk->available) != 0) {
        std::cerr << "packed TFHE-to-CKKS sparse rotation-key schedule "
                     "is not minimal"
                  << std::endl;
        return 1;
    }
    std::vector<TFHEpp::TLWE<TFHE>> eval_lwes(2);
    TFHEpp::tlweSymEncrypt<TFHE>(eval_lwes[0], std::uint64_t{1} << 60, 0.0,
                                  tfhe_key);
    TFHEpp::tlweSymEncrypt<TFHE>(eval_lwes[1],
                                  std::uint64_t{0} - (std::uint64_t{1} << 59),
                                  0.0, tfhe_key);

    std::vector<TFHEpp::TLWE<TFHE>> dense_lwes(EvalCKKSParam::n / 2);
    for (std::uint32_t i = 0; i < dense_lwes.size(); i++) {
        const std::int64_t message =
            static_cast<std::int64_t>(static_cast<int>(i % 5) - 2) *
            (std::int64_t{1} << 59);
        TFHEpp::tlweSymEncrypt<TFHE>(
            dense_lwes[i], static_cast<std::uint64_t>(message), 0.0,
            tfhe_key);
    }

    // OpenFHE-style dense rectangular partial decryption packs phases into
    // consecutive slots from one encrypted copy of the TFHE secret.
    TFHEpp::TFHEToCKKSSlotPackedSwitchKey<EvalCKKSParam, TFHE, eval_log_q,
                                          eval_log_delta>
        dense_switch_key{};
    TFHEpp::TFHEToCKKSSlotPackedKeyGen<EvalCKKSParam, TFHE, eval_log_q,
                                       eval_log_delta>(
        dense_switch_key, eval_ckks_key, tfhe_key, {0.0, 0});
    auto dense_input_gk = std::make_unique<
        TFHEpp::CKKSSparseGaloisKey<EvalCKKSParam, eval_log_q>>();
    auto dense_output_gk = std::make_unique<
        TFHEpp::CKKSSparseGaloisKey<EvalCKKSParam, eval_phase_log_q>>();
    TFHEpp::TFHEToCKKSSlotPackedPhaseSparseGaloisKeyGen<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta>(*dense_input_gk, *dense_output_gk,
                                    eval_ckks_key, dense_lwes.size(), 8,
                                    {0.0, 0});
    using DensePhase = TFHEpp::CKKSPlainMulResult<
        EvalCKKSParam, eval_log_q, eval_log_delta, eval_phase_plain_log_delta>;
    DensePhase dense_phase{};
    TFHEpp::TFHEToCKKSSlotPackedPhase<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta>(
        dense_phase, dense_lwes, dense_switch_key, 8, *dense_input_gk,
        *dense_output_gk);
    TFHEpp::CKKSSlotVector<EvalCKKSParam> dense_plain{};
    TFHEpp::ckksSlotDecrypt<EvalCKKSParam, DensePhase::log_q,
                            DensePhase::log_delta>(dense_plain, dense_phase,
                                                   eval_ckks_key);
    for (std::uint32_t i = 0; i < dense_lwes.size(); i++) {
        double expected =
            TFHEpp::ckks_scheme_switch_detail::tfheTorusToUnitInterval<TFHE>(
                dense_lwes[i][TFHE::n]);
        for (std::uint32_t j = 0; j < TFHE::n; j++)
            expected -=
                TFHEpp::ckks_scheme_switch_detail::tfheTorusToUnitInterval<
                    TFHE>(dense_lwes[i][j]) *
                TFHEpp::ckks_scheme_switch_detail::tfheSecretCoefficientToDouble<
                    TFHE>(tfhe_key[j]);
        if (std::abs(dense_plain[i].real() - expected) > 0.05) {
            std::cerr << "dense TFHE-to-CKKS phase mismatch slot=" << i
                      << " got=" << dense_plain[i].real()
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    // A non-dividing LWE dimension and a non-power-of-two batch must still
    // map each public mask to its own consecutive CKKS slot.  This mirrors a
    // boundary case in OpenFHE's rectangular partial-decryption path.
    TFHEpp::Key<NonDivTFHEParam> nondiv_tfhe_key{};
    for (std::uint32_t i = 0; i < NonDivTFHEParam::n; i++)
        nondiv_tfhe_key[i] = static_cast<NonDivTFHEParam::T>(
            static_cast<int>(i % 3) - 1);
    std::vector<TFHEpp::TLWE<NonDivTFHEParam>> nondiv_lwes(5);
    for (std::uint32_t i = 0; i < nondiv_lwes.size(); i++) {
        const std::int64_t message =
            static_cast<std::int64_t>(static_cast<int>(i) - 2) *
            (std::int64_t{1} << 59);
        TFHEpp::tlweSymEncrypt<NonDivTFHEParam>(
            nondiv_lwes[i], static_cast<std::uint64_t>(message), 0.0,
            nondiv_tfhe_key);
    }
    TFHEpp::TFHEToCKKSSlotPackedSwitchKey<EvalCKKSParam, NonDivTFHEParam,
                                          eval_log_q, eval_log_delta>
        nondiv_switch_key{};
    TFHEpp::TFHEToCKKSSlotPackedKeyGen<
        EvalCKKSParam, NonDivTFHEParam, eval_log_q, eval_log_delta>(
        nondiv_switch_key, eval_ckks_key, nondiv_tfhe_key, {0.0, 0});
    auto nondiv_input_gk = std::make_unique<
        TFHEpp::CKKSSparseGaloisKey<EvalCKKSParam, eval_log_q>>();
    auto nondiv_output_gk = std::make_unique<
        TFHEpp::CKKSSparseGaloisKey<EvalCKKSParam, eval_phase_log_q>>();
    TFHEpp::TFHEToCKKSSlotPackedPhaseSparseGaloisKeyGen<
        EvalCKKSParam, NonDivTFHEParam, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta>(*nondiv_input_gk, *nondiv_output_gk,
                                    eval_ckks_key, nondiv_lwes.size(), 4,
                                    {0.0, 0});
    DensePhase nondiv_phase{};
    TFHEpp::TFHEToCKKSSlotPackedPhase<
        EvalCKKSParam, NonDivTFHEParam, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta>(nondiv_phase, nondiv_lwes,
                                    nondiv_switch_key, 4, *nondiv_input_gk,
                                    *nondiv_output_gk);
    TFHEpp::CKKSSlotVector<EvalCKKSParam> nondiv_plain{};
    TFHEpp::ckksSlotDecrypt<EvalCKKSParam, DensePhase::log_q,
                            DensePhase::log_delta>(nondiv_plain, nondiv_phase,
                                                   eval_ckks_key);
    for (std::uint32_t i = 0; i < nondiv_lwes.size(); i++) {
        double expected = TFHEpp::ckks_scheme_switch_detail::
            tfheTorusToUnitInterval<NonDivTFHEParam>(
                nondiv_lwes[i][NonDivTFHEParam::n]);
        for (std::uint32_t j = 0; j < NonDivTFHEParam::n; j++)
            expected -= TFHEpp::ckks_scheme_switch_detail::
                tfheTorusToUnitInterval<NonDivTFHEParam>(nondiv_lwes[i][j]) *
                TFHEpp::ckks_scheme_switch_detail::
                    tfheSecretCoefficientToDouble<NonDivTFHEParam>(
                        nondiv_tfhe_key[j]);
        if (std::abs(nondiv_plain[i].real() - expected) > 0.05) {
            std::cerr << "non-dividing dense TFHE-to-CKKS phase mismatch "
                      << "slot=" << i << " got=" << nondiv_plain[i].real()
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    using DenseEvalKey = TFHEpp::TFHEToCKKSSlotPackedEvalModKey<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta, eval_coeff_log_delta, eval_degree,
        eval_double_angle>;
    auto dense_eval_key = std::make_unique<DenseEvalKey>();
    TFHEpp::TFHEToCKKSSlotPackedEvalModKeyGen<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta, eval_coeff_log_delta, eval_degree,
        eval_double_angle>(*dense_eval_key, eval_ckks_key, tfhe_key,
                           eval_phase_bound, eval_log_message_ratio,
                           {0.0, 0});
    EvalOut dense_eval_out{};
    TFHEpp::TFHEToCKKSSlotPackedEvalMod<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta, eval_coeff_log_delta, eval_degree,
        eval_double_angle>(dense_eval_out, dense_lwes, *dense_eval_key, 8,
                           *dense_input_gk, *dense_output_gk);
    TFHEpp::CKKSSlotVector<EvalCKKSParam> dense_eval_plain{};
    TFHEpp::ckksSlotDecrypt<EvalCKKSParam, EvalOut::log_q,
                            EvalOut::log_delta>(dense_eval_plain,
                                                dense_eval_out, eval_ckks_key);
    for (std::uint32_t i = 0; i < dense_lwes.size(); i++) {
        const double expected = std::remainder(
            std::ldexp(
                static_cast<double>(TFHEpp::tlweSymPhase<TFHE>(dense_lwes[i],
                                                                tfhe_key)),
                -64),
            1.0);
        if (std::abs(dense_eval_plain[i].real() - expected) > 0.12) {
            std::cerr << "dense TFHE-to-CKKS EvalMod mismatch slot=" << i
                      << " got=" << dense_eval_plain[i].real()
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    // The public post-scale/post-bias maps the canonical residue into the
    // application-facing CKKS range, matching EvalFHEWtoCKKS's last stage.
    constexpr std::uint32_t affine_plain_log_delta = 10;
    using DenseAffineOut = TFHEpp::CKKSPlainMulResult<
        EvalCKKSParam, EvalOut::log_q, eval_log_delta,
        affine_plain_log_delta>;
    DenseAffineOut dense_affine_out{};
    TFHEpp::TFHEToCKKSSlotPackedEvalModAffine<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta, eval_coeff_log_delta, eval_degree,
        eval_double_angle, affine_plain_log_delta>(
        dense_affine_out, dense_lwes, *dense_eval_key, 8, *dense_input_gk,
        *dense_output_gk, 2.0, 0.25);
    TFHEpp::CKKSSlotVector<EvalCKKSParam> dense_affine_plain{};
    TFHEpp::ckksSlotDecrypt<EvalCKKSParam, DenseAffineOut::log_q,
                            DenseAffineOut::log_delta>(
        dense_affine_plain, dense_affine_out, eval_ckks_key);
    for (std::uint32_t i = 0; i < dense_lwes.size(); i++) {
        const double phase = std::remainder(
            std::ldexp(static_cast<double>(TFHEpp::tlweSymPhase<TFHE>(
                           dense_lwes[i], tfhe_key)),
                       -64),
            1.0);
        const double expected = 2.0 * phase + 0.25;
        if (std::abs(dense_affine_plain[i].real() - expected) > 0.25) {
            std::cerr << "dense TFHE-to-CKKS affine mismatch slot=" << i
                      << " got=" << dense_affine_plain[i].real()
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    // Reproduce EvalFHEWtoCKKS's public p/pmin/pmax post-processing rule.
    // For p=8 and [-1, 1], OpenFHE applies 4*phase + 0.5.
    DenseAffineOut dense_openfhe_out{};
    TFHEpp::TFHEToCKKSSlotPackedEvalModOpenFHE<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta, eval_coeff_log_delta, eval_degree,
        eval_double_angle, affine_plain_log_delta>(
        dense_openfhe_out, dense_lwes, *dense_eval_key, 8,
        *dense_input_gk, *dense_output_gk, 8, -1.0, 1.0);
    TFHEpp::CKKSSlotVector<EvalCKKSParam> dense_openfhe_plain{};
    TFHEpp::ckksSlotDecrypt<EvalCKKSParam, DenseAffineOut::log_q,
                            DenseAffineOut::log_delta>(
        dense_openfhe_plain, dense_openfhe_out, eval_ckks_key);
    for (std::uint32_t i = 0; i < dense_lwes.size(); i++) {
        const double phase = std::remainder(
            std::ldexp(static_cast<double>(TFHEpp::tlweSymPhase<TFHE>(
                           dense_lwes[i], tfhe_key)),
                       -64),
            1.0);
        const double expected = 4.0 * phase + 0.5;
        if (std::abs(dense_openfhe_plain[i].real() - expected) > 0.5) {
            std::cerr << "OpenFHE TFHE-to-CKKS post-map mismatch slot=" << i
                      << " got=" << dense_openfhe_plain[i].real()
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    // Match OpenFHE's sparse-output convention: a smaller switched batch can
    // be repeated over the logical CKKS slot count after periodic reduction.
    auto dense_replication_gk = std::make_unique<
        TFHEpp::CKKSSparseGaloisKey<EvalCKKSParam, EvalOut::log_q>>();
    TFHEpp::CKKSSchemeSwitchReplicationSparseGaloisKeyGen<
        EvalCKKSParam, EvalOut::log_q>(*dense_replication_gk, eval_ckks_key,
                                       2, EvalCKKSParam::n / 2, {0.0, 0});
    if (TFHEpp::CKKSRotationKeyIndexSetCount<EvalCKKSParam>(
            dense_replication_gk->available) != 2) {
        std::cerr << "dense TFHE-to-CKKS replication rotation-key schedule "
                     "is not minimal"
                  << std::endl;
        return 1;
    }
    const std::vector<TFHEpp::TLWE<TFHE>> dense_subset{
        dense_lwes[0], dense_lwes[1]};
    EvalOut repeated_dense_eval_out{};
    TFHEpp::TFHEToCKKSSlotPackedEvalModToSlots<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta,
        eval_phase_plain_log_delta, eval_coeff_log_delta, eval_degree,
        eval_double_angle>(repeated_dense_eval_out, dense_subset,
                           *dense_eval_key, 8, EvalCKKSParam::n / 2,
                           *dense_input_gk, *dense_output_gk,
                           *dense_replication_gk);
    TFHEpp::CKKSSlotVector<EvalCKKSParam> repeated_dense_plain{};
    TFHEpp::ckksSlotDecrypt<EvalCKKSParam, EvalOut::log_q,
                            EvalOut::log_delta>(repeated_dense_plain,
                                                repeated_dense_eval_out,
                                                eval_ckks_key);
    for (std::uint32_t i = 0; i < EvalCKKSParam::n / 2; i++) {
        const auto &lwe = dense_subset[i % dense_subset.size()];
        const double expected = std::remainder(
            std::ldexp(static_cast<double>(TFHEpp::tlweSymPhase<TFHE>(
                           lwe, tfhe_key)),
                       -64),
            1.0);
        if (std::abs(repeated_dense_plain[i].real() - expected) > 0.12) {
            std::cerr << "dense TFHE-to-CKKS sparse output mismatch slot="
                      << i << " got=" << repeated_dense_plain[i].real()
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    EvalOut eval_out{};
    TFHEpp::TFHEToCKKSPackedEvalMod<
        EvalCKKSParam, TFHE, eval_log_q, eval_log_delta, eval_phase_plain_log_delta,
        eval_coeff_log_delta, eval_degree, eval_double_angle>(
        eval_out, eval_lwes, *eval_key, *eval_input_gk, *eval_output_gk);
    TFHEpp::CKKSSlotVector<EvalCKKSParam> eval_plain{};
    TFHEpp::ckksSlotDecrypt<EvalCKKSParam, EvalOut::log_q, EvalOut::log_delta>(
        eval_plain, eval_out, eval_ckks_key);
    for (std::uint32_t i = 0; i < eval_lwes.size(); i++) {
        const double expected = std::remainder(
            std::ldexp(
                static_cast<double>(TFHEpp::tlweSymPhase<TFHE>(eval_lwes[i],
                                                                tfhe_key)),
                -64),
            1.0);
        const std::uint32_t slot = i * TFHE::n;
        if (std::abs(eval_plain[slot].real() - expected) > 0.12) {
            std::cerr << "TFHE-to-CKKS EvalMod mismatch slot=" << i
                      << " got=" << eval_plain[slot].real()
                      << " expected=" << expected << std::endl;
            return 1;
        }
    }

    std::cout << "Passed" << std::endl;
}
