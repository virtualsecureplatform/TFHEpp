#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>
#include <memory>
#include <string_view>
#include <tfhe++.hpp>

// Functional test only: independently randomized keys and nonzero Gaussian
// CKKS noise. These parameter combinations are not security certifications.
namespace scheme_switch_test {

template <class P, class L, std::uint32_t reverse_q,
          std::uint32_t reverse_scale, std::uint32_t reverse_plain,
          std::size_t degree, std::uint32_t double_angle>
int run()
{
    constexpr std::uint32_t slots = 8;
    constexpr std::uint32_t forward_q = 86, forward_scale = 40, forward_plain = 40;
    using ReverseKey = TFHEpp::TFHEToCKKSSlotPackedEvalModKey<
        P, L, reverse_q, reverse_scale, reverse_plain, reverse_scale, degree,
        double_angle>;
    using Output = TFHEpp::CKKSEvalModBoundedCosResult<
        P, reverse_q - reverse_plain, reverse_scale, reverse_scale, degree,
        double_angle>;
    std::cout << "Reverse EvalMod key bytes=" << sizeof(ReverseKey) << std::endl;
    std::cout << "N=" << P::n << " LWE dimension=" << L::n
              << " logical slots=" << slots << std::endl;
    auto ckks_key = std::make_unique<TFHEpp::Key<P>>();
    TFHEpp::Key<L> tfhe_key{};
    TFHEpp::keyGen<P>(*ckks_key);
    TFHEpp::keyGen<L>(tfhe_key);
    // Centered Gaussian noise at every active level. Both test parameter
    // types supply a 3.2-unit floor; half-open four-bit uniform noise has a
    // -0.5 mean whose bias accumulates in large rotation schedules.
    const TFHEpp::CKKSNoise noise{std::ldexp(3.2, -static_cast<int>(reverse_q)), 0};
    const std::array<double, slots> values{-1.0, -0.75, -0.5, 0.0,
                                          0.25, 0.5, 0.75, 1.0};
    std::vector<TFHEpp::TLWE<L>> lwes;
    {
        auto packed = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
        for (std::uint32_t i = 0; i < P::n / 2; i++)
            (*packed)[i] = {values[i % slots], 0.0};
        auto ct = std::make_unique<
            TFHEpp::CKKSCiphertext<P, forward_q, forward_scale>>();
        TFHEpp::ckksSlotEncrypt<P, forward_q, forward_scale>(
            *ct, *packed, *ckks_key, noise);
        auto input_keys = std::make_unique<
            TFHEpp::CKKSSeededSparseGaloisKey<P, forward_q>>();
        auto output_keys = std::make_unique<
            TFHEpp::CKKSSeededSparseGaloisKey<P, forward_q - forward_plain>>();
        TFHEpp::CKKSSlotsToCoeffsSparseGaloisKeyGen<
            P, forward_q, forward_scale, forward_plain>(
            *input_keys, *output_keys, *ckks_key, 4, noise, slots);
        auto ring_key = std::make_unique<TFHEpp::CKKSToTFHERingSwitchKey<
            P, L, forward_q - forward_plain>>();
        TFHEpp::CKKSToTFHERingSwitchKeyGen<P, L, forward_q - forward_plain>(
            *ring_key, *ckks_key, tfhe_key, noise);
        std::cout << "Forward keys ready" << std::endl;
        TFHEpp::CKKSSparseSlotsToTFHEViaRingSwitch<
            P, L, forward_q, forward_scale, forward_plain>(
            lwes, *ct, slots, slots, 4, *input_keys, *output_keys, *ring_key);
    }
    double forward_error = 0;
    for (std::uint32_t i = 0; i < slots; i++) {
        const double phase = std::remainder(std::ldexp(static_cast<double>(
            TFHEpp::tlweSymPhase<L>(lwes[i], tfhe_key)), -std::numeric_limits<typename L::T>::digits), 1.0);
        forward_error = std::max(forward_error, std::abs(64 * phase - values[i]));
    }
    std::cout << "Forward max error=" << forward_error << std::endl;
    if (!std::isfinite(forward_error) || forward_error > 0.01)
        throw std::runtime_error("forward switch failed");

    auto reverse_key = std::make_unique<ReverseKey>();
    TFHEpp::TFHEToCKKSSlotPackedEvalModKeyGen<
        P, L, reverse_q, reverse_scale, reverse_plain, reverse_scale, degree,
        double_angle>(*reverse_key, *ckks_key, tfhe_key, L::k * L::n + 1, 0, noise);
    auto input_keys = std::make_unique<
        TFHEpp::CKKSSeededSparseGaloisKey<P, reverse_q>>();
    auto output_keys = std::make_unique<
        TFHEpp::CKKSSeededSparseGaloisKey<P, reverse_q - reverse_plain>>();
    TFHEpp::TFHEToCKKSSlotPackedPhaseSparseGaloisKeyGen<
        P, L, reverse_q, reverse_scale, reverse_plain>(
        *input_keys, *output_keys, *ckks_key, slots, 32, noise);
    std::cout << "Reverse keys ready" << std::endl;
    auto result = std::make_unique<Output>();
    using Phase = TFHEpp::CKKSPlainMulResult<P, reverse_q, reverse_scale, reverse_plain>;
    auto phase_ct = std::make_unique<Phase>();
    TFHEpp::TFHEToCKKSSlotPackedPhase<P, L, reverse_q, reverse_scale, reverse_plain>(
        *phase_ct, lwes, reverse_key->switch_key, 32, *input_keys, *output_keys,
        reverse_key->phase_scale);
    auto phase_plain = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Phase::log_q, Phase::log_delta>(
        *phase_plain, *phase_ct, *ckks_key);
    for (std::uint32_t i = 0; i < slots; i++) {
        double expected = TFHEpp::ckks_scheme_switch_detail::tfheTorusToUnitInterval<L>(lwes[i][L::n]);
        for (std::uint32_t j = 0; j < L::n; j++)
            expected -= TFHEpp::ckks_scheme_switch_detail::tfheTorusToUnitInterval<L>(lwes[i][j]) * static_cast<double>(tfhe_key[j]);
        expected *= reverse_key->phase_scale;
        const double error = std::abs((*phase_plain)[i] - expected);
        std::cout << "phase " << i << " error=" << error << std::endl;
        if (!std::isfinite(error) || error > 0.05 / (64 * (L::k * L::n + 1)))
            throw std::runtime_error("reverse partial decryption failed");
    }
    TFHEpp::CKKSEvalModBoundedCosNormalized<P, Phase::log_q, reverse_scale,
        reverse_scale, degree, double_angle>(*result, *phase_ct,
            reverse_key->polynomial, reverse_key->evalmod_key);
    // Public post-scale restores the original CKKS units while encrypted.
    TFHEpp::CKKSMulIntegerInPlace<P, Output::log_q, Output::log_delta>(*result, 64);
    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Output::log_q, Output::log_delta>(
        *decoded, *result, *ckks_key);
    double max_error = 0;
    for (std::uint32_t i = 0; i < slots; i++) {
        const auto got = (*decoded)[i];
        std::cout << "slot " << i << ": " << values[i] << " -> " << got << '\n';
        if (!std::isfinite(std::abs(got))) throw std::runtime_error("non-finite output");
        max_error = std::max(max_error, std::abs(got - values[i]));
    }
    std::cout << "Roundtrip max error=" << max_error << std::endl;
    if (max_error > 0.05) throw std::runtime_error("roundtrip accuracy failed");
    std::cout << "Passed" << std::endl;
    return 0;
}

}  // namespace scheme_switch_test
