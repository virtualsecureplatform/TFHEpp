#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>

#include <tfhe++.hpp>

namespace {

using CKKS = TFHEpp::lvl5param;
using TFHE = TFHEpp::lvl2param;
constexpr std::uint32_t log_q = 128;
constexpr std::uint32_t log_delta = 88;

void fill_ckks_key(TFHEpp::Key<CKKS> &key)
{
    for (std::uint32_t i = 0; i < CKKS::n; i++) {
        if (i % 3 == 0)
            key[i] = CKKS::T{0} - CKKS::T{1};
        else
            key[i] = CKKS::T{i % 3 == 1 ? 0U : 1U};
    }
}

void fill_tfhe_key(TFHEpp::Key<TFHE> &key)
{
    for (std::uint32_t i = 0; i < TFHE::k * TFHE::n; i++)
        key[i] = static_cast<TFHE::T>(static_cast<int>(i % 3) - 1);
}

std::int64_t centred(std::uint64_t value)
{
    return static_cast<std::int64_t>(value);
}

}  // namespace

// This is deliberately separate from the compact regression test: it checks
// the seeded RLWE-to-padded-LWE path at TFHEpp's production-sized lvl5 CKKS
// ring degree (2^14) and a standard 64-bit TFHE level.
int main()
{
    auto ckks_key = std::make_unique<TFHEpp::Key<CKKS>>();
    TFHEpp::Key<TFHE> tfhe_key{};
    fill_ckks_key(*ckks_key);
    fill_tfhe_key(tfhe_key);

    auto unit_lhs = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    auto unit_rhs = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    auto unit_product = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    (*unit_lhs)[0] = CKKS::T{1};
    (*unit_rhs)[0] = CKKS::T{1};
    TFHEpp::ckks_scheme_switch_detail::exactActiveTorusMulBySmall<CKKS,
                                                                    log_q>(
        *unit_product, *unit_lhs, *unit_rhs);
    if (static_cast<std::uint64_t>((*unit_product)[0]) != 1) {
        std::cerr << "full-size exact NTT unit product failed" << std::endl;
        return 1;
    }

    std::array<double, CKKS::n> values{};
    values.fill(0.0);
    values[0] = 0.25;
    auto input = std::make_unique<
        TFHEpp::CKKSCiphertext<CKKS, log_q, log_delta>>();
    // Build this full-ring fixture directly to isolate the cross-key
    // invariant; public encryption is checked independently below.
    auto encoded = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    encoded->fill(CKKS::T{0});
    (*encoded)[0] = TFHEpp::ckksEncodeCoeff<CKKS, log_q, log_delta>(values[0]);
    if (std::abs(TFHEpp::ckksDecodeCoeff<CKKS, log_q, log_delta>(
                     (*encoded)[0]) - values[0]) > 1e-12) {
        std::cerr << "full-size multi-limb CKKS encoding failed" << std::endl;
        return 1;
    }
    const auto negative_encoded =
        TFHEpp::ckksEncodeCoeff<CKKS, log_q, log_delta>(-0.5);
    if (std::abs(TFHEpp::ckksDecodeCoeff<CKKS, log_q, log_delta>(
                     negative_encoded) + 0.5) > 1e-12) {
        std::cerr << "full-size signed multi-limb CKKS encoding failed"
                  << std::endl;
        return 1;
    }
    auto source_key_poly = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    auto mask_phase = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    for (std::uint32_t i = 0; i < CKKS::n; i++) {
        input->ct[0][i] = TFHEpp::ckks_detail::uniformAtLevel<CKKS, log_q>();
        (*source_key_poly)[i] = (*ckks_key)[i];
    }
    TFHEpp::ckks_scheme_switch_detail::exactActiveTorusMulBySmall<CKKS,
                                                                    log_q>(
        *mask_phase, input->ct[0], *source_key_poly);
    for (std::uint32_t i = 0; i < CKKS::n; i++)
        input->ct[1][i] = TFHEpp::ckks_detail::reduceToLevel<CKKS, log_q>(
            (*encoded)[i] + (*mask_phase)[i]);
    constexpr auto output_scale =
        TFHEpp::CKKSToTFHEOutputScaleBits<CKKS, TFHE, log_q, log_delta>;
    static_assert(output_scale == 24);
    const auto expected = static_cast<std::int64_t>(
        std::llround(values[0] * std::ldexp(1.0, output_scale)));

    auto switch_key = std::make_unique<
        TFHEpp::CKKSToTFHERingSwitchKey<CKKS, TFHE, log_q>>();
    TFHEpp::CKKSToTFHERingSwitchKeyGen<CKKS, TFHE, log_q>(
        *switch_key, *ckks_key, tfhe_key, {0.0, 0});

    std::vector<TFHEpp::TLWE<TFHE>> output;
    TFHEpp::CKKSCoeffsToTFHEViaRingSwitch<CKKS, TFHE, log_q, log_delta>(
        output, *input, 1, *switch_key);
    const auto phase = centred(TFHEpp::tlweSymPhase<TFHE>(output[0], tfhe_key));
    if (std::llabs(phase - expected) > 64) {
        std::cerr << "full-size compact CKKS-to-TFHE mismatch got=" << phase
                  << " expected=" << expected << std::endl;
        return 1;
    }

    // Exercise the public encryption entry point as well as the isolated
    // exact fixture above.
    auto public_input = std::make_unique<
        TFHEpp::CKKSCiphertext<CKKS, log_q, log_delta>>();
    TFHEpp::ckksEncrypt<CKKS, log_q, log_delta>(
        *public_input, values, *ckks_key, {0.0, 0});
    TFHEpp::CKKSCoeffsToTFHEViaRingSwitch<CKKS, TFHE, log_q, log_delta>(
        output, *public_input, 1, *switch_key);
    const auto public_phase =
        centred(TFHEpp::tlweSymPhase<TFHE>(output[0], tfhe_key));
    if (std::llabs(public_phase - expected) > 64) {
        std::cerr << "public full-size compact CKKS-to-TFHE mismatch got="
                  << public_phase << " expected=" << expected << std::endl;
        return 1;
    }

    // Exercise the reverse core at the same production-sized ring, including
    // one nonzero LWE-mask coordinate against the encrypted TFHE secret.
    constexpr std::uint32_t reverse_plain_log_delta = 16;
    using ReverseCt = TFHEpp::CKKSPlainMulResult<
        CKKS, log_q, log_delta, reverse_plain_log_delta>;
    auto reverse_switch_key = std::make_unique<
        TFHEpp::TFHEToCKKSSwitchKey<CKKS, TFHE, log_q, log_delta>>();
    TFHEpp::TFHEToCKKSKeyGen<CKKS, TFHE, log_q, log_delta>(
        *reverse_switch_key, *ckks_key, tfhe_key, {0.0, 0});
    TFHEpp::TLWE<TFHE> reverse_input{};
    reverse_input[0] = std::uint64_t{1} << 62;
    reverse_input[TFHE::k * TFHE::n] = std::uint64_t{1} << 62;
    auto reverse_phase = std::make_unique<ReverseCt>();
    TFHEpp::TFHEToCKKSPhase<CKKS, TFHE, log_q, log_delta,
                             reverse_plain_log_delta>(
        *reverse_phase, reverse_input, *reverse_switch_key);
    auto reverse_key_poly = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    auto reverse_mask_phase = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    auto reverse_plain = std::make_unique<TFHEpp::Polynomial<CKKS>>();
    for (std::uint32_t i = 0; i < CKKS::n; i++)
        (*reverse_key_poly)[i] = (*ckks_key)[i];
    TFHEpp::ckks_detail::exactActiveTorusMulBySmall<CKKS, ReverseCt::log_q>(
        *reverse_mask_phase, reverse_phase->ct[0], *reverse_key_poly);
    for (std::uint32_t i = 0; i < CKKS::n; i++)
        (*reverse_plain)[i] = TFHEpp::ckks_detail::reduceToLevel<
            CKKS, ReverseCt::log_q>(
            reverse_phase->ct[CKKS::k][i] - (*reverse_mask_phase)[i]);
    const auto reverse_value =
        TFHEpp::ckksDecodeCoeff<CKKS, ReverseCt::log_q, log_delta>(
            (*reverse_plain)[0]);
    if (std::abs(reverse_value - 0.5) > 1e-12) {
        std::cerr << "full-size TFHE-to-CKKS phase mismatch got="
                  << reverse_value << " expected=0.5" << std::endl;
        return 1;
    }

    std::cout << "Passed" << std::endl;
}
