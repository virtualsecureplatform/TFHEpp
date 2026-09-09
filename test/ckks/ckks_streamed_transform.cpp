#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>
#include <string_view>
#include "scheme_switch_test_params.hpp"

// A test-only automorphism provider for transparent ciphertexts. This checks
// the full-size fused arithmetic without generating hundreds of MiB of keys.
namespace TFHEpp {
struct TransparentTestRotations {};
template <class P, std::uint32_t LogQ>
void CKKSRotateSlots(TRLWE<P> &out, const TRLWE<P> &in, int steps,
                     const TransparentTestRotations &)
{
    const int half = static_cast<int>(P::n / 2);
    steps = ((steps % half) + half) % half;
    std::uint64_t exponent = 1;
    for (int i = 0; i < steps; i++) exponent = exponent * 5 % (2 * P::n);
    for (std::uint32_t c = 0; c <= P::k; c++)
        Automorphism<P>(out[c], in[c], static_cast<int>(exponent));
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(out);
}
}

namespace {
struct P : TFHEpp::lvl5param {
    static constexpr std::uint32_t nbit = 4, n = 1 << nbit;
};
constexpr std::uint32_t q = 140, scale = 40, plain_scale = 30;
using Input = TFHEpp::CKKSCiphertext<P, q, scale>;
using Output = TFHEpp::CKKSPlainMulResult<P, q, scale, plain_scale>;
using Slots = TFHEpp::CKKSSlotVector<P>;

void require(bool ok, const char *message)
{
    if (!ok) throw std::runtime_error(message);
}
void full_arithmetic()
{
    using Large = scheme_switch_test::FullCKKS;
    constexpr std::uint32_t q = 1000, scale = 58, plain_scale = 60;
    using Out = TFHEpp::CKKSPlainMulResult<Large, q, scale, plain_scale>;
    auto input = std::make_unique<TFHEpp::CKKSSlotVector<Large>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<Large>>();
    for (std::uint32_t i = 0; i < Large::n / 2; i++)
        (*input)[i] = {0.1 * (static_cast<int>(i % 8) - 4), 0.03 * (i % 7)};
    auto ct = std::make_unique<TFHEpp::CKKSCiphertext<Large, q, scale>>();
    TFHEpp::ckksSlotEncode<Large, q, scale>(ct->ct[Large::k], *input);
    const std::vector<int> offsets{0, 3, 35};
    auto fill = [](int offset, TFHEpp::CKKSSlotVector<Large> &diagonal) {
        for (std::uint32_t i = 0; i < Large::n / 2; i++)
            diagonal[i] = {0.01 * (i % 5 + offset), -0.02};
    };
    auto diagonal = std::make_unique<TFHEpp::CKKSSlotVector<Large>>();
    for (int offset : offsets) {
        fill(offset, *diagonal);
        for (std::uint32_t i = 0; i < Large::n / 2; i++)
            (*expected)[i] += (*diagonal)[i] * (*input)[(i + offset) % (Large::n / 2)];
    }
    auto out = std::make_unique<Out>();
    const TFHEpp::TransparentTestRotations rotations;
    TFHEpp::CKKSLinearTransformBSGSStreamed<Large, q, scale, plain_scale>(
        *out, *ct, offsets, 32, fill, rotations, rotations);
    auto zero_key = std::make_unique<TFHEpp::Key<Large>>();
    auto got = std::make_unique<TFHEpp::CKKSSlotVector<Large>>();
    TFHEpp::ckksSlotDecrypt<Large, Out::log_q, scale>(*got, *out, *zero_key);
    double max_error = 0;
    for (std::uint32_t i = 0; i < Large::n / 2; i++) {
        const double error = std::abs((*got)[i] - (*expected)[i]);
        require(std::isfinite(error), "non-finite full-size transform output");
        max_error = std::max(max_error, error);
    }
    std::cout << "Full-size fused transform max error=" << max_error << std::endl;
    require(max_error < 1e-9, "full-size fused transform failed");
}

}

int main(int argc, char **argv)
try {
    if (argc == 2 && std::string_view(argv[1]) == "--full-arithmetic") {
        full_arithmetic();
        std::cout << "Passed" << std::endl;
        return 0;
    }
    if (argc != 1) throw std::invalid_argument("usage: ckks_streamed_transform [--full-arithmetic]");
    TFHEpp::Key<P> key{};
    TFHEpp::keyGen<P>(key);
    const TFHEpp::CKKSNoise noise{0.0, 4};
    // Non-dividing BSGS step exercises the partial final giant group.
    for (int step : {1, 3, 8}) {
        TFHEpp::CKKSRotationKeyIndexSet<P> input_indices{}, output_indices{};
        TFHEpp::CKKSSlotsToCoeffsRotationKeyIndices<P, q, scale, plain_scale>(
            input_indices, output_indices, step);
        TFHEpp::CKKSRealLinearTransformPlan<P, q, scale, plain_scale> reference;
        TFHEpp::CKKSBuildPackedSlotToCoeffPlan<P, q, scale, plain_scale>(reference, step);
        TFHEpp::CKKSRotationKeyIndexSet<P> ref_input{}, ref_output{};
        TFHEpp::CKKSMarkConjugationKeyIndex<P>(ref_input);
        TFHEpp::CKKSCollectLinearTransformPlanRotationKeyIndices<P, q, scale, plain_scale>(
            ref_input, ref_output, reference.direct);
        TFHEpp::CKKSCollectLinearTransformPlanRotationKeyIndices<P, q, scale, plain_scale>(
            ref_input, ref_output, reference.conjugate);
        require(input_indices == ref_input && output_indices == ref_output,
                "streamed and materialized schedules differ");
        auto input_keys = std::make_unique<TFHEpp::CKKSSeededSparseGaloisKey<P, q>>();
        auto output_keys = std::make_unique<
            TFHEpp::CKKSSeededSparseGaloisKey<P, q - plain_scale>>();
        TFHEpp::CKKSSlotsToCoeffsSparseGaloisKeyGen<P, q, scale, plain_scale>(
            *input_keys, *output_keys, key, step, noise);
        Slots input{};
        for (std::uint32_t i = 0; i < P::n / 2; i++)
            input[i] = {0.1 * (static_cast<int>(i) - 4), 0.03 * (i + 1)};
        Input ct{};
        TFHEpp::ckksSlotEncrypt<P, q, scale>(ct, input, key, noise);
        // Arbitrary complex diagonals, unsorted and negative offsets, and
        // an analytically computed matrix-vector product as the oracle.
        const std::vector<int> offsets{3, -1, 0};
        Slots expected{};
        const auto fill = [](int offset, Slots &diag) {
            for (std::uint32_t i = 0; i < P::n / 2; i++)
                diag[i] = {0.07 * (i + offset + 1), -0.02 * (offset + 1)};
        };
        for (int offset : offsets) {
            offset = (offset + P::n / 2) % (P::n / 2);
            Slots diag{};
            fill(offset, diag);
            for (std::uint32_t i = 0; i < P::n / 2; i++)
                expected[i] += diag[i] * input[(i + offset) % (P::n / 2)];
        }
        Output result{};
        TFHEpp::CKKSLinearTransformBSGSStreamed<P, q, scale, plain_scale>(
            result, ct, offsets, step, fill, *input_keys, *output_keys);
        Slots got{};
        TFHEpp::ckksSlotDecrypt<P, Output::log_q, scale>(got, result, key);
        for (std::uint32_t i = 0; i < P::n / 2; i++) {
            require(std::isfinite(std::abs(got[i])) &&
                        std::abs(got[i] - expected[i]) < 1e-4,
                    "streamed matrix-vector product mismatch");
        }
        bool rejected = false;
        try {
            TFHEpp::CKKSLinearTransformBSGSStreamed<P, q, scale, plain_scale>(
                result, ct, {0, static_cast<int>(P::n / 2)}, step, fill,
                *input_keys, *output_keys);
        }
        catch (const std::invalid_argument &) { rejected = true; }
        require(rejected, "duplicate normalized offset accepted");

        for (std::uint32_t logical_slots : {1U, 2U, 4U, 8U}) {
            for (std::uint32_t i = 0; i < P::n / 2; i++)
                input[i] = {0.1 * (static_cast<int>(i % logical_slots) - 2),
                            0.03 * (i % logical_slots + 1)};
            TFHEpp::ckksSlotEncrypt<P, q, scale>(ct, input, key, noise);
            // Provision only the logical sparse schedule, not the full set.
            TFHEpp::CKKSSlotsToCoeffsSparseGaloisKeyGen<P, q, scale, plain_scale>(
                *input_keys, *output_keys, key, step, noise, logical_slots);
            TFHEpp::CKKSSlotsToCoeffs<P, q, scale, plain_scale>(
                result, ct, step, *input_keys, *output_keys, logical_slots);
            std::array<double, P::n> coeffs{};
            TFHEpp::ckksDecrypt<P>(coeffs, result, key);
            const auto gap = P::n / (2 * logical_slots);
            for (std::uint32_t i = 0; i < P::n; i++) {
                double want = 0;
                if (i % gap == 0)
                    want = i < P::n / 2 ? input[i / gap].real()
                                            : input[(i - P::n / 2) / gap].imag();
                require(std::isfinite(coeffs[i]) && std::abs(coeffs[i] - want) < 1e-4,
                        "sparse SlotToCoeff layout mismatch");
            }
        }
    }
    std::cout << "Passed" << std::endl;
}
catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
}
