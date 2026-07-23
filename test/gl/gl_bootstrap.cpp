#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>

namespace {

struct BootstrapTestBaseParameter {
    static constexpr std::int32_t key_value_max = 1;
    static constexpr std::int32_t key_value_min = -1;
    static constexpr std::uint32_t n = 16;
    static constexpr std::uint32_t k = 1;
    using T = __uint128_t;
    static constexpr std::uint32_t Bgbit = 8;
    static constexpr std::uint32_t B̅gbit = 8;
};

using GLP = TFHEpp::GLParameter<BootstrapTestBaseParameter, 2, 5, 20>;

struct BootstrapMultiLimbTestBaseParameter {
    static constexpr std::int32_t key_value_max = 1;
    static constexpr std::int32_t key_value_min = -1;
    static constexpr std::uint32_t n = 16;
    static constexpr std::uint32_t k = 1;
    using T = TFHEpp::MultiLimbUInt<4>;
    static constexpr std::uint32_t Bgbit = 8;
    static constexpr std::uint32_t B̅gbit = 8;
};

using MultiLimbGLP =
    TFHEpp::GLParameter<BootstrapMultiLimbTestBaseParameter, 2, 5, 34>;
using Schedule =
    TFHEpp::GLSHIPBootstrapSchedule<GLP, 24, 4, 8, 8, 90, 18, 3, 1, 2, 2, 8, 8>;

static_assert(Schedule::input_log_q == 40);
static_assert(Schedule::input_log_delta == 20);
static_assert(Schedule::tree_depth == 2);
static_assert(Schedule::output_log_q == 54);
static_assert(TFHEpp::GLSHIPPaperProfileFitsStorage<TFHEpp::GL256p17Parameter>);
static_assert(
    std::numeric_limits<typename TFHEpp::GL256p17Parameter::T>::digits >=
    TFHEpp::GLSHIPPaperParameterProfile<TFHEpp::GL256p17Parameter>::log_pq);

TFHEpp::GLMatrixBatch<GLP> makeInput()
{
    TFHEpp::GLMatrixBatch<GLP> input;
    for (std::uint32_t batch = 0; batch < GLP::phi; batch++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                input(batch, row, column) = {
                    0.025 * (1 + batch) + 0.015 * row - 0.01 * column,
                    -0.018 * (1 + batch) + 0.008 * row + 0.006 * column};
    return input;
}

double matrixError(const TFHEpp::GLMatrixBatch<GLP> &lhs,
                   const TFHEpp::GLMatrixBatch<GLP> &rhs)
{
    double error = 0;
    for (std::uint32_t batch = 0; batch < GLP::phi; batch++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                error = std::max(error, std::abs(lhs(batch, row, column) -
                                                 rhs(batch, row, column)));
    return error;
}

std::complex<double> sineRefresh(const std::complex<double> value)
{
    constexpr double pi = 3.14159265358979323846;
    constexpr double gap = 16.0;
    return {gap / (2 * pi) * std::sin(2 * pi * value.real() / gap),
            gap / (2 * pi) * std::sin(2 * pi * value.imag() / gap)};
}

}  // namespace

int main()
{
    std::cout << "GL SHIP bootstrap regression" << std::endl;

    TFHEpp::Key<BootstrapTestBaseParameter> dense_key{};
    dense_key[0] = 1;
    dense_key[1] = static_cast<BootstrapTestBaseParameter::T>(-1);
    dense_key[4] = 1;

    TFHEpp::Key<BootstrapTestBaseParameter> sparse_key{};
    sparse_key[0] = 1;
    sparse_key[3] = static_cast<BootstrapTestBaseParameter::T>(-1);
    sparse_key[10] = 1;
    const std::array<TFHEpp::GLSHIPSupportInterval,
                     Schedule::sparse_hamming_weight>
        intervals{{{0, 2}, {2, 3}, {9, 3}}};

    auto bootstrap_key =
        std::make_unique<TFHEpp::GLSHIPBootstrapKey<Schedule>>();
    std::cout << "  generating evaluation key" << std::endl;
    TFHEpp::GLSHIPBootstrapKeyGen(*bootstrap_key, dense_key, sparse_key,
                                  intervals);

    // Check the slice canonical encoder independently; masked-column tables
    // and the refreshed Gaussian channels both use this representation.
    TFHEpp::GLBaseSlotTable<GLP> base_values;
    for (std::uint32_t batch = 0; batch < GLP::phi; batch++)
        for (std::uint32_t x = 0; x < GLP::matrix_dimension; x++)
            base_values(batch, x) = std::complex<double>(
                0.07 * batch + 0.03 * x, -0.02 * batch + 0.01 * x);
    TFHEpp::GLBasePlaintext<GLP, 40, 20> base_plaintext;
    TFHEpp::GLBaseEncode(base_plaintext, base_values);
    TFHEpp::GLBaseSlotTable<GLP> base_roundtrip;
    TFHEpp::GLBaseDecode(base_roundtrip, base_plaintext);
    double base_error = 0;
    for (std::uint32_t batch = 0; batch < GLP::phi; batch++)
        for (std::uint32_t x = 0; x < GLP::matrix_dimension; x++)
            base_error = std::max(
                base_error,
                std::abs(base_roundtrip(batch, x) - base_values(batch, x)));
    std::cout << "  base encode/decode max error=" << base_error << std::endl;
    if (base_error > 2e-5) return 1;

    TFHEpp::GLBaseSlotTable<MultiLimbGLP> high_scale_values;
    for (std::uint32_t batch = 0; batch < MultiLimbGLP::phi; batch++)
        for (std::uint32_t x = 0; x < MultiLimbGLP::matrix_dimension; x++)
            high_scale_values(batch, x) = std::complex<double>(
                0.04 * batch + 0.01 * x, -0.015 * batch + 0.005 * x);
    TFHEpp::GLBasePlaintext<MultiLimbGLP, 250, 220> high_scale_plaintext;
    TFHEpp::GLBaseEncode(high_scale_plaintext, high_scale_values);
    TFHEpp::GLBaseSlotTable<MultiLimbGLP> high_scale_roundtrip;
    TFHEpp::GLBaseDecode(high_scale_roundtrip, high_scale_plaintext);
    double high_scale_error = 0;
    for (std::uint32_t batch = 0; batch < MultiLimbGLP::phi; batch++)
        for (std::uint32_t x = 0; x < MultiLimbGLP::matrix_dimension; x++)
            high_scale_error = std::max(
                high_scale_error, std::abs(high_scale_roundtrip(batch, x) -
                                           high_scale_values(batch, x)));
    std::cout << "  multi-limb high-scale encode error=" << high_scale_error
              << std::endl;
    if (high_scale_error > 1e-15) return 1;

    // Half-bootstrap a deliberately coefficient-ordered message.  This tests
    // dense-to-sparse switching, both Gaussian channels, MaskedColumn, HMux,
    // the balanced product tree, and dense-key output without relying on StC.
    Schedule::CoefficientCiphertext coefficient_ciphertext;
    TFHEpp::GLPlaintext<GLP, Schedule::q0_log_q, Schedule::input_log_delta>
        coefficient_plaintext;
    std::array<TFHEpp::GLBaseSlotTable<GLP>, GLP::matrix_dimension>
        expected_slices;
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
        for (std::uint32_t w = 0; w < GLP::phi; w++) {
            for (std::uint32_t x = 0; x < GLP::matrix_dimension; x++) {
                const std::complex<double> value(
                    0.018 * (1 + w) + 0.011 * x - 0.006 * y,
                    -0.014 * (1 + w) + 0.007 * x + 0.005 * y);
                expected_slices[y](w, x) = sineRefresh(value);
                coefficient_plaintext
                    .poly[y][TFHEpp::gl_detail::baseIndex<GLP>(0, x, w)] =
                    TFHEpp::ckks_detail::signedToLevel<
                        BootstrapTestBaseParameter, Schedule::q0_log_q>(
                        static_cast<__int128_t>(std::llround(std::ldexp(
                            value.real(), Schedule::input_log_delta))));
                coefficient_plaintext
                    .poly[y][TFHEpp::gl_detail::baseIndex<GLP>(1, x, w)] =
                    TFHEpp::ckks_detail::signedToLevel<
                        BootstrapTestBaseParameter, Schedule::q0_log_q>(
                        static_cast<__int128_t>(std::llround(std::ldexp(
                            value.imag(), Schedule::input_log_delta))));
            }
        }
    }
    TFHEpp::GLEncrypt(coefficient_ciphertext, coefficient_plaintext, dense_key);
    Schedule::OutputCiphertext half_output;
    TFHEpp::GLSHIPHalfBootstrap<Schedule>(half_output, coefficient_ciphertext,
                                          *bootstrap_key);
    TFHEpp::GLPlaintext<GLP, Schedule::output_log_q, Schedule::output_log_delta>
        half_phase;
    TFHEpp::GLDecrypt(half_phase, half_output, dense_key);
    double half_error = 0;
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
        TFHEpp::GLBasePlaintext<GLP, Schedule::output_log_q,
                                Schedule::output_log_delta>
            slice_plaintext;
        slice_plaintext.poly = half_phase.poly[y];
        TFHEpp::GLBaseSlotTable<GLP> got;
        TFHEpp::GLBaseDecode(got, slice_plaintext);
        for (std::uint32_t w = 0; w < GLP::phi; w++)
            for (std::uint32_t x = 0; x < GLP::matrix_dimension; x++)
                half_error = std::max(
                    half_error, std::abs(got(w, x) - expected_slices[y](w, x)));
    }
    std::cout << "  half bootstrap max error=" << half_error << std::endl;
    if (half_error > 0.08) return 1;

    const auto input = makeInput();
    Schedule::InputCiphertext input_ciphertext;
    TFHEpp::GLEncrypt(input_ciphertext, input, dense_key);

    Schedule::CoefficientCiphertext stc_output;
    TFHEpp::GLSHIPSlotsToCoefficients<Schedule>(stc_output, input_ciphertext,
                                                bootstrap_key->stc_key);
    std::cout << "  grouped StC completed" << std::endl;

    Schedule::OutputCiphertext output;
    TFHEpp::GLSHIPBootstrap<Schedule>(output, input_ciphertext, *bootstrap_key);
    TFHEpp::GLMatrixBatch<GLP> decoded;
    TFHEpp::GLDecrypt(decoded, output, dense_key);
    const double full_error = matrixError(decoded, input);
    std::cout << "  full bootstrap max error=" << full_error << std::endl;
    if (full_error > 0.12) return 1;

    std::cout << "PASS" << std::endl;
    return 0;
}
