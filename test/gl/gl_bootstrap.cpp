#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>
#include <vector>

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
using MultiLimbSchedule =
    TFHEpp::GLSHIPBootstrapSchedule<MultiLimbGLP, 24, 4, 8, 8, 90, 18, 3, 1, 2,
                                    2, 8, 16>;
using N512Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;

// Compile the production key-generation and end-to-end evaluation templates.
// Running them still requires the optimized Rader/NTT backend discussed in
// docs/GL.md; the regression below exercises the same path at toy dimensions.
[[maybe_unused]] void instantiateN512Bootstrap(
    TFHEpp::GLSHIPBootstrapKey<N512Schedule> &bootstrap_key,
    const TFHEpp::Key<typename N512Schedule::Parameter::baseP> &dense_key,
    const TFHEpp::Key<typename N512Schedule::Parameter::baseP> &sparse_key,
    const std::array<TFHEpp::GLSHIPSupportInterval,
                     N512Schedule::sparse_hamming_weight> &intervals,
    N512Schedule::OutputCiphertext &output,
    const N512Schedule::InputCiphertext &input)
{
    TFHEpp::GLSHIPBootstrapKeyGen(bootstrap_key, dense_key, sparse_key,
                                  intervals);
    TFHEpp::GLSHIPBootstrap<N512Schedule>(
        output, input, bootstrap_key,
        TFHEpp::GLSHIPBootstrapExecutionOptions{16, 256, true});
}

static_assert(Schedule::input_log_q == 40);
static_assert(Schedule::input_log_delta == 20);
static_assert(Schedule::tree_depth == 2);
static_assert(Schedule::output_log_q == 54);
using FirstProductRelinKey =
    std::tuple_element_t<0,
                         TFHEpp::GLSHIPProductRelinKeyChain<Schedule>::Tuple>;
using SecondProductRelinKey =
    std::tuple_element_t<1,
                         TFHEpp::GLSHIPProductRelinKeyChain<Schedule>::Tuple>;
static_assert(FirstProductRelinKey::log_q == 90);
static_assert(SecondProductRelinKey::log_q == 72);
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

    // The term-major evaluator uses an online tree so transient NTT caches
    // can be released after each sparse term.  Its carry order must remain
    // coefficient-for-coefficient identical to the original balanced tree.
    std::vector<TFHEpp::GLBaseCiphertextData<GLP>> tree_factors(
        Schedule::factor_count);
    for (std::size_t factor = 0; factor < tree_factors.size(); factor++)
        for (std::size_t component = 0; component < 2; component++)
            for (std::size_t coefficient = 0;
                 coefficient < BootstrapTestBaseParameter::n; coefficient++)
                tree_factors[factor][component][coefficient] =
                    static_cast<BootstrapTestBaseParameter::T>(
                        17 * factor + 7 * component + coefficient + 1);
    TFHEpp::GLBaseCiphertextData<GLP> reference_tree_product{};
    TFHEpp::gl_ship_detail::productTreeLevel<0, Schedule>(
        reference_tree_product, tree_factors,
        bootstrap_key->product_relin_keys);
    std::array<std::unique_ptr<std::vector<TFHEpp::GLBaseCiphertextData<GLP>>>,
               Schedule::tree_depth + 1>
        online_tree{};
    for (auto &factor : tree_factors) {
        std::vector<TFHEpp::GLBaseCiphertextData<GLP>> batch;
        batch.push_back(std::move(factor));
        TFHEpp::gl_ship_detail::insertProductTreeBatch<0, Schedule>(
            online_tree, std::move(batch), bootstrap_key->product_relin_keys);
    }
    if (!online_tree[Schedule::tree_depth] ||
        online_tree[Schedule::tree_depth]->size() != 1 ||
        (*online_tree[Schedule::tree_depth])[0] != reference_tree_product)
        return 1;

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

    TFHEpp::GLBasePolynomial<MultiLimbGLP> unpacked_digits{};
    TFHEpp::GLBasePolynomial<MultiLimbGLP> packed_digit_source{};
    packed_digit_source[0] =
        BootstrapMultiLimbTestBaseParameter::T::from_signed_i64(-32768);
    packed_digit_source[1] =
        BootstrapMultiLimbTestBaseParameter::T::from_signed_i64(32767);
    TFHEpp::GLPackedBasePolynomial<MultiLimbGLP, 16> packed_digits{};
    TFHEpp::gl_detail::packDigitPolynomial<MultiLimbGLP, 16>(
        packed_digits, packed_digit_source);
    TFHEpp::gl_detail::unpackDigitPolynomial<MultiLimbGLP, 16>(unpacked_digits,
                                                               packed_digits);
    if (unpacked_digits != packed_digit_source) return 1;

    // Exercise the actual n512 four-row, 85-bit primary decomposition without
    // invoking the intentionally slow production-size multiplication kernel.
    using N512GLP = TFHEpp::GL512p17Parameter;
    auto n512_source = std::make_unique<TFHEpp::GLBasePolynomial<N512GLP>>();
    auto n512_expected = std::make_unique<TFHEpp::GLBasePolynomial<N512GLP>>();
    auto n512_recombined =
        std::make_unique<TFHEpp::GLBasePolynomial<N512GLP>>();
    (*n512_source)[0] = typename N512GLP::T{1} << 337;
    (*n512_source)[1] = typename N512GLP::T{0} - typename N512GLP::T{1};
    (*n512_source)[2] = (typename N512GLP::T{1} << 337) +
                        (typename N512GLP::T{1} << 169) +
                        typename N512GLP::T{12345};
    *n512_expected = *n512_source;
    TFHEpp::gl_detail::reduce<N512GLP, 338>(*n512_expected);
    const auto n512_primary_rows =
        TFHEpp::gl_detail::activeDecompose<N512GLP, 338, 85>(*n512_source);
    if (n512_primary_rows.size() != 4) return 1;
    TFHEpp::gl_detail::activeRecombine<N512GLP, 338, 85>(*n512_recombined,
                                                         n512_primary_rows);
    if (*n512_recombined != *n512_expected) return 1;

    // Exercise the production-width relinearize-before-rescale kernel.  The
    // product itself is reduced modulo its input Q, so only the low Q bits of
    // the multi-limb multiplication are needed before DD key switching.
    TFHEpp::Key<BootstrapMultiLimbTestBaseParameter> multi_limb_key{};
    multi_limb_key[0] = 1;
    multi_limb_key[1] = static_cast<BootstrapMultiLimbTestBaseParameter::T>(-1);
    multi_limb_key[4] = 1;

    // Exercise the complete unseeded packed-key archive with a multi-limb
    // torus, matching the coefficient type used by the n512 profile.
    TFHEpp::Key<BootstrapMultiLimbTestBaseParameter> multi_limb_sparse_key{};
    multi_limb_sparse_key[0] = 1;
    multi_limb_sparse_key[3] =
        static_cast<BootstrapMultiLimbTestBaseParameter::T>(-1);
    multi_limb_sparse_key[10] = 1;
    auto archived_key =
        std::make_unique<TFHEpp::GLSHIPBootstrapKey<MultiLimbSchedule>>();
    TFHEpp::GLSHIPBootstrapKeyGen(*archived_key, multi_limb_key,
                                  multi_limb_sparse_key, intervals);
    const std::filesystem::path archive_path =
        std::filesystem::temp_directory_path() /
        "tfhepp_gl_ship_unseeded_packed_test.bin";
    std::error_code archive_error;
    std::filesystem::remove(archive_path, archive_error);
    TFHEpp::GLSHIPSaveBootstrapKey<MultiLimbSchedule>(archive_path,
                                                      *archived_key);
    if (!std::filesystem::is_regular_file(archive_path) ||
        std::filesystem::file_size(archive_path) == 0)
        return 1;
    auto restored_key =
        std::make_unique<TFHEpp::GLSHIPBootstrapKey<MultiLimbSchedule>>();
    TFHEpp::GLSHIPLoadBootstrapKey(*restored_key, archive_path);
    std::filesystem::remove(archive_path, archive_error);
    if (TFHEpp::GLSHIPBootstrapKeyPackedPayloadBytes(*restored_key) !=
            TFHEpp::GLSHIPBootstrapKeyPackedPayloadBytes(*archived_key) ||
        restored_key->stc_key.w_rotation_keys.size() !=
            archived_key->stc_key.w_rotation_keys.size() ||
        restored_key->hmux_keys[0].stages.size() !=
            archived_key->hmux_keys[0].stages.size() ||
        restored_key->dense_to_sparse_key.data.size() !=
            archived_key->dense_to_sparse_key.data.size() ||
        restored_key->dense_to_sparse_key.data.front()[0][0] !=
            archived_key->dense_to_sparse_key.data.front()[0][0] ||
        restored_key->masked_column_keys[0].encrypted_masks.front()[0][0] !=
            archived_key->masked_column_keys[0].encrypted_masks.front()[0][0])
        return 1;
    archived_key.reset();
    restored_key.reset();

    TFHEpp::GLBaseSlotTable<MultiLimbGLP> product_lhs_values;
    TFHEpp::GLBaseSlotTable<MultiLimbGLP> product_rhs_values;
    for (std::uint32_t batch = 0; batch < MultiLimbGLP::phi; batch++) {
        for (std::uint32_t x = 0; x < MultiLimbGLP::matrix_dimension; x++) {
            product_lhs_values(batch, x) = {0.025 * (1 + batch) + 0.01 * x,
                                            -0.012 * (1 + x)};
            product_rhs_values(batch, x) = {-0.018 * (1 + x),
                                            0.009 * (1 + batch)};
        }
    }
    TFHEpp::GLBasePlaintext<MultiLimbGLP, 180, 30> product_lhs_plain;
    TFHEpp::GLBasePlaintext<MultiLimbGLP, 180, 30> product_rhs_plain;
    TFHEpp::GLBaseEncode(product_lhs_plain, product_lhs_values);
    TFHEpp::GLBaseEncode(product_rhs_plain, product_rhs_values);
    TFHEpp::GLBaseCiphertext<MultiLimbGLP, 180, 30> product_lhs;
    TFHEpp::GLBaseCiphertext<MultiLimbGLP, 180, 30> product_rhs;
    TFHEpp::GLBaseEncrypt(product_lhs, product_lhs_plain, multi_limb_key);
    TFHEpp::GLBaseEncrypt(product_rhs, product_rhs_plain, multi_limb_key);
    TFHEpp::GLHadamardRelinKey<MultiLimbGLP, 180, 8, 8> product_relin_key;
    TFHEpp::GLHadamardRelinKeyGen(product_relin_key, multi_limb_key);
    TFHEpp::GLBaseCiphertext<MultiLimbGLP, 150, 30> product_ciphertext;
    TFHEpp::GLBaseHadamardMultiply<MultiLimbGLP, 180, 30, 180, 30, 30, 8, 8>(
        product_ciphertext, product_lhs, product_rhs, product_relin_key);
    TFHEpp::GLBasePlaintext<MultiLimbGLP, 150, 30> product_plaintext;
    TFHEpp::GLBaseDecrypt(product_plaintext, product_ciphertext,
                          multi_limb_key);
    TFHEpp::GLBaseSlotTable<MultiLimbGLP> product_values;
    TFHEpp::GLBaseDecode(product_values, product_plaintext);
    double multi_limb_product_error = 0;
    for (std::uint32_t batch = 0; batch < MultiLimbGLP::phi; batch++)
        for (std::uint32_t x = 0; x < MultiLimbGLP::matrix_dimension; x++)
            multi_limb_product_error =
                std::max(multi_limb_product_error,
                         std::abs(product_values(batch, x) -
                                  product_lhs_values(batch, x) *
                                      product_rhs_values(batch, x)));
    std::cout << "  multi-limb fused product error=" << multi_limb_product_error
              << std::endl;
    if (multi_limb_product_error > 2e-5) return 1;

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
    Schedule::OutputCiphertext staged_half_output;
    TFHEpp::GLSHIPHalfBootstrap<Schedule>(
        staged_half_output, coefficient_ciphertext, *bootstrap_key,
        TFHEpp::GLSHIPBootstrapExecutionOptions{1, 3, true});
    for (std::size_t component = 0; component < 2; component++)
        for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++)
            if (staged_half_output[component][y] !=
                half_output[component][y])
                return 1;
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

    // The production StC path uses the compact W^0/signed-small trace.  At
    // toy dimensions it must agree coefficient-for-coefficient with the
    // generic GL trace followed by the matrix-dimension normalization.
    using AfterX =
        TFHEpp::GLRawProductCiphertext<GLP, Schedule::input_log_q,
                                       Schedule::input_log_delta,
                                       Schedule::x_transform_log_scale>;
    TFHEpp::GLPlaintext<GLP, Schedule::input_log_q,
                        Schedule::x_transform_log_scale>
        x_plaintext;
    TFHEpp::gl_ship_detail::buildXTransformPlaintext<Schedule>(x_plaintext);
    AfterX reference_x_product;
    AfterX compact_x_product;
    TFHEpp::GLPlaintextMatrixMultiplyRaw(reference_x_product, input_ciphertext,
                                         x_plaintext);
    TFHEpp::GLXTransformMatrixMultiplyRaw<Schedule>(compact_x_product,
                                                    input_ciphertext);
    for (std::size_t component = 0; component < 2; component++)
        for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++)
            for (std::size_t coefficient = 0; coefficient < GLP::baseP::n;
                 coefficient++)
                if (reference_x_product[component][y][coefficient] !=
                    compact_x_product[component][y][coefficient])
                    return 1;

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
