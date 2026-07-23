#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <tfhe++.hpp>

namespace {

struct GLTestBaseParameter {
    static constexpr std::uint32_t n = 16;
    static constexpr std::uint32_t k = 1;
    using T = __uint128_t;
    static constexpr std::uint32_t B̅gbit = 8;
};

using GLP = TFHEpp::GLParameter<GLTestBaseParameter, 2, 5>;
static_assert(TFHEpp::GL256p17Parameter::slice_dimension == 8192);
constexpr std::uint32_t input_log_q = 56;
constexpr std::uint32_t log_delta = 12;
constexpr std::uint32_t output_log_q = input_log_q - log_delta;
static_assert(TFHEpp::GLDDBigKeySwitchKey<GLP, output_log_q, 8, 8>::key_log_q ==
              output_log_q + GLP::auxiliary_log_q);

TFHEpp::GLMatrixBatch<GLP> makeBatch(const double offset)
{
    TFHEpp::GLMatrixBatch<GLP> batch;
    for (std::uint32_t b = 0; b < GLP::phi; b++) {
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++) {
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++) {
                const double real =
                    offset + 0.35 * b + 0.2 * row - 0.15 * column;
                const double imag =
                    -0.1 * offset + 0.08 * b - 0.12 * row + 0.09 * column;
                batch(b, row, column) = {real, imag};
            }
        }
    }
    return batch;
}

TFHEpp::GLMatrixBatch<GLP> adjointProduct(const TFHEpp::GLMatrixBatch<GLP> &lhs,
                                          const TFHEpp::GLMatrixBatch<GLP> &rhs)
{
    TFHEpp::GLMatrixBatch<GLP> result;
    for (std::uint32_t b = 0; b < GLP::phi; b++) {
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++) {
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++) {
                std::complex<double> value = 0;
                for (std::uint32_t inner = 0; inner < GLP::matrix_dimension;
                     inner++)
                    value +=
                        lhs(b, row, inner) * std::conj(rhs(b, column, inner));
                result(b, row, column) = value;
            }
        }
    }
    return result;
}

TFHEpp::GLMatrixBatch<GLP> hadamardProduct(
    const TFHEpp::GLMatrixBatch<GLP> &lhs,
    const TFHEpp::GLMatrixBatch<GLP> &rhs)
{
    TFHEpp::GLMatrixBatch<GLP> result;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                result(b, row, column) =
                    lhs(b, row, column) * rhs(b, row, column);
    return result;
}

TFHEpp::GLMatrixBatch<GLP> conjugateBatch(
    const TFHEpp::GLMatrixBatch<GLP> &input)
{
    TFHEpp::GLMatrixBatch<GLP> result;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                result(b, row, column) = std::conj(input(b, row, column));
    return result;
}

TFHEpp::GLMatrixBatch<GLP> transposeBatch(
    const TFHEpp::GLMatrixBatch<GLP> &input, const bool conjugate)
{
    TFHEpp::GLMatrixBatch<GLP> result;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++) {
                const auto value = input(b, column, row);
                result(b, row, column) = conjugate ? std::conj(value) : value;
            }
    return result;
}

TFHEpp::GLMatrixBatch<GLP> rotateRowsBatch(
    const TFHEpp::GLMatrixBatch<GLP> &input, const std::uint32_t amount)
{
    TFHEpp::GLMatrixBatch<GLP> result;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                result(b, row, column) =
                    input(b, (row + amount) % GLP::matrix_dimension, column);
    return result;
}

TFHEpp::GLMatrixBatch<GLP> rotateColumnsBatch(
    const TFHEpp::GLMatrixBatch<GLP> &input, const std::uint32_t amount)
{
    TFHEpp::GLMatrixBatch<GLP> result;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                result(b, row, column) =
                    input(b, row, (column + amount) % GLP::matrix_dimension);
    return result;
}

TFHEpp::GLMatrixBatch<GLP> rotateBatches(
    const TFHEpp::GLMatrixBatch<GLP> &input, const std::uint32_t amount)
{
    TFHEpp::GLMatrixBatch<GLP> result;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                result(b, row, column) =
                    input((b + amount) % GLP::phi, row, column);
    return result;
}

double maxError(const TFHEpp::GLMatrixBatch<GLP> &got,
                const TFHEpp::GLMatrixBatch<GLP> &want)
{
    double error = 0;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                error = std::max(error, std::abs(got(b, row, column) -
                                                 want(b, row, column)));
    return error;
}

bool requireError(const char *label, const TFHEpp::GLMatrixBatch<GLP> &got,
                  const TFHEpp::GLMatrixBatch<GLP> &want,
                  const double tolerance)
{
    const double error = maxError(got, want);
    std::cout << "  " << label << " max error=" << error << std::endl;
    if (error <= tolerance) return true;
    std::cerr << "FAIL: " << label << " exceeded tolerance " << tolerance
              << std::endl;
    return false;
}

}  // namespace

int main()
{
    std::cout << "GL scheme / Double Decomposition regression" << std::endl;

    const auto lhs = makeBatch(0.4);
    const auto rhs = makeBatch(-0.25);

    TFHEpp::GLPlaintext<GLP, input_log_q, log_delta> lhs_plain;
    TFHEpp::GLPlaintext<GLP, input_log_q, log_delta> rhs_plain;
    TFHEpp::GLEncode(lhs_plain, lhs);
    TFHEpp::GLEncode(rhs_plain, rhs);

    TFHEpp::GLMatrixBatch<GLP> encoded_roundtrip;
    TFHEpp::GLDecode(encoded_roundtrip, lhs_plain);
    if (!requireError("encode/decode", encoded_roundtrip, lhs, 0.002)) return 1;

    TFHEpp::GLPolynomial<GLP> trace_poly;
    TFHEpp::GLTraceProduct(trace_poly, lhs_plain.poly, rhs_plain.poly);
    TFHEpp::GLPlaintext<GLP, input_log_q, 2 * log_delta> trace_plain;
    trace_plain.poly = std::move(trace_poly);
    TFHEpp::GLMatrixBatch<GLP> trace_decoded;
    TFHEpp::GLDecode(trace_decoded, trace_plain);
    auto expected = adjointProduct(lhs, rhs);
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                expected(b, row, column) /= GLP::matrix_dimension;
    if (!requireError("paper trace product", trace_decoded, expected, 0.004))
        return 1;

    TFHEpp::Key<GLTestBaseParameter> key{};
    key[0] = 1;
    key[1] = static_cast<GLTestBaseParameter::T>(-1);
    key[4] = 1;

    TFHEpp::GLCiphertext<GLP, input_log_q, log_delta> lhs_ciphertext;
    TFHEpp::GLCiphertext<GLP, input_log_q, log_delta> rhs_ciphertext;
    TFHEpp::GLEncrypt(lhs_ciphertext, lhs_plain, key);
    TFHEpp::GLEncrypt(rhs_ciphertext, rhs_plain, key);

    TFHEpp::GLMatrixBatch<GLP> encrypted_roundtrip;
    TFHEpp::GLDecrypt(encrypted_roundtrip, lhs_ciphertext, key);
    if (!requireError("encrypt/decrypt", encrypted_roundtrip, lhs, 0.02))
        return 1;

    TFHEpp::GLMatrixRelinKey<GLP, output_log_q, 8, 8> relin_key;
    TFHEpp::GLMatrixRelinKeyGen(relin_key, key);

    TFHEpp::GLMatrixMultResult<GLP, input_log_q, log_delta, input_log_q,
                               log_delta>
        product;
    TFHEpp::GLMatrixMultiply(product, lhs_ciphertext, rhs_ciphertext,
                             relin_key);

    TFHEpp::GLMatrixBatch<GLP> product_decrypted;
    TFHEpp::GLDecrypt(product_decrypted, product, key);
    expected = adjointProduct(lhs, rhs);
    if (!requireError("encrypted matrix product", product_decrypted, expected,
                      0.04))
        return 1;

    decltype(product) unnormalized_product;
    TFHEpp::GLMatrixMultiply<false>(unnormalized_product, lhs_ciphertext,
                                    rhs_ciphertext, relin_key);
    TFHEpp::GLMatrixBatch<GLP> unnormalized_decrypted;
    TFHEpp::GLDecrypt(unnormalized_decrypted, unnormalized_product, key);
    auto expected_unnormalized = expected;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                expected_unnormalized(b, row, column) /= GLP::matrix_dimension;
    if (!requireError("unnormalized encrypted trace", unnormalized_decrypted,
                      expected_unnormalized, 0.03))
        return 1;

    TFHEpp::GLMatrixMultResult<GLP, input_log_q, log_delta, input_log_q,
                               log_delta>
        plaintext_product;
    TFHEpp::GLPlaintextMatrixMultiply(plaintext_product, lhs_ciphertext,
                                      rhs_plain);
    TFHEpp::GLMatrixBatch<GLP> plaintext_product_decrypted;
    TFHEpp::GLDecrypt(plaintext_product_decrypted, plaintext_product, key);
    if (!requireError("plaintext/ciphertext matrix product",
                      plaintext_product_decrypted, expected, 0.04))
        return 1;

    TFHEpp::GLPlaintext<GLP, input_log_q - 4, log_delta - 2> rhs_low_plain;
    TFHEpp::GLEncode(rhs_low_plain, rhs);
    TFHEpp::GLMatrixMultResult<GLP, input_log_q, log_delta, input_log_q - 4,
                               log_delta - 2>
        mixed_level_product;
    TFHEpp::GLPlaintextMatrixMultiply(mixed_level_product, lhs_ciphertext,
                                      rhs_low_plain);
    TFHEpp::GLMatrixBatch<GLP> mixed_level_decrypted;
    TFHEpp::GLDecrypt(mixed_level_decrypted, mixed_level_product, key);
    if (!requireError("mixed-level plaintext/ciphertext product",
                      mixed_level_decrypted, expected, 0.05))
        return 1;

    TFHEpp::GLHadamardRelinKey<GLP, output_log_q, 8, 8> hadamard_relin_key;
    TFHEpp::GLHadamardRelinKeyGen(hadamard_relin_key, key);

    TFHEpp::GLMatrixMultResult<GLP, input_log_q, log_delta, input_log_q,
                               log_delta>
        hadamard;
    TFHEpp::GLHadamardMultiply(hadamard, lhs_ciphertext, rhs_ciphertext,
                               hadamard_relin_key);

    TFHEpp::GLMatrixBatch<GLP> hadamard_decrypted;
    TFHEpp::GLDecrypt(hadamard_decrypted, hadamard, key);
    const auto expected_hadamard = hadamardProduct(lhs, rhs);
    if (!requireError("encrypted Hadamard product", hadamard_decrypted,
                      expected_hadamard, 0.04))
        return 1;

    TFHEpp::GLConjugationKey<GLP, input_log_q, 8, 8> conjugation_key;
    TFHEpp::GLConjugationKeyGen(conjugation_key, key);
    TFHEpp::GLCiphertext<GLP, input_log_q, log_delta> transformed;
    TFHEpp::GLConjugate(transformed, lhs_ciphertext, conjugation_key);
    TFHEpp::GLMatrixBatch<GLP> transformed_decrypted;
    TFHEpp::GLDecrypt(transformed_decrypted, transformed, key);
    if (!requireError("conjugation", transformed_decrypted, conjugateBatch(lhs),
                      0.02))
        return 1;

    TFHEpp::GLRowRotationKey<GLP, input_log_q, 8, 8> row_rotation_key;
    TFHEpp::GLRowRotationKeyGen(row_rotation_key, key, 1);
    TFHEpp::GLRotateRows(transformed, lhs_ciphertext, row_rotation_key);
    TFHEpp::GLDecrypt(transformed_decrypted, transformed, key);
    if (!requireError("row rotation", transformed_decrypted,
                      rotateRowsBatch(lhs, 1), 0.02))
        return 1;

    TFHEpp::GLRotateColumns(transformed, lhs_ciphertext, 1);
    TFHEpp::GLDecrypt(transformed_decrypted, transformed, key);
    if (!requireError("column rotation", transformed_decrypted,
                      rotateColumnsBatch(lhs, 1), 0.02))
        return 1;

    TFHEpp::GLBatchRotationKey<GLP, input_log_q, 8, 8> batch_rotation_key;
    TFHEpp::GLBatchRotationKeyGen(batch_rotation_key, key, 1);
    TFHEpp::GLRotateBatches(transformed, lhs_ciphertext, batch_rotation_key);
    TFHEpp::GLDecrypt(transformed_decrypted, transformed, key);
    if (!requireError("batch rotation", transformed_decrypted,
                      rotateBatches(lhs, 1), 0.02))
        return 1;

    TFHEpp::GLTransposeKey<GLP, input_log_q, 8, 8> transpose_key;
    TFHEpp::GLTransposeKeyGen(transpose_key, key);
    TFHEpp::GLTranspose(transformed, lhs_ciphertext, transpose_key);
    TFHEpp::GLDecrypt(transformed_decrypted, transformed, key);
    if (!requireError("transpose", transformed_decrypted,
                      transposeBatch(lhs, false), 0.02))
        return 1;

    TFHEpp::GLConjugateTransposeKey<GLP, input_log_q, 8, 8>
        conjugate_transpose_key;
    TFHEpp::GLConjugateTransposeKeyGen(conjugate_transpose_key, key);
    TFHEpp::GLConjugateTranspose(transformed, lhs_ciphertext,
                                 conjugate_transpose_key);
    TFHEpp::GLDecrypt(transformed_decrypted, transformed, key);
    if (!requireError("conjugate transpose", transformed_decrypted,
                      transposeBatch(lhs, true), 0.02))
        return 1;

    TFHEpp::GLCiphertext<GLP, input_log_q, log_delta> sum;
    TFHEpp::GLAdd(sum, lhs_ciphertext, rhs_ciphertext);
    TFHEpp::GLMatrixBatch<GLP> sum_decrypted;
    TFHEpp::GLDecrypt(sum_decrypted, sum, key);
    auto expected_sum = lhs;
    for (std::uint32_t b = 0; b < GLP::phi; b++)
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            for (std::uint32_t column = 0; column < GLP::matrix_dimension;
                 column++)
                expected_sum(b, row, column) += rhs(b, row, column);
    if (!requireError("encrypted addition", sum_decrypted, expected_sum, 0.03))
        return 1;

    std::cout << "PASS" << std::endl;
    return 0;
}
