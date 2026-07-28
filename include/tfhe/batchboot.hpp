#pragma once

#include <array>
#include <bit>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>

#include "evalkeygens.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace TFHEpp {

// BatchBoot, USENIX Security 2026, Algorithms 2--4.
//
// The implementation deliberately uses TFHEpp's native TRLWE/TRGSW types and
// FFT backend.  BatchBoot's two-bit monomial multiplier is represented by a
// freshly encrypted one-hot selector for every radix-4 digit.  In particular,
// selector ciphertexts are not shared between digits: sharing them would make
// equality of secret digits visible from the evaluation-key representation.

template <class SourceP, std::uint32_t Slots>
struct BatchRingSwitchParameter {
    static_assert(std::has_single_bit(Slots));
    static_assert(Slots <= SourceP::n && SourceP::n % Slots == 0);

    using T = typename SourceP::T;
    static constexpr std::uint32_t n = Slots;
    static constexpr std::uint32_t nbit = std::countr_zero(Slots);
    static constexpr std::uint32_t k = SourceP::n / Slots;
};

template <class SourceP, std::uint32_t Slots>
using BatchRingSwitchP = BatchRingSwitchParameter<SourceP, Slots>;

// Decompose s(X) = sum_r X^r s_r(X^k) into the module secret
// (s_0, ..., s_{k-1}).
template <class SourceP, std::uint32_t Slots>
void BatchRingSwitchSecret(Key<BatchRingSwitchP<SourceP, Slots>> &result,
                           const Key<SourceP> &source_key)
{
    static_assert(SourceP::k == 1,
                  "RLWE-to-MLWE extraction expects an RLWE source key");
    using ModuleP = BatchRingSwitchP<SourceP, Slots>;
    for (std::uint32_t r = 0; r < ModuleP::k; r++)
        for (std::uint32_t i = 0; i < Slots; i++)
            result[r * Slots + i] = source_key[i * ModuleP::k + r];
}

template <class SourceP, std::uint32_t Slots>
Key<BatchRingSwitchP<SourceP, Slots>> BatchRingSwitchSecret(
    const Key<SourceP> &source_key)
{
    Key<BatchRingSwitchP<SourceP, Slots>> result;
    BatchRingSwitchSecret<SourceP, Slots>(result, source_key);
    return result;
}

// Extract the coefficient-0 coset of an RLWE ciphertext as an MLWE
// ciphertext.  If a(X) = sum_r X^r a_r(Y), Y=X^k, then
//
//   Coef_0(a*s) = a_0*s_0 + Y * sum_{r=1}^{k-1} a_{k-r}*s_r.
//
// The module mask below applies exactly that permutation and multiplication
// by Y.  The body is the coefficient-0 coset of b.
template <class SourceP, std::uint32_t Slots>
void BatchRingSwitch(
    TRLWE<BatchRingSwitchP<SourceP, Slots>> &result,
    const TRLWE<SourceP> &source)
{
    static_assert(SourceP::k == 1,
                  "RLWE-to-MLWE extraction expects an RLWE ciphertext");
    using ModuleP = BatchRingSwitchP<SourceP, Slots>;
    constexpr std::uint32_t module_rank = ModuleP::k;

    for (std::uint32_t i = 0; i < Slots; i++) {
        result[0][i] = source[0][i * module_rank];
        result[module_rank][i] = source[1][i * module_rank];
    }

    for (std::uint32_t r = 1; r < module_rank; r++) {
        const std::uint32_t source_coset = module_rank - r;
        result[r][0] = -source[0][(Slots - 1) * module_rank + source_coset];
        for (std::uint32_t i = 1; i < Slots; i++)
            result[r][i] = source[0][(i - 1) * module_rank + source_coset];
    }
}

template <class P>
struct BatchEMPDigitKey {
    // Rows 0,1,2 select shifts 1,2,3; row 3 selects shift 0.
    std::array<TRGSWFFT<P>, 4> direct;
    // Parametrized external-product keys for coefficients that cross X^n.
    std::array<TRGSWFFT<P>, 3> inverse;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(direct, inverse);
    }
};

template <class P>
struct BatchEMPKey {
    BatchEMPKey() = default;
    BatchEMPKey(const BatchEMPKey &) = delete;
    BatchEMPKey &operator=(const BatchEMPKey &) = delete;
    BatchEMPKey(BatchEMPKey &&) noexcept = default;
    BatchEMPKey &operator=(BatchEMPKey &&) noexcept = default;

    // Digits are consumed sequentially by BatchEMP.  Store them inline so the
    // hot evaluation-key rows do not require one allocation and pointer chase
    // per radix-4 digit.
    std::vector<BatchEMPDigitKey<P>> digits;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(digits);
    }
};

template <class P>
struct BatchBootComponentKey {
    BatchBootComponentKey() = default;
    BatchBootComponentKey(const BatchBootComponentKey &) = delete;
    BatchBootComponentKey &operator=(const BatchBootComponentKey &) = delete;
    BatchBootComponentKey(BatchBootComponentKey &&) noexcept = default;
    BatchBootComponentKey &operator=(BatchBootComponentKey &&) noexcept =
        default;

    // Sparse positions are also visited sequentially during accumulation.
    std::vector<BatchEMPKey<P>> negative_gaps;
    std::unique_ptr<BatchEMPKey<P>> final_positive;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(negative_gaps, final_positive);
    }
};

template <class DomainP, class TargetP>
struct BatchBootKey {
    BatchBootKey() = default;
    BatchBootKey(const BatchBootKey &) = delete;
    BatchBootKey &operator=(const BatchBootKey &) = delete;
    BatchBootKey(BatchBootKey &&) noexcept = default;
    BatchBootKey &operator=(BatchBootKey &&) noexcept = default;

    std::array<BatchBootComponentKey<TargetP>, DomainP::k> components;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(components);
    }
};

template <class P>
void BatchTauInverseTRGSWEncrypt(TRGSWFFT<P> &result, const bool message,
                                 const Key<P> &key)
{
    static_assert(P::l̅ == 1 && P::l̅ₐ == 1,
                  "BatchBoot EMP currently requires standard decomposition");

    constexpr auto nonce_gadget = hgen<P, true>();
    constexpr auto body_gadget = hgen<P, false>();
    auto coefficient_key = std::make_unique<TRGSW<P>>();
    for (auto &row : *coefficient_key) trlweSymEncryptZero<P>(row, key);

    if (message) {
        for (std::uint32_t component = 0; component < P::k; component++) {
            Polynomial<P> part_key;
            for (std::uint32_t i = 0; i < P::n; i++)
                part_key[i] = key[component * P::n + i];
            Polynomial<P> inverse_key;
            Automorphism<P>(inverse_key, part_key, 2 * P::n - 1);

            for (std::uint32_t level = 0; level < P::lₐ; level++) {
                auto &body =
                    (*coefficient_key)[component * P::lₐ + level][P::k];
                for (std::uint32_t i = 0; i < P::n; i++)
                    body[i] -= inverse_key[i] * nonce_gadget[level];
            }
        }

        for (std::uint32_t level = 0; level < P::l; level++)
            (*coefficient_key)[P::k * P::lₐ + level][P::k][0] +=
                body_gadget[level];
    }

    ApplyFFT2trgsw<P>(result, *coefficient_key);
}

template <class P>
void BatchEMPKeyGen(BatchEMPKey<P> &result, const std::uint32_t exponent,
                    const std::uint32_t slots, const Key<P> &target_key)
{
    if (!std::has_single_bit(slots) || slots > P::n)
        throw std::invalid_argument(
            "BatchEMPKeyGen: slots must be a power of two no larger than N");
    if (exponent >= slots)
        throw std::invalid_argument(
            "BatchEMPKeyGen: exponent must be smaller than slots");

    const std::uint32_t digit_count =
        (std::bit_width(slots - 1U) + 1U) / 2U;
    result.digits.clear();
    result.digits.reserve(digit_count);

    for (std::uint32_t level = 0; level < digit_count; level++) {
        result.digits.emplace_back();
        auto &digit_key = result.digits.back();
        const std::uint32_t digit = (exponent >> (2 * level)) & 3U;
        const std::uint32_t selected_row = digit == 0 ? 3 : digit - 1;

        for (std::uint32_t row = 0; row < 4; row++) {
            Polynomial<P> selector{};
            selector[0] = row == selected_row ? 1 : 0;
            trgswSymEncrypt<P>(digit_key.direct[row], selector, target_key);
        }
        for (std::uint32_t row = 0; row < 3; row++)
            BatchTauInverseTRGSWEncrypt<P>(digit_key.inverse[row],
                                           row == selected_row, target_key);
    }
}

template <class DomainP, class TargetP>
void BatchBootKeyGen(BatchBootKey<DomainP, TargetP> &result,
                     const Key<DomainP> &domain_key,
                     const Key<TargetP> &target_key)
{
    static_assert(std::is_same_v<typename DomainP::T, typename TargetP::T>,
                  "BatchBoot currently expects equal domain/target torus types");
    static_assert(TargetP::l̅ == 1 && TargetP::l̅ₐ == 1,
                  "BatchBoot EMP currently requires standard decomposition");
    static_assert(DomainP::n <= TargetP::n);
    static_assert(std::has_single_bit(DomainP::n));

    for (std::uint32_t component = 0; component < DomainP::k; component++) {
        auto &component_key = result.components[component];
        component_key.negative_gaps.clear();
        component_key.final_positive.reset();

        std::vector<std::uint32_t> positions;
        for (std::uint32_t i = 0; i < DomainP::n; i++) {
            const auto value = domain_key[component * DomainP::n + i];
            if (value == 1)
                positions.push_back(i);
            else if (value != 0)
                throw std::invalid_argument(
                    "BatchBootKeyGen: the domain secret must be sparse binary");
        }
        if (positions.empty()) continue;

        std::uint32_t previous = 0;
        component_key.negative_gaps.reserve(positions.size());
        for (const std::uint32_t position : positions) {
            component_key.negative_gaps.emplace_back();
            BatchEMPKeyGen<TargetP>(component_key.negative_gaps.back(),
                                    position - previous,
                                    DomainP::n, target_key);
            previous = position;
        }

        component_key.final_positive =
            std::make_unique<BatchEMPKey<TargetP>>();
        BatchEMPKeyGen<TargetP>(*component_key.final_positive, previous,
                                DomainP::n, target_key);
    }
}

template <class P>
struct BatchHoistedTRLWEFFT {
    static constexpr std::uint32_t rows = P::k * P::lₐ + P::l;
    aligned_array<PolynomialInFD<P>, rows> decomposition;
};

template <class P>
struct BatchEMPWorkspace {
    std::vector<BatchHoistedTRLWEFFT<P>> hoisted;
    std::vector<TRLWEInFD<P>> transformed;

    void prepare(const std::uint32_t slots)
    {
        if (hoisted.size() != slots) hoisted.resize(slots);
        if (transformed.size() != slots) transformed.resize(slots);
    }
};

template <class P>
void BatchConjugateInFD(PolynomialInFD<P> &result,
                        const PolynomialInFD<P> &input)
{
#ifdef USE_INTERLEAVED_FORMAT
    for (std::uint32_t i = 0; i < P::n; i += 2) {
        result[i] = input[i];
        result[i + 1] = -input[i + 1];
    }
#else
    for (std::uint32_t i = 0; i < P::n / 2; i++) {
        result[i] = input[i];
        result[i + P::n / 2] = -input[i + P::n / 2];
    }
#endif
}

template <class P>
void BatchHoist(BatchHoistedTRLWEFFT<P> &result, const TRLWE<P> &input)
{
    static_assert(P::l̅ == 1 && P::l̅ₐ == 1,
                  "BatchBoot EMP currently requires standard decomposition");

    for (std::uint32_t component = 0; component < P::k; component++) {
        DecomposedNoncePolynomial<P> decomposed;
        NonceDecomposition<P>(decomposed, input[component]);
        for (std::uint32_t level = 0; level < P::lₐ; level++)
            TwistIFFT<P>(result.decomposition[component * P::lₐ + level],
                         decomposed[level]);
    }

    DecomposedPolynomial<P> decomposed;
    Decomposition<P>(decomposed, input[P::k]);
    for (std::uint32_t level = 0; level < P::l; level++)
        TwistIFFT<P>(result.decomposition[P::k * P::lₐ + level],
                     decomposed[level]);
}

template <class P>
void BatchExternalProductAddFD(TRLWEInFD<P> &result,
                               const BatchHoistedTRLWEFFT<P> &input,
                               const TRGSWFFT<P> &key,
                               const bool inverse,
                               const bool overwrite = false)
{
    constexpr std::uint32_t rows = P::k * P::lₐ + P::l;
#if defined(__AVX512F__) && defined(USE_INTERLEAVED_FORMAT)
    if constexpr (rows == 2) {
        for (std::uint32_t i = 0; i < P::n; i += 8) {
            const __m512d value0 =
                _mm512_load_pd(input.decomposition[0].data() + i);
            const __m512d value1 =
                _mm512_load_pd(input.decomposition[1].data() + i);
            const __m512d value0_re =
                _mm512_unpacklo_pd(value0, value0);
            const __m512d value1_re =
                _mm512_unpacklo_pd(value1, value1);
            __m512d value0_im = _mm512_unpackhi_pd(value0, value0);
            __m512d value1_im = _mm512_unpackhi_pd(value1, value1);
            if (inverse) {
                const __m512d zero = _mm512_setzero_pd();
                value0_im = _mm512_sub_pd(zero, value0_im);
                value1_im = _mm512_sub_pd(zero, value1_im);
            }
            for (std::uint32_t component = 0; component <= P::k;
                 component++) {
                const __m512d key0 =
                    _mm512_load_pd(key[0][component].data() + i);
                const __m512d key1 =
                    _mm512_load_pd(key[1][component].data() + i);
                const __m512d product0 = _mm512_fmaddsub_pd(
                    value0_re, key0,
                    _mm512_mul_pd(value0_im,
                                  _mm512_permute_pd(key0, 0x55)));
                const __m512d product1 = _mm512_fmaddsub_pd(
                    value1_re, key1,
                    _mm512_mul_pd(value1_im,
                                  _mm512_permute_pd(key1, 0x55)));
                __m512d sum =
                    overwrite
                        ? product0
                        : _mm512_add_pd(
                              _mm512_load_pd(result[component].data() + i),
                              product0);
                sum = _mm512_add_pd(sum, product1);
                _mm512_store_pd(result[component].data() + i, sum);
            }
        }
        return;
    }
#endif

    alignas(64) PolynomialInFD<P> inverse_digit;
    for (std::uint32_t row = 0; row < rows; row++) {
        if (inverse) {
            BatchConjugateInFD<P>(inverse_digit, input.decomposition[row]);
            if (overwrite && row == 0)
                MulInFD_Multi<P::n, P::k + 1>(result, inverse_digit,
                                               key[row]);
            else
                FMAInFD_Multi<P::n, P::k + 1>(result, inverse_digit,
                                               key[row]);
        }
        else {
            if (overwrite && row == 0)
                MulInFD_Multi<P::n, P::k + 1>(
                    result, input.decomposition[row], key[row]);
            else
                FMAInFD_Multi<P::n, P::k + 1>(
                    result, input.decomposition[row], key[row]);
        }
    }
}

// Algorithm 2 (and Algorithm 7 when Positive=false).  The expensive gadget
// decomposition and forward Fourier transforms are performed once per input
// ciphertext and digit; the four selector products stay in the Fourier domain
// until the output sum is complete.
template <class P>
void BatchEMP(std::vector<TRLWE<P>> &ciphertexts,
              const BatchEMPKey<P> &key, const bool positive,
              BatchEMPWorkspace<P> &workspace)
{
    const std::uint32_t slots = ciphertexts.size();
    if (!std::has_single_bit(slots) || slots > P::n)
        throw std::invalid_argument(
            "BatchEMP: ciphertext count must be a power of two no larger than N");
    const std::uint32_t expected_digits =
        (std::bit_width(slots - 1U) + 1U) / 2U;
    if (key.digits.size() != expected_digits)
        throw std::invalid_argument("BatchEMP: key has the wrong digit count");

    workspace.prepare(slots);
    auto &hoisted = workspace.hoisted;
    auto &transformed = workspace.transformed;

    for (std::uint32_t level = 0; level < key.digits.size(); level++) {
        for (std::uint32_t i = 0; i < slots; i++)
            BatchHoist<P>(hoisted[i], ciphertexts[i]);

        const std::uint32_t step = 1U << (2 * level);
        const auto &digit_key = key.digits[level];
        for (std::uint32_t output = 0; output < slots; output++) {
            for (std::uint32_t row = 0; row < 3; row++) {
                const std::uint32_t shift = (row + 1) * step;
                std::uint32_t source;
                bool inverse;
                if (positive) {
                    inverse = output < shift;
                    source = (output + slots - (shift % slots)) % slots;
                }
                else {
                    inverse = output + shift >= slots;
                    source = (output + (shift % slots)) % slots;
                }
                BatchExternalProductAddFD<P>(
                    transformed[output], hoisted[source],
                    inverse ? digit_key.inverse[row]
                            : digit_key.direct[row],
                    inverse, row == 0);
            }

            BatchExternalProductAddFD<P>(transformed[output],
                                         hoisted[output], digit_key.direct[3],
                                         false);
        }

        for (std::uint32_t i = 0; i < slots; i++)
            for (std::uint32_t component = 0; component <= P::k; component++)
                TwistFFT<P>(ciphertexts[i][component],
                            transformed[i][component]);
    }
}

template <class P>
void BatchEMP(std::vector<TRLWE<P>> &ciphertexts,
              const BatchEMPKey<P> &key, const bool positive)
{
    BatchEMPWorkspace<P> workspace;
    BatchEMP<P>(ciphertexts, key, positive, workspace);
}

template <class T>
std::uint32_t BatchModSwitchTorus(const T value,
                                  const std::uint32_t modulus)
{
    if (modulus < 2 || !std::has_single_bit(modulus))
        throw std::invalid_argument(
            "BatchModSwitchTorus: modulus must be a power of two of at least 2");
    const std::uint32_t modulus_bits = std::countr_zero(modulus);
    constexpr std::uint32_t width = std::numeric_limits<T>::digits;
    if (modulus_bits >= width)
        throw std::invalid_argument(
            "BatchModSwitchTorus: modulus does not fit the torus");
    const std::uint32_t shift = width - modulus_bits;
    const T rounding = static_cast<T>(1) << (shift - 1);
    return static_cast<std::uint32_t>((value + rounding) >> shift) &
           (modulus - 1);
}

template <class DomainP>
using BatchModSwitchResult =
    std::array<std::vector<std::uint32_t>, DomainP::k + 1>;

template <class DomainP>
BatchModSwitchResult<DomainP> BatchModSwitch(
    const TRLWE<DomainP> &input, const std::uint32_t modulus)
{
    BatchModSwitchResult<DomainP> result;
    for (std::uint32_t component = 0; component <= DomainP::k; component++) {
        result[component].resize(DomainP::n);
        for (std::uint32_t i = 0; i < DomainP::n; i++)
            result[component][i] =
                BatchModSwitchTorus(input[component][i], modulus);
    }
    return result;
}

template <class P>
void BatchMultiplyByMonomial(TRLWE<P> &result, const TRLWE<P> &input,
                             const std::uint32_t exponent)
{
    for (std::uint32_t component = 0; component <= P::k; component++)
        PolynomialMulByXai<P>(result[component], input[component], exponent);
}

template <class DomainP, class TargetP>
void BatchBootAccumulateModSwitched(
    std::vector<TRLWE<TargetP>> &result,
    const BatchModSwitchResult<DomainP> &modswitched,
    const std::span<const Polynomial<TargetP>> test_vectors,
    const BatchBootKey<DomainP, TargetP> &boot_key)
{
    const std::uint32_t slots = DomainP::n;
    if (test_vectors.size() != 1 && test_vectors.size() != slots)
        throw std::invalid_argument(
            "BatchBootAccumulate: provide either one or one-per-slot test vector");

    result.resize(slots);
    for (std::uint32_t i = 0; i < slots; i++) {
        result[i] = {};
        const auto &test = test_vectors.size() == 1 ? test_vectors[0]
                                                    : test_vectors[i];
        PolynomialMulByXai<TargetP>(
            result[i][TargetP::k], test,
            modswitched[DomainP::k][i] % (2 * TargetP::n));
    }

    TRLWE<TargetP> rotated;
    BatchEMPWorkspace<TargetP> workspace;
    for (std::uint32_t component = 0; component < DomainP::k; component++) {
        const auto &component_key = boot_key.components[component];
        if (!component_key.final_positive) continue;

        for (const auto &gap_key : component_key.negative_gaps) {
            BatchEMP<TargetP>(result, gap_key, false, workspace);
            for (std::uint32_t i = 0; i < slots; i++) {
                const std::uint32_t a = modswitched[component][i];
                const std::uint32_t exponent =
                    a == 0 ? 0 : 2 * TargetP::n - a;
                BatchMultiplyByMonomial<TargetP>(rotated, result[i], exponent);
                result[i] = rotated;
            }
        }
        BatchEMP<TargetP>(result, *component_key.final_positive, true,
                          workspace);
    }
}

template <class DomainP, class TargetP>
void BatchBootAccumulate(
    std::vector<TRLWE<TargetP>> &result, const TRLWE<DomainP> &input,
    const std::span<const Polynomial<TargetP>> test_vectors,
    const BatchBootKey<DomainP, TargetP> &boot_key,
    const std::uint32_t exponent_bias = 0)
{
    auto modswitched = BatchModSwitch<DomainP>(input, 2 * TargetP::n);
    for (auto &exponent : modswitched[DomainP::k])
        exponent = (exponent + exponent_bias) % (2 * TargetP::n);
    BatchBootAccumulateModSwitched<DomainP, TargetP>(
        result, modswitched, test_vectors, boot_key);
}

template <class DomainP, class TargetP>
void BatchBootAccumulate(
    std::vector<TRLWE<TargetP>> &result, const TRLWE<DomainP> &input,
    const Polynomial<TargetP> &test_vector,
    const BatchBootKey<DomainP, TargetP> &boot_key,
    const std::uint32_t exponent_bias = 0)
{
    BatchBootAccumulate<DomainP, TargetP>(
        result, input, std::span<const Polynomial<TargetP>>(&test_vector, 1),
        boot_key, exponent_bias);
}

template <class P>
void BatchHomomorphicTrace(TRLWE<P> &result, const TRLWE<P> &input,
                           const std::uint32_t subring_degree,
                           const AnnihilateKey<P> &automorphism_keys)
{
    if (!std::has_single_bit(subring_degree) || subring_degree > P::n)
        throw std::invalid_argument(
            "BatchHomomorphicTrace: invalid subring degree");

    result = input;
    const std::uint32_t first_level =
        std::countr_zero(subring_degree) + 1;
    for (std::uint32_t level = first_level; level <= P::nbit; level++) {
        for (auto &poly : result)
            for (auto &coefficient : poly) coefficient /= 2;
        TRLWE<P> automorphism;
        EvalAuto<P>(automorphism, result, (1U << level) + 1,
                    automorphism_keys[level - 1]);
        TRLWEAdd<P>(result, result, automorphism);
    }
}

// Algorithm 6, including the sparse trace when ciphertexts.size() < N.
template <class P>
void BatchRepack(TRLWE<P> &result,
                 const std::span<const TRLWE<P>> ciphertexts,
                 const AnnihilateKey<P> &automorphism_keys)
{
    const std::uint32_t slots = ciphertexts.size();
    if (!std::has_single_bit(slots) || slots > P::n)
        throw std::invalid_argument(
            "BatchRepack: ciphertext count must be a power of two no larger than N");

    std::vector<TRLWE<P>> work(ciphertexts.begin(), ciphertexts.end());
    const std::uint32_t levels = std::countr_zero(slots);
    for (std::uint32_t level = 1; level <= levels; level++) {
        const std::uint32_t remaining = slots >> level;
        const std::uint32_t exponent = P::n >> level;
        for (std::uint32_t j = 0; j < remaining; j++) {
            TRLWE<P> rotated;
            BatchMultiplyByMonomial<P>(rotated, work[j + remaining], exponent);
            TRLWE<P> difference;
            TRLWESub<P>(difference, work[j], rotated);
            for (auto &poly : difference)
                for (auto &coefficient : poly) coefficient /= 2;

            TRLWE<P> automorphism;
            EvalAuto<P>(automorphism, difference, (1U << level) + 1,
                        automorphism_keys[level - 1]);
            TRLWEAdd<P>(work[j], difference, automorphism, rotated);
        }
    }

    if (slots == P::n)
        result = work[0];
    else
        BatchHomomorphicTrace<P>(result, work[0], slots,
                                 automorphism_keys);
}

template <class DomainP, class TargetP>
void BatchBootstrap(
    TRLWE<TargetP> &result, const TRLWE<DomainP> &input,
    const std::span<const Polynomial<TargetP>> test_vectors,
    const BatchBootKey<DomainP, TargetP> &boot_key,
    const AnnihilateKey<TargetP> &automorphism_keys,
    const std::uint32_t exponent_bias = 0)
{
    std::vector<TRLWE<TargetP>> accumulators;
    BatchBootAccumulate<DomainP, TargetP>(
        accumulators, input, test_vectors, boot_key, exponent_bias);
    BatchRepack<TargetP>(result, accumulators, automorphism_keys);
}

template <class DomainP, class TargetP>
void BatchBootstrap(
    TRLWE<TargetP> &result, const TRLWE<DomainP> &input,
    const Polynomial<TargetP> &test_vector,
    const BatchBootKey<DomainP, TargetP> &boot_key,
    const AnnihilateKey<TargetP> &automorphism_keys,
    const std::uint32_t exponent_bias = 0)
{
    BatchBootstrap<DomainP, TargetP>(
        result, input,
        std::span<const Polynomial<TargetP>>(&test_vector, 1), boot_key,
        automorphism_keys, exponent_bias);
}

template <class P, std::uint32_t InputBits, std::uint32_t OutputBits>
struct BatchBootTestVector {
    Polynomial<P> polynomial;
    static constexpr std::uint32_t exponent_bias = P::n >> InputBits;
};

// Build Algorithm 1's anti-cyclic test polynomial.  BatchBoot reserves one
// input padding bit, so the supported message range is [0, 2^(InputBits-1)).
template <class P, std::uint32_t InputBits,
          std::uint32_t OutputBits = InputBits, class Function>
BatchBootTestVector<P, InputBits, OutputBits> MakeBatchBootTestVector(
    Function &&function)
{
    static_assert(InputBits > 0 && OutputBits > 0);
    static_assert(InputBits <= P::nbit);
    static_assert(InputBits < std::numeric_limits<std::uint32_t>::digits);
    static_assert(OutputBits < std::numeric_limits<std::uint32_t>::digits);
    static_assert(OutputBits < std::numeric_limits<typename P::T>::digits);
    constexpr std::uint32_t input_modulus = 1U << InputBits;
    constexpr std::uint32_t input_domain = input_modulus / 2;
    constexpr std::uint32_t output_modulus = 1U << OutputBits;
    constexpr std::uint32_t width =
        std::numeric_limits<typename P::T>::digits;

    BatchBootTestVector<P, InputBits, OutputBits> result{};
    for (std::uint32_t j = 0; j < P::n; j++) {
        // The bias moves each valid input to the center of a width-2N/p
        // interval.  Use the interval index here (rather than rounding to the
        // closest integer), or a one-unit modulus-switch error at the center
        // would select the adjacent LUT entry.
        const std::uint32_t interval = static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(j + 1) * input_modulus) /
            (2 * P::n));
        const std::uint32_t rounded_input =
            interval < input_domain ? interval : input_domain - 1;
        const std::uint32_t output =
            static_cast<std::uint32_t>(function(rounded_input)) %
            output_modulus;
        const typename P::T encoded =
            static_cast<typename P::T>(output) << (width - OutputBits);
        result.polynomial[P::n - j - 1] = -encoded;
    }
    return result;
}

template <std::uint32_t PlainBits, class T>
std::uint32_t BatchTorusDecode(const T value)
{
    static_assert(PlainBits > 0);
    static_assert(PlainBits < std::numeric_limits<std::uint32_t>::digits);
    constexpr std::uint32_t width = std::numeric_limits<T>::digits;
    static_assert(PlainBits < width);
    const T rounding = static_cast<T>(1) << (width - PlainBits - 1);
    return static_cast<std::uint32_t>((value + rounding) >>
                                      (width - PlainBits)) &
           ((1U << PlainBits) - 1);
}

// Optional ring key switch used by Algorithm 3 when the accumulator/repacking
// key differs from the sparse key used by the surrounding computation.
template <class P>
using BatchRLWEKeySwitchKey = std::array<HalfTRGSWFFT<P>, P::k>;

template <class P>
void BatchRLWEKeySwitchKeyGen(BatchRLWEKeySwitchKey<P> &result,
                              const Key<P> &source_key,
                              const Key<P> &target_key)
{
    for (std::uint32_t component = 0; component < P::k; component++) {
        Polynomial<P> part_key;
        for (std::uint32_t i = 0; i < P::n; i++)
            part_key[i] = source_key[component * P::n + i];
        halftrgswSymEncrypt<P>(result[component], part_key, target_key);
    }
}

template <class P>
void BatchRLWEKeySwitch(TRLWE<P> &result, const TRLWE<P> &input,
                        const BatchRLWEKeySwitchKey<P> &key)
{
    result = {};
    result[P::k] = input[P::k];
    for (std::uint32_t component = 0; component < P::k; component++) {
        TRLWE<P> product;
        ExternalProduct<P>(product, input[component], key[component]);
        TRLWESub<P>(result, result, product);
    }
}

// Scheme-switching key and operation used by Algorithm 4.  This is the k=1
// EP-CBS conversion RLWE_s(m) -> RLWE_s(-s*m).
template <class P>
using BatchSchemeSwitchKey = HalfTRGSWFFT<P>;

template <class P>
void BatchSchemeSwitchKeyGen(BatchSchemeSwitchKey<P> &result,
                             const Key<P> &key)
{
    static_assert(P::k == 1,
                  "BatchBoot scheme switching currently supports k=1");
    Polynomial<P> part_key;
    for (std::uint32_t i = 0; i < P::n; i++) part_key[i] = key[i];
    Polynomial<P> key_square;
    PolyMulNaive<P>(key_square, part_key, part_key);
    halftrgswSymEncrypt<P>(result, key_square, key);
}

template <class P>
void BatchSchemeSwitch(TRLWE<P> &result, const TRLWE<P> &input,
                       const BatchSchemeSwitchKey<P> &key)
{
    static_assert(P::k == 1,
                  "BatchBoot scheme switching currently supports k=1");
    ExternalProduct<P>(result, input[0], key);
    for (std::uint32_t i = 0; i < P::n; i++) result[0][i] += input[1][i];
}

// Algorithm 4, first stage (Half.BatchCBS).  The returned TRGSWs encrypt
// X^m for every active input coefficient.  TargetP::l is the output gadget
// length d and therefore must be a power of two.
template <class DomainP, class TargetP>
void BatchHalfCircuitBootstrap(
    std::vector<TRGSW<TargetP>> &result, const TRLWE<DomainP> &input,
    const BatchBootKey<DomainP, TargetP> &boot_key,
    const AnnihilateKey<TargetP> &automorphism_keys,
    const BatchSchemeSwitchKey<TargetP> &scheme_switch_key)
{
    static_assert(TargetP::k == 1);
    static_assert(TargetP::lₐ == TargetP::l,
                  "BatchHalfCircuitBootstrap requires matching gadget lengths");
    static_assert(std::has_single_bit(TargetP::l),
                  "the circuit-bootstrap gadget length must be a power of two");
    static_assert(TargetP::l <= TargetP::n,
                  "the circuit-bootstrap gadget length must not exceed N");

    constexpr std::uint32_t gadget_length = TargetP::l;
    constexpr auto gadget = hgen<TargetP, false>();
    auto modswitched =
        BatchModSwitch<DomainP>(input, 2 * TargetP::n / gadget_length);
    for (auto &component : modswitched)
        for (auto &exponent : component)
            exponent = (exponent * gadget_length) % (2 * TargetP::n);

    Polynomial<TargetP> initial_accumulator{};
    for (std::uint32_t i = 0; i < gadget_length; i++)
        initial_accumulator[i] = gadget[i];
    std::vector<TRLWE<TargetP>> accumulators;
    BatchBootAccumulateModSwitched<DomainP, TargetP>(
        accumulators, modswitched,
        std::span<const Polynomial<TargetP>>(&initial_accumulator, 1),
        boot_key);

    result.resize(DomainP::n);
    for (std::uint32_t slot = 0; slot < DomainP::n; slot++) {
        for (std::uint32_t level = 0; level < gadget_length; level++) {
            TRLWE<TargetP> rotated;
            BatchMultiplyByMonomial<TargetP>(
                rotated, accumulators[slot],
                level == 0 ? 0 : 2 * TargetP::n - level);
            auto &main_row =
                result[slot][TargetP::k * TargetP::lₐ + level];
            BatchHomomorphicTrace<TargetP>(
                main_row, rotated, TargetP::n / gadget_length,
                automorphism_keys);
            BatchSchemeSwitch<TargetP>(result[slot][level], main_row,
                                       scheme_switch_key);
        }
    }
}

template <class DomainP, class TargetP>
void BatchHalfCircuitBootstrapFFT(
    std::vector<TRGSWFFT<TargetP>> &result, const TRLWE<DomainP> &input,
    const BatchBootKey<DomainP, TargetP> &boot_key,
    const AnnihilateKey<TargetP> &automorphism_keys,
    const BatchSchemeSwitchKey<TargetP> &scheme_switch_key)
{
    std::vector<TRGSW<TargetP>> coefficient_result;
    BatchHalfCircuitBootstrap<DomainP, TargetP>(
        coefficient_result, input, boot_key, automorphism_keys,
        scheme_switch_key);
    result.resize(coefficient_result.size());
    for (std::uint32_t i = 0; i < coefficient_result.size(); i++)
        ApplyFFT2trgsw<TargetP>(result[i], coefficient_result[i]);
}

template <class P>
void BatchTRLWEMulByPolynomial(TRLWE<P> &result, const TRLWE<P> &input,
                               const Polynomial<P> &polynomial)
{
    for (std::uint32_t component = 0; component <= P::k; component++)
        PolyMul<P>(result[component], input[component], polynomial);
}

// Algorithm 8.  `test_polynomials` are ordered in radix groups exactly as in
// the paper and must contain radix^(controls.size()-1) entries.
template <class P>
void BatchExternalProductTree(
    TRLWE<P> &result,
    const std::span<const TRGSWFFT<P>> monomial_controls,
    const std::span<const Polynomial<P>> test_polynomials,
    const std::uint32_t radix,
    const AnnihilateKey<P> &automorphism_keys)
{
    if (monomial_controls.empty())
        throw std::invalid_argument(
            "BatchExternalProductTree: at least one control is required");
    if (!std::has_single_bit(radix) || radix < 2 || radix > P::n)
        throw std::invalid_argument(
            "BatchExternalProductTree: radix must be a power of two in [2, N]");

    std::size_t expected = 1;
    for (std::size_t i = 1; i < monomial_controls.size(); i++) {
        if (expected > std::numeric_limits<std::size_t>::max() / radix)
            throw std::invalid_argument(
                "BatchExternalProductTree: test-polynomial count overflows");
        expected *= radix;
    }
    if (test_polynomials.size() != expected)
        throw std::invalid_argument(
            "BatchExternalProductTree: incorrect test-polynomial count");

    std::vector<TRLWE<P>> current(test_polynomials.size());
    for (std::size_t i = 0; i < current.size(); i++) {
        current[i] = {};
        current[i][P::k] = test_polynomials[i];
    }

    Polynomial<P> normalization{};
    for (std::uint32_t i = 0; i < 2 * P::n / radix; i++)
        normalization[i] = -static_cast<typename P::T>(1);

    for (std::size_t level = 0; level < monomial_controls.size(); level++) {
        for (auto &ciphertext : current) {
            TRLWE<P> product;
            ExternalProduct<P>(product, ciphertext,
                               monomial_controls[level]);
            ciphertext = product;
        }
        if (level + 1 == monomial_controls.size()) break;

        const std::size_t groups = current.size() / radix;
        std::vector<TRLWE<P>> packed(groups);
        for (std::size_t group = 0; group < groups; group++) {
            BatchRepack<P>(
                packed[group],
                std::span<const TRLWE<P>>(current.data() + group * radix,
                                          radix),
                automorphism_keys);
            TRLWE<P> normalized;
            BatchTRLWEMulByPolynomial<P>(normalized, packed[group],
                                         normalization);
            packed[group] = normalized;
        }
        current = std::move(packed);
    }
    result = current[0];
}

// Algorithm 4 for one output sub-function.  The first `input_arity`
// coefficient slots become the controls of Algorithm 8; callers evaluating
// multiple output sub-functions can reuse BatchHalfCircuitBootstrapFFT and
// invoke BatchExternalProductTree once per output instead.
template <class DomainP, class TargetP>
void BatchCircuitBootstrap(
    TRLWE<TargetP> &result, const TRLWE<DomainP> &input,
    const BatchBootKey<DomainP, TargetP> &boot_key,
    const AnnihilateKey<TargetP> &automorphism_keys,
    const BatchSchemeSwitchKey<TargetP> &scheme_switch_key,
    const std::uint32_t input_arity,
    const std::span<const Polynomial<TargetP>> test_polynomials,
    const std::uint32_t radix)
{
    if (input_arity == 0 || input_arity > DomainP::n)
        throw std::invalid_argument(
            "BatchCircuitBootstrap: invalid input arity");

    std::vector<TRGSWFFT<TargetP>> controls;
    BatchHalfCircuitBootstrapFFT<DomainP, TargetP>(
        controls, input, boot_key, automorphism_keys, scheme_switch_key);
    BatchExternalProductTree<TargetP>(
        result,
        std::span<const TRGSWFFT<TargetP>>(controls.data(), input_arity),
        test_polynomials, radix, automorphism_keys);
}

}  // namespace TFHEpp
