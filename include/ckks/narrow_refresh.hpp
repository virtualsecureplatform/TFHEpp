#pragma once

#include "narrow_bootstrapping.hpp"
#include "narrow_frequency_basis.hpp"
#include "narrow_approximation.hpp"
#include "compact_linear.hpp"
#include <cereal/types/optional.hpp>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// SHIP Section 3.2 and ePrint 2026/2119 Section 2.3: a single low-level
// matrix S2C, sparse-secret encapsulation, and narrow decryption in the
// exponent. The gap gamma=q0/input_scale controls SHIP's sine approximation.
template <class NarrowSchedule, std::uint32_t InputLogDelta,
          std::uint32_t S2CPlainLogDelta, int BSGSStep = 0, bool CompactS2C = false,
          std::uint32_t S2CAuxiliaryLogScale = NarrowSchedule::auxiliary_log_scale>
struct CKKSNarrowRefreshSchedule {
    using Narrow = NarrowSchedule;
    using Param = typename Narrow::Param;
    static constexpr std::uint32_t input_log_delta = InputLogDelta;
    static constexpr std::uint32_t coefficient_log_q = Narrow::input_log_q;
    static constexpr std::uint32_t input_log_q = coefficient_log_q + S2CPlainLogDelta;
    static constexpr std::uint32_t output_log_q = Narrow::output_log_q;
    static constexpr std::uint32_t output_log_delta = Narrow::log_delta;
    static constexpr std::uint32_t s2c_plain_log_delta = S2CPlainLogDelta;
    static_assert(InputLogDelta < coefficient_log_q);
    static_assert(input_log_q <= std::numeric_limits<typename Param::T>::digits);
    static_assert(output_log_q > input_log_q);
    static_assert(S2CPlainLogDelta > 0 && S2CPlainLogDelta <= 50);
    static constexpr std::uint32_t log_gap = coefficient_log_q - InputLogDelta;
    static constexpr int bsgs_step = BSGSStep > 0 ? BSGSStep :
        (1 << ((Param::nbit - 1) / 2));
    using Input = CKKSCiphertext<Param, input_log_q, InputLogDelta>;
    using Coefficients = CKKSCiphertext<Param, coefficient_log_q, InputLogDelta>;
    using Output = typename Narrow::Output;
    // Lift the entire linear transform and encapsulation, then divide once
    // under the sparse secret. Fresh input error is preserved by this lift.
    // S2C and blind-rotation switches consume different scales and levels.
    // Keep the historical shared setting by default, but allow independent
    // budgets without lifting S2C keys solely for blind-rotation precision.
    static constexpr auto auxiliary_log_scale = S2CAuxiliaryLogScale;
    static constexpr auto evaluation_input_log_q = input_log_q + auxiliary_log_scale;
    static constexpr auto evaluation_coefficient_log_q = coefficient_log_q + auxiliary_log_scale;
    static constexpr auto evaluation_input_log_delta = InputLogDelta + auxiliary_log_scale;
    static_assert(evaluation_input_log_q <= std::numeric_limits<typename Param::T>::digits);
    using EvaluationInput = CKKSCiphertext<Param, evaluation_input_log_q, evaluation_input_log_delta>;
    using EvaluationCoefficients = CKKSCiphertext<Param, evaluation_coefficient_log_q, evaluation_input_log_delta>;
    static constexpr bool compact_s2c = CompactS2C;
    using S2CPlan = std::conditional_t<CompactS2C,
        CKKSCompactRealLinearTransformPlan<Param, evaluation_input_log_q, evaluation_input_log_delta, S2CPlainLogDelta>,
        CKKSRealLinearTransformPlan<Param, evaluation_input_log_q, evaluation_input_log_delta, S2CPlainLogDelta>>;

    static double gap() { return std::ldexp(1.0, log_gap); }

    // |gamma/(2*pi)*sin(2*pi*x/gamma)-x| <= (2*pi)^2*M^3/(6*gamma^2).
    // This bounds analytic sine error only; encryption and evaluation error
    // must be measured or bounded separately for an instantiated parameter set.
    static double sine_error_bound(double message_bound)
    {
        if (!std::isfinite(message_bound) || message_bound < 0)
            throw std::invalid_argument("invalid narrow refresh message bound");
        const double ratio = 2 * std::acos(-1.0) * message_bound / gap();
        return message_bound * ratio * ratio / 6;
    }
};

template <class RefreshSchedule>
struct CKKSNarrowRefreshKey {
    using P = typename RefreshSchedule::Param;
    typename RefreshSchedule::S2CPlan s2c_plan;
    CKKSGaloisKey<P, RefreshSchedule::evaluation_input_log_q> s2c_input_galois;
    CKKSGaloisKey<P, RefreshSchedule::evaluation_coefficient_log_q> s2c_output_galois;
    CKKSSecretKeySwitchKey<P, RefreshSchedule::evaluation_coefficient_log_q> to_sparse;
    CKKSNarrowBootstrapKey<typename RefreshSchedule::Narrow> narrow;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(s2c_plan.has_direct, s2c_plan.has_conjugate,
            s2c_plan.direct, s2c_plan.conjugate, s2c_input_galois,
            s2c_output_galois, to_sparse, narrow);
    }
};

template <class NarrowSchedule>
inline void CKKSNarrowStructuredKeyGen(Key<typename NarrowSchedule::Param> &key)
{
    using P = typename NarrowSchedule::Param;
    key.fill(typename P::T{0});
    std::uniform_int_distribution<int> sign(0, 1);
    for (std::size_t x = 0; x < NarrowSchedule::blocks; x++) {
        std::uniform_int_distribution<std::size_t> position(0, NarrowSchedule::block_size(x) - 1);
        key[x * NarrowSchedule::block_width + position(generator)] =
            sign(generator) ? typename P::T{1} : -typename P::T{1};
    }
}

namespace narrow_detail {
template <class RefreshSchedule>
inline void generateRefreshLinearKeys(CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const Key<typename RefreshSchedule::Param> &input_secret,
    const Key<typename RefreshSchedule::Param> &structured_secret, CKKSNoise noise)
{
    using P = typename RefreshSchedule::Param;
    narrow_detail::structuredPositions<typename RefreshSchedule::Narrow>(structured_secret);
    CKKSBuildPackedSlotToCoeffPlan<P, RefreshSchedule::evaluation_input_log_q,
        RefreshSchedule::evaluation_input_log_delta, RefreshSchedule::s2c_plain_log_delta>(
        key.s2c_plan, RefreshSchedule::bsgs_step);
    CKKSGaloisKeyGen<P, RefreshSchedule::evaluation_input_log_q>(key.s2c_input_galois, input_secret, noise);
    CKKSGaloisKeyGen<P, RefreshSchedule::evaluation_coefficient_log_q>(key.s2c_output_galois, input_secret, noise);
    CKKSSecretKeySwitchKeyGen<P, RefreshSchedule::evaluation_coefficient_log_q>(key.to_sparse,
        input_secret, structured_secret, noise);
}
}  // namespace narrow_detail

// Reusable S2C, encapsulation, and narrow evaluation setup for a compressed
// selector key. Reconstruct key.narrow from a public CKKSNarrowPackedKey
// before refreshing. Expanded selectors are never generated by this function.
template <class RefreshSchedule>
inline void CKKSNarrowRefreshSetupKeyGen(
    CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const Key<typename RefreshSchedule::Param> &input_secret,
    const Key<typename RefreshSchedule::Param> &structured_secret,
    const Key<typename RefreshSchedule::Param> &evaluation_secret,
    CKKSNoise noise = {RefreshSchedule::Param::α, 0},
    CKKSNarrowCacheMode cache_mode = CKKSNarrowCacheMode::Prepared)
{
    narrow_detail::generateRefreshLinearKeys<RefreshSchedule>(key, input_secret, structured_secret, noise);
    CKKSNarrowSetupKeyGen<typename RefreshSchedule::Narrow>(key.narrow,
        evaluation_secret, noise, cache_mode);
}

// input/evaluation keys may be independent. Reusing evaluation as the input
// key gives a repeatable refresh cycle under the same dense secret.
template <class RefreshSchedule>
inline void CKKSNarrowRefreshKeyGen(
    CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const Key<typename RefreshSchedule::Param> &input_secret,
    const Key<typename RefreshSchedule::Param> &structured_secret,
    const Key<typename RefreshSchedule::Param> &evaluation_secret,
    CKKSNoise noise = {RefreshSchedule::Param::α, 0},
    CKKSNarrowCacheMode cache_mode = CKKSNarrowCacheMode::Prepared)
{
    narrow_detail::generateRefreshLinearKeys<RefreshSchedule>(key, input_secret, structured_secret, noise);
    CKKSNarrowBootstrapKeyGen<typename RefreshSchedule::Narrow>(key.narrow,
        structured_secret, evaluation_secret, noise, cache_mode);
}

namespace narrow_detail {

// Keep derived rotation/conjugation spectra local to one S2C evaluation.
// BSGS releases them through clear_cached_keys at its existing boundaries.
template <class P, std::uint32_t Q, std::uint32_t Delta,
          std::uint32_t PlainDelta, class Result, class Input, class Plan>
inline void realTransformWithPreparedGalois(
    Result &result, const Input &input, const Plan &plan,
    const CKKSGaloisKey<P, Q> &input_key,
    const CKKSGaloisKey<P, Q - PlainDelta> &output_key)
{
    constexpr auto secondary_bits = CKKSPreparedAutoSecondaryBits<P, Q>();
    constexpr auto rows = CKKSKeySwitchRowCountForLevel<P, Q>();
    if constexpr (use_multilimb_digit_fft_v<P> &&
        P::B̅gbit + secondary_bits + P::nbit + std::bit_width(rows) + 3 <
            std::numeric_limits<double>::digits) {
        CKKSLazyPreparedGaloisKeyFFT<P, Q> prepared_input(input_key);
        CKKSLazyPreparedGaloisKeyFFT<P, Q - PlainDelta> prepared_output(output_key);
        CKKSRealLinearTransform<P, Q, Delta, PlainDelta>(
            result, input, plan, prepared_input, prepared_output);
    }
    else {
        CKKSRealLinearTransform<P, Q, Delta, PlainDelta>(
            result, input, plan, input_key, output_key);
    }
}
} // namespace narrow_detail

template <class RefreshSchedule>
inline void CKKSNarrowPrepareCoefficients(
    typename RefreshSchedule::Coefficients &coefficients,
    const typename RefreshSchedule::Input &input,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key)
{
    using P = typename RefreshSchedule::Param;
    if constexpr (RefreshSchedule::auxiliary_log_scale != 0) {
        auto lifted = std::make_unique<typename RefreshSchedule::EvaluationInput>();
        for (std::size_t c = 0; c <= P::k; c++)
            for (std::size_t i = 0; i < P::n; i++)
                lifted->ct[c][i] = ckks_detail::reduceToLevel<P, RefreshSchedule::input_log_q>(input.ct[c][i])
                    << RefreshSchedule::auxiliary_log_scale;
        auto dense = std::make_unique<typename RefreshSchedule::EvaluationCoefficients>();
        auto sparse = std::make_unique<typename RefreshSchedule::EvaluationCoefficients>();
        narrow_detail::realTransformWithPreparedGalois<P, RefreshSchedule::evaluation_input_log_q,
            RefreshSchedule::evaluation_input_log_delta, RefreshSchedule::s2c_plain_log_delta>(
            *dense, *lifted, key.s2c_plan, key.s2c_input_galois, key.s2c_output_galois);
        CKKSSecretKeySwitch<P>(*sparse, *dense, key.to_sparse);
        CKKSScaleDown<P, RefreshSchedule::evaluation_coefficient_log_q,
            RefreshSchedule::evaluation_input_log_delta, RefreshSchedule::auxiliary_log_scale>(coefficients, *sparse);
    }
    else {
        auto dense_storage = std::make_unique<typename RefreshSchedule::Coefficients>();
        auto &dense = *dense_storage;
        narrow_detail::realTransformWithPreparedGalois<P, RefreshSchedule::input_log_q,
            RefreshSchedule::input_log_delta, RefreshSchedule::s2c_plain_log_delta>(
            dense, input, key.s2c_plan, key.s2c_input_galois, key.s2c_output_galois);
        CKKSSecretKeySwitch<P>(coefficients, dense, key.to_sparse);
    }
    // A materialized ModRaise is unnecessary: narrow exponentiation uses
    // omega^q0=1 and therefore discards exactly the q0*I introduced by lifting
    // the public coefficients. The evaluator never computes the secret phase.
}

namespace narrow_detail {

// X^-N/2 moves the imaginary coefficient half into the first half. Unlike
// slot rotation, this multiplies the public ciphertext polynomials by a
// monomial and needs no automorphism or key switch.
template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void imaginaryCoefficientHalf(
    CKKSCiphertext<P, LogQ, LogDelta> &result,
    const CKKSCiphertext<P, LogQ, LogDelta> &input)
{
    for (std::size_t c = 0; c <= P::k; c++)
        for (std::size_t i = 0; i < P::n / 2; i++) {
            result.ct[c][i] = input.ct[c][i + P::n / 2];
            result.ct[c][i + P::n / 2] = -input.ct[c][i];
        }
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(result.ct);
}

// X^(N/2) evaluates to i at every canonical CKKS slot (5^j == 1 mod 4).
// Thus real+i*imag can be combined exactly with no plaintext rescale level.
template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void combineRealImaginary(CKKSCiphertext<P, LogQ, LogDelta> &result,
    const CKKSCiphertext<P, LogQ, LogDelta> &real,
    const CKKSCiphertext<P, LogQ, LogDelta> &imaginary)
{
    for (std::size_t c = 0; c <= P::k; c++)
        for (std::size_t i = 0; i < P::n / 2; i++) {
            // Read both halves before writing, including when result aliases
            // either input. Avoid a full ciphertext temporary at large N.
            const auto low = real.ct[c][i] - imaginary.ct[c][i + P::n / 2];
            const auto high = real.ct[c][i + P::n / 2] + imaginary.ct[c][i];
            result.ct[c][i] = low;
            result.ct[c][i + P::n / 2] = high;
        }
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(result.ct);
}

}  // namespace narrow_detail

// SHIP's standard small-angle refresh. Both coefficient halves use the same
// narrow selectors and DD product keys, preserving full complex slot order.
template <class RefreshSchedule>
inline void CKKSNarrowRefreshPrepared(
    typename RefreshSchedule::Output &result,
    const typename RefreshSchedule::Coefficients &coefficients,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key)
{
    using P = typename RefreshSchedule::Param;
    using Narrow = typename RefreshSchedule::Narrow;
    auto imaginary_storage = std::make_unique<typename RefreshSchedule::Coefficients>();
    auto &imaginary_coefficients = *imaginary_storage;
    narrow_detail::imaginaryCoefficientHalf<P>(imaginary_coefficients, coefficients);
    const std::complex<double> weight{0, -RefreshSchedule::gap() / (4 * std::acos(-1.0))};
    using Request = CKKSNarrowBlindRotationRequest<Narrow, RefreshSchedule::input_log_delta>;
    const std::array<Request, 2> requests = {
        Request{&coefficients, 1, weight}, Request{&imaginary_coefficients, 1, weight}};
    std::vector<typename RefreshSchedule::Output> roots;
    CKKSNarrowBlindRotateBatch<Narrow, RefreshSchedule::input_log_delta>(roots, requests, key.narrow);
    auto &real = roots[0];
    auto &imaginary = roots[1];
    auto conjugated_storage = std::make_unique<typename RefreshSchedule::Output>();
    auto &conjugated = *conjugated_storage;
    narrow_detail::conjugate<Narrow, false>(conjugated.ct, real.ct, key.narrow);
    CKKSAddInPlace<P>(real, conjugated);
    narrow_detail::conjugate<Narrow, false>(conjugated.ct, imaginary.ct, key.narrow);
    CKKSAddInPlace<P>(imaginary, conjugated);
    narrow_detail::combineRealImaginary<P>(result, real, imaginary);
}

// Evaluate from slots, preparing S2C and sparse encapsulation once. The
// prepared entry point above accepts the same coefficient ciphertext when
// callers also need coefficient diagnostics or multiple evaluations.
template <class RefreshSchedule>
inline void CKKSNarrowRefresh(
    typename RefreshSchedule::Output &result,
    const typename RefreshSchedule::Input &input,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key)
{
    auto prepared = std::make_unique<typename RefreshSchedule::Coefficients>();
    CKKSNarrowPrepareCoefficients<RefreshSchedule>(*prepared, input, key);
    CKKSNarrowRefreshPrepared<RefreshSchedule>(result, *prepared, key);
}

// Section 7's full FBS pipeline. The offline series represents f(gamma*x)
// over coefficient phases x=(message/gamma) modulo one. A complex input is
// mapped componentwise to f(Re(input)) + i*f(Im(input)). The tile size
// controls how many input/frequency columns share each streamed selector
// transform; larger tiles reuse more setup work while retaining more factors.
template <class RefreshSchedule>
inline void CKKSNarrowFunctionalRefreshPrepared(
    typename RefreshSchedule::Output &result,
    const typename RefreshSchedule::Coefficients &prepared,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const std::vector<std::complex<double>> &coefficients,
    bool complex_input = true, std::size_t tile_size = 8,
    const CKKSNarrowRootObserver<typename RefreshSchedule::Narrow,
        RefreshSchedule::input_log_delta> *observer = nullptr)
{
    if (tile_size == 0) throw std::invalid_argument("functional refresh tile size must be positive");
    using P = typename RefreshSchedule::Param;
    using Narrow = typename RefreshSchedule::Narrow;
    if (!complex_input) {
        CKKSNarrowEvalSeries<Narrow>(result, prepared, key.narrow, coefficients, tile_size, observer);
        return;
    }
    std::vector<typename RefreshSchedule::Coefficients> inputs(2);
    inputs[0] = prepared;
    narrow_detail::imaginaryCoefficientHalf<P>(inputs[1], inputs[0]);
    std::vector<typename RefreshSchedule::Output> results;
    CKKSNarrowEvalSeriesBatch<Narrow, RefreshSchedule::input_log_delta>(results, inputs, key.narrow, coefficients, tile_size, observer);
    narrow_detail::combineRealImaginary<P>(result, results[0], results[1]);
}

template <class RefreshSchedule>
inline void CKKSNarrowFunctionalRefresh(
    typename RefreshSchedule::Output &result,
    const typename RefreshSchedule::Input &input,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const std::vector<std::complex<double>> &coefficients,
    bool complex_input = true, std::size_t tile_size = 8,
    const CKKSNarrowRootObserver<typename RefreshSchedule::Narrow,
        RefreshSchedule::input_log_delta> *observer = nullptr)
{
    if (tile_size == 0) throw std::invalid_argument("functional refresh tile size must be positive");
    auto prepared = std::make_unique<typename RefreshSchedule::Coefficients>();
    CKKSNarrowPrepareCoefficients<RefreshSchedule>(*prepared, input, key);
    CKKSNarrowFunctionalRefreshPrepared<RefreshSchedule>(result, *prepared, key,
        coefficients, complex_input, tile_size, observer);
}

template <class RefreshSchedule>
inline void CKKSNarrowFunctionalRefreshPrepared(
    typename RefreshSchedule::Output &result,
    const typename RefreshSchedule::Coefficients &prepared,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const CKKSNarrowTrigonometricApproximation &approximation,
    bool complex_input = true, std::size_t tile_size = 8,
    const CKKSNarrowRootObserver<typename RefreshSchedule::Narrow,
        RefreshSchedule::input_log_delta> *observer = nullptr)
{
    if (approximation.gap != RefreshSchedule::gap() || approximation.message_bound <= 0 ||
        approximation.message_bound >= approximation.gap / 2)
        throw std::invalid_argument("functional approximation does not match refresh gap");
    CKKSNarrowFunctionalRefreshPrepared<RefreshSchedule>(result, prepared, key,
        approximation.coefficients, complex_input, tile_size, observer);
}

template <class RefreshSchedule>
inline void CKKSNarrowFunctionalRefresh(
    typename RefreshSchedule::Output &result,
    const typename RefreshSchedule::Input &input,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const CKKSNarrowTrigonometricApproximation &approximation,
    bool complex_input = true, std::size_t tile_size = 8,
    const CKKSNarrowRootObserver<typename RefreshSchedule::Narrow,
        RefreshSchedule::input_log_delta> *observer = nullptr)
{
    if (approximation.gap != RefreshSchedule::gap() || approximation.message_bound <= 0 ||
        approximation.message_bound >= approximation.gap / 2)
        throw std::invalid_argument("functional approximation does not match refresh gap");
    CKKSNarrowFunctionalRefresh<RefreshSchedule>(result, input, key,
        approximation.coefficients, complex_input, tile_size, observer);
}

// Sections 5.2 and 7: slot-level batching shares W_x across ciphertexts as
// well as frequencies. S2C uses the same public plan and switching keys for
// every input. The tile size bounds retained product-tree factors.
template <class RefreshSchedule, std::uint32_t CoefficientBits = 32>
using CKKSNarrowBasisRefreshSchedule =
    CKKSNarrowFrequencyBasisSchedule<typename RefreshSchedule::Narrow, CoefficientBits>;

// Optional Section7 basis reconstruction after the same S2C/encapsulation
// pipeline. The changed output modulus and coefficient scale are explicit in
// the result type, including when the real and imaginary halves are combined.
template <class RefreshSchedule, std::uint32_t CoefficientBits = 32>
inline void CKKSNarrowFunctionalRefreshWithBasisBatch(
    std::vector<typename CKKSNarrowBasisRefreshSchedule<RefreshSchedule, CoefficientBits>::Output> &results,
    std::span<const typename RefreshSchedule::Input> inputs,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const CKKSNarrowFrequencyBasisKey<CKKSNarrowBasisRefreshSchedule<RefreshSchedule, CoefficientBits>> &basis_key,
    const CKKSNarrowFrequencyBasis &plan,
    const std::vector<std::complex<double>> &coefficients,
    bool complex_input = true, std::size_t tile_size = 8, bool use_fft = true)
{
    using P = typename RefreshSchedule::Param;
    using Basis = CKKSNarrowBasisRefreshSchedule<RefreshSchedule, CoefficientBits>;
    if (tile_size == 0) throw std::invalid_argument("functional basis refresh tile size must be positive");
    const std::size_t halves = complex_input ? 2 : 1;
    if (inputs.size() > std::numeric_limits<std::size_t>::max() / halves)
        throw std::length_error("functional basis refresh batch is too large");
    std::vector<typename RefreshSchedule::Coefficients> prepared(inputs.size() * halves);
    for (std::size_t i = 0; i < inputs.size(); i++) {
        CKKSNarrowPrepareCoefficients<RefreshSchedule>(prepared[halves * i], inputs[i], key);
        if (complex_input)
            narrow_detail::imaginaryCoefficientHalf<P>(prepared[halves * i + 1], prepared[halves * i]);
    }
    std::vector<typename Basis::Output> evaluated;
    CKKSNarrowEvalSeriesWithBasisBatch<Basis, RefreshSchedule::input_log_delta>(
        evaluated, prepared, key.narrow, basis_key, plan, coefficients, tile_size, use_fft);
    results.resize(inputs.size());
    for (std::size_t i = 0; i < inputs.size(); i++) {
        if (complex_input)
            narrow_detail::combineRealImaginary<P>(results[i], evaluated[2 * i], evaluated[2 * i + 1]);
        else results[i] = std::move(evaluated[i]);
    }
}

template <class RefreshSchedule, std::uint32_t CoefficientBits = 32>
inline void CKKSNarrowFunctionalRefreshWithBasis(
    typename CKKSNarrowBasisRefreshSchedule<RefreshSchedule, CoefficientBits>::Output &result,
    const typename RefreshSchedule::Input &input,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const CKKSNarrowFrequencyBasisKey<CKKSNarrowBasisRefreshSchedule<RefreshSchedule, CoefficientBits>> &basis_key,
    const CKKSNarrowFrequencyBasis &plan,
    const std::vector<std::complex<double>> &coefficients,
    bool complex_input = true, std::size_t tile_size = 8, bool use_fft = true)
{
    std::vector<typename CKKSNarrowBasisRefreshSchedule<RefreshSchedule, CoefficientBits>::Output> results;
    CKKSNarrowFunctionalRefreshWithBasisBatch<RefreshSchedule, CoefficientBits>(results, {&input, 1},
        key, basis_key, plan, coefficients, complex_input, tile_size, use_fft);
    result = std::move(results.front());
}

template <class RefreshSchedule, std::uint32_t CoefficientBits = 32>
inline void CKKSNarrowFunctionalRefreshWithBasis(
    typename CKKSNarrowBasisRefreshSchedule<RefreshSchedule, CoefficientBits>::Output &result,
    const typename RefreshSchedule::Input &input,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const CKKSNarrowFrequencyBasisKey<CKKSNarrowBasisRefreshSchedule<RefreshSchedule, CoefficientBits>> &basis_key,
    const CKKSNarrowFrequencyBasis &plan,
    const CKKSNarrowTrigonometricApproximation &approximation,
    bool complex_input = true, std::size_t tile_size = 8, bool use_fft = true)
{
    if (approximation.gap != RefreshSchedule::gap() || approximation.message_bound <= 0 ||
        approximation.message_bound >= approximation.gap / 2)
        throw std::invalid_argument("functional approximation does not match basis refresh gap");
    CKKSNarrowFunctionalRefreshWithBasis<RefreshSchedule, CoefficientBits>(result, input, key,
        basis_key, plan, approximation.coefficients, complex_input, tile_size, use_fft);
}

// Direct-frequency batching retains its original output level and scale.
template <class RefreshSchedule>
inline void CKKSNarrowFunctionalRefreshBatch(
    std::vector<typename RefreshSchedule::Output> &results,
    std::span<const typename RefreshSchedule::Input> inputs,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const std::vector<std::complex<double>> &coefficients,
    bool complex_input = true, std::size_t tile_size = 8)
{
    using P = typename RefreshSchedule::Param;
    using Narrow = typename RefreshSchedule::Narrow;
    if (tile_size == 0) throw std::invalid_argument("functional refresh tile size must be positive");
    const std::size_t halves = complex_input ? 2 : 1;
    if (inputs.size() > std::numeric_limits<std::size_t>::max() / halves)
        throw std::length_error("functional refresh batch is too large");
    std::vector<typename RefreshSchedule::Coefficients> prepared(inputs.size() * halves);
    for (std::size_t i = 0; i < inputs.size(); i++) {
        CKKSNarrowPrepareCoefficients<RefreshSchedule>(prepared[halves * i], inputs[i], key);
        if (complex_input)
            narrow_detail::imaginaryCoefficientHalf<P>(prepared[halves * i + 1], prepared[halves * i]);
    }
    std::vector<typename RefreshSchedule::Output> evaluated;
    CKKSNarrowEvalSeriesBatch<Narrow, RefreshSchedule::input_log_delta>(
        evaluated, prepared, key.narrow, coefficients, tile_size);
    results.resize(inputs.size());
    for (std::size_t i = 0; i < inputs.size(); i++) {
        if (complex_input)
            narrow_detail::combineRealImaginary<P>(results[i], evaluated[2 * i], evaluated[2 * i + 1]);
        else results[i] = std::move(evaluated[i]);
    }
}

template <class RefreshSchedule>
inline void CKKSNarrowFunctionalRefreshBatch(
    std::vector<typename RefreshSchedule::Output> &results,
    std::span<const typename RefreshSchedule::Input> inputs,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    const CKKSNarrowTrigonometricApproximation &approximation,
    bool complex_input = true, std::size_t tile_size = 8)
{
    if (approximation.gap != RefreshSchedule::gap() || approximation.message_bound <= 0 ||
        approximation.message_bound >= approximation.gap / 2)
        throw std::invalid_argument("functional approximation does not match refresh gap");
    CKKSNarrowFunctionalRefreshBatch<RefreshSchedule>(results, inputs, key,
        approximation.coefficients, complex_input, tile_size);
}

// The standard sine refresh is the one-frequency real series with coefficient
// gamma/(4*pi*i). This retains the unary refresh's exact DD arithmetic.
template <class RefreshSchedule>
inline void CKKSNarrowRefreshBatch(
    std::vector<typename RefreshSchedule::Output> &results,
    std::span<const typename RefreshSchedule::Input> inputs,
    const CKKSNarrowRefreshKey<RefreshSchedule> &key,
    std::size_t tile_size = 8)
{
    const std::vector<std::complex<double>> coefficients = {
        {0, 0}, {0, -RefreshSchedule::gap() / (4 * std::acos(-1.0))}};
    CKKSNarrowFunctionalRefreshBatch<RefreshSchedule>(results, inputs, key,
        coefficients, true, tile_size);
}

}  // namespace TFHEpp
#endif
