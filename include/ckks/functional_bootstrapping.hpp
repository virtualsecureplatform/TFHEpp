#pragma once

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "ckks/ckks.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS

namespace TFHEpp {

// First-order trigonometric Hermite interpolation from Theorem 1 / Corollary 1
// of "General Functional Bootstrapping using CKKS" (IACR 2024/1623).
// The input representatives are x = m / p (mod 1), and coefficients[k] is the
// coefficient of exp(2*pi*i*k*x) in Equation (4).
struct CKKSFunctionalBootstrapLUT {
    std::vector<double> values{};
    std::vector<std::complex<double>> coefficients{};

    std::size_t plaintextModulus() const { return values.size(); }
};

inline CKKSFunctionalBootstrapLUT CKKSBuildFunctionalBootstrapLUT(
    const std::vector<double> &values)
{
    if (values.size() < 2)
        throw std::invalid_argument(
            "a functional-bootstrap LUT needs at least two values");
    for (const double value : values)
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "functional-bootstrap LUT values must be finite");

    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const std::size_t p = values.size();
    const long double p_ld = static_cast<long double>(p);

    CKKSFunctionalBootstrapLUT lut;
    lut.values = values;
    lut.coefficients.assign(p, {0.0, 0.0});

    long double mean = 0.0L;
    for (const double value : values) mean += static_cast<long double>(value);
    lut.coefficients[0] = {static_cast<double>(mean / p_ld), 0.0};

    for (std::size_t k = 1; k < p; k++) {
        std::complex<long double> dft{0.0L, 0.0L};
        for (std::size_t l = 0; l < p; l++) {
            const long double angle =
                -2.0L * pi * static_cast<long double>(k * l) / p_ld;
            dft += static_cast<long double>(values[l]) *
                   std::complex<long double>{std::cos(angle),
                                             std::sin(angle)};
        }
        const long double scale =
            2.0L * static_cast<long double>(p - k) / (p_ld * p_ld);
        lut.coefficients[k] = {
            static_cast<double>(scale * dft.real()),
            static_cast<double>(scale * dft.imag())};
    }
    return lut;
}

inline std::complex<double> CKKSEvaluateFunctionalBootstrapLUTComplex(
    const CKKSFunctionalBootstrapLUT &lut, double x)
{
    if (lut.values.size() < 2 ||
        lut.coefficients.size() != lut.values.size())
        throw std::invalid_argument("invalid functional-bootstrap LUT");
    constexpr double two_pi =
        6.283185307179586476925286766559005768;
    // Reducing x before constructing the unit-circle point improves the
    // reference evaluator for the integer overflow representatives produced by
    // modulus raising.
    const double reduced = std::remainder(x, 1.0);
    const std::complex<double> z =
        std::polar(1.0, two_pi * reduced);
    std::complex<double> value = lut.coefficients.back();
    for (std::size_t k = lut.coefficients.size() - 1; k-- > 0;)
        value = value * z + lut.coefficients[k];
    return value;
}

inline double CKKSEvaluateFunctionalBootstrapLUT(
    const CKKSFunctionalBootstrapLUT &lut, double x)
{
    return CKKSEvaluateFunctionalBootstrapLUTComplex(lut, x).real();
}

inline double CKKSEvaluateFunctionalBootstrapLUTDerivative(
    const CKKSFunctionalBootstrapLUT &lut, double x,
    std::uint32_t derivative_order = 1)
{
    if (derivative_order == 0)
        return CKKSEvaluateFunctionalBootstrapLUT(lut, x);
    if (lut.values.size() < 2 ||
        lut.coefficients.size() != lut.values.size())
        throw std::invalid_argument("invalid functional-bootstrap LUT");

    constexpr double two_pi =
        6.283185307179586476925286766559005768;
    const double reduced = std::remainder(x, 1.0);
    const std::complex<double> z = std::polar(1.0, two_pi * reduced);
    std::complex<double> z_to_k{1.0, 0.0};
    std::complex<double> derivative{0.0, 0.0};
    for (std::size_t k = 1; k < lut.coefficients.size(); k++) {
        z_to_k *= z;
        std::complex<double> factor{1.0, 0.0};
        const std::complex<double> ik{0.0, two_pi * static_cast<double>(k)};
        for (std::uint32_t j = 0; j < derivative_order; j++) factor *= ik;
        derivative += lut.coefficients[k] * factor * z_to_k;
    }
    return derivative.real();
}

// A real Chebyshev approximation of the exact periodic interpolation.  This is
// a TFHEpp-native EvalLUT representation: homomorphic evaluation reuses the
// existing balanced Chebyshev basis and relinearization-key chain.  The input
// to EvalLUT is t = x / input_bound, so it must lie in [-1, 1].
struct CKKSFunctionalBootstrapChebyshevPolynomial {
    std::size_t plaintext_modulus = 0;
    std::size_t degree = 0;
    double input_bound = 1.0;
    std::vector<double> coefficients{};
};

inline CKKSFunctionalBootstrapChebyshevPolynomial
CKKSBuildFunctionalBootstrapChebyshevPolynomial(
    const CKKSFunctionalBootstrapLUT &lut, std::size_t degree,
    double input_bound = 1.0)
{
    if (lut.values.size() < 2 ||
        lut.coefficients.size() != lut.values.size())
        throw std::invalid_argument("invalid functional-bootstrap LUT");
    if (degree < 2)
        throw std::invalid_argument(
            "functional-bootstrap Chebyshev degree must be at least two");
    if (!(input_bound > 0.0) || !std::isfinite(input_bound))
        throw std::invalid_argument(
            "functional-bootstrap input bound must be positive and finite");

    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const std::size_t nodes = degree + 1;
    std::vector<long double> samples(nodes);
    for (std::size_t j = 0; j < nodes; j++) {
        const long double theta =
            pi * (static_cast<long double>(j) + 0.5L) /
            static_cast<long double>(nodes);
        const double t = static_cast<double>(std::cos(theta));
        samples[j] = static_cast<long double>(
            CKKSEvaluateFunctionalBootstrapLUT(lut, input_bound * t));
    }

    CKKSFunctionalBootstrapChebyshevPolynomial polynomial;
    polynomial.plaintext_modulus = lut.plaintextModulus();
    polynomial.degree = degree;
    polynomial.input_bound = input_bound;
    polynomial.coefficients.resize(nodes);
    for (std::size_t k = 0; k < nodes; k++) {
        long double coefficient = 0.0L;
        for (std::size_t j = 0; j < nodes; j++) {
            const long double theta =
                pi * static_cast<long double>(k) *
                (static_cast<long double>(j) + 0.5L) /
                static_cast<long double>(nodes);
            coefficient += samples[j] * std::cos(theta);
        }
        coefficient *= 2.0L / static_cast<long double>(nodes);
        if (k == 0) coefficient *= 0.5L;
        polynomial.coefficients[k] = static_cast<double>(coefficient);
    }
    return polynomial;
}

inline double CKKSEvaluateFunctionalBootstrapChebyshevPolynomialNormalized(
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial,
    double normalized_input)
{
    if (polynomial.coefficients.size() != polynomial.degree + 1)
        throw std::invalid_argument(
            "invalid functional-bootstrap Chebyshev polynomial");
    return CKKSEvaluateChebyshevUnit(polynomial.coefficients,
                                     normalized_input);
}

inline double CKKSEvaluateFunctionalBootstrapChebyshevPolynomial(
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial, double x)
{
    return CKKSEvaluateFunctionalBootstrapChebyshevPolynomialNormalized(
        polynomial, x / polynomial.input_bound);
}

// FHE-friendly complex exponential from Section 3.2.  The Chebyshev stage
// approximates exp(2*pi*i*input_bound*t/2^r) on t in [-1,1], after which r
// squarings recover E(x) = exp(2*pi*i*x).  Evaluating Equation (4) in powers of
// E(x) avoids directly approximating every oscillation of R(x).
struct CKKSFunctionalBootstrapComplexExponentialPolynomial {
    std::size_t degree = 0;
    std::uint32_t double_angle = 0;
    double input_bound = 1.0;
    std::vector<std::complex<double>> coefficients{};
};

inline CKKSFunctionalBootstrapComplexExponentialPolynomial
CKKSBuildFunctionalBootstrapComplexExponentialPolynomial(
    std::size_t degree, std::uint32_t double_angle,
    double input_bound = 1.0)
{
    if (degree < 2)
        throw std::invalid_argument(
            "functional-bootstrap exponential degree must be at least two");
    if (!(input_bound > 0.0) || !std::isfinite(input_bound))
        throw std::invalid_argument(
            "functional-bootstrap exponential input bound must be positive");

    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const std::size_t nodes = degree + 1;
    const long double divisor =
        std::ldexp(1.0L, static_cast<int>(double_angle));
    std::vector<std::complex<long double>> samples(nodes);
    for (std::size_t j = 0; j < nodes; j++) {
        const long double theta =
            pi * (static_cast<long double>(j) + 0.5L) /
            static_cast<long double>(nodes);
        const long double t = std::cos(theta);
        const long double phase =
            2.0L * pi * static_cast<long double>(input_bound) * t / divisor;
        samples[j] = {std::cos(phase), std::sin(phase)};
    }

    CKKSFunctionalBootstrapComplexExponentialPolynomial polynomial;
    polynomial.degree = degree;
    polynomial.double_angle = double_angle;
    polynomial.input_bound = input_bound;
    polynomial.coefficients.resize(nodes);
    for (std::size_t k = 0; k < nodes; k++) {
        std::complex<long double> coefficient{0.0L, 0.0L};
        for (std::size_t j = 0; j < nodes; j++) {
            const long double theta =
                pi * static_cast<long double>(k) *
                (static_cast<long double>(j) + 0.5L) /
                static_cast<long double>(nodes);
            coefficient += samples[j] * std::cos(theta);
        }
        coefficient *= 2.0L / static_cast<long double>(nodes);
        if (k == 0) coefficient *= 0.5L;
        polynomial.coefficients[k] =
            static_cast<std::complex<double>>(coefficient);
    }
    return polynomial;
}

inline std::complex<double>
CKKSEvaluateFunctionalBootstrapComplexExponentialNormalized(
    const CKKSFunctionalBootstrapComplexExponentialPolynomial &polynomial,
    double normalized_input)
{
    if (polynomial.coefficients.size() != polynomial.degree + 1)
        throw std::invalid_argument(
            "invalid functional-bootstrap exponential polynomial");
    std::complex<double> b_next{0.0, 0.0};
    std::complex<double> b_next_next{0.0, 0.0};
    for (std::size_t i = polynomial.coefficients.size() - 1; i > 0; i--) {
        const auto b = 2.0 * normalized_input * b_next - b_next_next +
                       polynomial.coefficients[i];
        b_next_next = b_next;
        b_next = b;
    }
    std::complex<double> value = polynomial.coefficients[0] +
                                 normalized_input * b_next - b_next_next;
    for (std::uint32_t i = 0; i < polynomial.double_angle; i++) value *= value;
    return value;
}

inline double CKKSEvaluateFunctionalBootstrapLUTFHEFriendly(
    const CKKSFunctionalBootstrapLUT &lut,
    const CKKSFunctionalBootstrapComplexExponentialPolynomial &exponential,
    double x)
{
    const std::complex<double> z =
        CKKSEvaluateFunctionalBootstrapComplexExponentialNormalized(
            exponential, x / exponential.input_bound);
    std::complex<double> value = lut.coefficients.back();
    for (std::size_t k = lut.coefficients.size() - 1; k-- > 0;)
        value = value * z + lut.coefficients[k];
    return value.real();
}

namespace ckks_detail {

template <std::size_t I, class Traits, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t CoeffLogDelta,
          std::size_t Degree>
inline void CKKSAddComplexPolynomialTermsImpl(
    typename Traits::Ciphertext &res, const typename Traits::PowerBasis &basis,
    const std::vector<std::complex<double>> &coefficients)
{
    if constexpr (I <= Degree) {
        if (I < coefficients.size() && coefficients[I] !=
                                           std::complex<double>{0.0, 0.0}) {
            using PowerPtr =
                std::tuple_element_t<I - 1, typename Traits::PowerBasis>;
            using PowerCt = typename PowerPtr::element_type;
            const auto &power = std::get<I - 1>(basis);
            assert(power != nullptr);
            auto reduced = std::make_unique<
                CKKSCiphertext<P, Traits::term_input_log_q, LogDelta>>();
            CKKSLevelReduce<P, PowerCt::log_q, Traits::term_input_log_q,
                            LogDelta>(*reduced, *power);
            auto term = std::make_unique<CKKSPlainMulResult<
                P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>>();
            CKKSPlainMulByConstant<P, Traits::term_input_log_q, LogDelta,
                                   CoeffLogDelta>(*term, *reduced,
                                                  coefficients[I]);
            CKKSAddInPlace<P, Traits::log_q, LogDelta>(res, *term);
        }
        CKKSAddComplexPolynomialTermsImpl<
            I + 1, Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
            res, basis, coefficients);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t DoubleAngle,
          class RelinKeyProvider>
inline void CKKSFunctionalBootstrapDoubleAngleImpl(
    CKKSCiphertext<P, StartLogQ - DoubleAngle * LogDelta, LogDelta> &res,
    const CKKSCiphertext<P, StartLogQ - I * LogDelta, LogDelta> &ct,
    const RelinKeyProvider &keys)
{
    if constexpr (I == DoubleAngle) {
        res = ct;
    }
    else {
        using Next = CKKSMultResult<
            P, StartLogQ - I * LogDelta, LogDelta,
            StartLogQ - I * LogDelta, LogDelta>;
        auto next = std::make_unique<Next>();
        CKKSSquare<P>(*next, ct, keys.template get<I>());
        maybe_release_key<I>(keys);
        CKKSFunctionalBootstrapDoubleAngleImpl<
            I + 1, P, StartLogQ, LogDelta, DoubleAngle>(res, *next, keys);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          class RelinKeyProvider>
inline void CKKSEvalComplexChebyshevPolynomialWithKeyProvider(
    CKKSPowerPolynomialResult<P, StartLogQ, LogDelta, CoeffLogDelta, Degree>
        &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const std::vector<std::complex<double>> &coefficients,
    const RelinKeyProvider &keys)
{
    using Traits = CKKSPowerPolynomialEvaluatorTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>;
    assert(coefficients.size() <= Degree + 1);
    typename Traits::PowerBasis basis;
    ckks_detail::CKKSBuildChebyshevBasisImpl<
        1, Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(basis, ct,
                                                                  keys);
    CKKSSetTransparentConstant<P, Traits::log_q, LogDelta>(
        res, coefficients.empty() ? std::complex<double>{0.0, 0.0}
                                  : coefficients[0]);
    ckks_detail::CKKSAddComplexPolynomialTermsImpl<
        1, Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
        res, basis, coefficients);
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          class RelinKeyProvider>
inline void CKKSEvalComplexPowerPolynomialWithKeyProvider(
    CKKSPowerPolynomialResult<P, StartLogQ, LogDelta, CoeffLogDelta, Degree>
        &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const std::vector<std::complex<double>> &coefficients,
    const RelinKeyProvider &keys)
{
    using Traits = CKKSPowerPolynomialEvaluatorTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>;
    assert(coefficients.size() <= Degree + 1);
    typename Traits::PowerBasis basis;
    ckks_detail::CKKSBuildPowerBasisImpl<
        1, Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(basis, ct,
                                                                  keys);
    CKKSSetTransparentConstant<P, Traits::log_q, LogDelta>(
        res, coefficients.empty() ? std::complex<double>{0.0, 0.0}
                                  : coefficients[0]);
    ckks_detail::CKKSAddComplexPolynomialTermsImpl<
        1, Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
        res, basis, coefficients);
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, class PolynomialRelinKeyProvider,
          class DoubleAngleRelinKeyProvider>
inline void CKKSFunctionalBootstrapComplexExponential(
    CKKSCiphertext<
        P, CKKSPowerPolynomialEvaluatorTraits<
                   P, StartLogQ, LogDelta, CoeffLogDelta, Degree>::log_q -
               DoubleAngle * LogDelta,
        LogDelta> &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &normalized_input,
    const CKKSFunctionalBootstrapComplexExponentialPolynomial &polynomial,
    const PolynomialRelinKeyProvider &polynomial_keys,
    const DoubleAngleRelinKeyProvider &double_angle_keys)
{
    using PolynomialTraits = CKKSPowerPolynomialEvaluatorTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>;
    assert(polynomial.degree <= Degree);
    assert(polynomial.double_angle == DoubleAngle);
    auto base = std::make_unique<typename PolynomialTraits::Ciphertext>();
    CKKSEvalComplexChebyshevPolynomialWithKeyProvider<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
        *base, normalized_input, polynomial.coefficients, polynomial_keys);
    ckks_detail::CKKSFunctionalBootstrapDoubleAngleImpl<
        0, P, PolynomialTraits::log_q, LogDelta, DoubleAngle>(
        res, *base, double_angle_keys);
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
using CKKSFunctionalBootstrapEvalLUTTraits =
    CKKSPowerPolynomialEvaluatorTraits<P, StartLogQ, LogDelta,
                                       CoeffLogDelta, Degree>;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
using CKKSFunctionalBootstrapEvalLUTResult =
    typename CKKSFunctionalBootstrapEvalLUTTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>::Ciphertext;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
using CKKSFunctionalBootstrapEvalLUTRelinKey =
    typename CKKSFunctionalBootstrapEvalLUTTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>::RelinKeyChain;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
inline void CKKSFunctionalBootstrapEvalLUTKeyGen(
    CKKSFunctionalBootstrapEvalLUTRelinKey<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree> &keys,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    using Traits = CKKSFunctionalBootstrapEvalLUTTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>;
    CKKSRelinKeyChainGen<P, StartLogQ, LogDelta, Traits::relin_depth>(
        keys, key, noise);
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          class RelinKeyProvider>
inline void CKKSFunctionalBootstrapEvalLUTNormalizedWithKeyProvider(
    CKKSFunctionalBootstrapEvalLUTResult<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree> &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &normalized_input,
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial,
    const RelinKeyProvider &keys)
{
    assert(polynomial.degree <= Degree);
    assert(polynomial.coefficients.size() <= Degree + 1);
    CKKSEvalChebyshevPolynomialWithKeyProvider<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
        res, normalized_input, polynomial.coefficients, keys);
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
inline void CKKSFunctionalBootstrapEvalLUTNormalized(
    CKKSFunctionalBootstrapEvalLUTResult<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree> &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &normalized_input,
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial,
    const CKKSFunctionalBootstrapEvalLUTRelinKey<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree> &keys)
{
    CKKSFunctionalBootstrapEvalLUTNormalizedWithKeyProvider<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
        res, normalized_input, polynomial, keys);
}

// Build the direct bounded EvalLUT used by the dense functional-bootstrap
// wrapper below.  Dense C2S emits
//
//   t = (m / p + message_ratio * overflow) /
//       (evalmod_k * message_ratio).
//
// The periodic LUT is therefore sampled at x = input_bound * t, while the
// result is divided by message_ratio to compensate for the scaling applied by
// the existing dense StC plan.
template <class Schedule>
inline CKKSFunctionalBootstrapChebyshevPolynomial
CKKSBuildDenseFunctionalBootstrapChebyshevPolynomial(
    const CKKSFunctionalBootstrapLUT &lut)
{
    std::vector<double> scaled_values = lut.values;
    for (double &value : scaled_values) value /= Schedule::message_ratio;
    const auto scaled_lut =
        CKKSBuildFunctionalBootstrapLUT(scaled_values);
    return CKKSBuildFunctionalBootstrapChebyshevPolynomial(
        scaled_lut, Schedule::evalmod_degree,
        static_cast<double>(Schedule::evalmod_k) * Schedule::message_ratio);
}

namespace ckks_detail {

template <class KeyProvider, bool Retain>
struct CKKSDenseFunctionalBootstrapPolynomialKeyChain {
    const KeyProvider &provider;

    template <std::size_t I>
    decltype(auto) get() const
    {
        return provider.template polynomial_relin<I>();
    }

    template <std::size_t I>
    void release() const
    {
        if constexpr (!Retain && requires {
                          provider.template release_polynomial_relin<I>();
                      }) {
            provider.template release_polynomial_relin<I>();
        }
    }
};

template <class Schedule>
using CKKSDenseFunctionalBootstrapEvalTraits =
    CKKSFunctionalBootstrapEvalLUTTraits<
        typename Schedule::Param, Schedule::after_component_split_log_q,
        Schedule::log_delta, Schedule::evalmod_log_scale,
        Schedule::evalmod_degree>;

template <class Schedule, class KeyProvider>
inline void CKKSDenseFunctionalBootstrapWithKeyProviderImpl(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial,
    const KeyProvider &key_provider, CKKSDenseBootstrapTimings *timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    using P = typename Schedule::Param;
    using EvalTraits = CKKSDenseFunctionalBootstrapEvalTraits<Schedule>;
    static_assert(EvalTraits::log_q >= Schedule::after_evalmod_log_q,
                  "dense EvalLUT consumes more levels than the schedule");
    assert(polynomial.degree == Schedule::evalmod_degree);
    assert(polynomial.plaintext_modulus >= 2);

    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan =
        key_provider.linear_plan();

    auto raised = std::make_unique<typename Schedule::BootstrapCiphertext>();
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->modraise_ms,
        [&] {
            CKKSModRaiseBoundedPhaseRandomized<
                P, Schedule::input_log_q, Schedule::boot_log_q,
                Schedule::log_delta, Schedule::modraise_mask_bound>(*raised,
                                                                     ct);
        },
        progress, "functional_modraise");

    auto coeff_to_slot =
        std::make_unique<typename Schedule::CoeffToSlotCiphertext>();
    const CKKSDenseBootstrapLinearKeyProviderChain<KeyProvider, true>
        coeff_to_slot_galois{key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->coeff_to_slot_ms,
        [&] {
            CKKSDenseBootstrapCoeffToSlotStagesBSGS<Schedule>(
                *coeff_to_slot, *raised, linear_plan, coeff_to_slot_galois);
        },
        progress, "functional_coeff_to_slot");
    raised.reset();

    auto real_component =
        std::make_unique<typename Schedule::ComponentCiphertext>();
    auto imag_component =
        std::make_unique<typename Schedule::ComponentCiphertext>();
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->split_ms,
        [&] {
            CKKSExtractRealImagSlots<
                P, Schedule::after_coeff_to_slot_log_q, Schedule::log_delta,
                Schedule::component_split_plain_log_delta>(
                *real_component, *imag_component, *coeff_to_slot,
                key_provider.packed_conjugate_galois());
        },
        progress, "functional_split");
    coeff_to_slot.reset();
    if constexpr (requires { key_provider.release_packed_conjugate_galois(); })
        key_provider.release_packed_conjugate_galois();

    auto real_lut = std::make_unique<typename EvalTraits::Ciphertext>();
    const CKKSDenseFunctionalBootstrapPolynomialKeyChain<KeyProvider, true>
        retained_polynomial_keys{key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->real_evalmod_ms,
        [&] {
            CKKSFunctionalBootstrapEvalLUTNormalizedWithKeyProvider<
                P, Schedule::after_component_split_log_q,
                Schedule::log_delta, Schedule::evalmod_log_scale,
                Schedule::evalmod_degree>(*real_lut, *real_component,
                                           polynomial,
                                           retained_polynomial_keys);
        },
        progress, "functional_real_eval_lut");
    real_component.reset();

    auto imag_lut = std::make_unique<typename EvalTraits::Ciphertext>();
    const CKKSDenseFunctionalBootstrapPolynomialKeyChain<KeyProvider, false>
        polynomial_keys{key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->imag_evalmod_ms,
        [&] {
            CKKSFunctionalBootstrapEvalLUTNormalizedWithKeyProvider<
                P, Schedule::after_component_split_log_q,
                Schedule::log_delta, Schedule::evalmod_log_scale,
                Schedule::evalmod_degree>(*imag_lut, *imag_component,
                                           polynomial, polynomial_keys);
        },
        progress, "functional_imag_eval_lut");
    imag_component.reset();

    auto real_stc_input =
        std::make_unique<typename Schedule::EvalModCiphertext>();
    auto imag_stc_input =
        std::make_unique<typename Schedule::EvalModCiphertext>();
    CKKSLevelReduce<P, EvalTraits::log_q, Schedule::after_evalmod_log_q,
                    Schedule::log_delta>(*real_stc_input, *real_lut);
    CKKSLevelReduce<P, EvalTraits::log_q, Schedule::after_evalmod_log_q,
                    Schedule::log_delta>(*imag_stc_input, *imag_lut);
    real_lut.reset();
    imag_lut.reset();

    const CKKSDenseBootstrapLinearKeyProviderChain<KeyProvider, false>
        slot_to_coeff_galois{key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->slot_to_coeff_ms,
        [&] {
            CKKSDenseBootstrapSlotToCoeffStagesBSGSDualInputSharedTail<
                Schedule>(res, *real_stc_input, *imag_stc_input, linear_plan,
                          slot_to_coeff_galois);
        },
        progress, "functional_slot_to_coeff");
}

}  // namespace ckks_detail

template <class Schedule, class KeyProvider>
inline void CKKSDenseFunctionalBootstrapWithKeyProvider(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial,
    const KeyProvider &key_provider)
{
    ckks_detail::CKKSDenseFunctionalBootstrapWithKeyProviderImpl<Schedule>(
        res, ct, polynomial, key_provider, nullptr);
}

template <class Schedule, class KeyProvider>
inline void CKKSDenseFunctionalBootstrapWithKeyProviderTimed(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial,
    const KeyProvider &key_provider, CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    timings = {};
    ckks_detail::CKKSDenseFunctionalBootstrapWithKeyProviderImpl<Schedule>(
        res, ct, polynomial, key_provider, &timings, progress);
}

template <class Schedule>
inline void CKKSDenseFunctionalBootstrap(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial,
    const CKKSDenseBootstrapKey<Schedule> &bootstrap_key)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<Schedule> provider(
        bootstrap_key);
    CKKSDenseFunctionalBootstrapWithKeyProvider<Schedule>(
        res, ct, polynomial, provider);
}

template <class Schedule>
inline void CKKSDenseFunctionalBootstrapTimed(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapChebyshevPolynomial &polynomial,
    const CKKSDenseBootstrapKey<Schedule> &bootstrap_key,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<Schedule> provider(
        bootstrap_key);
    CKKSDenseFunctionalBootstrapWithKeyProviderTimed<Schedule>(
        res, ct, polynomial, provider, timings, progress);
}

// Schedule for Algorithm 1's FHE-friendly EvalLUT.  The expensive functional
// stage may deliberately use a smaller scale than the surrounding dense CKKS
// transforms.  This is the key parameter-selection trick that lets the p <= 8
// circuit fit under lvl6's already security-checked 896-bit largest modulus.
template <class BaseSchedule, std::uint32_t EvalLogDelta,
          std::uint32_t ExponentialCoeffLogDelta,
          std::size_t ExponentialDegree, std::uint32_t DoubleAngle,
          std::uint32_t LUTCoeffLogDelta, std::size_t LUTDegree>
struct CKKSDenseFHEFriendlyFunctionalBootstrapSchedule : BaseSchedule {
    using Param = typename BaseSchedule::Param;
    static_assert(EvalLogDelta < BaseSchedule::log_delta);
    static_assert(ExponentialDegree > 1);
    static_assert(LUTDegree > 1);

    static constexpr std::uint32_t eval_log_delta = EvalLogDelta;
    static constexpr std::uint32_t eval_scale_down_bits =
        BaseSchedule::log_delta - EvalLogDelta;
    static constexpr std::uint32_t exponential_coeff_log_delta =
        ExponentialCoeffLogDelta;
    static constexpr std::size_t exponential_degree = ExponentialDegree;
    static constexpr std::uint32_t functional_double_angle = DoubleAngle;
    static constexpr std::uint32_t lut_coeff_log_delta = LUTCoeffLogDelta;
    static constexpr std::size_t lut_degree = LUTDegree;
    static constexpr std::uint32_t functional_start_log_q =
        BaseSchedule::after_component_split_log_q - eval_scale_down_bits;
    static constexpr double functional_input_bound =
        static_cast<double>(BaseSchedule::evalmod_k) *
        BaseSchedule::message_ratio;

    using ExponentialTraits = CKKSPowerPolynomialEvaluatorTraits<
        Param, functional_start_log_q, EvalLogDelta,
        ExponentialCoeffLogDelta, ExponentialDegree>;
    static constexpr std::uint32_t after_exponential_log_q =
        ExponentialTraits::log_q - DoubleAngle * EvalLogDelta;
    using LUTTraits = CKKSPowerPolynomialEvaluatorTraits<
        Param, after_exponential_log_q, EvalLogDelta, LUTCoeffLogDelta,
        LUTDegree>;

    static constexpr std::uint32_t after_evalmod_log_q = LUTTraits::log_q;
    static constexpr std::uint32_t output_log_q =
        after_evalmod_log_q -
        BaseSchedule::slot_to_coeff_level_count *
            BaseSchedule::slot_to_coeff_plain_log_delta;
    static_assert(ExponentialTraits::log_q >
                  DoubleAngle * EvalLogDelta);
    static_assert(after_evalmod_log_q > BaseSchedule::log_delta);
    static_assert(output_log_q > BaseSchedule::input_log_q);

    using EvalModCiphertext =
        CKKSCiphertext<Param, after_evalmod_log_q, BaseSchedule::log_delta>;
    using OutputCiphertext = CKKSStagedPlainMulResult<
        Param, after_evalmod_log_q, BaseSchedule::log_delta,
        BaseSchedule::slot_to_coeff_plain_log_delta,
        BaseSchedule::slot_to_coeff_level_count>;
};

// Functional bootstrapping only needs the minimum legal message ratio 2: the
// LUT is periodic modulo one.  Spending that input headroom on a 51-bit
// EvalLUT scale keeps TFHEpp's torus/FFT approximation error below the
// double-angle stability threshold without increasing the 896-bit top level.
struct lvl6CKKSDenseFunctionalBootstrapBaseSchedule
    : CKKSDenseBootstrapSchedule<lvl6param, 52, 1, 896, 44, 7, 63, 18, 2, 7,
                                  24, 128, 0, 44, 44, 14, 7, 7, 2, 1> {
    template <std::size_t I>
    static consteval int coeff_to_slot_bsgs_step()
    {
        return I == 0 ? 2048 : 16;
    }

    template <std::size_t I>
    static consteval int slot_to_coeff_bsgs_step()
    {
        return I + 1 == slot_to_coeff_level_count ? 1024 : 64;
    }
};

using lvl6CKKSDenseFunctionalBootstrapP8Schedule =
    CKKSDenseFHEFriendlyFunctionalBootstrapSchedule<
        lvl6CKKSDenseFunctionalBootstrapBaseSchedule, 51, 34, 58, 3, 34, 7>;

namespace ckks_detail {

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta, class Seq>
struct CKKSFunctionalSeededRelinKeyTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t... Is>
struct CKKSFunctionalSeededRelinKeyTuple<
    P, StartLogQ, LogDelta, std::index_sequence<Is...>> {
    using type = std::tuple<
        CKKSSeededRelinKey<P, StartLogQ - (Is + 1) * LogDelta>...>;
};

template <std::size_t I, class Chain, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::size_t Depth>
inline void CKKSFunctionalSeededRelinKeyChainGenImpl(
    Chain &chain, const Key<P> &key, CKKSNoise noise)
{
    if constexpr (I < Depth) {
        constexpr std::uint32_t log_q =
            StartLogQ - (I + 1) * LogDelta;
        chain.template get<I>() =
            *makeCKKSSeededRelinKey<P, log_q>(key, noise);
        CKKSFunctionalSeededRelinKeyChainGenImpl<
            I + 1, Chain, P, StartLogQ, LogDelta, Depth>(chain, key, noise);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t Depth>
struct CKKSFunctionalSeededRelinKeyChain {
    using Tuple = typename ckks_detail::CKKSFunctionalSeededRelinKeyTuple<
        P, StartLogQ, LogDelta,
        std::make_index_sequence<Depth>>::type;
    Tuple keys{};

    template <std::size_t I>
    auto &get()
    {
        return std::get<I>(keys);
    }
    template <std::size_t I>
    const auto &get() const
    {
        return std::get<I>(keys);
    }
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t Depth>
inline void CKKSFunctionalSeededRelinKeyChainGen(
    CKKSFunctionalSeededRelinKeyChain<P, StartLogQ, LogDelta, Depth> &chain,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    ckks_detail::CKKSFunctionalSeededRelinKeyChainGenImpl<
        0, decltype(chain), P, StartLogQ, LogDelta, Depth>(chain, key, noise);
}

template <class Schedule>
struct CKKSDenseFHEFriendlyFunctionalBootstrapKey {
    using P = typename Schedule::Param;
    using CoeffToSlotGaloisKeyChain = CKKSSparseGaloisKeyChain<
        P, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>;
    using SlotToCoeffGaloisKeyChain = CKKSSparseGaloisKeyChain<
        P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>;
    using ExponentialRelinKeyChain =
        typename Schedule::ExponentialTraits::RelinKeyChain;
    using DoubleAngleRelinKeyChain = CKKSRelinKeyChain<
        P, Schedule::ExponentialTraits::log_q, Schedule::eval_log_delta,
        Schedule::functional_double_angle>;
    using LUTRelinKeyChain = typename Schedule::LUTTraits::RelinKeyChain;

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan{};
    CoeffToSlotGaloisKeyChain coeff_to_slot_galois{};
    CKKSSparseGaloisKey<P, Schedule::after_coeff_to_slot_log_q>
        packed_conjugate_galois{};
    ExponentialRelinKeyChain exponential_relin{};
    DoubleAngleRelinKeyChain double_angle_relin{};
    LUTRelinKeyChain lut_relin{};
    CKKSSparseGaloisKey<P, Schedule::after_evalmod_log_q>
        output_conjugate_galois{};
    SlotToCoeffGaloisKeyChain slot_to_coeff_galois{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(linear_plan, coeff_to_slot_galois,
                packed_conjugate_galois, exponential_relin,
                double_angle_relin, lut_relin, output_conjugate_galois,
                slot_to_coeff_galois);
    }
};

template <class Schedule>
inline void CKKSDenseFHEFriendlyFunctionalBootstrapKeyGen(
    CKKSDenseFHEFriendlyFunctionalBootstrapKey<Schedule> &bootstrap_key,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(bootstrap_key.linear_plan);
    CKKSDenseBootstrapRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapRotationKeyUsage<Schedule>(
        rotation_usage, bootstrap_key.linear_plan);
    CKKSDenseBootstrapCoeffToSlotGaloisKeyChainGen<Schedule>(
        bootstrap_key.coeff_to_slot_galois, rotation_usage, key, noise);
    CKKSDenseBootstrapPackedConjugateGaloisKeyGen<Schedule>(
        bootstrap_key.packed_conjugate_galois, rotation_usage, key, noise);
    CKKSRelinKeyChainGen<P, Schedule::functional_start_log_q,
                         Schedule::eval_log_delta,
                         Schedule::ExponentialTraits::relin_depth>(
        bootstrap_key.exponential_relin, key, noise);
    CKKSRelinKeyChainGen<P, Schedule::ExponentialTraits::log_q,
                         Schedule::eval_log_delta,
                         Schedule::functional_double_angle>(
        bootstrap_key.double_angle_relin, key, noise);
    CKKSRelinKeyChainGen<P, Schedule::after_exponential_log_q,
                         Schedule::eval_log_delta,
                         Schedule::LUTTraits::relin_depth>(
        bootstrap_key.lut_relin, key, noise);
    CKKSRotationKeyIndexSet<P> conjugation_usage{};
    CKKSMarkConjugationKeyIndex<P>(conjugation_usage);
    CKKSSparseGaloisKeyGen<P, Schedule::after_evalmod_log_q>(
        bootstrap_key.output_conjugate_galois, key, conjugation_usage, noise);
    CKKSDenseBootstrapSlotToCoeffGaloisKeyChainGen<Schedule>(
        bootstrap_key.slot_to_coeff_galois, rotation_usage, key, noise);
}

// Practical in-memory lvl6 key: hybrid linear-transform keys avoid the dense
// rotation-key footprint, while seeded relinearization keys are expanded only
// while their functional phase is active.
template <class Schedule>
struct CKKSDenseFHEFriendlyFunctionalBootstrapHybridGiantSeededKey {
    using P = typename Schedule::Param;
    using CoeffToSlotGaloisKeyChain = CKKSHybridSparseGaloisKeyChain<
        P, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>;
    using SlotToCoeffGaloisKeyChain = CKKSHybridSparseGaloisKeyChain<
        P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>;
    using ExponentialRelinKeyChain = CKKSFunctionalSeededRelinKeyChain<
        P, Schedule::functional_start_log_q, Schedule::eval_log_delta,
        Schedule::ExponentialTraits::relin_depth>;
    using DoubleAngleRelinKeyChain = CKKSFunctionalSeededRelinKeyChain<
        P, Schedule::ExponentialTraits::log_q, Schedule::eval_log_delta,
        Schedule::functional_double_angle>;
    using LUTRelinKeyChain = CKKSFunctionalSeededRelinKeyChain<
        P, Schedule::after_exponential_log_q, Schedule::eval_log_delta,
        Schedule::LUTTraits::relin_depth>;

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan{};
    mutable std::unique_ptr<CoeffToSlotGaloisKeyChain>
        coeff_to_slot_galois =
            std::make_unique<CoeffToSlotGaloisKeyChain>();
    mutable std::unique_ptr<
        CKKSSparseGaloisKey<P, Schedule::after_coeff_to_slot_log_q>>
        packed_conjugate_galois = std::make_unique<
            CKKSSparseGaloisKey<P, Schedule::after_coeff_to_slot_log_q>>();
    mutable std::unique_ptr<ExponentialRelinKeyChain> exponential_relin =
        std::make_unique<ExponentialRelinKeyChain>();
    mutable std::unique_ptr<DoubleAngleRelinKeyChain> double_angle_relin =
        std::make_unique<DoubleAngleRelinKeyChain>();
    mutable std::unique_ptr<LUTRelinKeyChain> lut_relin =
        std::make_unique<LUTRelinKeyChain>();
    mutable std::unique_ptr<
        CKKSSparseGaloisKey<P, Schedule::after_evalmod_log_q>>
        output_conjugate_galois = std::make_unique<
            CKKSSparseGaloisKey<P, Schedule::after_evalmod_log_q>>();
    mutable std::unique_ptr<SlotToCoeffGaloisKeyChain>
        slot_to_coeff_galois =
            std::make_unique<SlotToCoeffGaloisKeyChain>();

    void release_coeff_to_slot() const { coeff_to_slot_galois.reset(); }
    void release_packed_conjugate() const
    {
        packed_conjugate_galois.reset();
    }
    void release_eval_lut() const
    {
        exponential_relin.reset();
        double_angle_relin.reset();
        lut_relin.reset();
        output_conjugate_galois.reset();
    }
    void release_slot_to_coeff() const { slot_to_coeff_galois.reset(); }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(linear_plan, *coeff_to_slot_galois,
                *packed_conjugate_galois, *exponential_relin,
                *double_angle_relin, *lut_relin,
                *output_conjugate_galois, *slot_to_coeff_galois);
    }
};

template <class Schedule>
inline void CKKSDenseFHEFriendlyFunctionalBootstrapHybridGiantSeededKeyGen(
    CKKSDenseFHEFriendlyFunctionalBootstrapHybridGiantSeededKey<Schedule>
        &bootstrap_key,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(bootstrap_key.linear_plan);
    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        rotation_usage, bootstrap_key.linear_plan);
    CKKSHybridSparseGaloisKeyChainGen<
        P, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>(
        *bootstrap_key.coeff_to_slot_galois,
        rotation_usage.coeff_to_slot_binary,
        rotation_usage.coeff_to_slot_direct, key, noise);
    CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
        *bootstrap_key.packed_conjugate_galois, key,
        rotation_usage.packed_conjugate, noise);
    CKKSFunctionalSeededRelinKeyChainGen<
        P, Schedule::functional_start_log_q, Schedule::eval_log_delta,
        Schedule::ExponentialTraits::relin_depth>(
        *bootstrap_key.exponential_relin, key, noise);
    CKKSFunctionalSeededRelinKeyChainGen<
        P, Schedule::ExponentialTraits::log_q, Schedule::eval_log_delta,
        Schedule::functional_double_angle>(*bootstrap_key.double_angle_relin,
                                           key, noise);
    CKKSFunctionalSeededRelinKeyChainGen<
        P, Schedule::after_exponential_log_q, Schedule::eval_log_delta,
        Schedule::LUTTraits::relin_depth>(*bootstrap_key.lut_relin, key,
                                          noise);
    CKKSRotationKeyIndexSet<P> conjugation_usage{};
    CKKSMarkConjugationKeyIndex<P>(conjugation_usage);
    CKKSSparseGaloisKeyGen<P, Schedule::after_evalmod_log_q>(
        *bootstrap_key.output_conjugate_galois, key, conjugation_usage, noise);
    CKKSHybridSparseGaloisKeyChainGen<
        P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *bootstrap_key.slot_to_coeff_galois,
        rotation_usage.slot_to_coeff_binary,
        rotation_usage.slot_to_coeff_direct, key, noise);
}

template <class Schedule, class BootstrapKey>
struct CKKSDenseFHEFriendlyFunctionalBootstrapInMemoryKeyProvider {
    const BootstrapKey *key;

    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan() const
    {
        return key->linear_plan;
    }
    template <std::size_t I>
    const auto &coeff_to_slot_galois() const
    {
        if constexpr (requires { *key->coeff_to_slot_galois; })
            return key->coeff_to_slot_galois->template get<I>();
        else
            return key->coeff_to_slot_galois.template get<I>();
    }
    const auto &packed_conjugate_galois() const
    {
        if constexpr (requires { *key->packed_conjugate_galois; })
            return *key->packed_conjugate_galois;
        else
            return key->packed_conjugate_galois;
    }
    template <std::size_t I>
    const auto &slot_to_coeff_galois() const
    {
        if constexpr (requires { *key->slot_to_coeff_galois; })
            return key->slot_to_coeff_galois->template get<I>();
        else
            return key->slot_to_coeff_galois.template get<I>();
    }
};

namespace ckks_detail {

template <class FunctionalKey>
inline decltype(auto) CKKSFunctionalExponentialRelinKey(
    const FunctionalKey &key)
{
    if constexpr (requires { *key.exponential_relin; })
        return *key.exponential_relin;
    else
        return (key.exponential_relin);
}

template <class FunctionalKey>
inline decltype(auto) CKKSFunctionalDoubleAngleRelinKey(
    const FunctionalKey &key)
{
    if constexpr (requires { *key.double_angle_relin; })
        return *key.double_angle_relin;
    else
        return (key.double_angle_relin);
}

template <class FunctionalKey>
inline decltype(auto) CKKSFunctionalLUTRelinKey(const FunctionalKey &key)
{
    if constexpr (requires { *key.lut_relin; })
        return *key.lut_relin;
    else
        return (key.lut_relin);
}

template <class FunctionalKey>
inline decltype(auto) CKKSFunctionalOutputConjugateKey(
    const FunctionalKey &key)
{
    if constexpr (requires { *key.output_conjugate_galois; })
        return *key.output_conjugate_galois;
    else
        return (key.output_conjugate_galois);
}

template <class Schedule, class FunctionalKey>
inline void CKKSDenseFHEFriendlyFunctionalBootstrapEvalComponent(
    typename Schedule::EvalModCiphertext &res,
    const typename Schedule::ComponentCiphertext &component,
    const CKKSFunctionalBootstrapLUT &lut,
    const CKKSFunctionalBootstrapComplexExponentialPolynomial &exponential,
    const FunctionalKey &keys)
{
    using P = typename Schedule::Param;
    assert(lut.plaintextModulus() >= 2);
    assert(lut.plaintextModulus() <= Schedule::lut_degree + 1);
    assert(exponential.degree == Schedule::exponential_degree);
    assert(exponential.double_angle == Schedule::functional_double_angle);
    assert(std::abs(exponential.input_bound -
                    Schedule::functional_input_bound) < 1e-12);

    auto scaled = std::make_unique<CKKSScaleDownResult<
        P, Schedule::after_component_split_log_q, Schedule::log_delta,
        Schedule::eval_scale_down_bits>>();
    CKKSScaleDown<P, Schedule::after_component_split_log_q,
                  Schedule::log_delta, Schedule::eval_scale_down_bits>(
        *scaled, component);

    auto exponential_value = std::make_unique<CKKSCiphertext<
        P, Schedule::after_exponential_log_q, Schedule::eval_log_delta>>();
    CKKSFunctionalBootstrapComplexExponential<
        P, Schedule::functional_start_log_q, Schedule::eval_log_delta,
        Schedule::exponential_coeff_log_delta, Schedule::exponential_degree,
        Schedule::functional_double_angle>(
        *exponential_value, *scaled, exponential,
        CKKSFunctionalExponentialRelinKey(keys),
        CKKSFunctionalDoubleAngleRelinKey(keys));
    scaled.reset();

    std::vector<std::complex<double>> coefficients = lut.coefficients;
    const double output_scale = 0.5 / Schedule::message_ratio;
    for (auto &coefficient : coefficients) coefficient *= output_scale;
    auto complex_lut =
        std::make_unique<typename Schedule::LUTTraits::Ciphertext>();
    CKKSEvalComplexPowerPolynomialWithKeyProvider<
        P, Schedule::after_exponential_log_q, Schedule::eval_log_delta,
        Schedule::lut_coeff_log_delta, Schedule::lut_degree>(
        *complex_lut, *exponential_value, coefficients,
        CKKSFunctionalLUTRelinKey(keys));
    exponential_value.reset();

    auto conjugated =
        std::make_unique<typename Schedule::LUTTraits::Ciphertext>();
    CKKSConjugateSlots<P, Schedule::after_evalmod_log_q>(
        conjugated->ct, complex_lut->ct,
        CKKSFunctionalOutputConjugateKey(keys));
    CKKSAddInPlace<P, Schedule::after_evalmod_log_q,
                   Schedule::eval_log_delta>(*complex_lut, *conjugated);
    CKKSScaleUpAtSameLevel<P, Schedule::after_evalmod_log_q,
                           Schedule::eval_log_delta, Schedule::log_delta>(
        res, *complex_lut);
}

template <class Schedule, bool ConsumeKey, class FunctionalKey,
          class KeyProvider>
inline void CKKSDenseFHEFriendlyFunctionalBootstrapImpl(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapLUT &lut,
    const CKKSFunctionalBootstrapComplexExponentialPolynomial &exponential,
    const FunctionalKey &functional_key,
    const KeyProvider &key_provider, CKKSDenseBootstrapTimings *timings,
    const CKKSDenseBootstrapProgress *progress)
{
    using P = typename Schedule::Param;
    const auto &linear_plan = key_provider.linear_plan();
    auto raised = std::make_unique<typename Schedule::BootstrapCiphertext>();
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->modraise_ms,
        [&] {
            CKKSModRaiseBoundedPhaseRandomized<
                P, Schedule::input_log_q, Schedule::boot_log_q,
                Schedule::log_delta, Schedule::modraise_mask_bound>(*raised,
                                                                     ct);
        },
        progress, "functional_fhe_friendly_modraise");
    auto slots = std::make_unique<typename Schedule::CoeffToSlotCiphertext>();
    const CKKSDenseBootstrapLinearKeyProviderChain<KeyProvider, true> c2s_keys{
        key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->coeff_to_slot_ms,
        [&] {
            CKKSDenseBootstrapCoeffToSlotStagesBSGS<Schedule>(
                *slots, *raised, linear_plan, c2s_keys);
        },
        progress, "functional_fhe_friendly_coeff_to_slot");
    raised.reset();
    if constexpr (ConsumeKey && requires {
                      functional_key.release_coeff_to_slot();
                  })
        functional_key.release_coeff_to_slot();

    auto real = std::make_unique<typename Schedule::ComponentCiphertext>();
    auto imag = std::make_unique<typename Schedule::ComponentCiphertext>();
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->split_ms,
        [&] {
            CKKSExtractRealImagSlots<
                P, Schedule::after_coeff_to_slot_log_q, Schedule::log_delta,
                Schedule::component_split_plain_log_delta>(
                *real, *imag, *slots,
                key_provider.packed_conjugate_galois());
        },
        progress, "functional_fhe_friendly_split");
    slots.reset();
    if constexpr (ConsumeKey && requires {
                      functional_key.release_packed_conjugate();
                  })
        functional_key.release_packed_conjugate();

    auto real_lut = std::make_unique<typename Schedule::EvalModCiphertext>();
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->real_evalmod_ms,
        [&] {
            CKKSDenseFHEFriendlyFunctionalBootstrapEvalComponent<Schedule>(
                *real_lut, *real, lut, exponential, functional_key);
        },
        progress, "functional_fhe_friendly_real_eval_lut");
    real.reset();
    auto imag_lut = std::make_unique<typename Schedule::EvalModCiphertext>();
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->imag_evalmod_ms,
        [&] {
            CKKSDenseFHEFriendlyFunctionalBootstrapEvalComponent<Schedule>(
                *imag_lut, *imag, lut, exponential, functional_key);
        },
        progress, "functional_fhe_friendly_imag_eval_lut");
    imag.reset();
    if constexpr (ConsumeKey && requires {
                      functional_key.release_eval_lut();
                  })
        functional_key.release_eval_lut();

    const CKKSDenseBootstrapLinearKeyProviderChain<KeyProvider, false> stc_keys{
        key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->slot_to_coeff_ms,
        [&] {
            CKKSDenseBootstrapSlotToCoeffStagesBSGSDualInputSharedTail<
                Schedule>(res, *real_lut, *imag_lut, linear_plan, stc_keys);
        },
        progress, "functional_fhe_friendly_slot_to_coeff");
    if constexpr (ConsumeKey && requires {
                      functional_key.release_slot_to_coeff();
                  })
        functional_key.release_slot_to_coeff();
}

}  // namespace ckks_detail

template <class Schedule, class FunctionalKey>
inline void CKKSDenseFHEFriendlyFunctionalBootstrapTimed(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapLUT &lut,
    const CKKSFunctionalBootstrapComplexExponentialPolynomial &exponential,
    const FunctionalKey &key,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    timings = {};
    const CKKSDenseFHEFriendlyFunctionalBootstrapInMemoryKeyProvider<
        Schedule, FunctionalKey>
        provider{&key};
    ckks_detail::CKKSDenseFHEFriendlyFunctionalBootstrapImpl<Schedule, false>(
        res, ct, lut, exponential, key, provider, &timings, progress);
}

// Releases phase-owned keys when FunctionalKey provides release methods. Such
// a key cannot be reused or serialized after this call.
template <class Schedule, class FunctionalKey>
inline void CKKSDenseFHEFriendlyFunctionalBootstrapConsumeTimed(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapLUT &lut,
    const CKKSFunctionalBootstrapComplexExponentialPolynomial &exponential,
    const FunctionalKey &key, CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    timings = {};
    const CKKSDenseFHEFriendlyFunctionalBootstrapInMemoryKeyProvider<
        Schedule, FunctionalKey>
        provider{&key};
    ckks_detail::CKKSDenseFHEFriendlyFunctionalBootstrapImpl<Schedule, true>(
        res, ct, lut, exponential, key, provider, &timings, progress);
}

template <class Schedule, class FunctionalKey>
inline void CKKSDenseFHEFriendlyFunctionalBootstrap(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSFunctionalBootstrapLUT &lut,
    const CKKSFunctionalBootstrapComplexExponentialPolynomial &exponential,
    const FunctionalKey &key)
{
    const CKKSDenseFHEFriendlyFunctionalBootstrapInMemoryKeyProvider<
        Schedule, FunctionalKey>
        provider{&key};
    ckks_detail::CKKSDenseFHEFriendlyFunctionalBootstrapImpl<Schedule, false>(
        res, ct, lut, exponential, key, provider, nullptr, nullptr);
}

}  // namespace TFHEpp

#endif  // TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
