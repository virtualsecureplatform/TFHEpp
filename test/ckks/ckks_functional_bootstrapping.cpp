#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tfhe++.hpp>
#include <vector>

namespace {

void require(bool condition, const char *label)
{
    if (!condition) {
        std::cerr << label << std::endl;
        std::exit(1);
    }
}

void require_close(double got, double want, double tolerance,
                   const char *label)
{
    if (std::abs(got - want) > tolerance) {
        std::cerr << label << " got=" << got << " want=" << want
                  << " tolerance=" << tolerance << std::endl;
        std::exit(1);
    }
}

void test_exact_trigonometric_hermite_lut()
{
    const std::vector<double> values{3.0, -2.0, 7.0, 1.0, -5.0, 4.0, 9.0};
    const auto lut = TFHEpp::CKKSBuildFunctionalBootstrapLUT(values);
    require(lut.plaintextModulus() == values.size(), "LUT plaintext modulus");

    for (std::size_t k = 0; k < values.size(); k++) {
        const double x = static_cast<double>(k) /
                         static_cast<double>(values.size());
        require_close(TFHEpp::CKKSEvaluateFunctionalBootstrapLUT(lut, x),
                      values[k], 2e-12, "Hermite interpolation value");
        require_close(
            TFHEpp::CKKSEvaluateFunctionalBootstrapLUTDerivative(lut, x),
            0.0, 2e-10, "Hermite interpolation zero first derivative");
        require_close(TFHEpp::CKKSEvaluateFunctionalBootstrapLUT(lut, x + 3.0),
                      values[k], 2e-12, "Hermite interpolation periodicity");
    }

    // The zero first derivative makes the perturbation error quadratic.
    const double x = 3.0 / static_cast<double>(values.size());
    const double epsilon = 1e-5;
    const double error = std::abs(
        TFHEpp::CKKSEvaluateFunctionalBootstrapLUT(lut, x + epsilon) -
        values[3]);
    require(error < 1e4 * epsilon * epsilon,
            "first-order Hermite LUT quadratically cleans noise");

    bool rejected = false;
    try {
        (void)TFHEpp::CKKSBuildFunctionalBootstrapLUT({1.0});
    }
    catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "reject one-entry functional-bootstrap LUT");
}

void test_chebyshev_eval_lut_reference()
{
    // For p=2 and outputs 0,1, Equation (4) is sin^2(pi*x).
    const auto lut =
        TFHEpp::CKKSBuildFunctionalBootstrapLUT({0.0, 1.0});
    const auto polynomial =
        TFHEpp::CKKSBuildFunctionalBootstrapChebyshevPolynomial(lut, 24);
    for (int i = -100; i <= 100; i++) {
        const double x = static_cast<double>(i) / 100.0;
        const double want = std::sin(3.14159265358979323846 * x) *
                            std::sin(3.14159265358979323846 * x);
        require_close(
            TFHEpp::CKKSEvaluateFunctionalBootstrapChebyshevPolynomial(
                polynomial, x),
            want, 2e-10, "Chebyshev EvalLUT reference");
    }
}

void test_fhe_friendly_complex_exponential_reference()
{
    const std::vector<double> values{3.0, -2.0, 7.0, 1.0,
                                     -5.0, 4.0,  9.0};
    const auto lut = TFHEpp::CKKSBuildFunctionalBootstrapLUT(values);
    constexpr double input_bound = 32.0;
    const auto exponential =
        TFHEpp::CKKSBuildFunctionalBootstrapComplexExponentialPolynomial(
            15, 8, input_bound);
    for (int overflow = -2; overflow <= 2; overflow++) {
        for (std::size_t message = 0; message < values.size(); message++) {
            const double x = static_cast<double>(overflow * 7) +
                             static_cast<double>(message) / values.size();
            require_close(
                TFHEpp::CKKSEvaluateFunctionalBootstrapLUTFHEFriendly(
                    lut, exponential, x),
                values[message], 2e-8,
                "FHE-friendly complex-exponential EvalLUT reference");
        }
    }

    using Schedule = TFHEpp::lvl6CKKSDenseFunctionalBootstrapP8Schedule;
    const auto lvl6_lut = TFHEpp::CKKSBuildFunctionalBootstrapLUT(
        {-0.25, 0.125, 0.5, -0.375, 0.25, 0.0, -0.125, 0.375});
    const auto lvl6_exponential =
        TFHEpp::CKKSBuildFunctionalBootstrapComplexExponentialPolynomial(
            Schedule::exponential_degree, Schedule::functional_double_angle,
            Schedule::functional_input_bound);
    static_assert(Schedule::after_evalmod_log_q == 83);
    static_assert(Schedule::output_log_q == 55);
    for (int overflow = -18; overflow <= 18; overflow++) {
        for (std::size_t message = 0; message < lvl6_lut.values.size();
             message++) {
            const double x =
                Schedule::message_ratio * static_cast<double>(overflow) +
                static_cast<double>(message) / lvl6_lut.values.size();
            require_close(
                TFHEpp::CKKSEvaluateFunctionalBootstrapLUTFHEFriendly(
                    lvl6_lut, lvl6_exponential, x),
                lvl6_lut.values[message], 2e-8,
                "lvl6 p8 FHE-friendly EvalLUT approximation");
        }
    }
}

struct TinyFunctionalBootstrapParam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 4;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 1;
    static constexpr std::uint32_t l_a = 1;
    static constexpr std::uint32_t lₐ = l_a;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t Bg_abit = 16;
    static constexpr std::uint32_t Bgₐbit = Bg_abit;
    using T = TFHEpp::MultiLimbUInt<5>;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bg_a = T{1} << Bg_abit;
    static constexpr T Bgₐ = Bg_a;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = 0.0;
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 20;
    static constexpr std::uint32_t l̅ₐ = 20;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;
};

void test_homomorphic_eval_lut()
{
    using P = TinyFunctionalBootstrapParam;
    constexpr std::uint32_t log_q = 190;
    constexpr std::uint32_t log_delta = 24;
    constexpr std::uint32_t coeff_log_delta = 24;
    constexpr std::size_t degree = 15;
    using Input = TFHEpp::CKKSCiphertext<P, log_q, log_delta>;
    using Output = TFHEpp::CKKSFunctionalBootstrapEvalLUTResult<
        P, log_q, log_delta, coeff_log_delta, degree>;
    using EvalKey = TFHEpp::CKKSFunctionalBootstrapEvalLUTRelinKey<
        P, log_q, log_delta, coeff_log_delta, degree>;

    const auto lut =
        TFHEpp::CKKSBuildFunctionalBootstrapLUT({0.0, 1.0});
    const auto polynomial =
        TFHEpp::CKKSBuildFunctionalBootstrapChebyshevPolynomial(lut, degree);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    for (std::size_t i = 0; i < P::n; i++)
        (*key)[i] = static_cast<P::T>(static_cast<int>(i % 3) - 1);
    auto eval_key = std::make_unique<EvalKey>();
    TFHEpp::CKKSFunctionalBootstrapEvalLUTKeyGen<
        P, log_q, log_delta, coeff_log_delta, degree>(*eval_key, *key,
                                                       {0.0, 0});

    TFHEpp::CKKSSlotVector<P> input_slots{};
    for (std::size_t i = 0; i < input_slots.size(); i++)
        input_slots[i] = {static_cast<double>(static_cast<int>(i % 5) - 2) /
                              2.0,
                          0.0};
    auto input = std::make_unique<Input>();
    TFHEpp::ckksSlotEncrypt<P, log_q, log_delta>(*input, input_slots, *key,
                                                 {0.0, 0});

    auto output = std::make_unique<Output>();
    TFHEpp::CKKSFunctionalBootstrapEvalLUTNormalized<
        P, log_q, log_delta, coeff_log_delta, degree>(
        *output, *input, polynomial, *eval_key);

    TFHEpp::CKKSSlotVector<P> decoded{};
    TFHEpp::ckksSlotDecrypt<P>(decoded, *output, *key);
    for (std::size_t i = 0; i < decoded.size(); i++) {
        const double x = input_slots[i].real();
        const double want = TFHEpp::CKKSEvaluateFunctionalBootstrapLUT(lut, x);
        require_close(decoded[i].real(), want, 0.015,
                      "homomorphic functional-bootstrap EvalLUT");
        require_close(decoded[i].imag(), 0.0, 0.015,
                      "homomorphic functional-bootstrap imaginary error");
    }
}

void test_homomorphic_fhe_friendly_eval_lut()
{
    using P = TinyFunctionalBootstrapParam;
    constexpr std::uint32_t start_log_q = 310;
    constexpr std::uint32_t log_delta = 20;
    constexpr std::uint32_t coeff_log_delta = 20;
    constexpr std::size_t exponential_degree = 15;
    constexpr std::uint32_t double_angle = 4;
    constexpr std::size_t lut_degree = 3;
    constexpr double input_bound = 4.0;
    using ExpTraits = TFHEpp::CKKSPowerPolynomialEvaluatorTraits<
        P, start_log_q, log_delta, coeff_log_delta, exponential_degree>;
    constexpr std::uint32_t after_exponential_log_q =
        ExpTraits::log_q - double_angle * log_delta;
    using LUTTraits = TFHEpp::CKKSPowerPolynomialEvaluatorTraits<
        P, after_exponential_log_q, log_delta, coeff_log_delta, lut_degree>;

    const auto lut =
        TFHEpp::CKKSBuildFunctionalBootstrapLUT({-0.2, 0.1, 0.3, -0.1});
    const auto exponential =
        TFHEpp::CKKSBuildFunctionalBootstrapComplexExponentialPolynomial(
            exponential_degree, double_angle, input_bound);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    for (std::size_t i = 0; i < P::n; i++)
        (*key)[i] = static_cast<P::T>(static_cast<int>(i % 3) - 1);
    auto exponential_keys =
        std::make_unique<typename ExpTraits::RelinKeyChain>();
    TFHEpp::CKKSRelinKeyChainGen<P, start_log_q, log_delta,
                                 ExpTraits::relin_depth>(
        *exponential_keys, *key, {0.0, 0});
    auto double_angle_keys =
        std::make_unique<TFHEpp::CKKSRelinKeyChain<
            P, ExpTraits::log_q, log_delta, double_angle>>();
    TFHEpp::CKKSRelinKeyChainGen<P, ExpTraits::log_q, log_delta,
                                 double_angle>(*double_angle_keys, *key,
                                               {0.0, 0});
    auto lut_keys = std::make_unique<typename LUTTraits::RelinKeyChain>();
    TFHEpp::CKKSRelinKeyChainGen<P, after_exponential_log_q, log_delta,
                                 LUTTraits::relin_depth>(*lut_keys, *key,
                                                        {0.0, 0});
    TFHEpp::CKKSRotationKeyIndexSet<P> conjugation_usage{};
    TFHEpp::CKKSMarkConjugationKeyIndex<P>(conjugation_usage);
    auto conjugation_key =
        std::make_unique<TFHEpp::CKKSSparseGaloisKey<P, LUTTraits::log_q>>();
    TFHEpp::CKKSSparseGaloisKeyGen<P, LUTTraits::log_q>(
        *conjugation_key, *key, conjugation_usage, {0.0, 0});

    TFHEpp::CKKSSlotVector<P> inputs{};
    TFHEpp::CKKSSlotVector<P> expected{};
    for (std::size_t i = 0; i < inputs.size(); i++) {
        const std::size_t message = (3 * i + 1) % lut.values.size();
        const int overflow = static_cast<int>(i % 3) - 1;
        const double x = static_cast<double>(message) / lut.values.size() +
                         static_cast<double>(overflow);
        inputs[i] = {x / input_bound, 0.0};
        expected[i] = {lut.values[message], 0.0};
    }
    auto input = std::make_unique<
        TFHEpp::CKKSCiphertext<P, start_log_q, log_delta>>();
    TFHEpp::ckksSlotEncrypt<P, start_log_q, log_delta>(
        *input, inputs, *key, {0.0, 0});
    auto exponential_value = std::make_unique<TFHEpp::CKKSCiphertext<
        P, after_exponential_log_q, log_delta>>();
    TFHEpp::CKKSFunctionalBootstrapComplexExponential<
        P, start_log_q, log_delta, coeff_log_delta, exponential_degree,
        double_angle>(*exponential_value, *input, exponential,
                      *exponential_keys, *double_angle_keys);

    auto half_coefficients = lut.coefficients;
    for (auto &coefficient : half_coefficients) coefficient *= 0.5;
    auto complex_lut = std::make_unique<typename LUTTraits::Ciphertext>();
    TFHEpp::CKKSEvalComplexPowerPolynomialWithKeyProvider<
        P, after_exponential_log_q, log_delta, coeff_log_delta, lut_degree>(
        *complex_lut, *exponential_value, half_coefficients, *lut_keys);
    auto conjugated = std::make_unique<typename LUTTraits::Ciphertext>();
    TFHEpp::CKKSConjugateSlots<P, LUTTraits::log_q>(
        conjugated->ct, complex_lut->ct, *conjugation_key);
    TFHEpp::CKKSAddInPlace<P, LUTTraits::log_q, log_delta>(*complex_lut,
                                                           *conjugated);

    TFHEpp::CKKSSlotVector<P> decoded{};
    TFHEpp::ckksSlotDecrypt<P>(decoded, *complex_lut, *key);
    for (std::size_t i = 0; i < decoded.size(); i++) {
        require_close(decoded[i].real(), expected[i].real(), 0.015,
                      "homomorphic FHE-friendly EvalLUT");
        require_close(decoded[i].imag(), 0.0, 0.015,
                      "homomorphic FHE-friendly EvalLUT imaginary error");
    }
}

void test_dense_functional_bootstrap_e2e()
{
    using P = TinyFunctionalBootstrapParam;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        P, 16, 1, 310, 16, 3, 63, 2, 1, 0, 16, 2>;
    static_assert(Schedule::coeff_to_slot_level_count == 1);
    static_assert(Schedule::slot_to_coeff_level_count == 1);

    const auto lut =
        TFHEpp::CKKSBuildFunctionalBootstrapLUT({-0.125, 0.25});
    const auto polynomial =
        TFHEpp::CKKSBuildDenseFunctionalBootstrapChebyshevPolynomial<
            Schedule>(lut);

    // Check all representatives that a weight-one modulus raise can produce.
    for (int overflow = -1; overflow <= 1; overflow++) {
        for (int message = 0; message < 2; message++) {
            const double x = static_cast<double>(message) / 2.0 +
                             Schedule::message_ratio * overflow;
            const double normalized = x / polynomial.input_bound;
            const double want = lut.values[message] / Schedule::message_ratio;
            require_close(
                TFHEpp::CKKSEvaluateFunctionalBootstrapChebyshevPolynomialNormalized(
                    polynomial, normalized),
                want, 2e-10, "dense functional EvalLUT approximation");
        }
    }

    auto key = std::make_unique<TFHEpp::Key<P>>();
    key->fill(P::T{0});
    (*key)[0] = P::T{1};

    auto bootstrap_key =
        std::make_unique<TFHEpp::CKKSDenseBootstrapKey<Schedule>>();
    TFHEpp::CKKSDenseBootstrapKeyGen<Schedule>(*bootstrap_key, *key,
                                               {0.0, 0});

    std::array<double, P::n> messages{};
    std::array<double, P::n> expected{};
    for (std::size_t i = 0; i < P::n; i++) {
        const std::size_t message = (3 * i + 1) % 2;
        messages[i] = static_cast<double>(message) / 2.0;
        expected[i] = lut.values[message];
    }

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    TFHEpp::ckksEncrypt<P, Schedule::input_log_q, Schedule::log_delta>(
        *input, messages, *key, {0.0, 0});
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapTimings timings;
    TFHEpp::CKKSDenseFunctionalBootstrapTimed<Schedule>(
        *output, *input, polynomial, *bootstrap_key, timings);

    std::array<double, P::n> decoded{};
    TFHEpp::ckksDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        decoded, *output, *key);
    for (std::size_t i = 0; i < P::n; i++)
        require_close(decoded[i], expected[i], 0.01,
                      "dense functional bootstrap end to end");
    require(timings.total_ms() > 0.0,
            "dense functional bootstrap records timings");
    std::cout << "dense_functional_bootstrap_ms=" << timings.total_ms()
              << " amortized_ms=" << timings.total_ms() / P::n << std::endl;
}

void test_dense_fhe_friendly_functional_bootstrap_e2e()
{
    using P = TinyFunctionalBootstrapParam;
    using BaseSchedule = TFHEpp::CKKSDenseBootstrapSchedule<
        P, 20, 1, 310, 16, 3, 63, 2, 1, 0, 16, 2>;
    using Schedule =
        TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapSchedule<
            BaseSchedule, 16, 16, 15, 4, 16, 3>;
    static_assert(Schedule::output_log_q == 66);

    const auto lut = TFHEpp::CKKSBuildFunctionalBootstrapLUT(
        {-0.2, 0.1, 0.3, -0.1});
    const auto exponential =
        TFHEpp::CKKSBuildFunctionalBootstrapComplexExponentialPolynomial(
            Schedule::exponential_degree, Schedule::functional_double_angle,
            Schedule::functional_input_bound);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    key->fill(P::T{0});
    (*key)[0] = P::T{1};
    auto bootstrap_key = std::make_unique<
        TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapKey<Schedule>>();
    TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapKeyGen<Schedule>(
        *bootstrap_key, *key, {0.0, 0});

    std::array<double, P::n> messages{};
    std::array<double, P::n> expected{};
    for (std::size_t i = 0; i < P::n; i++) {
        const std::size_t message = (3 * i + 1) % lut.values.size();
        messages[i] = static_cast<double>(message) / lut.values.size();
        expected[i] = lut.values[message];
    }
    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    TFHEpp::ckksEncrypt<P, Schedule::input_log_q, Schedule::log_delta>(
        *input, messages, *key, {0.0, 0});
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapTimings timings;
    TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapTimed<Schedule>(
        *output, *input, lut, exponential, *bootstrap_key, timings);

    std::array<double, P::n> decoded{};
    TFHEpp::ckksDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        decoded, *output, *key);
    for (std::size_t i = 0; i < P::n; i++)
        require_close(decoded[i], expected[i], 0.04,
                      "dense FHE-friendly functional bootstrap end to end");
    require(timings.total_ms() > 0.0,
            "dense FHE-friendly bootstrap records timings");
    std::cout << "dense_fhe_friendly_functional_bootstrap_ms="
              << timings.total_ms()
              << " amortized_ms=" << timings.total_ms() / P::n << std::endl;

    using HybridKey =
        TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapHybridGiantSeededKey<
            Schedule>;
    auto hybrid_key = std::make_unique<HybridKey>();
    TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapHybridGiantSeededKeyGen<
        Schedule>(*hybrid_key, *key, {0.0, 0});
    TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapConsumeTimed<Schedule>(
        *output, *input, lut, exponential, *hybrid_key, timings);
    TFHEpp::ckksDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        decoded, *output, *key);
    for (std::size_t i = 0; i < P::n; i++)
        require_close(decoded[i], expected[i], 0.04,
                      "consuming hybrid functional bootstrap end to end");
    require(!hybrid_key->coeff_to_slot_galois &&
                !hybrid_key->packed_conjugate_galois &&
                !hybrid_key->exponential_relin &&
                !hybrid_key->double_angle_relin && !hybrid_key->lut_relin &&
                !hybrid_key->output_conjugate_galois &&
                !hybrid_key->slot_to_coeff_galois,
            "consuming functional bootstrap releases phase keys");
}

template <class F>
double average_runtime_ms(std::size_t repeats, F &&operation)
{
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; i++) operation();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count() /
           static_cast<double>(repeats);
}

void benchmark_functional_bootstrap()
{
    using P = TinyFunctionalBootstrapParam;
    using DirectSchedule = TFHEpp::CKKSDenseBootstrapSchedule<
        P, 20, 1, 310, 16, 3, 63, 2, 1, 0, 16, 2>;
    using FHEFriendlySchedule =
        TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapSchedule<
            DirectSchedule, 16, 16, 15, 4, 16, 3>;
    constexpr std::size_t repeats = 5;

    const auto lut = TFHEpp::CKKSBuildFunctionalBootstrapLUT(
        {-0.2, 0.1, 0.3, -0.1});
    const auto direct_polynomial =
        TFHEpp::CKKSBuildDenseFunctionalBootstrapChebyshevPolynomial<
            DirectSchedule>(lut);
    const auto exponential =
        TFHEpp::CKKSBuildFunctionalBootstrapComplexExponentialPolynomial(
            FHEFriendlySchedule::exponential_degree,
            FHEFriendlySchedule::functional_double_angle,
            FHEFriendlySchedule::functional_input_bound);
    auto key = std::make_unique<TFHEpp::Key<P>>();
    key->fill(P::T{0});
    (*key)[0] = P::T{1};
    auto direct_key =
        std::make_unique<TFHEpp::CKKSDenseBootstrapKey<DirectSchedule>>();
    TFHEpp::CKKSDenseBootstrapKeyGen<DirectSchedule>(*direct_key, *key,
                                                     {0.0, 0});
    auto fhe_friendly_key = std::make_unique<
        TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapKey<
            FHEFriendlySchedule>>();
    TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapKeyGen<
        FHEFriendlySchedule>(*fhe_friendly_key, *key, {0.0, 0});
    std::array<double, P::n> messages{};
    for (std::size_t i = 0; i < P::n; i++)
        messages[i] = static_cast<double>(i % lut.values.size()) /
                      lut.values.size();
    auto input =
        std::make_unique<typename DirectSchedule::InputCiphertext>();
    auto direct_output =
        std::make_unique<typename DirectSchedule::OutputCiphertext>();
    auto fhe_friendly_output =
        std::make_unique<typename FHEFriendlySchedule::OutputCiphertext>();
    TFHEpp::ckksEncrypt<P, DirectSchedule::input_log_q,
                       DirectSchedule::log_delta>(
        *input, messages, *key, {0.0, 0});
    TFHEpp::CKKSDenseFunctionalBootstrap<DirectSchedule>(
        *direct_output, *input, direct_polynomial, *direct_key);
    const double direct_total_ms = average_runtime_ms(repeats, [&] {
        TFHEpp::CKKSDenseFunctionalBootstrap<DirectSchedule>(
            *direct_output, *input, direct_polynomial, *direct_key);
    });
    TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrap<FHEFriendlySchedule>(
        *fhe_friendly_output, *input, lut, exponential, *fhe_friendly_key);
    const double fhe_friendly_total_ms = average_runtime_ms(repeats, [&] {
        TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrap<FHEFriendlySchedule>(
            *fhe_friendly_output, *input, lut, exponential,
            *fhe_friendly_key);
    });
    const double direct_amortized_ms = direct_total_ms / P::n;
    const double fhe_friendly_amortized_ms = fhe_friendly_total_ms / P::n;

    using BootstrapP = TFHEpp::lvl02param;
    using IKSP = TFHEpp::lvl20param;
    auto tfhe_secret = std::make_unique<TFHEpp::SecretKey>();
    auto tfhe_eval = std::make_unique<TFHEpp::EvalKey>();
    tfhe_eval->emplacebkfft<BootstrapP>(*tfhe_secret);
    tfhe_eval->emplaceiksk<IKSP>(*tfhe_secret);
    std::array<TFHEpp::TLWE<typename IKSP::domainP>, P::n> tfhe_inputs{};
    std::array<TFHEpp::TLWE<typename BootstrapP::targetP>, P::n>
        tfhe_outputs{};
    for (std::size_t i = 0; i < P::n; i++)
        TFHEpp::tlweSymEncrypt<typename IKSP::domainP>(
            tfhe_inputs[i], (i & 1) ? IKSP::domainP::μ : -IKSP::domainP::μ,
            tfhe_secret->key.getSubset<typename IKSP::domainP>());
    const double tfhe_total_ms = average_runtime_ms(2, [&] {
        for (std::size_t i = 0; i < P::n; i++)
            TFHEpp::GateBootstrapping<IKSP, BootstrapP,
                                      BootstrapP::targetP::μ>(
                tfhe_outputs[i], tfhe_inputs[i], *tfhe_eval);
    });
    const double tfhe_amortized_ms = tfhe_total_ms / P::n;

    std::cout << "benchmark,security,packed_values,total_ms,amortized_ms\n";
    std::cout << "ckks_direct_functional_bootstrap,toy," << P::n << ','
              << direct_total_ms << ',' << direct_amortized_ms << '\n';
    std::cout << "ckks_fhe_friendly_functional_bootstrap,toy," << P::n << ','
              << fhe_friendly_total_ms << ',' << fhe_friendly_amortized_ms
              << '\n';
    std::cout << "tfhe_gate_bootstrap,default," << P::n << ','
              << tfhe_total_ms << ',' << tfhe_amortized_ms << '\n';
    std::cout << "fhe_friendly_over_direct_speedup="
              << direct_amortized_ms / fhe_friendly_amortized_ms << '\n';
    std::cout << "fhe_friendly_over_tfhe_speedup="
              << tfhe_amortized_ms / fhe_friendly_amortized_ms << '\n';
    std::cout << "warning=parameter sets have different security; speedup is "
                 "a throughput smoke result, not a secure comparison\n";
}

void lvl6_stage_begin(const char *stage, const void *)
{
    std::cout << "lvl6_stage_begin=" << stage << std::endl;
}

void lvl6_stage_end(const char *stage, double elapsed_ms, const void *)
{
    std::cout << "lvl6_stage_end=" << stage
              << " elapsed_ms=" << elapsed_ms << std::endl;
}

void benchmark_lvl6_functional_bootstrap(
    const std::filesystem::path &key_file = {})
{
    using Schedule = TFHEpp::lvl6CKKSDenseFunctionalBootstrapP8Schedule;
    using P = typename Schedule::Param;
    using FunctionalKey =
        TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapHybridGiantSeededKey<
            Schedule>;
    constexpr std::size_t sparse_weight = 16;
    constexpr std::size_t stride = 7919;

    std::cout << "lvl6_n=" << P::n
              << " boot_log_q=" << Schedule::boot_log_q
              << " input_log_q=" << Schedule::input_log_q
              << " output_log_q=" << Schedule::output_log_q
              << " sparse_weight=" << sparse_weight << std::endl;
    auto key = std::make_unique<TFHEpp::Key<P>>();
    key->fill(typename P::T{0});
    for (std::size_t i = 0; i < sparse_weight; i++) {
        const std::size_t index = (i * stride) % P::n;
        const int sign = (i & 1) == 0 ? 1 : -1;
        (*key)[index] = static_cast<typename P::T>(sign);
    }

    std::unique_ptr<FunctionalKey> functional_key;
    const auto key_start = std::chrono::steady_clock::now();
    if (!key_file.empty() && std::filesystem::exists(key_file)) {
        std::cout << "lvl6_key_load=" << key_file << std::endl;
        functional_key = std::make_unique<FunctionalKey>();
        TFHEpp::CKKSLoadPortableBinary(*functional_key, key_file);
    }
    else {
        std::cout << "lvl6_keygen_begin=1" << std::endl;
        auto linear_plan =
            std::make_unique<TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule>>();
        TFHEpp::CKKSBuildDenseBootstrapLinearPlan<Schedule>(*linear_plan);
        std::cout << "lvl6_keygen_phase=linear_plan" << std::endl;
        functional_key = std::make_unique<FunctionalKey>();
        functional_key->linear_plan = std::move(*linear_plan);
        linear_plan.reset();
        auto rotation_usage = std::make_unique<
            TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>>();
        TFHEpp::CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
            *rotation_usage, functional_key->linear_plan);
        std::cout << "lvl6_keygen_phase=rotation_usage" << std::endl;
        TFHEpp::CKKSHybridSparseGaloisKeyChainGen<
            P, Schedule::boot_log_q,
            Schedule::coeff_to_slot_plain_log_delta,
            Schedule::coeff_to_slot_level_count>(
            *functional_key->coeff_to_slot_galois,
            rotation_usage->coeff_to_slot_binary,
            rotation_usage->coeff_to_slot_direct, *key);
        std::cout << "lvl6_keygen_phase=coeff_to_slot" << std::endl;
        TFHEpp::CKKSSparseGaloisKeyGen<
            P, Schedule::after_coeff_to_slot_log_q>(
            *functional_key->packed_conjugate_galois, *key,
            rotation_usage->packed_conjugate);
        std::cout << "lvl6_keygen_phase=packed_conjugate" << std::endl;
        TFHEpp::CKKSFunctionalSeededRelinKeyChainGen<
            P, Schedule::functional_start_log_q, Schedule::eval_log_delta,
            Schedule::ExponentialTraits::relin_depth>(
            *functional_key->exponential_relin, *key);
        std::cout << "lvl6_keygen_phase=exponential_relin" << std::endl;
        TFHEpp::CKKSFunctionalSeededRelinKeyChainGen<
            P, Schedule::ExponentialTraits::log_q,
            Schedule::eval_log_delta, Schedule::functional_double_angle>(
            *functional_key->double_angle_relin, *key);
        std::cout << "lvl6_keygen_phase=double_angle_relin" << std::endl;
        TFHEpp::CKKSFunctionalSeededRelinKeyChainGen<
            P, Schedule::after_exponential_log_q,
            Schedule::eval_log_delta, Schedule::LUTTraits::relin_depth>(
            *functional_key->lut_relin, *key);
        std::cout << "lvl6_keygen_phase=lut_relin" << std::endl;
        TFHEpp::CKKSRotationKeyIndexSet<P> output_conjugation_usage{};
        TFHEpp::CKKSMarkConjugationKeyIndex<P>(output_conjugation_usage);
        TFHEpp::CKKSSparseGaloisKeyGen<P, Schedule::after_evalmod_log_q>(
            *functional_key->output_conjugate_galois, *key,
            output_conjugation_usage);
        std::cout << "lvl6_keygen_phase=output_conjugate" << std::endl;
        TFHEpp::CKKSHybridSparseGaloisKeyChainGen<
            P, Schedule::after_evalmod_log_q,
            Schedule::slot_to_coeff_plain_log_delta,
            Schedule::slot_to_coeff_level_count>(
            *functional_key->slot_to_coeff_galois,
            rotation_usage->slot_to_coeff_binary,
            rotation_usage->slot_to_coeff_direct, *key);
        std::cout << "lvl6_keygen_phase=slot_to_coeff" << std::endl;
        rotation_usage.reset();
        std::cout << "lvl6_keygen_complete=1" << std::endl;
        if (!key_file.empty()) {
            std::cout << "lvl6_key_save=" << key_file << std::endl;
            TFHEpp::CKKSSavePortableBinaryAtomic(key_file, *functional_key);
        }
    }
    const double key_ms = std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - key_start)
                              .count();
    std::cout << "lvl6_key_ready_ms=" << key_ms << std::endl;

    const auto lut = TFHEpp::CKKSBuildFunctionalBootstrapLUT(
        {-0.25, 0.125, 0.5, -0.375, 0.25, 0.0, -0.125, 0.375});
    const auto exponential =
        TFHEpp::CKKSBuildFunctionalBootstrapComplexExponentialPolynomial(
            Schedule::exponential_degree, Schedule::functional_double_angle,
            Schedule::functional_input_bound);
    auto messages = std::make_unique<std::array<double, P::n>>();
    auto expected = std::make_unique<std::array<double, P::n>>();
    for (std::size_t i = 0; i < P::n; i++) {
        const std::size_t message = (5 * i + 3) % lut.values.size();
        (*messages)[i] = static_cast<double>(message) / lut.values.size();
        (*expected)[i] = lut.values[message];
    }
    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    TFHEpp::ckksEncrypt<P, Schedule::input_log_q, Schedule::log_delta>(
        *input, *messages, *key);
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapTimings timings;
    const TFHEpp::CKKSDenseBootstrapProgress progress{
        lvl6_stage_begin, lvl6_stage_end, nullptr};
    TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapConsumeTimed<Schedule>(
        *output, *input, lut, exponential, *functional_key, timings,
        &progress);

    auto decoded = std::make_unique<std::array<double, P::n>>();
    TFHEpp::ckksDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *output, *key);
    double max_error = 0.0;
    for (std::size_t i = 0; i < P::n; i++)
        max_error =
            std::max(max_error, std::abs((*decoded)[i] - (*expected)[i]));
    std::cout << "lvl6_functional_bootstrap_ms=" << timings.total_ms()
              << " amortized_ms=" << timings.total_ms() / P::n
              << " max_error=" << max_error << std::endl;
    require(max_error < 0.2, "lvl6 functional bootstrap correctness");

    // Free the multi-gigabyte CKKS key before constructing the TFHE baseline.
    functional_key.reset();
    using BootstrapP = TFHEpp::lvl02param;
    using IKSP = TFHEpp::lvl20param;
    auto tfhe_secret = std::make_unique<TFHEpp::SecretKey>();
    auto tfhe_eval = std::make_unique<TFHEpp::EvalKey>();
    tfhe_eval->emplacebkfft<BootstrapP>(*tfhe_secret);
    tfhe_eval->emplaceiksk<IKSP>(*tfhe_secret);
    TFHEpp::TLWE<typename IKSP::domainP> tfhe_input{};
    TFHEpp::TLWE<typename BootstrapP::targetP> tfhe_output{};
    TFHEpp::tlweSymEncrypt<typename IKSP::domainP>(
        tfhe_input, IKSP::domainP::μ,
        tfhe_secret->key.getSubset<typename IKSP::domainP>());
    const double tfhe_ms = average_runtime_ms(16, [&] {
        TFHEpp::GateBootstrapping<IKSP, BootstrapP, BootstrapP::targetP::μ>(
            tfhe_output, tfhe_input, *tfhe_eval);
    });
    std::cout << "tfhe_gate_bootstrap_ms=" << tfhe_ms
              << " lvl6_over_tfhe_amortized_speedup="
              << tfhe_ms / (timings.total_ms() / P::n) << std::endl;
}

}  // namespace

int main(int argc, char **argv)
{
    if (argc == 2 && std::string(argv[1]) == "--benchmark") {
        benchmark_functional_bootstrap();
        return 0;
    }
    if (argc >= 2 && std::string(argv[1]) == "--benchmark-lvl6") {
        benchmark_lvl6_functional_bootstrap(argc >= 3 ? argv[2] : "");
        return 0;
    }
    test_exact_trigonometric_hermite_lut();
    test_chebyshev_eval_lut_reference();
    test_fhe_friendly_complex_exponential_reference();
    test_homomorphic_eval_lut();
    test_homomorphic_fhe_friendly_eval_lut();
    test_dense_functional_bootstrap_e2e();
    test_dense_fhe_friendly_functional_bootstrap_e2e();
    std::cout << "Passed" << std::endl;
}
