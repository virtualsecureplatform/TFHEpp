#include <ckks/narrow_bootstrapping.hpp>
#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <ckks/narrow_frequency_basis.hpp>
#include <iostream>
#include <sstream>
#include <chrono>
#include <fstream>

namespace {

template <std::uint32_t Nbit>
struct TestParam : TFHEpp::lvl5param {
    static constexpr std::uint32_t nbit = Nbit;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t Bgbit = 12;
    static constexpr T Bg = T{1} << Bgbit;
    static const inline double α = 0.0;
};

void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

template <class P, std::uint32_t LogQ, std::size_t DigitBits>
std::size_t testSelectorDigits()
{
    using namespace TFHEpp;
    std::mt19937_64 random(20262119);
    std::size_t checked = 0;
    for (std::size_t trial = 0; trial < 2048; ++trial) {
        typename P::T value{};
        if constexpr (requires { value.limb; })
            for (auto &limb : value.limb) limb = random();
        else value = random();
        if (trial == 0) value = std::numeric_limits<typename P::T>::max();
        if (trial == 1) value = typename P::T{};
        if (trial >= 2 && trial < 2 + ckks_detail::torus_width_v<P>)
            value = typename P::T{1} << (trial - 2);
        const auto reduced = ckks_detail::reduceToLevel<P, LogQ>(value);
        for (std::size_t digit = 0; digit < (LogQ + DigitBits - 1) / DigitBits; ++digit) {
            const auto expected = static_cast<std::uint64_t>(reduced >> (digit * DigitBits)) &
                ((std::uint64_t{1} << DigitBits) - 1);
            require(narrow_detail::selectorDigit<P, LogQ, DigitBits>(value, digit) == expected,
                    "native selector digit differs from complete torus reduction/shift");
            ++checked;
        }
        require(narrow_detail::selectorDigit<P, LogQ, DigitBits>(value,
                    (LogQ + DigitBits - 1) / DigitBits) == 0,
                "selector digit beyond ciphertext modulus is nonzero");
    }
    return checked;
}

void testSelectorDigitExtraction()
{
    struct Wide { using T = TFHEpp::MultiLimbUInt<9>; };
    struct Scalar { using T = __uint128_t; };
    const auto checked = testSelectorDigits<Wide, 518, 56>() +
        testSelectorDigits<Wide, 518, 48>() + testSelectorDigits<Wide, 518, 16>() +
        testSelectorDigits<Wide, 576, 56>() + testSelectorDigits<Wide, 64, 48>() +
        testSelectorDigits<Scalar, 61, 16>() + testSelectorDigits<Scalar, 64, 56>();
    std::cout << "{\"selector_digits_checked\":" << checked
              << ",\"all_digits_exact\":true}" << std::endl;
}

template <std::uint32_t Nbit>
void testDigitFFT()
{
    using namespace TFHEpp;
    using P = TestParam<Nbit>;
    std::array<std::int64_t, P::n> a, b;
    Polynomial<P> lhs, rhs, product;
    for (std::size_t i = 0; i < P::n; i++) {
        a[i] = static_cast<std::int64_t>((1031 * i + 53) % 65536) - 32768;
        b[i] = static_cast<std::int64_t>((3137 * i + 97) % 65536) - 32768;
        lhs[i] = P::T::from_signed_i64(a[i]);
        rhs[i] = P::T::from_signed_i64(b[i]);
    }
    PolyMulDigit<P>(product, lhs, rhs);
    for (std::size_t i = 0; i < P::n; i++) {
        std::int64_t expected = 0;
        for (std::size_t j = 0; j <= i; j++) expected += a[j] * b[i - j];
        for (std::size_t j = i + 1; j < P::n; j++) expected -= a[j] * b[P::n + i - j];
        require(product[i] == P::T::from_signed_i64(expected),
                "ring-sized DD FFT differs from integer convolution");
    }
}

template <std::uint32_t Nbit>
struct ActiveDigitRingParam : TestParam<Nbit> {
    using T = typename TestParam<Nbit>::T;
    static constexpr std::uint32_t B̅gbit = 4, B̅gₐbit = 4;
    static constexpr std::uint32_t l̅ = std::numeric_limits<T>::digits / B̅gbit;
    static constexpr std::uint32_t l̅ₐ = l̅;
};
using ActiveDigitParam = ActiveDigitRingParam<9>;
template <std::uint32_t Nbit>
struct GaussianActiveDigitRingParam : ActiveDigitRingParam<Nbit> {
    using T = typename ActiveDigitRingParam<Nbit>::T;
    static const inline double α = std::ldexp(1.0, -std::numeric_limits<T>::digits);
    static constexpr double ckks_min_noise_stddev = 3.2;
};
using GaussianActiveDigitParam = GaussianActiveDigitRingParam<9>;

template <std::uint32_t Nbit, std::size_t HammingWeight>
void testFrequencyBasis()
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitRingParam<Nbit>;
    using Base = CKKSNarrowBootstrapSchedule<P, 20, 440, 40, 34, HammingWeight, true, 16, true>;
    using Schedule = CKKSNarrowFrequencyBasisSchedule<Base, 32>;
    Key<P> source{}, evaluation{};
    std::mt19937_64 engine(20262119);
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : evaluation) value = P::T::from_signed_i64(ternary(engine));
    for (std::size_t x = 0; x < Base::blocks; x++)
        source[x * Base::block_width + (3 * x + 1) % Base::block_size(x)] =
            P::T::from_signed_i64(x % 2 ? -1 : 1);
    const CKKSNoise noise{P::α, 0};
    auto blind_key = std::make_unique<CKKSNarrowBootstrapKey<Base>>();
    CKKSNarrowBootstrapKeyGen<Base>(*blind_key, source, evaluation, noise);
    if constexpr (use_multilimb_digit_fft_v<P>) CKKSNarrowUseStagedFFT(*blind_key, 7);
    else {
        blind_key->cache_mode = CKKSNarrowCacheMode::Streamed;
        CKKSNarrowPrepareKey(*blind_key);
    }
    auto basis_key = std::make_unique<CKKSNarrowFrequencyBasisKey<Schedule>>();
    CKKSNarrowFrequencyBasisKeyGen<Schedule>(*basis_key, evaluation, noise);
    std::stringstream archive_bytes;
    { cereal::PortableBinaryOutputArchive archive(archive_bytes); archive(*basis_key); }
    auto loaded = std::make_unique<CKKSNarrowFrequencyBasisKey<Schedule>>();
    { cereal::PortableBinaryInputArchive archive(archive_bytes); archive(*loaded); }
    require(loaded->product.data == basis_key->product.data &&
            loaded->conjugation.data == basis_key->conjugation.data, "basis key archive changed public keys");
    using Input = CKKSCiphertext<P, 20, 12>;
    std::array<Input, 2> inputs;
    std::array<std::array<double, P::n>, 2> messages;
    for (std::size_t j = 0; j < inputs.size(); j++) {
        for (std::size_t i = 0; i < P::n; i++)
            messages[j][i] = (int((7 * i + 3 * j) % 17) - 8) / 32.0;
        ckksEncrypt<P>(inputs[j], messages[j], source, noise);
    }
    for (const auto degree : {12U, 40U}) {
        const CKKSNarrowFrequencyBasis plan(degree, degree == 12
            ? std::vector<std::uint32_t>{1, 3, 5, 6}
            : std::vector<std::uint32_t>{1, 3, 4, 9, 11, 16, 17, 19, 20});
        std::vector<std::complex<double>> coefficients(degree + 1);
        coefficients[0] = 0.125;
        for (std::size_t k = 1; k <= degree; k++)
            coefficients[k] = {0.03 / (k * k), -0.02 / (k * (k + 1))};
        std::vector<typename Schedule::Output> got;
        CKKSNarrowEvalSeriesWithBasisBatch<Schedule, 12>(got, inputs, *blind_key, *basis_key,
                                                       plan, coefficients, 3);
        if constexpr (Nbit == 6) {
            std::vector<typename Schedule::Output> raw;
            CKKSNarrowEvalSeriesWithBasisBatch<Schedule, 12>(raw, inputs, *blind_key, *loaded,
                                                           plan, coefficients, 2, false);
            for (std::size_t j = 0; j < got.size(); j++)
                require(raw[j].ct == got[j].ct, "raw/FFT basis or tile size changed exact ciphertexts");
        }
        if constexpr (use_multilimb_digit_fft_v<P>) {
            if (degree == 12) {
                auto raw = std::make_unique<typename Schedule::Output>();
                CKKSNarrowEvalSeriesWithBasis<Schedule, 12>(*raw, inputs[0], *blind_key,
                    *loaded, plan, coefficients, 2, false);
                require(raw->ct == got[0].ct, "raw and FFT extra DD/conjugation changed exact basis output");
            }
        }
        double error = 0;
        for (std::size_t j = 0; j < got.size(); j++) {
            const auto phase = trlwePhase<P>(inputs[j].ct, source);
            CKKSSlotVector<P> decoded;
            ckksSlotDecrypt<P>(decoded, got[j], evaluation);
            for (std::size_t i = 0; i < decoded.size(); i++) {
                const double cycles = std::ldexp(
                    ckks_detail::levelToLongDouble<P, Base::input_log_q>(phase[i]),
                    -int(Base::input_log_q));
                double expected = coefficients[0].real();
                for (std::size_t k = 1; k <= degree; k++)
                    expected += 2 * (coefficients[k] * std::polar(1.0,
                        2 * std::acos(-1.0) * k * cycles)).real();
                error = std::max(error, std::abs(decoded[i] - expected));
            }
        }
        std::cout << "basis N=" << P::n << " h=" << HammingWeight << " degree=" << degree
                  << " roots=" << plan.frequencies().size() << " error=" << error << std::endl;
        require(error < std::ldexp(1.0, -23), "basis series failed unchanged 23-bit gate");
        auto retained_key = makeCKKSCompactDDRelinKey<P, Schedule::output_log_q>(evaluation, noise);
        auto square = std::make_unique<CKKSMultResult<P, Schedule::output_log_q, Schedule::log_delta,
            Schedule::output_log_q, Schedule::log_delta>>();
        if constexpr (use_multilimb_digit_fft_v<P>) {
            auto retained_fft = makeCKKSDDRelinKeyFFT(*retained_key);
            CKKSMultRelinThenRescale<P>(*square, got[0], got[0], *retained_fft);
        }
        else CKKSMultRelinThenRescale<P>(*square, got[0], got[0], *retained_key);
        CKKSSlotVector<P> decoded, squared;
        ckksSlotDecrypt<P>(decoded, got[0], evaluation);
        ckksSlotDecrypt<P>(squared, *square, evaluation);
        for (std::size_t i = 0; i < squared.size(); i++)
            require(std::abs(squared[i] - decoded[i] * decoded[i]) < std::ldexp(1.0, -23),
                    "basis output cannot support its retained DD square");
    }
    bool rejected = false;
    try { CKKSNarrowFrequencyBasis invalid(12, {1, 3, 5}); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "basis plan accepted an uncovered frequency");
    require(!blind_key->relin_fft && !blind_key->factor_conjugation_fft &&
            !blind_key->tail_rotations_fft && !blind_key->output_conjugation_fft,
            "basis evaluator retained staged blind-rotation caches");
}

void testGaussianBasisFBS(const std::string &report_path = {})
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitParam;
    using Narrow = CKKSNarrowBootstrapSchedule<P, 48, 440, 40, 34, 31, true, 16, true>;
    using Refresh = CKKSNarrowRefreshSchedule<Narrow, 44, 37, 0, true>;
    using Basis = CKKSNarrowBasisRefreshSchedule<Refresh>;
    static_assert(Basis::output_log_q == 166 && Basis::log_delta == 72);
    static_assert(Basis::output_log_q - 2 * Basis::log_delta > 8);
    Key<P> sparse{}, dense{};
    std::mt19937_64 engine(20262119);
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : dense) value = P::T::from_signed_i64(ternary(engine));
    CKKSNarrowStructuredKeyGen<Narrow>(sparse);
    const CKKSNoise noise{P::α, 0};
    auto key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    CKKSNarrowRefreshKeyGen<Refresh>(*key, dense, sparse, dense, noise);
    CKKSNarrowUseStagedFFT(key->narrow, 7);
    auto basis_key = std::make_unique<CKKSNarrowFrequencyBasisKey<Basis>>();
    CKKSNarrowFrequencyBasisKeyGen<Basis>(*basis_key, dense, noise);
    CKKSSlotVector<P> slots;
    std::uniform_real_distribution<double> message(-4, 4);
    for (auto &slot : slots) slot = {message(engine), message(engine)};
    slots[0] = {-4, 4}; slots[1] = {4, -4}; slots[2] = {4, 4}; slots[3] = {0, 0};
    auto input = std::make_unique<typename Refresh::Input>();
    ckksSlotEncrypt<P>(*input, slots, dense, noise);
    auto retained = makeCKKSCompactDDRelinKey<P, Basis::output_log_q>(dense, noise);
    auto retained_fft = makeCKKSDDRelinKeyFFT(*retained);
    std::ostringstream report;
    report.precision(17);
    report << "{\n  \"n\": 512, \"h\": 31, \"gaussian_sigma\": 3.2,\n"
           << "  \"source_log_q\": 48, \"input_log_q\": 85, \"input_log_delta\": 44,\n"
           << "  \"selector_log_q\": 440, \"root_plain_log_delta\": 34, \"root_log_delta\": 40,\n"
           << "  \"coefficient_bits\": 32, \"output_log_q\": 166, \"output_log_delta\": 72,\n"
           << "  \"message_bound\": 4, \"target_bits\": 23, \"security_validated\": false,\n"
           << "  \"full_ring_validated\": false, \"results\": [\n";
    bool first = true;
    for (const std::string function : {"sigmoid", "gelu", "tanh"}) {
        const auto target = [&](double x) {
            if (function == "sigmoid") return 1 / (1 + std::exp(-x));
            if (function == "gelu") return 0.5 * x * (1 + std::erf(x / std::sqrt(2.0)));
            return std::tanh(x);
        };
        const auto curvature = [&](double x) {
            const double y = target(x);
            if (function == "sigmoid") return y * (1 - y) * (1 - 2 * y);
            if (function == "gelu")
                return (2 - x * x) * std::exp(-x * x / 2) / std::sqrt(2 * std::acos(-1.0));
            return -2 * y * (1 - y * y);
        };
        const bool windowed = function == "tanh";
        const std::uint32_t degree = windowed ? 40 : 12;
        const auto approximation = windowed ? CKKSBuildNarrowWindowedTrigonometricApproximation(
            target, [](double x) { const auto y = std::tanh(x); return 1 - y * y; },
            curvature, {1, 1, 2, 4, 40}, 4, Refresh::gap(), degree) :
            CKKSBuildNarrowTrigonometricApproximation(target, 4, Refresh::gap(), degree,
                function == "gelu" ? 2.0 : 1.0, 65536, 40, 2.0, curvature, 160.0);
        const CKKSNarrowFrequencyBasis plan(degree, windowed
            ? std::vector<std::uint32_t>{1, 3, 4, 9, 11, 16, 17, 19, 20}
            : std::vector<std::uint32_t>{1, 3, 5, 6});
        const auto start = std::chrono::steady_clock::now();
        auto output = std::make_unique<typename Basis::Output>();
        CKKSNarrowFunctionalRefreshWithBasis<Refresh>(*output, *input, *key,
            *basis_key, plan, approximation, true, 3);
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        CKKSSlotVector<P> decoded, expected, approximated;
        ckksSlotDecrypt<P>(decoded, *output, dense);
        double error = 0, evaluation_error = 0, amplitude = 0, approximation_amplitude = 0;
        for (std::size_t i = 0; i < decoded.size(); i++) {
            expected[i] = {target(slots[i].real()), target(slots[i].imag())};
            approximated[i] = {CKKSEvaluateNarrowApproximation(approximation, slots[i].real()),
                CKKSEvaluateNarrowApproximation(approximation, slots[i].imag())};
            error = std::max(error, std::abs(decoded[i] - expected[i]));
            evaluation_error = std::max(evaluation_error, std::abs(decoded[i] - approximated[i]));
            amplitude = std::max(amplitude, std::abs(expected[i]));
            approximation_amplitude = std::max(approximation_amplitude, std::abs(approximated[i]));
        }
        auto square = std::make_unique<CKKSMultResult<P, Basis::output_log_q, Basis::log_delta,
            Basis::output_log_q, Basis::log_delta>>();
        CKKSMultRelinThenRescale<P>(*square, *output, *output, *retained_fft);
        ckksSlotDecrypt<P>(decoded, *square, dense);
        double product_error = 0, product_evaluation_error = 0;
        for (std::size_t i = 0; i < decoded.size(); i++) {
            product_error = std::max(product_error, std::abs(decoded[i] - expected[i] * expected[i]));
            product_evaluation_error = std::max(product_evaluation_error,
                std::abs(decoded[i] - approximated[i] * approximated[i]));
        }
        const double tolerance = std::ldexp(1.0, -23);
        const double approximation_bound = std::sqrt(2.0) * *approximation.uniform_error_bound;
        const double square_tolerance = 2 * amplitude * tolerance + tolerance * tolerance + tolerance;
        const double evaluation_square_tolerance =
            2 * approximation_amplitude * tolerance + tolerance * tolerance + tolerance;
        const bool passed = error < tolerance && evaluation_error < tolerance && approximation_bound < tolerance &&
            product_error < square_tolerance && product_evaluation_error < evaluation_square_tolerance;
        if (!first) report << ",\n";
        first = false;
        report << "    {\"function\": \"" << function << "\", \"degree\": " << degree
               << ", \"basis_size\": " << plan.frequencies().size() << ", \"target_tolerance\": " << tolerance
               << ", \"approximation_uniform_bound\": " << approximation_bound
               << ", \"max_error\": " << error << ", \"evaluation_max_error\": " << evaluation_error
               << ", \"post_refresh_product_max_error\": " << product_error
               << ", \"post_refresh_product_evaluation_max_error\": " << product_evaluation_error
               << ", \"product_tolerance\": " << square_tolerance
               << ", \"product_evaluation_tolerance\": " << evaluation_square_tolerance
               << ", \"refresh_seconds\": " << elapsed << ", \"accuracy_passed\": " << (passed ? "true" : "false") << "}";
        std::cout << "Gaussian basis FBS " << function << " degree=" << degree << " roots="
                  << plan.frequencies().size() << " target=" << error << " evaluation=" << evaluation_error
                  << " square=" << product_error << " passed=" << passed << std::endl;
        require(passed, "Gaussian nonlinear basis FBS failed unchanged precision/amplitude gates");
    }
    report << "\n  ]\n}\n";
    if (report_path.empty()) std::cout << report.str();
    else {
        std::ofstream file(report_path);
        if (!file) throw std::runtime_error("cannot write Gaussian basis FBS report");
        file << report.str();
    }
}

template <std::uint32_t Q, std::uint32_t Bbar = 4>
void testCompactDD()
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitParam;
    using Compact = CKKSCompactDDRelinKey<P, Q, 4, Bbar>;
    static_assert(sizeof(Compact) * sizeof(P::T) == sizeof(typename Compact::Raw) * sizeof(typename Compact::Digit));
    auto secret = std::make_unique<Key<P>>();
    std::mt19937_64 engine(20262119 + Q);
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : *secret) value = P::T::from_signed_i64(ternary(engine));
    const CKKSNoise noise{P::α, 0};
    auto raw = makeCKKSDDRelinKey<P, Q, 4, Bbar>(*secret, noise);
    auto compact = makeCKKSCompactDDRelinKey(*raw);
    auto expanded = expandCKKSCompactDDRelinKey(*compact);
    require(raw->data == expanded->data, "compact DD failed its exact coefficient round trip");
    expanded.reset();
    auto raw_fft = makeCKKSDDRelinKeyFFT(*raw);
    auto compact_fft = makeCKKSDDRelinKeyFFT(*compact);
    require(raw_fft->data == compact_fft->data, "compact DD changed a prepared key spectrum");
    auto secondary16 = makeCKKSCompactDDRelinKeyWithSecondaryBase<16>(*compact);
    auto restored = makeCKKSCompactDDRelinKeyWithSecondaryBase<Bbar>(*secondary16);
    require(restored->data == compact->data, "secondary DD conversion changed encrypted primary rows");
    auto secondary16_fft = makeCKKSDDRelinKeyFFT(*secondary16);
    auto input = std::make_unique<Polynomial<P>>();
    auto reference = std::make_unique<TRLWE<P>>();
    auto actual = std::make_unique<TRLWE<P>>();
    for (std::size_t sample = 0; sample < 3; sample++) {
        for (auto &value : *input) {
            for (auto &limb : value.limb) limb = engine();
            value = ckks_detail::reduceToLevel<P, Q>(value);
        }
        CKKSRelinKeySwitch<P>(*reference, *input, *raw);
        CKKSRelinKeySwitch<P>(*actual, *input, *compact);
        require(*reference == *actual, "compact DD switching changed exact coefficients");
        CKKSRelinKeySwitch<P>(*actual, *input, *compact_fft);
        require(*reference == *actual, "compact prepared DD switching changed coefficients");
        CKKSRelinKeySwitch<P>(*actual, *input, *secondary16);
        require(*reference == *actual, "secondary DD conversion changed scalar switching");
        CKKSRelinKeySwitch<P>(*actual, *input, *secondary16_fft);
        require(*reference == *actual, "secondary DD conversion changed FFT switching");
    }
    std::stringstream payload;
    { cereal::PortableBinaryOutputArchive archive(payload); archive(*compact); }
    auto loaded = std::make_unique<Compact>();
    { cereal::PortableBinaryInputArchive archive(payload); archive(*loaded); }
    require(loaded->data == compact->data, "compact DD archive changed digits");
    std::stringstream invalid;
    { cereal::PortableBinaryOutputArchive archive(invalid); std::array<std::uint64_t, 7> metadata{}; archive(metadata); }
    bool rejected_metadata = false;
    try { cereal::PortableBinaryInputArchive archive(invalid); archive(*loaded); }
    catch (const std::runtime_error &) { rejected_metadata = true; }
    require(rejected_metadata, "compact DD accepted foreign archive metadata");
    if constexpr (Bbar < sizeof(typename Compact::Digit) * 8) {
        const auto saved_digit = compact->data[0][0][0][0];
        compact->data[0][0][0][0] = std::int64_t{1} << (Bbar - 1);
        bool rejected_digit = false;
        try { makeCKKSCompactDDRelinKeyWithSecondaryBase<16>(*compact); }
        catch (const std::invalid_argument &) { rejected_digit = true; }
        compact->data[0][0][0][0] = saved_digit;
        require(rejected_digit, "secondary DD conversion accepted an invalid source digit");
    }
    constexpr auto half = std::int64_t{1} << (Bbar - 1);
    require(Compact::pack(P::T::from_signed_i64(-half)) == -half &&
                Compact::pack(P::T::from_signed_i64(half - 1)) == half - 1,
            "compact DD lost a signed digit boundary");
    for (const auto value : {P::T::from_signed_i64(half), P::T{1} << 64}) {
        bool rejected = false;
        try { Compact::pack(value); }
        catch (const std::invalid_argument &) { rejected = true; }
        require(rejected, "compact DD accepted a non-digit coefficient");
    }
    if constexpr (Q > 40) {
        auto direct = makeCKKSCompactDDRelinKey<P, Q, 4, Bbar>(*secret, noise);
        auto lhs = std::make_unique<CKKSCiphertext<P, Q, 20>>();
        auto rhs = std::make_unique<CKKSCiphertext<P, Q, 20>>();
        CKKSSlotVector<P> a, b, before_a, before_b, got;
        for (std::size_t i = 0; i < a.size(); i++) {
            a[i] = {0.125 * std::sin(double(i)), 0.125 * std::cos(double(3 * i))};
            b[i] = {0.125 * std::cos(double(5 * i)), -0.125 * std::sin(double(2 * i))};
        }
        ckksSlotEncrypt<P>(*lhs, a, *secret, noise);
        ckksSlotEncrypt<P>(*rhs, b, *secret, noise);
        ckksSlotDecrypt<P>(before_a, *lhs, *secret);
        ckksSlotDecrypt<P>(before_b, *rhs, *secret);
        using Result = CKKSMultResult<P, Q, 20, Q, 20>;
        auto result = std::make_unique<Result>();
        auto reference_product = std::make_unique<Result>();
        CKKSMultRelinThenRescale<P>(*reference_product, *lhs, *rhs, *raw);
        CKKSMultRelinThenRescale<P>(*result, *lhs, *rhs, *compact);
        require(result->ct == reference_product->ct, "compact DD changed a late-rescale product");
        CKKSMultRelinThenRescale<P>(*result, *lhs, *rhs, *direct);
        ckksSlotDecrypt<P>(got, *result, *secret);
        for (std::size_t i = 0; i < got.size(); i++)
            require(std::abs(got[i] - before_a[i] * before_b[i]) < 0.01,
                    "direct compact DD generation has excessive multiplication error");
    }
    std::cout << "compact DD Q=" << Q << " raw/compact bytes=" << sizeof(*raw)
              << '/' << sizeof(*compact) << " Bbar=" << Bbar
              << " exact codec/switch/FFT/archive checks passed\n";
}

void testCompactS2C()
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitParam;
    using Compact = CKKSCompactRealLinearTransformPlan<P, 85, 44, 37>;
    using Raw = CKKSRealLinearTransformPlan<P, 85, 44, 37>;
    auto raw = std::make_unique<Raw>();
    auto compact = std::make_unique<Compact>();
    Key<P> secret{};
    std::mt19937_64 engine(20262119);
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &v : secret) v = P::T::from_signed_i64(ternary(engine));
    auto input_gk = std::make_unique<CKKSGaloisKey<P, 85>>();
    auto output_gk = std::make_unique<CKKSGaloisKey<P, 48>>();
    const CKKSNoise noise{P::α, 0};
    CKKSGaloisKeyGen<P, 85>(*input_gk, secret, noise);
    CKKSGaloisKeyGen<P, 48>(*output_gk, secret, noise);
    CKKSCiphertext<P, 85, 44> ct;
    CKKSSlotVector<P> slots;
    for (std::size_t i = 0; i < slots.size(); i++)
        slots[i] = {0.1 * std::sin(double(i)), 0.1 * std::cos(double(3 * i))};
    ckksSlotEncrypt<P>(ct, slots, secret, noise);
    auto expected = std::make_unique<CKKSPlainMulResult<P, 85, 44, 37>>();
    auto got = std::make_unique<CKKSPlainMulResult<P, 85, 44, 37>>();
    for (const int step : {16, 19}) {
        CKKSBuildPackedSlotToCoeffPlan<P, 85, 44, 37>(*raw, step);
        CKKSBuildPackedSlotToCoeffPlan<P, 85, 44, 37>(*compact, step);
        auto check = [&](const auto &before, const auto &after) {
            require(before.k_step == after.k_step && before.groups.size() == after.groups.size(), "compact S2C changed group order");
            for (std::size_t g = 0; g < before.groups.size(); g++) {
                require(before.groups[g].giant_step == after.groups[g].giant_step && before.groups[g].terms.size() == after.groups[g].terms.size(), "compact S2C changed giant step");
                for (std::size_t t = 0; t < before.groups[g].terms.size(); t++) {
                    require(before.groups[g].terms[t].baby_step == after.groups[g].terms[t].baby_step, "compact S2C changed baby step");
                    for (std::size_t n = 0; n < P::n; n++)
                        require(before.groups[g].terms[t].plain.poly[n] == ckks_detail::reduceToLevel<P, 85>(P::T::from_signed_i64(after.groups[g].terms[t].coefficients[n])), "compact S2C changed encoded coefficient");
                }
            }
        };
        check(raw->direct, compact->direct);
        check(raw->conjugate, compact->conjugate);
        CKKSRealLinearTransform<P, 85, 44, 37>(*expected, ct, *raw, *input_gk, *output_gk);
        CKKSRealLinearTransform<P, 85, 44, 37>(*got, ct, *compact, *input_gk, *output_gk);
        require(expected->ct == got->ct, "compact S2C changed evaluated ciphertext");
        std::stringstream payload;
        { cereal::PortableBinaryOutputArchive archive(payload); archive(compact->direct, compact->conjugate); }
        Compact loaded;
        loaded.has_direct = loaded.has_conjugate = true;
        { cereal::PortableBinaryInputArchive archive(payload); archive(loaded.direct, loaded.conjugate); }
        CKKSRealLinearTransform<P, 85, 44, 37>(*got, ct, loaded, *input_gk, *output_gk);
        require(expected->ct == got->ct, "compact S2C archive changed evaluation");
    }
    using Plan = CKKSCompactLinearTransformPlan<P, 85, 44, 37>;
    require(Plan::pack(P::T::from_signed_i64(INT64_MIN)) == INT64_MIN && Plan::pack(P::T::from_signed_i64(INT64_MAX)) == INT64_MAX, "compact S2C lost int64 boundaries");
    bool rejected = false;
    try { Plan::pack(P::T{1} << 64); } catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "compact S2C truncated a wide coefficient");
    using Small = CKKSCompactLinearTransformPlan<P, 17, 12, 8>;
    require(Small::pack(ckks_detail::reduceToLevel<P, 17>(P::T::from_signed_i64(-13))) == -13, "compact S2C lost active-level negative sign");
    std::stringstream invalid;
    { cereal::PortableBinaryOutputArchive archive(invalid); std::array<std::uint64_t, 7> metadata{}; archive(metadata); }
    rejected = false;
    try { cereal::PortableBinaryInputArchive archive(invalid); archive(compact->direct); } catch (const std::runtime_error &) { rejected = true; }
    require(rejected, "compact S2C accepted foreign metadata");
    std::cout << "compact S2C exact coefficients, evaluation, archive, range checks passed\n";
}

template <std::uint32_t PlainDelta, std::size_t DigitBits = 0>
void testSelectorChunks()
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitParam;
    using Schedule = CKKSNarrowBootstrapSchedule<P, 20, 400, 40, PlainDelta, 31>;
    using Cache = CKKSNarrowNTTSelector<Schedule, DigitBits>;
    static_assert(Cache::prime_count == (DigitBits == 0 && PlainDelta == 18 ? 1 : 2));
    Key<P> secret{};
    std::mt19937_64 engine(20262119 + PlainDelta);
    for (auto &v : secret) v = P::T::from_signed_i64(static_cast<int>(engine() % 3) - 1);
    auto key = std::make_unique<CKKSNarrowBootstrapKey<Schedule>>();
    CKKSAutoKeyGen<P, Schedule::factor_key_log_q>(key->factor_galois[P::nbit],
        2 * P::n - 1, secret, {P::α, 0});
    key->factor_conjugation_fft = makeCKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q>(key->factor_galois[P::nbit]);
    // Valid Gaussian ciphertexts exercise all active selector decomposition
    // digits. File tests below additionally use the actual signed indicators.
    std::vector<typename Schedule::Selector> selectors(Schedule::block_width);
    CKKSSlotVector<P> slots;
    for (std::size_t y = 0; y < selectors.size(); y++) {
        for (std::size_t i = 0; i < slots.size(); i++)
            slots[i] = {0.125 * std::sin(double(i + y)), 0.125 * std::cos(double(3 * i + y))};
        ckksSlotEncrypt<P>(selectors[y], slots, secret, {P::α, 0});
    }
    std::vector<Cache> prepared;
    narrow_detail::prepareSelectors<Schedule, DigitBits>(prepared, selectors);
    std::vector<CKKSSlotVector<P>> tables(3);
    for (std::size_t column = 0; column < tables.size(); column++)
        for (std::size_t i = 0; i < slots.size(); i++)
            tables[column][i] = std::polar(1.0, 0.03 * double((column + 1) * (i + 1)));
    for (const std::size_t columns : {1U, 3U}) {
        const auto used = std::span<const CKKSSlotVector<P>>(tables).first(columns);
        const auto baseline = narrow_detail::packedFactorsBatch<Schedule, DigitBits>(used, selectors, prepared, 0, *key);
        for (std::size_t i = 0; i < columns; i++)
            require(baseline[i].ct == narrow_detail::packedFactorReference<Schedule>(
                used[i], selectors, 0, *key).ct,
                "NTT factor recomposition differs from coefficient reference");
        for (const std::size_t chunk : {1U, 3U, 9U, 17U, 32U}) {
            const auto got = narrow_detail::packedFactorsBatch<Schedule, DigitBits>(used, selectors, {}, 0, *key, chunk);
            for (std::size_t i = 0; i < columns; i++)
                require(got[i].ct == baseline[i].ct, "selector chunks changed exact factor ciphertexts");
        }
    }
    bool rejected = false;
    try { CKKSNarrowUseSelectorChunks<Schedule>(*key, 0); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "selector chunk policy accepted zero size");
    std::cout << "selector chunks exact, primes=" << Cache::prime_count << " full group=" << Schedule::block_width << '\n';
}

template <std::uint32_t AuxiliaryScale>
void testBoundedRotations()
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitParam;
    constexpr auto q = 85;
    constexpr auto key_q = q + AuxiliaryScale;
    Key<P> secret{};
    std::mt19937_64 engine(20262119 + AuxiliaryScale);
    for (auto &v : secret) v = P::T::from_signed_i64(static_cast<int>(engine() % 3) - 1);
    auto raw = std::make_unique<CKKSGaloisKey<P, key_q>>();
    CKKSGaloisKeyGen<P, key_q>(*raw, secret, {P::α, 0});
    const std::vector<int> steps{0, 1, 32, 224, -1, 2 * int(P::n / 2) + 87, -int(P::n / 2)};
    std::vector<TRLWE<P>> input(steps.size()), expected(steps.size()), got;
    for (auto &ct : input)
        for (auto &component : ct)
            for (auto &v : component)
                for (auto &limb : v.limb) limb = engine();
    std::uint32_t mask = 0;
    for (const int step : steps) mask |= ((step % int(P::n / 2)) + int(P::n / 2)) % int(P::n / 2);
    auto full = makeCKKSPreparedRotationKeyFFT<P, key_q>(*raw, mask);
    for (std::size_t i = 0; i < steps.size(); i++) {
        if constexpr (AuxiliaryScale != 0)
            ckks_detail::applyWithAuxiliaryScale<P, q, AuxiliaryScale>(expected[i], input[i],
                [&](auto &out, const auto &in) { CKKSRotateSlots<P, key_q>(out, in, steps[i], *full); });
        else CKKSRotateSlots<P, q>(expected[i], input[i], steps[i], *full);
    }
    got = input;
    const auto peak = CKKSRotateSlotsBatchWithBoundedFFT<P, q, AuxiliaryScale>(steps,
        [&](std::size_t i) -> auto & { return got[i]; }, *raw);
    require(got == expected, "bounded FFT rotation changed exact ciphertext arithmetic");
    require(peak == CKKSPreparedAutoKeyFFT<P, key_q>::spectrum_count * sizeof(PolynomialInFD<P>) && peak < full->byteSize(), "bounded FFT rotation retained more than one key");
    const auto unchanged = got;
    require(CKKSRotateSlotsBatchWithBoundedFFT<P, q, AuxiliaryScale>({},
        [&](std::size_t i) -> auto & { return got[i]; }, *raw) == 0 && got == unchanged, "empty bounded rotation batch changed ciphertexts");
    const std::array<int, 1> zero{0};
    got = input;
    require(CKKSRotateSlotsBatchWithBoundedFFT<P, q, AuxiliaryScale>(zero,
        [&](std::size_t i) -> auto & { return got[i]; }, *raw) == 0 && got[0] == expected[0], "zero bounded rotation allocated a key or changed rounding");
    std::cout << "bounded rotations auxiliary=" << AuxiliaryScale << " exact, peak FFT bytes=" << peak << " full bytes=" << full->byteSize() << '\n';
}

void testFileSelectors()
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitRingParam<9>;
    using Schedule = CKKSNarrowBootstrapSchedule<P, 20, 400, 40, 40, 7, true, 16>;
    const auto suffix = std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    const auto directory = std::filesystem::temp_directory_path() / ("narrow-selector-test-" + suffix);
    const auto generated_directory = std::filesystem::temp_directory_path() / ("narrow-selector-generated-" + suffix);
    Key<P> source{}, evaluation{};
    std::mt19937_64 engine(20262119);
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : evaluation) value = P::T::from_signed_i64(ternary(engine));
    for (std::size_t x = 0; x < Schedule::blocks; x++)
        source[x * Schedule::block_width + Schedule::block_size(x) - 1] =
            P::T::from_signed_i64(x % 2 ? -1 : 1);
    auto key = std::make_unique<CKKSNarrowBootstrapKey<Schedule>>();
    const CKKSNoise noise{P::α, 0};
    CKKSNarrowBootstrapKeyGen<Schedule>(*key, source, evaluation, noise);
    using Input = CKKSCiphertext<P, Schedule::input_log_q, 12>;
    Input input;
    CKKSSlotVector<P> slots;
    for (std::size_t i = 0; i < slots.size(); i++)
        slots[i] = {0.125 * std::sin(double(i)), 0.125 * std::cos(double(3 * i))};
    ckksSlotEncrypt<P>(input, slots, source, noise);
    using Request = CKKSNarrowBlindRotationRequest<Schedule, 12>;
    const std::array<Request, 2> requests{{{&input, 1, {0.2, -0.1}}, {&input, 3, {-0.15, 0.25}}}};
    std::vector<typename Schedule::Output> baseline, got;
    CKKSNarrowBlindRotateBatch<Schedule, 12>(baseline, requests, *key, 2);
    const std::vector<typename Schedule::Selector> original_group(
        key->head.begin(), key->head.begin() + Schedule::block_width);
    const std::array<std::uint64_t, 9> file_header{
        0x4e4152524f573031ULL, P::n, Schedule::selector_log_q, Schedule::log_delta,
        Schedule::plain_log_delta, Schedule::block_width,
        std::numeric_limits<typename P::T>::digits, 0, 0};
    std::stringstream legacy_bytes;
    {
        cereal::PortableBinaryOutputArchive archive(legacy_bytes);
        archive(file_header);
        for (const auto &selector : original_group) archive(selector);
    }
    CKKSNarrowStoreSelectors<Schedule>(*key, directory);
    {
        std::ifstream file(directory / "head-0.bin", std::ios::binary);
        const std::string bulk_bytes((std::istreambuf_iterator<char>(file)), {});
        require(bulk_bytes == legacy_bytes.str(), "bulk selector writer changed legacy archive bytes");
    }
    // Exercise the old per-coefficient format with opposite byte order.
    {
        std::ofstream file(directory / "head-0.bin", std::ios::binary | std::ios::trunc);
        cereal::PortableBinaryOutputArchive archive(file,
            cereal::PortableBinaryOutputArchive::Options::BigEndian());
        archive(file_header);
        for (const auto &selector : original_group) archive(selector);
    }
    {
        std::vector<typename Schedule::Selector> loaded_group;
        key->selector_source->load(false, 0, loaded_group);
        for (std::size_t i = 0; i < original_group.size(); i++)
            require(loaded_group[i].ct == original_group[i].ct,
                    "bulk selector reader changed legacy big-endian ciphertexts");
    }
    require(key->head.capacity() == 0 && key->tail.capacity() == 0 && key->selector_source,
            "external selector cache retained expanded torus storage");
    require(key->selector_source->supportsRawChunkPrefetch(),
            "file selector source did not enable bounded raw read-ahead");
    for (const auto mode : {CKKSNarrowCacheMode::Streamed, CKKSNarrowCacheMode::Reference}) {
        key->cache_mode = mode;
        CKKSNarrowBlindRotateBatch<Schedule, 12>(got, requests, *key, 2);
        for (std::size_t i = 0; i < got.size(); i++)
            require(got[i].ct == baseline[i].ct, "file selector evaluation changed exact ciphertexts");
    }
    key->cache_mode = CKKSNarrowCacheMode::Prepared;
    CKKSNarrowPrepareAutomorphisms<Schedule>(*key);
    CKKSNarrowPrepareTailRotations<Schedule>(*key);
    CKKSNarrowPrepareRelinearization<Schedule>(*key);
    CKKSNarrowBlindRotateBatch<Schedule, 12>(got, requests, *key, 1);
    for (std::size_t i = 0; i < got.size(); i++)
        require(got[i].ct == baseline[i].ct, "prepared file selector evaluation changed ciphertexts");
    CKKSNarrowUseBoundedTailRotations<Schedule>(*key);
    require(!key->tail_rotations_fft && key->bounded_tail_rotation_fft, "bounded tail policy retained a full FFT cache");
    for (const auto tile : {1U, 2U}) {
        CKKSNarrowBlindRotateBatch<Schedule, 12>(got, requests, *key, tile);
        for (std::size_t i = 0; i < got.size(); i++)
            require(got[i].ct == baseline[i].ct, "bounded file selector evaluation changed ciphertexts");
    }
    for (const auto chunk : {1U, 7U, 73U, 100U}) {
        CKKSNarrowUseSelectorChunks<Schedule>(*key, chunk);
        for (const auto tile : {1U, 2U}) {
            CKKSNarrowBlindRotateBatch<Schedule, 12>(got, requests, *key, tile);
            for (std::size_t i = 0; i < got.size(); i++)
                require(got[i].ct == baseline[i].ct, "chunked file selector evaluation changed ciphertexts");
        }
    }
    const std::array<Input, 1> series_inputs{input};
    const std::vector<std::complex<double>> series{{0.1, 0}, {0.05, -0.025}, {-0.02, 0.03}};
    std::vector<typename Schedule::Output> series_baseline, series_got;
    CKKSNarrowEvalSeriesBatch<Schedule, 12>(series_baseline, series_inputs, *key, series, 2);
    CKKSNarrowUseStagedFFT<Schedule>(*key, 7);
    require(key->staged_fft_caches && !key->factor_conjugation_fft && !key->output_conjugation_fft &&
            !key->relin_fft && !key->tail_rotations_fft, "staged policy retained permanent FFT caches");
    std::stringstream public_before;
    { cereal::PortableBinaryOutputArchive archive(public_before); archive(*key); }
    for (const auto tile : {1U, 2U}) {
        CKKSNarrowBlindRotateBatch<Schedule, 12>(got, requests, *key, tile);
        for (std::size_t i = 0; i < got.size(); i++)
            require(got[i].ct == baseline[i].ct, "staged FFT evaluation changed ciphertexts");
        CKKSNarrowEvalSeriesBatch<Schedule, 12>(series_got, series_inputs, *key, series, tile);
        require(series_got[0].ct == series_baseline[0].ct, "staged FFT functional series changed ciphertexts");
    }
    std::stringstream public_after;
    { cereal::PortableBinaryOutputArchive archive(public_after); archive(*key); }
    require(public_before.str() == public_after.str() && !key->factor_conjugation_fft &&
            !key->output_conjugation_fft && !key->relin_fft && !key->tail_rotations_fft,
            "staged FFT evaluation mutated public setup or retained caches");

    std::stringstream payload;
    { cereal::PortableBinaryOutputArchive archive(payload); archive(*key); }
    auto loaded = std::make_unique<CKKSNarrowBootstrapKey<Schedule>>();
    loaded->selector_source = key->selector_source;
    loaded->bounded_tail_rotation_fft = true;
    loaded->selector_ntt_chunk_size = 7;
    loaded->staged_fft_caches = true;
    loaded->factor_observer = [](bool, std::size_t, const CKKSNarrowFactorTimings &) {};
    loaded->root_stage_observer = [](const CKKSNarrowRootTimings &) {};
    { cereal::PortableBinaryInputArchive archive(payload); archive(*loaded); }
    require(!loaded->selector_source, "key archive retained an external selector source");
    require(!loaded->bounded_tail_rotation_fft, "key archive retained a derived rotation policy");
    require(loaded->selector_ntt_chunk_size == 0, "key archive retained a derived selector chunk policy");
    require(!loaded->staged_fft_caches, "key archive retained a staged FFT policy");
    require(!loaded->factor_observer && !loaded->root_stage_observer,
            "key archive retained a caller-owned timing observer");

    auto generated = std::make_unique<CKKSNarrowBootstrapKey<Schedule>>();
    CKKSNarrowSetupKeyGen<Schedule>(*generated, evaluation, noise, CKKSNarrowCacheMode::Streamed);
    CKKSNarrowEncryptSelectorsToFile<Schedule>(*generated, source, evaluation, generated_directory, noise);
    CKKSNarrowBlindRotateBatch<Schedule, 12>(got, requests, *generated, 2);
    CKKSSlotVector<P> expected_slots, actual_slots;
    for (std::size_t i = 0; i < got.size(); i++) {
        ckksSlotDecrypt<P>(expected_slots, baseline[i], evaluation);
        ckksSlotDecrypt<P>(actual_slots, got[i], evaluation);
        for (std::size_t j = 0; j < actual_slots.size(); j++)
            require(std::abs(actual_slots[j] - expected_slots[j]) < 1e-7,
                    "direct-to-file selector generation has excessive error");
    }
    generated.reset();
    require(!std::filesystem::exists(generated_directory), "selector cache did not clean its owned directory");

    bool rejected_existing = false;
    try { CKKSNarrowSelectorFile<Schedule> duplicate(directory); }
    catch (const std::runtime_error &) { rejected_existing = true; }
    require(rejected_existing, "selector cache overwrote an existing directory");
    std::filesystem::resize_file(directory / "head-0.bin", 4);
    bool rejected_truncated = false;
    try {
        std::vector<typename Schedule::Selector> group;
        key->selector_source->load(false, 0, group);
    }
    catch (const std::exception &) { rejected_truncated = true; }
    require(rejected_truncated, "selector cache accepted a truncated group file");
    bool rejected_prefetched_truncated = false;
    try {
        CKKSNarrowBlindRotateBatch<Schedule, 12>(got, requests, *key, 2);
    }
    catch (const std::exception &) { rejected_prefetched_truncated = true; }
    require(rejected_prefetched_truncated,
            "prefetched selector evaluation accepted a truncated group file");
    { std::ofstream unrelated(directory / "keep.txt"); unrelated << "keep"; }
    key->selector_source.reset();
    require(std::filesystem::exists(directory / "keep.txt"), "selector cleanup removed unrelated data");
    std::filesystem::remove(directory / "keep.txt");
    std::filesystem::remove(directory);
    std::cout << "file selector batch/reference/FFT/generation/archive/cleanup checks passed\n";
}
struct NarrowDDDigitParam : ActiveDigitParam {
    static constexpr std::uint32_t Bgbit = 4;
    static constexpr T Bg = T{1} << Bgbit;
};

template <std::uint32_t LogQ, std::uint32_t Nbit = 9>
void testPreparedAuto()
{
    using namespace TFHEpp;
    using P = ActiveDigitRingParam<Nbit>;
    auto secret = std::make_unique<Key<P>>();
    for (std::size_t i = 0; i < P::n; i++)
        (*secret)[i] = i % 3 == 0 ? -typename P::T{1} : typename P::T{i % 3 == 1 ? 1U : 0U};
    auto raw = std::make_unique<CKKSAutoKey<P, LogQ>>();
    CKKSAutoKeyGen<P, LogQ>(*raw, 2 * P::n - 1, *secret, {0, 3});
    // Include unused high bits in both key and input; decomposition carries
    // must match the original full-torus offsets at unaligned moduli.
    if constexpr (LogQ + 2 < std::numeric_limits<typename P::T>::digits)
        (*raw)[0][0][0][0] += typename P::T{1} << (LogQ + 2);
    const auto setup_start = std::chrono::steady_clock::now();
    auto prepared = makeCKKSPreparedAutoKeyFFT<P, LogQ>(*raw);
    const double setup_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - setup_start).count();
    auto input = std::make_unique<TRLWE<P>>();
    std::mt19937_64 engine(20262119 + LogQ);
    for (auto &component : *input)
        for (auto &value : component) {
            value = typename P::T{};
            for (auto &limb : value.limb) limb = engine();
        }
    auto reference = std::make_unique<TRLWE<P>>();
    auto actual = std::make_unique<TRLWE<P>>();
    const auto raw_start = std::chrono::steady_clock::now();
    CKKSEvalAuto<P, LogQ>(*reference, *input, 2 * P::n - 1, *raw);
    const double raw_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - raw_start).count();
    const auto prepared_start = std::chrono::steady_clock::now();
    CKKSEvalAuto<P, LogQ>(*actual, *input, 2 * P::n - 1, *prepared);
    const double prepared_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - prepared_start).count();
    require(*actual == *reference, "prepared auto switch changed exact full-torus arithmetic");
    auto secondary16 = makeCKKSPreparedAutoKeyFFT<P, LogQ, 16>(*raw);
    CKKSEvalAuto<P, LogQ>(*actual, *input, 2 * P::n - 1, *secondary16);
    require(*actual == *reference, "secondary16 auto switch changed exact full-torus arithmetic");
    CKKSEvalAutoStreamedFFT<P, LogQ, P::B̅gbit>(*actual, *input, 2 * P::n - 1, *raw);
    require(*actual == *reference, "streamed original-radix auto switch changed exact arithmetic");
    CKKSEvalAutoStreamedFFT<P, LogQ>(*actual, *input, 2 * P::n - 1, *raw);
    require(*actual == *reference, "streamed adaptive-radix auto switch changed exact arithmetic");
    std::cout << "prepared auto N=" << P::n << " Q=" << LogQ << " setup=" << setup_seconds
              << " raw=" << raw_seconds << " prepared=" << prepared_seconds
              << " cache_bytes=" << prepared->spectra.size() * sizeof(PolynomialInFD<P>) << std::endl;
}

template <bool LateRescale = false>
void testPreparedDDTree()
{
    using namespace TFHEpp;
    using P = NarrowDDDigitParam;
    using Schedule = CKKSNarrowBootstrapSchedule<P, 20, 189, 20, 20, 7, LateRescale>;
    auto secret = std::make_unique<Key<P>>();
    for (std::size_t i = 0; i < P::n; i++)
        (*secret)[i] = P::T::from_signed_i64(static_cast<int>(i % 3) - 1);
    auto raw = std::make_unique<typename Schedule::RelinKeys>();
    CKKSDDRelinKeyChainGen<P, Schedule::relin_chain_start_log_q, Schedule::log_delta,
        Schedule::tree_depth>(*raw, *secret, {0, 3});
    auto cache = makeCKKSPreparedDDRelinKeyChainFFT(*raw);
    require(cache->byteSize() > 0, "prepared DD chain has no spectra");
    std::vector<typename Schedule::Factor> factors(Schedule::blocks + 1);
    CKKSSlotVector<P> slots;
    for (std::size_t leaf = 0; leaf < factors.size(); leaf++) {
        for (std::size_t i = 0; i < slots.size(); i++)
            slots[i] = {0.25 + static_cast<double>((i + leaf) % 7) / 32,
                        static_cast<double>((i * 3 + leaf) % 5) / 64};
        ckksSlotEncrypt<P>(factors[leaf], slots, *secret, {0, 3});
    }
    auto reference = std::make_unique<typename Schedule::TreeResult>();
    auto actual = std::make_unique<typename Schedule::TreeResult>();
    narrow_detail::productTree<Schedule>(*reference, factors, *raw);
    narrow_detail::productTree<Schedule>(*actual, std::move(factors), *raw, cache.get());
    require(actual->ct == reference->ct,
            "prepared DD product tree changed the exact ciphertext");
    std::cout << "prepared DD three-level encrypted product tree matched exactly\n";
}

void testAuxiliaryScaleAuto()
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitParam;
    constexpr std::uint32_t q = 90, delta = 20, auxiliary = 26;
    auto secret = std::make_unique<Key<P>>();
    std::mt19937_64 engine(20262119);
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : *secret) value = P::T::from_signed_i64(ternary(engine));
    auto raw_key = std::make_unique<CKKSAutoKey<P, q + auxiliary>>();
    CKKSAutoKeyGen<P, q + auxiliary>(*raw_key, 2 * P::n - 1, *secret, {P::α, 0});
    const auto fft_key = makeCKKSPreparedAutoKeyFFT<P, q + auxiliary>(*raw_key);
    auto input = std::make_unique<CKKSCiphertext<P, q, delta>>();
    auto raw = std::make_unique<CKKSCiphertext<P, q, delta>>();
    auto fft = std::make_unique<CKKSCiphertext<P, q, delta>>();
    CKKSSlotVector<P> slots, before, got;
    for (std::size_t i = 0; i < slots.size(); i++)
        slots[i] = {0.125 * std::sin(double(i)), 0.125 * std::cos(double(3 * i))};
    ckksSlotEncrypt<P>(*input, slots, *secret, {P::α, 0});
    input->ct[0][0] += P::T{1} << (q + 1); // noncanonical upper bits are ignored
    ckksSlotDecrypt<P>(before, *input, *secret);
    CKKSEvalAutoWithAuxiliaryScale<P, q, delta, auxiliary>(*raw, *input, 2 * P::n - 1, *raw_key);
    CKKSEvalAutoWithAuxiliaryScale<P, q, delta, auxiliary>(*fft, *input, 2 * P::n - 1, *fft_key);
    require(raw->ct == fft->ct, "auxiliary scale raw/FFT switching differs");
    *fft = *input;
    CKKSEvalAutoWithAuxiliaryScale<P, q, delta, auxiliary>(*fft, *fft, 2 * P::n - 1, *fft_key);
    require(raw->ct == fft->ct, "auxiliary scale switching failed aliasing");
    ckksSlotDecrypt<P>(got, *raw, *secret);
    for (std::size_t i = 0; i < got.size(); i++)
        require(std::abs(got[i] - std::conj(before[i])) < 0.005,
                "auxiliary scale switch failed to preserve the input phase");
    std::cout << "auxiliary scale DD auto raw/FFT/alias phase checks passed\n";
}

void testRelinThenRescale()
{
    using namespace TFHEpp;
    using P = GaussianActiveDigitParam;
    auto secret = std::make_unique<Key<P>>();
    // The accuracy screen assumes an iid dense ternary evaluation key.
    // A fixed seed gives that distribution without a periodic Fourier spike.
    std::mt19937_64 secret_engine(20262119);
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : *secret) value = P::T::from_signed_i64(ternary(secret_engine));
    auto lhs = std::make_unique<CKKSCiphertext<P, 90, 20>>();
    auto rhs = std::make_unique<CKKSCiphertext<P, 85, 17>>();
    CKKSSlotVector<P> a, b, before_a, before_b, got;
    for (std::size_t i = 0; i < a.size(); i++) {
        a[i] = {0.125 * std::sin(double(i)), 0.125 * std::cos(double(3 * i))};
        b[i] = {0.125 * std::cos(double(5 * i)), -0.125 * std::sin(double(2 * i))};
    }
    for (const CKKSNoise noise : {CKKSNoise{0, 0}, CKKSNoise{P::α, 0}}) {
        ckksSlotEncrypt<P>(*lhs, a, *secret, noise);
        ckksSlotEncrypt<P>(*rhs, b, *secret, noise);
        ckksSlotDecrypt<P>(before_a, *lhs, *secret);
        ckksSlotDecrypt<P>(before_b, *rhs, *secret);
        const auto raw_key = makeCKKSDDRelinKey<P, 85, 4, 4>(*secret, noise);
        const auto fft_key = makeCKKSDDRelinKeyFFT(*raw_key);
        auto raw = std::make_unique<CKKSMultResult<P, 90, 20, 85, 17>>();
        auto fft = std::make_unique<CKKSMultResult<P, 90, 20, 85, 17>>();
        static_assert(CKKSMultTraits<P, 90, 20, 85, 17>::log_q == 65 &&
                      CKKSMultTraits<P, 90, 20, 85, 17>::log_delta == 17);
        CKKSMultRelinThenRescale<P>(*raw, *lhs, *rhs, *raw_key);
        CKKSMultRelinThenRescale<P>(*fft, *lhs, *rhs, *fft_key);
        require(raw->ct == fft->ct, "late rescale raw/FFT DD arithmetic differs");
        ckksSlotDecrypt<P>(got, *raw, *secret);
        double maximum_error = 0;
        for (std::size_t i = 0; i < got.size(); i++)
            maximum_error = std::max(maximum_error, std::abs(got[i] - before_a[i] * before_b[i]));
        std::cout << "late rescale mixed-scale max error=" << maximum_error << std::endl;
        // Output scale 2^17 includes final ciphertext rounding. This API
        // accuracy check is independent of the unchanged bootstrap gates.
        require(maximum_error < 0.01,
                "late rescale mixed-level/mixed-scale DD multiplication has excessive error");
    }
    std::cout << "DD relinearize-then-rescale mixed levels/scales passed\n";
}

void packedWorkspaceFixture(const std::string &directory, bool write)
{
    using namespace TFHEpp;
    using P = ActiveDigitParam;
    using Schedule = CKKSNarrowBootstrapSchedule<P, 20, 71, 8, 10, 31>;
    const auto root = std::filesystem::path(directory);
    auto packed = std::make_unique<CKKSNarrowPackedKey<Schedule>>();
    if (write) {
        require(std::filesystem::create_directory(root), "packed fixture directory already exists");
        std::filesystem::permissions(root, std::filesystem::perms::owner_all,
            std::filesystem::perm_options::replace);
        auto source = std::make_unique<Key<P>>();
        auto dense = std::make_unique<Key<P>>();
        CKKSNarrowStructuredKeyGen<Schedule>(*source);
        for (std::size_t i = 0; i < P::n; ++i)
            (*dense)[i] = P::T::from_signed_i64(static_cast<int>(i % 3) - 1);
        CKKSNarrowPackedKeyGen<Schedule>(*packed, *source, *dense, {0, 3});
        std::ofstream file(root / "packed.bin", std::ios::binary);
        { cereal::PortableBinaryOutputArchive archive(file); archive(*packed); }
        file.close(); require(bool(file), "cannot save packed fixture key");
    }
    else {
        std::ifstream file(root / "packed.bin", std::ios::binary);
        require(bool(file), "cannot open packed fixture key");
        cereal::PortableBinaryInputArchive archive(file); archive(*packed);
        require(file.peek() == std::char_traits<char>::eof(), "trailing packed fixture key data");
    }
    for (std::size_t levels : {0, 2}) {
        CKKSNarrowPreparePackedSubstitutions<Schedule>(*packed, levels);
        const auto path = root / ("terms-" + std::to_string(levels) + ".bin");
        std::fstream file(path, std::ios::binary | (write ? std::ios::out | std::ios::trunc : std::ios::in));
        require(bool(file), "cannot open packed fixture terms");
        const auto start = std::chrono::steady_clock::now();
        std::size_t terms = 0;
        if (write) {
            cereal::PortableBinaryOutputArchive archive(file);
            CKKSNarrowExtractSelectorTerms<Schedule>(*packed,
                [&](bool tail, std::size_t index, const auto &term) {
                    const auto position = static_cast<std::uint64_t>(index);
                    archive(tail, position, term); ++terms;
                });
        }
        else {
            cereal::PortableBinaryInputArchive archive(file);
            auto expected = std::make_unique<typename Schedule::Selector>();
            CKKSNarrowExtractSelectorTerms<Schedule>(*packed,
                [&](bool tail, std::size_t index, const auto &term) {
                    bool expected_tail{}; std::uint64_t position{};
                    archive(expected_tail, position, *expected);
                    require(tail == expected_tail && index == position && term.ct == expected->ct,
                        "packed extraction differs from saved pre-change terms");
                    ++terms;
                });
            require(file.peek() == std::char_traits<char>::eof(), "trailing packed fixture terms");
        }
        require(terms == 4 * P::n, "incomplete packed fixture terms");
        const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::cout << "packed fixture " << (write ? "written" : "matched") << " N=" << P::n
                  << " cached_levels=" << levels << " terms=" << terms << " seconds=" << elapsed << '\n';
    }
}

void testPreparedPackedSubstitutions()
{
    using namespace TFHEpp;
    using P = ActiveDigitParam;
    using Schedule = CKKSNarrowBootstrapSchedule<P, 20, 71, 8, 10, 31>;
    static_assert(CKKSNarrowPackedKey<Schedule>::log_q == 90);
    auto source = std::make_unique<Key<P>>();
    auto dense = std::make_unique<Key<P>>();
    CKKSNarrowStructuredKeyGen<Schedule>(*source);
    for (std::size_t i = 0; i < P::n; i++)
        (*dense)[i] = P::T::from_signed_i64(static_cast<int>(i % 3) - 1);
    auto packed = std::make_unique<CKKSNarrowPackedKey<Schedule>>();
    CKKSNarrowPackedKeyGen<Schedule>(*packed, *source, *dense, {0, 3});
    auto key = std::make_unique<CKKSNarrowBootstrapKey<Schedule>>();
    key->cache_mode = CKKSNarrowCacheMode::Streamed;
    CKKSNarrowReconstructSelectors<Schedule>(*key, *packed);
    const auto raw_head = key->head, raw_tail = key->tail;
    const auto directory = std::filesystem::temp_directory_path() /
        ("tfhepp-packed-selector-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    auto check_file_reconstruction = [&] {
        CKKSNarrowReconstructSelectorsToFile<Schedule>(*key, *packed, directory);
        require(key->head.empty() && key->tail.empty() && key->head_ntt.empty() && key->tail_ntt.empty(),
                "packed file reconstruction retained expanded selectors or caches");
        require(!std::filesystem::exists(directory / "packed-accumulator.bin"),
                "packed file reconstruction retained its temporary accumulator");
        std::vector<typename Schedule::Selector> group;
        for (bool tail : {false, true}) {
            const auto &expected = tail ? raw_tail : raw_head;
            for (std::size_t index = 0; index < (tail ? Schedule::tail_packs : Schedule::blocks); index++) {
                key->selector_source->load(tail, index, group);
                for (std::size_t y = 0; y < Schedule::block_width; y++)
                    require(group[y].ct == expected[index * Schedule::block_width + y].ct,
                            "packed disk reconstruction changed a selector");
            }
        }
        auto prior = key->selector_source;
        bool rejected_existing = false;
        try { CKKSNarrowReconstructSelectorsToFile<Schedule>(*key, *packed, directory); }
        catch (const std::runtime_error &) { rejected_existing = true; }
        require(rejected_existing && key->selector_source == prior,
                "packed reconstruction overwrote its existing directory or source");
        prior.reset();
        key->selector_source.reset();
        require(!std::filesystem::exists(directory), "packed disk reconstruction did not clean its files");
    };
    check_file_reconstruction();
    CKKSNarrowPreparePackedSubstitutions<Schedule>(*packed, 2);
    for (std::size_t r = 0; r < P::nbit; r++)
        require(bool(packed->substitutions_fft[r]) == (r >= P::nbit - 2),
                "packed substitution cache selected the wrong tree depths");
    CKKSNarrowReconstructSelectors<Schedule>(*key, *packed);
    for (std::size_t i = 0; i < raw_head.size(); i++)
        require(key->head[i].ct == raw_head[i].ct, "prepared packed reconstruction changed a head selector");
    for (std::size_t i = 0; i < raw_tail.size(); i++)
        require(key->tail[i].ct == raw_tail[i].ct, "prepared packed reconstruction changed a tail selector");
    check_file_reconstruction();
    bool rejected = false;
    try { CKKSNarrowPreparePackedSubstitutions<Schedule>(*packed, P::nbit + 1); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected && packed->substitutions_fft.back(),
            "invalid packed cache depth was accepted or destroyed the prior cache");
    CKKSNarrowPreparePackedSubstitutions<Schedule>(*packed, 0);
    for (const auto &cache : packed->substitutions_fft)
        require(!cache, "zero-depth packed preparation retained a cache");
    std::cout << "prepared packed extraction matched every h=31 selector exactly\n";
}

void testPreparedTailRotations()
{
    using namespace TFHEpp;
    using P = ActiveDigitParam;
    using Schedule = CKKSNarrowBootstrapSchedule<P, 20, 119, 10, 34, 31>;
    auto secret = std::make_unique<Key<P>>();
    for (std::size_t i = 0; i < P::n; i++)
        (*secret)[i] = P::T::from_signed_i64(static_cast<int>(i % 3) - 1);
    auto key = std::make_unique<CKKSNarrowBootstrapKey<Schedule>>();
    CKKSGaloisKeyGen<P, Schedule::factor_log_q>(key->factor_galois, *secret, {0, 3});
    CKKSNarrowPrepareTailRotations<Schedule>(*key);
    for (std::size_t bit = 0; bit < key->tail_rotations_fft->bits.size(); bit++)
        require(bool(key->tail_rotations_fft->bits[bit]) == (bit >= 5),
                "tail cache prepared an unused rotation bit or omitted a required bit");
    auto input = std::make_unique<TRLWE<P>>();
    std::mt19937_64 engine(20262119);
    for (auto &component : *input)
        for (auto &value : component)
            for (auto &limb : value.limb) limb = engine();
    auto raw = std::make_unique<TRLWE<P>>();
    auto prepared = std::make_unique<TRLWE<P>>();
    for (std::size_t ell = 0; ell < Schedule::tail_blocks_per_pack; ell++) {
        const int steps = -static_cast<int>(Schedule::slots - 2 * (ell + 1) * Schedule::tail_width);
        CKKSRotateSlots<P, Schedule::factor_log_q>(*raw, *input, steps, key->factor_galois);
        CKKSRotateSlots<P, Schedule::factor_log_q>(*prepared, *input, steps, *key->tail_rotations_fft);
        require(*raw == *prepared, "prepared tail rotation changed exact ciphertext arithmetic");
        *prepared = *input;
        CKKSRotateSlots<P, Schedule::factor_log_q>(*prepared, *prepared,
            steps + 2 * Schedule::slots, *key->tail_rotations_fft);
        require(*raw == *prepared, "prepared rotation failed aliasing or step normalization");
    }
    bool missing_rejected = false;
    try {
        CKKSRotateSlots<P, Schedule::factor_log_q>(*prepared, *input, 1, *key->tail_rotations_fft);
    }
    catch (const std::invalid_argument &) { missing_rejected = true; }
    require(missing_rejected, "partial rotation cache accepted a missing key");
    std::cout << "prepared tail rotation chain matched exactly, cache_bytes="
              << key->tail_rotations_fft->byteSize() << '\n';
}

template <std::uint32_t LogQ>
void testActiveDigitProduct()
{
    using namespace TFHEpp;
    using P = ActiveDigitParam;
    Polynomial<P> torus, digit, product, full;
    std::array<int, P::n> signed_digits;
    std::mt19937_64 engine(20262119 + LogQ);
    for (std::size_t i = 0; i < P::n; i++) {
        // Include nonzero bits above the active level and negative digits.
        for (auto &limb : torus[i].limb) limb = engine();
        signed_digits[i] = static_cast<int>(engine() % 16) - 8;
        digit[i] = P::T{signed_digits[i]};
    }
    PolyMulTorusByDigit<P, LogQ>(product, torus, digit);
    PolyMulTorusByDigit<P>(full, torus, digit);
    for (std::size_t i = 0; i < P::n; i++) {
        P::T expected{};
        for (std::size_t j = 0; j < P::n; j++) {
            const auto index = (P::n + i - j) % P::n;
            const auto term = torus[j] * signed_digits[index];
            if (j <= i) expected += term;
            else expected -= term;
        }
        expected = ckks_detail::reduceToLevel<P, LogQ>(expected);
        require(ckks_detail::reduceToLevel<P, LogQ>(product[i]) == expected,
                "active-level DD product differs from integer convolution");
        require(ckks_detail::reduceToLevel<P, LogQ>(full[i]) == expected,
                "full-width DD product differs from integer convolution");
    }
}

template <class Schedule>
void testTables()
{
    using namespace TFHEpp;
    using P = typename Schedule::Param;
    std::array<double, P::n> a;
    for (std::size_t i = 0; i < P::n; i++)
        a[i] = static_cast<double>((37 * i + 11) % 257) / 257.0 - 0.5;
    for (std::size_t x = 0; x < Schedule::blocks; x++) {
        const auto head = narrow_detail::headTable<Schedule>(a, x, 3);
        const auto tail = narrow_detail::tailTable<Schedule>(a, x / Schedule::tail_blocks_per_pack, 3);
        for (std::size_t y = 0; y < Schedule::block_width; y++) {
            CKKSSlotVector<P> head_rotated, tail_rotated;
            rotateCKKSSlotVector<P>(head_rotated, head, -static_cast<int>(y));
            rotateCKKSSlotVector<P>(tail_rotated, tail, -static_cast<int>(y));
            for (std::size_t i = 0; i < Schedule::slots; i++) {
                const auto j = x * Schedule::block_width + y;
                const auto expected = narrow_detail::root<P>(
                    a, static_cast<std::int64_t>(i) - j, j <= i ? 1 : -1, 3);
                const std::size_t ell = x % Schedule::tail_blocks_per_pack;
                const auto actual = i < Schedule::head_slots ? head_rotated[i] :
                    tail_rotated[(2 * ell + 1) * Schedule::tail_width + i - Schedule::head_slots];
                require(std::abs(actual - expected) < 1e-12, "narrow table does not match negacyclic convolution");
            }
        }
    }
    // Lemma 5.1, checked against coefficient automorphisms for both primes.
    const auto table = narrow_detail::headTable<Schedule>(a, Schedule::blocks - 1, 3);
    Polynomial<P> polynomial;
    ckksSlotEncode<P, Schedule::selector_log_q, Schedule::plain_log_delta>(polynomial, table);
    for (std::size_t prime = 0; prime < 2; prime++) {
        const auto modulus = modular_ntt::degree65536_primes[prime].value;
        const modular_ntt::NegacyclicNTTPlan plan(P::n, modular_ntt::degree65536_primes[prime]);
        auto to_ntt = [&](const Polynomial<P> &poly) {
            std::vector<std::uint64_t> result(P::n);
            for (std::size_t i = 0; i < P::n; i++) {
                const auto [negative, magnitude] =
                    ckks_detail::smallSignedMagnitude<P, Schedule::selector_log_q>(poly[i]);
                const auto reduced = magnitude % modulus;
                result[i] = negative && reduced != 0 ? modulus - reduced : reduced;
            }
            plan.forward(result);
            return result;
        };
        const auto base = to_ntt(polynomial);
        for (std::size_t y = 0; y < Schedule::block_width; y++) {
            const auto exponent = modular_ntt::power(5,
                (Schedule::slots - y) % Schedule::slots, 2 * P::n);
            Polynomial<P> rotated;
            Automorphism<P>(rotated, polynomial, static_cast<int>(exponent));
            const auto actual = to_ntt(rotated);
            for (std::size_t i = 0; i < P::n; i++) {
                const auto index = ((exponent * (2 * i + 1)) % (2 * P::n) - 1) / 2;
                require(actual[i] == base[index], "NTT slot permutation does not match automorphism");
            }
        }
    }
}

template <class Schedule>
void testEncrypted()
{
    using namespace TFHEpp;
    using P = typename Schedule::Param;
    Key<P> source{}, evaluation{};
    for (std::size_t x = 0; x < Schedule::blocks; x++) {
        const auto y = (3 * x + 1) % Schedule::block_size(x);
        source[x * Schedule::block_width + y] = x % 2 ? -typename P::T{1} : typename P::T{1};
    }
    for (std::size_t i = 0; i < P::n; i++)
        evaluation[i] = i % 3 == 0 ? -typename P::T{1} : typename P::T{i % 3 == 1 ? 1U : 0U};
    auto key = std::make_unique<CKKSNarrowBootstrapKey<Schedule>>();
    CKKSNarrowBootstrapKeyGen<Schedule>(*key, source, evaluation, {0, 3});
    require(key->head.size() + key->tail.size() == Schedule::selector_count,
            "incorrect complex-packed selector count");
    using Input = CKKSCiphertext<P, Schedule::input_log_q, 12>;
    Input input;
    std::array<double, P::n> coefficients;
    for (std::size_t i = 0; i < P::n; i++)
        coefficients[i] = static_cast<double>(static_cast<int>(i % 11) - 5) / 32;
    ckksEncrypt<P>(input, coefficients, source, {0, 3});
    const auto phase = trlwePhase<P>(input.ct, source);
    for (const auto frequency : {1U, 3U}) {
        typename Schedule::Output result;
        const std::complex<double> weight{0.7, -0.2};
        CKKSNarrowBlindRotate<Schedule>(result, input, *key, frequency, weight);
        CKKSSlotVector<P> got;
        ckksSlotDecrypt<P>(got, result, evaluation);
        double error = 0;
        for (std::size_t i = 0; i < Schedule::slots; i++) {
            const double cycles = static_cast<double>(std::ldexp(
                ckks_detail::levelToLongDouble<P, Schedule::input_log_q>(phase[i]),
                -static_cast<int>(Schedule::input_log_q)));
            const auto expected = weight * std::polar(1.0, 2 * std::acos(-1.0) * frequency * cycles);
            error = std::max(error, std::abs(got[i] - expected));
        }
        std::cout << "N=" << P::n << " h=" << Schedule::blocks
                  << " frequency=" << frequency << " error=" << error << '\n';
        require(error < 2e-4, "encrypted narrow rotation has excessive error");
    }
    Input different_input;
    for (auto &value : coefficients) value = -value + 0.07;
    ckksEncrypt<P>(different_input, coefficients, source, {0, 3});
    using Request = CKKSNarrowBlindRotationRequest<Schedule, 12>;
    const std::array<Request, 3> requests = {
        Request{&input, 1, {1, 0}}, Request{&different_input, 3, {0.4, 0.3}},
        Request{&input, 2, {-0.2, 0.1}}};
    std::vector<typename Schedule::Output> batched;
    CKKSNarrowBlindRotateBatch<Schedule, 12>(batched, requests, *key, 2);
    key->cache_mode = CKKSNarrowCacheMode::Streamed;
    CKKSNarrowPrepareKey<Schedule>(*key);
    require(key->head_ntt.empty() && key->tail_ntt.empty(),
            "streamed key retained the complete NTT cache");
    std::vector<typename Schedule::Output> streamed;
    CKKSNarrowBlindRotateBatch<Schedule, 12>(streamed, requests, *key, 2);
    for (std::size_t column = 0; column < requests.size(); column++)
        require(streamed[column].ct == batched[column].ct,
                "streamed NTT batch differs from the fully prepared evaluation");
    if constexpr (P::n == 16) {
        key->cache_mode = CKKSNarrowCacheMode::Reference;
        CKKSNarrowPrepareKey<Schedule>(*key);
        std::vector<typename Schedule::Output> reference;
        CKKSNarrowBlindRotateBatch<Schedule, 12>(reference, requests, *key, 2);
        for (std::size_t column = 0; column < requests.size(); column++) {
            CKKSSlotVector<P> cached_got, reference_got;
            ckksSlotDecrypt<P>(cached_got, batched[column], evaluation);
            ckksSlotDecrypt<P>(reference_got, reference[column], evaluation);
            for (std::size_t i = 0; i < Schedule::slots; i++)
                require(std::abs(reference_got[i] - cached_got[i]) < 3e-4,
                        "coefficient reference evaluation differs from the NTT batch");
        }
    }
    key->cache_mode = CKKSNarrowCacheMode::Prepared;
    CKKSNarrowPrepareKey<Schedule>(*key);
    for (std::size_t column = 0; column < requests.size(); column++) {
        const auto batch_phase = trlwePhase<P>(requests[column].input->ct, source);
        CKKSSlotVector<P> got;
        ckksSlotDecrypt<P>(got, batched[column], evaluation);
        for (std::size_t i = 0; i < Schedule::slots; i++) {
            const double cycles = static_cast<double>(std::ldexp(
                ckks_detail::levelToLongDouble<P, Schedule::input_log_q>(batch_phase[i]),
                -static_cast<int>(Schedule::input_log_q)));
            const auto expected = requests[column].weight * std::polar(1.0,
                2 * std::acos(-1.0) * requests[column].frequency * cycles);
            require(std::abs(got[i] - expected) < 3e-4,
                    "shared NTT input/frequency batch has excessive error");
        }
    }
    const std::vector<std::complex<double>> series = {
        {0.2, 0}, {0.3, -0.1}, {}, {-0.05, 0.04}};
    typename Schedule::Output series_result;
    CKKSNarrowEvalSeries<Schedule>(series_result, input, *key, series);
    CKKSSlotVector<P> series_got;
    ckksSlotDecrypt<P>(series_got, series_result, evaluation);
    for (std::size_t i = 0; i < Schedule::slots; i++) {
        const double cycles = static_cast<double>(std::ldexp(
            ckks_detail::levelToLongDouble<P, Schedule::input_log_q>(phase[i]),
            -static_cast<int>(Schedule::input_log_q)));
        double expected = series[0].real();
        for (std::size_t k = 1; k < series.size(); k++)
            expected += 2 * (series[k] * std::polar(1.0, 2 * std::acos(-1.0) * k * cycles)).real();
        require(std::abs(series_got[i] - expected) < 3e-4,
                "direct-frequency functional series has excessive error");
    }
    // Exercise accumulation across tiles and independent batch owners, plus
    // a constant-only series that must not need an output conjugation key.
    const std::array<Input, 2> series_inputs{input, different_input};
    std::vector<typename Schedule::Output> series_batch;
    CKKSNarrowEvalSeriesBatch<Schedule, 12>(series_batch, series_inputs, *key, series, 1);
    for (std::size_t column = 0; column < series_inputs.size(); column++) {
        const auto batch_phase = trlwePhase<P>(series_inputs[column].ct, source);
        ckksSlotDecrypt<P>(series_got, series_batch[column], evaluation);
        for (std::size_t i = 0; i < Schedule::slots; i++) {
            const double cycles = static_cast<double>(std::ldexp(
                ckks_detail::levelToLongDouble<P, Schedule::input_log_q>(batch_phase[i]),
                -static_cast<int>(Schedule::input_log_q)));
            double expected = series[0].real();
            for (std::size_t k = 1; k < series.size(); k++)
                expected += 2 * (series[k] * std::polar(1.0, 2 * std::acos(-1.0) * k * cycles)).real();
            require(std::abs(series_got[i] - expected) < 3e-4,
                    "tiled functional series mixed owners or doubled its constant");
        }
    }
    CKKSNarrowEvalSeriesBatch<Schedule, 12>(series_batch, series_inputs, *key, {{-0.35, 0}}, 1);
    for (const auto &constant : series_batch) {
        ckksSlotDecrypt<P>(series_got, constant, evaluation);
        for (const auto value : series_got)
            require(std::abs(value - std::complex<double>{-0.35, 0}) < 1e-5,
                    "constant-only functional series has excessive error");
    }
    auto packed_key = std::make_unique<CKKSNarrowPackedKey<Schedule>>();
    CKKSNarrowPackedKeyGen<Schedule>(*packed_key, source, evaluation, {0, 3});
    if constexpr (P::n == 16) {
        std::stringstream payload;
        { cereal::PortableBinaryOutputArchive archive(payload); archive(*packed_key); }
        auto loaded = std::make_unique<CKKSNarrowPackedKey<Schedule>>();
        loaded->substitutions_fft.back() =
            std::make_unique<CKKSPreparedAutoKeyFFT<P, CKKSNarrowPackedKey<Schedule>::log_q>>();
        payload.seekg(0);
        { cereal::PortableBinaryInputArchive archive(payload); archive(*loaded); }
        require(!loaded->substitutions_fft.back(), "loaded packed key retained a stale substitution cache");
        packed_key.swap(loaded);
    }
    // Reconstruction needs only the four ciphertexts and public setup keys.
    // Verify every reconstructed selector against independently encrypted masks.
    auto head_before = key->head;
    auto tail_before = key->tail;
    CKKSNarrowReconstructSelectors<Schedule>(*key, *packed_key);
    for (std::size_t group = 0; group < 2; group++) {
        const auto &before = group == 0 ? head_before : tail_before;
        const auto &after = group == 0 ? key->head : key->tail;
        for (std::size_t j = 0; j < before.size(); j++) {
            CKKSSlotVector<P> original, reconstructed;
            ckksSlotDecrypt<P>(original, before[j], evaluation);
            ckksSlotDecrypt<P>(reconstructed, after[j], evaluation);
            for (std::size_t i = 0; i < Schedule::slots; i++)
                require(std::abs(original[i] - reconstructed[i]) < 2e-4,
                        "coefficient-packed key reconstructed an incorrect selector");
        }
    }
    typename Schedule::Output reconstructed_result;
    CKKSNarrowBlindRotate<Schedule>(reconstructed_result, input, *key);
    CKKSSlotVector<P> reconstructed_got;
    ckksSlotDecrypt<P>(reconstructed_got, reconstructed_result, evaluation);
    for (std::size_t i = 0; i < Schedule::slots; i++) {
        const double cycles = static_cast<double>(std::ldexp(
            ckks_detail::levelToLongDouble<P, Schedule::input_log_q>(phase[i]),
            -static_cast<int>(Schedule::input_log_q)));
        require(std::abs(reconstructed_got[i] - std::polar(1.0, 2 * std::acos(-1.0) * cycles)) < 3e-4,
                "blind rotation with reconstructed key has excessive error");
    }
    auto invalid = source;
    invalid[0] = invalid[1] = typename P::T{0};
    bool rejected = false;
    try { narrow_detail::structuredPositions<Schedule>(invalid); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "invalid structured key accepted");
}

struct RefreshTestParam : TestParam<4> {
    static constexpr std::uint32_t B̅gbit = 8, B̅gₐbit = 8;
    static constexpr std::uint32_t l̅ = std::numeric_limits<T>::digits / B̅gbit;
    static constexpr std::uint32_t l̅ₐ = l̅;
};

template <bool Packed = false, bool LateRescale = false, std::uint32_t AuxiliaryScale = 0,
          bool CompactDD = false, bool CompactS2C = false>
void testFullRefresh()
{
    using namespace TFHEpp;
    using P = RefreshTestParam;
    using Narrow = CKKSNarrowBootstrapSchedule<P, 36, 400, 40, 40, 4, LateRescale, AuxiliaryScale, CompactDD>;
    using Refresh = CKKSNarrowRefreshSchedule<Narrow, 32, 40, 0, CompactS2C>;
    static_assert(Refresh::output_log_q > Refresh::input_log_q + Refresh::output_log_delta);
    Key<P> sparse, dense;
    CKKSNarrowStructuredKeyGen<Narrow>(sparse);
    for (std::size_t i = 0; i < P::n; i++)
        dense[i] = i % 3 == 0 ? -P::T{1} : P::T{i % 3 == 1 ? 1U : 0U};
    auto key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    if constexpr (Packed) {
        key->narrow.head.resize(1);
        key->narrow.tail.resize(1);
        key->narrow.head_ntt.resize(1);
        key->narrow.tail_ntt.resize(1);
        key->narrow.factor_conjugation_fft =
            std::make_unique<CKKSPreparedAutoKeyFFT<P, Narrow::factor_key_log_q>>();
        CKKSNarrowRefreshSetupKeyGen<Refresh>(*key, dense, sparse, dense, {0, 3});
        require(key->narrow.head.empty() && key->narrow.tail.empty() &&
                key->narrow.head_ntt.empty() && key->narrow.tail_ntt.empty() &&
                !key->narrow.factor_conjugation_fft,
                "compressed setup retained old selectors or derived caches");
        auto packed = std::make_unique<CKKSNarrowPackedKey<Narrow>>();
        CKKSNarrowPackedKeyGen<Narrow>(*packed, sparse, dense, {0, 3});
        CKKSNarrowReconstructSelectors<Narrow>(key->narrow, *packed);
    }
    else CKKSNarrowRefreshKeyGen<Refresh>(*key, dense, sparse, dense, {0, 3});
    CKKSSlotVector<P> slots;
    for (std::size_t i = 0; i < P::n / 2; i++)
        slots[i] = {static_cast<double>(static_cast<int>(i) - 4) / 32,
                    static_cast<double>(3 - static_cast<int>(i)) / 40};
    typename Refresh::Input input;
    ckksSlotEncrypt<P>(input, slots, dense, {0, 3});
    typename Refresh::Coefficients prepared;
    CKKSNarrowPrepareCoefficients<Refresh>(prepared, input, *key);
    std::array<double, P::n> decoded_coefficients;
    ckksDecrypt<P>(decoded_coefficients, prepared, sparse);
    for (std::size_t i = 0; i < P::n / 2; i++) {
        require(std::abs(decoded_coefficients[i] - slots[i].real()) < 2e-5,
                "low-modulus S2C failed to encode real coefficients in natural order");
        require(std::abs(decoded_coefficients[i + P::n / 2] - slots[i].imag()) < 2e-5,
                "low-modulus S2C failed to encode imaginary coefficients in natural order");
    }
    typename Refresh::Output result;
    CKKSNarrowRefresh<Refresh>(result, input, *key);
    typename Refresh::Output prepared_result;
    const auto original_prepared = prepared.ct;
    CKKSNarrowRefreshPrepared<Refresh>(prepared_result, prepared, *key);
    require(prepared_result.ct == result.ct, "prepared refresh changed its ciphertext");
    require(prepared.ct == original_prepared, "prepared refresh mutated its input");
    CKKSSlotVector<P> got;
    ckksSlotDecrypt<P>(got, result, dense);
    double error = 0;
    for (std::size_t i = 0; i < P::n / 2; i++) error = std::max(error, std::abs(got[i] - slots[i]));
    std::cout << "full complex refresh error=" << error << " sine bound="
              << Refresh::sine_error_bound(0.125) << '\n';
    require(error < 3e-4, "full narrow complex refresh has excessive error");

    std::array<typename Refresh::Input, 2> refresh_inputs{input, input};
    auto different_slots = slots;
    for (auto &slot : different_slots) slot = {0.8 * slot.imag(), -0.6 * slot.real()};
    ckksSlotEncrypt<P>(refresh_inputs[1], different_slots, dense, {0, 3});
    if constexpr (Packed && LateRescale && AuxiliaryScale == 16 && CompactDD && CompactS2C) {
        using Basis = CKKSNarrowBasisRefreshSchedule<Refresh>;
        auto basis_key = std::make_unique<CKKSNarrowFrequencyBasisKey<Basis>>();
        CKKSNarrowFrequencyBasisKeyGen<Basis>(*basis_key, dense, {0, 3});
        const CKKSNarrowFrequencyBasis plan(12, {1, 3, 5, 6});
        std::vector<std::complex<double>> series(13);
        series[0] = 0.125;
        for (std::size_t k = 1; k < series.size(); k++)
            series[k] = {0.03 / (k * k), -0.02 / (k * (k + 1))};
        auto expected = [&](double message) {
            double value = series[0].real();
            for (std::size_t k = 1; k < series.size(); k++)
                value += 2 * (series[k] * std::polar(1.0,
                    2 * std::acos(-1.0) * k * message / Refresh::gap())).real();
            return value;
        };
        for (const bool complex_input : {false, true}) {
            std::vector<typename Basis::Output> basis_batch;
            CKKSNarrowFunctionalRefreshWithBasisBatch<Refresh>(basis_batch, refresh_inputs,
                *key, *basis_key, plan, series, complex_input, 3, false);
            for (std::size_t j = 0; j < basis_batch.size(); j++) {
                typename Basis::Output unary;
                CKKSNarrowFunctionalRefreshWithBasis<Refresh>(unary, refresh_inputs[j],
                    *key, *basis_key, plan, series, complex_input, 1, false);
                require(unary.ct == basis_batch[j].ct, "basis FBS batch mixed input owners or halves");
                CKKSSlotVector<P> decoded;
                ckksSlotDecrypt<P>(decoded, unary, dense);
                const auto &original = j == 0 ? slots : different_slots;
                for (std::size_t i = 0; i < decoded.size(); i++) {
                    const std::complex<double> reference{expected(original[i].real()),
                        complex_input ? expected(original[i].imag()) : 0.0};
                    require(std::abs(decoded[i] - reference) < std::ldexp(1.0, -23),
                            "basis full FBS failed unchanged 23-bit gate");
                }
            }
            CKKSNarrowFunctionalRefreshWithBasisBatch<Refresh>(basis_batch, {},
                *key, *basis_key, plan, series, complex_input);
            require(basis_batch.empty(), "empty basis FBS batch retained outputs");
        }
        std::cout << "basis full real/complex FBS and batch passed" << std::endl;
    }
    std::vector<typename Refresh::Output> refresh_batch;
    CKKSNarrowRefreshBatch<Refresh>(refresh_batch, refresh_inputs, *key, 3);
    require(refresh_batch[0].ct == result.ct, "slot refresh batch changed unary evaluation");
    typename Refresh::Output independent_refresh;
    CKKSNarrowRefresh<Refresh>(independent_refresh, refresh_inputs[1], *key);
    require(refresh_batch[1].ct == independent_refresh.ct, "slot refresh batch mixed independent inputs");
    CKKSNarrowRefreshBatch<Refresh>(refresh_batch, {}, *key);
    require(refresh_batch.empty(), "empty slot refresh batch retained previous outputs");

    // A saved key contains no derived NTT cache. Exercise both evaluation
    // immediately after loading and exact reconstruction of that cache.
    const auto original_ciphertext = result.ct;
    std::stringstream key_payload;
    { cereal::PortableBinaryOutputArchive archive(key_payload); archive(*key); }
    auto loaded_key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    loaded_key->narrow.head_ntt = key->narrow.head_ntt;
    loaded_key->narrow.tail_ntt = key->narrow.tail_ntt;
    loaded_key->narrow.factor_conjugation_fft =
        std::make_unique<CKKSPreparedAutoKeyFFT<P, Narrow::factor_key_log_q>>();
    loaded_key->narrow.output_conjugation_fft =
        std::make_unique<CKKSPreparedAutoKeyFFT<P, Narrow::output_key_log_q>>();
    loaded_key->narrow.relin_fft = std::make_unique<typename Narrow::PreparedRelinKeys>();
    loaded_key->narrow.tail_rotations_fft =
        std::make_unique<CKKSPreparedRotationKeyFFT<P, Narrow::factor_key_log_q>>();
    key_payload.seekg(0);
    { cereal::PortableBinaryInputArchive archive(key_payload); archive(*loaded_key); }
    require(!loaded_key->narrow.factor_conjugation_fft && !loaded_key->narrow.output_conjugation_fft,
            "loaded key retained a stale prepared automorphism cache");
    require(!loaded_key->narrow.relin_fft, "loaded key retained a stale prepared DD cache");
    require(!loaded_key->narrow.tail_rotations_fft, "loaded key retained a stale rotation cache");
    require(loaded_key->narrow.head_ntt.empty() && loaded_key->narrow.tail_ntt.empty(),
            "deserialized narrow key retained a stale derived NTT cache");
    require(loaded_key->narrow.cache_mode == CKKSNarrowCacheMode::Prepared,
            "deserialized narrow key lost its cache mode");
    CKKSNarrowRefresh<Refresh>(result, input, *loaded_key);
    ckksSlotDecrypt<P>(got, result, dense);
    for (std::size_t i = 0; i < P::n / 2; i++)
        require(std::abs(got[i] - slots[i]) < 3e-4,
                "deserialized refresh key failed coefficient reference evaluation");
    CKKSNarrowPrepareKey<Narrow>(loaded_key->narrow);
    CKKSNarrowRefresh<Refresh>(result, input, *loaded_key);
    require(result.ct == original_ciphertext,
            "deserialized refresh key changed evaluation after rebuilding its NTT cache");
    key.swap(loaded_key);

    const auto relin = makeCKKSDDRelinKey<P,
        Refresh::output_log_q - Refresh::output_log_delta>(dense, {0, 3});
    CKKSMultResult<P, Refresh::output_log_q, Refresh::output_log_delta,
        Refresh::output_log_q, Refresh::output_log_delta> product;
    CKKSMult<P>(product, result, result, *relin);
    ckksSlotDecrypt<P>(got, product, dense);
    for (std::size_t i = 0; i < P::n / 2; i++)
        require(std::abs(got[i] - slots[i] * slots[i]) < 3e-4,
                "refreshed ciphertext failed a subsequent DD multiplication");

    CKKSScaleDownResult<P, Refresh::output_log_q, Refresh::output_log_delta,
        Refresh::output_log_delta - Refresh::input_log_delta> lowered;
    CKKSScaleDown<P, Refresh::output_log_q, Refresh::output_log_delta,
        Refresh::output_log_delta - Refresh::input_log_delta>(lowered, result);
    typename Refresh::Input repeated_input;
    CKKSLevelReduce<P>(repeated_input, lowered);
    CKKSNarrowRefresh<Refresh>(result, repeated_input, *key);
    ckksSlotDecrypt<P>(got, result, dense);
    for (std::size_t i = 0; i < P::n / 2; i++)
        require(std::abs(got[i] - slots[i]) < 5e-4,
                "narrow refresh did not preserve a repeated bootstrap cycle");

    const std::vector<std::complex<double>> coefficients = {{0.2, 0}, {0.15, 0}, {}, {0, -0.1}};
    CKKSNarrowFunctionalRefresh<Refresh>(result, input, *key, coefficients);
    ckksSlotDecrypt<P>(got, result, dense);
    auto target = [](double x) {
        const double angle = 2 * std::acos(-1.0) * x / Refresh::gap();
        return 0.2 + 0.3 * std::cos(angle) + 0.2 * std::sin(3 * angle);
    };
    for (std::size_t i = 0; i < P::n / 2; i++)
        require(std::abs(got[i] - std::complex<double>{target(slots[i].real()), target(slots[i].imag())}) < 3e-4,
                "full complex functional refresh has excessive error");

    const auto exact = CKKSBuildNarrowTrigonometricApproximation(target, 4, Refresh::gap(), 3);
    require(exact.max_sampled_error < 1e-10, "Remez failed to reproduce a trigonometric polynomial");
    for (std::size_t k = 0; k < coefficients.size(); k++)
        require(std::abs(exact.coefficients[k] - coefficients[k]) < 1e-9,
                "Remez used the wrong positive-frequency coefficient convention");
    auto tanh_approximation = CKKSBuildNarrowTrigonometricApproximation(
        [](double x) { return std::tanh(x); }, 4, Refresh::gap(), 6, 1.0);
    const auto curvature_approximation = CKKSBuildNarrowTrigonometricApproximation(
        [](double x) { return std::tanh(x); }, 4, Refresh::gap(), 6, 1.0, 8192, 40, 2.0);
    require(curvature_approximation.coefficients == tanh_approximation.coefficients &&
                curvature_approximation.uniform_error_bound &&
                *curvature_approximation.uniform_error_bound < *tanh_approximation.uniform_error_bound &&
                *curvature_approximation.uniform_error_bound >= curvature_approximation.max_sampled_error,
            "curvature validation changed coefficients or failed to tighten the bound");
    for (std::size_t i = 0; i <= 32768; i++) {
        const double x = -4 + 8.0 * i / 32768;
        require(std::abs(CKKSEvaluateNarrowApproximation(curvature_approximation, x) - std::tanh(x)) <=
                    *curvature_approximation.uniform_error_bound,
                "curvature uniform bound failed an independent grid");
    }
    const auto tanh_curvature = [](double x) {
        const double y = std::tanh(x);
        return -2 * y * (1 - y * y);
    };
    const auto refined_tanh = CKKSBuildNarrowTrigonometricApproximation(
        [](double x) { return std::tanh(x); }, 4, 16, 20, 1.0, 65536, 40, 2.0,
        tanh_curvature, 40.0);
    require(refined_tanh.uniform_error_bound &&
                std::sqrt(2.0) * *refined_tanh.uniform_error_bound < std::ldexp(1.0, -23),
            "residual curvature failed to establish the degree-20 offline tanh budget");
    for (std::size_t i = 0; i <= 32768; i++) {
        const double x = -4 + 8.0 * i / 32768;
        require(std::abs(CKKSEvaluateNarrowApproximation(refined_tanh, x) - std::tanh(x)) <=
                    *refined_tanh.uniform_error_bound,
                "refined residual bound failed an independent grid");
    }
    const auto windowed_tanh = CKKSBuildNarrowWindowedTrigonometricApproximation(
        [](double x) { return std::tanh(x); },
        [](double x) { const double y = std::tanh(x); return 1 - y * y; },
        tanh_curvature, {1, 1, 2, 4, 40}, 4, 16, 40);
    double windowed_norm = std::abs(windowed_tanh.coefficients[0]);
    for (std::size_t k = 1; k < windowed_tanh.coefficients.size(); k++)
        windowed_norm += 2 * std::abs(windowed_tanh.coefficients[k]);
    require(windowed_norm < 2 && windowed_tanh.message_bound == 4 &&
                std::sqrt(2.0) * *windowed_tanh.uniform_error_bound < std::ldexp(1.0, -23),
            "windowed Remez failed the protected domain or coefficient budget");
    for (std::size_t i = 0; i <= 32768; i++) {
        const double x = -4 + 8.0 * i / 32768;
        require(std::abs(CKKSEvaluateNarrowApproximation(windowed_tanh, x) - std::tanh(x)) <=
                    *windowed_tanh.uniform_error_bound,
                "windowed Remez changed tanh within the protected domain");
    }
    bool rejected_incomplete_curvature = false;
    try {
        CKKSBuildNarrowTrigonometricApproximation(
            [](double x) { return x; }, 1, 4, 2, {}, 8192, 40, {}, tanh_curvature);
    }
    catch (const std::invalid_argument &) { rejected_incomplete_curvature = true; }
    require(rejected_incomplete_curvature, "residual curvature accepted an incomplete derivative contract");
    std::stringstream approximation_payload;
    { cereal::PortableBinaryOutputArchive archive(approximation_payload); archive(tanh_approximation); }
    CKKSNarrowTrigonometricApproximation loaded_approximation;
    approximation_payload.seekg(0);
    { cereal::PortableBinaryInputArchive archive(approximation_payload); archive(loaded_approximation); }
    require(loaded_approximation.coefficients == tanh_approximation.coefficients &&
            loaded_approximation.message_bound == tanh_approximation.message_bound &&
            loaded_approximation.gap == tanh_approximation.gap &&
            loaded_approximation.max_sampled_error == tanh_approximation.max_sampled_error &&
            loaded_approximation.uniform_error_bound == tanh_approximation.uniform_error_bound &&
            loaded_approximation.iterations == tanh_approximation.iterations,
            "serialized approximation lost coefficients or validation metadata");
    tanh_approximation = std::move(loaded_approximation);
    require(tanh_approximation.uniform_error_bound.has_value(), "Remez omitted Lipschitz validation");
    require(tanh_approximation.max_sampled_error < 0.005, "tanh Remez accuracy is insufficient");
    for (std::size_t i = 0; i < 10001; i++) {
        const double x = -4 + 8.0 * i / 10000;
        require(std::abs(CKKSEvaluateNarrowApproximation(tanh_approximation, x) - std::tanh(x)) <=
                    *tanh_approximation.uniform_error_bound,
                "Remez validation bound failed on an independent grid");
    }
    // Exercise the whole declared approximation domain in ciphertexts,
    // including both endpoints, rather than only the refresh's small inputs.
    for (std::size_t i = 0; i < P::n / 2; i++) {
        const double x = -4 + 8.0 * i / (P::n / 2 - 1);
        slots[i] = {x, -x};
    }
    ckksSlotEncrypt<P>(input, slots, dense, {0, 3});
    CKKSNarrowFunctionalRefresh<Refresh>(result, input, *key, tanh_approximation);
    CKKSNarrowPrepareCoefficients<Refresh>(prepared, input, *key);
    const auto functional_prepared = prepared.ct;
    CKKSNarrowFunctionalRefreshPrepared<Refresh>(prepared_result, prepared, *key, tanh_approximation);
    require(prepared_result.ct == result.ct, "prepared complex FBS changed its ciphertext");
    require(prepared.ct == functional_prepared, "prepared complex FBS mutated its input");
    typename Refresh::Output wide_result;
    CKKSNarrowFunctionalRefresh<Refresh>(wide_result, input, *key, tanh_approximation, true, 12);
    require(wide_result.ct == result.ct, "wide functional refresh changed its ciphertext");
    key->narrow.cache_mode = CKKSNarrowCacheMode::Streamed;
    CKKSNarrowPrepareKey<Narrow>(key->narrow);
    CKKSNarrowFunctionalRefresh<Refresh>(wide_result, input, *key, tanh_approximation, true, 12);
    require(wide_result.ct == result.ct, "streamed wide functional refresh changed its ciphertext");
    key->narrow.cache_mode = CKKSNarrowCacheMode::Prepared;
    CKKSNarrowPrepareKey<Narrow>(key->narrow);
    typename Refresh::Output real_only, real_wide;
    CKKSNarrowFunctionalRefresh<Refresh>(real_only, input, *key, tanh_approximation, false, 1);
    CKKSNarrowFunctionalRefresh<Refresh>(real_wide, input, *key, tanh_approximation, false, 6);
    require(real_only.ct == real_wide.ct, "real-only functional refresh ignored its tiling convention");
    CKKSNarrowFunctionalRefreshPrepared<Refresh>(prepared_result, prepared, *key, tanh_approximation, false, 6);
    require(prepared_result.ct == real_only.ct, "prepared real FBS changed its ciphertext");
    require(prepared.ct == functional_prepared, "prepared real FBS mutated its input");
    ckksSlotDecrypt<P>(got, real_only, dense);
    for (std::size_t i = 0; i < P::n / 2; i++)
        require(std::abs(got[i] - CKKSEvaluateNarrowApproximation(tanh_approximation, slots[i].real())) < 6e-4,
                "real-only functional refresh has excessive error");
    const std::array<typename Refresh::Input, 2> functional_inputs{input, refresh_inputs[1]};
    std::vector<typename Refresh::Output> functional_batch;
    CKKSNarrowFunctionalRefreshBatch<Refresh>(functional_batch, functional_inputs, *key,
        tanh_approximation, true, 5);
    require(functional_batch[0].ct == result.ct, "functional batch changed unary evaluation");
    CKKSNarrowFunctionalRefresh<Refresh>(wide_result, functional_inputs[1], *key, tanh_approximation);
    require(functional_batch[1].ct == wide_result.ct, "functional batch mixed independent inputs");
    CKKSNarrowFunctionalRefreshBatch<Refresh>(functional_batch, functional_inputs, *key,
        tanh_approximation, false, 5);
    require(functional_batch[0].ct == real_only.ct, "real-only functional batch changed unary evaluation");
    CKKSNarrowFunctionalRefresh<Refresh>(wide_result, functional_inputs[1], *key, tanh_approximation, false);
    require(functional_batch[1].ct == wide_result.ct, "real-only functional batch mixed inputs");
    key->narrow.cache_mode = CKKSNarrowCacheMode::Streamed;
    CKKSNarrowPrepareKey<Narrow>(key->narrow);
    CKKSNarrowFunctionalRefreshBatch<Refresh>(functional_batch, functional_inputs, *key,
        tanh_approximation, true, 5);
    require(functional_batch[0].ct == result.ct, "streamed functional batch changed unary evaluation");
    key->narrow.cache_mode = CKKSNarrowCacheMode::Prepared;
    CKKSNarrowPrepareKey<Narrow>(key->narrow);
    bool rejected_zero_tile = false;
    try { CKKSNarrowFunctionalRefresh<Refresh>(wide_result, input, *key, tanh_approximation, true, 0); }
    catch (const std::invalid_argument &) { rejected_zero_tile = true; }
    require(rejected_zero_tile, "functional refresh accepted a zero tile size");
    ckksSlotDecrypt<P>(got, result, dense);
    for (std::size_t i = 0; i < P::n / 2; i++) {
        const std::complex<double> approximate = {
            CKKSEvaluateNarrowApproximation(tanh_approximation, slots[i].real()),
            CKKSEvaluateNarrowApproximation(tanh_approximation, slots[i].imag())};
        require(std::abs(got[i] - approximate) < 6e-4,
                "full Remez functional refresh has excessive encrypted evaluation error");
        require(std::abs(got[i] - std::complex<double>{std::tanh(slots[i].real()), std::tanh(slots[i].imag())}) <
                    2 * *tanh_approximation.uniform_error_bound + 6e-4,
                "full Remez functional refresh failed its declared approximation domain");
    }
    std::cout << "tanh Remez sampled error=" << tanh_approximation.max_sampled_error
              << " uniform bound=" << *tanh_approximation.uniform_error_bound << '\n';

    const std::array<std::function<double(double)>, 2> targets = {
        [](double x) { return 1 / (1 + std::exp(-x)); },
        [](double x) { return x * (1 + std::erf(x / std::sqrt(2.0))) / 2; }};
    const std::array<double, 2> lipschitz_bounds = {0.25, 1.13};
    const std::array<const char *, 2> names = {"sigmoid", "GeLU"};
    for (std::size_t target_index = 0; target_index < targets.size(); target_index++) {
        const auto approximation = CKKSBuildNarrowTrigonometricApproximation(
            targets[target_index], 4, Refresh::gap(), 6, lipschitz_bounds[target_index]);
        require(approximation.max_sampled_error < std::ldexp(1.0, -7),
                "nonlinear FBS approximation missed the seven-bit accuracy target");
        CKKSNarrowFunctionalRefresh<Refresh>(result, input, *key, approximation);
        ckksSlotDecrypt<P>(got, result, dense);
        for (std::size_t i = 0; i < P::n / 2; i++) {
            const std::complex<double> approximate{
                CKKSEvaluateNarrowApproximation(approximation, slots[i].real()),
                CKKSEvaluateNarrowApproximation(approximation, slots[i].imag())};
            const std::complex<double> target{
                targets[target_index](slots[i].real()), targets[target_index](slots[i].imag())};
            require(std::abs(got[i] - approximate) < 6e-4,
                    "nonlinear FBS encrypted evaluation missed its error budget");
            require(std::abs(got[i] - target) <
                        2 * *approximation.uniform_error_bound + 6e-4,
                    "nonlinear FBS failed the combined approximation/evaluation bound");
        }
        std::cout << names[target_index] << " Remez sampled error="
                  << approximation.max_sampled_error << " uniform bound="
                  << *approximation.uniform_error_bound << '\n';
    }
}

}  // namespace

int main(int argc, char **argv)
{
    if (argc > 1) {
        if (argc == 2 && std::string(argv[1]) == "--prepared-refresh") {
            testFullRefresh<>();
            testFullRefresh<true, true, 16, true, true>();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--selector-digits") {
            testSelectorDigitExtraction(); return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--prepared-auto-parallel") {
            testPreparedAuto<85, 13>();
            testPreparedAuto<379, 13>();
            return 0;
        }
        if (argc == 3 && std::string(argv[1]) == "--packed-workspace-fixture-write") {
            packedWorkspaceFixture(argv[2], true); return 0;
        }
        if (argc == 3 && std::string(argv[1]) == "--packed-workspace-fixture-check") {
            packedWorkspaceFixture(argv[2], false); return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--packed-selectors") {
            testPreparedPackedSubstitutions();
            return 0;
        }

        if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--frequency-basis-fbs") {
            testGaussianBasisFBS(argc == 3 ? argv[2] : "");
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--frequency-basis") {
            testFrequencyBasis<6, 7>();
            testFrequencyBasis<9, 31>();
            testFullRefresh<true, true, 16, true, true>();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--staged-fft") {
            testFileSelectors();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--selector-chunks") {
            testSelectorChunks<18>();
            testSelectorChunks<40>();
            testSelectorChunks<18,48>();
            testSelectorChunks<40,48>();
            testFileSelectors();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--bounded-rotations") {
            testBoundedRotations<0>();
            testBoundedRotations<28>();
            testFileSelectors();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--compact-s2c") {
            testCompactS2C();
            testFullRefresh<true, true, 16, true, true>();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--compact-dd-wide") {
            testCompactDD<97, 12>();
            testCompactDD<97, 24>();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--compact-dd") {
            testCompactDD<17>();
            testCompactDD<85>();
            testCompactDD<379>();
            testFullRefresh<true, true, 16, true>();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--file-selectors") {
            testFileSelectors();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--auxiliary-auto") {
            testAuxiliaryScaleAuto();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--relin-first") {
            testRelinThenRescale();
            return 0;
        }
        return 2;
    }
    using namespace TFHEpp;
    testFrequencyBasis<6, 7>();
    testFrequencyBasis<9, 31>();
    testGaussianBasisFBS();
    testDigitFFT<9>();
    testDigitFFT<13>();
    testPreparedAuto<17>();
    testPreparedAuto<85>();
    testPreparedAuto<379>();
    testPreparedAuto<448>();
    testPreparedAuto<85, 14>();
    testPreparedDDTree<>();
    testPreparedDDTree<true>();
    testPreparedTailRotations();
    testSelectorChunks<18>();
    testSelectorChunks<40>();
    testSelectorChunks<18,48>();
    testSelectorChunks<40,48>();
    testBoundedRotations<0>();
    testBoundedRotations<28>();
    testPreparedPackedSubstitutions();
    testRelinThenRescale();
    testAuxiliaryScaleAuto();
    testCompactDD<17>();
    testCompactDD<85>();
    testCompactDD<379>();
    testCompactDD<97, 12>();
    testCompactDD<97, 24>();
    testFileSelectors();
    testCompactS2C();
    testActiveDigitProduct<17>();
    testActiveDigitProduct<85>();
    testActiveDigitProduct<379>();
    testActiveDigitProduct<448>();
    using S1 = CKKSNarrowBootstrapSchedule<TestParam<4>, 20, 400, 40, 40, 4>;
    using S2 = CKKSNarrowBootstrapSchedule<TestParam<5>, 20, 400, 40, 40, 8>;
    using S3 = CKKSNarrowBootstrapSchedule<TestParam<6>, 20, 400, 40, 40, 16>;
    using S4 = CKKSNarrowBootstrapSchedule<TestParam<6>, 20, 400, 40, 40, 7>;
    using SinglePrime = CKKSNarrowBootstrapSchedule<TestParam<4>, 20, 400, 40, 30, 4>;
    static_assert(CKKSNarrowNTTSelector<SinglePrime>::prime_count == 1);
    static_assert(CKKSNarrowNTTSelector<S1>::prime_count == 2);
    using PaperWeight = CKKSNarrowBootstrapSchedule<TestParam<9>, 20, 400, 40, 40, 31>;
    testTables<S1>();
    testTables<S2>();
    testTables<S3>();  // B=5 and h=16: includes a partial final tail pack.
    testTables<S4>();  // Non-dividing h with a shorter final block.
    testTables<PaperWeight>();  // The paper's actual h=31, including odd tree topology.
    testEncrypted<S1>();
    testEncrypted<S2>();
    testEncrypted<S3>();
    testEncrypted<S4>();
    using LateOdd = CKKSNarrowBootstrapSchedule<TestParam<6>, 20, 400, 40, 40, 7, true>;
    static_assert(std::remove_reference_t<decltype(std::declval<typename LateOdd::RelinKeys>().template get<0>())>::log_q == LateOdd::factor_log_q);
    testEncrypted<LateOdd>();
    testEncrypted<SinglePrime>();
    testFullRefresh<>();
    testFullRefresh<true>();
    testFullRefresh<true, true>();
    testFullRefresh<true, true, 16>();
    testFullRefresh<true, true, 16, true>();
    std::cout << "CKKS narrow DD blind rotation passed\n";
}
