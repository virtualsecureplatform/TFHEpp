#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>

namespace {
template <std::uint32_t Nbit>
struct ProbeParam {
    static constexpr std::int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = Nbit, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t l = 1, l_a = 1, lₐ = 1;
    static constexpr std::uint32_t Bgbit = 4, Bg_abit = 4, Bgₐbit = 4;
    using T = TFHEpp::MultiLimbUInt<7>;
    static constexpr T Bg = T{1} << Bgbit, Bg_a = Bg, Bgₐ = Bg;
    static constexpr TFHEpp::ErrorDistribution errordist = TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::ldexp(1.0, -448);
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << 253;
    static constexpr std::uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 112, l̅ₐ = l̅, B̅gbit = 4, B̅gₐbit = 4;
};
using Clock = std::chrono::steady_clock;

double maximum(double current, double value)
{
    if (!std::isfinite(value)) throw std::runtime_error("nonfinite error measurement");
    return std::max(current, value);
}

template <std::uint32_t Nbit>
int run(const std::filesystem::path &root, bool inject_failure)
{
    using namespace TFHEpp;
    using P = ProbeParam<Nbit>;
    using Narrow = CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, true, 8, true>;
    using Refresh = CKKSNarrowRefreshSchedule<Narrow, 44, 37, 0, true>;
    static_assert(Refresh::input_log_q == 85 && Refresh::output_log_q == 108);
    static_assert(Refresh::evaluation_input_log_q == 93);
    if (!std::filesystem::create_directory(root)) throw std::runtime_error("probe directory must be new");
    std::filesystem::permissions(root, std::filesystem::perms::owner_all,
        std::filesystem::perm_options::replace);
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    std::mt19937_64 random(20262119);
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : *dense) value = P::T::from_signed_i64(ternary(random));
    CKKSNarrowStructuredKeyGen<Narrow>(*sparse);
    const CKKSNoise noise{P::α, 0};
    auto key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    std::cout << "generating N" << P::n << " auxiliary8 probe keys" << std::endl;
    CKKSNarrowRefreshSetupKeyGen<Refresh>(*key, *dense, *sparse, *dense,
        noise, CKKSNarrowCacheMode::Streamed);
    CKKSNarrowEncryptSelectorsToFile<Narrow>(key->narrow, *sparse, *dense, root / "selectors", noise);
    key->narrow.selector_ntt_chunk_size = std::min<std::size_t>(40, Narrow::block_width);
    key->narrow.bounded_tail_rotation_fft = true;
    key->narrow.staged_fft_caches = true;
    auto square_key = makeCKKSCompactDDRelinKey<P, 108>(*dense, noise);
    auto input = std::make_unique<typename Refresh::Input>();
    auto output = std::make_unique<typename Refresh::Output>();
    auto square = std::make_unique<CKKSMultResult<P, 108, 48, 108, 48>>();
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto decoded = std::make_unique<CKKSSlotVector<P>>();
    auto target_slots = std::make_unique<CKKSSlotVector<P>>();
    auto approximate_slots = std::make_unique<CKKSSlotVector<P>>();
    bool all_passed = true;
    for (const std::string function : {"sigmoid", "gelu", "tanh"}) {
        const auto target = [&](double x) {
            if (function == "sigmoid") return 1 / (1 + std::exp(-x));
            if (function == "gelu") return 0.5 * x * (1 + std::erf(x / std::sqrt(2.0)));
            return std::tanh(x);
        };
        const auto curvature = [&](double x) {
            const double y = target(x);
            if (function == "sigmoid") return y * (1-y) * (1-2*y);
            if (function == "gelu") return (2-x*x) * std::exp(-x*x/2) / std::sqrt(2*std::acos(-1.0));
            return -2*y*(1-y*y);
        };
        const auto approximation = function == "tanh" ?
            CKKSBuildNarrowWindowedTrigonometricApproximation(target,
                [](double x) { const double y = std::tanh(x); return 1-y*y; },
                curvature, {1, 1, 2, 4, 40}, 4, Refresh::gap(), 40) :
            CKKSBuildNarrowTrigonometricApproximation(target, 4, Refresh::gap(), 12,
                function == "gelu" ? 2.0 : 1.0, 65536, 40, 2.0, curvature, 160.0);
        std::uniform_real_distribution<double> message(-4, 4);
        for (auto &slot : *slots) slot = {message(random), message(random)};
        slots->front() = {-4, 4}; slots->back() = {4, -4};
        (*slots)[1] = {4, 4}; (*slots)[2] = {0, 0};
        ckksSlotEncrypt<P>(*input, *slots, *dense, noise);
        ckksSlotDecrypt<P>(*decoded, *input, *dense);
        double fresh = 0;
        for (std::size_t i = 0; i < P::n/2; ++i)
            fresh = maximum(fresh, std::abs((*decoded)[i] - (*slots)[i]));
        bool refresh_finished = false;
        bool capture_attempted = false;
        auto capture = [&] {
            if (capture_attempted) return;
            capture_attempted = true;
            const auto source = std::dynamic_pointer_cast<const CKKSNarrowSelectorFile<Narrow>>(key->narrow.selector_source);
            if (!source) throw std::runtime_error("probe selectors are not file-backed");
            source->preserve();
            std::ofstream file(root / (function + "-failure.bin"), std::ios::binary);
            if (!file) throw std::runtime_error("cannot create probe failure state");
            { cereal::PortableBinaryOutputArchive archive(file);
              const std::array<std::uint64_t, 10> metadata{
                  0x4155583850523031ULL, P::n, 85, 44, 385, 48, 37, 108, 8, 448};
              archive(metadata, function, refresh_finished, *dense, *sparse, *slots,
                      *input, *output, *key, *square_key); }
            file.close();
            if (!file) throw std::runtime_error("cannot close probe failure state");
        };
        try {
            std::cout << "evaluating " << function << std::endl;
            const auto start = Clock::now();
            CKKSNarrowFunctionalRefresh<Refresh>(*output, *input, *key, approximation, true, 12);
            refresh_finished = true;
            const double seconds = std::chrono::duration<double>(Clock::now()-start).count();
            ckksSlotDecrypt<P>(*decoded, *output, *dense);
            double error = 0, evaluation = 0, squared_error = 0, squared_evaluation = 0;
            double sum_squared = 0, amplitude = 0, approximate_amplitude = 0;
            for (std::size_t i = 0; i < P::n/2; ++i) {
                const auto x = (*slots)[i];
                (*target_slots)[i] = {target(x.real()), target(x.imag())};
                (*approximate_slots)[i] = {CKKSEvaluateNarrowApproximation(approximation, x.real()),
                    CKKSEvaluateNarrowApproximation(approximation, x.imag())};
                const double e = std::abs((*decoded)[i]-(*target_slots)[i]);
                error = maximum(error, e); sum_squared += e*e;
                evaluation = maximum(evaluation, std::abs((*decoded)[i]-(*approximate_slots)[i]));
                amplitude = maximum(amplitude, std::abs((*target_slots)[i]));
                approximate_amplitude = maximum(approximate_amplitude, std::abs((*approximate_slots)[i]));
            }
            CKKSMultRelinThenRescale<P>(*square, *output, *output, *square_key);
            ckksSlotDecrypt<P>(*decoded, *square, *dense);
            for (std::size_t i = 0; i < P::n/2; ++i) {
                squared_error = maximum(squared_error, std::abs((*decoded)[i]-(*target_slots)[i]*(*target_slots)[i]));
                squared_evaluation = maximum(squared_evaluation,
                    std::abs((*decoded)[i]-(*approximate_slots)[i]*(*approximate_slots)[i]));
            }
            const double tolerance = std::ldexp(1.0, -23);
            const double uniform = std::sqrt(2.0) * approximation.uniform_error_bound.value();
            const double rms = std::sqrt(sum_squared/(P::n/2));
            const double product_tolerance = 2*amplitude*tolerance+tolerance*tolerance+tolerance;
            const double evaluation_product_tolerance = 2*approximate_amplitude*tolerance+tolerance*tolerance+tolerance;
            const bool passed = !inject_failure && fresh < tolerance && error < tolerance && evaluation < tolerance &&
                rms < tolerance && uniform < tolerance && squared_error < product_tolerance &&
                squared_evaluation < evaluation_product_tolerance;
            all_passed &= passed;
            std::ofstream file(root / (function + ".json"));
            file << std::setprecision(17) << "{\"profile\":\"auxiliary8-"
                 << (Nbit == 9 ? "small" : "full") << "-probe\",\"n\":" << P::n << ",\"h\":31,"
                 << "\"function\":\"" << function << "\",\"auxiliary_bits\":8,\"gaussian_sigma\":3.2,"
                 << "\"source_log_q\":48,\"selector_log_q\":385,\"structured_switch_key_log_q\":56,"
                 << "\"factor_key_log_q\":356,\"s2c_input_key_log_q\":93,\"output_conjugation_key_log_q\":116,"
                 << "\"input_log_q\":85,\"input_log_delta\":44,\"output_log_q\":108,\"output_log_delta\":48,"
                 << "\"root_plain_log_delta\":37,\"s2c_plain_log_delta\":37,\"storage_bits\":448,"
                 << "\"series_tile_size\":12,\"degree\":" << (function == "tanh" ? 40 : 12)
                 << ",\"compact_dd\":true,\"compact_s2c\":true,\"relinearize_before_rescale\":true,"
                 << "\"packed_selectors\":false,\"file_backed_selectors\":true,\"input_pattern\":\"random-with-endpoints\","
                 << "\"message_bound\":4,\"target_precision_bits\":23,\"max_error\":" << error
                 << ",\"evaluation_max_error\":" << evaluation << ",\"fresh_input_max_error\":" << fresh
                 << ",\"rms_error\":" << rms << ",\"complex_uniform_bound\":" << uniform
                 << ",\"post_refresh_product_max_error\":" << squared_error
                 << ",\"post_refresh_product_evaluation_max_error\":" << squared_evaluation
                 << ",\"product_tolerance\":" << product_tolerance
                 << ",\"evaluation_product_tolerance\":" << evaluation_product_tolerance
                 << ",\"refresh_seconds\":" << seconds << ",\"accuracy_passed\":" << (passed ? "true" : "false")
                 << ",\"injected_failure\":" << (inject_failure ? "true" : "false") << "}\n";
            file.close(); if (!file) throw std::runtime_error("cannot write probe report");
            if (!passed) { capture(); return 1; }
        }
        catch (...) { capture(); throw; }
    }
    return all_passed ? 0 : 1;
}
} // namespace

int main(int argc, char **argv)
{
    try {
        if (argc == 2) return run<9>(argv[1], false);
        if (argc == 3 && std::string(argv[2]) == "--inject-failure") return run<9>(argv[1], true);
        if (argc == 3 && std::string(argv[2]) == "--full") return run<14>(argv[1], false);
        if (argc == 4 && std::string(argv[2]) == "--full" && std::string(argv[3]) == "--inject-failure")
            return run<14>(argv[1], true);
        std::cerr << "Usage: ckks_narrow_auxiliary_probe NEW_DIRECTORY [--full] [--inject-failure]\n";
        return 2;
    }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
}
