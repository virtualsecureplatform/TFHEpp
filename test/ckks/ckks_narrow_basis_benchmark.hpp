#pragma once

// Shared CPU benchmark implementation: the dedicated retained replay runner
// instantiates only the two lean-basis profiles, avoiding every other benchmark.
#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <ckks/narrow_selector_ntt_file.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>

namespace narrow_benchmark {

template <class T>
T finiteMaximum(T previous, T candidate)
{
    if (!std::isfinite(previous) || !std::isfinite(candidate))
        throw std::runtime_error("nonfinite narrow validation metric");
    return std::max(previous, candidate);
}


template <std::uint32_t Nbit, std::size_t Limbs = 4>
struct BenchmarkParam {
    static constexpr int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = Nbit, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t l = 1, l_a = 1, lₐ = 1;
    static constexpr std::uint32_t Bgbit = 4, Bg_abit = 4, Bgₐbit = 4;
    using T = TFHEpp::MultiLimbUInt<Limbs>;
    static constexpr T Bg = T{1} << Bgbit, Bg_a = T{1} << Bg_abit, Bgₐ = Bg_a;
    static constexpr TFHEpp::ErrorDistribution errordist = TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::ldexp(1.0, -static_cast<int>(64 * Limbs));
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << 253;
    static constexpr std::uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 16 * Limbs, l̅ₐ = l̅, B̅gbit = 4, B̅gₐbit = 4;
};

// Separate DD candidate: wider primary/secondary digits keep the FFT
// accumulation within its existing mantissa guard while shrinking caches.
template <std::uint32_t Nbit, std::size_t Limbs = 9>
struct LeanBasisParam : BenchmarkParam<Nbit, Limbs> {
    using T = typename BenchmarkParam<Nbit, Limbs>::T;
    static constexpr std::uint32_t Bgbit = 18;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr std::uint32_t B̅gbit = 8, B̅gₐbit = 8;
    static constexpr std::uint32_t l̅ = 8 * Limbs, l̅ₐ = l̅;
};


// Optional replay policy: the serialized coefficient keys remain unchanged.
// Only the derived product-tree FFT cache uses the tested secondary base.
template <class Base>
struct Secondary11BasisSchedule : Base {
    using P = typename Base::Param;
    static constexpr std::uint32_t prepared_relin_bbarbit = 11;
    static constexpr auto first_q = Base::relin_chain_start_log_q - Base::log_delta;
    static_assert(P::Bgbit + prepared_relin_bbarbit + P::nbit +
        std::bit_width(TFHEpp::CKKSDDRelinPrimaryRowCountForLevel<P, first_q>()) + 3 <
        std::numeric_limits<double>::digits);
    static_assert(((first_q + prepared_relin_bbarbit - 1) / prepared_relin_bbarbit) *
        prepared_relin_bbarbit <= std::numeric_limits<typename P::T>::digits);
    using PreparedRelinKeys = TFHEpp::CKKSPreparedDDRelinKeyChainFFT<P,
        Base::relin_chain_start_log_q, Base::log_delta, Base::tree_depth,
        P::Bgbit, prepared_relin_bbarbit>;
};

template <class Base>
struct Selector56BasisSchedule : Base {
    static constexpr std::size_t selector_matrix_digit_bits = 56;
    static_assert(TFHEpp::CKKSNarrowNTTSelector<Base,56>::fits(56));
};

template <std::uint32_t Nbit, bool Lean = false, bool Secondary11 = false,
          bool Selector56 = false>
int basisBenchmark(const std::string &function_option, const std::string &directory,
                   const std::string &output_path, bool storage_only = false,
                   const std::string &checkpoint_directory = "", bool replay_checkpoint = false,
                   const std::string &prepared_selector_directory = "")
try {
    using namespace TFHEpp;
    using P = std::conditional_t<Lean, LeanBasisParam<Nbit>, BenchmarkParam<Nbit, 10>>;
    constexpr auto selector_bits = Lean ? 518 : 536;
    constexpr auto root_delta = Lean ? 47 : 49;
    constexpr auto coefficient_bits = Lean ? 31 : 32;
    using OriginalNarrow = CKKSNarrowBootstrapSchedule<P, 48, selector_bits, root_delta, 39, 31, true, 44, true>;
    using DDSchedule = std::conditional_t<Secondary11, Secondary11BasisSchedule<OriginalNarrow>, OriginalNarrow>;
    using Narrow = std::conditional_t<Selector56, Selector56BasisSchedule<DDSchedule>, DDSchedule>;
    using Refresh = CKKSNarrowRefreshSchedule<Narrow, 44, 37, 0, true>;
    using Basis = CKKSNarrowBasisRefreshSchedule<Refresh, coefficient_bits>;
    using S2CTerm = typename CKKSCompactLinearTransformPlan<P, Refresh::evaluation_input_log_q,
        Refresh::evaluation_input_log_delta, Refresh::s2c_plain_log_delta>::Term;
    using Square = CKKSMultResult<P, Basis::output_log_q, Basis::log_delta,
        Basis::output_log_q, Basis::log_delta>;
    static_assert(Basis::output_log_q == (Lean ? 197 : 203) && Basis::log_delta == (Lean ? 78 : 81));
    static_assert(Square::log_q - (Basis::log_delta - Refresh::input_log_delta) == Refresh::input_log_q);
    if constexpr (Secondary11)
        if (!replay_checkpoint || storage_only)
            throw std::invalid_argument("secondary11 mode requires retained-key replay");
    if constexpr (Selector56)
        if (!replay_checkpoint || storage_only)
            throw std::invalid_argument("selector56 mode requires retained-key replay");
    std::ofstream file;
    if (!output_path.empty()) {
        file.open(output_path);
        if (!file) throw std::runtime_error("cannot open frequency basis report");
    }
    auto &report = file.is_open() ? file : std::cout;
    report.precision(17);
    if (storage_only) {
        constexpr auto factor_fft = CKKSPreparedAutoKeyFFT<P, Narrow::factor_key_log_q>::spectrum_count * sizeof(PolynomialInFD<P>);
        constexpr auto tree_fft = []<std::size_t... I>(std::index_sequence<I...>) {
            return (std::size_t{} + ... + sizeof(typename std::tuple_element_t<I,
                typename Narrow::PreparedRelinKeys::Tuple>::element_type));
        }(std::make_index_sequence<Narrow::tree_depth>{});
        using Cache = CKKSNarrowNTTSelector<Narrow>;
        report << "{\n  \"scope\": \"Compiled public payloads; excludes workspaces, allocator overhead and process reserve\",\n"
               << "  \"n\": " << P::n << ", \"storage_bits\": " << std::numeric_limits<typename P::T>::digits
               << ", \"selector_bits\": " << selector_bits << ", \"primary_bits\": " << P::Bgbit << ", \"secondary_bits\": " << P::B̅gbit << ",\n"
               << "  \"output_bits\": " << Basis::output_log_q << ", \"output_scale_bits\": " << Basis::log_delta << ",\n"
               << "  \"compact_setup_bytes\": " << sizeof(CKKSNarrowRefreshKey<Refresh>)
               << ",\n  \"compact_s2c_polynomial_bytes\": " << P::n * sizeof(S2CTerm)
               << ",\n  \"basis_public_key_bytes\": " << sizeof(CKKSNarrowFrequencyBasisKey<Basis>)
               << ",\n  \"retained_dd_key_bytes\": " << sizeof(CKKSCompactDDRelinKey<P, Basis::output_log_q>)
               << ",\n  \"expanded_selector_bytes\": " << Narrow::selector_count * sizeof(typename Narrow::Selector)
               << ",\n  \"one_selector_group_bytes\": " << Narrow::block_width * sizeof(typename Narrow::Selector)
               << ",\n  \"selector_ntt_chunk_bytes\": " << std::min<std::size_t>(40, Narrow::block_width) *
                    Cache::prime_count * (P::k + 1) * Cache::digits * P::n * sizeof(std::uint64_t)
               << ",\n  \"factor_fft_bytes\": " << factor_fft
               << ",\n  \"bounded_rotation_fft_bytes\": " << factor_fft
               << ",\n  \"dd_tree_fft_bytes\": " << tree_fft
               << ",\n  \"basis_product_fft_bytes\": " << sizeof(CKKSDDRelinKeyFFT<P, Basis::root_log_q>)
               << ",\n  \"basis_output_fft_bytes\": " << CKKSPreparedAutoKeyFFT<P, Basis::output_log_q>::spectrum_count * sizeof(PolynomialInFD<P>)
               << ",\n  \"retained_dd_fft_bytes\": " << sizeof(CKKSDDRelinKeyFFT<P, Basis::output_log_q>)
               << ",\n  \"maximum_frequency_columns\": 22"
               << ",\n  \"factor_digit_workspace_bytes\": " << Cache::digits * 22 * 2 * (P::k + 1) * P::n * sizeof(typename P::T)
               << ",\n  \"factor_modular_sum_bytes\": " << 2 * (P::k + 1) * Cache::digits * Cache::prime_count * 22 * P::n * sizeof(std::uint64_t)
               << ",\n  \"head_tail_factor_capacity_bytes\": " << 2 * (Narrow::blocks + 1) * 22 * sizeof(typename Narrow::Factor)
               << ",\n  \"full_validation_passed\": false\n}\n";
        return 0;
    }
    if (function_option != "sigmoid" && function_option != "gelu" &&
        function_option != "tanh" && function_option != "all")
        throw std::invalid_argument("unsupported frequency basis function");
    if (!checkpoint_directory.empty() && output_path.empty())
        throw std::invalid_argument("basis checkpoint evaluation requires a report path");
    if (directory.empty()) throw std::invalid_argument("basis benchmark requires a new selector directory");
    std::mt19937_64 engine(20262119);
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    const CKKSNoise noise{P::α, 0};
    const auto start = std::chrono::steady_clock::now();
    auto key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    auto basis_key = std::make_unique<CKKSNarrowFrequencyBasisKey<Basis>>();
    auto retained = std::make_unique<CKKSCompactDDRelinKey<P, Basis::output_log_q>>();
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto input = std::make_unique<typename Refresh::Input>();
    const std::array<std::uint64_t, 12> checkpoint_metadata{
        0x4e41525242433031ULL, P::n, std::numeric_limits<typename P::T>::digits,
        selector_bits, root_delta, coefficient_bits, Refresh::input_log_q,
        Refresh::input_log_delta, P::Bgbit, P::B̅gbit, Narrow::auxiliary_log_scale,
        Refresh::s2c_plain_log_delta};
    const auto checkpoint = std::filesystem::path(checkpoint_directory) / "checkpoint.bin";
    if (replay_checkpoint) {
        if (checkpoint_directory.empty()) throw std::invalid_argument("basis replay requires a checkpoint");
        std::ifstream saved(checkpoint, std::ios::binary);
        if (!saved) throw std::runtime_error("cannot open basis checkpoint");
        cereal::PortableBinaryInputArchive archive(saved);
        std::array<std::uint64_t, 12> metadata;
        archive(metadata);
        if (metadata != checkpoint_metadata) throw std::runtime_error("basis checkpoint profile mismatch");
        archive(*dense, *sparse, *slots, *input, *key, *basis_key, *retained);
        if (saved.peek() != std::char_traits<char>::eof())
            throw std::runtime_error("basis checkpoint has trailing data");
        key->narrow.selector_source = std::make_shared<CKKSNarrowSelectorFile<Narrow>>(
            directory, typename CKKSNarrowSelectorFile<Narrow>::Existing{});
        narrow_detail::structuredPositions<Narrow>(*sparse);
        std::cout << "basis n=" << P::n << " reused checkpoint keys and input" << std::endl;
    }
    else {
        if (!checkpoint_directory.empty()) {
            if (!std::filesystem::create_directory(checkpoint_directory))
                throw std::runtime_error("basis checkpoint directory must not exist");
            std::filesystem::permissions(checkpoint_directory, std::filesystem::perms::owner_all,
                std::filesystem::perm_options::replace);
        }
        // Public benchmark inputs may be reproducible; secret keys must use
        // the library's key sampler rather than the fixed input seed.
        keyGen<P>(*dense);
        CKKSNarrowStructuredKeyGen<Narrow>(*sparse);
        std::cout << "basis n=" << P::n << " generating public setup and file selectors" << std::endl;
        CKKSNarrowRefreshSetupKeyGen<Refresh>(*key, *dense, *sparse, *dense, noise, CKKSNarrowCacheMode::Streamed);
        CKKSNarrowEncryptSelectorsToFile<Narrow>(key->narrow, *sparse, *dense, directory, noise);
        CKKSNarrowFrequencyBasisKeyGen<Basis>(*basis_key, *dense, noise);
        // All functions use the same secret and output level. Generate this
        // validation key once, and retain it with the reusable public setup.
        CKKSCompactDDRelinKeyGen(*retained, *dense, noise);
        std::uniform_real_distribution<double> message(-4, 4);
        for (auto &slot : *slots) slot = {message(engine), message(engine)};
        (*slots)[0] = {-4, 4}; (*slots)[1] = {4, -4}; (*slots)[2] = {4, 4}; (*slots)[3] = {0, 0};
        ckksSlotEncrypt<P>(*input, *slots, *dense, noise);
        if (!checkpoint_directory.empty()) {
            const auto pending = std::filesystem::path(checkpoint_directory) / "checkpoint.pending";
            {
                std::ofstream saved(pending, std::ios::binary);
                if (!saved) throw std::runtime_error("cannot create basis checkpoint");
                { cereal::PortableBinaryOutputArchive archive(saved);
                  archive(checkpoint_metadata, *dense, *sparse, *slots, *input, *key, *basis_key, *retained); }
                saved.flush();
                if (!saved) throw std::runtime_error("cannot flush basis checkpoint");
            }
            auto source = std::dynamic_pointer_cast<const CKKSNarrowSelectorFile<Narrow>>(key->narrow.selector_source);
            if (!source) throw std::logic_error("basis checkpoint requires file selectors");
            source->preserve();
            std::filesystem::rename(pending, checkpoint);
        }
    }
    if (!prepared_selector_directory.empty()) {
#if defined(USE_BLAKE3)
        auto raw = std::dynamic_pointer_cast<const CKKSNarrowSelectorFile<Narrow>>(key->narrow.selector_source);
        if (!raw) throw std::invalid_argument("prepared cache requires file-backed raw selectors");
        key->narrow.selector_source = std::make_shared<CKKSNarrowNTTSelectorFile<Narrow>>(
            prepared_selector_directory, std::move(raw), typename CKKSNarrowNTTSelectorFile<Narrow>::Existing{});
#else
        throw std::invalid_argument("prepared selector files require BLAKE3 support");
#endif
    }
    CKKSNarrowUseStagedFFT(key->narrow, 40);
    const auto setup_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "basis setup seconds=" << setup_seconds << "; evaluating full complex FBS" << std::endl;
    auto retained_fft = makeCKKSDDRelinKeyFFT(*retained);
    const std::vector<std::string> functions = function_option == "all"
        ? std::vector<std::string>{"sigmoid", "gelu", "tanh"}
        : std::vector<std::string>{function_option};
    std::unique_ptr<std::array<typename Refresh::Coefficients, 2>> prepared_halves;
    std::array<std::map<std::uint32_t, typename Basis::Root>, 2> root_cache;
    const bool share_roots = functions.size() > 1;
    const auto make_plan = [](const std::string &function) {
        return function == "tanh"
            ? CKKSNarrowFrequencyBasis(40, {1, 3, 4, 9, 11, 16, 17, 19, 20})
            : CKKSNarrowFrequencyBasis(12, {1, 3, 5, 6});
    };
    std::vector<CKKSNarrowFrequencyBasis> shared_plans;
    if (share_roots)
        for (const auto &function : functions) shared_plans.push_back(make_plan(function));
    const auto shared_columns = 2 * CKKSNarrowFrequencyBasisUnion(shared_plans).size();
    double coefficient_prepare_seconds = 0;
    std::size_t shared_s2c_term_payload_bytes_released = 0;
    if (share_roots) {
        prepared_halves = std::make_unique<std::array<typename Refresh::Coefficients, 2>>();
        const auto prepare_start = std::chrono::steady_clock::now();
        CKKSNarrowPrepareCoefficients<Refresh>((*prepared_halves)[0], *input, *key);
        narrow_detail::imaginaryCoefficientHalf<P>((*prepared_halves)[1], (*prepared_halves)[0]);
        coefficient_prepare_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-prepare_start).count();
        // All functions consume the prepared coefficient halves below. The
        // retained checkpoint already contains the reusable S2C plan.
        for (const auto &group : key->s2c_plan.direct.groups)
            shared_s2c_term_payload_bytes_released += group.terms.size() * sizeof(group.terms[0]);
        for (const auto &group : key->s2c_plan.conjugate.groups)
            shared_s2c_term_payload_bytes_released += group.terms.size() * sizeof(group.terms[0]);
        key->s2c_plan = {};
    }
    bool all_passed = true;
    bool first_result = true;
    if (function_option == "all")
        report << "{\n  \"shared_public_setup\": true, \"shared_encrypted_input\": true,\n"
               << "  \"setup_generation_count\": " << (replay_checkpoint ? 0 : 1) << ", \"results\": [\n";
    for (const auto &function : functions) {
        const auto target = [&](double x) {
            if (function == "sigmoid") return 1 / (1 + std::exp(-x));
            if (function == "gelu") return 0.5 * x * (1 + std::erf(x / std::sqrt(2.0)));
            return std::tanh(x);
        };
        const auto curvature = [&](double x) {
            const auto y = target(x);
            if (function == "sigmoid") return y * (1 - y) * (1 - 2 * y);
            if (function == "gelu") return (2 - x * x) * std::exp(-x * x / 2) / std::sqrt(2 * std::acos(-1.0));
            return -2 * y * (1 - y * y);
        };
        const bool windowed = function == "tanh";
        const std::uint32_t degree = windowed ? 40 : 12;
        const auto approximation = windowed ? CKKSBuildNarrowWindowedTrigonometricApproximation(
            target, [](double x) { const auto y = std::tanh(x); return 1 - y * y; }, curvature,
            {1, 1, 2, 4, 40}, 4, 16, degree) : CKKSBuildNarrowTrigonometricApproximation(
            target, 4, 16, degree, function == "gelu" ? 2.0 : 1.0, 65536, 40, 2.0, curvature, 160.0);
        const auto plan = make_plan(function);
        std::cout << "basis evaluating " << function << " with shared public setup" << std::endl;
        auto output = std::make_unique<typename Basis::Output>();
        const auto refresh_start = std::chrono::steady_clock::now();
        std::size_t newly_evaluated_roots = 0, reused_roots = 0;
        const auto series_tile_size = share_roots ? shared_columns : 2 * plan.frequencies().size();
        std::size_t factor_groups_completed = 0;
        CKKSNarrowFactorTimings factor_totals;
        CKKSNarrowRootTimings root_totals;
        double root_generation_seconds = 0, basis_reconstruction_seconds = 0;
        key->narrow.root_stage_observer = [&](const CKKSNarrowRootTimings &times) {
            root_totals.input_and_factor_cache_seconds += times.input_and_factor_cache_seconds;
            root_totals.leaf_factors_and_inline_alignment_seconds += times.leaf_factors_and_inline_alignment_seconds;
            root_totals.bounded_tail_alignment_seconds += times.bounded_tail_alignment_seconds;
            root_totals.dd_cache_preparation_seconds += times.dd_cache_preparation_seconds;
            root_totals.dd_trees_and_beta_seconds += times.dd_trees_and_beta_seconds;
        };
        key->narrow.factor_observer = [&](bool tail, std::size_t group, const CKKSNarrowFactorTimings &times) {
            ++factor_groups_completed;
            factor_totals.selector_load_seconds += times.selector_load_seconds;
            factor_totals.table_ntt_seconds += times.table_ntt_seconds;
            factor_totals.selector_ntt_seconds += times.selector_ntt_seconds;
            factor_totals.matrix_seconds += times.matrix_seconds;
            factor_totals.inverse_crt_seconds += times.inverse_crt_seconds;
            factor_totals.conjugation_seconds += times.conjugation_seconds;
            std::cout << "basis " << function << " completed " << (tail ? "tail" : "head")
                      << " group=" << group << " factor_groups=" << factor_groups_completed
                      << " elapsed_seconds=" << std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - refresh_start).count() << std::endl;
        };
        if (share_roots) {
            if (root_cache[0].empty()) {
                std::vector<std::map<std::uint32_t, typename Basis::Root>> shared_roots;
                const auto root_start = std::chrono::steady_clock::now();
                CKKSNarrowBlindRotateBasisPlansBatch<Narrow, Refresh::input_log_delta>(
                    shared_roots, *prepared_halves, shared_plans, key->narrow, series_tile_size);
                root_generation_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - root_start).count();
                for (std::size_t half = 0; half < 2; ++half) {
                    newly_evaluated_roots += shared_roots[half].size();
                    root_cache[half] = std::move(shared_roots[half]);
                }
            }
            else reused_roots = 2 * plan.frequencies().size();
            const auto reconstruction_start = std::chrono::steady_clock::now();
            auto halves = std::make_unique<std::array<typename Basis::Output, 2>>();
            for (std::size_t half = 0; half < 2; half++) {
                std::vector<typename Basis::Root> roots;
                for (const auto frequency : plan.frequencies()) roots.push_back(root_cache[half].at(frequency));
                CKKSNarrowEvalSeriesFromBasisRoots<Basis>((*halves)[half], roots, *basis_key,
                    plan, approximation.coefficients);
            }
            narrow_detail::combineRealImaginary<P>(*output, (*halves)[0], (*halves)[1]);
            basis_reconstruction_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - reconstruction_start).count();
        }
        else CKKSNarrowFunctionalRefreshWithBasis<Refresh, coefficient_bits>(*output, *input, *key, *basis_key, plan, approximation, true, series_tile_size);
        key->narrow.factor_observer = {};
        key->narrow.root_stage_observer = {};
        const auto refresh_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - refresh_start).count();
        auto decoded = std::make_unique<CKKSSlotVector<P>>(), expected = std::make_unique<CKKSSlotVector<P>>(), approximated = std::make_unique<CKKSSlotVector<P>>();
        ckksSlotDecrypt<P>(*decoded, *output, *dense);
        double error = 0, evaluation_error = 0, amplitude = 0, approximation_amplitude = 0;
        for (std::size_t i = 0; i < P::n / 2; i++) {
            (*expected)[i] = {target((*slots)[i].real()), target((*slots)[i].imag())};
            (*approximated)[i] = {CKKSEvaluateNarrowApproximation(approximation, (*slots)[i].real()),
                CKKSEvaluateNarrowApproximation(approximation, (*slots)[i].imag())};
            error = finiteMaximum(error, std::abs((*decoded)[i] - (*expected)[i]));
            evaluation_error = finiteMaximum(evaluation_error, std::abs((*decoded)[i] - (*approximated)[i]));
            amplitude = finiteMaximum(amplitude, std::abs((*expected)[i]));
            approximation_amplitude = finiteMaximum(approximation_amplitude, std::abs((*approximated)[i]));
        }
        std::cout << "basis refresh seconds=" << refresh_seconds << " target error=" << error << std::endl;
        auto square = std::make_unique<Square>();
        CKKSMultRelinThenRescale<P>(*square, *output, *output, *retained_fft);
        ckksSlotDecrypt<P>(*decoded, *square, *dense);
        double product_error = 0, product_evaluation_error = 0;
        for (std::size_t i = 0; i < P::n / 2; i++) {
            product_error = finiteMaximum(product_error, std::abs((*decoded)[i] - (*expected)[i] * (*expected)[i]));
            product_evaluation_error = finiteMaximum(product_evaluation_error, std::abs((*decoded)[i] - (*approximated)[i] * (*approximated)[i]));
        }
        auto normalized = std::make_unique<typename Refresh::Input>();
        CKKSScaleDown<P, Square::log_q, Square::log_delta, Basis::log_delta - Refresh::input_log_delta>(*normalized, *square);
        auto normalized_slots = std::make_unique<CKKSSlotVector<P>>();
        ckksSlotDecrypt<P>(*normalized_slots, *normalized, *dense);
        double normalization_error = 0;
        for (std::size_t i = 0; i < P::n / 2; i++)
            normalization_error = finiteMaximum(normalization_error, std::abs((*normalized_slots)[i] - (*decoded)[i]));
        const double tolerance = std::ldexp(1.0, -23);
        const double approximation_bound = std::sqrt(2.0) * *approximation.uniform_error_bound;
        const double product_tolerance = 2 * amplitude * tolerance + tolerance * tolerance + tolerance;
        const double product_evaluation_tolerance = 2 * approximation_amplitude * tolerance + tolerance * tolerance + tolerance;
        const bool passed = error < tolerance && evaluation_error < tolerance && approximation_bound < tolerance &&
            product_error < product_tolerance && product_evaluation_error < product_evaluation_tolerance && normalization_error < tolerance;
        if (!checkpoint_directory.empty()) {
            const auto result_path = std::filesystem::path(output_path).parent_path() /
                (function + (replay_checkpoint ? "-replay.bin" : "-initial.bin"));
            std::ofstream saved(result_path, std::ios::binary);
            if (!saved) throw std::runtime_error("cannot save basis validation output");
            { cereal::PortableBinaryOutputArchive archive(saved);
              archive(checkpoint_metadata, function, *output, *square, *normalized); }
            saved.flush();
            if (!saved) throw std::runtime_error("cannot flush basis validation output");
        }
        if (function_option == "all" && !first_result) report << ",\n";
        first_result = false;
        report << "{\n  \"profile\": \"frequency-basis-staged-dd\", \"n\": " << P::n
               << ", \"h\": 31, \"gaussian_sigma\": 3.2, \"storage_bits\": " << std::numeric_limits<typename P::T>::digits << ",\n"
               << "  \"function\": \"" << function << "\", \"degree\": " << degree << ", \"basis_size\": " << plan.frequencies().size()
               << ",\n  \"source_log_q\": 48, \"input_log_q\": 85, \"input_log_delta\": 44, \"selector_log_q\": " << selector_bits << ",\n"
               << "  \"root_plain_log_delta\": 39, \"root_log_delta\": " << root_delta << ", \"coefficient_bits\": " << coefficient_bits << ",\n"
               << "  \"output_log_q\": " << Basis::output_log_q << ", \"output_log_delta\": " << Basis::log_delta
               << ", \"gadget_bits\": " << P::Bgbit << ", \"secondary_digit_bits\": " << P::B̅gbit << ", \"auxiliary_switch_bits\": 44,\n"
               << "  \"derived_tree_dd_secondary_bits\": " << Narrow::prepared_relin_bbarbit
               << ", \"secondary11_replay\": " << (Secondary11 ? "true" : "false") << ",\n"
               << "  \"selector_matrix_digit_bits\": " << CKKSNarrowNTTSelector<Narrow>::digit_bits
               << ", \"selector_matrix_digits\": " << CKKSNarrowNTTSelector<Narrow>::digits
               << ", \"selector56_replay\": " << (Selector56 ? "true" : "false") << ",\n"
               << "  \"staged_fft\": true, \"selector_chunk_size_requested\": 40, \"selector_chunk_size\": "
               << key->narrow.selector_ntt_chunk_size << ", \"series_tile_size\": " << series_tile_size << ", \"message_bound\": 4,\n"
               << "  \"target_tolerance\": " << tolerance << ", \"approximation_uniform_bound\": " << approximation_bound
               << ",\n  \"max_error\": " << error << ", \"evaluation_max_error\": " << evaluation_error
               << ",\n  \"post_refresh_product_max_error\": " << product_error << ", \"post_refresh_product_evaluation_max_error\": " << product_evaluation_error
               << ",\n  \"product_tolerance\": " << product_tolerance << ", \"product_evaluation_tolerance\": " << product_evaluation_tolerance
               << ",\n  \"square_normalization_max_error\": " << normalization_error
               << ",\n  \"square_normalization_scope\": \"Dense ciphertext scale/modulus adapter only; a repeated nonlinear refresh also requires its message-domain bound\",\n"
               << "  \"shared_frequency_roots\": " << (share_roots ? "true" : "false")
               << ", \"newly_evaluated_roots\": " << newly_evaluated_roots << ", \"reused_roots\": " << reused_roots
               << ", \"factor_groups_completed\": " << factor_groups_completed
               << ", \"root_stage_timings_available\": " << (share_roots ? "true" : "false")
               << ", \"root_generation_seconds\": " << root_generation_seconds
               << ", \"basis_reconstruction_seconds\": " << basis_reconstruction_seconds
               << ", \"root_tile_stage_seconds\": {\"input_and_factor_cache\": " << root_totals.input_and_factor_cache_seconds
               << ", \"leaf_factors_and_inline_alignment\": " << root_totals.leaf_factors_and_inline_alignment_seconds
               << ", \"bounded_tail_alignment\": " << root_totals.bounded_tail_alignment_seconds
               << ", \"dd_cache_preparation\": " << root_totals.dd_cache_preparation_seconds
               << ", \"dd_trees_and_beta\": " << root_totals.dd_trees_and_beta_seconds << "}"
               << ", \"factor_stage_seconds\": {\"selector_load\": " << factor_totals.selector_load_seconds
               << ", \"table_ntt\": " << factor_totals.table_ntt_seconds
               << ", \"selector_ntt\": " << factor_totals.selector_ntt_seconds
               << ", \"matrix\": " << factor_totals.matrix_seconds
               << ", \"inverse_crt\": " << factor_totals.inverse_crt_seconds
               << ", \"conjugation\": " << factor_totals.conjugation_seconds << "}"
               << ", \"shared_coefficient_prepare_seconds\": " << coefficient_prepare_seconds << ",\n"
               << "  \"shared_s2c_term_payload_bytes_released\": " << shared_s2c_term_payload_bytes_released << ",\n"
               << "  \"key_setup_seconds\": " << setup_seconds << ", \"refresh_seconds\": " << refresh_seconds
               << ",\n  \"fresh_keys_generated\": " << (replay_checkpoint ? "false" : "true")
               << ", \"input_reused\": " << (replay_checkpoint ? "true" : "false")
               << ", \"retained_square_key_shared\": true"
               << ", \"shared_basis_union\": " << (share_roots ? "true" : "false")
               << ", \"refresh_timing_scope\": \"" << (share_roots
                    ? (newly_evaluated_roots ? "shared all-function root generation and this function reconstruction"
                                             : "this function reconstruction from shared roots")
                    : "standalone function refresh") << "\""
               << ",\n  \"accuracy_passed\": " << (passed ? "true" : "false") << ", \"security_validated\": false\n}\n";
        report.flush();
        all_passed &= passed;
    }
    if (function_option == "all")
        report << "\n  ], \"accuracy_passed\": " << (all_passed ? "true" : "false")
               << ", \"security_validated\": false\n}\n";
    return all_passed ? 0 : 1;
}
catch (const std::exception &error) {
    std::cerr << "frequency basis benchmark: " << error.what() << '\n';
    return 2;
}


}  // namespace narrow_benchmark
