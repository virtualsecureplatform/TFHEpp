#include <algorithm>
#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tfhe++.hpp>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

constexpr const char *validation_test_key_metadata_filename =
    ".ckks_bootstrap_validation_key";

struct TinyDeepMultiLimbCKKSParam {
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
    using T = TFHEpp::MultiLimbUInt<9>;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bg_a = T{1} << Bg_abit;
    static constexpr T Bgₐ = Bg_a;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = 0.0;
    static constexpr T mu = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr T μ = mu;
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Delta = 0.0;
    static constexpr double Δ = Delta;
    static constexpr std::uint32_t lbar = 36;
    static constexpr std::uint32_t l̅ = lbar;
    static constexpr std::uint32_t lbar_a = 36;
    static constexpr std::uint32_t l̅ₐ = lbar_a;
    static constexpr std::uint32_t Bbarbit = 16;
    static constexpr std::uint32_t B̅gbit = Bbarbit;
    static constexpr std::uint32_t Bbar_abit = 16;
    static constexpr std::uint32_t B̅gₐbit = Bbar_abit;
};

template <class Param>
void fill_test_key(TFHEpp::Key<Param> &key)
{
    for (std::size_t i = 0; i < Param::n; i++) {
        const int v = static_cast<int>(i % 3) - 1;
        key[i] = static_cast<typename Param::T>(v);
    }
}

template <class Param>
void fill_sparse_test_key(TFHEpp::Key<Param> &key, std::size_t weight)
{
    key.fill(typename Param::T{0});
    if (weight == 0) return;

    constexpr std::size_t stride = 7919;
    for (std::size_t i = 0; i < weight; i++) {
        const std::size_t index = (i * stride) % Param::n;
        const int sign = (i % 2 == 0) ? 1 : -1;
        key[index] = static_cast<typename Param::T>(sign);
    }
}

template <class Param>
void fill_bootstrap_test_key(TFHEpp::Key<Param> &key, std::size_t sparse_weight)
{
    if (sparse_weight == 0)
        fill_test_key<Param>(key);
    else
        fill_sparse_test_key<Param>(key, sparse_weight);
}

template <class Schedule>
constexpr bool requires_bounded_modraise_test_key()
{
    using P = typename Schedule::Param;
    return P::n >= (1U << 15);
}

template <class Schedule>
constexpr std::size_t bounded_modraise_sparse_weight_max()
{
    return Schedule::bounded_sparse_secret_key_weight;
}

template <class Schedule>
bool sparse_weight_fits_bounded_modraise(std::size_t sparse_weight)
{
    return sparse_weight != 0 &&
           sparse_weight <= bounded_modraise_sparse_weight_max<Schedule>();
}

template <class Schedule>
int validate_bounded_modraise_test_key(std::size_t sparse_weight,
                                       const char *label)
{
    if constexpr (!requires_bounded_modraise_test_key<Schedule>()) {
        return 0;
    }
    else {
        constexpr std::size_t max_weight =
            bounded_modraise_sparse_weight_max<Schedule>();
        if (sparse_weight == 0) {
            std::cerr << label
                      << "_dense_key_unbounded_for_evalmod=1 evalmod_k="
                      << Schedule::evalmod_k
                      << " modraise_mask_bound="
                      << Schedule::modraise_mask_bound
                      << " bounded_sparse_weight_max=" << max_weight << '\n';
            return 2;
        }
        if (sparse_weight > max_weight) {
            std::cerr << label
                      << "_sparse_weight_out_of_evalmod_range sparse_weight="
                      << sparse_weight << " evalmod_k=" << Schedule::evalmod_k
                      << " evalmod_mask_bound="
                      << Schedule::evalmod_mask_bound
                      << " modraise_mask_bound="
                      << Schedule::modraise_mask_bound
                      << " bounded_sparse_weight_max=" << max_weight << '\n';
            return 2;
        }
        return 0;
    }
}

template <class Param>
void fill_test_slots(TFHEpp::CKKSSlotVector<Param> &slots)
{
    for (std::size_t i = 0; i < Param::n / 2; i++) {
        const double re =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 64.0;
        const double im =
            static_cast<double>(static_cast<int>(i % 7) - 3) / 128.0;
        slots[i] = {re, im};
    }
}

template <class Param>
void fill_alternate_test_slots(TFHEpp::CKKSSlotVector<Param> &slots)
{
    for (std::size_t i = 0; i < Param::n / 2; i++) {
        const double re =
            static_cast<double>(static_cast<int>((3 * i + 1) % 7) - 3) / 96.0;
        const double im =
            static_cast<double>(static_cast<int>((5 * i + 2) % 9) - 4) / 160.0;
        slots[i] = {re, im};
    }
}

template <class Param>
void multiply_slots(TFHEpp::CKKSSlotVector<Param> &out,
                    const TFHEpp::CKKSSlotVector<Param> &lhs,
                    const TFHEpp::CKKSSlotVector<Param> &rhs)
{
    for (std::size_t i = 0; i < Param::n / 2; i++) out[i] = lhs[i] * rhs[i];
}

template <class Param>
double max_error(const TFHEpp::CKKSSlotVector<Param> &got,
                 const TFHEpp::CKKSSlotVector<Param> &want)
{
    double err = 0.0;
    for (std::size_t i = 0; i < Param::n / 2; i++)
        err = std::max(err, std::abs(got[i] - want[i]));
    return err;
}

template <class Param>
double max_abs(const TFHEpp::CKKSSlotVector<Param> &slots)
{
    double value = 0.0;
    for (const auto &slot : slots) value = std::max(value, std::abs(slot));
    return value;
}

template <class Param>
void print_slot_diagnostic(const char *label,
                           const TFHEpp::CKKSSlotVector<Param> &got,
                           const TFHEpp::CKKSSlotVector<Param> *want = nullptr)
{
    std::cout << std::setprecision(17);
    std::cout << label << "_max_abs=" << max_abs<Param>(got);
    if (want != nullptr)
        std::cout << " " << label
                  << "_max_error=" << max_error<Param>(got, *want)
                  << " " << label << "_want_max_abs=" << max_abs<Param>(*want);
    std::cout << " " << label << "_slot0=(" << got[0].real() << ","
              << got[0].imag() << ")";
    if (want != nullptr)
        std::cout << " " << label << "_want0=(" << (*want)[0].real() << ","
                  << (*want)[0].imag() << ")";
    std::cout << '\n';
}

template <class Schedule, std::uint32_t LogQ>
void print_phase_coefficient_diagnostic(
    const char *label,
    const TFHEpp::CKKSCiphertext<typename Schedule::Param, LogQ,
                                 Schedule::log_delta> &ct,
    const TFHEpp::Key<typename Schedule::Param> &key)
{
    using P = typename Schedule::Param;
    TFHEpp::Polynomial<P> phase = TFHEpp::trlwePhase<P>(ct.ct, key);
    TFHEpp::ckks_detail::reducePolynomialToLevel<P, LogQ>(phase);

    double max_coeff = 0.0;
    double max_normalized = 0.0;
    double coeff0 = 0.0;
    double normalized0 = 0.0;
    for (std::uint32_t i = 0; i < P::n; i++) {
        const long double coeff = std::ldexp(
            TFHEpp::ckks_detail::levelToLongDouble<P, LogQ>(phase[i]),
            -static_cast<int>(Schedule::log_delta));
        const double coeff_double = static_cast<double>(coeff);
        const double normalized =
            coeff_double /
            (static_cast<double>(Schedule::evalmod_k) *
             Schedule::message_ratio);
        if (i == 0) {
            coeff0 = coeff_double;
            normalized0 = normalized;
        }
        max_coeff = std::max(max_coeff, std::abs(coeff_double));
        max_normalized = std::max(max_normalized, std::abs(normalized));
    }
    std::cout << std::setprecision(17);
    std::cout << label << "_phase_coeff_max_abs=" << max_coeff << " "
              << label << "_phase_normalized_max_abs=" << max_normalized
              << " " << label << "_phase_coeff0=" << coeff0 << " " << label
              << "_phase_normalized0=" << normalized0 << '\n';
}

template <class Param>
void apply_complex_linear(
    TFHEpp::CKKSSlotVector<Param> &out,
    const TFHEpp::CKKSSlotVector<Param> &in,
    const std::vector<TFHEpp::CKKSSlotVector<Param>> &diagonals,
    const std::vector<int> &offsets)
{
    constexpr int half = static_cast<int>(Param::n) / 2;
    out.fill({0.0, 0.0});
    for (std::size_t d = 0; d < diagonals.size(); d++) {
        const int offset = ((offsets[d] % half) + half) % half;
        for (int i = 0; i < half; i++)
            out[i] += diagonals[d][i] * in[(i + offset) % half];
    }
}

template <class Param>
void apply_complex_stages(
    TFHEpp::CKKSSlotVector<Param> &out,
    const TFHEpp::CKKSSlotVector<Param> &in,
    const TFHEpp::CKKSLinearTransformStages<Param> &stages)
{
    auto cur = std::make_unique<TFHEpp::CKKSSlotVector<Param>>(in);
    auto next = std::make_unique<TFHEpp::CKKSSlotVector<Param>>();
    for (const auto &stage : stages) {
        apply_complex_linear<Param>(*next, *cur, stage.diagonals,
                                    stage.rotation_offsets);
        std::swap(cur, next);
    }
    out = *cur;
}

template <class F>
double elapsed_ms(F &&fn)
{
    const auto begin = Clock::now();
    fn();
    const auto end = Clock::now();
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

void print_bootstrap_timings(const TFHEpp::CKKSDenseBootstrapTimings &timings)
{
    std::cout << "bootstrap_normalize_ms=" << timings.normalize_ms << '\n';
    std::cout << "bootstrap_modraise_ms=" << timings.modraise_ms << '\n';
    std::cout << "bootstrap_coeff_to_slot_ms=" << timings.coeff_to_slot_ms
              << '\n';
    std::cout << "bootstrap_split_ms=" << timings.split_ms << '\n';
    std::cout << "bootstrap_real_evalmod_ms=" << timings.real_evalmod_ms
              << '\n';
    std::cout << "bootstrap_imag_evalmod_ms=" << timings.imag_evalmod_ms
              << '\n';
    std::cout << "bootstrap_slot_to_coeff_ms=" << timings.slot_to_coeff_ms
              << '\n';
    std::cout << "bootstrap_timed_total_ms=" << timings.total_ms() << '\n';
}

std::filesystem::path
validation_test_key_metadata_file(const std::filesystem::path &key_dir)
{
    return key_dir / validation_test_key_metadata_filename;
}

bool is_validation_test_key_metadata_file(const std::filesystem::path &path)
{
    return path.filename() == validation_test_key_metadata_filename;
}

bool read_validation_test_key_sparse_weight(
    const std::filesystem::path &key_dir, std::size_t &sparse_weight)
{
    std::ifstream input(validation_test_key_metadata_file(key_dir));
    std::string name;
    return static_cast<bool>(input >> name >> sparse_weight) &&
           name == "sparse_weight";
}

int write_validation_test_key_sparse_weight(
    const std::filesystem::path &key_dir, std::size_t sparse_weight)
{
    std::filesystem::create_directories(key_dir);
    std::ofstream output(validation_test_key_metadata_file(key_dir),
                         std::ios::trunc);
    if (!output) {
        std::cerr << "test_key_metadata_write_failed="
                  << validation_test_key_metadata_file(key_dir).string()
                  << '\n';
        return 2;
    }
    output << "sparse_weight " << sparse_weight << '\n';
    return output ? 0 : 2;
}

std::uintmax_t directory_size_bytes(const std::filesystem::path &root)
{
    if (!std::filesystem::exists(root)) return 0;
    std::uintmax_t total = 0;
    for (const auto &entry : std::filesystem::directory_iterator(root)) {
        if (entry.is_regular_file()) total += entry.file_size();
    }
    return total;
}

std::size_t regular_file_count(const std::filesystem::path &root)
{
    if (!std::filesystem::exists(root)) return 0;
    std::size_t count = 0;
    for (const auto &entry : std::filesystem::directory_iterator(root)) {
        if (entry.is_regular_file() &&
            !is_validation_test_key_metadata_file(entry.path()))
            count++;
    }
    return count;
}

std::size_t temporary_file_count(const std::filesystem::path &root)
{
    if (!std::filesystem::exists(root)) return 0;
    std::size_t count = 0;
    for (const auto &entry : std::filesystem::directory_iterator(root)) {
        if (entry.is_regular_file() && entry.path().extension() == ".tmp")
            count++;
    }
    return count;
}

template <class Schedule>
int validate_or_create_validation_test_key_metadata(
    const std::filesystem::path &key_dir, std::size_t sparse_weight,
    bool allow_create)
{
    if constexpr (!requires_bounded_modraise_test_key<Schedule>()) {
        (void)key_dir;
        (void)sparse_weight;
        (void)allow_create;
        return 0;
    }
    else {
        const std::filesystem::path metadata =
            validation_test_key_metadata_file(key_dir);
        if (std::filesystem::exists(metadata)) {
            std::size_t stored_sparse_weight = 0;
            if (!read_validation_test_key_sparse_weight(
                    key_dir, stored_sparse_weight)) {
                std::cerr << "test_key_metadata_unreadable="
                          << metadata.string() << '\n';
                return 2;
            }
            if (stored_sparse_weight != sparse_weight) {
                std::cerr << "test_key_metadata_sparse_weight_mismatch="
                          << metadata.string()
                          << " stored_sparse_weight=" << stored_sparse_weight
                          << " requested_sparse_weight=" << sparse_weight
                          << '\n';
                return 2;
            }
            return 0;
        }

        if (!allow_create) {
            std::cerr << "test_key_metadata_missing=" << metadata.string()
                      << " requested_sparse_weight=" << sparse_weight
                      << '\n';
            return 2;
        }

        const std::size_t existing_files = regular_file_count(key_dir);
        if (existing_files != 0) {
            std::cerr << "test_key_metadata_missing=" << metadata.string()
                      << " existing_key_files=" << existing_files
                      << " requested_sparse_weight=" << sparse_weight << '\n';
            return 2;
        }

        return write_validation_test_key_sparse_weight(key_dir, sparse_weight);
    }
}

void print_missing_key_files(
    const std::vector<std::filesystem::path> &missing, std::size_t limit = 8)
{
    const std::size_t shown = std::min(limit, missing.size());
    for (std::size_t i = 0; i < shown; i++)
        std::cerr << "missing_key_file=" << missing[i].string() << '\n';
    if (missing.size() > shown)
        std::cerr << "missing_key_file_more=" << missing.size() - shown
                  << '\n';
}

bool path_list_contains(const std::vector<std::filesystem::path> &paths,
                        const std::filesystem::path &path)
{
    return std::find(paths.begin(), paths.end(), path) != paths.end();
}

void print_created_key_files(
    const std::vector<std::filesystem::path> &before_missing,
    const std::vector<std::filesystem::path> &after_missing,
    std::size_t limit = 8)
{
    std::size_t count = 0;
    for (const std::filesystem::path &path : before_missing) {
        if (!path_list_contains(after_missing, path)) count++;
    }
    std::cout << "created_files=" << count << '\n';

    std::size_t shown = 0;
    for (const std::filesystem::path &path : before_missing) {
        if (path_list_contains(after_missing, path)) continue;
        if (shown < limit) std::cout << "created_file=" << path.string() << '\n';
        shown++;
    }
    if (shown > limit) std::cout << "created_file_more=" << shown - limit << '\n';
}

template <class Schedule>
std::string manifest_status(const std::filesystem::path &root)
{
    const std::filesystem::path manifest_path =
        TFHEpp::CKKSDenseBootstrapKeyDirectoryManifestFile(root);
    if (!std::filesystem::exists(manifest_path)) return "missing";
    try {
        if (TFHEpp::CKKSDenseBootstrapKeyDirectoryManifestMatches<Schedule>(
                root))
            return "sparse-match";
        if (TFHEpp::CKKSDenseBootstrapHybridGiantKeyDirectoryManifestMatches<
                Schedule>(root))
            return "hybrid-match";
        return "mismatch";
    }
    catch (...) {
        return "unreadable";
    }
}

template <int HybridThreshold>
using Lvl6HybridThresholdSchedule = TFHEpp::CKKSDenseBootstrapSchedule<
    TFHEpp::lvl6param, 50, 8, 880, 50, 5, 30, 16, 3, 0, 50, 128, 0, 50, 50,
    35, 5, 5, HybridThreshold>;

template <int HybridThreshold>
using Lvl6RobustHybridThresholdSchedule =
    TFHEpp::lvl6CKKSDenseBootstrapRobustHybridSchedule<HybridThreshold>;

using Lvl6InverseSchedule = TFHEpp::lvl6CKKSDenseBootstrapInverseSchedule;
static_assert(Lvl6InverseSchedule::evalmod_inv_degree == 3);
static_assert(Lvl6InverseSchedule::evalmod_log_q_consumption == 560);
static_assert(Lvl6InverseSchedule::output_log_q == 60);

struct Lvl6InverseBudgetParams {
    std::uint32_t log_delta = 50;
    std::uint32_t log_message_ratio = 8;
    std::uint32_t boot_log_q = 880;
    std::uint32_t coeff_to_slot_plain_log_delta = 50;
    std::uint32_t component_split_plain_log_delta = 50;
    std::uint32_t slot_to_coeff_plain_log_delta = 25;
    std::uint32_t coeff_to_slot_fuse_radix = 5;
    std::uint32_t slot_to_coeff_fuse_radix = 5;
    std::uint32_t evalmod_degree = 34;
    std::uint32_t evalmod_k = 18;
    std::uint32_t evalmod_double_angle = 3;
    std::uint32_t evalmod_inv_degree = 5;
    std::uint32_t evalmod_log_scale = 50;
};

struct Lvl6InverseBudgetResult {
    std::uint32_t c2s_levels = 0;
    std::uint32_t stc_levels = 0;
    std::uint32_t polynomial_power_depth = 0;
    std::uint32_t inverse_power_depth = 0;
    std::uint32_t coeff_rescale_count = 0;
    std::int64_t input_log_q = 0;
    std::int64_t after_c2s_log_q = 0;
    std::int64_t after_split_log_q = 0;
    std::int64_t evalmod_log_q_loss = 0;
    std::int64_t after_evalmod_log_q = 0;
    std::int64_t output_log_q = 0;
    std::int64_t output_margin_bits = 0;
    bool fits = false;
};

struct EvalModApproximationMetrics {
    double max_message_error = 0.0;
    double max_basis_error = 0.0;
    double message_bits = 0.0;
    double basis_bits = 0.0;
};

std::uint32_t budget_ceil_div(std::uint32_t value, std::uint32_t divisor)
{
    return (value + divisor - 1) / divisor;
}

std::uint32_t budget_bit_width(std::uint32_t value)
{
    std::uint32_t width = 0;
    while (value != 0) {
        width++;
        value >>= 1;
    }
    return width;
}

double positive_log2_precision(double value)
{
    if (value <= 0.0) return std::numeric_limits<double>::infinity();
    return -std::log2(value);
}

Lvl6InverseBudgetResult compute_lvl6_inverse_budget(
    const Lvl6InverseBudgetParams &params)
{
    using P = TFHEpp::lvl6param;
    Lvl6InverseBudgetResult result;
    const std::uint32_t raw_linear_stage_count = P::nbit - 1;
    result.c2s_levels =
        budget_ceil_div(raw_linear_stage_count,
                        params.coeff_to_slot_fuse_radix);
    result.stc_levels =
        budget_ceil_div(raw_linear_stage_count,
                        params.slot_to_coeff_fuse_radix);
    const std::uint32_t evalmod_degree_bound = std::max(
        params.evalmod_degree, 2 * (params.evalmod_k - 1));
    result.polynomial_power_depth =
        budget_bit_width(evalmod_degree_bound);
    result.inverse_power_depth =
        params.evalmod_inv_degree == 0
            ? 0
            : budget_bit_width(params.evalmod_inv_degree - 1);
    result.coeff_rescale_count =
        1 + (params.evalmod_inv_degree == 0 ? 0 : 1);
    result.input_log_q =
        params.log_delta + params.log_message_ratio;
    result.after_c2s_log_q =
        static_cast<std::int64_t>(params.boot_log_q) -
        static_cast<std::int64_t>(result.c2s_levels) *
            params.coeff_to_slot_plain_log_delta;
    result.after_split_log_q =
        result.after_c2s_log_q - params.component_split_plain_log_delta;
    result.evalmod_log_q_loss =
        static_cast<std::int64_t>(result.polynomial_power_depth +
                                  params.evalmod_double_angle +
                                  result.inverse_power_depth) *
            params.log_delta +
        static_cast<std::int64_t>(result.coeff_rescale_count) *
            params.evalmod_log_scale;
    result.after_evalmod_log_q =
        result.after_split_log_q - result.evalmod_log_q_loss;
    result.output_log_q =
        result.after_evalmod_log_q -
        static_cast<std::int64_t>(result.stc_levels) *
            params.slot_to_coeff_plain_log_delta;
    result.output_margin_bits =
        result.output_log_q - static_cast<std::int64_t>(params.log_delta);
    result.fits =
        result.input_log_q < params.boot_log_q &&
        params.boot_log_q <=
            static_cast<std::uint32_t>(
                std::numeric_limits<typename P::T>::digits) &&
        result.after_c2s_log_q > params.component_split_plain_log_delta &&
        result.after_split_log_q > result.evalmod_log_q_loss &&
        result.output_log_q > params.log_delta;
    return result;
}

EvalModApproximationMetrics measure_evalmod_approximation(
    std::uint32_t k, std::uint32_t degree, std::uint32_t log_message_ratio,
    std::uint32_t double_angle, std::uint32_t inverse_degree)
{
    const auto poly = TFHEpp::CKKSBuildBoundedCosEvalModPolynomial(
        k, degree, log_message_ratio, double_angle, 1.0,
        inverse_degree != 0);
    EvalModApproximationMetrics metrics;

    for (int mask = -static_cast<int>(k) + 1; mask < static_cast<int>(k);
         mask++) {
        for (int sample = 0; sample <= 256; sample++) {
            const double message =
                -1.0 + 2.0 * static_cast<double>(sample) / 256.0;
            const double masked =
                static_cast<double>(mask) * poly.message_ratio + message;
            const double cheb =
                TFHEpp::CKKSPlainEvalModBoundedCos(poly, masked,
                                                   inverse_degree);
            const double power =
                TFHEpp::CKKSPlainEvalModBoundedCosPower(poly, masked,
                                                        inverse_degree);
            metrics.max_message_error =
                std::max(metrics.max_message_error, std::abs(power - message));
            metrics.max_basis_error =
                std::max(metrics.max_basis_error, std::abs(cheb - power));
        }
    }

    metrics.message_bits =
        positive_log2_precision(metrics.max_message_error);
    metrics.basis_bits = positive_log2_precision(metrics.max_basis_error);
    return metrics;
}

void print_lvl6_inverse_budget(const char *label,
                               const Lvl6InverseBudgetParams &params)
{
    const Lvl6InverseBudgetResult result =
        compute_lvl6_inverse_budget(params);

    std::cout << label << " budget log_delta=" << params.log_delta
              << " evalmod_log_scale=" << params.evalmod_log_scale
              << " evalmod_inv_degree=" << params.evalmod_inv_degree
              << " after_split_logQ=" << result.after_split_log_q
              << " evalmod_log_q_loss=" << result.evalmod_log_q_loss
              << " output_logQ=" << result.output_log_q
              << " fits=" << (result.fits ? "yes" : "no") << '\n';
}

template <class Schedule>
void print_evalmod_approximation_report(const char *label)
{
    const EvalModApproximationMetrics metrics =
        measure_evalmod_approximation(Schedule::evalmod_k,
                                      Schedule::evalmod_degree,
                                      Schedule::log_message_ratio,
                                      Schedule::evalmod_double_angle,
                                      Schedule::evalmod_inv_degree);
    const int output_margin_bits =
        static_cast<int>(Schedule::output_log_q) -
        static_cast<int>(Schedule::log_delta);

    std::cout << label
              << " evalmod_plain_approx max_message_error="
              << metrics.max_message_error
              << " message_bits=" << metrics.message_bits
              << " max_basis_error=" << metrics.max_basis_error
              << " basis_bits=" << metrics.basis_bits
              << " output_margin_bits=" << output_margin_bits
              << " sample_messages=257\n";
}

template <class Schedule>
void print_schedule_report(const char *label,
                           const std::filesystem::path *key_dir = nullptr)
{
    using P = typename Schedule::Param;
    TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    TFHEpp::CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    TFHEpp::CKKSDenseBootstrapRotationKeyUsage<Schedule> usage;
    TFHEpp::CKKSBuildDenseBootstrapRotationKeyUsage<Schedule>(usage,
                                                              linear_plan);
    TFHEpp::CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> direct_usage;
    TFHEpp::CKKSBuildDenseBootstrapDirectRotationKeyUsage<Schedule>(
        direct_usage, linear_plan);
    TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        hybrid_usage;
    TFHEpp::CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        hybrid_usage, linear_plan);

    const std::size_t sparse_rows =
        TFHEpp::CKKSDenseBootstrapSparseKeySwitchRowCount<Schedule>(usage);
    const std::size_t sparse_bytes =
        TFHEpp::CKKSDenseBootstrapSparseKeyByteEstimate<Schedule>(usage);
    const std::size_t streamed_peak_rows =
        TFHEpp::CKKSDenseBootstrapStreamedKeySwitchPeakRowCount<Schedule>(
            usage);
    const std::size_t streamed_peak_bytes =
        TFHEpp::CKKSDenseBootstrapStreamedPeakKeyByteEstimate<Schedule>(
            usage);
    const std::size_t direct_rows =
        TFHEpp::CKKSDenseBootstrapDirectKeySwitchRowCount<Schedule>(
            direct_usage);
    const std::size_t direct_bytes =
        TFHEpp::CKKSDenseBootstrapDirectKeyByteEstimate<Schedule>(
            direct_usage);
    const std::size_t direct_streamed_peak_rows =
        TFHEpp::CKKSDenseBootstrapDirectStreamedKeySwitchPeakRowCount<
            Schedule>(direct_usage);
    const std::size_t direct_streamed_peak_bytes =
        TFHEpp::CKKSDenseBootstrapDirectStreamedPeakKeyByteEstimate<Schedule>(
            direct_usage);
    const std::size_t hybrid_rows =
        TFHEpp::CKKSDenseBootstrapHybridGiantKeySwitchRowCount<Schedule>(
            hybrid_usage);
    const std::size_t hybrid_bytes =
        TFHEpp::CKKSDenseBootstrapHybridGiantKeyByteEstimate<Schedule>(
            hybrid_usage);
    const std::size_t hybrid_streamed_peak_rows =
        TFHEpp::CKKSDenseBootstrapHybridGiantStreamedKeySwitchPeakRowCount<
            Schedule>(hybrid_usage);
    const std::size_t hybrid_streamed_peak_bytes =
        TFHEpp::
            CKKSDenseBootstrapHybridGiantStreamedPeakKeyByteEstimate<Schedule>(
                hybrid_usage);
    const std::size_t c2s_current_evalautos =
        TFHEpp::CKKSLinearTransformStagesRotationEvalAutoCount<P>(
            linear_plan.coeff_to_slot_stages, 0,
            linear_plan.coeff_to_slot_stages.size(),
            Schedule::linear_bsgs_step);
    const std::size_t c2s_direct_evalautos =
        TFHEpp::CKKSLinearTransformStagesDirectRotationEvalAutoCount<P>(
            linear_plan.coeff_to_slot_stages, 0,
            linear_plan.coeff_to_slot_stages.size(),
            Schedule::linear_bsgs_step);
    const std::size_t c2s_hybrid_evalautos =
        TFHEpp::CKKSLinearTransformStagesHybridGiantRotationEvalAutoCount<P>(
            linear_plan.coeff_to_slot_stages, 0,
            linear_plan.coeff_to_slot_stages.size(),
            Schedule::linear_bsgs_step,
            Schedule::hybrid_giant_direct_popcount_threshold);
    const std::size_t stc_current_evalautos =
        TFHEpp::CKKSLinearTransformStagesDualInputSharedTailRotationEvalAutoCount<
            P>(linear_plan.slot_to_coeff_stages, 0,
               linear_plan.slot_to_coeff_stages.size(),
               Schedule::linear_bsgs_step);
    const std::size_t stc_direct_evalautos =
        TFHEpp::
            CKKSLinearTransformStagesDualInputSharedTailDirectRotationEvalAutoCount<
                P>(linear_plan.slot_to_coeff_stages, 0,
                   linear_plan.slot_to_coeff_stages.size(),
                   Schedule::linear_bsgs_step);
    const std::size_t stc_hybrid_evalautos =
        TFHEpp::
            CKKSLinearTransformStagesDualInputSharedTailHybridGiantRotationEvalAutoCount<
                P>(linear_plan.slot_to_coeff_stages, 0,
                   linear_plan.slot_to_coeff_stages.size(),
                   Schedule::linear_bsgs_step,
                   Schedule::hybrid_giant_direct_popcount_threshold);

    std::cout << label << " n=" << P::n << " logQ="
              << Schedule::boot_log_q << " input_logQ="
              << Schedule::input_log_q << " output_logQ="
              << Schedule::output_log_q << '\n';
    std::cout << label << " c2s_levels="
              << Schedule::coeff_to_slot_level_count << " stc_levels="
              << Schedule::slot_to_coeff_level_count << " evalmod_depth="
              << Schedule::evalmod_depth << '\n';
    std::cout << label << " evalmod degree/k/double_angle/inv/log_scale="
              << Schedule::evalmod_degree << "/" << Schedule::evalmod_k << "/"
              << Schedule::evalmod_double_angle << "/"
              << Schedule::evalmod_inv_degree << "/"
              << Schedule::evalmod_log_scale
              << " log_q_loss=" << Schedule::evalmod_log_q_consumption
              << " after_evalmod_logQ=" << Schedule::after_evalmod_log_q
              << '\n';
    std::cout << label << " modraise_mask_bound="
              << Schedule::modraise_mask_bound
              << " evalmod_mask_bound=" << Schedule::evalmod_mask_bound
              << " bounded_sparse_weight_max="
              << bounded_modraise_sparse_weight_max<Schedule>() << '\n';
    std::cout << label << " linear_plain_log_delta c2s/split/stc="
              << Schedule::coeff_to_slot_plain_log_delta << "/"
              << Schedule::component_split_plain_log_delta << "/"
              << Schedule::slot_to_coeff_plain_log_delta
              << " fuse c2s/stc=" << Schedule::coeff_to_slot_fuse_radix << "/"
              << Schedule::slot_to_coeff_fuse_radix << '\n';
    std::cout << label << " rotation_indices="
              << TFHEpp::CKKSDenseBootstrapRotationKeyUsageCount<Schedule>(
                     usage)
              << "/"
              << TFHEpp::CKKSDenseBootstrapFullGaloisKeyIndexCount<Schedule>()
              << " sparse_rows=" << sparse_rows
              << " streamed_peak_rows=" << streamed_peak_rows << '\n';
    std::cout << label << " direct_rotation_indices="
              << TFHEpp::CKKSDenseBootstrapDirectRotationKeyUsageCount<
                     Schedule>(direct_usage)
              << " direct_rows=" << direct_rows
              << " direct_streamed_peak_rows=" << direct_streamed_peak_rows
              << '\n';
    std::cout << label << " hybrid_rotation_indices="
              << TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsageCount<
                     Schedule>(hybrid_usage)
              << " hybrid_rows=" << hybrid_rows
              << " hybrid_streamed_peak_rows=" << hybrid_streamed_peak_rows
              << '\n';
    std::cout << label << " rotation_evalautos current/direct c2s="
              << c2s_current_evalautos << "/" << c2s_direct_evalautos
              << " stc=" << stc_current_evalautos << "/"
              << stc_direct_evalautos << '\n';
    std::cout << label << " rotation_evalautos hybrid c2s="
              << c2s_hybrid_evalautos << " stc=" << stc_hybrid_evalautos
              << " direct_popcount_threshold="
              << Schedule::hybrid_giant_direct_popcount_threshold
              << '\n';
    std::cout << label << " sparse_key_bytes=" << sparse_bytes
              << " streamed_peak_bytes=" << streamed_peak_bytes
              << " direct_key_bytes=" << direct_bytes
              << " direct_streamed_peak_bytes=" << direct_streamed_peak_bytes
              << " hybrid_key_bytes=" << hybrid_bytes
              << " hybrid_streamed_peak_bytes="
              << hybrid_streamed_peak_bytes
              << " full_key_bytes="
              << TFHEpp::CKKSDenseBootstrapFullKeyByteEstimate<Schedule>()
              << '\n';
    print_evalmod_approximation_report<Schedule>(label);
    if (key_dir != nullptr) {
        const auto expected =
            TFHEpp::CKKSDenseBootstrapKeyDirectoryFiles<Schedule>(*key_dir);
        const auto missing =
            TFHEpp::CKKSDenseBootstrapMissingKeyDirectoryFiles<Schedule>(
                *key_dir);
        std::cout << label << " key_dir=" << key_dir->string()
                  << " files=" << regular_file_count(*key_dir) << "/"
                  << expected.size() << " missing=" << missing.size()
                  << " manifest=" << manifest_status<Schedule>(*key_dir)
                  << " tmp_files=" << temporary_file_count(*key_dir)
                  << " disk_bytes=" << directory_size_bytes(*key_dir)
                  << '\n';
    }
}

void print_lvl6_hybrid_threshold_reports()
{
    print_schedule_report<Lvl6HybridThresholdSchedule<1>>("lvl6-hybrid-th1");
    print_schedule_report<Lvl6HybridThresholdSchedule<2>>("lvl6-hybrid-th2");
    print_schedule_report<Lvl6HybridThresholdSchedule<3>>("lvl6-hybrid-th3");
    print_schedule_report<Lvl6HybridThresholdSchedule<4>>("lvl6-hybrid-th4");
    print_schedule_report<Lvl6HybridThresholdSchedule<5>>("lvl6-hybrid-th5");
}

void print_lvl6_robust_reports()
{
    print_schedule_report<Lvl6RobustHybridThresholdSchedule<3>>(
        "lvl6-robust-th3");
    print_schedule_report<Lvl6RobustHybridThresholdSchedule<4>>(
        "lvl6-robust-th4");
}

void print_lvl6_inverse_reports()
{
    Lvl6InverseBudgetParams current_with_inverse;
    print_lvl6_inverse_budget("lvl6-inverse-current-log50",
                              current_with_inverse);

    Lvl6InverseBudgetParams selected;
    selected.log_delta = 40;
    selected.coeff_to_slot_plain_log_delta = 50;
    selected.component_split_plain_log_delta = 50;
    selected.slot_to_coeff_plain_log_delta = 20;
    selected.evalmod_degree = 52;
    selected.evalmod_double_angle = 4;
    selected.evalmod_inv_degree = 3;
    selected.evalmod_log_scale = 40;
    print_lvl6_inverse_budget("lvl6-inverse-selected-budget", selected);
    print_schedule_report<Lvl6InverseSchedule>("lvl6-inverse-selected");
}

struct Lvl6InverseSearchResult {
    Lvl6InverseBudgetParams params{};
    Lvl6InverseBudgetResult budget{};
    EvalModApproximationMetrics approx{};
};

void print_lvl6_inverse_search_report()
{
    std::vector<Lvl6InverseSearchResult> results;

    for (const std::uint32_t log_delta : {36U, 38U, 40U, 42U, 44U}) {
        for (const std::uint32_t degree : {34U, 40U, 46U, 52U, 58U, 63U}) {
            for (const std::uint32_t double_angle : {2U, 3U, 4U}) {
                for (const std::uint32_t inverse_degree : {3U, 5U, 7U, 9U}) {
                    Lvl6InverseBudgetParams params;
                    params.log_delta = log_delta;
                    params.coeff_to_slot_plain_log_delta = log_delta;
                    params.component_split_plain_log_delta = log_delta;
                    params.slot_to_coeff_plain_log_delta = 25;
                    params.evalmod_degree = degree;
                    params.evalmod_double_angle = double_angle;
                    params.evalmod_inv_degree = inverse_degree;
                    params.evalmod_log_scale = log_delta;

                    const Lvl6InverseBudgetResult budget =
                        compute_lvl6_inverse_budget(params);
                    if (!budget.fits || budget.output_margin_bits < 30)
                        continue;

                    const EvalModApproximationMetrics approx =
                        measure_evalmod_approximation(
                            params.evalmod_k, params.evalmod_degree,
                            params.log_message_ratio,
                            params.evalmod_double_angle,
                            params.evalmod_inv_degree);
                    if (approx.message_bits < 18.0) continue;

                    results.push_back({params, budget, approx});
                }
            }
        }
    }

    std::sort(results.begin(), results.end(),
              [](const Lvl6InverseSearchResult &lhs,
                 const Lvl6InverseSearchResult &rhs) {
                  if (lhs.approx.message_bits != rhs.approx.message_bits)
                      return lhs.approx.message_bits > rhs.approx.message_bits;
                  if (lhs.budget.output_margin_bits !=
                      rhs.budget.output_margin_bits)
                      return lhs.budget.output_margin_bits >
                             rhs.budget.output_margin_bits;
                  return lhs.params.log_delta > rhs.params.log_delta;
              });

    std::cout << "lvl6-inverse-search candidates=" << results.size()
              << " shown=" << std::min<std::size_t>(results.size(), 12)
              << '\n';
    for (std::size_t i = 0; i < std::min<std::size_t>(results.size(), 12);
         i++) {
        const Lvl6InverseSearchResult &r = results[i];
        std::cout << "lvl6-inverse-search rank=" << (i + 1)
                  << " log_delta=" << r.params.log_delta
                  << " degree=" << r.params.evalmod_degree
                  << " double_angle=" << r.params.evalmod_double_angle
                  << " inv_degree=" << r.params.evalmod_inv_degree
                  << " output_logQ=" << r.budget.output_log_q
                  << " output_margin_bits=" << r.budget.output_margin_bits
                  << " evalmod_log_q_loss=" << r.budget.evalmod_log_q_loss
                  << " message_bits=" << r.approx.message_bits
                  << " max_message_error=" << r.approx.max_message_error
                  << " basis_bits=" << r.approx.basis_bits << '\n';
    }
}

template <class Schedule>
int validate_filesystem_key_dir(const std::filesystem::path &key_dir)
{
    const auto missing =
        TFHEpp::CKKSDenseBootstrapMissingKeyDirectoryFiles<Schedule>(key_dir);
    if (!missing.empty()) {
        std::cerr << "key_dir_incomplete=" << key_dir.string()
                  << " missing=" << missing.size() << '\n';
        print_missing_key_files(missing);
        return 2;
    }
    try {
        if (!TFHEpp::CKKSDenseBootstrapKeyDirectoryManifestMatches<Schedule>(
                key_dir)) {
            std::cerr << "key_dir_manifest_mismatch=" << key_dir.string()
                      << '\n';
            return 2;
        }
    }
    catch (const std::exception &e) {
        std::cerr << "key_dir_manifest_unreadable=" << key_dir.string()
                  << " error=" << e.what() << '\n';
        return 2;
    }
    return 0;
}

template <class Schedule>
int run_filesystem_bootstrap(const std::filesystem::path &key_dir, double tol,
                             std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "bootstrap");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, false);
        status != 0)
        return status;
    if (const int status = validate_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);
    std::cout << "key_sparse_weight=" << sparse_weight << '\n';

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    fill_test_slots<P>(*slots);

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    const double encrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotEncrypt<P, Schedule::input_log_q,
                                Schedule::log_delta>(*input, *slots, *key);
    });

    TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapTimings bootstrap_timings;
    const double bootstrap_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseBootstrapWithKeyProviderTimed<Schedule>(
            *output, *input, provider, bootstrap_timings);
    });

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    const double decrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q,
                                Schedule::log_delta>(*decoded, *output, *key);
    });
    const double err = max_error<P>(*decoded, *slots);
    std::cout << "encrypt_ms=" << encrypt_ms << '\n';
    std::cout << "bootstrap_ms=" << bootstrap_ms << '\n';
    print_bootstrap_timings(bootstrap_timings);
    std::cout << "decrypt_ms=" << decrypt_ms << '\n';
    std::cout << "max_error=" << err << '\n';
    return err <= tol ? 0 : 1;
}

template <class Schedule>
int validate_hybrid_filesystem_key_dir(const std::filesystem::path &key_dir)
{
    const auto missing =
        TFHEpp::CKKSDenseBootstrapHybridGiantMissingKeyDirectoryFiles<
            Schedule>(key_dir);
    if (!missing.empty()) {
        std::cerr << "hybrid_key_dir_incomplete=" << key_dir.string()
                  << " missing=" << missing.size() << '\n';
        print_missing_key_files(missing);
        return 2;
    }
    try {
        if (!TFHEpp::CKKSDenseBootstrapHybridGiantKeyDirectoryManifestMatches<
                Schedule>(key_dir)) {
            std::cerr << "hybrid_key_dir_manifest_mismatch="
                      << key_dir.string() << '\n';
            return 2;
        }
    }
    catch (const std::exception &e) {
        std::cerr << "hybrid_key_dir_manifest_unreadable="
                  << key_dir.string() << " error=" << e.what() << '\n';
        return 2;
    }
    return 0;
}

template <class Schedule>
int run_hybrid_filesystem_bootstrap(const std::filesystem::path &key_dir,
                                    double tol,
                                    std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "bootstrap");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, false);
        status != 0)
        return status;
    if (const int status =
            validate_hybrid_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);
    std::cout << "key_sparse_weight=" << sparse_weight << '\n';

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    fill_test_slots<P>(*slots);

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    const double encrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotEncrypt<P, Schedule::input_log_q,
                                Schedule::log_delta>(*input, *slots, *key);
    });

    TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapTimings bootstrap_timings;
    const double bootstrap_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseBootstrapWithKeyProviderTimed<Schedule>(
            *output, *input, provider, bootstrap_timings);
    });

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    const double decrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q,
                                Schedule::log_delta>(*decoded, *output, *key);
    });
    const double err = max_error<P>(*decoded, *slots);
    std::cout << "encrypt_ms=" << encrypt_ms << '\n';
    std::cout << "bootstrap_ms=" << bootstrap_ms << '\n';
    print_bootstrap_timings(bootstrap_timings);
    std::cout << "decrypt_ms=" << decrypt_ms << '\n';
    std::cout << "max_error=" << err << '\n';
    return err <= tol ? 0 : 1;
}

template <class Schedule>
int run_modraise_plain_diagnostics(std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);
    std::cout << "diag_key_sparse_weight=" << sparse_weight << '\n';
    std::cout << "diag_modraise_bounded_sparse_weight_max="
              << bounded_modraise_sparse_weight_max<Schedule>()
              << " diag_modraise_sparse_weight_valid="
              << (sparse_weight_fits_bounded_modraise<Schedule>(sparse_weight)
                      ? 1
                      : 0)
              << '\n';

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    fill_test_slots<P>(*slots);

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    TFHEpp::ckksSlotEncrypt<P, Schedule::input_log_q, Schedule::log_delta>(
        *input, *slots, *key);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::input_log_q, Schedule::log_delta>(
        *decoded, *input, *key);
    print_slot_diagnostic<P>("diag_input", *decoded, slots.get());
    print_phase_coefficient_diagnostic<Schedule, Schedule::input_log_q>(
        "diag_input", *input, *key);

    auto raised = std::make_unique<typename Schedule::BootstrapCiphertext>();
    TFHEpp::CKKSModRaiseBoundedPhaseRandomized<
        P, Schedule::input_log_q, Schedule::boot_log_q, Schedule::log_delta,
        Schedule::modraise_mask_bound>(*raised, *input);
    auto raised_slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::boot_log_q, Schedule::log_delta>(
        *raised_slots, *raised, *key);
    print_slot_diagnostic<P>("diag_raised", *raised_slots);
    print_phase_coefficient_diagnostic<Schedule, Schedule::boot_log_q>(
        "diag_raised", *raised, *key);

    TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    TFHEpp::CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    auto expected_c2s = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    apply_complex_stages<P>(*expected_c2s, *raised_slots,
                            linear_plan.coeff_to_slot_stages);
    print_slot_diagnostic<P>("diag_plain_c2s", *expected_c2s);
    return 0;
}

template <class Schedule, class KeyProvider>
int run_key_provider_bootstrap_diagnostics(const std::filesystem::path &key_dir,
                                           bool full_pipeline,
                                           std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "diagnostic");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, false);
        status != 0)
        return status;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);
    std::cout << "diag_key_sparse_weight=" << sparse_weight << '\n';

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    fill_test_slots<P>(*slots);

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    const double encrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotEncrypt<P, Schedule::input_log_q,
                                Schedule::log_delta>(*input, *slots, *key);
    });
    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::input_log_q, Schedule::log_delta>(
        *decoded, *input, *key);
    print_slot_diagnostic<P>("diag_input", *decoded, slots.get());
    print_phase_coefficient_diagnostic<Schedule, Schedule::input_log_q>(
        "diag_input", *input, *key);

    KeyProvider provider(key_dir);
    const TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan =
        provider.linear_plan();

    auto raised = std::make_unique<typename Schedule::BootstrapCiphertext>();
    const double modraise_ms = elapsed_ms([&] {
        TFHEpp::CKKSModRaiseBoundedPhaseRandomized<
            P, Schedule::input_log_q, Schedule::boot_log_q,
            Schedule::log_delta, Schedule::modraise_mask_bound>(*raised, *input);
    });
    auto raised_slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::boot_log_q, Schedule::log_delta>(
        *raised_slots, *raised, *key);
    print_slot_diagnostic<P>("diag_raised", *raised_slots);
    print_phase_coefficient_diagnostic<Schedule, Schedule::boot_log_q>(
        "diag_raised", *raised, *key);

    auto expected_c2s = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    apply_complex_stages<P>(*expected_c2s, *raised_slots,
                            linear_plan.coeff_to_slot_stages);

    auto coeff_to_slot =
        std::make_unique<typename Schedule::CoeffToSlotCiphertext>();
    const double c2s_ms = elapsed_ms([&] {
        const TFHEpp::ckks_detail::CKKSDenseBootstrapLinearKeyProviderChain<
            KeyProvider, true>
            coeff_to_slot_galois{provider};
        TFHEpp::CKKSLinearTransformStagesBSGS<
            P, Schedule::boot_log_q, Schedule::log_delta,
            Schedule::coeff_to_slot_plain_log_delta,
            Schedule::coeff_to_slot_level_count>(
            *coeff_to_slot, *raised, linear_plan.coeff_to_slot_stages, 0,
            Schedule::linear_bsgs_step, coeff_to_slot_galois);
    });
    raised.reset();
    raised_slots.reset();

    TFHEpp::ckksSlotDecrypt<P, Schedule::after_coeff_to_slot_log_q,
                            Schedule::log_delta>(*decoded, *coeff_to_slot, *key);
    print_slot_diagnostic<P>("diag_c2s", *decoded, expected_c2s.get());

    std::cout << "diag_encrypt_ms=" << encrypt_ms << '\n';
    std::cout << "diag_modraise_ms=" << modraise_ms << '\n';
    std::cout << "diag_c2s_ms=" << c2s_ms << '\n';
    if (!full_pipeline)
        return max_error<P>(*decoded, *expected_c2s) <= 0.1 ? 0 : 1;

    auto expected_real = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected_imag = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    for (std::size_t i = 0; i < P::n / 2; i++) {
        (*expected_real)[i] = {(*expected_c2s)[i].real(), 0.0};
        (*expected_imag)[i] = {(*expected_c2s)[i].imag(), 0.0};
    }

    auto real_component =
        std::make_unique<typename Schedule::ComponentCiphertext>();
    auto imag_component =
        std::make_unique<typename Schedule::ComponentCiphertext>();
    const double split_ms = elapsed_ms([&] {
        TFHEpp::CKKSExtractRealSlots<P, Schedule::after_coeff_to_slot_log_q,
                                     Schedule::log_delta,
                                     Schedule::component_split_plain_log_delta>(
            *real_component, *coeff_to_slot,
            provider.packed_conjugate_galois());
        TFHEpp::CKKSExtractImagSlots<P, Schedule::after_coeff_to_slot_log_q,
                                     Schedule::log_delta,
                                     Schedule::component_split_plain_log_delta>(
            *imag_component, *coeff_to_slot,
            provider.packed_conjugate_galois());
    });
    coeff_to_slot.reset();
    if constexpr (requires { provider.release_packed_conjugate_galois(); }) {
        provider.release_packed_conjugate_galois();
    }

    TFHEpp::ckksSlotDecrypt<P, Schedule::after_component_split_log_q,
                            Schedule::log_delta>(*decoded, *real_component,
                                                 *key);
    print_slot_diagnostic<P>("diag_real_split", *decoded, expected_real.get());
    TFHEpp::ckksSlotDecrypt<P, Schedule::after_component_split_log_q,
                            Schedule::log_delta>(*decoded, *imag_component,
                                                 *key);
    print_slot_diagnostic<P>("diag_imag_split", *decoded, expected_imag.get());
    expected_c2s.reset();

    const TFHEpp::CKKSBoundedCosEvalModPolynomial &poly =
        provider.evalmod_polynomial();
    auto expected_real_eval = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected_imag_eval = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    for (std::size_t i = 0; i < P::n / 2; i++) {
        (*expected_real_eval)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, (*expected_real)[i].real(),
                 Schedule::evalmod_inv_degree) /
                 Schedule::message_ratio,
             0.0};
        (*expected_imag_eval)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, (*expected_imag)[i].real(),
                 Schedule::evalmod_inv_degree) /
                 Schedule::message_ratio,
             0.0};
    }
    expected_real.reset();
    expected_imag.reset();

    auto real_evalmod =
        std::make_unique<TFHEpp::CKKSDenseEvalModBoundedCosResult<Schedule>>();
    const double real_evalmod_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseEvalModBoundedCosNormalizedWithKeyProvider<Schedule>(
            *real_evalmod, *real_component, poly, provider);
    });
    real_component.reset();
    auto actual_real_eval = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::after_evalmod_log_q,
                            Schedule::log_delta>(*decoded, *real_evalmod, *key);
    print_slot_diagnostic<P>("diag_real_evalmod", *decoded,
                             expected_real_eval.get());
    *actual_real_eval = *decoded;

    auto imag_evalmod =
        std::make_unique<TFHEpp::CKKSDenseEvalModBoundedCosResult<Schedule>>();
    const double imag_evalmod_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseEvalModBoundedCosNormalizedWithKeyProvider<Schedule>(
            *imag_evalmod, *imag_component, poly, provider);
    });
    imag_component.reset();
    auto actual_imag_eval = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::after_evalmod_log_q,
                            Schedule::log_delta>(*decoded, *imag_evalmod, *key);
    print_slot_diagnostic<P>("diag_imag_evalmod", *decoded,
                             expected_imag_eval.get());
    *actual_imag_eval = *decoded;

    auto expected_real_out = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected_imag_out = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected_pipeline = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto actual_real_out = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto actual_imag_out = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto actual_pipeline = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    apply_complex_stages<P>(*expected_real_out, *expected_real_eval,
                            linear_plan.slot_to_coeff_stages);
    apply_complex_stages<P>(*expected_imag_out, *expected_imag_eval,
                            linear_plan.slot_to_coeff_imag_stages);
    apply_complex_stages<P>(*actual_real_out, *actual_real_eval,
                            linear_plan.slot_to_coeff_stages);
    apply_complex_stages<P>(*actual_imag_out, *actual_imag_eval,
                            linear_plan.slot_to_coeff_imag_stages);
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*expected_pipeline)[i] =
            (*expected_real_out)[i] + (*expected_imag_out)[i];
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*actual_pipeline)[i] = (*actual_real_out)[i] + (*actual_imag_out)[i];
    print_slot_diagnostic<P>("diag_actual_pipeline_vs_pipeline",
                             *actual_pipeline, expected_pipeline.get());
    print_slot_diagnostic<P>("diag_actual_pipeline_vs_message",
                             *actual_pipeline, slots.get());
    expected_real_eval.reset();
    expected_imag_eval.reset();
    expected_real_out.reset();
    expected_imag_out.reset();
    actual_real_eval.reset();
    actual_imag_eval.reset();
    actual_real_out.reset();
    actual_imag_out.reset();

    const TFHEpp::ckks_detail::CKKSDenseBootstrapLinearKeyProviderChain<
        KeyProvider, false>
        slot_to_coeff_galois{provider};
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    const double stc_ms = elapsed_ms([&] {
        TFHEpp::CKKSLinearTransformStagesBSGSDualInputSharedTail<
            P, Schedule::after_evalmod_log_q, Schedule::log_delta,
            Schedule::slot_to_coeff_plain_log_delta,
            Schedule::slot_to_coeff_level_count>(
            *output, *real_evalmod, *imag_evalmod,
            linear_plan.slot_to_coeff_stages,
            linear_plan.slot_to_coeff_imag_stages, 0,
            Schedule::linear_bsgs_step, slot_to_coeff_galois);
    });
    real_evalmod.reset();
    imag_evalmod.reset();

    TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *output, *key);
    print_slot_diagnostic<P>("diag_output_vs_actual_pipeline", *decoded,
                             actual_pipeline.get());
    print_slot_diagnostic<P>("diag_output_vs_pipeline", *decoded,
                             expected_pipeline.get());
    print_slot_diagnostic<P>("diag_output_vs_message", *decoded, slots.get());
    print_slot_diagnostic<P>("diag_pipeline_vs_message", *expected_pipeline,
                             slots.get());

    std::cout << "diag_split_ms=" << split_ms << '\n';
    std::cout << "diag_real_evalmod_ms=" << real_evalmod_ms << '\n';
    std::cout << "diag_imag_evalmod_ms=" << imag_evalmod_ms << '\n';
    std::cout << "diag_stc_ms=" << stc_ms << '\n';
    return max_error<P>(*decoded, *slots) <= 0.1 ? 0 : 1;
}

template <class Schedule>
int run_filesystem_bootstrap_diagnostics(const std::filesystem::path &key_dir,
                                         bool full_pipeline,
                                         std::size_t sparse_weight = 0)
{
    if (const int status = validate_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;
    return run_key_provider_bootstrap_diagnostics<
        Schedule, TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule>>(
        key_dir, full_pipeline, sparse_weight);
}

template <class Schedule>
int run_hybrid_filesystem_bootstrap_diagnostics(
    const std::filesystem::path &key_dir, bool full_pipeline,
    std::size_t sparse_weight = 0)
{
    if (const int status =
            validate_hybrid_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;
    return run_key_provider_bootstrap_diagnostics<
        Schedule,
        TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule>>(
        key_dir, full_pipeline, sparse_weight);
}

template <class Schedule, class KeyProvider>
int run_key_provider_evalmod_diagnostics(const std::filesystem::path &key_dir,
                                         std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "diagnostic");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, false);
        status != 0)
        return status;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);
    std::cout << "diag_key_sparse_weight=" << sparse_weight << '\n';

    KeyProvider provider(key_dir);
    const TFHEpp::CKKSBoundedCosEvalModPolynomial &poly =
        provider.evalmod_polynomial();

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    constexpr int k = static_cast<int>(Schedule::evalmod_k);
    for (std::size_t i = 0; i < P::n / 2; i++) {
        const int mask = static_cast<int>(i % (2 * k + 3)) - k - 1;
        const double message =
            static_cast<double>(static_cast<int>(i % 17) - 8) / 256.0;
        const double normalized =
            (static_cast<double>(mask) * Schedule::message_ratio + message) /
            (static_cast<double>(Schedule::evalmod_k) *
             Schedule::message_ratio);
        (*slots)[i] = {normalized, 0.0};
        (*expected)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, normalized, Schedule::evalmod_inv_degree) /
                 Schedule::message_ratio,
             0.0};
    }
    print_slot_diagnostic<P>("diag_evalmod_direct_input", *slots);
    print_slot_diagnostic<P>("diag_evalmod_direct_expected", *expected);

    auto input = std::make_unique<typename Schedule::ComponentCiphertext>();
    const double encrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotEncrypt<P, Schedule::after_component_split_log_q,
                                Schedule::log_delta>(*input, *slots, *key);
    });

    auto output =
        std::make_unique<TFHEpp::CKKSDenseEvalModBoundedCosResult<Schedule>>();
    const double evalmod_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseEvalModBoundedCosNormalizedWithKeyProvider<Schedule>(
            *output, *input, poly, provider);
    });
    input.reset();

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::after_evalmod_log_q,
                            Schedule::log_delta>(*decoded, *output, *key);
    print_slot_diagnostic<P>("diag_evalmod_direct", *decoded, expected.get());
    std::cout << "diag_encrypt_ms=" << encrypt_ms << '\n';
    std::cout << "diag_evalmod_ms=" << evalmod_ms << '\n';
    return max_error<P>(*decoded, *expected) <= 0.01 ? 0 : 1;
}

template <class Schedule>
int run_filesystem_evalmod_diagnostics(const std::filesystem::path &key_dir,
                                       std::size_t sparse_weight = 0)
{
    if (const int status = validate_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;
    return run_key_provider_evalmod_diagnostics<
        Schedule, TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule>>(
        key_dir, sparse_weight);
}

template <class Schedule>
int run_hybrid_filesystem_evalmod_diagnostics(
    const std::filesystem::path &key_dir, std::size_t sparse_weight = 0)
{
    if (const int status =
            validate_hybrid_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;
    return run_key_provider_evalmod_diagnostics<
        Schedule,
        TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule>>(
        key_dir, sparse_weight);
}

template <class Schedule, class KeyProvider>
int run_key_provider_stc_diagnostics(const std::filesystem::path &key_dir,
                                     std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "diagnostic");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, false);
        status != 0)
        return status;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);
    std::cout << "diag_key_sparse_weight=" << sparse_weight << '\n';

    KeyProvider provider(key_dir);
    const TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan =
        provider.linear_plan();

    const TFHEpp::CKKSBoundedCosEvalModPolynomial &poly =
        provider.evalmod_polynomial();

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    fill_test_slots<P>(*slots);
    print_slot_diagnostic<P>("diag_stc_message", *slots);

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    const double input_encrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotEncrypt<P, Schedule::input_log_q, Schedule::log_delta>(
            *input, *slots, *key);
    });

    auto raised = std::make_unique<typename Schedule::BootstrapCiphertext>();
    const double modraise_ms = elapsed_ms([&] {
        TFHEpp::CKKSModRaiseBoundedPhaseRandomized<
            P, Schedule::input_log_q, Schedule::boot_log_q, Schedule::log_delta,
            Schedule::modraise_mask_bound>(*raised, *input);
    });
    input.reset();

    auto raised_slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::boot_log_q, Schedule::log_delta>(
        *raised_slots, *raised, *key);
    raised.reset();

    auto coeff_to_slot = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    apply_complex_stages<P>(*coeff_to_slot, *raised_slots,
                            linear_plan.coeff_to_slot_stages);
    raised_slots.reset();
    print_slot_diagnostic<P>("diag_stc_plain_c2s", *coeff_to_slot);

    auto real_eval = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto imag_eval = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    for (std::size_t i = 0; i < P::n / 2; i++) {
        (*real_eval)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, (*coeff_to_slot)[i].real(),
                 Schedule::evalmod_inv_degree) /
                 Schedule::message_ratio,
             0.0};
        (*imag_eval)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, (*coeff_to_slot)[i].imag(),
                 Schedule::evalmod_inv_degree) /
                 Schedule::message_ratio,
             0.0};
    }
    coeff_to_slot.reset();

    auto expected_real_out = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected_imag_out = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    apply_complex_stages<P>(*expected_real_out, *real_eval,
                            linear_plan.slot_to_coeff_stages);
    apply_complex_stages<P>(*expected_imag_out, *imag_eval,
                            linear_plan.slot_to_coeff_imag_stages);
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*expected)[i] = (*expected_real_out)[i] + (*expected_imag_out)[i];
    print_slot_diagnostic<P>("diag_stc_real_input", *real_eval);
    print_slot_diagnostic<P>("diag_stc_imag_input", *imag_eval);
    print_slot_diagnostic<P>("diag_stc_expected", *expected, slots.get());

    auto real_input = std::make_unique<typename Schedule::EvalModCiphertext>();
    auto imag_input = std::make_unique<typename Schedule::EvalModCiphertext>();
    const double eval_encrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotEncrypt<P, Schedule::after_evalmod_log_q,
                                Schedule::log_delta>(*real_input, *real_eval,
                                                     *key);
        TFHEpp::ckksSlotEncrypt<P, Schedule::after_evalmod_log_q,
                                Schedule::log_delta>(*imag_input, *imag_eval,
                                                     *key);
    });
    real_eval.reset();
    imag_eval.reset();

    const TFHEpp::ckks_detail::CKKSDenseBootstrapLinearKeyProviderChain<
        KeyProvider, false>
        slot_to_coeff_galois{provider};
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    const double stc_ms = elapsed_ms([&] {
        TFHEpp::CKKSLinearTransformStagesBSGSDualInputSharedTail<
            P, Schedule::after_evalmod_log_q, Schedule::log_delta,
            Schedule::slot_to_coeff_plain_log_delta,
            Schedule::slot_to_coeff_level_count>(
            *output, *real_input, *imag_input,
            linear_plan.slot_to_coeff_stages,
            linear_plan.slot_to_coeff_imag_stages, 0,
            Schedule::linear_bsgs_step, slot_to_coeff_galois);
    });
    real_input.reset();
    imag_input.reset();

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *output, *key);
    print_slot_diagnostic<P>("diag_stc_output", *decoded, expected.get());
    print_slot_diagnostic<P>("diag_stc_output_vs_message", *decoded,
                             slots.get());
    std::cout << "diag_input_encrypt_ms=" << input_encrypt_ms << '\n';
    std::cout << "diag_modraise_ms=" << modraise_ms << '\n';
    std::cout << "diag_eval_encrypt_ms=" << eval_encrypt_ms << '\n';
    std::cout << "diag_stc_ms=" << stc_ms << '\n';
    return max_error<P>(*decoded, *expected) <= 0.01 ? 0 : 1;
}

template <class Schedule>
int run_filesystem_stc_diagnostics(const std::filesystem::path &key_dir,
                                   std::size_t sparse_weight = 0)
{
    if (const int status = validate_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;
    return run_key_provider_stc_diagnostics<
        Schedule, TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule>>(
        key_dir, sparse_weight);
}

template <class Schedule>
int run_hybrid_filesystem_stc_diagnostics(
    const std::filesystem::path &key_dir, std::size_t sparse_weight = 0)
{
    if (const int status =
            validate_hybrid_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;
    return run_key_provider_stc_diagnostics<
        Schedule,
        TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule>>(
        key_dir, sparse_weight);
}

template <class Schedule>
int check_manifest_mismatch_rejected(const std::filesystem::path &key_dir)
{
    const std::filesystem::path mismatch_dir =
        key_dir.parent_path() / (key_dir.filename().string() + "_mismatch");
    std::filesystem::remove_all(mismatch_dir);
    std::filesystem::copy(key_dir, mismatch_dir,
                          std::filesystem::copy_options::recursive);

    auto manifest =
        TFHEpp::CKKSDenseBootstrapBuildKeyDirectoryManifest<Schedule>(
            mismatch_dir);
    manifest.boot_log_q++;
    TFHEpp::CKKSSavePortableBinaryAtomic(
        TFHEpp::CKKSDenseBootstrapKeyDirectoryManifestFile(mismatch_dir),
        manifest);

    bool rejected = false;
    try {
        TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(
            mismatch_dir);
    }
    catch (const std::exception &) {
        rejected = true;
    }
    std::filesystem::remove_all(mismatch_dir);
    if (!rejected) {
        std::cerr << "manifest_mismatch_was_accepted\n";
        return 1;
    }
    return 0;
}

template <class Schedule>
int run_toy_schedule_validation(bool keep_dir, const char *label,
                                const char *directory_name, double tol)
{
    using P = typename Schedule::Param;
    const std::filesystem::path key_dir =
        std::filesystem::temp_directory_path() / directory_name;
    std::filesystem::remove_all(key_dir);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_test_key<P>(*key);

    std::size_t generated_slices = 0;
    const double keygen_ms = elapsed_ms([&] {
        while (TFHEpp::CKKSDenseBootstrapKeyGenNextMissingToDirectory<Schedule>(
            key_dir, *key, {0.0, 0})) {
            generated_slices++;
        }
    });
    TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions resume_options;
    resume_options.overwrite_existing = false;
    TFHEpp::CKKSDenseBootstrapKeyGenToDirectory<Schedule>(
        key_dir, *key, {0.0, 0}, resume_options);

    print_schedule_report<Schedule>(label, &key_dir);
    std::cout << label << " keygen_ms=" << keygen_ms << '\n';
    std::cout << label << " keygen_slices=" << generated_slices << '\n';
    if (check_manifest_mismatch_rejected<Schedule>(key_dir) != 0) return 1;
    const int result = run_filesystem_bootstrap<Schedule>(key_dir, tol);

    if (!keep_dir) std::filesystem::remove_all(key_dir);
    return result;
}

template <class Schedule>
int run_toy_schedule_product_bootstrap_validation(
    bool keep_dir, const char *label, const char *directory_name, double tol)
{
    using P = typename Schedule::Param;
    constexpr std::uint32_t fresh_log_q =
        Schedule::input_log_q + 2 * Schedule::log_delta;
    using FreshCt = TFHEpp::CKKSCiphertext<P, fresh_log_q, Schedule::log_delta>;
    using ProductCt =
        TFHEpp::CKKSMultResult<P, FreshCt::log_q, FreshCt::log_delta,
                               FreshCt::log_q, FreshCt::log_delta>;
    static_assert(ProductCt::log_q > Schedule::input_log_q);

    const std::filesystem::path key_dir =
        std::filesystem::temp_directory_path() / directory_name;
    std::filesystem::remove_all(key_dir);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_test_key<P>(*key);
    TFHEpp::CKKSDenseBootstrapKeyGenToDirectory<Schedule>(key_dir, *key,
                                                          {0.0, 0});

    auto lhs = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto rhs = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    fill_test_slots<P>(*lhs);
    fill_alternate_test_slots<P>(*rhs);
    multiply_slots<P>(*expected, *lhs, *rhs);

    auto lhs_ct = std::make_unique<FreshCt>();
    auto rhs_ct = std::make_unique<FreshCt>();
    TFHEpp::ckksSlotEncrypt<P, FreshCt::log_q, FreshCt::log_delta>(
        *lhs_ct, *lhs, *key, {0.0, 0});
    TFHEpp::ckksSlotEncrypt<P, FreshCt::log_q, FreshCt::log_delta>(
        *rhs_ct, *rhs, *key, {0.0, 0});

    auto relin = TFHEpp::makeCKKSRelinKey<P, ProductCt::log_q>(*key, {0.0, 0});
    auto product = std::make_unique<ProductCt>();
    TFHEpp::CKKSMult<P>(*product, *lhs_ct, *rhs_ct, *relin);

    TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapTimings timings;
    TFHEpp::CKKSDenseBootstrapFromLevelWithKeyProviderTimed<Schedule>(
        *output, *product, provider, timings);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *output, *key);
    const double err = max_error<P>(*decoded, *expected);

    print_schedule_report<Schedule>(label, &key_dir);
    std::cout << label << " product_fresh_logQ=" << FreshCt::log_q
              << " product_logQ=" << ProductCt::log_q
              << " normalized_logQ=" << Schedule::input_log_q << '\n';
    print_bootstrap_timings(timings);
    std::cout << label << "_max_error=" << err << '\n';

    if (!keep_dir) std::filesystem::remove_all(key_dir);
    return err <= tol ? 0 : 1;
}

int run_toy_validation(bool keep_dir)
{
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        TinyDeepMultiLimbCKKSParam, 40, 8, 560, 40, 3, 31, 4, 2, 0, 40, 2>;
    return run_toy_schedule_validation<Schedule>(
        keep_dir, "toy", "tfhepp_ckks_bootstrap_validation_toy", 0.02);
}

int run_toy_inverse_validation(bool keep_dir)
{
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        TinyDeepMultiLimbCKKSParam, 30, 8, 560, 30, 3, 31, 4, 2, 5, 30, 2>;
    static_assert(Schedule::evalmod_inv_degree == 5);
    if (run_toy_schedule_validation<Schedule>(
            keep_dir, "toy-inverse",
            "tfhepp_ckks_bootstrap_validation_toy_inverse", 0.05) != 0)
        return 1;
    return run_toy_schedule_product_bootstrap_validation<Schedule>(
        keep_dir, "toy-inverse-product",
        "tfhepp_ckks_bootstrap_validation_toy_inverse_product", 0.05);
}

template <class Schedule>
int run_keygen_next(const std::filesystem::path &key_dir,
                    std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;
    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "keygen");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, true);
        status != 0)
        return status;
    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);

    print_schedule_report<Schedule>("keygen-next-before", &key_dir);
    std::cout << "keygen_next_sparse_weight=" << sparse_weight << '\n';
    const auto before_missing =
        TFHEpp::CKKSDenseBootstrapMissingKeyDirectoryFiles<Schedule>(key_dir);
    bool generated = false;
    const double keygen_ms = elapsed_ms([&] {
        generated =
            TFHEpp::CKKSDenseBootstrapKeyGenNextMissingToDirectory<Schedule>(
                key_dir, *key, {P::α, 0});
    });
    const auto after_missing =
        TFHEpp::CKKSDenseBootstrapMissingKeyDirectoryFiles<Schedule>(key_dir);
    print_schedule_report<Schedule>("keygen-next-after", &key_dir);
    std::cout << "keygen_next_generated=" << (generated ? 1 : 0) << '\n';
    std::cout << "keygen_next_ms=" << keygen_ms << '\n';
    print_created_key_files(before_missing, after_missing);
    return 0;
}

template <class Schedule>
int run_keygen(const std::filesystem::path &key_dir, bool resume,
               std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;
    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "keygen");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, true);
        status != 0)
        return status;
    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);

    TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions options;
    options.overwrite_existing = !resume;
    print_schedule_report<Schedule>("keygen-before", &key_dir);
    std::cout << "keygen_sparse_weight=" << sparse_weight << '\n';
    const double keygen_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseBootstrapKeyGenToDirectory<Schedule>(key_dir, *key,
                                                              {P::α, 0},
                                                              options);
    });
    print_schedule_report<Schedule>("keygen-after", &key_dir);
    std::cout << "keygen_ms=" << keygen_ms << '\n';
    return 0;
}

template <class Schedule>
int run_hybrid_keygen_next(const std::filesystem::path &key_dir,
                           std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;
    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "keygen");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, true);
        status != 0)
        return status;
    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);

    print_schedule_report<Schedule>("hybrid-keygen-next-before", &key_dir);
    std::cout << "keygen_next_sparse_weight=" << sparse_weight << '\n';
    const auto before_missing =
        TFHEpp::CKKSDenseBootstrapHybridGiantMissingKeyDirectoryFiles<
            Schedule>(key_dir);
    bool generated = false;
    const double keygen_ms = elapsed_ms([&] {
        generated =
            TFHEpp::CKKSDenseBootstrapHybridGiantKeyGenNextMissingToDirectory<
                Schedule>(key_dir, *key, {P::α, 0});
    });
    const auto after_missing =
        TFHEpp::CKKSDenseBootstrapHybridGiantMissingKeyDirectoryFiles<
            Schedule>(key_dir);
    print_schedule_report<Schedule>("hybrid-keygen-next-after", &key_dir);
    std::cout << "keygen_next_generated=" << (generated ? 1 : 0) << '\n';
    std::cout << "keygen_next_ms=" << keygen_ms << '\n';
    print_created_key_files(before_missing, after_missing);
    return 0;
}

template <class Schedule>
int run_hybrid_keygen(const std::filesystem::path &key_dir, bool resume,
                      std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;
    if (const int status =
            validate_bounded_modraise_test_key<Schedule>(sparse_weight,
                                                         "keygen");
        status != 0)
        return status;
    if (const int status = validate_or_create_validation_test_key_metadata<
            Schedule>(key_dir, sparse_weight, true);
        status != 0)
        return status;
    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);

    TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions options;
    options.overwrite_existing = !resume;
    print_schedule_report<Schedule>("hybrid-keygen-before", &key_dir);
    std::cout << "keygen_sparse_weight=" << sparse_weight << '\n';
    const double keygen_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseBootstrapHybridGiantKeyGenToDirectory<Schedule>(
            key_dir, *key, {P::α, 0}, options);
    });
    print_schedule_report<Schedule>("hybrid-keygen-after", &key_dir);
    std::cout << "keygen_ms=" << keygen_ms << '\n';
    return 0;
}

void print_usage(const char *program)
{
    std::cerr << "Usage: " << program
              << " [--toy] [--toy-inverse] [--keep]"
                 " [--lvl6-sparse-key H|--lvl6-dense-key]"
                 " [--lvl6-plan] [--lvl6-keygen DIR]"
                 " [--lvl6-hybrid-thresholds-plan]"
                 " [--lvl6-inverse-plan]"
                 " [--lvl6-inverse-search]"
                 " [--lvl6-inverse-debug-modraise]"
                 " [--lvl6-inverse-debug-modraise-sparse H]"
                 " [--lvl6-inverse-hybrid-keygen DIR]"
                 " [--lvl6-inverse-hybrid-keygen-next DIR]"
                 " [--lvl6-inverse-hybrid-run DIR]"
                 " [--lvl6-inverse-hybrid-debug-c2s DIR]"
                 " [--lvl6-inverse-hybrid-debug-evalmod DIR]"
                 " [--lvl6-inverse-hybrid-debug-stc DIR]"
                 " [--lvl6-inverse-hybrid-debug DIR]"
                 " [--lvl6-keygen-next DIR] [--lvl6-run DIR]"
                 " [--lvl6-hybrid-keygen DIR]"
                 " [--lvl6-hybrid-keygen-next DIR]"
                 " [--lvl6-hybrid-run DIR]"
                 " [--lvl6-robust-plan]"
                 " [--lvl6-fast-hybrid-keygen DIR]"
                 " [--lvl6-fast-hybrid-keygen-next DIR]"
                 " [--lvl6-fast-hybrid-run DIR]"
                 " [--lvl6-robust-hybrid-th3-keygen DIR]"
                 " [--lvl6-robust-hybrid-th3-keygen-next DIR]"
                 " [--lvl6-robust-hybrid-th3-run DIR]"
                 " [--lvl6-robust-hybrid-th3-debug-c2s DIR]"
                 " [--lvl6-robust-hybrid-th3-debug-evalmod DIR]"
                 " [--lvl6-robust-hybrid-th3-debug-stc DIR]"
                 " [--lvl6-robust-hybrid-th3-debug DIR]"
                 " [--lvl6-compact-hybrid-keygen DIR]"
                 " [--lvl6-compact-hybrid-keygen-next DIR]"
                 " [--lvl6-compact-hybrid-run DIR]"
                 " [--lvl6-compact-hybrid-debug-c2s DIR]"
                 " [--lvl6-compact-hybrid-debug-evalmod DIR]"
                 " [--lvl6-compact-hybrid-debug-stc DIR]"
                 " [--lvl6-compact-hybrid-debug DIR]"
                 " [--lvl6-robust-hybrid-th4-keygen DIR]"
                 " [--lvl6-robust-hybrid-th4-keygen-next DIR]"
                 " [--lvl6-robust-hybrid-th4-run DIR]"
                 " [--lvl6-hybrid-th3-keygen DIR]"
                 " [--lvl6-hybrid-th3-keygen-next DIR]"
                 " [--lvl6-hybrid-th3-run DIR]"
                 " [--lvl6-hybrid-th3-debug-c2s DIR]"
                 " [--lvl6-hybrid-th3-debug-evalmod DIR]"
                 " [--lvl6-hybrid-th3-debug-stc DIR]"
                 " [--lvl6-hybrid-th3-debug DIR]"
                 " [--lvl6-debug-modraise]"
                 " [--lvl6-debug-modraise-sparse H]"
                 " [--lvl6-debug-c2s DIR] [--lvl6-debug-evalmod DIR]"
                 " [--lvl6-debug-stc DIR] [--lvl6-debug DIR]"
                 " [--lvl6-hybrid-debug-c2s DIR]"
                 " [--lvl6-hybrid-debug-evalmod DIR]"
                 " [--lvl6-hybrid-debug-stc DIR]"
                 " [--lvl6-hybrid-debug DIR]"
                 " [--lvl6-all DIR] [--resume]\n";
}

}  // namespace

int main(int argc, char **argv)
{
    using Lvl6Schedule = TFHEpp::lvl6CKKSDenseBootstrapSchedule;
    using Lvl6HybridTh3Schedule = Lvl6HybridThresholdSchedule<3>;
    using Lvl6FastSchedule = TFHEpp::lvl6CKKSDenseBootstrapFastSchedule;
    using Lvl6CompactSchedule = TFHEpp::lvl6CKKSDenseBootstrapCompactSchedule;
    using Lvl6RobustHybridTh3Schedule = Lvl6FastSchedule;
    using Lvl6RobustHybridTh4Schedule = Lvl6CompactSchedule;
    constexpr std::size_t default_lvl6_sparse_weight =
        Lvl6Schedule::bounded_sparse_secret_key_weight;
    static_assert(default_lvl6_sparse_weight == 16);
    static_assert(Lvl6InverseSchedule::bounded_sparse_secret_key_weight ==
                  default_lvl6_sparse_weight);

    bool saw_action = false;
    bool keep = false;
    bool resume = false;
    std::size_t lvl6_sparse_weight = default_lvl6_sparse_weight;
    std::vector<std::string> args(argv + 1, argv + argc);

    for (std::size_t i = 0; i < args.size(); i++) {
        const std::string &arg = args[i];
        if (arg == "--keep") {
            keep = true;
        }
        else if (arg == "--resume") {
            resume = true;
        }
        else if (arg == "--lvl6-dense-key") {
            lvl6_sparse_weight = 0;
        }
        else if (arg == "--lvl6-sparse-key") {
            if (i + 1 >= args.size()) {
                print_usage(argv[0]);
                return 2;
            }
            lvl6_sparse_weight =
                static_cast<std::size_t>(std::stoull(args[++i]));
        }
        else if (arg == "--toy") {
            saw_action = true;
            if (run_toy_validation(keep) != 0) return 1;
        }
        else if (arg == "--toy-inverse") {
            saw_action = true;
            if (run_toy_inverse_validation(keep) != 0) return 1;
        }
        else if (arg == "--lvl6-plan") {
            saw_action = true;
            print_schedule_report<Lvl6Schedule>("lvl6");
        }
        else if (arg == "--lvl6-hybrid-thresholds-plan") {
            saw_action = true;
            print_lvl6_hybrid_threshold_reports();
        }
        else if (arg == "--lvl6-robust-plan") {
            saw_action = true;
            print_lvl6_robust_reports();
        }
        else if (arg == "--lvl6-inverse-plan") {
            saw_action = true;
            print_lvl6_inverse_reports();
        }
        else if (arg == "--lvl6-inverse-search") {
            saw_action = true;
            print_lvl6_inverse_search_report();
        }
        else if (arg == "--lvl6-inverse-debug-modraise") {
            saw_action = true;
            print_schedule_report<Lvl6InverseSchedule>("lvl6-inverse");
            if (run_modraise_plain_diagnostics<Lvl6InverseSchedule>(
                    lvl6_sparse_weight) != 0)
                return 1;
        }
        else if (arg == "--lvl6-inverse-debug-modraise-sparse") {
            if (i + 1 >= args.size()) {
                print_usage(argv[0]);
                return 2;
            }
            saw_action = true;
            const std::size_t sparse_weight =
                static_cast<std::size_t>(std::stoull(args[++i]));
            print_schedule_report<Lvl6InverseSchedule>("lvl6-inverse");
            if (run_modraise_plain_diagnostics<Lvl6InverseSchedule>(
                    sparse_weight) != 0)
                return 1;
        }
        else if (arg == "--lvl6-debug-modraise") {
            saw_action = true;
            print_schedule_report<Lvl6Schedule>("lvl6");
            if (run_modraise_plain_diagnostics<Lvl6Schedule>(
                    lvl6_sparse_weight) != 0)
                return 1;
        }
        else if (arg == "--lvl6-debug-modraise-sparse") {
            if (i + 1 >= args.size()) {
                print_usage(argv[0]);
                return 2;
            }
            saw_action = true;
            const std::size_t sparse_weight =
                static_cast<std::size_t>(std::stoull(args[++i]));
            print_schedule_report<Lvl6Schedule>("lvl6");
            if (run_modraise_plain_diagnostics<Lvl6Schedule>(sparse_weight) != 0)
                return 1;
        }
        else if (arg == "--lvl6-keygen" || arg == "--lvl6-keygen-next" ||
                 arg == "--lvl6-run" || arg == "--lvl6-debug-c2s" ||
                 arg == "--lvl6-hybrid-keygen" ||
                 arg == "--lvl6-hybrid-keygen-next" ||
                 arg == "--lvl6-hybrid-run" ||
                 arg == "--lvl6-inverse-hybrid-keygen" ||
                 arg == "--lvl6-inverse-hybrid-keygen-next" ||
                 arg == "--lvl6-inverse-hybrid-run" ||
                 arg == "--lvl6-inverse-hybrid-debug-c2s" ||
                 arg == "--lvl6-inverse-hybrid-debug-evalmod" ||
                 arg == "--lvl6-inverse-hybrid-debug-stc" ||
                 arg == "--lvl6-inverse-hybrid-debug" ||
                 arg == "--lvl6-fast-hybrid-keygen" ||
                 arg == "--lvl6-fast-hybrid-keygen-next" ||
                 arg == "--lvl6-fast-hybrid-run" ||
                 arg == "--lvl6-robust-hybrid-th3-keygen" ||
                 arg == "--lvl6-robust-hybrid-th3-keygen-next" ||
                 arg == "--lvl6-robust-hybrid-th3-run" ||
                 arg == "--lvl6-robust-hybrid-th3-debug-c2s" ||
                 arg == "--lvl6-robust-hybrid-th3-debug-evalmod" ||
                 arg == "--lvl6-robust-hybrid-th3-debug-stc" ||
                 arg == "--lvl6-robust-hybrid-th3-debug" ||
                 arg == "--lvl6-compact-hybrid-keygen" ||
                 arg == "--lvl6-compact-hybrid-keygen-next" ||
                 arg == "--lvl6-compact-hybrid-run" ||
                 arg == "--lvl6-compact-hybrid-debug-c2s" ||
                 arg == "--lvl6-compact-hybrid-debug-evalmod" ||
                 arg == "--lvl6-compact-hybrid-debug-stc" ||
                 arg == "--lvl6-compact-hybrid-debug" ||
                 arg == "--lvl6-robust-hybrid-th4-keygen" ||
                 arg == "--lvl6-robust-hybrid-th4-keygen-next" ||
                 arg == "--lvl6-robust-hybrid-th4-run" ||
                 arg == "--lvl6-hybrid-th3-keygen" ||
                 arg == "--lvl6-hybrid-th3-keygen-next" ||
                 arg == "--lvl6-hybrid-th3-run" ||
                 arg == "--lvl6-hybrid-th3-debug-c2s" ||
                 arg == "--lvl6-hybrid-th3-debug-evalmod" ||
                 arg == "--lvl6-hybrid-th3-debug-stc" ||
                 arg == "--lvl6-hybrid-th3-debug" ||
                 arg == "--lvl6-hybrid-debug-c2s" ||
                 arg == "--lvl6-hybrid-debug-evalmod" ||
                 arg == "--lvl6-hybrid-debug-stc" ||
                 arg == "--lvl6-hybrid-debug" ||
                 arg == "--lvl6-debug-evalmod" || arg == "--lvl6-debug-stc" ||
                 arg == "--lvl6-debug" ||
                 arg == "--lvl6-all") {
            if (i + 1 >= args.size()) {
                print_usage(argv[0]);
                return 2;
            }
            saw_action = true;
            const std::filesystem::path key_dir = args[++i];
            if (arg == "--lvl6-keygen") {
                if (run_keygen<Lvl6Schedule>(key_dir, resume,
                                             lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-keygen-next") {
                if (run_keygen_next<Lvl6Schedule>(key_dir,
                                                  lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-run") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_filesystem_bootstrap<Lvl6Schedule>(
                        key_dir, 0.1, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-keygen") {
                if (run_hybrid_keygen<Lvl6Schedule>(key_dir, resume,
                                                    lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-keygen-next") {
                if (run_hybrid_keygen_next<Lvl6Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-run") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_hybrid_filesystem_bootstrap<Lvl6Schedule>(
                        key_dir, 0.1, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-inverse-hybrid-keygen") {
                if (run_hybrid_keygen<Lvl6InverseSchedule>(
                        key_dir, resume, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-inverse-hybrid-keygen-next") {
                if (run_hybrid_keygen_next<Lvl6InverseSchedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-inverse-hybrid-run") {
                print_schedule_report<Lvl6InverseSchedule>("lvl6-inverse",
                                                           &key_dir);
                if (run_hybrid_filesystem_bootstrap<Lvl6InverseSchedule>(
                        key_dir, 0.1, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-inverse-hybrid-debug-c2s") {
                print_schedule_report<Lvl6InverseSchedule>("lvl6-inverse",
                                                           &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<
                        Lvl6InverseSchedule>(
                        key_dir, false, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-inverse-hybrid-debug-evalmod") {
                print_schedule_report<Lvl6InverseSchedule>("lvl6-inverse",
                                                           &key_dir);
                if (run_hybrid_filesystem_evalmod_diagnostics<
                        Lvl6InverseSchedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-inverse-hybrid-debug-stc") {
                print_schedule_report<Lvl6InverseSchedule>("lvl6-inverse",
                                                           &key_dir);
                if (run_hybrid_filesystem_stc_diagnostics<
                        Lvl6InverseSchedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-inverse-hybrid-debug") {
                print_schedule_report<Lvl6InverseSchedule>("lvl6-inverse",
                                                           &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<
                        Lvl6InverseSchedule>(
                        key_dir, true, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-fast-hybrid-keygen" ||
                     arg == "--lvl6-robust-hybrid-th3-keygen") {
                if (run_hybrid_keygen<Lvl6RobustHybridTh3Schedule>(
                        key_dir, resume, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-fast-hybrid-keygen-next" ||
                     arg == "--lvl6-robust-hybrid-th3-keygen-next") {
                if (run_hybrid_keygen_next<Lvl6RobustHybridTh3Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-fast-hybrid-run" ||
                     arg == "--lvl6-robust-hybrid-th3-run") {
                print_schedule_report<Lvl6RobustHybridTh3Schedule>(
                    "lvl6-robust-th3", &key_dir);
                if (run_hybrid_filesystem_bootstrap<
                        Lvl6RobustHybridTh3Schedule>(
                        key_dir, 0.1, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-robust-hybrid-th3-debug-c2s") {
                print_schedule_report<Lvl6RobustHybridTh3Schedule>(
                    "lvl6-robust-th3", &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<
                        Lvl6RobustHybridTh3Schedule>(
                        key_dir, false, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-robust-hybrid-th3-debug-evalmod") {
                print_schedule_report<Lvl6RobustHybridTh3Schedule>(
                    "lvl6-robust-th3", &key_dir);
                if (run_hybrid_filesystem_evalmod_diagnostics<
                        Lvl6RobustHybridTh3Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-robust-hybrid-th3-debug-stc") {
                print_schedule_report<Lvl6RobustHybridTh3Schedule>(
                    "lvl6-robust-th3", &key_dir);
                if (run_hybrid_filesystem_stc_diagnostics<
                        Lvl6RobustHybridTh3Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-robust-hybrid-th3-debug") {
                print_schedule_report<Lvl6RobustHybridTh3Schedule>(
                    "lvl6-robust-th3", &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<
                        Lvl6RobustHybridTh3Schedule>(
                        key_dir, true, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-compact-hybrid-keygen" ||
                     arg == "--lvl6-robust-hybrid-th4-keygen") {
                if (run_hybrid_keygen<Lvl6RobustHybridTh4Schedule>(
                        key_dir, resume, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-compact-hybrid-keygen-next" ||
                     arg == "--lvl6-robust-hybrid-th4-keygen-next") {
                if (run_hybrid_keygen_next<Lvl6RobustHybridTh4Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-compact-hybrid-run" ||
                     arg == "--lvl6-robust-hybrid-th4-run") {
                print_schedule_report<Lvl6RobustHybridTh4Schedule>(
                    "lvl6-robust-th4", &key_dir);
                if (run_hybrid_filesystem_bootstrap<
                        Lvl6RobustHybridTh4Schedule>(
                        key_dir, 0.1, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-compact-hybrid-debug-c2s") {
                print_schedule_report<Lvl6RobustHybridTh4Schedule>(
                    "lvl6-robust-th4", &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<
                        Lvl6RobustHybridTh4Schedule>(
                        key_dir, false, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-compact-hybrid-debug-evalmod") {
                print_schedule_report<Lvl6RobustHybridTh4Schedule>(
                    "lvl6-robust-th4", &key_dir);
                if (run_hybrid_filesystem_evalmod_diagnostics<
                        Lvl6RobustHybridTh4Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-compact-hybrid-debug-stc") {
                print_schedule_report<Lvl6RobustHybridTh4Schedule>(
                    "lvl6-robust-th4", &key_dir);
                if (run_hybrid_filesystem_stc_diagnostics<
                        Lvl6RobustHybridTh4Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-compact-hybrid-debug") {
                print_schedule_report<Lvl6RobustHybridTh4Schedule>(
                    "lvl6-robust-th4", &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<
                        Lvl6RobustHybridTh4Schedule>(
                        key_dir, true, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-th3-keygen") {
                if (run_hybrid_keygen<Lvl6HybridTh3Schedule>(
                        key_dir, resume, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-th3-keygen-next") {
                if (run_hybrid_keygen_next<Lvl6HybridTh3Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-th3-run") {
                print_schedule_report<Lvl6HybridTh3Schedule>("lvl6-th3",
                                                             &key_dir);
                if (run_hybrid_filesystem_bootstrap<Lvl6HybridTh3Schedule>(
                        key_dir, 0.1, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-th3-debug-c2s") {
                print_schedule_report<Lvl6HybridTh3Schedule>("lvl6-th3",
                                                             &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<
                        Lvl6HybridTh3Schedule>(
                        key_dir, false, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-th3-debug-evalmod") {
                print_schedule_report<Lvl6HybridTh3Schedule>("lvl6-th3",
                                                             &key_dir);
                if (run_hybrid_filesystem_evalmod_diagnostics<
                        Lvl6HybridTh3Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-th3-debug-stc") {
                print_schedule_report<Lvl6HybridTh3Schedule>("lvl6-th3",
                                                             &key_dir);
                if (run_hybrid_filesystem_stc_diagnostics<
                        Lvl6HybridTh3Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-th3-debug") {
                print_schedule_report<Lvl6HybridTh3Schedule>("lvl6-th3",
                                                             &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<
                        Lvl6HybridTh3Schedule>(
                        key_dir, true, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-debug-c2s") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<Lvl6Schedule>(
                        key_dir, false, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-debug-evalmod") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_hybrid_filesystem_evalmod_diagnostics<Lvl6Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-debug-stc") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_hybrid_filesystem_stc_diagnostics<Lvl6Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-hybrid-debug") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_hybrid_filesystem_bootstrap_diagnostics<Lvl6Schedule>(
                        key_dir, true, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-debug-c2s") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_filesystem_bootstrap_diagnostics<Lvl6Schedule>(
                        key_dir, false, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-debug-evalmod") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_filesystem_evalmod_diagnostics<Lvl6Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-debug-stc") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_filesystem_stc_diagnostics<Lvl6Schedule>(
                        key_dir, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else if (arg == "--lvl6-debug") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_filesystem_bootstrap_diagnostics<Lvl6Schedule>(
                        key_dir, true, lvl6_sparse_weight) != 0)
                    return 1;
            }
            else {
                if (run_keygen<Lvl6Schedule>(key_dir, resume,
                                             lvl6_sparse_weight) != 0)
                    return 1;
                if (run_filesystem_bootstrap<Lvl6Schedule>(
                        key_dir, 0.1, lvl6_sparse_weight) != 0)
                    return 1;
            }
        }
        else {
            print_usage(argv[0]);
            return 2;
        }
    }

    if (!saw_action) {
        if (run_toy_validation(keep) != 0) return 1;
        print_schedule_report<Lvl6Schedule>("lvl6");
    }
    return 0;
}
