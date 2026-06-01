#include <algorithm>
#include <chrono>
#include <complex>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tfhe++.hpp>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

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
        if (entry.is_regular_file()) count++;
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
            Schedule::linear_bsgs_step);
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
                   Schedule::linear_bsgs_step);

    std::cout << label << " n=" << P::n << " logQ="
              << Schedule::boot_log_q << " input_logQ="
              << Schedule::input_log_q << " output_logQ="
              << Schedule::output_log_q << '\n';
    std::cout << label << " c2s_levels="
              << Schedule::coeff_to_slot_level_count << " stc_levels="
              << Schedule::slot_to_coeff_level_count << " evalmod_depth="
              << Schedule::evalmod_depth << '\n';
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
    const double bootstrap_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseBootstrapWithKeyProvider<Schedule>(*output, *input,
                                                            provider);
    });

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    const double decrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q,
                                Schedule::log_delta>(*decoded, *output, *key);
    });
    const double err = max_error<P>(*decoded, *slots);
    std::cout << "encrypt_ms=" << encrypt_ms << '\n';
    std::cout << "bootstrap_ms=" << bootstrap_ms << '\n';
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
    const double bootstrap_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseBootstrapWithKeyProvider<Schedule>(*output, *input,
                                                            provider);
    });

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    const double decrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q,
                                Schedule::log_delta>(*decoded, *output, *key);
    });
    const double err = max_error<P>(*decoded, *slots);
    std::cout << "encrypt_ms=" << encrypt_ms << '\n';
    std::cout << "bootstrap_ms=" << bootstrap_ms << '\n';
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

template <class Schedule>
int run_filesystem_bootstrap_diagnostics(const std::filesystem::path &key_dir,
                                         bool full_pipeline,
                                         std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    if (const int status = validate_filesystem_key_dir<Schedule>(key_dir);
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

    TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
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
            TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule>, true>
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
                 poly, (*expected_real)[i].real()) /
                 Schedule::message_ratio,
             0.0};
        (*expected_imag_eval)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, (*expected_imag)[i].real()) /
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
        TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule>, false>
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
int run_filesystem_evalmod_diagnostics(const std::filesystem::path &key_dir,
                                       std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    if (const int status = validate_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);
    std::cout << "diag_key_sparse_weight=" << sparse_weight << '\n';

    TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
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
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(poly,
                                                               normalized) /
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
int run_filesystem_stc_diagnostics(const std::filesystem::path &key_dir,
                                   std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;

    if (const int status = validate_filesystem_key_dir<Schedule>(key_dir);
        status != 0)
        return status;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_bootstrap_test_key<P>(*key, sparse_weight);
    std::cout << "diag_key_sparse_weight=" << sparse_weight << '\n';

    TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
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
                 poly, (*coeff_to_slot)[i].real()) /
                 Schedule::message_ratio,
             0.0};
        (*imag_eval)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, (*coeff_to_slot)[i].imag()) /
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
        TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule>, false>
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

int run_toy_validation(bool keep_dir)
{
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        TinyDeepMultiLimbCKKSParam, 40, 8, 560, 40, 3, 31, 4, 2, 0, 40, 2>;
    const std::filesystem::path key_dir =
        std::filesystem::temp_directory_path() /
        "tfhepp_ckks_bootstrap_validation_toy";
    std::filesystem::remove_all(key_dir);

    auto key = std::make_unique<TFHEpp::Key<TinyDeepMultiLimbCKKSParam>>();
    fill_test_key<TinyDeepMultiLimbCKKSParam>(*key);

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

    print_schedule_report<Schedule>("toy", &key_dir);
    std::cout << "toy keygen_ms=" << keygen_ms << '\n';
    std::cout << "toy keygen_slices=" << generated_slices << '\n';
    if (check_manifest_mismatch_rejected<Schedule>(key_dir) != 0) return 1;
    const int result = run_filesystem_bootstrap<Schedule>(key_dir, 0.02);

    if (!keep_dir) std::filesystem::remove_all(key_dir);
    return result;
}

template <class Schedule>
int run_keygen_next(const std::filesystem::path &key_dir,
                    std::size_t sparse_weight = 0)
{
    using P = typename Schedule::Param;
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
              << " [--toy] [--keep] [--lvl6-sparse-key H|--lvl6-dense-key]"
                 " [--lvl6-plan] [--lvl6-keygen DIR]"
                 " [--lvl6-keygen-next DIR] [--lvl6-run DIR]"
                 " [--lvl6-hybrid-keygen DIR]"
                 " [--lvl6-hybrid-keygen-next DIR]"
                 " [--lvl6-hybrid-run DIR]"
                 " [--lvl6-debug-modraise]"
                 " [--lvl6-debug-modraise-sparse H]"
                 " [--lvl6-debug-c2s DIR] [--lvl6-debug-evalmod DIR]"
                 " [--lvl6-debug-stc DIR]"
                 " [--lvl6-debug DIR]"
                 " [--lvl6-all DIR] [--resume]\n";
}

}  // namespace

int main(int argc, char **argv)
{
    using Lvl6Schedule = TFHEpp::CKKSDenseBootstrapSchedule<TFHEpp::lvl6param>;
    constexpr std::size_t default_lvl6_sparse_weight = 192;

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
        else if (arg == "--lvl6-plan") {
            saw_action = true;
            print_schedule_report<Lvl6Schedule>("lvl6");
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
