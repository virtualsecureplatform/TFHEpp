#include <algorithm>
#include <chrono>
#include <complex>
#include <filesystem>
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

template <class Schedule>
std::string manifest_status(const std::filesystem::path &root)
{
    const std::filesystem::path manifest_path =
        TFHEpp::CKKSDenseBootstrapKeyDirectoryManifestFile(root);
    if (!std::filesystem::exists(manifest_path)) return "missing";
    try {
        return TFHEpp::CKKSDenseBootstrapKeyDirectoryManifestMatches<Schedule>(
                   root)
                   ? "match"
                   : "mismatch";
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

    std::cout << label << " n=" << P::n << " logQ="
              << Schedule::boot_log_q << " input_logQ="
              << Schedule::input_log_q << " output_logQ="
              << Schedule::output_log_q << '\n';
    std::cout << label << " c2s_levels="
              << Schedule::coeff_to_slot_level_count << " stc_levels="
              << Schedule::slot_to_coeff_level_count << " evalmod_depth="
              << Schedule::evalmod_depth << '\n';
    std::cout << label << " rotation_indices="
              << TFHEpp::CKKSDenseBootstrapRotationKeyUsageCount<Schedule>(
                     usage)
              << "/"
              << TFHEpp::CKKSDenseBootstrapFullGaloisKeyIndexCount<Schedule>()
              << " sparse_rows=" << sparse_rows
              << " streamed_peak_rows=" << streamed_peak_rows << '\n';
    std::cout << label << " sparse_key_bytes=" << sparse_bytes
              << " streamed_peak_bytes=" << streamed_peak_bytes
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
int run_filesystem_bootstrap(const std::filesystem::path &key_dir, double tol)
{
    using P = typename Schedule::Param;

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

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_test_key<P>(*key);

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

    const double keygen_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseBootstrapKeyGenToDirectory<Schedule>(key_dir, *key,
                                                              {0.0, 0});
    });
    TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions resume_options;
    resume_options.overwrite_existing = false;
    TFHEpp::CKKSDenseBootstrapKeyGenToDirectory<Schedule>(
        key_dir, *key, {0.0, 0}, resume_options);

    print_schedule_report<Schedule>("toy", &key_dir);
    std::cout << "toy keygen_ms=" << keygen_ms << '\n';
    if (check_manifest_mismatch_rejected<Schedule>(key_dir) != 0) return 1;
    const int result = run_filesystem_bootstrap<Schedule>(key_dir, 0.02);

    if (!keep_dir) std::filesystem::remove_all(key_dir);
    return result;
}

template <class Schedule>
int run_keygen(const std::filesystem::path &key_dir, bool resume)
{
    using P = typename Schedule::Param;
    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_test_key<P>(*key);

    TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions options;
    options.overwrite_existing = !resume;
    print_schedule_report<Schedule>("keygen-before", &key_dir);
    const double keygen_ms = elapsed_ms([&] {
        TFHEpp::CKKSDenseBootstrapKeyGenToDirectory<Schedule>(key_dir, *key,
                                                              {P::α, 0},
                                                              options);
    });
    print_schedule_report<Schedule>("keygen-after", &key_dir);
    std::cout << "keygen_ms=" << keygen_ms << '\n';
    return 0;
}

void print_usage(const char *program)
{
    std::cerr << "Usage: " << program
              << " [--toy] [--keep] [--lvl6-plan] [--lvl6-keygen DIR]"
                 " [--lvl6-run DIR] [--lvl6-all DIR] [--resume]\n";
}

}  // namespace

int main(int argc, char **argv)
{
    using Lvl6Schedule = TFHEpp::CKKSDenseBootstrapSchedule<TFHEpp::lvl6param>;

    bool saw_action = false;
    bool keep = false;
    bool resume = false;
    std::vector<std::string> args(argv + 1, argv + argc);

    for (std::size_t i = 0; i < args.size(); i++) {
        const std::string &arg = args[i];
        if (arg == "--keep") {
            keep = true;
        }
        else if (arg == "--resume") {
            resume = true;
        }
        else if (arg == "--toy") {
            saw_action = true;
            if (run_toy_validation(keep) != 0) return 1;
        }
        else if (arg == "--lvl6-plan") {
            saw_action = true;
            print_schedule_report<Lvl6Schedule>("lvl6");
        }
        else if (arg == "--lvl6-keygen" || arg == "--lvl6-run" ||
                 arg == "--lvl6-all") {
            if (i + 1 >= args.size()) {
                print_usage(argv[0]);
                return 2;
            }
            saw_action = true;
            const std::filesystem::path key_dir = args[++i];
            if (arg == "--lvl6-keygen") {
                if (run_keygen<Lvl6Schedule>(key_dir, resume) != 0) return 1;
            }
            else if (arg == "--lvl6-run") {
                print_schedule_report<Lvl6Schedule>("lvl6", &key_dir);
                if (run_filesystem_bootstrap<Lvl6Schedule>(key_dir, 0.1) != 0)
                    return 1;
            }
            else {
                if (run_keygen<Lvl6Schedule>(key_dir, resume) != 0) return 1;
                if (run_filesystem_bootstrap<Lvl6Schedule>(key_dir, 0.1) != 0)
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
