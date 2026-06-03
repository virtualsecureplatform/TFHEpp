#include <algorithm>
#include <complex>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <tfhe++.hpp>

namespace {

struct TinyWorkflowCKKSParam {
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
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 36;
    static constexpr std::uint32_t l̅ₐ = 36;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;
};

using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
    TinyWorkflowCKKSParam, 30, 8, 560, 30, 3, 31, 4, 2, 5, 30, 2>;

template <class P>
void fill_dense_key(TFHEpp::Key<P> &key)
{
    for (std::size_t i = 0; i < P::n; i++) {
        const int v = static_cast<int>(i % 3) - 1;
        key[i] = static_cast<typename P::T>(v);
    }
}

template <class P>
void fill_sparse_key(TFHEpp::Key<P> &key, std::size_t weight)
{
    key.fill(typename P::T{0});
    for (std::size_t i = 0; i < weight; i++) {
        const std::size_t index = (i * 5) % P::n;
        const int sign = (i % 2 == 0) ? 1 : -1;
        key[index] = static_cast<typename P::T>(sign);
    }
}

template <class P>
void fill_lhs_slots(TFHEpp::CKKSSlotVector<P> &slots)
{
    for (std::size_t i = 0; i < P::n / 2; i++) {
        const double re =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 64.0;
        const double im =
            static_cast<double>(static_cast<int>(i % 7) - 3) / 128.0;
        slots[i] = {re, im};
    }
}

template <class P>
void fill_rhs_slots(TFHEpp::CKKSSlotVector<P> &slots)
{
    for (std::size_t i = 0; i < P::n / 2; i++) {
        const double re =
            static_cast<double>(static_cast<int>((3 * i + 1) % 7) - 3) / 96.0;
        const double im =
            static_cast<double>(static_cast<int>((5 * i + 2) % 9) - 4) / 160.0;
        slots[i] = {re, im};
    }
}

template <class P>
double max_error(const TFHEpp::CKKSSlotVector<P> &got,
                 const TFHEpp::CKKSSlotVector<P> &want)
{
    double err = 0.0;
    for (std::size_t i = 0; i < P::n / 2; i++)
        err = std::max(err, std::abs(got[i] - want[i]));
    return err;
}

std::uintmax_t directory_bytes(const std::filesystem::path &root)
{
    std::uintmax_t bytes = 0;
    for (const auto &entry :
         std::filesystem::recursive_directory_iterator(root)) {
        if (entry.is_regular_file()) bytes += entry.file_size();
    }
    return bytes;
}

template <class ScheduleT>
std::size_t seeded_hybrid_streamed_bootstrap_key_estimate_bytes()
{
    TFHEpp::CKKSDenseBootstrapLinearPlan<ScheduleT> linear_plan;
    TFHEpp::CKKSBuildDenseBootstrapLinearPlan<ScheduleT>(linear_plan);
    TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsage<ScheduleT> usage;
    TFHEpp::CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<ScheduleT>(
        usage, linear_plan);
    return TFHEpp::CKKSDenseBootstrapHybridGiantSeededKeyByteEstimate<
        ScheduleT>(usage);
}

}  // namespace

int main()
{
    using P = Schedule::Param;
    constexpr std::uint32_t fresh_log_q =
        Schedule::input_log_q + 2 * Schedule::log_delta;
    using FreshCt = TFHEpp::CKKSCiphertext<P, fresh_log_q, Schedule::log_delta>;
    using ProductCt =
        TFHEpp::CKKSMultResult<P, FreshCt::log_q, FreshCt::log_delta,
                               FreshCt::log_q, FreshCt::log_delta>;
    using PostBootstrapProductCt = TFHEpp::CKKSMultResult<
        P, Schedule::output_log_q, Schedule::log_delta,
        Schedule::output_log_q, Schedule::log_delta>;
    static_assert(Schedule::supports_post_bootstrap_product);
    static_assert(PostBootstrapProductCt::log_q ==
                  Schedule::post_bootstrap_product_log_q);

    const std::filesystem::path root =
        std::filesystem::temp_directory_path() /
        "tfhepp_ckks_bootstrap_workflow";
    const std::filesystem::path bootstrap_key_dir = root / "bootstrap_key";
    const std::filesystem::path eval_key_dir = root / "seeded_external_eval_key";
    const auto eval_key_files =
        TFHEpp::CKKSDenseBootstrapSeededProductEvalKeyFilesInDirectory(
            eval_key_dir);
    const std::filesystem::path encapsulation_key_file =
        eval_key_files.encapsulation_key;
    const std::filesystem::path relin_key_file =
        eval_key_files.product_relin_key;
    const std::filesystem::path post_bootstrap_relin_key_file =
        eval_key_files.post_bootstrap_product_relin_key;
    std::filesystem::remove_all(root);

    auto external_key = std::make_unique<TFHEpp::Key<P>>();
    auto bootstrap_key = std::make_unique<TFHEpp::Key<P>>();
    fill_dense_key<P>(*external_key);
    fill_sparse_key<P>(*bootstrap_key, 2);

    TFHEpp::CKKSDenseBootstrapSeededHybridGiantStreamedKeyGenToDirectory<
        Schedule>(bootstrap_key_dir, *bootstrap_key, {0.0, 0});
    TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions resume_options;
    resume_options.overwrite_existing = false;
    TFHEpp::CKKSDenseBootstrapSeededHybridGiantStreamedKeyGenToDirectory<
        Schedule>(bootstrap_key_dir, *bootstrap_key, {0.0, 0}, resume_options);
    if (!TFHEpp::
            CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryComplete<
                Schedule>(bootstrap_key_dir) ||
        !TFHEpp::
            CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryManifestMatches<
                Schedule>(bootstrap_key_dir)) {
        std::cerr << "workflow streamed seeded bootstrap key directory invalid\n";
        std::filesystem::remove_all(root);
        return 1;
    }
    if (TFHEpp::CKKSDenseBootstrapSeededProductEvalKeyDirectoryComplete(
            eval_key_dir)) {
        std::cerr << "workflow incomplete seeded eval key directory accepted\n";
        std::filesystem::remove_all(root);
        return 1;
    }

    TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions eval_key_options;
    eval_key_options.overwrite_existing = false;
    TFHEpp::CKKSDenseBootstrapSeededProductEvalKeyGenToDirectory<Schedule>(
        eval_key_dir, *external_key, *bootstrap_key, true, {0.0, 0},
        eval_key_options);
    TFHEpp::CKKSDenseBootstrapSeededProductEvalKeyGenToDirectory<Schedule>(
        eval_key_dir, *external_key, *bootstrap_key, true, {0.0, 0},
        eval_key_options);
    if (!TFHEpp::CKKSDenseBootstrapSeededProductEvalKeyDirectoryComplete(
            eval_key_dir)) {
        std::cerr << "workflow seeded eval key directory incomplete\n";
        std::filesystem::remove_all(root);
        return 1;
    }
    if (!TFHEpp::CKKSDenseBootstrapSeededProductEvalKeyDirectoryManifestMatches<
            Schedule>(eval_key_dir)) {
        std::cerr << "workflow seeded eval key manifest mismatch\n";
        std::filesystem::remove_all(root);
        return 1;
    }

    auto lhs = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto rhs = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto chained_expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    fill_lhs_slots<P>(*lhs);
    fill_rhs_slots<P>(*rhs);
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*expected)[i] = (*lhs)[i] * (*rhs)[i];
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*chained_expected)[i] = (*expected)[i] * (*expected)[i];

    auto lhs_ct = std::make_unique<FreshCt>();
    auto rhs_ct = std::make_unique<FreshCt>();
    TFHEpp::ckksSlotEncrypt<P, FreshCt::log_q, FreshCt::log_delta>(
        *lhs_ct, *lhs, *external_key, {0.0, 0});
    TFHEpp::ckksSlotEncrypt<P, FreshCt::log_q, FreshCt::log_delta>(
        *rhs_ct, *rhs, *external_key, {0.0, 0});

    auto bootstrapped = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapProductTimings timings;
    TFHEpp::
        CKKSDenseBootstrapProductWithSeededHybridGiantStreamedFilesystemSeededEvalKeyDirectoryTimed<
            Schedule>(
        *bootstrapped, *lhs_ct, *rhs_ct, bootstrap_key_dir, eval_key_dir,
        timings);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *bootstrapped, *external_key);

    const double err = max_error<P>(*decoded, *expected);

    auto chained_bootstrapped =
        std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapProductTimings chained_timings;
    TFHEpp::
        CKKSDenseBootstrapPostBootstrapProductWithSeededHybridGiantStreamedFilesystemSeededEvalKeyDirectoryTimed<
            Schedule>(
        *chained_bootstrapped, *bootstrapped, *bootstrapped, bootstrap_key_dir,
        eval_key_dir, chained_timings);
    TFHEpp::ckksSlotDecrypt<P, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *chained_bootstrapped, *external_key);
    const double chained_err = max_error<P>(*decoded, *chained_expected);

    std::cout << "workflow product_logQ=" << ProductCt::log_q
              << " normalized_logQ=" << Schedule::input_log_q
              << " output_logQ=" << Schedule::output_log_q << '\n';
    std::cout << "workflow post_bootstrap_product_logQ="
              << PostBootstrapProductCt::log_q
              << " product_bootstrap_slack="
              << Schedule::post_bootstrap_product_slack << '\n';
    std::cout << "workflow streamed_seeded_bootstrap_key_dir="
              << bootstrap_key_dir.string()
              << " bootstrap_key_bytes=" << directory_bytes(bootstrap_key_dir)
              << " estimated_streamed_seeded_hybrid_key_bytes="
              << seeded_hybrid_streamed_bootstrap_key_estimate_bytes<
                     Schedule>()
              << '\n';
    std::cout << "workflow seeded_eval_key_dir=" << eval_key_dir.string()
              << " seeded_encapsulation_key_bytes="
              << std::filesystem::file_size(encapsulation_key_file)
              << " estimated_seeded_encapsulation_key_bytes="
              << TFHEpp::CKKSDenseBootstrapEncapsulationSeededKeyByteEstimate<
                     Schedule>()
              << " seeded_relin_key_bytes="
              << std::filesystem::file_size(relin_key_file)
              << " estimated_seeded_relin_key_bytes="
              << TFHEpp::CKKSSeededRelinKeyByteEstimate<P, ProductCt::log_q>()
              << " seeded_post_bootstrap_relin_key_bytes="
              << std::filesystem::file_size(post_bootstrap_relin_key_file)
              << " estimated_seeded_post_bootstrap_relin_key_bytes="
              << TFHEpp::CKKSSeededRelinKeyByteEstimate<
                     P, PostBootstrapProductCt::log_q>()
              << '\n';
    std::cout << "workflow multiply_ms=" << timings.multiply_ms
              << " bootstrap_ms=" << timings.bootstrap.total_ms()
              << " max_error=" << err << '\n';
    std::cout << "workflow chained_multiply_ms="
              << chained_timings.multiply_ms
              << " chained_bootstrap_ms="
              << chained_timings.bootstrap.total_ms()
              << " chained_max_error=" << chained_err << '\n';

    std::filesystem::remove_all(root);
    return err <= 0.05 && chained_err <= 0.05 ? 0 : 1;
}
