#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <ckks/narrow_packed_file.hpp>
#include <chrono>
#include <fstream>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <string>

namespace {
using narrow_benchmark::finiteMaximum;
using narrow_benchmark::BenchmarkParam;
using narrow_benchmark::LeanBasisParam;
using narrow_benchmark::basisBenchmark;

template <class Refresh>
constexpr std::uint64_t plain37CheckpointMagic()
{
    static_assert(Refresh::auxiliary_log_scale == 28 || Refresh::auxiliary_log_scale == 16);
    // Different gadget-row counts must never be read as the original archive.
    return Refresh::auxiliary_log_scale == 16 ? 0x4e41525243503136ULL : 0x4e41525243503031ULL;
}

template <class Refresh>
std::filesystem::path capturePlain37Failure(
    const std::string &selectors, const std::string &function,
    const TFHEpp::Key<typename Refresh::Param> &dense,
    const TFHEpp::Key<typename Refresh::Param> &sparse,
    const TFHEpp::CKKSSlotVector<typename Refresh::Param> &slots,
    const typename Refresh::Input &input,
    const TFHEpp::CKKSNarrowRefreshKey<Refresh> &key,
    const typename Refresh::Output &output, bool refresh_finished, const char *reason)
{
    using P = typename Refresh::Param;
    using Narrow = typename Refresh::Narrow;
    if (!std::dynamic_pointer_cast<const TFHEpp::CKKSNarrowSelectorFile<Narrow>>(key.narrow.selector_source))
        throw std::runtime_error("failure snapshot requires file-backed selectors");
    // Reserve a new private directory. Never replace an earlier failed state.
    std::filesystem::path directory;
    for (std::size_t index = 0; index < 1024; ++index) {
        const auto candidate = std::filesystem::path(selectors + ".failed-state-" + std::to_string(index));
        if (std::filesystem::create_directory(candidate)) { directory = candidate; break; }
    }
    if (directory.empty()) throw std::runtime_error("cannot reserve failure snapshot directory");
    std::filesystem::permissions(directory, std::filesystem::perms::owner_all,
        std::filesystem::perm_options::replace);
    std::filesystem::create_directory(directory / "selectors");
    // Retain only the known public selector files. Hard links survive the
    // source cache's cleanup without moving the caller's directory or files.
    for (bool tail : {false, true})
        for (std::size_t group = 0; group < (tail ? Narrow::tail_packs : Narrow::blocks); ++group) {
            const auto name = (tail ? "tail-" : "head-") + std::to_string(group) + ".bin";
            std::error_code error;
            std::filesystem::create_hard_link(std::filesystem::path(selectors) / name,
                directory / "selectors" / name, error);
            if (error == std::errc::cross_device_link)
                std::filesystem::copy_file(std::filesystem::path(selectors) / name, directory / "selectors" / name);
            else if (error) throw std::filesystem::filesystem_error("cannot retain failed selectors", error);
        }
    std::ofstream file(directory / "checkpoint.pending", std::ios::binary);
    if (!file) throw std::runtime_error("cannot create failure checkpoint");
    {
        cereal::PortableBinaryOutputArchive archive(file);
        const std::array<std::uint64_t, 9> metadata{plain37CheckpointMagic<Refresh>(), P::n,
            Refresh::input_log_q, Refresh::input_log_delta, Narrow::selector_log_q,
            Narrow::log_delta, Narrow::plain_log_delta, Narrow::output_log_q,
            std::numeric_limits<typename P::T>::digits};
        archive(metadata, function, dense, sparse, slots, input, key);
    }
    file.flush(); if (!file) throw std::runtime_error("cannot flush failure checkpoint");
    file.close(); if (!file) throw std::runtime_error("cannot close failure checkpoint");
    std::ofstream ciphertext(directory / "output.ciphertext", std::ios::binary);
    if (!ciphertext) throw std::runtime_error("cannot create failed output ciphertext");
    { cereal::PortableBinaryOutputArchive archive(ciphertext); archive(output); }
    ciphertext.flush(); if (!ciphertext) throw std::runtime_error("cannot flush failed output ciphertext");
    ciphertext.close(); if (!ciphertext) throw std::runtime_error("cannot close failed output ciphertext");
    std::ofstream note(directory / "failure.txt");
    note << "reason: " << reason << "\nrefresh_completed: " << refresh_finished << '\n';
    note.close(); if (!note) throw std::runtime_error("cannot record failure snapshot reason");
    // A completed checkpoint appears only after all replay inputs are retained.
    std::filesystem::rename(directory / "checkpoint.pending", directory / "checkpoint.bin");
    return directory;
}

// std::max can silently keep a finite prior maximum when its next value is
// NaN. Validation must reject that slot rather than report an unchanged gate.

int validationFinitenessDiagnostic()
{
    if (finiteMaximum(0.0, 0.25) != 0.25 || finiteMaximum(0.25, 0.0) != 0.25)
        throw std::runtime_error("finite validation maximum changed");
    std::size_t rejected = 0;
    for (const auto bad : {std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::infinity(),
                          -std::numeric_limits<double>::infinity()})
        for (bool prior : {false, true}) {
            try { (void)finiteMaximum(prior ? bad : 0.0, prior ? 0.0 : bad); }
            catch (const std::runtime_error &) { rejected++; }
        }
    if (rejected != 6) throw std::runtime_error("nonfinite validation maximum accepted");
    // Several switching diagnostics accumulate long double errors.
    bool wide_rejected = false;
    try { (void)finiteMaximum(0.0L, std::numeric_limits<long double>::quiet_NaN()); }
    catch (const std::runtime_error &) { wide_rejected = true; }
    if (!wide_rejected) throw std::runtime_error("nonfinite long double validation metric accepted");
    std::cout << "finite validation maxima preserved; all seven nonfinite injections rejected\n";
    return 0;
}


struct RetainedUnitStep {
    std::uint32_t input_log_q, output_log_q, log_delta;
    double maximum_error, evaluation_error, rms_error;
};

template <class Refresh>
void reportStorageProfile(std::ostream &report, const char *name)
{
    using P = typename Refresh::Param;
    using Narrow = typename Refresh::Narrow;
    using CompactNarrow = TFHEpp::CKKSNarrowBootstrapSchedule<P, Narrow::input_log_q,
        Narrow::selector_log_q, Narrow::log_delta, Narrow::plain_log_delta, Narrow::blocks,
        Narrow::relinearize_before_rescale, Narrow::auxiliary_log_scale, true>;
    using CompactRefresh = TFHEpp::CKKSNarrowRefreshSchedule<CompactNarrow,
        Refresh::input_log_delta, Refresh::s2c_plain_log_delta, Refresh::bsgs_step>;
    constexpr auto term_bytes = sizeof(TFHEpp::CKKSLinearTransformTerm<P,
        Refresh::evaluation_input_log_q, Refresh::s2c_plain_log_delta>);
    constexpr auto compact_term_bytes = sizeof(typename TFHEpp::CKKSCompactLinearTransformPlan<P,
        Refresh::evaluation_input_log_q, Refresh::evaluation_input_log_delta, Refresh::s2c_plain_log_delta>::Term);
    constexpr auto factor_fft_bytes = TFHEpp::CKKSPreparedAutoKeyFFT<P, Narrow::factor_key_log_q>::spectrum_count * sizeof(TFHEpp::PolynomialInFD<P>);
    constexpr auto output_fft_bytes = TFHEpp::CKKSPreparedAutoKeyFFT<P, Narrow::output_key_log_q>::spectrum_count * sizeof(TFHEpp::PolynomialInFD<P>);
    constexpr auto dd_fft_bytes = []<std::size_t... I>(std::index_sequence<I...>) {
        return (std::size_t{} + ... + sizeof(typename std::tuple_element_t<I, typename Narrow::PreparedRelinKeys::Tuple>::element_type));
    }(std::make_index_sequence<Narrow::tree_depth>{});
    std::uint32_t rotation_mask = 0;
    for (std::size_t ell = 0; ell < std::min(Narrow::blocks, Narrow::tail_blocks_per_pack); ell++)
        rotation_mask |= (2 * (ell + 1) * Narrow::tail_width) % Narrow::slots;
    using NTT = TFHEpp::CKKSNarrowNTTSelector<Narrow>;
    constexpr auto ntt_group_bytes = Narrow::block_width * (P::k + 1) * NTT::digits * P::n * NTT::prime_count * sizeof(std::uint64_t);
    report << "    {\"profile\": \"" << name << "\", \"n\": " << P::n
           << ", \"selector_bits\": " << Narrow::selector_log_q
           << ", \"root_scale_bits\": " << Narrow::log_delta
           << ", \"root_plain_bits\": " << Narrow::plain_log_delta
           << ", \"primary_bits\": " << P::Bgbit << ", \"secondary_bits\": " << P::B̅gbit
           << ", \"raw_setup_bytes\": " << sizeof(TFHEpp::CKKSNarrowRefreshKey<Refresh>)
           << ", \"raw_dd_chain_bytes\": " << sizeof(typename Narrow::RelinKeys)
           << ", \"compact_dd_setup_bytes\": " << sizeof(TFHEpp::CKKSNarrowRefreshKey<CompactRefresh>)
           << ", \"compact_dd_chain_bytes\": " << sizeof(typename CompactNarrow::RelinKeys)
           << ", \"dense_s2c_term_bytes\": " << P::n * term_bytes
           << ", \"compact_s2c_term_bytes\": " << P::n * compact_term_bytes
           << ", \"prepared_factor_auto_bytes\": " << factor_fft_bytes
           << ", \"prepared_output_auto_bytes\": " << output_fft_bytes
           << ", \"prepared_dd_chain_bytes\": " << dd_fft_bytes
           << ", \"tail_rotation_binary_keys\": " << std::popcount(rotation_mask)
           << ", \"prepared_tail_rotation_bytes\": " << std::popcount(rotation_mask) * factor_fft_bytes
           << ", \"bounded_tail_rotation_peak_bytes\": " << (rotation_mask ? factor_fft_bytes : 0)
           << ", \"one_selector_group_ntt_bytes\": " << ntt_group_bytes
           << ", \"expanded_selector_bytes\": " << Narrow::selector_count * sizeof(typename Narrow::Selector)
           << ", \"one_selector_group_bytes\": " << Narrow::block_width * sizeof(typename Narrow::Selector)
           << ", \"packed_key_object_bytes\": " << sizeof(TFHEpp::CKKSNarrowPackedKey<Narrow>)
           << ", \"packed_archive_payload_bytes\": " <<
                1 + sizeof(TFHEpp::CKKSNarrowPackedKey<Narrow>::indicators) +
                    sizeof(TFHEpp::CKKSNarrowPackedKey<Narrow>::substitutions)
           << ", \"packed_two_depth_fft_bytes\": " <<
                2 * TFHEpp::CKKSPreparedAutoKeyFFT<P, TFHEpp::CKKSNarrowPackedKey<Narrow>::log_q>::spectrum_count * sizeof(TFHEpp::PolynomialInFD<P>)
           << ", \"packed_extraction_ciphertext_scratch_bytes\": " <<
                (2 * P::nbit + 4) * sizeof(typename TFHEpp::CKKSNarrowPackedKey<Narrow>::Packed)
           << '}';
}

int storageDiagnostic(const std::string &output_path)
{
    using namespace TFHEpp;
    std::ofstream file;
    if (!output_path.empty()) {
        file.open(output_path);
        if (!file) throw std::runtime_error("cannot open storage diagnostic report");
    }
    std::ostream &report = file.is_open() ? file : std::cout;
    using LL13FBS = CKKSNarrowBootstrapSchedule<BenchmarkParam<13, 4>, 25, 217, 27, 24, 31, true, 24>;
    using LL14 = CKKSNarrowBootstrapSchedule<BenchmarkParam<14, 7>, 48, 387, 49, 34, 31, true, 28>;
    using LL14Plain37 = CKKSNarrowBootstrapSchedule<BenchmarkParam<14, 7>, 48, 385, 48, 37, 31, true, 28>;
    using HT14 = CKKSNarrowBootstrapSchedule<BenchmarkParam<14, 7>, 28, 383, 25, 23, 31, true, 58>;
    using HT15 = CKKSNarrowBootstrapSchedule<BenchmarkParam<15, 14>, 50, 669, 46, 39, 31, true, 212>;
    using LeanHT14 = CKKSNarrowBootstrapSchedule<LeanBasisParam<14, 7>, 28, 383, 25, 23, 31, true, 58>;
    using LeanHT15 = CKKSNarrowBootstrapSchedule<LeanBasisParam<15, 14>, 50, 669, 46, 39, 31, true, 212>;
    report << "{\n  \"scope\": \"Static key, selector, encoded-plan and derived FFT/NTT payloads; excludes workspaces and allocator overhead\",\n  \"profiles\": [\n";
    reportStorageProfile<CKKSNarrowRefreshSchedule<LL13FBS, 21, 20>>(report, "ll13-fbs-key217-dd");
    report << ",\n";
    reportStorageProfile<CKKSNarrowRefreshSchedule<LL14, 44, 37>>(report, "ll14-auxiliary-fbs");
    report << ",\n";
    reportStorageProfile<CKKSNarrowRefreshSchedule<LL14Plain37, 44, 37>>(report, "ll14-fbs-plain37-dd");
    report << ",\n";
    reportStorageProfile<CKKSNarrowRefreshSchedule<HT14, 26, 23>>(report, "ht14-auxiliary");
    report << ",\n";
    reportStorageProfile<CKKSNarrowRefreshSchedule<HT15, 43, 39>>(report, "ht15-auxiliary");
    report << ",\n";
    reportStorageProfile<CKKSNarrowRefreshSchedule<LeanHT14, 26, 23>>(report, "ht14-lean-dd");
    report << ",\n";
    reportStorageProfile<CKKSNarrowRefreshSchedule<LeanHT15, 43, 39>>(report, "ht15-lean-dd");
    report << ",\n";
    using PackedBounded = CKKSNarrowBootstrapSchedule<BenchmarkParam<14, 7>,
        48, 379, 49, 34, 31, false, 0, true>;
    reportStorageProfile<CKKSNarrowRefreshSchedule<PackedBounded, 41, 37, 0, true>>(
        report, "ll14-packed-bounded-file");
    report << "\n  ]\n}\n";
    return 0;
}


// Fresh encrypted ones preserve the refreshed signal through every level.
// Compare raw/prepared DD exactly and measure both target and added error.
template <bool CompactDD, class P, std::uint32_t Q, std::uint32_t Delta>
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#elif defined(_MSC_VER)
__declspec(noinline)
#endif
std::unique_ptr<TFHEpp::CKKSMultResult<P, Q, Delta, Q, Delta>>
retainedUnitProduct(const TFHEpp::CKKSCiphertext<P, Q, Delta> &input,
    const TFHEpp::Key<P> &secret, TFHEpp::CKKSNoise noise,
    const TFHEpp::CKKSSlotVector<P> &target,
    const TFHEpp::CKKSSlotVector<P> &refreshed,
    std::vector<RetainedUnitStep> &steps)
{
    using namespace TFHEpp;
    auto unit = std::make_unique<CKKSCiphertext<P, Q, Delta>>();
    auto ones_storage = std::make_unique<CKKSSlotVector<P>>();
    auto &ones = *ones_storage;
    ones.fill({1, 0});
    ckksSlotEncrypt<P>(*unit, ones, secret, noise);
    auto raw_key = [&] {
        if constexpr (CompactDD) return makeCKKSCompactDDRelinKey<P, Q>(secret, noise);
        else return makeCKKSDDRelinKey<P, Q>(secret, noise);
    }();
    auto prepared = makeCKKSDDRelinKeyFFT(*raw_key);
    using Result = CKKSMultResult<P, Q, Delta, Q, Delta>;
    auto result = std::make_unique<Result>();
    auto fft_result = std::make_unique<Result>();
    CKKSMultRelinThenRescale<P>(*result, input, *unit, *raw_key);
    CKKSMultRelinThenRescale<P>(*fft_result, input, *unit, *prepared);
    if (result->ct != fft_result->ct)
        throw std::runtime_error("retained unit-chain raw/prepared DD mismatch");
    auto got_storage = std::make_unique<CKKSSlotVector<P>>();
    auto &got = *got_storage;
    ckksSlotDecrypt<P>(got, *result, secret);
    double maximum = 0, evaluation = 0, squared = 0;
    for (std::size_t i = 0; i < got.size(); i++) {
        const double error = std::abs(got[i] - target[i]);
        maximum = finiteMaximum(maximum, error);
        evaluation = finiteMaximum(evaluation, std::abs(got[i] - refreshed[i]));
        squared += error * error;
    }
    steps.push_back({Q, Q - Delta, Delta, maximum, evaluation,
                     std::sqrt(squared / got.size())});
    std::cout << "retained unit product " << steps.size() << " Q=" << Q - Delta
              << " delta=" << Delta << " max error=" << maximum << std::endl;
    return result;
}

// Keep each level's arithmetic frame outside the recursive driver. Inlining
// several levels can otherwise retain tens of MiB of stack until all nine
// checks finish, even though the ciphertexts and keys themselves use the heap.
template <std::size_t Remaining, bool CompactDD = false, class P, std::uint32_t Q, std::uint32_t Delta>
void retainedUnitChain(const TFHEpp::CKKSCiphertext<P, Q, Delta> &input,
    const TFHEpp::Key<P> &secret, TFHEpp::CKKSNoise noise,
    const TFHEpp::CKKSSlotVector<P> &target,
    const TFHEpp::CKKSSlotVector<P> &refreshed,
    std::vector<RetainedUnitStep> &steps)
{
    auto result = retainedUnitProduct<CompactDD>(input, secret, noise, target, refreshed, steps);
    if constexpr (Remaining > 1)
        retainedUnitChain<Remaining - 1, CompactDD>(*result, secret, noise, target, refreshed, steps);
}

template <std::uint32_t Nbit, bool HighPrecision = false, bool Functional = false, bool FullDomain = false,
          bool StreamedSelectors = HighPrecision, bool LateRescale = false, bool AuxiliaryProfile = false,
          std::uint32_t ThroughputProfile = 0, bool CompactDD = false, bool CompactS2C = false,
          bool LeanDD = false, bool LL14PlainPrecision = false, bool LL13Functional = false,
          std::uint32_t Plain37AuxiliaryBits = 28>
int benchmark(const std::string &output_path, bool zero_noise = false,
              const std::string &function = "tanh", bool precise = false, bool constant_input = false, bool windowed = false, bool prepared_autos = false, bool prepared_dd = false, bool prepared_rotations = false, bool packed_selectors = false, std::size_t packed_fft_levels = 0,
              const std::string &selector_directory = "", bool bounded_rotations = false,
              std::size_t selector_chunk_size = 0, bool staged_fft = false,
              bool trace_roots = false, const std::string &checkpoint_directory = "",
              bool replay_checkpoint = false, bool failure_capture_test = false)
{
    using namespace TFHEpp;
    static_assert(!FullDomain || Functional);
    static_assert(ThroughputProfile == 0 || ThroughputProfile == 14 || ThroughputProfile == 15);
    static_assert(ThroughputProfile == 0 || (AuxiliaryProfile && !Functional));
    static_assert(!LL14PlainPrecision ||
        (HighPrecision && Functional && FullDomain && AuxiliaryProfile && ThroughputProfile == 0));
    static_assert(!LL13Functional || (!HighPrecision && Functional && FullDomain && AuxiliaryProfile));
    static_assert(Plain37AuxiliaryBits == 28 || (LL14PlainPrecision && Plain37AuxiliaryBits == 16));
    constexpr bool PackedCheckpoint = HighPrecision && !Functional && !AuxiliaryProfile &&
        ThroughputProfile == 0 && StreamedSelectors && CompactDD && CompactS2C && !LeanDD;
    if (LL13Functional && packed_selectors)
        throw std::invalid_argument("LL13 FBS packed setup exceeds the 218-bit key budget");
    if (function != "tanh" && function != "sigmoid" && function != "gelu")
        throw std::invalid_argument("unsupported functional benchmark target");
    if (!checkpoint_directory.empty()) {
        if ((!LL14PlainPrecision && !(PackedCheckpoint && packed_selectors)) || selector_directory.empty())
            throw std::invalid_argument("checkpoints require a supported file-selector profile");
        if (!replay_checkpoint && !std::filesystem::create_directory(checkpoint_directory))
            throw std::runtime_error("checkpoint directory must not exist");
        if (!replay_checkpoint)
            std::filesystem::permissions(checkpoint_directory, std::filesystem::perms::owner_all,
                std::filesystem::perm_options::replace);
    }
    const auto target = [&](double x) {
        if (function == "sigmoid") return 1 / (1 + std::exp(-x));
        if (function == "gelu") return 0.5 * x * (1 + std::erf(x / std::sqrt(2.0)));
        return std::tanh(x);
    };
    const auto target_curvature = [&](double x) {
        if (function == "sigmoid") {
            const double y = target(x);
            return y * (1 - y) * (1 - 2 * y);
        }
        if (function == "gelu")
            return (2 - x * x) * std::exp(-x * x / 2) / std::sqrt(2 * std::acos(-1.0));
        const double y = std::tanh(x);
        return -2 * y * (1 - y * y);
    };
    const std::size_t degree = windowed ? 40 : precise ? (function == "tanh" ? 20 : 12) : 6;
    const std::size_t series_tile_size = Functional ?
        std::min(2 * degree, (bounded_rotations || staged_fft) &&
            (LL14PlainPrecision || Nbit >= (LL13Functional ? 13 : 14)) ?
            (LL14PlainPrecision || LL13Functional ? std::size_t{12} : std::size_t{24}) : 2 * degree) : 8;
    constexpr double message_bound = FullDomain ? 4.0 : Functional ? 1.0 : 0.125;
    constexpr auto storage_limbs = ThroughputProfile == 15 ? 14 :
                                  ThroughputProfile == 14 || HighPrecision ? 7 : 4;
    using P = std::conditional_t<LeanDD, LeanBasisParam<Nbit, storage_limbs>,
                                BenchmarkParam<Nbit, storage_limbs>>;
    // LL13-inspired DD port: q0=2^25, one 20-bit S2C level, h=31.
    // Fixed 20-bit scales use five 20-bit tree rescalings. Q=2^184 leaves
    // Qout=2^64 and the optional packed key lives at 2^(204+nbit), <=2^217.
    // These are port parameters, not aliases for SEAL's Table-2 RNS chain.
    // High-precision schedule selected by Parameter-Selection/narrowCKKSnoise.py:
    // N=16384, q0=48, input scale=41, S2C=37, selector=379, tree scale=49,
    // root plaintext scale=34; output Q=100 and packed-key Q=427 <= 437.
    using OriginalNarrow = std::conditional_t<ThroughputProfile == 14,
        CKKSNarrowBootstrapSchedule<P, 28, 383, 25, 23, 31, true, 58, CompactDD>,
        std::conditional_t<ThroughputProfile == 15,
        CKKSNarrowBootstrapSchedule<P, 50, 669, 46, 39, 31, true, 212, CompactDD>,
        std::conditional_t<AuxiliaryProfile,
        std::conditional_t<HighPrecision,
            CKKSNarrowBootstrapSchedule<P, 48, Functional ? 387 : 379, 49, 34, 31, true, 28, CompactDD>,
            CKKSNarrowBootstrapSchedule<P, 25, 187, 24, 17, 31, true, 28, CompactDD>>, std::conditional_t<HighPrecision,
        CKKSNarrowBootstrapSchedule<P, 48, Functional ? 387 : 379, 49, 34, 31, LateRescale, 0, CompactDD>,
        CKKSNarrowBootstrapSchedule<P, 25, 184, 20, 20, 31, LateRescale, 0, CompactDD>>>>>;
    // Test higher table precision against the full GeLU negative control.
    // Raise its scale while staying within LL14's packed
    // key budget: 385+37+14=436<=437, output385-37-5*48=108.
    using Narrow = std::conditional_t<LL13Functional,
        CKKSNarrowBootstrapSchedule<P, 25, 217, 27, 24, 31, true, 24, CompactDD>,
        std::conditional_t<LL14PlainPrecision,
        CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, true, Plain37AuxiliaryBits, CompactDD>, OriginalNarrow>>;
    // Parameter-Selection retained_square_range selects Qout=108 (selector
    // 387, packed 435 at N=16384) for component outputs bounded by 4.01,
    // with three guard bits after a DD square. Standard refresh keeps Qout=100.
    // Functional inputs use gamma=4 on [-1,1] or gamma=16 on [-4,4].
    // This changes the input scale, leaving selector and DD tree levels intact.
    using Refresh = CKKSNarrowRefreshSchedule<Narrow,
        LL13Functional ? 21 : ThroughputProfile == 14 ? 26 : ThroughputProfile == 15 ? 43 :
        AuxiliaryProfile && !HighPrecision ? 24 : FullDomain ? 44 : Functional ? 46 : HighPrecision ? 41 : 21,
        LL13Functional ? 20 : ThroughputProfile == 14 ? 23 : ThroughputProfile == 15 ? 39 :
        AuxiliaryProfile && !HighPrecision ? 18 : HighPrecision ? 37 : 20, 0, CompactS2C>;
    using Clock = std::chrono::steady_clock;
    auto seconds_since = [](const auto start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    auto dense_storage = std::make_unique<Key<P>>();
    auto sparse_storage = std::make_unique<Key<P>>();
    auto &dense = *dense_storage;
    auto &sparse = *sparse_storage;
    auto input_storage = std::make_unique<typename Refresh::Input>();
    auto &input = *input_storage;
    auto slots_storage = std::make_unique<CKKSSlotVector<P>>();
    auto &slots = *slots_storage;
    if (!replay_checkpoint) {
        std::uniform_int_distribution<int> ternary(-1, 1);
        for (auto &coefficient : dense) coefficient = static_cast<P::T>(ternary(generator));
        CKKSNarrowStructuredKeyGen<Narrow>(sparse);
    }
    std::cout << "n=" << P::n << " h=31 block_width=" << Narrow::block_width
              << " final_block=" << Narrow::block_size(30)
              << " selectors=" << Narrow::selector_count << " auxiliary_primes="
              << CKKSNarrowNTTSelector<Narrow>::prime_count << std::endl;
    std::cout << (replay_checkpoint ? "loading retained refresh setup" :
        "generating S2C, encapsulation, selector, DD, and NTT keys") << std::endl;
    auto key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    const auto key_start = Clock::now();
    const CKKSNoise noise = zero_noise ? CKKSNoise{0, 0} : CKKSNoise{P::α, 0};
    const auto cache_mode = StreamedSelectors ? CKKSNarrowCacheMode::Streamed : CKKSNarrowCacheMode::Prepared;
    const auto checkpoint_magic = [] {
        if constexpr (LL14PlainPrecision) return plain37CheckpointMagic<Refresh>();
        else return std::uint64_t{0x4e415252504b3031ULL};
    };
    if (replay_checkpoint) {
        if constexpr (LL14PlainPrecision || PackedCheckpoint) {
            std::ifstream file(std::filesystem::path(checkpoint_directory) / "checkpoint.bin", std::ios::binary);
            if (!file) throw std::runtime_error("cannot open completed refresh checkpoint");
            cereal::PortableBinaryInputArchive archive(file);
            std::array<std::uint64_t, 9> metadata;
            std::string saved_function;
            archive(metadata, saved_function);
            const std::array<std::uint64_t, 9> expected{checkpoint_magic(), P::n,
                Refresh::input_log_q, Refresh::input_log_delta, Narrow::selector_log_q,
                Narrow::log_delta, Narrow::plain_log_delta, Narrow::output_log_q,
                std::numeric_limits<typename P::T>::digits};
            if (metadata != expected || saved_function != function)
                throw std::runtime_error("refresh checkpoint profile/function mismatch");
            archive(dense, sparse, slots, input, *key);
            if (file.peek() != std::char_traits<char>::eof())
                throw std::runtime_error("refresh checkpoint has trailing data");
            key->narrow.selector_source = std::make_shared<CKKSNarrowSelectorFile<Narrow>>(
                selector_directory, typename CKKSNarrowSelectorFile<Narrow>::Existing{});
            std::cout << "loaded exact checkpoint keys, secrets, ciphertext and selectors" << std::endl;
        }
    }
    else {
        if (packed_selectors || !selector_directory.empty())
            CKKSNarrowRefreshSetupKeyGen<Refresh>(*key, dense, sparse, dense, noise, cache_mode);
        else CKKSNarrowRefreshKeyGen<Refresh>(*key, dense, sparse, dense, noise, cache_mode);
        if (!packed_selectors && !selector_directory.empty())
            CKKSNarrowEncryptSelectorsToFile<Narrow>(key->narrow, sparse, dense, selector_directory, noise);
    }
    double packed_setup_seconds = 0, selector_reconstruction_seconds = 0, selector_error = 0;
    double selector_validation_seconds = 0;
    std::size_t packed_serialized_key_bytes = 0, packed_substitution_cache_bytes = 0;
    if (packed_selectors) {
        if (!replay_checkpoint) {
            std::cout << "generating and round-tripping packed selector key" << std::endl;
            const auto packed_start = Clock::now();
            auto packed = std::make_unique<CKKSNarrowPackedKey<Narrow>>();
            CKKSNarrowPackedKeyGen<Narrow>(*packed, sparse, dense, noise);
            std::stringstream payload;
            { cereal::PortableBinaryOutputArchive archive(payload); archive(*packed); }
            packed_serialized_key_bytes = static_cast<std::size_t>(payload.tellp());
            auto loaded = std::make_unique<CKKSNarrowPackedKey<Narrow>>();
            payload.seekg(0);
            { cereal::PortableBinaryInputArchive archive(payload); archive(*loaded); }
            packed.swap(loaded);
            loaded.reset();
            if (!checkpoint_directory.empty()) {
                CKKSNarrowSavePackedKey<Narrow>(*packed,
                    std::filesystem::path(checkpoint_directory) / "packed.bin");
                std::cout << "retained raw packed key before selector reconstruction" << std::endl;
            }
            CKKSNarrowPreparePackedSubstitutions<Narrow>(*packed, packed_fft_levels);
            for (const auto &cache : packed->substitutions_fft)
                if (cache) packed_substitution_cache_bytes += cache->spectra.size() * sizeof(PolynomialInFD<P>);
            packed_setup_seconds = seconds_since(packed_start);
            std::cout << "reconstructing selectors from public packed key" << std::endl;
            const auto reconstruct_start = Clock::now();
            if (selector_directory.empty())
                CKKSNarrowReconstructSelectors<Narrow>(key->narrow, *packed);
            else CKKSNarrowReconstructSelectorsToFile<Narrow>(key->narrow, *packed, selector_directory);
            selector_reconstruction_seconds = seconds_since(reconstruct_start);
        }
        // Secrets are used only by this validator, never by reconstruction.
        auto wanted = [&](std::size_t x, std::size_t y) -> std::complex<double> {
            if (x >= Narrow::blocks || y >= Narrow::block_size(x)) return {};
            const auto value = sparse[x * Narrow::block_width + y];
            if (value == typename P::T{1}) return {1, 0};
            if (value == -typename P::T{1}) return {0, 1};
            return {};
        };
        std::cout << "validating reconstructed selector groups" << std::endl;
        const auto validation_start = Clock::now();
        std::vector<typename Narrow::Selector> loaded_group;
        for (bool tail : {false, true})
            for (std::size_t group = 0; group < (tail ? Narrow::tail_packs : Narrow::blocks); group++) {
                if (key->narrow.selector_source)
                    key->narrow.selector_source->load(tail, group, loaded_group);
                std::vector<double> errors(Narrow::block_width, 0);
                std::vector<std::exception_ptr> failures(Narrow::block_width);
                // The FFT backend owns thread-local scratch. Each worker
                // decrypts into its own buffer and checks the same full slot set.
#ifdef _OPENMP
#pragma omp parallel if(P::n >= 512)
#endif
                {
                    auto selector_got = std::make_unique<CKKSSlotVector<P>>();
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
                    for (std::size_t y = 0; y < Narrow::block_width; y++) {
                        try {
                            const auto &selector = key->narrow.selector_source ? loaded_group[y] :
                                (tail ? key->narrow.tail : key->narrow.head)[group * Narrow::block_width + y];
                            ckksSlotDecrypt<P>(*selector_got, selector, dense);
                            for (std::size_t i = 0; i < Narrow::slots; i++) {
                                std::complex<double> expected{};
                                if (!tail) {
                                    if (i < Narrow::head_slots) expected = wanted(group, y);
                                }
                                else {
                                    const auto ell = i / (2 * Narrow::tail_width);
                                    if (ell < Narrow::tail_blocks_per_pack &&
                                        i % (2 * Narrow::tail_width) >= Narrow::tail_width)
                                        expected = wanted(group * Narrow::tail_blocks_per_pack + ell, y);
                                }
                                errors[y] = finiteMaximum(errors[y], std::abs((*selector_got)[i] - expected));
                            }
                        }
                        catch (...) { failures[y] = std::current_exception(); }
                    }
                }
                for (std::size_t y = 0; y < Narrow::block_width; y++) {
                    if (failures[y]) std::rethrow_exception(failures[y]);
                    selector_error = finiteMaximum(selector_error, errors[y]);
                }
                const auto completed = (tail ? Narrow::blocks : 0) + group + 1;
                std::cout << "selector validation groups=" << completed << "/"
                          << Narrow::blocks + Narrow::tail_packs
                          << " seconds=" << seconds_since(validation_start) << std::endl;
            }
        selector_validation_seconds = seconds_since(validation_start);
        std::cout << "selector reconstruction seconds=" << selector_reconstruction_seconds
                  << " max error=" << selector_error << std::endl;
    }
    if (staged_fft) CKKSNarrowUseStagedFFT<Narrow>(key->narrow, selector_chunk_size ? selector_chunk_size : 40);
    else {
        if (prepared_autos) CKKSNarrowPrepareAutomorphisms<Narrow>(key->narrow);
        if (prepared_dd) CKKSNarrowPrepareRelinearization<Narrow>(key->narrow);
        if (bounded_rotations) CKKSNarrowUseBoundedTailRotations<Narrow>(key->narrow);
        else if (prepared_rotations) CKKSNarrowPrepareTailRotations<Narrow>(key->narrow);
        if (selector_chunk_size != 0) CKKSNarrowUseSelectorChunks<Narrow>(key->narrow, selector_chunk_size);
    }
    const double key_seconds = seconds_since(key_start);
    std::size_t cache_bytes = 0;
    for (const auto *cache : {&key->narrow.head_ntt, &key->narrow.tail_ntt})
        for (const auto &selector : *cache)
            for (const auto &spectrum : selector.spectra) cache_bytes += spectrum.size() * sizeof(std::uint64_t);
    const std::size_t auto_cache_bytes =
        ((key->narrow.factor_conjugation_fft ? key->narrow.factor_conjugation_fft->spectra.size() : 0) +
         (key->narrow.output_conjugation_fft ? key->narrow.output_conjugation_fft->spectra.size() : 0)) * sizeof(PolynomialInFD<P>);
    std::cout << "key setup seconds=" << key_seconds << " NTT cache bytes=" << cache_bytes << std::endl;
    std::mt19937_64 engine(20262119);
    std::uniform_real_distribution<double> message(-message_bound, message_bound);
    if (!replay_checkpoint) for (auto &slot : slots) slot = {message(engine), message(engine)};
    if (!replay_checkpoint && constant_input) slots.fill({message_bound, message_bound});
    if (!replay_checkpoint && !constant_input) {
        slots.front() = {-message_bound, message_bound};
        slots.back() = {message_bound, -message_bound};
        slots[1] = {message_bound, message_bound};
        slots[2] = {0, 0};
    }
    std::cout << "encrypting benchmark input" << std::endl;
    if (!replay_checkpoint) ckksSlotEncrypt<P>(input, slots, dense, noise);
    if constexpr (LL14PlainPrecision || PackedCheckpoint) {
        if (!checkpoint_directory.empty() && !replay_checkpoint) {
            const auto directory = std::filesystem::path(checkpoint_directory);
            std::ofstream file(directory / "checkpoint.pending", std::ios::binary | std::ios::trunc);
            if (!file) throw std::runtime_error("cannot create refresh checkpoint");
            {
                cereal::PortableBinaryOutputArchive archive(file);
                const std::array<std::uint64_t, 9> metadata{checkpoint_magic(), P::n,
                    Refresh::input_log_q, Refresh::input_log_delta, Narrow::selector_log_q,
                    Narrow::log_delta, Narrow::plain_log_delta, Narrow::output_log_q,
                    std::numeric_limits<typename P::T>::digits};
                archive(metadata, function, dense, sparse, slots, input, *key);
            }
            file.flush(); if (!file) throw std::runtime_error("cannot flush refresh checkpoint");
            file.close(); if (!file) throw std::runtime_error("cannot close refresh checkpoint");
            const auto source = std::dynamic_pointer_cast<const CKKSNarrowSelectorFile<Narrow>>(key->narrow.selector_source);
            if (!source) throw std::runtime_error("checkpoint selector source is not a file cache");
            source->preserve();
            std::filesystem::rename(directory / "checkpoint.pending", directory / "checkpoint.bin");
            std::cout << "published replay checkpoint before refresh evaluation" << std::endl;
        }
    }
    auto decrypted_input_storage = std::make_unique<CKKSSlotVector<P>>();
    auto &decrypted_input = *decrypted_input_storage;
    ckksSlotDecrypt<P>(decrypted_input, input, dense);
    double input_error = 0;
    for (std::size_t i = 0; i < P::n / 2; i++)
        input_error = finiteMaximum(input_error, std::abs(decrypted_input[i] - slots[i]));
    std::cout << "fresh input max error=" << input_error << std::endl;
    auto prepared_storage = std::make_unique<typename Refresh::Coefficients>();
    auto &prepared = *prepared_storage;
    std::cout << "evaluating S2C and sparse encapsulation" << std::endl;
    const auto preparation_start = Clock::now();
    CKKSNarrowPrepareCoefficients<Refresh>(prepared, input, *key);
    const double coefficient_preparation_seconds = seconds_since(preparation_start);
    std::cout << "S2C complete; decrypting coefficient diagnostics" << std::endl;
    auto coefficient_storage = std::make_unique<std::array<std::array<double, P::n>, 2>>();
    auto &decoded_coefficients = (*coefficient_storage)[0];
    auto &normalized_mask = (*coefficient_storage)[1];
    ckksDecrypt<P>(decoded_coefficients, prepared, sparse);
    double s2c_error = 0;
    for (std::size_t i = 0; i < P::n / 2; i++) {
        s2c_error = finiteMaximum(s2c_error, std::abs(decoded_coefficients[i] - slots[i].real()));
        s2c_error = finiteMaximum(s2c_error, std::abs(decoded_coefficients[i + P::n / 2] - slots[i].imag()));
    }
    for (std::size_t i = 0; i < P::n; i++)
        normalized_mask[i] = -static_cast<double>(std::ldexp(
            ckks_detail::levelToLongDouble<P, Narrow::input_log_q>(prepared.ct[0][i]),
            -static_cast<int>(Narrow::input_log_q)));
    auto table_storage = std::make_unique<CKKSSlotVector<P>>(
        narrow_detail::headTable<Narrow>(normalized_mask, 0, 1));
    const auto &table = *table_storage;
    std::cout << "evaluating first packed factor" << std::endl;
    auto factor = std::make_unique<typename Narrow::Factor>();
    {
        std::vector<typename Narrow::Selector> first_group;
        std::vector<CKKSNarrowNTTSelector<Narrow>> first_cache;
        if (key->narrow.selector_source) {
            key->narrow.selector_source->load(false, 0, first_group);
            if (key->narrow.selector_ntt_chunk_size == 0) narrow_detail::prepareSelectors<Narrow>(first_cache, first_group);
        }
        const auto selectors = key->narrow.selector_source ? std::span<const typename Narrow::Selector>(first_group) :
            std::span<const typename Narrow::Selector>(key->narrow.head);
        const auto cache = key->narrow.selector_source ? std::span<const CKKSNarrowNTTSelector<Narrow>>(first_cache) :
            std::span<const CKKSNarrowNTTSelector<Narrow>>(key->narrow.head_ntt);
        auto factors = narrow_detail::packedFactorsBatch<Narrow>({&table, 1}, selectors, cache, 0, key->narrow, key->narrow.selector_ntt_chunk_size);
        *factor = std::move(factors.front());
    }
    const auto positions = narrow_detail::structuredPositions<Narrow>(sparse);
    auto factor_slot_storage = std::make_unique<std::array<CKKSSlotVector<P>, 2>>();
    auto &factor_got = (*factor_slot_storage)[0];
    auto &factor_expected = (*factor_slot_storage)[1];
    ckksSlotDecrypt<P>(factor_got, *factor, dense);
    rotateCKKSSlotVector<P>(factor_expected, table, -positions[0]);
    const bool negative = sparse[positions[0]] == -typename P::T{1};
    double factor_error = 0;
    for (std::size_t i = 0; i < P::n / 2; i++) {
        const auto expected = i >= Narrow::head_slots ? std::complex<double>{} :
            negative ? std::conj(factor_expected[i]) : factor_expected[i];
        factor_error = finiteMaximum(factor_error, std::abs(factor_got[i] - expected));
    }
    std::cout << "S2C max error=" << s2c_error << " first packed factor max error=" << factor_error << std::endl;
    auto output_storage = std::make_unique<typename Refresh::Output>();
    auto &output = *output_storage;
    bool refresh_finished = false, failure_snapshot_attempted = false, failure_snapshot_written = false;
    auto capture_failure = [&](const char *reason) {
        if constexpr (LL14PlainPrecision) {
            if (!checkpoint_directory.empty() || selector_directory.empty() || failure_snapshot_attempted) return;
            failure_snapshot_attempted = true;
            try {
                const auto directory = capturePlain37Failure<Refresh>(selector_directory, function,
                    dense, sparse, slots, input, *key, output, refresh_finished, reason);
                failure_snapshot_written = true;
                std::cerr << "retained failed plain37 state at " << directory << std::endl;
            }
            catch (const std::exception &error) {
                std::cerr << "failed to retain plain37 failure state: " << error.what() << std::endl;
            }
        }
    };
    try {
        std::cout << "evaluating full complex refresh" << std::endl;
        const auto approximation_start = Clock::now();
        CKKSNarrowTrigonometricApproximation approximation;
        if constexpr (Functional) {
            if (windowed) {
                approximation = CKKSBuildNarrowWindowedTrigonometricApproximation(target,
                    [](double x) { const double y = std::tanh(x); return 1 - y * y; },
                    target_curvature, {1, 1, 2, 4, 40}, message_bound, Refresh::gap(), degree);
            }
            else approximation = CKKSBuildNarrowTrigonometricApproximation(
                target, message_bound, Refresh::gap(), degree, function == "gelu" ? 2.0 : 1.0,
                precise ? 65536 : 8192, 40, 2.0, target_curvature, 160.0);
        }
        const double approximation_seconds = seconds_since(approximation_start);
        CKKSNarrowRootObserver<Narrow, Refresh::input_log_delta> root_observer;
        if (trace_roots) {
            root_observer = [&](const auto &request, const auto &root) {
                auto phase = std::make_unique<std::array<double, P::n>>();
                auto decoded = std::make_unique<CKKSSlotVector<P>>();
                ckksDecrypt<P>(*phase, *request.input, sparse);
                ckksSlotDecrypt<P>(*decoded, root, dense);
                double maximum = 0, error = 0;
                for (std::size_t i = 0; i < P::n / 2; i++) {
                    const auto expected = request.weight * std::polar(1.0,
                        2 * std::acos(-1.0) * request.frequency * (*phase)[i] / Refresh::gap());
                    maximum = finiteMaximum(maximum, std::abs((*decoded)[i]));
                    error = finiteMaximum(error, std::abs((*decoded)[i] - expected));
                }
                std::cout << "root diagnostic frequency=" << request.frequency
                          << " weight_abs=" << std::abs(request.weight)
                          << " max_abs=" << maximum << " max_error=" << error << std::endl;
            };
        }
        const auto refresh_start = Clock::now();
        if constexpr (Functional) {
            CKKSNarrowFunctionalRefreshPrepared<Refresh>(output, prepared, *key, approximation, true, series_tile_size,
                trace_roots ? &root_observer : nullptr);
        }
        else CKKSNarrowRefreshPrepared<Refresh>(output, prepared, *key);
        refresh_finished = true;
        const double prepared_evaluation_seconds = seconds_since(refresh_start);
        // Report the complete encrypted pipeline, excluding diagnostic decryption.
        const double refresh_seconds = coefficient_preparation_seconds + prepared_evaluation_seconds;
        if (!checkpoint_directory.empty() && !output_path.empty()) {
            std::ofstream file(output_path + ".ciphertext", std::ios::binary | std::ios::trunc);
            if (!file) throw std::runtime_error("cannot write checkpoint evaluation ciphertext");
            { cereal::PortableBinaryOutputArchive archive(file); archive(output); }
            file.flush(); if (!file) throw std::runtime_error("cannot flush checkpoint evaluation ciphertext");
        }
        auto got_storage = std::make_unique<CKKSSlotVector<P>>();
        auto &got = *got_storage;
        ckksSlotDecrypt<P>(got, output, dense);
        auto refreshed_storage = std::make_unique<CKKSSlotVector<P>>(got);
        const auto &refreshed = *refreshed_storage;
        double maximum_error = 0, squared_error = 0, evaluation_error = 0;
        auto expected_storage = std::make_unique<CKKSSlotVector<P>>(slots);
        auto approximated_storage = std::make_unique<CKKSSlotVector<P>>(slots);
        auto &expected = *expected_storage;
        auto &approximated = *approximated_storage;
        double maximum_target_magnitude = 0, maximum_approximant_magnitude = 0;
        for (std::size_t i = 0; i < P::n / 2; i++) {
            if constexpr (Functional) {
                expected[i] = {target(slots[i].real()), target(slots[i].imag())};
                approximated[i] = {
                    CKKSEvaluateNarrowApproximation(approximation, slots[i].real()),
                    CKKSEvaluateNarrowApproximation(approximation, slots[i].imag())};
                evaluation_error = finiteMaximum(evaluation_error, std::abs(got[i] - approximated[i]));
            }
            maximum_target_magnitude = finiteMaximum(maximum_target_magnitude, std::abs(expected[i]));
            maximum_approximant_magnitude = finiteMaximum(maximum_approximant_magnitude, std::abs(approximated[i]));
            const double error = std::abs(got[i] - expected[i]);
            maximum_error = finiteMaximum(maximum_error, error);
            squared_error += error * error;
        }
        const double precision_bits = -std::log2(maximum_error);
        std::cout << "refresh seconds=" << refresh_seconds << " max error=" << maximum_error
                  << " precision bits=" << precision_bits << std::endl;
        bool post_refresh_square_key_reused = false;
        if constexpr (LL13Functional) {
            // Normalize before squaring: GeLU needs six bits for its complex
            // product range; sigmoid and tanh need three. Both leave Q=31.
            const auto square_normalized = [&]<std::uint32_t Drop>() {
                constexpr auto q = Refresh::output_log_q - Drop;
                constexpr auto delta = Refresh::output_log_delta - Drop;
                auto reduced = std::make_unique<CKKSScaleDownResult<P,
                    Refresh::output_log_q, Refresh::output_log_delta, Drop>>();
                CKKSScaleDown<P, Refresh::output_log_q, Refresh::output_log_delta, Drop>(*reduced, output);
                const auto relin = makeCKKSCompactDDRelinKey<P, q>(dense, noise);
                auto squared = std::make_unique<CKKSMultResult<P, q, delta, q, delta>>();
                CKKSMultRelinThenRescale<P>(*squared, *reduced, *reduced, *relin);
                ckksSlotDecrypt<P>(got, *squared, dense);
            };
            if (function == "gelu") square_normalized.template operator()<6>();
            else square_normalized.template operator()<3>();
        }
        else {
            const auto relin = [&] {
                constexpr auto q = Refresh::output_log_q - (AuxiliaryProfile ? 0 : Refresh::output_log_delta);
                if constexpr (CompactDD) {
                    if constexpr (PackedCheckpoint) {
                        if (!checkpoint_directory.empty()) {
                            const auto path = std::filesystem::path(checkpoint_directory) / "square-key.bin";
                            const std::array<std::uint64_t, 4> expected{
                                0x4e41525253513031ULL, P::n, q,
                                std::numeric_limits<typename P::T>::digits};
                            if (std::filesystem::exists(path)) {
                                auto retained = std::make_unique<CKKSCompactDDRelinKey<P, q>>();
                                auto saved_secret = std::make_unique<Key<P>>();
                                std::ifstream file(path, std::ios::binary);
                                if (!file) throw std::runtime_error("cannot open retained square key");
                                cereal::PortableBinaryInputArchive archive(file);
                                std::array<std::uint64_t, 4> metadata;
                                archive(metadata, *saved_secret);
                                if (metadata != expected || *saved_secret != dense)
                                    throw std::runtime_error("square key profile/secret mismatch");
                                archive(*retained);
                                if (file.peek() != std::char_traits<char>::eof())
                                    throw std::runtime_error("square key has trailing data");
                                post_refresh_square_key_reused = true;
                                return retained;
                            }
                            auto retained = makeCKKSCompactDDRelinKey<P, q>(dense, noise);
                            const auto pending = std::filesystem::path(checkpoint_directory) / "square-key.pending";
                            std::ofstream file(pending, std::ios::binary | std::ios::trunc);
                            if (!file) throw std::runtime_error("cannot create retained square key");
                            { cereal::PortableBinaryOutputArchive archive(file); archive(expected, dense, *retained); }
                            file.flush(); if (!file) throw std::runtime_error("cannot flush retained square key");
                            file.close(); if (!file) throw std::runtime_error("cannot close retained square key");
                            std::filesystem::rename(pending, path);
                            return retained;
                        }
                    }
                    return makeCKKSCompactDDRelinKey<P, q>(dense, noise);
                }
                else return makeCKKSDDRelinKey<P, q>(dense, noise);
            }();
            auto squared_storage = std::make_unique<CKKSMultResult<P, Refresh::output_log_q,
                Refresh::output_log_delta, Refresh::output_log_q, Refresh::output_log_delta>>();
            auto &squared = *squared_storage;
            if constexpr (AuxiliaryProfile) CKKSMultRelinThenRescale<P>(squared, output, output, *relin);
            else CKKSMult<P>(squared, output, output, *relin);
            ckksSlotDecrypt<P>(got, squared, dense);
            if constexpr (PackedCheckpoint) {
                if (!checkpoint_directory.empty() && !output_path.empty()) {
                    std::ofstream file(output_path + ".square.ciphertext", std::ios::binary | std::ios::trunc);
                    if (!file) throw std::runtime_error("cannot write checkpoint square ciphertext");
                    { cereal::PortableBinaryOutputArchive archive(file); archive(squared); }
                    file.flush(); if (!file) throw std::runtime_error("cannot flush checkpoint square ciphertext");
                }
            }
        }
        double product_error = 0, product_evaluation_error = 0;
        for (std::size_t i = 0; i < P::n / 2; i++) {
            product_error = finiteMaximum(product_error, std::abs(got[i] - expected[i] * expected[i]));
            product_evaluation_error = finiteMaximum(product_evaluation_error,
                std::abs(got[i] - approximated[i] * approximated[i]));
        }
        const double tolerance = LL13Functional ? std::exp2(-(function == "gelu" ? 7.17 : function == "tanh" ? 7.89 : 9.04)) : precise ? std::ldexp(1.0, -23) : HighPrecision ? std::ldexp(1.0, -16) : 0.0625;
        const double approximation_bound = Functional ? std::sqrt(2.0) * *approximation.uniform_error_bound : 0;
        const double total_tolerance = precise ? tolerance : tolerance + approximation_bound;
        const double target_square_tolerance = Functional ?
            2 * maximum_target_magnitude * total_tolerance + total_tolerance * total_tolerance + tolerance : tolerance;
        const double evaluation_square_tolerance = Functional ?
            2 * maximum_approximant_magnitude * tolerance + tolerance * tolerance + tolerance : tolerance;
        std::vector<RetainedUnitStep> retained_steps;
        double retained_maximum = 0;
        if constexpr (ThroughputProfile != 0) {
            constexpr std::uint32_t product_delta = ThroughputProfile == 14 ? 23 : 39;
            constexpr auto drop = Refresh::output_log_delta - product_delta;
            auto reduced = std::make_unique<CKKSScaleDownResult<P, Refresh::output_log_q,
                Refresh::output_log_delta, drop>>();
            CKKSScaleDown<P, Refresh::output_log_q, Refresh::output_log_delta, drop>(*reduced, output);
            retainedUnitChain<9, CompactDD>(*reduced, dense, noise, expected, refreshed, retained_steps);
            for (const auto &step : retained_steps)
                retained_maximum = finiteMaximum(retained_maximum, step.maximum_error);
        }
        const bool accuracy_passed = !failure_capture_test && input_error < tolerance && maximum_error < total_tolerance && evaluation_error < tolerance &&
            product_error < target_square_tolerance && product_evaluation_error < evaluation_square_tolerance &&
            retained_maximum < total_tolerance &&
            (!precise || approximation_bound < tolerance) &&
            (!packed_selectors || selector_error < std::ldexp(1.0, -20));
        if (!accuracy_passed) capture_failure(failure_capture_test ? "injected accuracy failure for capture diagnostic" : "accuracy gate failed");
        std::ofstream file;
        if (!output_path.empty()) {
            file.open(output_path);
            if (!file) throw std::runtime_error("cannot open narrow benchmark output");
        }
        std::ostream &report = file.is_open() ? file : std::cout;
        double coefficient_norm = 0;
        for (const auto value : approximation.coefficients) coefficient_norm += 2 * std::abs(value);
        if (!approximation.coefficients.empty()) coefficient_norm -= std::abs(approximation.coefficients[0]);
        report.precision(17);
        report << "{\n  \"accuracy_passed\": " << (accuracy_passed ? "true" : "false") << ",\n"
               << "  \"accuracy_failure_injected\": " << (failure_capture_test ? "true" : "false") << ",\n"
               << "  \"failure_snapshot_attempted\": " << (failure_snapshot_attempted ? "true" : "false") << ",\n"
               << "  \"failure_snapshot_written\": " << (failure_snapshot_written ? "true" : "false") << ",\n"
               << "  \"evaluation_tolerance\": " << tolerance << ",\n"
               << "  \"target_error_tolerance\": " << total_tolerance << ",\n"
               << "  \"function\": \"" << (Functional ? function : "identity") << "\",\n"
               << "  \"packed_selectors\": " << (packed_selectors ? "true" : "false") << ",\n"
               << "  \"file_backed_selectors\": " << (key->narrow.selector_source ? "true" : "false") << ",\n"
               << "  \"packed_payload_bytes\": " << (packed_selectors ? 4 * sizeof(typename CKKSNarrowPackedKey<Narrow>::Packed) : 0) << ",\n"
               << "  \"packed_serialized_key_bytes\": " << packed_serialized_key_bytes << ",\n"
               << "  \"packed_fft_levels\": " << packed_fft_levels << ",\n"
               << "  \"packed_substitution_cache_bytes\": " << packed_substitution_cache_bytes << ",\n"
               << "  \"packed_setup_seconds\": " << packed_setup_seconds << ",\n"
               << "  \"selector_reconstruction_seconds\": " << selector_reconstruction_seconds << ",\n"
               << "  \"selector_validation_seconds\": " << selector_validation_seconds << ",\n"
               << "  \"selector_reconstruction_max_error\": " << selector_error << ",\n"
               << "  \"selector_reconstruction_tolerance\": " << std::ldexp(1.0, -20) << ",\n"
               << "  \"prepared_auto_cache_bytes\": " << auto_cache_bytes << ",\n"
               << "  \"prepared_automorphisms\": " << (key->narrow.factor_conjugation_fft && key->narrow.output_conjugation_fft ? "true" : "false") << ",\n"
               << "  \"prepared_dd\": " << (key->narrow.relin_fft ? "true" : "false") << ",\n"
               << "  \"compact_dd\": " << (CompactDD ? "true" : "false") << ",\n"
               << "  \"compact_s2c\": " << (CompactS2C ? "true" : "false") << ",\n"
               << "  \"bounded_rotation_fft\": " << (key->narrow.bounded_tail_rotation_fft ? "true" : "false") << ",\n"
               << "  \"staged_fft_caches\": " << (staged_fft ? "true" : "false") << ",\n"
               << "  \"selector_ntt_chunk_size\": " << key->narrow.selector_ntt_chunk_size << ",\n"
               << "  \"dd_chain_storage_bytes\": " << sizeof(typename Narrow::RelinKeys) << ",\n"
               << "  \"relinearize_before_rescale\": " << (Narrow::relinearize_before_rescale ? "true" : "false") << ",\n"
               << "  \"post_refresh_product_relinearize_before_rescale\": " << (AuxiliaryProfile ? "true" : "false") << ",\n"
               << "  \"auxiliary_switch_bits\": " << Narrow::auxiliary_log_scale << ",\n"
               << "  \"factor_key_log_q\": " << Narrow::factor_key_log_q << ",\n"
               << "  \"structured_switch_key_log_q\": " << Refresh::evaluation_coefficient_log_q << ",\n"
               << "  \"s2c_input_key_log_q\": " << Refresh::evaluation_input_log_q << ",\n"
               << "  \"prepared_dd_cache_bytes\": " << (key->narrow.relin_fft ? key->narrow.relin_fft->byteSize() : 0) << ",\n"
               << "  \"prepared_rotations\": " << (key->narrow.tail_rotations_fft ? "true" : "false") << ",\n"
               << "  \"prepared_rotation_cache_bytes\": " << (key->narrow.tail_rotations_fft ? key->narrow.tail_rotations_fft->byteSize() : 0) << ",\n"
               << "  \"approximation_method\": \"" << (windowed ? "windowed-remez" : "direct-remez") << "\",\n"
               << "  \"input_pattern\": \"" << (constant_input ? "constant-endpoint" : "random-with-endpoints") << "\",\n"
               << "  \"series_tile_size\": " << series_tile_size << ",\n"
               << "  \"root_diagnostics\": " << (trace_roots ? "true" : "false") << ",\n"
               << "  \"checkpoint_written\": " << (!checkpoint_directory.empty() && !replay_checkpoint ? "true" : "false") << ",\n"
               << "  \"checkpoint_replayed\": " << (replay_checkpoint ? "true" : "false") << ",\n"
               << "  \"refresh_keys_reused\": " << (replay_checkpoint ? "true" : "false") << ",\n"
               << "  \"post_refresh_square_key_reused\": " << (post_refresh_square_key_reused ? "true" : "false") << ",\n"
               << "  \"degree\": " << (Functional ? degree : 0) << ",\n"
               << "  \"coefficient_l1_norm\": " << coefficient_norm << ",\n"
               << "  \"approximation_sampled_error\": " << approximation.max_sampled_error << ",\n"
               << "  \"evaluation_max_error\": " << evaluation_error << ",\n"
               << "  \"approximation_uniform_bound\": " << approximation.uniform_error_bound.value_or(0) << ",\n"
               << "  \"message_bound\": " << message_bound << ",\n"
               << "  \"paper\": \"2026/2119\",\n  \"profile\": \""
               << (LL13Functional ? "ll13-fbs-key217-dd" : LL14PlainPrecision ?
                   (Plain37AuxiliaryBits == 16 ? "ll14-fbs-plain37-aux16-dd" : "ll14-fbs-plain37-dd") :
                   ThroughputProfile == 14 ? "ht14-inspired-dd" : ThroughputProfile == 15 ? "ht15-inspired-dd" : HighPrecision ? (Functional ? "ll14-fbs-range-dd-streamed" : "ll14-inspired-dd-streamed") : "ll13-inspired-dd") << "\",\n"
               << "  \"security_validated\": false,\n  \"n\": " << P::n
               << ",\n  \"h\": 31,\n  \"source_log_q\": " << Narrow::input_log_q
               << ",\n  \"selector_log_q\": " << Narrow::selector_log_q << ",\n"
               << "  \"packed_log_q\": " << CKKSNarrowPackedKey<Narrow>::log_q
               << ",\n  \"input_log_q\": " << Refresh::input_log_q
               << ",\n  \"input_log_delta\": " << Refresh::input_log_delta
               << ",\n  \"output_log_q\": " << Refresh::output_log_q
               << ",\n  \"output_log_delta\": " << Refresh::output_log_delta
               << ",\n  \"root_plain_log_delta\": " << Narrow::plain_log_delta
               << ",\n  \"s2c_plain_log_delta\": " << Refresh::s2c_plain_log_delta
               << ",\n  \"storage_bits\": " << std::numeric_limits<typename P::T>::digits
               << ",\n  \"cache_mode\": \"" << (StreamedSelectors ? "streamed" : "prepared") << "\""
               << ",\n  \"gaussian_sigma\": " << (zero_noise ? 0.0 : 3.2)
               << ",\n  \"gadget_bits\": " << P::Bgbit
               << ",\n  \"secondary_digit_bits\": " << P::B̅gbit
               << ",\n  \"lean_dd\": " << (LeanDD ? "true" : "false") << ",\n"
               << "  \"tree_depth\": " << Narrow::tree_depth << ",\n  \"ntt_cache_bytes\": " << cache_bytes
               << ",\n  \"key_setup_seconds\": " << key_seconds << ",\n  \"refresh_seconds\": " << refresh_seconds
               << ",\n  \"coefficient_preparation_seconds\": " << coefficient_preparation_seconds
               << ",\n  \"prepared_evaluation_seconds\": " << prepared_evaluation_seconds
               << ",\n  \"diagnostic_preparation_reused\": true"
               << ",\n  \"approximation_seconds\": " << approximation_seconds
               << ",\n  \"fresh_input_max_error\": " << input_error
               << ",\n  \"s2c_max_error\": " << s2c_error << ",\n  \"first_packed_factor_max_error\": " << factor_error
               << ",\n  \"max_error\": " << maximum_error << ",\n  \"rms_error\": " << std::sqrt(squared_error / (P::n / 2))
               << ",\n  \"precision_bits\": " << precision_bits << ",\n  \"post_refresh_product_max_error\": " << product_error
               << ",\n  \"post_refresh_product_evaluation_max_error\": " << product_evaluation_error
               << ",\n  \"retained_unit_product_levels\": " << retained_steps.size()
               << ",\n  \"retained_unit_product_max_error\": " << retained_maximum
               << ",\n  \"retained_unit_products\": [";
        for (std::size_t i = 0; i < retained_steps.size(); i++) {
            const auto &step = retained_steps[i];
            if (i) report << ',';
            report << "\n    {\"input_log_q\": " << step.input_log_q
                   << ", \"output_log_q\": " << step.output_log_q
                   << ", \"log_delta\": " << step.log_delta
                   << ", \"max_error\": " << step.maximum_error
                   << ", \"evaluation_max_error\": " << step.evaluation_error
                   << ", \"rms_error\": " << step.rms_error << '}';
        }
        report << "\n  ]";
        if constexpr (LL13Functional) {
            report << ",\n  \"maximum_used_evaluation_key_bits\": 217"
                   << ",\n  \"optional_packed_setup_within_key_budget\": false"
                   << ",\n  \"retained_square_normalization_bits\": " << (function == "gelu" ? 6 : 3)
                   << ",\n  \"retained_square_log_q\": 31"
                   << ",\n  \"retained_square_log_delta\": " << (function == "gelu" ? 21 : 24);
        }
        report << "\n}\n";
        return accuracy_passed ? 0 : 1;
    }
    catch (...) {
        capture_failure("exception during functional refresh or validation; output may be partial");
        throw;
    }
}

int auxiliarySwitchDiagnostic(const std::string &output_path)
{
    using namespace TFHEpp;
    using P = BenchmarkParam<13>;
    constexpr std::uint32_t q = 45, delta = 21, auxiliary = 28, samples = 32;
    auto secret = std::make_unique<Key<P>>();
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : *secret) value = static_cast<P::T>(ternary(generator));
    auto ordinary_raw = std::make_unique<CKKSAutoKey<P, q>>();
    auto auxiliary_raw = std::make_unique<CKKSAutoKey<P, q + auxiliary>>();
    CKKSAutoKeyGen<P, q>(*ordinary_raw, 2 * P::n - 1, *secret, {P::α, 0});
    CKKSAutoKeyGen<P, q + auxiliary>(*auxiliary_raw, 2 * P::n - 1, *secret, {P::α, 0});
    auto ordinary_key = makeCKKSPreparedAutoKeyFFT<P, q>(*ordinary_raw);
    auto auxiliary_key = makeCKKSPreparedAutoKeyFFT<P, q + auxiliary>(*auxiliary_raw);
    auto input = std::make_unique<CKKSCiphertext<P, q, delta>>();
    auto ordinary = std::make_unique<CKKSCiphertext<P, q, delta>>();
    auto scaled = std::make_unique<CKKSCiphertext<P, q, delta>>();
    auto reference = std::make_unique<CKKSCiphertext<P, q, delta>>();
    long double ordinary_energy = 0, auxiliary_energy = 0, ordinary_max = 0, auxiliary_max = 0;
    bool exact = true;
    std::uniform_real_distribution<double> message(-0.125, 0.125);
    for (std::size_t sample = 0; sample < samples; sample++) {
        CKKSSlotVector<P> slots, before, got_ordinary, got_auxiliary;
        for (auto &value : slots) value = {message(generator), message(generator)};
        ckksSlotEncrypt<P>(*input, slots, *secret, {P::α, 0});
        ckksSlotDecrypt<P>(before, *input, *secret);
        CKKSEvalAuto<P, q>(ordinary->ct, input->ct, 2 * P::n - 1, *ordinary_key);
        CKKSEvalAutoWithAuxiliaryScale<P, q, delta, auxiliary>(*scaled, *input, 2 * P::n - 1, *auxiliary_key);
        CKKSEvalAutoWithAuxiliaryScale<P, q, delta, auxiliary>(*reference, *input, 2 * P::n - 1, *auxiliary_raw);
        exact &= reference->ct == scaled->ct;
        ckksSlotDecrypt<P>(got_ordinary, *ordinary, *secret);
        ckksSlotDecrypt<P>(got_auxiliary, *scaled, *secret);
        for (std::size_t i = 0; i < before.size(); i++) {
            const auto first = got_ordinary[i] - std::conj(before[i]);
            const auto second = got_auxiliary[i] - std::conj(before[i]);
            ordinary_energy += std::norm(first);
            auxiliary_energy += std::norm(second);
            ordinary_max = finiteMaximum(ordinary_max, static_cast<long double>(std::abs(first)));
            auxiliary_max = finiteMaximum(auxiliary_max, static_cast<long double>(std::abs(second)));
        }
    }
    std::ostringstream report;
    report.precision(17);
    report << "{\n  \"diagnostic\": \"dd-auxiliary-scale-switch\",\n"
           << "  \"n\": " << P::n << ",\n  \"input_output_log_q\": " << q
           << ",\n  \"log_delta\": " << delta << ",\n  \"auxiliary_scale_bits\": " << auxiliary
           << ",\n  \"evaluation_key_log_q\": " << q + auxiliary
           << ",\n  \"sigma\": 3.2,\n  \"samples\": " << samples
           << ",\n  \"all_raw_prepared_equal\": " << (exact ? "true" : "false")
           << ",\n  \"ordinary_switch_rms\": " << std::sqrt(ordinary_energy / (samples * (P::n / 2)))
           << ",\n  \"ordinary_switch_max\": " << ordinary_max
           << ",\n  \"auxiliary_switch_rms\": " << std::sqrt(auxiliary_energy / (samples * (P::n / 2)))
           << ",\n  \"auxiliary_switch_max\": " << auxiliary_max << "\n}\n";
    std::cout << report.str();
    if (!output_path.empty()) {
        std::ofstream file(output_path);
        if (!file) throw std::runtime_error("cannot open auxiliary switch diagnostic output");
        file << report.str();
    }
    return exact && std::isfinite(static_cast<double>(auxiliary_energy)) ? 0 : 1;
}

// Compare DD multiplication orders at identical input/output levels.
// The experimental path keeps tensor products modulo Q, switches s^2 at Q,
// then rescales the two-component ciphertext. Production CKKSMult is unchanged.
int multiplicationOrderDiagnostic(const std::string &output_path)
{
    using namespace TFHEpp;
    using P = BenchmarkParam<13>;
    constexpr std::uint32_t q = 64, delta = 20, samples = 8;
    using Input = CKKSCiphertext<P, q, delta>;
    using Output = CKKSMultResult<P, q, delta, q, delta>;
    auto secret = std::make_unique<Key<P>>();
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : *secret) value = static_cast<P::T>(ternary(generator));
    std::uniform_real_distribution<double> message(-0.125, 0.125);
    struct Statistics { long double early_energy{}, late_energy{}, early_max{}, late_max{}; };
    std::array<Statistics, 2> stats;
    for (std::size_t mode = 0; mode < 2; mode++) {
        const CKKSNoise noise = mode ? CKKSNoise{P::α, 0} : CKKSNoise{0, 0};
        const auto early_raw = makeCKKSDDRelinKey<P, q - delta>(*secret, noise);
        const auto late_raw = makeCKKSDDRelinKey<P, q>(*secret, noise);
        const auto early_key = makeCKKSDDRelinKeyFFT(*early_raw);
        const auto late_key = makeCKKSDDRelinKeyFFT(*late_raw);
        auto lhs = std::make_unique<Input>(), rhs = std::make_unique<Input>();
        auto early = std::make_unique<Output>(), late = std::make_unique<Output>();
        for (std::size_t sample = 0; sample < samples; sample++) {
            CKKSSlotVector<P> a, b, before_a, before_b, got_early, got_late;
            for (auto &value : a) value = {message(generator), message(generator)};
            for (auto &value : b) value = {message(generator), message(generator)};
            ckksSlotEncrypt<P>(*lhs, a, *secret, noise);
            ckksSlotEncrypt<P>(*rhs, b, *secret, noise);
            ckksSlotDecrypt<P>(before_a, *lhs, *secret);
            ckksSlotDecrypt<P>(before_b, *rhs, *secret);
            CKKSMult<P>(*early, *lhs, *rhs, *early_key);
            CKKSMultRelinThenRescale<P>(*late, *lhs, *rhs, *late_key);
            ckksSlotDecrypt<P>(got_early, *early, *secret);
            ckksSlotDecrypt<P>(got_late, *late, *secret);
            for (std::size_t i = 0; i < P::n / 2; i++) {
                const auto expected = before_a[i] * before_b[i];
                const auto early_error = got_early[i] - expected;
                const auto late_error = got_late[i] - expected;
                stats[mode].early_energy += std::norm(early_error);
                stats[mode].late_energy += std::norm(late_error);
                stats[mode].early_max = finiteMaximum(stats[mode].early_max, static_cast<long double>(std::abs(early_error)));
                stats[mode].late_max = finiteMaximum(stats[mode].late_max, static_cast<long double>(std::abs(late_error)));
            }
        }
    }
    std::ostringstream report;
    report.precision(17);
    report << "{\n  \"diagnostic\": \"dd-relinearization-rescale-order\",\n"
           << "  \"n\": " << P::n << ",\n  \"input_log_q\": " << q
           << ",\n  \"output_log_q\": " << q - delta << ",\n  \"log_delta\": " << delta
           << ",\n  \"samples_per_noise_mode\": " << samples << ",\n  \"modes\": [\n";
    for (std::size_t mode = 0; mode < stats.size(); mode++) {
        if (mode) report << ",\n";
        report << "    {\"sigma\": " << (mode ? 3.2 : 0)
               << ", \"early_rescale_rms\": " << std::sqrt(stats[mode].early_energy / (samples * (P::n / 2)))
               << ", \"early_rescale_max\": " << stats[mode].early_max
               << ", \"late_rescale_rms\": " << std::sqrt(stats[mode].late_energy / (samples * (P::n / 2)))
               << ", \"late_rescale_max\": " << stats[mode].late_max << "}";
    }
    report << "\n  ]\n}\n";
    std::cout << report.str();
    if (!output_path.empty()) {
        std::ofstream file(output_path);
        if (!file) throw std::runtime_error("cannot open multiplication diagnostic output");
        file << report.str();
    }
    return std::isfinite(static_cast<double>(stats[0].late_energy + stats[1].late_energy)) ? 0 : 1;
}

// Isolate one public automorphism switch from input encryption/encoding error.
// Reuse one noisy key over independent ciphertext masks to expose persistent
// slot error that an average-case variance screen can hide.
int switchNoiseDiagnostic(const std::string &output_path)
{
    using namespace TFHEpp;
    using P = BenchmarkParam<13>;
    constexpr std::uint32_t q = 45, delta = 21, samples = 32;
    constexpr auto rows = CKKSKeySwitchRowCountForLevel<P, q>();
    constexpr auto first = CKKSKeySwitchFirstRowForLevel<P, q>();
    auto secret = std::make_unique<Key<P>>();
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &value : *secret) value = static_cast<P::T>(ternary(generator));
    auto raw_key = std::make_unique<CKKSAutoKey<P, q>>();
    CKKSAutoKeyGen<P, q>(*raw_key, 2 * P::n - 1, *secret, {P::α, 0});
    auto prepared = makeCKKSPreparedAutoKeyFFT<P, q>(*raw_key);
    auto input = std::make_unique<CKKSCiphertext<P, q, delta>>();
    auto output = std::make_unique<CKKSCiphertext<P, q, delta>>();
    auto reference = std::make_unique<CKKSCiphertext<P, q, delta>>();
    auto mask = std::make_unique<Polynomial<P>>();
    auto digits = std::make_unique<std::array<Polynomial<P>, rows>>();
    std::vector<std::complex<double>> error_sum(P::n / 2), first_error(P::n / 2);
    std::array<long double, rows> digit_sum{}, balanced_digit_sum{};
    std::vector<std::complex<double>> balanced_error_sum(P::n / 2);
    auto balanced_digits = std::make_unique<std::array<Polynomial<P>, rows>>();
    auto switched = std::make_unique<TRLWE<P>>();
    auto product = std::make_unique<Polynomial<P>>();
    long double balanced_energy = 0, balanced_maximum = 0;
    bool balanced_recomposition_exact = true;
    long double energy = 0, maximum = 0, correlation_numerator = 0;
    long double first_energy = 0, next_energy = 0;
    std::size_t maximum_slot = 0;
    bool raw_prepared_equal = true;
    std::uniform_real_distribution<double> message(-0.125, 0.125);
    for (std::size_t sample = 0; sample < samples; sample++) {
        CKKSSlotVector<P> slots, before, after;
        for (auto &value : slots) value = {message(generator), message(generator)};
        ckksSlotEncrypt<P>(*input, slots, *secret, {P::α, 0});
        ckksSlotDecrypt<P>(before, *input, *secret);
        CKKSEvalAuto<P, q>(output->ct, input->ct, 2 * P::n - 1, *prepared);
        if (sample == 0) {
            CKKSEvalAuto<P, q>(reference->ct, input->ct, 2 * P::n - 1, *raw_key);
            raw_prepared_equal = reference->ct == output->ct;
        }
        ckksSlotDecrypt<P>(after, *output, *secret);
        for (std::size_t i = 0; i < P::n / 2; i++) {
            const auto error = after[i] - std::conj(before[i]);
            error_sum[i] += error;
            energy += std::norm(error);
            if (std::abs(error) > maximum) { maximum = std::abs(error); maximum_slot = i; }
            if (sample == 0) { first_error[i] = error; first_energy += std::norm(error); }
            if (sample == 1) {
                correlation_numerator += std::real(error * std::conj(first_error[i]));
                next_energy += std::norm(error);
            }
        }
        Automorphism<P>(*mask, input->ct[0], 2 * P::n - 1);
        auto centered = std::make_unique<Polynomial<P>>();
        ckks_detail::centeredPolynomialAtLevel<P, q>(*centered, *mask);
        ckks_detail::baseBbarDecomposePolynomialRows<P, first>(*digits, *centered);
        for (std::size_t row = 0; row < rows; row++)
            for (const auto value : (*digits)[row])
                digit_sum[row] += ckks_detail::levelToLongDouble<P, 256>(value);
        // Diagnostic intervention: distribute +/-half ties using a public
        // deterministic hash. Carry from lower digits preserves the exact
        // torus value. No production decomposition path uses this experiment.
        for (std::size_t i = 0; i < P::n; i++) {
            auto remaining = (*centered)[i];
            P::T recomposed{};
            for (std::size_t row = rows; row-- > 0;) {
                const auto shift = 256 - (first + row + 1) * P::B̅gbit;
                int digit = static_cast<std::uint64_t>(remaining >> shift) & 15;
                std::uint64_t hash = i + row * P::n + 0x9e3779b97f4a7c15ULL;
                hash = (hash ^ (hash >> 30)) * 0xbf58476d1ce4e5b9ULL;
                hash = (hash ^ (hash >> 27)) * 0x94d049bb133111ebULL;
                hash ^= hash >> 31;
                if (digit > 8 || (digit == 8 && (hash & 1))) digit -= 16;
                (*balanced_digits)[row][i] = P::T{digit};
                const auto term = P::T{digit} << shift;
                recomposed += term;
                remaining -= term;
                balanced_digit_sum[row] += digit;
            }
            balanced_recomposition_exact &= recomposed == (*centered)[i];
        }
        for (auto &component : *switched) component.fill(P::T{});
        for (std::size_t row = 0; row < rows; row++)
            for (std::size_t component = 0; component <= P::k; component++) {
                ckks_detail::polyMulTorusByBbarDigit<P, q>(*product,
                    (*raw_key)[0][row][component], (*balanced_digits)[row]);
                for (std::size_t i = 0; i < P::n; i++)
                    (*switched)[component][i] += (*product)[i];
            }
        for (auto &component : reference->ct) component.fill(P::T{});
        Automorphism<P>(reference->ct[P::k], input->ct[P::k], 2 * P::n - 1);
        CKKSSubTRLWEInPlace<P, q>(reference->ct, *switched);
        ckksSlotDecrypt<P>(after, *reference, *secret);
        for (std::size_t i = 0; i < P::n / 2; i++) {
            const auto error = after[i] - std::conj(before[i]);
            balanced_error_sum[i] += error;
            balanced_energy += std::norm(error);
            balanced_maximum = finiteMaximum(balanced_maximum, static_cast<long double>(std::abs(error)));
        }
    }
    long double persistent_energy = 0, maximum_mean = 0;
    for (const auto value : error_sum) {
        persistent_energy += std::norm(value / double(samples));
        maximum_mean = finiteMaximum(maximum_mean, static_cast<long double>(std::abs(value / double(samples))));
    }
    long double balanced_persistent_energy = 0, balanced_maximum_mean = 0;
    for (const auto value : balanced_error_sum) {
        balanced_persistent_energy += std::norm(value / double(samples));
        balanced_maximum_mean = finiteMaximum(balanced_maximum_mean,
            static_cast<long double>(std::abs(value / double(samples))));
    }
    std::ostringstream report;
    report.precision(17);
    report << "{\n  \"diagnostic\": \"reused-noisy-conjugation-key\",\n"
           << "  \"n\": " << P::n << ",\n  \"log_q\": " << q
           << ",\n  \"log_delta\": " << delta << ",\n  \"samples\": " << samples
           << ",\n  \"sigma\": 3.2,\n  \"raw_prepared_equal\": " << (raw_prepared_equal ? "true" : "false")
           << ",\n  \"switch_rms_error\": " << std::sqrt(energy / (samples * (P::n / 2)))
           << ",\n  \"switch_max_error\": " << maximum
           << ",\n  \"maximum_slot\": " << maximum_slot
           << ",\n  \"mean_error_rms\": " << std::sqrt(persistent_energy / (P::n / 2))
           << ",\n  \"mean_error_max\": " << maximum_mean
           << ",\n  \"first_two_error_correlation\": " << correlation_numerator / std::sqrt(first_energy * next_energy)
           << ",\n  \"primary_digit_means\": [";
    for (std::size_t row = 0; row < rows; row++) {
        if (row) report << ", ";
        report << digit_sum[row] / (samples * P::n);
    }
    report << "],\n  \"experimental_balanced_recomposition_exact\": " << (balanced_recomposition_exact ? "true" : "false")
           << ",\n  \"experimental_balanced_switch_rms\": " << std::sqrt(balanced_energy / (samples * (P::n / 2)))
           << ",\n  \"experimental_balanced_switch_max\": " << balanced_maximum
           << ",\n  \"experimental_balanced_mean_rms\": " << std::sqrt(balanced_persistent_energy / (P::n / 2))
           << ",\n  \"experimental_balanced_mean_max\": " << balanced_maximum_mean
           << ",\n  \"experimental_balanced_digit_means\": [";
    for (std::size_t row = 0; row < rows; row++) {
        if (row) report << ", ";
        report << balanced_digit_sum[row] / (samples * P::n);
    }
    report << "]\n}\n";
    std::cout << report.str();
    if (!output_path.empty()) {
        std::ofstream output_file(output_path);
        if (!output_file) throw std::runtime_error("cannot open switch diagnostic output");
        output_file << report.str();
    }
    return raw_prepared_equal && balanced_recomposition_exact &&
        std::isfinite(static_cast<double>(energy + balanced_energy)) ? 0 : 1;
}

// Constant imaginary slots expose the range needed after a functional square.
// No keys/noise are needed: encoding itself must represent the desired value.
int retainedRangeDiagnostic(const std::string &output_path)
{
    using namespace TFHEpp;
    using P = BenchmarkParam<9, 7>;
    CKKSSlotVector<P> slots, decoded;
    slots.fill({0, 32});
    Polynomial<P> encoded;
    ckksSlotEncode<P, 51, 49>(encoded, slots);
    ckksSlotDecode<P, 51, 49>(decoded, encoded);
    double old_error = 0, new_error = 0;
    for (std::size_t i = 0; i < slots.size(); i++)
        old_error = finiteMaximum(old_error, std::abs(decoded[i] - slots[i]));
    ckksSlotEncode<P, 59, 49>(encoded, slots);
    ckksSlotDecode<P, 59, 49>(decoded, encoded);
    for (std::size_t i = 0; i < slots.size(); i++)
        new_error = finiteMaximum(new_error, std::abs(decoded[i] - slots[i]));
    std::ofstream file;
    if (!output_path.empty()) {
        file.open(output_path);
        if (!file) throw std::runtime_error("cannot open range diagnostic output");
    }
    auto &report = file.is_open() ? static_cast<std::ostream &>(file) : std::cout;
    report.precision(17);
    report << "{\n  \"old_square_log_q\": 51, \"new_square_log_q\": 59, \"log_delta\": 49,\n"
           << "  \"constant_imaginary_value\": 32, \"old_max_error\": " << old_error
           << ", \"new_max_error\": " << new_error << "\n}\n";
    return old_error > 1 && new_error < 1e-9 ? 0 : 1;
}

// Reproduce offline accuracy/coefficient-growth evidence without encrypted keys.
int approximationSweep(const std::string &output_path)
{
    using namespace TFHEpp;
    std::ofstream file;
    if (!output_path.empty()) {
        file.open(output_path);
        if (!file) throw std::runtime_error("cannot open approximation sweep output");
    }
    auto &report = file.is_open() ? static_cast<std::ostream &>(file) : std::cout;
    report.precision(17);
    report << "{\n  \"message_bound\": 4, \"gap\": 16, \"grid_intervals\": 65536,\n"
           << "  \"bound_method\": \"residual-curvature\",\n"
           << "  \"results\": [\n";
    bool first = true;
    for (const std::string name : {"tanh", "sigmoid", "gelu"}) {
        const auto target = [&](double x) {
            if (name == "sigmoid") return 1 / (1 + std::exp(-x));
            if (name == "gelu") return 0.5 * x * (1 + std::erf(x / std::sqrt(2.0)));
            return std::tanh(x);
        };
        const auto target_curvature = [&](double x) {
            if (name == "sigmoid") {
                const double y = target(x);
                return y * (1 - y) * (1 - 2 * y);
            }
            if (name == "gelu")
                return (2 - x * x) * std::exp(-x * x / 2) / std::sqrt(2 * std::acos(-1.0));
            const double y = std::tanh(x);
            return -2 * y * (1 - y * y);
        };
        for (std::size_t degree : {6, 8, 10, 12, 14, 16, 18, 20}) {
            const auto approximation = CKKSBuildNarrowTrigonometricApproximation(
                target, 4, 16, degree, 2.0, 65536, 40, 2.0, target_curvature, 160.0);
            double norm = -std::abs(approximation.coefficients[0]);
            for (auto value : approximation.coefficients) norm += 2 * std::abs(value);
            if (!first) report << ",\n";
            first = false;
            report << "    {\"function\": \"" << name << "\", \"degree\": " << degree
                   << ", \"sampled_error\": " << approximation.max_sampled_error
                   << ", \"uniform_bound\": " << *approximation.uniform_error_bound
                   << ", \"coefficient_l1_norm\": " << norm
                   << ", \"complex_23_bit_bound_passed\": "
                   << (std::sqrt(2.0) * *approximation.uniform_error_bound < std::ldexp(1.0, -23) ? "true" : "false")
                   << "}";
        }
    }
    report << "\n  ]\n}\n";
    return 0;
}

// Isolate the LL14 S2C path without allocating selectors or product-tree keys.
// It uses the same selected profile and public transform as the full refresh.
int s2cDiagnostic()
{
    using namespace TFHEpp;
    using P = BenchmarkParam<14, 7>;
    using Narrow = CKKSNarrowBootstrapSchedule<P, 48, 379, 49, 34, 31>;
    using Refresh = CKKSNarrowRefreshSchedule<Narrow, 41, 37>;
    auto dense = std::make_unique<Key<P>>();
    auto sparse = std::make_unique<Key<P>>();
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &coefficient : *dense) coefficient = static_cast<P::T>(ternary(generator));
    CKKSNarrowStructuredKeyGen<Narrow>(*sparse);
    CKKSSlotVector<P> slots;
    std::mt19937_64 engine(20262119);
    std::uniform_real_distribution<double> message(-0.125, 0.125);
    for (auto &slot : slots) slot = {message(engine), message(engine)};
    auto input = std::make_unique<typename Refresh::Input>();
    ckksSlotEncrypt<P>(*input, slots, *dense);
    std::cout << "LL14 input encrypted; building S2C plan" << std::endl;
    auto plan = std::make_unique<typename Refresh::S2CPlan>();
    CKKSBuildPackedSlotToCoeffPlan<P, Refresh::input_log_q, Refresh::input_log_delta,
        Refresh::s2c_plain_log_delta>(*plan, Refresh::bsgs_step);
    auto input_galois = std::make_unique<CKKSGaloisKey<P, Refresh::input_log_q>>();
    auto output_galois = std::make_unique<CKKSGaloisKey<P, Refresh::coefficient_log_q>>();
    auto to_sparse = std::make_unique<CKKSSecretKeySwitchKey<P, Refresh::coefficient_log_q>>();
    CKKSGaloisKeyGen<P, Refresh::input_log_q>(*input_galois, *dense);
    CKKSGaloisKeyGen<P, Refresh::coefficient_log_q>(*output_galois, *dense);
    CKKSSecretKeySwitchKeyGen<P, Refresh::coefficient_log_q>(*to_sparse, *dense, *sparse);
    std::cout << "LL14 S2C keys ready; evaluating transform" << std::endl;
    auto coefficients = std::make_unique<typename Refresh::Coefficients>();
    auto switched = std::make_unique<typename Refresh::Coefficients>();
    CKKSRealLinearTransform<P, Refresh::input_log_q, Refresh::input_log_delta,
        Refresh::s2c_plain_log_delta>(*coefficients, *input, *plan, *input_galois, *output_galois);
    std::cout << "LL14 S2C transform done; switching to structured key" << std::endl;
    CKKSSecretKeySwitch<P>(*switched, *coefficients, *to_sparse);
    std::array<double, P::n> decoded;
    ckksDecrypt<P>(decoded, *switched, *sparse);
    double maximum_error = 0;
    for (std::size_t i = 0; i < P::n / 2; i++) {
        maximum_error = finiteMaximum(maximum_error, std::abs(decoded[i] - slots[i].real()));
        maximum_error = finiteMaximum(maximum_error, std::abs(decoded[i + P::n / 2] - slots[i].imag()));
    }
    std::cout << "LL14 S2C max error=" << maximum_error << std::endl;
    return maximum_error < std::ldexp(1.0, -16) ? 0 : 1;
}

int factorDiagnostic()
{
    using namespace TFHEpp;
    using P = BenchmarkParam<14, 7>;
    using Narrow = CKKSNarrowBootstrapSchedule<P, 48, 379, 49, 34, 31>;
    auto dense = std::make_unique<Key<P>>();
    auto sparse = std::make_unique<Key<P>>();
    std::uniform_int_distribution<int> ternary(-1, 1);
    for (auto &coefficient : *dense) coefficient = static_cast<P::T>(ternary(generator));
    CKKSNarrowStructuredKeyGen<Narrow>(*sparse);
    const auto positions = narrow_detail::structuredPositions<Narrow>(*sparse);
    const bool negative = (*sparse)[positions[0]] == -typename P::T{1};
    std::cout << "LL14 factor: generating one selector block and conjugation keys" << std::endl;
    auto key = std::make_unique<CKKSNarrowBootstrapKey<Narrow>>();
    key->cache_mode = CKKSNarrowCacheMode::Streamed;
    key->head.resize(Narrow::block_width);
#pragma omp parallel for schedule(static)
    for (std::size_t y = 0; y < Narrow::block_width; y++) {
        CKKSSlotVector<P> mask{};
        if (y == static_cast<std::size_t>(positions[0]))
            std::fill_n(mask.begin(), Narrow::head_slots,
                negative ? std::complex<double>{0, 1} : std::complex<double>{1, 0});
        ckksSlotEncrypt<P>(key->head[y], mask, *dense);
    }
    CKKSGaloisKeyGen<P, Narrow::factor_log_q>(key->factor_galois, *dense);
    auto a = std::make_unique<std::array<double, P::n>>();
    std::mt19937_64 engine(20262119);
    std::uniform_real_distribution<double> mask(-0.5, 0.5);
    for (auto &value : *a) value = mask(engine);
    const auto table = narrow_detail::headTable<Narrow>(*a, 0, 1);
    std::cout << "LL14 factor: evaluating streamed NTT matrix" << std::endl;
    const auto factor = std::make_unique<typename Narrow::Factor>(
        narrow_detail::packedFactor<Narrow>(table, key->head, key->head_ntt, 0, *key));
    std::cout << "LL14 factor: decrypting result" << std::endl;
    CKKSSlotVector<P> got, expected;
    ckksSlotDecrypt<P>(got, *factor, *dense);
    rotateCKKSSlotVector<P>(expected, table, -positions[0]);
    double maximum_error = 0;
    for (std::size_t i = 0; i < Narrow::slots; i++) {
        const auto wanted = i >= Narrow::head_slots ? std::complex<double>{} :
            negative ? std::conj(expected[i]) : expected[i];
        maximum_error = finiteMaximum(maximum_error, std::abs(got[i] - wanted));
    }
    std::cout << "LL14 first-factor max error=" << maximum_error << std::endl;
    return maximum_error < std::ldexp(1.0, -16) ? 0 : 1;
}

}  // namespace

int main(int argc, char **argv)
{
    if (argc >= 3 && argc <= 4) {
        const std::string mode = argv[1];
        const bool small = mode.starts_with("--small-high-packed-bounded-");
        const std::string prefix = small ? "--small-high-packed-bounded-" : "--ll14-packed-bounded-";
        if (mode == prefix + "checkpoint" || mode == prefix + "replay") {
            const bool replay = mode == prefix + "replay";
            const auto selectors = (std::filesystem::path(argv[2]) / "selectors").string();
            if (small) return benchmark<9, true, false, false, true, false, false, 0, true, true>(
                argc == 4 ? argv[3] : "", false, "tanh", false, false, false,
                true, true, true, true, 2, selectors, true, 40, true, false, argv[2], replay);
            return benchmark<14, true, false, false, true, false, false, 0, true, true>(
                argc == 4 ? argv[3] : "", false, "tanh", false, false, false,
                true, true, true, true, 2, selectors, true, 40, true, false, argv[2], replay);
        }
    }
    if (argc == 5 && std::string(argv[1]) == "--small-ll14-plain37-aux16-fbs-failure-capture-test")
        return benchmark<9, true, true, true, true, true, true, 0, true, true, false, true, false, 16>(
            argv[4], false, argv[2], true, false, std::string(argv[2]) == "tanh",
            false, false, false, false, 0, argv[3], false, 7, true, false, "", false, true);
    if (argc == 5 && std::string(argv[1]) == "--small-ll14-plain37-fbs-failure-capture-test")
        return benchmark<9, true, true, true, true, true, true, 0, true, true, false, true>(
            argv[4], false, argv[2], true, false, std::string(argv[2]) == "tanh",
            false, false, false, false, 0, argv[3], false, 7, true, false, "", false, true);
    if (argc >= 4 && argc <= 5) {
        const std::string mode = argv[1];
        const bool small = mode.starts_with("--small-");
        const bool aux16 = mode.starts_with(small ? "--small-ll14-plain37-aux16-fbs-" : "--ll14-plain37-aux16-fbs-");
        const auto prefix = aux16 ? (small ? "--small-ll14-plain37-aux16-fbs-" : "--ll14-plain37-aux16-fbs-") :
                                   (small ? "--small-ll14-plain37-fbs-" : "--ll14-plain37-fbs-");
        if (mode == std::string(prefix) + "checkpoint" || mode == std::string(prefix) + "replay" ||
            mode == std::string(prefix) + "replay-trace") {
            const bool replay = mode != std::string(prefix) + "checkpoint";
            const bool trace = mode == std::string(prefix) + "replay-trace";
            const auto selectors = (std::filesystem::path(argv[3]) / "selectors").string();
            if (aux16) {
                if (small) return benchmark<9, true, true, true, true, true, true, 0, true, true, false, true, false, 16>(
                    argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh",
                    false, false, false, false, 0, selectors, false, 7, true, trace, argv[3], replay);
                return benchmark<14, true, true, true, true, true, true, 0, true, true, false, true, false, 16>(
                    argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh",
                    false, false, false, false, 0, selectors, false, 40, true, trace, argv[3], replay);
            }
            if (small) return benchmark<9, true, true, true, true, true, true, 0, true, true, false, true>(
                argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh",
                false, false, false, false, 0, selectors, false, 7, true, trace, argv[3], replay);
            return benchmark<14, true, true, true, true, true, true, 0, true, true, false, true>(
                argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh",
                false, false, false, false, 0, selectors, false, 40, true, trace, argv[3], replay);
        }
    }
    if (argc == 2 && std::string(argv[1]) == "--validation-finiteness")
        return validationFinitenessDiagnostic();

    if (argc == 5) {
        const std::string mode = argv[1];
        for (const bool small : {false, true})
            for (const bool lean : {false, true}) {
                const std::string prefix = std::string(small ? "--small-" : "--") +
                    (lean ? "lean-basis-fbs-" : "basis-fbs-");
                if (mode != prefix + "checkpoint" && mode != prefix + "replay") continue;
                const bool replay = mode == prefix + "replay";
                const auto selectors = (std::filesystem::path(argv[3]) / "selectors").string();
                if (small && lean) return basisBenchmark<9, true>(argv[2], selectors, argv[4], false, argv[3], replay);
                if (small) return basisBenchmark<9>(argv[2], selectors, argv[4], false, argv[3], replay);
                if (lean) return basisBenchmark<15, true>(argv[2], selectors, argv[4], false, argv[3], replay);
                return basisBenchmark<15>(argv[2], selectors, argv[4], false, argv[3], replay);
            }
    }

    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--basis-storage")
        return basisBenchmark<15>("", "", argc == 3 ? argv[2] : "", true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--lean-basis-storage")
        return basisBenchmark<15, true>("", "", argc == 3 ? argv[2] : "", true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-lean-basis-fbs-file")
        return basisBenchmark<9, true>(argv[2], argv[3], argc == 5 ? argv[4] : "");
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--lean-basis-fbs-file")
        return basisBenchmark<15, true>(argv[2], argv[3], argc == 5 ? argv[4] : "");
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-basis-fbs-file")
        return basisBenchmark<9>(argv[2], argv[3], argc == 5 ? argv[4] : "");
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--basis-fbs-file")
        return basisBenchmark<15>(argv[2], argv[3], argc == 5 ? argv[4] : "");
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-ll13-fbs-file")
        return benchmark<9, false, true, true, true, true, true, 0, true, true, false, false, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 7, true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--ll13-fbs-file")
        return benchmark<13, false, true, true, true, true, true, 0, true, true, false, false, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 40, true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-ll14-plain37-fbs-file")
        return benchmark<9, true, true, true, true, true, true, 0, true, true, false, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 7, true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--ll14-plain37-fbs-file")
        return benchmark<14, true, true, true, true, true, true, 0, true, true, false, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 40, true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-ll14-plain37-aux16-fbs-file")
        return benchmark<9, true, true, true, true, true, true, 0, true, true, false, true, false, 16>(argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 7, true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--ll14-plain37-aux16-fbs-file")
        return benchmark<14, true, true, true, true, true, true, 0, true, true, false, true, false, 16>(argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 40, true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-ll14-plain37-fbs-trace-file")
        return benchmark<9, true, true, true, true, true, true, 0, true, true, false, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 7, true, true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--ll14-plain37-fbs-trace-file")
        return benchmark<14, true, true, true, true, true, true, 0, true, true, false, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false, std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 40, true, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-ht14-lean-staged-file")
        return benchmark<9, false, false, false, true, true, true, 14, true, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, false, false, false, false, 0, argv[2], false, 7, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-ht15-lean-staged-file")
        return benchmark<9, true, false, false, true, true, true, 15, true, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, false, false, false, false, 0, argv[2], false, 7, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--ht14-lean-staged-file")
        return benchmark<14, false, false, false, true, true, true, 14, true, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, false, false, false, false, 0, argv[2], false, 40, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--ht15-lean-staged-file")
        return benchmark<15, true, false, false, true, true, true, 15, true, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, false, false, false, false, 0, argv[2], false, 40, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-ht15-staged-file")
        return benchmark<9, true, false, false, true, true, true, 15, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, false, false, false, false, 0, argv[2], false, 7, true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-high-staged-fbs-file")
        return benchmark<9, true, true, true, true, true, true, 0, true, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false,
            std::string(argv[2]) == "tanh", false, false, false, false, 0, argv[3], false, 7, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--ht15-staged-file")
        return benchmark<15, true, false, false, true, true, true, 15, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, false, false, false, false, 0, argv[2], false, 40, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-ht15-chunked-file")
        return benchmark<9, true, false, false, true, true, true, 15, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, true, true, false, false, 0, argv[2], true, 7);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-high-chunked-fbs-file")
        return benchmark<9, true, true, true, true, true, true, 0, true, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false,
            std::string(argv[2]) == "tanh", true, true, false, false, 0, argv[3], true, 7);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--ll14-bounded-fbs-file")
        return benchmark<14, true, true, true, true, true, true, 0, true, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false,
            std::string(argv[2]) == "tanh", true, true, false, false, 0, argv[3], true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-ht15-bounded-file")
        return benchmark<9, true, false, false, true, true, true, 15, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, true, true, false, false, 0, argv[2], true);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-high-bounded-fbs-file")
        return benchmark<9, true, true, true, true, true, true, 0, true, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false,
            std::string(argv[2]) == "tanh", true, true, false, false, 0, argv[3], true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-ht15-compact-all-file")
        return benchmark<9, true, false, false, true, true, true, 15, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, true, true, true, false, 0, argv[2]);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-high-compact-all-fbs-file")
        return benchmark<9, true, true, true, true, true, true, 0, true, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false,
            std::string(argv[2]) == "tanh", true, true, true, false, 0, argv[3]);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-ht15-compact-file")
        return benchmark<9, true, false, false, true, true, true, 15, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, true, true, true, false, 0, argv[2]);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--small-high-compact-fbs-file")
        return benchmark<9, true, true, true, true, true, true, 0, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false,
            std::string(argv[2]) == "tanh", true, true, true, false, 0, argv[3]);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--storage")
        return storageDiagnostic(argc == 3 ? argv[2] : "");
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-auxiliary-file")
        return benchmark<9, false, false, false, true, true, true>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, true, true, true, false, 0, argv[2]);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-ht15-file")
        return benchmark<9, true, false, false, true, true, true, 15>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, true, true, true, false, 0, argv[2]);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--ht14-file")
        return benchmark<14, false, false, false, true, true, true, 14>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, true, true, true, false, 0, argv[2]);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--ht15-file")
        return benchmark<15, true, false, false, true, true, true, 15>(argc == 4 ? argv[3] : "", false, "tanh", false, false, false, false, false, false, false, 0, argv[2]);
    if (argc >= 4 && argc <= 5 && std::string(argv[1]) == "--ll14-auxiliary-fbs-file")
        return benchmark<14, true, true, true, true, true, true>(argc == 5 ? argv[4] : "", false, argv[2], true, false,
            std::string(argv[2]) == "tanh", true, true, true, false, 0, argv[3]);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-ht14")
        return benchmark<9, false, false, false, true, true, true, 14>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-ht15")
        return benchmark<9, true, false, false, true, true, true, 15>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ht14")
        return benchmark<14, false, false, false, true, true, true, 14>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ht15")
        return benchmark<15, true, false, false, true, true, true, 15>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-auxiliary")
        return benchmark<9, false, false, false, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll13-auxiliary")
        return benchmark<13, false, false, false, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 2 && argc <= 4 && (std::string(argv[1]) == "--small-high-auxiliary-fbs" ||
                                 std::string(argv[1]) == "--ll14-auxiliary-fbs")) {
        const std::string function = argc >= 3 ? argv[2] : "sigmoid";
        const std::string report = argc == 4 ? argv[3] : "";
        if (std::string(argv[1]) == "--small-high-auxiliary-fbs")
            return benchmark<9, true, true, true, true, true, true>(report, false, function, true, false,
                function == "tanh", true, true, true);
        return benchmark<14, true, true, true, true, true, true>(report, false, function, true, false,
            function == "tanh", true, true, true);
    }
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll13-auxiliary-switch")
        return auxiliarySwitchDiagnostic(argc == 3 ? argv[2] : "");
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-late-rescale")
        return benchmark<9, true, false, false, true, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-sigmoid-precise-late-rescale")
        return benchmark<9, true, true, true, true, true>(argc == 3 ? argv[2] : "", false, "sigmoid", true, false, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll13-multiplication-order")
        return multiplicationOrderDiagnostic(argc == 3 ? argv[2] : "");
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll13-switch-noise")
        return switchNoiseDiagnostic(argc == 3 ? argv[2] : "");
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-streamed-fft-all")
        return benchmark<9, false, false, false, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll13-streamed-fft-all")
        return benchmark<13, false, false, false, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--small-high-packed-bounded-file")
        return benchmark<9, true, false, false, true, false, false, 0, true, true>(
            argc == 4 ? argv[3] : "", false, "tanh", false, false, false,
            true, true, true, true, 2, argv[2], true, 40, true);
    if (argc >= 3 && argc <= 4 && std::string(argv[1]) == "--ll14-packed-bounded-file")
        return benchmark<14, true, false, false, true, false, false, 0, true, true>(
            argc == 4 ? argv[3] : "", false, "tanh", false, false, false,
            true, true, true, true, 2, argv[2], true, 40, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-packed-fft")
        return benchmark<9, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true, true, 2);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll14-packed-fft")
        return benchmark<14, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true, true, 2);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-packed")
        return benchmark<9, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll14-packed")
        return benchmark<14, true>(argc == 3 ? argv[2] : "", false, "tanh", false, false, false, true, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-tanh-windowed-fft-all")
        return benchmark<9, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", true, false, true, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll14-tanh-windowed-fft-all")
        return benchmark<14, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", true, false, true, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-tanh-windowed-fft-dd")
        return benchmark<9, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", true, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll14-tanh-windowed-fft-dd")
        return benchmark<14, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", true, false, true, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-tanh-windowed-fft")
        return benchmark<9, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", true, false, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll14-tanh-windowed-fft")
        return benchmark<14, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", true, false, true, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-tanh-windowed")
        return benchmark<9, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", true, false, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll14-tanh-windowed")
        return benchmark<14, true, true, true>(argc == 3 ? argv[2] : "", false, "tanh", true, false, true);
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--retained-range")
        return retainedRangeDiagnostic(argc == 3 ? argv[2] : "");
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--approximations")
        return approximationSweep(argc == 3 ? argv[2] : "");
    if (argc >= 2 && argc <= 3) {
        const std::string mode = argv[1];
        for (const std::string function : {"tanh", "sigmoid", "gelu"}) {
            if (mode == "--small-high-" + function + "-constant")
                return benchmark<9, true, true, true>(argc == 3 ? argv[2] : "", false, function, true, true);
            if (mode == "--ll14-" + function + "-constant")
                return benchmark<14, true, true, true>(argc == 3 ? argv[2] : "", false, function, true, true);
            if (mode == "--small-high-" + function + "-precise")
                return benchmark<9, true, true, true>(argc == 3 ? argv[2] : "", false, function, true);
            if (mode == "--ll14-" + function + "-precise")
                return benchmark<14, true, true, true>(argc == 3 ? argv[2] : "", false, function, true);
            if (mode == "--small-high-" + function + "-full")
                return benchmark<9, true, true, true>(argc == 3 ? argv[2] : "", false, function);
            if (mode == "--ll14-" + function + "-full")
                return benchmark<14, true, true, true>(argc == 3 ? argv[2] : "", false, function);
        }
    }
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--small-high-tanh")
        return benchmark<9, true, true>(argc == 3 ? argv[2] : "");
    if (argc >= 2 && argc <= 3 && std::string(argv[1]) == "--ll14-tanh")
        return benchmark<14, true, true>(argc == 3 ? argv[2] : "");
    if (argc == 2 && std::string(argv[1]) == "--ll14-s2c") return s2cDiagnostic();
    if (argc == 2 && std::string(argv[1]) == "--ll14-factor") return factorDiagnostic();
    if (argc < 2 || (std::string(argv[1]) != "--small" && std::string(argv[1]) != "--ll13" &&
                    std::string(argv[1]) != "--small-zero-noise" &&
                    std::string(argv[1]) != "--small-high" && std::string(argv[1]) != "--ll14") || argc > 3) {
        std::cerr << "Usage: ckks_narrow_benchmark --small|--ll13|--small-zero-noise|--small-high|--ll14 [report.json]\n"
                  << "       ckks_narrow_benchmark --small-streamed-fft-all|--ll13-streamed-fft-all [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-packed-bounded-file|--ll14-packed-bounded-file selector-dir [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-packed-bounded-checkpoint|--ll14-packed-bounded-checkpoint new-checkpoint-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-packed-bounded-replay|--ll14-packed-bounded-replay checkpoint-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-packed|--ll14-packed [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-packed-fft|--ll14-packed-fft [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-tanh|--ll14-tanh [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-{tanh,sigmoid,gelu}-full|--ll14-{tanh,sigmoid,gelu}-full [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-{tanh,sigmoid,gelu}-precise|--ll14-{tanh,sigmoid,gelu}-precise [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-{tanh,sigmoid,gelu}-constant|--ll14-{tanh,sigmoid,gelu}-constant [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-tanh-windowed|--ll14-tanh-windowed [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-tanh-windowed-fft|--ll14-tanh-windowed-fft [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-tanh-windowed-fft-dd|--ll14-tanh-windowed-fft-dd [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-tanh-windowed-fft-all|--ll14-tanh-windowed-fft-all [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-late-rescale|--small-high-sigmoid-precise-late-rescale [report.json]\n"
                  << "       ckks_narrow_benchmark --small-auxiliary|--ll13-auxiliary [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ht14|--small-ht15|--ht14|--ht15 [report.json]\n"
                  << "       ckks_narrow_benchmark --small-auxiliary-file|--small-ht15-file|--ht14-file|--ht15-file cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --ll14-auxiliary-fbs-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ht15-compact-file cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ht15-compact-all-file cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ht15-bounded-file cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ht15-staged-file cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-staged-fbs-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --ht15-staged-file cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ht14-lean-staged-file|--small-ht15-lean-staged-file cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --ht14-lean-staged-file|--ht15-lean-staged-file cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ll13-fbs-file|--ll13-fbs-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ll14-plain37-fbs-file|--ll14-plain37-fbs-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ll14-plain37-fbs-trace-file|--ll14-plain37-fbs-trace-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ll14-plain37-fbs-checkpoint|--ll14-plain37-fbs-checkpoint sigmoid|gelu|tanh new-checkpoint-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-ll14-plain37-fbs-replay[-trace]|--ll14-plain37-fbs-replay[-trace] sigmoid|gelu|tanh checkpoint-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-bounded-fbs-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --ll14-bounded-fbs-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-compact-all-fbs-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-compact-fbs-file sigmoid|gelu|tanh cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-high-auxiliary-fbs|--ll14-auxiliary-fbs [sigmoid|gelu|tanh] [report.json]\n"
                  << "       ckks_narrow_benchmark --ll13-auxiliary-switch [report.json]\n"
                  << "       ckks_narrow_benchmark --ll13-multiplication-order [report.json]\n"
                  << "       ckks_narrow_benchmark --ll13-switch-noise [report.json]\n"
                  << "       ckks_narrow_benchmark --retained-range [report.json]\n"
                  << "       ckks_narrow_benchmark --approximations [report.json]\n"
                  << "       ckks_narrow_benchmark --ll14-s2c|--ll14-factor\n"
                  << "       ckks_narrow_benchmark --storage [report.json]\n"
                  << "       ckks_narrow_benchmark --basis-storage [report.json]\n"
                  << "       ckks_narrow_benchmark --lean-basis-storage [report.json]\n"
                  << "       ckks_narrow_benchmark --small-lean-basis-fbs-file|--lean-basis-fbs-file sigmoid|gelu|tanh|all cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --small-basis-fbs-file|--basis-fbs-file sigmoid|gelu|tanh|all cache-directory [report.json]\n"
                  << "       ckks_narrow_benchmark --[small-][lean-]basis-fbs-checkpoint sigmoid|gelu|tanh|all new-checkpoint-directory report.json\n"
                  << "       ckks_narrow_benchmark --[small-][lean-]basis-fbs-replay sigmoid|gelu|tanh|all checkpoint-directory report.json\n";
        return argc == 1 ? 0 : 2;
    }
    const std::string output = argc == 3 ? argv[2] : "";
    if (std::string(argv[1]) == "--ll14") return benchmark<14, true>(output);
    if (std::string(argv[1]) == "--small-high") return benchmark<9, true>(output);
    return std::string(argv[1]) == "--ll13" ? benchmark<13>(output) :
        benchmark<9>(output, std::string(argv[1]) == "--small-zero-noise");
}
