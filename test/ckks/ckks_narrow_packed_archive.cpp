#include <ckks/narrow_packed_file.hpp>
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
template <std::uint32_t Nbit, std::uint32_t BarBits = 4>
struct ArchiveParam {
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
    static constexpr std::uint32_t l̅ = 448 / BarBits, l̅ₐ = l̅, B̅gbit = BarBits, B̅gₐbit = BarBits;
};

template <std::uint32_t Nbit, std::uint32_t BarBits = 4, bool HighPrecision = false,
          std::uint32_t MaskBits = 34>
int run(const std::string &mode, const std::filesystem::path &input,
        const std::filesystem::path &output, std::size_t depths,
        const std::filesystem::path &scratch)
{
    using P = ArchiveParam<Nbit, BarBits>;
    using Schedule = TFHEpp::CKKSNarrowBootstrapSchedule<P, 48, HighPrecision ? 387 : 379,
        HighPrecision ? 57 : 49, 34, 31>;
    using Reconstruction = std::conditional_t<MaskBits == 34, Schedule,
        TFHEpp::CKKSNarrowPackedReconstructionSchedule<Schedule, MaskBits>>;
    if (depths > P::nbit) throw std::invalid_argument("cache depth exceeds ring bit count");
    auto key = TFHEpp::CKKSNarrowLoadPackedKey<Schedule>(input);
    for (const auto &cache : key->substitutions_fft)
        if (cache) throw std::runtime_error("loaded archive contains a derived cache");
    double preparation = 0, reconstruction = 0;
    if (mode == "--copy") TFHEpp::CKKSNarrowSavePackedKey<Schedule>(*key, output);
    if (mode == "--reconstruct" || mode == "--reconstruct-mask42") {
        std::cerr << "Preparing " << depths << " substitution depths ("
            << depths * TFHEpp::CKKSNarrowPackedSubstitutionCacheBytesPerDepth<Schedule>()
            << " FFT spectrum bytes)\n";
        auto start = std::chrono::steady_clock::now();
        TFHEpp::CKKSNarrowPreparePackedSubstitutions<Schedule>(*key, depths);
        preparation = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
        auto selectors = TFHEpp::CKKSNarrowReconstructSelectorFile<Reconstruction>(*key, output,
            [&](std::size_t completed) {
                if (completed % 256 == 0 || completed == 4 * P::n)
                    std::cerr << "Reconstructed terms " << completed << "/" << 4 * P::n
                        << " elapsed_seconds=" << std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - start).count() << std::endl;
            }, scratch);
        reconstruction = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        selectors->preserve();
    }
    std::cout << std::setprecision(17)
        << "{\"n\":" << P::n << ",\"fresh_keys_generated\":false,\"archive_validated\":true"
        << ",\"packed_log_q\":" << TFHEpp::CKKSNarrowPackedKey<Schedule>::log_q
        << ",\"selector_log_delta\":" << Schedule::log_delta
        << ",\"selector_log_q\":" << Reconstruction::selector_log_q
        << ",\"mask_log_delta\":" << Reconstruction::plain_log_delta
        << ",\"primary_bits\":" << P::B̅gbit
        << ",\"reconstructed\":" << ((mode == "--reconstruct" || mode == "--reconstruct-mask42") ? "true" : "false")
        << ",\"cached_depths\":" << ((mode == "--reconstruct" || mode == "--reconstruct-mask42") ? depths : 0)
        << ",\"prepare_seconds\":" << preparation
        << ",\"reconstruction_seconds\":" << reconstruction << "}\n";
    return 0;
}
} // namespace
int main(int argc, char **argv)
{
    try {
        if (argc < 3 || argc > 6) throw std::runtime_error("Usage: ckks_narrow_packed_archive --validate ARCHIVE | --copy ARCHIVE NEW_ARCHIVE | --reconstruct[\u002dmask42] ARCHIVE NEW_SELECTOR_DIRECTORY [CACHED_DEPTHS [NEW_SCRATCH_DIRECTORY]]");
        const std::string mode = argv[1];
        if (!((mode == "--validate" && argc == 3) || (mode == "--copy" && argc == 4) ||
              ((mode == "--reconstruct" || mode == "--reconstruct-mask42") && (argc == 4 || argc == 5 || argc == 6))))
            throw std::runtime_error("invalid mode or arguments");
        const auto depths = argc >= 5 ? std::stoul(argv[4]) : 2;
        const std::filesystem::path scratch = argc == 6 ? argv[5] : "";
        const std::filesystem::path output = argc >= 4 ? argv[3] : "";
        if (!output.empty() && std::filesystem::exists(output)) throw std::runtime_error("destination already exists");
        std::ifstream input(argv[2], std::ios::binary);
        if (!input) throw std::runtime_error("cannot open packed archive");
        std::array<std::uint64_t, 9> metadata{};
        { cereal::PortableBinaryInputArchive archive(input); archive(metadata); }
        if (metadata[1] == 16384 && metadata[2] == 387 && metadata[3] == 57 &&
            metadata[4] == 34 && metadata[7] == 4 && metadata[8] == 8)
        {
            if (mode == "--reconstruct-mask42")
                return run<14, 8, true, 42>(mode, argv[2], output, depths, scratch);
            return run<14, 8, true>(mode, argv[2], output, depths, scratch);
        }
        if (mode == "--reconstruct-mask42")
            throw std::runtime_error("mask42 reconstruction requires the saved packed435/57 primary8 archive");
        if (metadata[1] == 512) return run<9>(mode, argv[2], output, depths, scratch);
        if (metadata[1] == 16384) return run<14>(mode, argv[2], output, depths, scratch);
        throw std::runtime_error("unsupported packed archive ring");
    }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
}
