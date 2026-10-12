#include <ckks/narrow_refresh.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>

namespace {
struct P {
    static constexpr std::int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = 9, n = 1 << nbit, k = 1;
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
struct FullP : P {
    static constexpr std::uint32_t nbit = 14, n = 1 << nbit;
};
using Schedule = TFHEpp::CKKSNarrowBootstrapSchedule<P, 48, 379, 49, 34, 31>;
using Packed = TFHEpp::CKKSNarrowPackedKey<Schedule>;
constexpr std::array<std::uint64_t, 9> metadata{
    0x5041434b43463031ULL, P::n, 379, 49, 34, 31, 448, 4, 4};
static_assert(Packed::log_q == 422);
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

void create(const std::filesystem::path &root)
{
    using namespace TFHEpp;
    if (!std::filesystem::create_directory(root)) throw std::runtime_error("fixture directory already exists");
    std::filesystem::permissions(root, std::filesystem::perms::owner_all,
        std::filesystem::perm_options::replace);
    auto source = std::make_unique<Key<P>>(), dense = std::make_unique<Key<P>>();
    CKKSNarrowStructuredKeyGen<Schedule>(*source);
    std::mt19937_64 random(20262119);
    for (auto &value : *dense) value = P::T::from_signed_i64(static_cast<int>(random() % 3) - 1);
    auto packed = std::make_unique<Packed>();
    CKKSNarrowPackedKeyGen<Schedule>(*packed, *source, *dense, {P::α, 0});
    std::ofstream file(root / "packed.bin", std::ios::binary);
    { cereal::PortableBinaryOutputArchive archive(file); archive(metadata, *packed); }
    file.close();
    if (!file) throw std::runtime_error("cannot write packed key fixture");
    std::cout << "saved public packed key fixture; no secrets serialized\n";
}

void run(const std::filesystem::path &root, std::size_t levels)
{
    using namespace TFHEpp;
    if (levels > P::nbit) throw std::runtime_error("cache depth exceeds extraction depth");
    auto packed = std::make_unique<Packed>();
    std::ifstream keyfile(root / "packed.bin", std::ios::binary);
    if (!keyfile) throw std::runtime_error("cannot open packed key fixture");
    std::array<std::uint64_t, 9> actual{};
    { cereal::PortableBinaryInputArchive archive(keyfile); archive(actual);
      if (actual != metadata) throw std::runtime_error("incompatible packed key fixture");
      archive(*packed); }
    if (keyfile.peek() != std::char_traits<char>::eof()) throw std::runtime_error("trailing packed key data");
    const auto prepare_start = Clock::now();
    constexpr auto per_depth = CKKSNarrowPackedSubstitutionCacheBytesPerDepth<Schedule>();
    static_assert(CKKSNarrowPackedSubstitutionCacheDepthForBudget<Schedule>(0) == 0);
    static_assert(CKKSNarrowPackedSubstitutionCacheDepthForBudget<Schedule>(per_depth - 1) == 0);
    static_assert(CKKSNarrowPackedSubstitutionCacheDepthForBudget<Schedule>(2 * per_depth - 1) == 1);
    const auto selected = CKKSNarrowPreparePackedSubstitutionsWithinBudget<Schedule>(
        *packed, levels * per_depth);
    if (selected != levels) throw std::runtime_error("cache budget selected incorrect depth");
    const auto prepare_seconds = elapsed(prepare_start);
    std::size_t cache_bytes = 0;
    for (const auto &cache : packed->substitutions_fft)
        if (cache) cache_bytes += cache->spectra.size() * sizeof(PolynomialInFD<P>);
    if (cache_bytes != levels * per_depth) throw std::runtime_error("cache exceeded spectrum budget");
    const auto terms_path = root / "terms.bin";
    const bool write = !std::filesystem::exists(terms_path);
    if (write && levels != 2) throw std::runtime_error("create the depth-two reference first");
    std::fstream file(terms_path, std::ios::binary | (write ? std::ios::out : std::ios::in));
    if (!file) throw std::runtime_error("cannot open ordered term reference");
    std::size_t terms = 0;
    const auto extract_start = Clock::now();
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
        auto expected = std::make_unique<Schedule::Selector>();
        CKKSNarrowExtractSelectorTerms<Schedule>(*packed,
            [&](bool tail, std::size_t index, const auto &term) {
                bool expected_tail{}; std::uint64_t position{};
                archive(expected_tail, position, *expected);
                if (tail != expected_tail || index != position || term.ct != expected->ct)
                    throw std::runtime_error("cache depth changed an ordered ciphertext term");
                ++terms;
            });
        if (file.peek() != std::char_traits<char>::eof()) throw std::runtime_error("trailing ordered term reference");
    }
    const auto extraction_seconds = elapsed(extract_start);
    if (terms != 4 * P::n) throw std::runtime_error("incomplete packed extraction");
    file.close();
    if (!file) throw std::runtime_error("ordered term file close failed");
    std::cout << std::setprecision(17) << "{\"n\":" << P::n << ",\"packed_log_q\":" << Packed::log_q
              << ",\"cached_depths\":" << levels << ",\"cache_bytes\":" << cache_bytes
        << ",\"cache_bytes_per_depth\":" << per_depth
              << ",\"prepare_seconds\":" << prepare_seconds << ",\"extraction_seconds\":" << extraction_seconds
              << ",\"terms\":" << terms << ",\"coefficients\":" << terms * (P::k + 1) * P::n
              << ",\"reference_written\":" << (write ? "true" : "false")
              << ",\"ordered_terms_identical\":" << (write ? "null" : "true") << "}" << std::endl;
}
} // namespace

int main(int argc, char **argv)
{
    try {
        if (argc == 2 && std::string(argv[1]) == "--full-budget-estimate") {
            using FullSchedule = TFHEpp::CKKSNarrowBootstrapSchedule<FullP, 48, 379, 49, 34, 31>;
            constexpr std::size_t budget = std::size_t{16} << 30;
            constexpr auto per_depth = TFHEpp::CKKSNarrowPackedSubstitutionCacheBytesPerDepth<FullSchedule>();
            constexpr auto levels = TFHEpp::CKKSNarrowPackedSubstitutionCacheDepthForBudget<FullSchedule>(budget);
            static_assert(levels == FullP::nbit && levels * per_depth <= budget);
            std::cout << "{\"n\":" << FullP::n << ",\"budget_bytes\":" << budget
                      << ",\"cache_bytes_per_depth\":" << per_depth
                      << ",\"cached_depths\":" << levels << ",\"cache_bytes\":" << levels * per_depth
                      << ",\"full_reconstruction_executed\":false}" << std::endl;
            return 0;
        }
        if (argc == 3 && std::string(argv[1]) == "--create") { create(argv[2]); return 0; }
        if (argc == 4 && std::string(argv[1]) == "--run") { run(argv[2], std::stoul(argv[3])); return 0; }
        std::cerr << "Usage: ckks_narrow_packed_cache_profile --create DIRECTORY | --run DIRECTORY DEPTH | --full-budget-estimate\n";
        return 2;
    }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
}
