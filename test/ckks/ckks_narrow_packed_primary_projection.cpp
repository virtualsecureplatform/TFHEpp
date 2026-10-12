#include <ckks/narrow_refresh.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>

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


// Diagnostic projection only: no encryption or key generation is performed.
struct WideP : P {
    static constexpr std::uint32_t l̅ = 56, l̅ₐ = l̅, B̅gbit = 8, B̅gₐbit = 8;
};
using WideSchedule = TFHEpp::CKKSNarrowBootstrapSchedule<WideP, 48, 379, 49, 34, 31>;
using WidePacked = TFHEpp::CKKSNarrowPackedKey<WideSchedule>;
using Integer = boost::multiprecision::cpp_int;
Integer integer(const P::T &x) {
    Integer out = 0;
    for (std::size_t i = P::T::limbs; i-- > 0;) { out <<= 64; out += x.limb[i]; }
    return out;
}
void run(const std::filesystem::path &root, bool projected) {
    using namespace TFHEpp;
    auto old = std::make_unique<Packed>();
    std::ifstream file(root / "packed.bin", std::ios::binary);
    if (!file) throw std::runtime_error("missing fixture");
    {
        cereal::PortableBinaryInputArchive ar(file);
        std::array<std::uint64_t, 9> actual{}; ar(actual);
        if (actual != metadata) throw std::runtime_error("wrong fixture profile");
        ar(*old);
    }
    if (file.peek() != std::char_traits<char>::eof()) throw std::runtime_error("trailing key data");
    auto wide = std::make_unique<WidePacked>();
    for (std::size_t i = 0; i < 4; ++i) wide->indicators[i].ct = old->indicators[i].ct;
    constexpr auto first = CKKSKeySwitchFirstRowForLevel<P, Packed::log_q>();
    constexpr auto wide_first = CKKSKeySwitchFirstRowForLevel<WideP, Packed::log_q>();
    constexpr auto rows = CKKSAutoKey<WideP, Packed::log_q>::key_switch_rows;
    for (std::size_t d = 0; d < P::nbit; ++d)
        for (std::size_t r = 0; r < rows; ++r) {
            const auto source = (wide_first + r + 1) * 2 - 1 - first;
            if (source >= old->substitutions[d].key_switch_rows) throw std::runtime_error("row projection out of bounds");
            wide->substitutions[d][0][r] = old->substitutions[d][0][source];
        }
    if (projected) old.reset();
    else wide.reset();
    auto start = Clock::now();
    if (projected) CKKSNarrowPreparePackedSubstitutions<WideSchedule>(*wide, 2);
    else CKKSNarrowPreparePackedSubstitutions<Schedule>(*old, 2);
    const double preparation = elapsed(start);
    Key<WideP> dense{};
    std::array<int, P::n> signed_key{};
    std::mt19937_64 random(20262119);
    for (std::size_t i = 0; i < P::n; ++i) {
        signed_key[i] = static_cast<int>(random() % 3) - 1;
        dense[i] = P::T::from_signed_i64(signed_key[i]);
    }
    std::ifstream termsfile(root / "terms.bin", std::ios::binary);
    if (!termsfile) throw std::runtime_error("missing reference");
    cereal::PortableBinaryInputArchive ar(termsfile);
    auto baseline = std::make_unique<Schedule::Selector>();
    std::array<std::complex<double>, P::n / 2> delta_slots{};
    std::size_t count = 0, oracle_count = 0;
    double max_error = 0, oracle_error = 0, validation_seconds = 0;
    bool identical = true;
    start = Clock::now();
    const auto validate = [&](bool tail, std::size_t index, const auto &term) {
        const auto validation_start = Clock::now();
        bool expected_tail{}; std::uint64_t position{}; ar(expected_tail, position, *baseline);
        if (tail != expected_tail || index != position) throw std::runtime_error("term ordering changed");
        identical = identical && term.ct == baseline->ct;
        auto difference = std::make_unique<WideSchedule::Selector>();
        for (std::size_t component = 0; component < 2; ++component)
            for (std::size_t i = 0; i < P::n; ++i)
                difference->ct[component][i] = term.ct[component][i] - baseline->ct[component][i];
        ckksSlotDecrypt<WideP, 379, 49>(delta_slots, *difference, dense);
        for (const auto value : delta_slots) {
            if (!std::isfinite(std::abs(value))) throw std::runtime_error("nonfinite decoded difference");
            max_error = std::max(max_error, std::abs(value));
        }
        // Independent direct integer phase for sixteen evenly spaced terms.
        if (count % 128 == 0) {
            std::array<Integer, P::n> phase{};
            for (std::size_t i = 0; i < P::n; ++i) phase[i] = integer(difference->ct[1][i]);
            for (std::size_t i = 0; i < P::n; ++i) {
                const Integer a = integer(difference->ct[0][i]);
                for (std::size_t j = 0; j < P::n; ++j) {
                    const auto pos = i + j;
                    phase[pos % P::n] -= (pos < P::n ? 1 : -1) * a * signed_key[j];
                }
            }
            const Integer modulus = Integer{1} << 379, mask = modulus - 1;
            std::array<double, P::n> actual{};
            ckksDecrypt<WideP, 379, 49>(actual, *difference, dense);
            for (std::size_t i = 0; i < P::n; ++i) {
                phase[i] &= mask;
                if (phase[i] >= (modulus >> 1)) phase[i] -= modulus;
                const double expected = std::ldexp(phase[i].convert_to<double>(), -49);
                oracle_error = std::max(oracle_error, std::abs(expected - actual[i]));
            }
            ++oracle_count;
        }
        ++count;
        validation_seconds += elapsed(validation_start);
    };
    if (projected) CKKSNarrowExtractSelectorTerms<WideSchedule>(*wide, validate);
    else CKKSNarrowExtractSelectorTerms<Schedule>(*old, validate);
    const double seconds = elapsed(start);
    if (termsfile.peek() != std::char_traits<char>::eof() || count != 4 * P::n || oracle_count != 16)
        throw std::runtime_error("incomplete extraction/reference");
    const bool passed = max_error < std::ldexp(1.0, -23) && oracle_error < 1e-15 && (projected || identical);
    std::cout << std::setprecision(17)
        << "{\"n\":512,\"primary_bits\":" << (projected ? 8 : 4)
        << ",\"rows\":" << (projected ? rows : CKKSAutoKey<P, Packed::log_q>::key_switch_rows)
        << ",\"fresh_keys_generated\":false,\"terms\":" << count
        << ",\"prepare_seconds\":" << preparation << ",\"extraction_and_validation_seconds\":" << seconds
        << ",\"validation_seconds\":" << validation_seconds
        << ",\"extraction_seconds_excluding_callbacks\":" << seconds - validation_seconds
        << ",\"ciphertexts_identical\":" << (identical ? "true" : "false")
        << ",\"max_slot_difference\":" << max_error << ",\"direct_phase_terms\":" << oracle_count
        << ",\"direct_phase_max_difference\":" << oracle_error
        << ",\"accuracy_passed\":" << (passed ? "true" : "false") << "}\n";
    if (!passed) throw std::runtime_error("projected primary decomposition accuracy failed");
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 2 && argc != 3) throw std::runtime_error("Usage: ckks_narrow_packed_primary_projection FIXTURE [--baseline]");
        if (argc == 3 && std::string(argv[2]) != "--baseline") throw std::runtime_error("unknown mode");
        run(argv[1], argc == 2); return 0;
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 2; }
}
