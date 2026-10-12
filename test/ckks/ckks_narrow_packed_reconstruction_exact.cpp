#include <ckks/narrow_selector_file.hpp>
#include <boost/multiprecision/cpp_int.hpp>
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
using Schedule = TFHEpp::CKKSNarrowBootstrapSchedule<P, 48, 379, 49, 34, 31>;

using Integer = boost::multiprecision::cpp_int;
Integer integer(const P::T &value)
{
    Integer result = 0;
    for (std::size_t i = P::T::limbs; i-- > 0;) { result <<= 64; result += value.limb[i]; }
    return result;
}
void check(const std::filesystem::path &directory, const std::filesystem::path &terms)
{
    constexpr std::size_t head_count = Schedule::blocks * Schedule::block_width;
    constexpr std::size_t count = Schedule::selector_count;
    std::vector<std::array<Integer, 2 * P::n>> sums(count);
    std::ifstream file(terms, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open ordered reference terms");
    cereal::PortableBinaryInputArchive archive(file);
    auto term = std::make_unique<Schedule::Selector>();
    for (std::size_t t = 0; t < 4 * P::n; ++t) {
        bool tail{}; std::uint64_t index{}; archive(tail, index, *term);
        if (index >= (tail ? count - head_count : head_count)) throw std::runtime_error("reference index out of range");
        auto &sum = sums[(tail ? head_count : 0) + index];
        for (std::size_t c = 0; c < 2; ++c)
            for (std::size_t i = 0; i < P::n; ++i) sum[c * P::n + i] += integer(term->ct[c][i]);
    }
    if (file.peek() != std::char_traits<char>::eof()) throw std::runtime_error("trailing reference data");
    // Reference sums use unbounded integers and one final modular reduction;
    // no production ciphertext addition, extraction, or FFT helper is used.
    const Integer mask = (Integer{1} << 379) - 1;
    TFHEpp::CKKSNarrowSelectorFile<Schedule> actual(directory,
        TFHEpp::CKKSNarrowSelectorFile<Schedule>::Existing{});
    std::vector<Schedule::Selector> group;
    std::size_t checked = 0;
    for (bool tail : {false, true})
        for (std::size_t x = 0; x < (tail ? Schedule::tail_packs : Schedule::blocks); ++x) {
            actual.load(tail, x, group);
            for (std::size_t y = 0; y < Schedule::block_width; ++y) {
                const auto &sum = sums[(tail ? head_count : 0) + x * Schedule::block_width + y];
                for (std::size_t c = 0; c < 2; ++c)
                    for (std::size_t i = 0; i < P::n; ++i) {
                        if (integer(group[y].ct[c][i]) != (sum[c * P::n + i] & mask))
                            throw std::runtime_error("reconstructed selector differs from independent sum");
                        ++checked;
                    }
            }
        }
    std::cout << "{\"n\":512,\"reference_terms\":2048,\"selectors\":" << count
              << ",\"coefficients_checked\":" << checked << ",\"independent_integer_sums_passed\":true}\n";
}
} // namespace
int main(int argc, char **argv)
{
    try {
        if (argc != 3) throw std::runtime_error("Usage: ckks_narrow_packed_reconstruction_exact SELECTOR_DIRECTORY ORDERED_TERMS");
        check(argv[1], argv[2]); return 0;
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 2; }
}
