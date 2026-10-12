#include <ckks/narrow_bootstrapping.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <random>

namespace {
using namespace TFHEpp;
using boost::multiprecision::cpp_int;

// Public arithmetic profiles: these controls allocate no encryption keys.
template <std::size_t N, std::size_t Limbs, std::size_t Q, std::size_t Scale>
struct PublicSchedule {
    struct Param {
        static constexpr std::size_t n = N, k = 1;
        using T = MultiLimbUInt<Limbs>;
    };
    static constexpr std::size_t selector_log_q = Q, plain_log_delta = Scale;
    static constexpr std::size_t log_delta = 48, blocks = 31;
    static constexpr std::size_t block_width = (N + 30) / 31;
    static constexpr std::size_t tail_packs =
        (blocks + N / (4 * (block_width - 1)) - 1) / (N / (4 * (block_width - 1)));
    struct Selector { std::array<Polynomial<Param>, 2> ct; };
};

template <class T>
cpp_int integer(const T &value)
{
    cpp_int result = 0;
    for (std::size_t i = 0; i < T::digits; i += 64)
        result += cpp_int(static_cast<std::uint64_t>(value >> i)) << i;
    return result;
}

template <class S, std::size_t Bits>
void check()
{
    using P = typename S::Param;
    using T = typename P::T;
    using Cache = CKKSNarrowNTTSelector<S, Bits>;
    static_assert(Cache::prime_count == 2);
    static_assert(!Cache::fits(0) && !Cache::fits(64));
    const cpp_int mask = (cpp_int(1) << S::selector_log_q) - 1;
    const std::int64_t bound = (std::int64_t{1} << S::plain_log_delta) + 1;
    const modular_ntt::TwoPrimeCRT crt(modular_ntt::degree65536_primes[0],
                                      modular_ntt::degree65536_primes[1]);
    const std::array<modular_ntt::NegacyclicNTTPlan, 2> plans{
        modular_ntt::NegacyclicNTTPlan(P::n, modular_ntt::degree65536_primes[0]),
        modular_ntt::NegacyclicNTTPlan(P::n, modular_ntt::degree65536_primes[1])};
    std::mt19937_64 random(20262119 + P::n);
    std::size_t checked = 0;
    for (const bool dense_boundary : {false, true}) {
        std::vector<typename S::Selector> selectors(1);
        for (std::size_t c = 0; c < 2; ++c)
            for (std::size_t i = 0; i < P::n; ++i) {
                auto &value = selectors[0].ct[c][i];
                if (dense_boundary) value = c == 0 ? T{} - T{1} : T{1};
                else {
                    for (std::size_t bit = 0; bit < T::digits; bit += 64)
                        value |= T{random()} << bit;
                    // Ensure native-width garbage above the active Q.
                    value |= T{1} << S::selector_log_q;
                }
            }
        std::vector<Cache> prepared;
        narrow_detail::prepareSelectors<S, Bits>(prepared, selectors);
        std::vector<std::int64_t> table(P::n);
        const std::array<std::size_t, 5> positions{0, 1, P::n / 2, P::n - 2, P::n - 1};
        if (dense_boundary) std::fill(table.begin(), table.end(), -bound);
        else for (std::size_t j = 0; j < positions.size(); ++j)
            table[positions[j]] = j % 2 ? bound : -bound;
        std::array<std::vector<std::uint64_t>, 2> transformed;
        for (std::size_t prime = 0; prime < 2; ++prime) {
            const auto modulus = modular_ntt::degree65536_primes[prime].value;
            transformed[prime].resize(P::n);
            for (std::size_t i = 0; i < P::n; ++i)
                transformed[prime][i] = table[i] < 0 ? modulus - (-table[i]) : table[i];
            plans[prime].forward(transformed[prime]);
        }
        for (std::size_t c = 0; c < 2; ++c) {
            auto actual = std::make_unique<Polynomial<P>>();
            for (std::size_t digit = 0; digit < Cache::digits; ++digit) {
                std::array<std::vector<std::uint64_t>, 2> product;
                for (std::size_t prime = 0; prime < 2; ++prime) {
                    const auto modulus = modular_ntt::degree65536_primes[prime].value;
                    product[prime].resize(P::n);
                    const auto row = (c * Cache::digits + digit) * P::n;
                    for (std::size_t i = 0; i < P::n; ++i) {
                        const auto multiplied = modular_ntt::multiply(
                            prepared[0].spectra[prime][row + i], transformed[prime][i], modulus);
                        // Repeat the selector contribution block_width times;
                        // dense boundary case approaches the stated CRT bound.
                        product[prime][i] = modular_ntt::multiply(multiplied, S::block_width, modulus);
                    }
                    plans[prime].inverse(product[prime]);
                }
                for (std::size_t i = 0; i < P::n; ++i)
                    (*actual)[i] += ckks_detail::signedToLevel<P, S::selector_log_q>(
                        crt.reconstructSigned(product[0][i], product[1][i])) << (digit * Bits);
            }
            for (std::size_t i = 0; i < P::n; ++i) {
                cpp_int expected = 0;
                if (dense_boundary) {
                    const cpp_int coefficient = c == 0 ? mask : cpp_int(1);
                    expected = coefficient * (-bound) *
                        (2 * static_cast<std::int64_t>(i) + 2 - static_cast<std::int64_t>(P::n));
                }
                else for (const auto position : positions) {
                    const auto source = (i + P::n - position) % P::n;
                    const cpp_int coefficient = integer(selectors[0].ct[c][source]) & mask;
                    expected += coefficient * table[position] * (i < position ? -1 : 1);
                }
                expected = (expected * S::block_width) & mask;
                if ((integer((*actual)[i]) & mask) != expected)
                    throw std::runtime_error("wide selector CRT differs from independent integer convolution");
                ++checked;
            }
        }
    }
    std::cout << "{\"n\":" << P::n << ",\"q\":" << S::selector_log_q
              << ",\"plain_scale\":" << S::plain_log_delta << ",\"digit_bits\":" << Bits
              << ",\"digits\":" << Cache::digits << ",\"coefficients_checked\":" << checked
              << ",\"dense_signed_boundary_and_sparse_wrap_exact\":true}" << std::endl;
}

template <class S, std::size_t Bits>
void retainedTiming(std::span<const typename S::Selector> selectors)
{
    using Cache = CKKSNarrowNTTSelector<S, Bits>;
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    std::vector<Cache> prepared;
    // Bound the derived cache independently of the complete retained group.
    for (std::size_t first = 0; first < selectors.size(); first += 8)
        narrow_detail::prepareSelectors<S, Bits>(prepared,
            selectors.subspan(first, std::min(std::size_t{8}, selectors.size() - first)));
    const auto seconds = std::chrono::duration<double>(Clock::now() - start).count();
    std::cout << "{\"retained_selector_count\":" << selectors.size()
              << ",\"digit_bits\":" << Bits << ",\"digits\":" << Cache::digits
              << ",\"selector_ntt_seconds\":" << seconds << "}" << std::endl;
}
}

int main(int argc, char **argv)
{
    using Aux = PublicSchedule<16384, 7, 385, 37>;
    using Packed = PublicSchedule<16384, 7, 379, 34>;
    using Basis = PublicSchedule<32768, 9, 518, 39>;
    static_assert(CKKSNarrowNTTSelector<Aux>::fits(60));
    static_assert(CKKSNarrowNTTSelector<Packed>::fits(60));
    static_assert(CKKSNarrowNTTSelector<Basis>::fits(56));
    static_assert(!CKKSNarrowNTTSelector<Basis>::fits(60));
    try {
        if (argc == 3 && std::string(argv[1]) == "--retained-aux8-timing") {
            CKKSNarrowSelectorFile<Aux> source(argv[2], CKKSNarrowSelectorFile<Aux>::Existing{});
            std::vector<Aux::Selector> selectors;
            source.load(false, 30, selectors);
            // Warm both immutable plans, then alternate complete group trials.
            retainedTiming<Aux, 48>(std::span<const Aux::Selector>(selectors).first(8));
            retainedTiming<Aux, 60>(std::span<const Aux::Selector>(selectors).first(8));
            retainedTiming<Aux, 48>(selectors);
            retainedTiming<Aux, 60>(selectors);
            retainedTiming<Aux, 60>(selectors);
            retainedTiming<Aux, 48>(selectors);
            return 0;
        }
        if (argc != 1) throw std::invalid_argument("unsupported wide-selector control arguments");
        check<Aux, 60>();
        check<Packed, 60>();
        check<Basis, 56>();
    }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
}
