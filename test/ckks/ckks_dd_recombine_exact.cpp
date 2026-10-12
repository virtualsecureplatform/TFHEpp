#include <ckks/ckks.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <bit>
#include <iostream>
#include <memory>
#include <random>

namespace {
using namespace TFHEpp;
using boost::multiprecision::cpp_int;

template <unsigned Nbit>
struct Param : lvl5param {
    static constexpr unsigned nbit = Nbit, n = 1U << Nbit;
    using T = MultiLimbUInt<9>;
};

template <unsigned Nbit, unsigned Q, unsigned Base>
void check()
{
    using P = Param<Nbit>;
    constexpr unsigned count = (Q + Base - 1) / Base;
    auto rows = std::make_unique<std::array<TRLWE<P>, count>>();
    auto output = std::make_unique<TRLWE<P>>();
    std::mt19937_64 rng(2026211908);
    const cpp_int mask = (cpp_int{1} << Q) - 1;
    for (unsigned sample = 0; sample < 3; ++sample) {
        for (auto &row : *rows)
            for (auto &component : row)
                for (auto &value : component)
                    for (auto &word : value.limb)
                        word = sample == 0 ? 0 : sample == 1 ? ~std::uint64_t{} : rng();
        ckks_detail::activeBaseRecombineTRLWERows<P, Q, Base, count>(*output, *rows);
        for (unsigned c = 0; c <= P::k; ++c)
            for (unsigned n = 0; n < P::n; ++n) {
                cpp_int expected = 0;
                for (unsigned row = 0; row < count; ++row) {
                    cpp_int value = 0;
                    for (unsigned word = P::T::limbs; word-- > 0;)
                        value = (value << 64) + (*rows)[row][c][n].limb[word];
                    expected += value << ((count - row - 1) * Base);
                }
                expected &= mask;
                for (const auto word : (*output)[c][n].limb) {
                    if (word != static_cast<std::uint64_t>(expected & ((cpp_int{1} << 64) - 1)))
                        throw std::runtime_error("DD recomposition differs from integer oracle");
                    expected >>= 64;
                }
            }
    }
    std::cout << "n=" << P::n << " q=" << Q << " base=" << Base
              << " oracle_coefficients=" << 3 * (P::k + 1) * P::n << '\n';
}

template <unsigned Nbit, unsigned Q, unsigned Base>
void checkNative()
{
    using P = Param<Nbit>;
    constexpr unsigned count = (Q + Base - 1) / Base;
    using Rows = std::array<std::array<std::array<std::uint64_t, P::n>, P::k + 1>, count>;
    auto rows = std::make_unique<Rows>();
    auto output = std::make_unique<TRLWE<P>>();
    std::mt19937_64 rng(2026211912);
    const cpp_int mask = (cpp_int{1} << Q) - 1;
    for (unsigned sample = 0; sample < 3; ++sample) {
        for (auto &row : *rows)
            for (auto &component : row)
                for (unsigned n = 0; n < P::n; ++n)
                    component[n] = sample == 0 ? (n % 2 ? ~std::uint64_t{} : 0)
                        : sample == 1 ? (std::uint64_t{1} << 63) - (n % 2) : rng();
        ckks_detail::activeBaseRecombineDigitRows<P, Q, Base, count>(*output, *rows);
        for (unsigned c = 0; c <= P::k; ++c)
            for (unsigned n = 0; n < P::n; ++n) {
                cpp_int expected = 0;
                for (unsigned row = 0; row < count; ++row)
                    expected += cpp_int{std::bit_cast<std::int64_t>((*rows)[row][c][n])}
                        << ((count - row - 1) * Base);
                expected &= mask;
                for (const auto word : (*output)[c][n].limb) {
                    if (word != static_cast<std::uint64_t>(expected & ((cpp_int{1} << 64) - 1)))
                        throw std::runtime_error("native DD recomposition differs from signed integer oracle");
                    expected >>= 64;
                }
            }
    }
    std::cout << "native n=" << P::n << " q=" << Q << " base=" << Base
              << " oracle_coefficients=" << 3 * (P::k + 1) * P::n << '\n';
}
}

int main() try
{
    check<9, 479, 8>();
    check<13, 479, 8>();
    check<13, 197, 11>();
    checkNative<9, 479, 8>();
    checkNative<13, 479, 8>();
    checkNative<13, 197, 11>();
}
catch (const std::exception &error)
{
    std::cerr << error.what() << '\n';
    return 1;
}
