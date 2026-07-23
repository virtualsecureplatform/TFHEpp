#include <chrono>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

namespace {

struct NTTTestBaseParameter {
    static constexpr std::int32_t key_value_max = 1;
    static constexpr std::int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 6;
    static constexpr std::uint32_t n = 1U << nbit;
    static constexpr std::uint32_t k = 1;
    using T = TFHEpp::MultiLimbUInt<2>;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t B̅gbit = 16;
};

using SmallGLP = TFHEpp::GLParameter<NTTTestBaseParameter, 2, 17, 32>;

__int128 randomSigned(std::mt19937_64 &rng, const std::uint32_t bits)
{
    unsigned __int128 magnitude = rng();
    if (bits > 64) magnitude |= static_cast<unsigned __int128>(rng()) << 64;
    if (bits < 128)
        magnitude &= (static_cast<unsigned __int128>(1) << (bits - 1)) - 1;
    if ((rng() & 1U) == 0) return static_cast<__int128>(magnitude);
    return -static_cast<__int128>(magnitude);
}

template <class GLP>
void fillSigned(TFHEpp::GLBasePolynomial<GLP> &polynomial, std::mt19937_64 &rng,
                const std::uint32_t bits)
{
    for (auto &coefficient : polynomial)
        coefficient = TFHEpp::gl_detail::signedI128ToTorus<typename GLP::T>(
            randomSigned(rng, bits));
}

bool checkSmallDense()
{
    std::mt19937_64 rng(0x474c4e5454ULL);
    TFHEpp::GLBasePolynomial<SmallGLP> lhs{};
    TFHEpp::GLBasePolynomial<SmallGLP> rhs{};
    TFHEpp::GLBasePolynomial<SmallGLP> expected{};
    TFHEpp::GLBasePolynomial<SmallGLP> got{};

    fillSigned<SmallGLP>(lhs, rng, 85);
    fillSigned<SmallGLP>(rhs, rng, 16);
    if (TFHEpp::gl_detail::baseMultiplyNTTPrimeCount<SmallGLP>(lhs, rhs) != 2)
        return false;
    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(expected, lhs, rhs);
    TFHEpp::gl_detail::baseMultiply<SmallGLP>(got, lhs, rhs);
    if (got != expected) return false;

    fillSigned<SmallGLP>(lhs, rng, 16);
    fillSigned<SmallGLP>(rhs, rng, 16);
    if (TFHEpp::gl_detail::baseMultiplyNTTPrimeCount<SmallGLP>(lhs, rhs) != 1)
        return false;
    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(expected, lhs, rhs);
    TFHEpp::gl_detail::baseMultiply<SmallGLP>(got, lhs, rhs);
    if (got != expected) return false;

    // Wide-by-small multiplication takes the unsigned chunked NTT path.
    for (auto &coefficient : lhs) {
        coefficient.limb[0] = rng();
        coefficient.limb[1] = rng();
    }
    for (std::size_t i = 0; i < rhs.size(); i++)
        rhs[i] = typename SmallGLP::T(i % 3 == 0 ? -1 : (i % 3 == 1 ? 0 : 1));
    if (!TFHEpp::gl_detail::baseMultiplyNTT<SmallGLP>(got, lhs, rhs))
        return false;
    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(expected, lhs, rhs);
    if (got != expected) return false;

    // Two unrestricted torus operands exceed the two-prime reconstruction
    // bound and must retain the exact coefficient-domain fallback.
    for (auto &coefficient : rhs) {
        coefficient.limb[0] = rng();
        coefficient.limb[1] = rng();
    }
    if (TFHEpp::gl_detail::baseMultiplyNTT<SmallGLP>(got, lhs, rhs))
        return false;
    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(expected, lhs, rhs);
    TFHEpp::gl_detail::baseMultiply<SmallGLP>(got, lhs, rhs);
    if (got != expected) return false;

    lhs.fill({});
    rhs.fill({});
    lhs[0] = TFHEpp::gl_detail::signedI128ToTorus<typename SmallGLP::T>(
        static_cast<__int128>(1) << 120);
    rhs[0] = typename SmallGLP::T{3};
    if (TFHEpp::gl_detail::baseMultiplyNTTPrimeCount<SmallGLP>(lhs, rhs) != 0)
        return false;
    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(expected, lhs, rhs);
    TFHEpp::gl_detail::baseMultiply<SmallGLP>(got, lhs, rhs);
    return got == expected;
}

bool checkProductionShape()
{
    using GLP = TFHEpp::GL512p17Parameter;
    std::mt19937_64 rng(0x4e353132ULL);
    auto lhs = std::make_unique<TFHEpp::GLBasePolynomial<GLP>>();
    auto rhs = std::make_unique<TFHEpp::GLBasePolynomial<GLP>>();
    auto expected = std::make_unique<TFHEpp::GLBasePolynomial<GLP>>();
    auto got = std::make_unique<TFHEpp::GLBasePolynomial<GLP>>();

    // A sparse production-shape comparison validates the I/X/W layout while
    // keeping the coefficient-domain reference test inexpensive.
    for (std::uint32_t i = 0; i < 6; i++) {
        const std::size_t lhs_index = (137 * i + 11) % lhs->size();
        const std::size_t rhs_index = (211 * i + 7) % rhs->size();
        (*lhs)[lhs_index] =
            TFHEpp::gl_detail::signedI128ToTorus<typename GLP::T>(
                randomSigned(rng, 85));
        (*rhs)[rhs_index] =
            TFHEpp::gl_detail::signedI128ToTorus<typename GLP::T>(
                randomSigned(rng, 16));
    }
    if (TFHEpp::gl_detail::baseMultiplyNTTPrimeCount<GLP>(*lhs, *rhs) != 2)
        return false;
    TFHEpp::gl_detail::baseMultiplyReference<GLP>(*expected, *lhs, *rhs);
    TFHEpp::gl_detail::baseMultiply<GLP>(*got, *lhs, *rhs);
    if (*got != *expected) return false;

    fillSigned<GLP>(*lhs, rng, 85);
    fillSigned<GLP>(*rhs, rng, 16);
    const auto start = std::chrono::steady_clock::now();
    TFHEpp::gl_detail::baseMultiply<GLP>(*got, *lhs, *rhs);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cout << "n512 p17 dense 85x16 base multiply: " << seconds << " s"
              << std::endl;

    // This is the encryption/key-generation shape: an unrestricted torus
    // mask multiplied by a sparse ternary secret.
    for (auto &coefficient : *lhs)
        for (auto &limb : coefficient.limb) limb = rng();
    rhs->fill({});
    for (std::uint32_t i = 0; i < 64; i++)
        (*rhs)[(257 * i + 19) % rhs->size()] = typename GLP::T(i & 1U ? -1 : 1);
    const auto wide_start = std::chrono::steady_clock::now();
    if (!TFHEpp::gl_detail::baseMultiplyNTT<GLP>(*got, *lhs, *rhs))
        return false;
    const double wide_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      wide_start)
            .count();
    std::cout << "n512 p17 dense torus x ternary base multiply: "
              << wide_seconds << " s" << std::endl;
    return true;
}

}  // namespace

int main()
{
    if (!checkSmallDense()) {
        std::cerr << "small GL NTT multiplication regression failed"
                  << std::endl;
        return 1;
    }
    if (!checkProductionShape()) {
        std::cerr << "n512 GL NTT multiplication regression failed"
                  << std::endl;
        return 1;
    }
    std::cout << "GL NTT multiplication regression passed" << std::endl;
    return 0;
}
