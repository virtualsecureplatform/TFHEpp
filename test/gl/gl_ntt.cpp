#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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
constexpr std::uint32_t dd_log_q = 96;
constexpr std::uint32_t dd_primary_bit = 48;
constexpr std::uint32_t dd_bbar_bit = 16;
using SmallDDSwitchKey =
    TFHEpp::GLDDSmallKeySwitchKey<SmallGLP, dd_log_q, dd_primary_bit,
                                  dd_bbar_bit>;
using SmallDDBigSwitchKey =
    TFHEpp::GLDDBigKeySwitchKey<SmallGLP, dd_log_q, dd_primary_bit,
                                dd_bbar_bit>;
using SmallWTransformSchedule =
    TFHEpp::GLSHIPBootstrapSchedule<SmallGLP, 32, 4, 8, 8, 96, 18, 3, 1, 2, 2,
                                    48, 16>;
static_assert(TFHEpp::gl_detail::smallKeySwitchAccumulationNTTPrimeCount<
                  SmallGLP, SmallDDSwitchKey> == 2);

constexpr std::uint32_t one_prime_log_q = 32;
constexpr std::uint32_t one_prime_primary_bit = 16;
constexpr std::uint32_t one_prime_bbar_bit = 8;
using SmallOnePrimeDDSwitchKey =
    TFHEpp::GLDDSmallKeySwitchKey<SmallGLP, one_prime_log_q,
                                  one_prime_primary_bit, one_prime_bbar_bit>;
static_assert(TFHEpp::gl_detail::smallKeySwitchAccumulationNTTPrimeCount<
                  SmallGLP, SmallOnePrimeDDSwitchKey> == 1);

using N512DDSwitchKey =
    TFHEpp::GLDDSmallKeySwitchKey<TFHEpp::GL512p17Parameter, 338, 85, 16>;
static_assert(TFHEpp::gl_detail::smallKeySwitchAccumulationNTTPrimeCount<
                  TFHEpp::GL512p17Parameter, N512DDSwitchKey> == 2);
static_assert(TFHEpp::gl_detail::smallKeySwitchNTTKeyCacheBytes<
                  TFHEpp::GL512p17Parameter, N512DDSwitchKey> == 56623104);

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

template <class GLP>
void fillSigned(TFHEpp::GLPolynomial<GLP> &polynomial, std::mt19937_64 &rng,
                const std::uint32_t bits)
{
    for (auto &slice : polynomial) fillSigned<GLP>(slice, rng, bits);
}

bool checkDirectBaseAutomorphism()
{
    std::mt19937_64 rng(0x424153454155544fULL);
    TFHEpp::GLBasePolynomial<SmallGLP> input{};
    fillSigned<SmallGLP>(input, rng, 60);
    constexpr std::array<std::uint32_t, 4> x_multipliers{1, 3, 5, 7};
    constexpr std::array<std::uint32_t, 4> w_multipliers{1, 3, 7, 16};
    for (const std::uint32_t x_multiplier : x_multipliers)
        for (const std::uint32_t w_multiplier : w_multipliers) {
            TFHEpp::GLBasePolynomial<SmallGLP> direct{};
            TFHEpp::gl_detail::baseAutomorphism<SmallGLP>(
                direct, input, x_multiplier, w_multiplier);
            TFHEpp::GLPolynomial<SmallGLP> lifted{};
            TFHEpp::GLPolynomial<SmallGLP> transformed{};
            TFHEpp::gl_detail::liftBase<SmallGLP>(lifted, input);
            TFHEpp::gl_detail::polynomialAutomorphism<SmallGLP>(
                transformed, lifted, x_multiplier & 3U, x_multiplier, 1,
                w_multiplier);
            if (direct != transformed[0]) return false;
        }

    for (const std::uint32_t x_multiplier : x_multipliers) {
        TFHEpp::GLBasePolynomial<SmallGLP> direct{};
        TFHEpp::gl_detail::baseAutomorphism<SmallGLP>(direct, input,
                                                      x_multiplier, 1);
        const auto z_map =
            TFHEpp::gl_detail::baseXAutomorphismSpectrumMap<SmallGLP>(
                x_multiplier);
        const auto check_prime = [&](const auto &plan) {
            std::vector<std::uint64_t> source, expected,
                got(TFHEpp::gl_detail::GLBaseNTTPlan<
                    SmallGLP>::coefficient_count);
            plan.forward(source, input);
            plan.forward(expected, direct);
            TFHEpp::gl_detail::applyBaseXAutomorphismSpectrum<SmallGLP>(
                std::span<std::uint64_t>(got),
                std::span<const std::uint64_t>(source), z_map);
            return got == expected;
        };
        if (!check_prime(TFHEpp::gl_detail::baseNTTPlan<SmallGLP, 0>()) ||
            !check_prime(TFHEpp::gl_detail::baseNTTPlan<SmallGLP, 1>()))
            return false;
    }
    return true;
}

bool checkSymmetricDecompositionSpectrumHoist()
{
    using P = typename SmallGLP::baseP;
    constexpr std::uint32_t log_q = dd_log_q;
    constexpr std::uint32_t base_bit = dd_primary_bit;
    constexpr std::size_t rows = (log_q + base_bit - 1) / base_bit;
    using T = typename P::T;
    std::mt19937_64 rng(0x53594d4d45545259ULL);
    TFHEpp::GLBasePolynomial<SmallGLP> input{};
    for (auto &coefficient : input)
        coefficient = TFHEpp::ckks_detail::signedToLevel<P, log_q>(
            randomSigned(rng, log_q - 2));
    input[0] = T{1} << (base_bit - 1);
    input[1] = T{0} - (T{1} << (base_bit - 1));
    input[2] = T{1} << (log_q - 1);
    input[3] = (T{1} << (base_bit - 1)) << base_bit;

    std::array<TFHEpp::Polynomial<P>, rows> digits{};
    TFHEpp::gl_detail::activeSymmetricBaseDecomposePolynomialRows<
        P, log_q, base_bit, rows>(digits, input);
    for (std::size_t coefficient = 0; coefficient < input.size();
         coefficient++) {
        T recombined{};
        for (std::size_t row = 0; row < rows; row++)
            recombined += digits[row][coefficient]
                          << ((rows - row - 1) * base_bit);
        if (TFHEpp::ckks_detail::reduceToLevel<P, log_q>(recombined) !=
            TFHEpp::ckks_detail::reduceToLevel<P, log_q>(input[coefficient]))
            return false;
    }

    TFHEpp::GLBasePolynomial<SmallGLP> negative{};
    for (std::size_t coefficient = 0; coefficient < input.size(); coefficient++)
        negative[coefficient] = T{0} - input[coefficient];
    std::array<TFHEpp::Polynomial<P>, rows> negative_digits{};
    TFHEpp::gl_detail::activeSymmetricBaseDecomposePolynomialRows<
        P, log_q, base_bit, rows>(negative_digits, negative);
    const T minimum = T{1} << (log_q - 1);
    for (std::size_t coefficient = 0; coefficient < input.size();
         coefficient++) {
        if (TFHEpp::ckks_detail::reduceToLevel<P, log_q>(input[coefficient]) ==
            minimum)
            continue;
        for (std::size_t row = 0; row < rows; row++)
            if (negative_digits[row][coefficient] !=
                T{0} - digits[row][coefficient])
                return false;
    }

    constexpr std::array<std::uint32_t, 4> multipliers{1, 3, 5, 7};
    constexpr std::size_t z_dimension =
        TFHEpp::gl_detail::GLBaseNTTPlan<SmallGLP>::z_dimension;
    constexpr std::size_t coefficient_count =
        TFHEpp::gl_detail::GLBaseNTTPlan<SmallGLP>::coefficient_count;
    for (const std::uint32_t multiplier : multipliers) {
        TFHEpp::GLBasePolynomial<SmallGLP> automorphed{};
        TFHEpp::gl_detail::baseAutomorphism<SmallGLP>(automorphed, digits[0],
                                                      multiplier, 1);
        for (std::size_t prime = 0; prime < 2; prime++) {
            const auto &plan =
                prime == 0 ? TFHEpp::gl_detail::baseNTTPlan<SmallGLP, 0>()
                           : TFHEpp::gl_detail::baseNTTPlan<SmallGLP, 1>();
            std::vector<std::uint64_t> source_spectrum;
            std::vector<std::uint64_t> automorphed_spectrum;
            plan.forward(source_spectrum, digits[0]);
            plan.forward(automorphed_spectrum, automorphed);
            for (std::size_t w = 0; w < coefficient_count / z_dimension;
                 w++) {
                for (std::size_t z = 0; z < z_dimension; z++) {
                    const std::uint32_t mapped_root =
                        static_cast<std::uint32_t>(
                            (static_cast<std::uint64_t>(multiplier) *
                             (2 * z + 1)) %
                            (4 * SmallGLP::matrix_dimension));
                    const std::size_t source_z = (mapped_root - 1) / 2;
                    const std::size_t destination_index =
                        TFHEpp::modular_ntt::NegacyclicNTTPlan::backendIndex(
                            z_dimension, z);
                    const std::size_t source_index =
                        TFHEpp::modular_ntt::NegacyclicNTTPlan::backendIndex(
                            z_dimension, source_z);
                    if (automorphed_spectrum[w * z_dimension +
                                             destination_index] !=
                        source_spectrum[w * z_dimension + source_index])
                        return false;
                }
            }
        }
    }
    return true;
}

bool checkFusedCiphertextTensorMultiply()
{
    constexpr std::uint32_t log_q = dd_log_q;
    std::mt19937_64 rng(0x54454e534f524e54ULL);
    TFHEpp::GLBaseCiphertextData<SmallGLP> lhs{};
    TFHEpp::GLBaseCiphertextData<SmallGLP> rhs{};
    for (std::size_t component = 0; component < 2; component++) {
        fillSigned<SmallGLP>(lhs[component], rng, log_q - 1);
        fillSigned<SmallGLP>(rhs[component], rng, log_q - 1);
    }
    lhs[0][0] = TFHEpp::gl_detail::signedI128ToTorus<typename SmallGLP::T>(
        static_cast<__int128>(1) << (log_q - 2));
    lhs[1][1] = typename SmallGLP::T{0} - lhs[0][0];
    rhs[0][2] =
        (typename SmallGLP::T{1} << (log_q - 1)) - typename SmallGLP::T{1};
    rhs[1][3] = typename SmallGLP::T{1} << (log_q - 1);

    std::array<TFHEpp::GLBasePolynomial<SmallGLP>, 3> expected{};
    TFHEpp::GLBasePolynomial<SmallGLP> second_cross{};
    TFHEpp::gl_detail::baseMultiplyAtLevel<SmallGLP, log_q>(expected[0], lhs[0],
                                                            rhs[0]);
    TFHEpp::gl_detail::baseMultiplyAtLevel<SmallGLP, log_q>(expected[1], lhs[0],
                                                            rhs[1]);
    TFHEpp::gl_detail::baseMultiplyAtLevel<SmallGLP, log_q>(second_cross,
                                                            lhs[1], rhs[0]);
    TFHEpp::gl_detail::addInPlace<SmallGLP>(expected[1], second_cross);
    TFHEpp::gl_detail::reduce<SmallGLP, log_q>(expected[1]);
    TFHEpp::gl_detail::baseMultiplyAtLevel<SmallGLP, log_q>(expected[2], lhs[1],
                                                            rhs[1]);

    std::array<TFHEpp::GLBasePolynomial<SmallGLP>, 3> got{};
    TFHEpp::gl_detail::BaseCiphertextTensorNTTWorkspace<SmallGLP> workspace;
    if (!TFHEpp::gl_detail::baseCiphertextTensorMultiplyNTT<SmallGLP, log_q>(
            got, lhs, rhs, &workspace))
        return false;
    return got == expected;
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

    // Two unrestricted torus operands use the exact double-chunk NTT path.
    for (auto &coefficient : rhs) {
        coefficient.limb[0] = rng();
        coefficient.limb[1] = rng();
    }
    if (!TFHEpp::gl_detail::baseMultiplyNTT<SmallGLP>(got, lhs, rhs))
        return false;
    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(expected, lhs, rhs);
    TFHEpp::gl_detail::baseMultiply<SmallGLP>(got, lhs, rhs);
    if (got != expected) return false;

    // Level-aware double chunking may omit diagonals whose shifts are already
    // zero modulo 2^LogQ.
    TFHEpp::gl_detail::reduce<SmallGLP, dd_log_q>(lhs);
    TFHEpp::gl_detail::reduce<SmallGLP, dd_log_q>(rhs);
    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(expected, lhs, rhs);
    TFHEpp::gl_detail::reduce<SmallGLP, dd_log_q>(expected);
    TFHEpp::gl_detail::baseMultiplyAtLevel<SmallGLP, dd_log_q>(got, lhs, rhs);
    if (got != expected) return false;

    const typename SmallGLP::T level_mask =
        (typename SmallGLP::T{1} << dd_log_q) - typename SmallGLP::T{1};
    lhs.fill(level_mask);
    rhs.fill(level_mask);
    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(expected, lhs, rhs);
    TFHEpp::gl_detail::reduce<SmallGLP, dd_log_q>(expected);
    TFHEpp::gl_detail::baseMultiplyAtLevel<SmallGLP, dd_log_q>(got, lhs, rhs);
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

bool checkFullPolynomialDense()
{
    std::mt19937_64 rng(0x46554c4c4e5454ULL);
    TFHEpp::GLPolynomial<SmallGLP> lhs;
    TFHEpp::GLPolynomial<SmallGLP> rhs;
    TFHEpp::GLPolynomial<SmallGLP> expected;
    TFHEpp::GLPolynomial<SmallGLP> got;
    const auto equal = [&] {
        for (std::size_t y = 0; y < SmallGLP::matrix_dimension; y++)
            if (got[y] != expected[y]) return false;
        return true;
    };

    fillSigned<SmallGLP>(lhs, rng, 16);
    fillSigned<SmallGLP>(rhs, rng, 16);
    if (TFHEpp::gl_detail::polynomialMultiplyNTTPrimeCount<SmallGLP>(lhs,
                                                                     rhs) != 1)
        return false;
    TFHEpp::gl_detail::polynomialMultiplyReference<SmallGLP>(expected, lhs,
                                                             rhs);
    TFHEpp::gl_detail::polynomialMultiply<SmallGLP>(got, lhs, rhs);
    if (!equal()) return false;

    fillSigned<SmallGLP>(lhs, rng, 85);
    fillSigned<SmallGLP>(rhs, rng, 16);
    if (TFHEpp::gl_detail::polynomialMultiplyNTTPrimeCount<SmallGLP>(lhs,
                                                                     rhs) != 2)
        return false;
    TFHEpp::gl_detail::polynomialMultiplyReference<SmallGLP>(expected, lhs,
                                                             rhs);
    TFHEpp::gl_detail::polynomialMultiply<SmallGLP>(got, lhs, rhs);
    if (!equal()) return false;

    for (auto &slice : lhs)
        for (auto &coefficient : slice) {
            coefficient.limb[0] = rng();
            coefficient.limb[1] = rng();
        }
    for (std::size_t y = 0; y < SmallGLP::matrix_dimension; y++)
        for (std::size_t i = 0; i < SmallGLP::baseP::n; i++)
            rhs[y][i] = typename SmallGLP::T(
                (i + y) % 3 == 0 ? -1 : ((i + y) % 3 == 1 ? 0 : 1));
    if (TFHEpp::gl_detail::polynomialMultiplyNTTPrimeCount<SmallGLP>(lhs,
                                                                     rhs) != 0)
        return false;
    TFHEpp::gl_detail::polynomialMultiplyReference<SmallGLP>(expected, lhs,
                                                             rhs);
    if (!TFHEpp::gl_detail::polynomialMultiplyNTT<SmallGLP>(got, lhs, rhs))
        return false;
    return equal();
}

void referenceSmallSwitch(TFHEpp::GLCiphertextData<SmallGLP> &result,
                          const TFHEpp::GLPolynomial<SmallGLP> &input,
                          const SmallDDSwitchKey &switch_key)
{
    auto input_digits =
        TFHEpp::gl_detail::activeDecompose<SmallGLP, dd_log_q, dd_primary_bit>(
            input);
    std::array<std::vector<TFHEpp::GLPolynomial<SmallGLP>>, 2> digit_rows{
        std::vector<TFHEpp::GLPolynomial<SmallGLP>>(
            SmallDDSwitchKey::bbar_rows),
        std::vector<TFHEpp::GLPolynomial<SmallGLP>>(
            SmallDDSwitchKey::bbar_rows)};
    TFHEpp::GLBasePolynomial<SmallGLP> product{};
    TFHEpp::GLBasePolynomial<SmallGLP> key_row{};
    for (std::uint32_t primary = 0; primary < SmallDDSwitchKey::primary_rows;
         primary++) {
        for (std::uint32_t bbar = 0; bbar < SmallDDSwitchKey::bbar_rows;
             bbar++) {
            for (std::size_t component = 0; component < 2; component++) {
                TFHEpp::gl_detail::unpackDigitPolynomial<SmallGLP, dd_bbar_bit>(
                    key_row, switch_key.at(primary, bbar)[component]);
                for (std::uint32_t y = 0; y < SmallGLP::matrix_dimension; y++) {
                    TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(
                        product, input_digits[primary][y], key_row);
                    TFHEpp::gl_detail::addInPlace<SmallGLP>(
                        digit_rows[component][bbar][y], product);
                }
            }
        }
    }
    for (std::size_t component = 0; component < 2; component++) {
        TFHEpp::gl_detail::activeRecombine<
            SmallGLP, SmallDDSwitchKey::key_log_q, dd_bbar_bit>(
            result[component], digit_rows[component]);
        TFHEpp::gl_detail::divideRoundLevel<SmallGLP,
                                            SmallDDSwitchKey::key_log_q,
                                            SmallDDSwitchKey::auxiliary_log_q>(
            result[component]);
    }
}

void referenceBigSwitch(TFHEpp::GLCiphertextData<SmallGLP> &result,
                        const TFHEpp::GLPolynomial<SmallGLP> &input,
                        const SmallDDBigSwitchKey &switch_key)
{
    auto input_digits =
        TFHEpp::gl_detail::activeDecompose<SmallGLP, dd_log_q, dd_primary_bit>(
            input);
    std::array<std::vector<TFHEpp::GLPolynomial<SmallGLP>>, 2> digit_rows{
        std::vector<TFHEpp::GLPolynomial<SmallGLP>>(
            SmallDDBigSwitchKey::bbar_rows),
        std::vector<TFHEpp::GLPolynomial<SmallGLP>>(
            SmallDDBigSwitchKey::bbar_rows)};
    TFHEpp::GLPolynomial<SmallGLP> product;
    TFHEpp::GLPolynomial<SmallGLP> key_row;
    for (std::uint32_t primary = 0; primary < SmallDDBigSwitchKey::primary_rows;
         primary++) {
        for (std::uint32_t bbar = 0; bbar < SmallDDBigSwitchKey::bbar_rows;
             bbar++) {
            for (std::size_t component = 0; component < 2; component++) {
                TFHEpp::gl_detail::unpackDigitPolynomial<SmallGLP, dd_bbar_bit>(
                    key_row, switch_key.at(primary, bbar)[component]);
                TFHEpp::gl_detail::polynomialMultiplyReference<SmallGLP>(
                    product, input_digits[primary], key_row);
                TFHEpp::gl_detail::addInPlace<SmallGLP>(
                    digit_rows[component][bbar], product);
            }
        }
    }
    for (std::size_t component = 0; component < 2; component++) {
        TFHEpp::gl_detail::activeRecombine<
            SmallGLP, SmallDDBigSwitchKey::key_log_q, dd_bbar_bit>(
            result[component], digit_rows[component]);
        TFHEpp::gl_detail::divideRoundLevel<
            SmallGLP, SmallDDBigSwitchKey::key_log_q,
            SmallDDBigSwitchKey::auxiliary_log_q>(result[component]);
    }
}

bool checkBigDDSwitchAccumulation()
{
    using P = typename SmallGLP::baseP;
    std::mt19937_64 rng(0x42494744444e5454ULL);
    SmallDDBigSwitchKey switch_key;
    switch_key.allocate();
    for (auto &ciphertext : switch_key.data)
        for (std::size_t component_index = 0; component_index < 2;
             component_index++) {
            auto &component = ciphertext[component_index];
            for (std::size_t y = 0; y < component.size(); y++)
                for (auto &coefficient : component[y]) {
                    const std::int32_t value =
                        static_cast<std::int32_t>(rng() & 0xffffU) - 32768;
                    coefficient = static_cast<std::int16_t>(value);
                }
        }

    TFHEpp::GLPolynomial<SmallGLP> input;
    for (auto &slice : input)
        for (auto &coefficient : slice)
            coefficient = TFHEpp::ckks_detail::signedToLevel<P, dd_log_q>(
                randomSigned(rng, 90));

    TFHEpp::GLCiphertextData<SmallGLP> expected;
    TFHEpp::GLCiphertextData<SmallGLP> got;
    referenceBigSwitch(expected, input, switch_key);
    TFHEpp::GLDDBigKeySwitch(got, input, switch_key);
    for (std::size_t component = 0; component < 2; component++)
        for (std::size_t y = 0; y < SmallGLP::matrix_dimension; y++)
            if (got[component][y] != expected[component][y]) return false;
    return true;
}

template <std::uint32_t LogQ, std::uint32_t PrimaryBit, std::uint32_t BbarBit>
void referenceSmallSwitchBaseRaw(
    TFHEpp::GLBaseCiphertextData<SmallGLP> &result,
    const TFHEpp::GLBasePolynomial<SmallGLP> &input,
    const TFHEpp::GLDDSmallKeySwitchKey<SmallGLP, LogQ, PrimaryBit, BbarBit>
        &switch_key)
{
    using SwitchKey =
        TFHEpp::GLDDSmallKeySwitchKey<SmallGLP, LogQ, PrimaryBit, BbarBit>;
    auto input_digits =
        TFHEpp::gl_detail::activeDecompose<SmallGLP, LogQ, PrimaryBit>(input);
    std::array<std::vector<TFHEpp::GLBasePolynomial<SmallGLP>>, 2> digit_rows{
        std::vector<TFHEpp::GLBasePolynomial<SmallGLP>>(SwitchKey::bbar_rows),
        std::vector<TFHEpp::GLBasePolynomial<SmallGLP>>(SwitchKey::bbar_rows)};
    TFHEpp::GLBasePolynomial<SmallGLP> product{};
    TFHEpp::GLBasePolynomial<SmallGLP> key_row{};
    for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
         primary++) {
        for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++) {
            for (std::size_t component = 0; component < 2; component++) {
                TFHEpp::gl_detail::unpackDigitPolynomial<SmallGLP, BbarBit>(
                    key_row, switch_key.at(primary, bbar)[component]);
                TFHEpp::gl_detail::baseMultiplyReference<SmallGLP>(
                    product, input_digits[primary], key_row);
                TFHEpp::gl_detail::addInPlace<SmallGLP>(
                    digit_rows[component][bbar], product);
            }
        }
    }
    for (std::size_t component = 0; component < 2; component++)
        TFHEpp::gl_detail::activeRecombine<SmallGLP, SwitchKey::key_log_q,
                                           BbarBit>(result[component],
                                                    digit_rows[component]);
}

bool checkDDSwitchAccumulation()
{
    using P = typename SmallGLP::baseP;
    std::mt19937_64 rng(0x44444e5454ULL);
    SmallDDSwitchKey switch_key;
    switch_key.allocate();
    for (auto &ciphertext : switch_key.data)
        for (auto &component : ciphertext)
            for (auto &coefficient : component) {
                const std::int32_t value =
                    static_cast<std::int32_t>(rng() & 0xffffU) - 32768;
                coefficient = static_cast<std::int16_t>(value);
            }

    TFHEpp::GLPolynomial<SmallGLP> input;
    for (auto &slice : input)
        for (auto &coefficient : slice)
            coefficient = TFHEpp::ckks_detail::signedToLevel<P, dd_log_q>(
                randomSigned(rng, 90));

    TFHEpp::GLCiphertextData<SmallGLP> expected;
    TFHEpp::GLCiphertextData<SmallGLP> got;
    referenceSmallSwitch(expected, input, switch_key);
    TFHEpp::GLDDSmallKeySwitch(got, input, switch_key);
    for (std::size_t component = 0; component < 2; component++)
        for (std::size_t y = 0; y < SmallGLP::matrix_dimension; y++)
            if (got[component][y] != expected[component][y]) return false;

    TFHEpp::GLBaseCiphertextData<SmallGLP> expected_raw{};
    TFHEpp::GLBaseCiphertextData<SmallGLP> got_raw{};
    referenceSmallSwitchBaseRaw<dd_log_q, dd_primary_bit, dd_bbar_bit>(
        expected_raw, input[0], switch_key);
    TFHEpp::GLDDSmallKeySwitchBaseRaw(got_raw, input[0], switch_key);
    if (got_raw != expected_raw) return false;

    TFHEpp::GLBaseCiphertextData<SmallGLP> second_raw{};
    TFHEpp::GLBaseCiphertextData<SmallGLP> expected_sum = expected_raw;
    referenceSmallSwitchBaseRaw<dd_log_q, dd_primary_bit, dd_bbar_bit>(
        second_raw, input[1], switch_key);
    TFHEpp::gl_detail::addInPlace<SmallGLP>(expected_sum[0], second_raw[0]);
    TFHEpp::gl_detail::addInPlace<SmallGLP>(expected_sum[1], second_raw[1]);
    TFHEpp::gl_detail::reduce<SmallGLP, SmallDDSwitchKey::key_log_q>(
        expected_sum[0]);
    TFHEpp::gl_detail::reduce<SmallGLP, SmallDDSwitchKey::key_log_q>(
        expected_sum[1]);
    const std::array<const TFHEpp::GLBasePolynomial<SmallGLP> *, 2> sum_inputs{
        &input[0], &input[1]};
    const std::array<const SmallDDSwitchKey *, 2> sum_keys{&switch_key,
                                                           &switch_key};
    TFHEpp::GLBaseCiphertextData<SmallGLP> got_sum{};
    if (!TFHEpp::gl_detail::accumulateSmallKeySwitchSumNTT<SmallGLP,
                                                           SmallDDSwitchKey>(
            got_sum, sum_inputs, sum_keys) ||
        got_sum != expected_sum)
        return false;
#ifdef USE_HEXL
    TFHEpp::GLBaseCiphertextData<SmallGLP> blocked_sum{};
    TFHEpp::gl_detail::SmallKeySwitchSumNTTWorkspace<SmallGLP,
                                                     SmallDDSwitchKey, 2>
        blocked_workspace;
    blocked_workspace.use_batched_vector_mac = true;
    blocked_workspace.use_coefficient_blocked_key_layout = true;
    if (!TFHEpp::gl_detail::accumulateSmallKeySwitchSumNTT<SmallGLP,
                                                           SmallDDSwitchKey>(
            blocked_sum, sum_inputs, sum_keys, &blocked_workspace) ||
        blocked_sum != expected_sum)
        return false;
#endif

    std::cout << "n512 DD transient NTT key cache: "
              << TFHEpp::gl_detail::smallKeySwitchNTTKeyCacheBytes<
                     TFHEpp::GL512p17Parameter, N512DDSwitchKey> /
                     (1024.0 * 1024.0)
              << " MiB" << std::endl;
    return true;
}

bool checkOnePrimeDDSwitchAccumulation()
{
    using P = typename SmallGLP::baseP;
    std::mt19937_64 rng(0x3144444e5454ULL);
    SmallOnePrimeDDSwitchKey switch_key;
    switch_key.allocate();
    for (auto &ciphertext : switch_key.data)
        for (auto &component : ciphertext)
            for (auto &coefficient : component) {
                const std::int32_t value =
                    static_cast<std::int32_t>(rng() & 0xffU) - 128;
                coefficient = static_cast<std::int8_t>(value);
            }

    TFHEpp::GLBasePolynomial<SmallGLP> input{};
    for (auto &coefficient : input)
        coefficient = TFHEpp::ckks_detail::signedToLevel<P, one_prime_log_q>(
            randomSigned(rng, 30));

    TFHEpp::GLBaseCiphertextData<SmallGLP> expected{};
    TFHEpp::GLBaseCiphertextData<SmallGLP> got{};
    referenceSmallSwitchBaseRaw<one_prime_log_q, one_prime_primary_bit,
                                one_prime_bbar_bit>(expected, input,
                                                    switch_key);
    TFHEpp::GLDDSmallKeySwitchBaseRaw(got, input, switch_key);
    return got == expected;
}

bool checkFusedWTransform()
{
    using Schedule = SmallWTransformSchedule;
    using P = typename SmallGLP::baseP;
    using AfterX =
        TFHEpp::GLRawProductCiphertext<SmallGLP, Schedule::input_log_q,
                                       Schedule::input_log_delta,
                                       Schedule::x_transform_log_scale>;
    using BeforeRescale =
        TFHEpp::GLRawProductCiphertext<SmallGLP, Schedule::input_log_q,
                                       Schedule::input_log_delta +
                                           Schedule::x_transform_log_scale,
                                       Schedule::w_transform_log_scale>;
    std::mt19937_64 rng(0x575452414e53464dULL);
    std::vector<AfterX> baby_rotations(Schedule::w_baby_step);
    for (auto &ciphertext : baby_rotations)
        for (std::size_t component = 0; component < 2; component++)
            for (auto &slice : ciphertext[component])
                for (auto &coefficient : slice)
                    coefficient = TFHEpp::ckks_detail::signedToLevel<
                        P, Schedule::input_log_q>(randomSigned(rng, 42));

    TFHEpp::gl_ship_detail::FusedWTransformNTTState<Schedule> state;
    if (!TFHEpp::gl_ship_detail::prepareFusedWTransformNTT<Schedule>(
            state, baby_rotations))
        return false;
    for (std::uint32_t b = 0; b < Schedule::w_giant_steps; b++) {
        const std::uint32_t begin = b * Schedule::w_baby_step;
        const std::uint32_t end = std::min<std::uint32_t>(
            SmallGLP::phi, begin + Schedule::w_baby_step);
        BeforeRescale expected;
        bool initialized = false;
        for (std::uint32_t t = begin; t < end; t++) {
            TFHEpp::GLBasePlaintext<SmallGLP, Schedule::input_log_q,
                                    Schedule::w_transform_log_scale>
                diagonal;
            TFHEpp::gl_ship_detail::buildAdjustedWDiagonalPlaintext<Schedule>(
                diagonal, t, begin);
            BeforeRescale term;
            TFHEpp::gl_ship_detail::basePlaintextHadamardMultiplyRaw(
                term, baby_rotations[t - begin], diagonal);
            if (!initialized) {
                expected = std::move(term);
                initialized = true;
            }
            else {
                TFHEpp::gl_ship_detail::addInPlace(expected, term);
            }
        }
        BeforeRescale got;
        TFHEpp::gl_ship_detail::fusedWTransformGroup<Schedule>(got, state, b);
        for (std::size_t component = 0; component < 2; component++)
            for (std::size_t y = 0; y < SmallGLP::matrix_dimension; y++)
                if (got[component][y] != expected[component][y]) return false;
    }
    return true;
}

bool checkFusedMaskedColumn()
{
    using Schedule = SmallWTransformSchedule;
    using P = typename SmallGLP::baseP;
    constexpr std::uint32_t key_log_q =
        Schedule::half_bootstrap_log_q + SmallGLP::auxiliary_log_q;
    std::mt19937_64 rng(0x4d41534b4e545446ULL);
    TFHEpp::GLSHIPMaskedColumnKey<Schedule> key;
    key.candidates = {{0, 0, 0}, {1, 3, 1}, {0, 7, 3}};
    key.encrypted_masks.resize(key.candidates.size());
    for (auto &ciphertext : key.encrypted_masks)
        for (auto &component : ciphertext)
            for (auto &coefficient : component) {
                coefficient.limb[0] = rng();
                coefficient.limb[1] = rng();
            }
    TFHEpp::GLBasePolynomial<SmallGLP> mask{};
    for (auto &coefficient : mask)
        coefficient = TFHEpp::ckks_detail::signedToLevel<P, Schedule::q0_log_q>(
            randomSigned(rng, 28));

    std::vector<std::complex<long double>> phase_roots;
    TFHEpp::gl_ship_detail::prepareCandidatePhaseRoots<Schedule>(phase_roots,
                                                                 mask);
    for (std::uint32_t fine_x = 0; fine_x < Schedule::theta; fine_x++)
        for (std::uint32_t w = 0; w < SmallGLP::phi; w++)
            for (std::uint32_t gaussian_phase = 0; gaussian_phase < 4;
                 gaussian_phase++)
                for (std::uint32_t channel = 0; channel < 2; channel++) {
                    const TFHEpp::GLSHIPCandidate descriptor{fine_x, w,
                                                             gaussian_phase};
                    TFHEpp::GLBasePlaintext<SmallGLP, key_log_q,
                                            Schedule::tree_log_delta>
                        direct, hoisted;
                    TFHEpp::gl_ship_detail::buildCandidatePlaintext<Schedule>(
                        direct, mask, descriptor, channel);
                    TFHEpp::gl_ship_detail::
                        buildCandidatePlaintextFromPhaseRoots<Schedule>(
                            hoisted, phase_roots, descriptor, channel);
                    if (hoisted.poly != direct.poly) return false;

                    TFHEpp::GLBasePlaintext<SmallGLP, key_log_q,
                                            Schedule::tree_log_delta>
                        unshifted;
                    TFHEpp::gl_ship_detail::
                        buildCandidatePlaintextFromPhaseRoots<Schedule>(
                            unshifted, phase_roots, {0, w, gaussian_phase},
                            channel);
                    TFHEpp::GLBasePolynomial<SmallGLP> shifted;
                    const std::uint32_t multiplier = TFHEpp::gl_detail::powMod(
                        5,
                        (SmallGLP::matrix_dimension - fine_x) %
                            SmallGLP::matrix_dimension,
                        4 * SmallGLP::matrix_dimension);
                    TFHEpp::gl_detail::baseAutomorphism<SmallGLP>(
                        shifted, unshifted.poly, multiplier, 1);
                    if (shifted != direct.poly) return false;
                }

    TFHEpp::GLBaseCiphertextData<SmallGLP> reference_raw{};
    for (std::size_t i = 0; i < key.candidates.size(); i++) {
        TFHEpp::GLBasePlaintext<SmallGLP, key_log_q, Schedule::tree_log_delta>
            candidate;
        TFHEpp::gl_ship_detail::buildCandidatePlaintext<Schedule>(
            candidate, mask, key.candidates[i], 1);
        TFHEpp::GLBaseCiphertextData<SmallGLP> term{};
        for (std::size_t component = 0; component < 2; component++)
            TFHEpp::gl_detail::baseMultiply<SmallGLP>(
                term[component], key.encrypted_masks[i][component],
                candidate.poly);
        TFHEpp::gl_ship_detail::addInPlace<SmallGLP>(reference_raw, term);
    }
    TFHEpp::gl_ship_detail::reduce<SmallGLP, key_log_q>(reference_raw);
    TFHEpp::GLBaseCiphertext<SmallGLP, Schedule::half_bootstrap_log_q,
                             Schedule::tree_log_delta>
        expected;
    TFHEpp::gl_ship_detail::rescaleBase<SmallGLP, key_log_q,
                                        SmallGLP::auxiliary_log_q>(
        expected.ct, reference_raw);

    decltype(expected) got;
    if (!TFHEpp::gl_ship_detail::maskedColumnNTT<Schedule>(got, mask, 1, key))
        return false;
    if (got.ct != expected.ct) return false;
    TFHEpp::gl_ship_detail::releaseMaskedColumnNTTCache<Schedule>(key);
    return true;
}

bool benchmarkProductionDDBaseSwitch()
{
    using GLP = TFHEpp::GL512p17Parameter;
    using P = typename GLP::baseP;
    using SwitchKey = N512DDSwitchKey;
    std::mt19937_64 rng(0x4e3531324444ULL);

    auto switch_key = std::make_unique<SwitchKey>();
    switch_key->allocate();
    for (auto &ciphertext : switch_key->data)
        for (auto &component : ciphertext)
            for (auto &coefficient : component) {
                const std::int32_t value =
                    static_cast<std::int32_t>(rng() & 0xffffU) - 32768;
                coefficient = static_cast<std::int16_t>(value);
            }

    auto input = std::make_unique<TFHEpp::GLBasePolynomial<GLP>>();
    for (auto &coefficient : *input) {
        for (auto &limb : coefficient.limb) limb = rng();
        coefficient = TFHEpp::ckks_detail::reduceToLevel<P, SwitchKey::log_q>(
            coefficient);
    }
    auto output = std::make_unique<TFHEpp::GLBaseCiphertextData<GLP>>();
    auto benchmark_digits = std::make_unique<
        std::array<TFHEpp::GLBasePolynomial<GLP>, SwitchKey::primary_rows>>();
    auto benchmark_output =
        std::make_unique<TFHEpp::GLBaseCiphertextData<GLP>>();
    const auto cache_start = std::chrono::steady_clock::now();
    const bool used_ntt =
        TFHEpp::gl_detail::accumulateSmallKeySwitchProductsNTT<GLP, SwitchKey>(
            0, *switch_key, [&](const std::size_t) {},
            [&](const std::uint32_t,
                const std::size_t) -> const TFHEpp::GLBasePolynomial<GLP> & {
                return (*benchmark_digits)[0];
            },
            [&](const std::size_t, const std::uint32_t, const std::size_t,
                const std::size_t, const typename GLP::T) {});
    const double cache_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      cache_start)
            .count();
    if (!used_ntt) return false;

    constexpr std::size_t measured_slices = 8;
    const auto slices_start = std::chrono::steady_clock::now();
    if (!TFHEpp::gl_detail::accumulateSmallKeySwitchProductsNTT<GLP, SwitchKey>(
            measured_slices, *switch_key,
            [&](const std::size_t) {
                TFHEpp::ckks_detail::activeBaseDecomposePolynomialRows<
                    P, SwitchKey::log_q, SwitchKey::primary_bit,
                    SwitchKey::primary_rows>(*benchmark_digits, *input);
                TFHEpp::gl_detail::clear<GLP>((*benchmark_output)[0]);
                TFHEpp::gl_detail::clear<GLP>((*benchmark_output)[1]);
            },
            [&](const std::uint32_t primary,
                const std::size_t) -> const TFHEpp::GLBasePolynomial<GLP> & {
                return (*benchmark_digits)[primary];
            },
            [&](const std::size_t component, const std::uint32_t bbar,
                const std::size_t, const std::size_t coefficient,
                const typename GLP::T value) {
                const std::uint32_t shift =
                    (SwitchKey::bbar_rows - bbar - 1) * SwitchKey::bbar_bit;
                (*benchmark_output)[component][coefficient] += value << shift;
            }))
        return false;
    const double slices_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      slices_start)
            .count();
    const double per_slice_seconds = slices_seconds / measured_slices;
    const double projected_full_seconds =
        cache_seconds + GLP::matrix_dimension * per_slice_seconds;

    const auto start = std::chrono::steady_clock::now();
    TFHEpp::GLDDSmallKeySwitchBaseRaw(*output, *input, *switch_key);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cout << "n512 dense raw DD base switch: " << seconds
              << " s (key transforms " << cache_seconds
              << " s; projected full-Y switch " << projected_full_seconds
              << " s)" << std::endl;
    if (std::getenv("TFHEPP_GL_N512_SMALL_FULL_BENCH") != nullptr) {
        auto full_input = std::make_unique<TFHEpp::GLPolynomial<GLP>>();
        for (std::size_t y = 0; y < GLP::matrix_dimension; y++)
            (*full_input)[y] = *input;
        auto full_output = std::make_unique<TFHEpp::GLCiphertextData<GLP>>();
        const auto full_start = std::chrono::steady_clock::now();
        TFHEpp::GLDDSmallKeySwitch(*full_output, *full_input, *switch_key);
        const double full_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          full_start)
                .count();
        std::cout << "n512 dense full-Y DD small switch: " << full_seconds
                  << " s" << std::endl;
        if ((*full_output)[0][0][0] == typename GLP::T{} &&
            (*full_output)[1][0][0] == typename GLP::T{})
            return false;
    }
    return (*output)[0][0] != typename GLP::T{} ||
           (*output)[1][0] != typename GLP::T{};
}

bool benchmarkProductionFullTransform()
{
    using GLP = TFHEpp::GL512p17Parameter;
    std::mt19937_64 rng(0x4e35313246554c4cULL);
    auto packed =
        std::make_unique<TFHEpp::GLPackedPolynomial<GLP, dd_bbar_bit>>();
    for (std::size_t y = 0; y < packed->size(); y++)
        for (auto &coefficient : (*packed)[y])
            coefficient = static_cast<std::int16_t>(rng());

    const auto &plan = TFHEpp::gl_detail::polynomialNTTPlan<GLP, 0>();
    std::vector<std::uint64_t> spectrum;
    std::vector<std::uint64_t> coefficients;
    const auto forward_start = std::chrono::steady_clock::now();
    plan.forwardPacked(spectrum, *packed);
    const double forward_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      forward_start)
            .count();
    const auto inverse_start = std::chrono::steady_clock::now();
    plan.inverse(coefficients, spectrum);
    const double inverse_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      inverse_start)
            .count();

    for (std::size_t sample = 0; sample < 64; sample++) {
        const std::size_t index = (sample * 104729 + 17) % coefficients.size();
        const std::int64_t value =
            (*packed)[index / GLP::baseP::n][index % GLP::baseP::n];
        const std::uint64_t expected =
            value >= 0 ? static_cast<std::uint64_t>(value)
                       : plan.modulus() - static_cast<std::uint64_t>(-value);
        if (coefficients[index] != expected) return false;
    }
    std::cout << "n512 p17 full-ring packed forward/inverse: "
              << forward_seconds << " s / " << inverse_seconds << " s"
              << std::endl;
    return true;
}

bool benchmarkProductionHMuxSwitchFusion()
{
    if (std::getenv("TFHEPP_GL_N512_HMUX_BENCH") == nullptr) return true;
    using GLP = TFHEpp::GL512p17Parameter;
    using P = typename GLP::baseP;
    using SwitchKey = N512DDSwitchKey;
    constexpr std::size_t switch_count = 8;
    std::mt19937_64 rng(0x4e353132484d5558ULL);
    std::array<std::unique_ptr<SwitchKey>, switch_count> keys;
    std::array<std::unique_ptr<TFHEpp::GLBasePolynomial<GLP>>, switch_count>
        inputs;
    std::array<const SwitchKey *, switch_count> key_pointers{};
    std::array<const TFHEpp::GLBasePolynomial<GLP> *, switch_count>
        input_pointers{};
    for (std::size_t term = 0; term < switch_count; term++) {
        keys[term] = std::make_unique<SwitchKey>();
        keys[term]->allocate();
        for (auto &ciphertext : keys[term]->data)
            for (auto &component : ciphertext)
                for (auto &coefficient : component)
                    coefficient = static_cast<std::int16_t>(rng());
        inputs[term] = std::make_unique<TFHEpp::GLBasePolynomial<GLP>>();
        for (auto &coefficient : *inputs[term]) {
            for (auto &limb : coefficient.limb) limb = rng();
            coefficient =
                TFHEpp::ckks_detail::reduceToLevel<P, SwitchKey::log_q>(
                    coefficient);
        }
        key_pointers[term] = keys[term].get();
        input_pointers[term] = inputs[term].get();
    }

    const auto cache_start = std::chrono::steady_clock::now();
    for (const auto &key : keys)
        TFHEpp::gl_detail::prepareSmallKeySwitchNTTCache<GLP, SwitchKey>(*key);
    const double cache_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      cache_start)
            .count();
#ifdef USE_HEXL
    const auto blocked_cache_start = std::chrono::steady_clock::now();
    for (const auto &key : keys)
        TFHEpp::gl_detail::prepareSmallKeySwitchBlockedNTTCache<GLP,
                                                                SwitchKey>(
            *key);
    const double blocked_cache_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      blocked_cache_start)
            .count();
#endif

    TFHEpp::GLBaseCiphertextData<GLP> fused{};
    const auto fused_start = std::chrono::steady_clock::now();
    if (!TFHEpp::gl_detail::accumulateSmallKeySwitchSumNTT<GLP, SwitchKey>(
            fused, input_pointers, key_pointers))
        return false;
    const double fused_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      fused_start)
            .count();

    std::size_t batch_count = 1;
#ifdef _OPENMP
    batch_count = 8 * static_cast<std::size_t>(omp_get_max_threads());
#endif
    std::vector<TFHEpp::GLBaseCiphertextData<GLP>> fused_batch(batch_count);
    const auto batch_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::gl_detail::SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey,
                                                         switch_count>
            workspace;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            if (!TFHEpp::gl_detail::accumulateSmallKeySwitchSumNTT<GLP,
                                                                   SwitchKey>(
                    fused_batch[batch], input_pointers, key_pointers,
                    &workspace))
                std::terminate();
    }
    const double batch_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      batch_start)
            .count();
    for (const auto &batch : fused_batch)
        if (batch != fused) return false;
#ifdef USE_HEXL
    std::vector<TFHEpp::GLBaseCiphertextData<GLP>> batched_mac_batch(
        batch_count);
    const auto batched_mac_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::gl_detail::SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey,
                                                         switch_count>
            workspace;
        workspace.use_batched_vector_mac = true;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            if (!TFHEpp::gl_detail::accumulateSmallKeySwitchSumNTT<GLP,
                                                                   SwitchKey>(
                    batched_mac_batch[batch], input_pointers, key_pointers,
                    &workspace))
                std::terminate();
    }
    const double batched_mac_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      batched_mac_start)
            .count();
    for (const auto &batch : batched_mac_batch)
        if (batch != fused) return false;
    std::vector<TFHEpp::GLBaseCiphertextData<GLP>> blocked_mac_batch(
        batch_count);
    const auto blocked_mac_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::gl_detail::SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey,
                                                         switch_count>
            workspace;
        workspace.use_batched_vector_mac = true;
        workspace.use_coefficient_blocked_key_layout = true;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            if (!TFHEpp::gl_detail::accumulateSmallKeySwitchSumNTT<GLP,
                                                                   SwitchKey>(
                    blocked_mac_batch[batch], input_pointers, key_pointers,
                    &workspace))
                std::terminate();
    }
    const double blocked_mac_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      blocked_mac_start)
            .count();
    for (const auto &batch : blocked_mac_batch)
        if (batch != fused) return false;
#endif

    TFHEpp::GLBaseCiphertextData<GLP> separate{};
    const auto separate_start = std::chrono::steady_clock::now();
    for (std::size_t term = 0; term < switch_count; term++) {
        TFHEpp::GLBaseCiphertextData<GLP> product{};
        TFHEpp::GLDDSmallKeySwitchBaseRaw(product, *inputs[term], *keys[term]);
        TFHEpp::gl_detail::addInPlace<GLP>(separate[0], product[0]);
        TFHEpp::gl_detail::addInPlace<GLP>(separate[1], product[1]);
    }
    TFHEpp::gl_detail::reduce<GLP, SwitchKey::key_log_q>(separate[0]);
    TFHEpp::gl_detail::reduce<GLP, SwitchKey::key_log_q>(separate[1]);
    const double separate_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      separate_start)
            .count();
    std::cout << "n512 HMux eight-switch fusion: " << fused_seconds << " s vs "
              << separate_seconds << " s separate (" << cache_seconds
              << " s cache preparation); " << batch_count
              << " warm calls with reused per-thread scratch in "
              << batch_seconds << " s";
#ifdef USE_HEXL
    std::cout << " vs " << batched_mac_seconds << " s batched-MAC";
    std::cout << " vs " << blocked_mac_seconds << " s blocked-MAC ("
              << blocked_cache_seconds << " s blocked cache preparation)";
#endif
    std::cout << std::endl;
    if (fused != separate) return false;

    using Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;
    using Ciphertext =
        TFHEpp::GLBaseCiphertext<GLP, Schedule::half_bootstrap_log_q,
                                 Schedule::tree_log_delta>;
    TFHEpp::GLSHIPHMuxKey<Schedule> hmux_key;
    hmux_key.stages.resize(1);
    hmux_key.stages[0].step = 1;
    hmux_key.stages[0].branches.resize(Schedule::hmux_radix);
    for (std::size_t digit = 0; digit < Schedule::hmux_radix; digit++) {
        hmux_key.stages[0].branches[digit].body_key =
            std::move(*keys[2 * digit]);
        hmux_key.stages[0].branches[digit].mask_key =
            std::move(*keys[2 * digit + 1]);
    }
    Ciphertext hmux_input;
    hmux_input[0] = *inputs[0];
    hmux_input[1] = *inputs[1];
    Ciphertext hmux_output;
    TFHEpp::gl_ship_detail::HMuxNTTWorkspace<Schedule> hmux_workspace;
    const auto hmux_start = std::chrono::steady_clock::now();
    TFHEpp::gl_ship_detail::hmux<Schedule>(hmux_output, hmux_input, hmux_key,
                                           &hmux_workspace);
    const double hmux_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      hmux_start)
            .count();
    std::vector<Ciphertext> legacy_hmux_batch(batch_count);
    const auto legacy_hmux_batch_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        auto rotated =
            std::make_unique<std::array<TFHEpp::GLBaseCiphertextData<GLP>,
                                        Schedule::hmux_radix>>();
        auto accumulated =
            std::make_unique<TFHEpp::GLBaseCiphertextData<GLP>>();
        TFHEpp::gl_detail::SmallKeySwitchSumNTTWorkspace<
            GLP, SwitchKey, 2 * Schedule::hmux_radix>
            workspace;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++) {
            std::array<const TFHEpp::GLBasePolynomial<GLP> *,
                       2 * Schedule::hmux_radix>
                switch_inputs{};
            std::array<const SwitchKey *, 2 * Schedule::hmux_radix>
                switch_keys{};
            const auto &stage = hmux_key.stages.front();
            for (std::uint32_t digit = 0; digit < Schedule::hmux_radix;
                 digit++) {
                const std::uint32_t displacement =
                    (digit * stage.step) % GLP::matrix_dimension;
                const std::uint32_t amount =
                    (GLP::matrix_dimension - displacement) %
                    GLP::matrix_dimension;
                TFHEpp::gl_ship_detail::rotateX<GLP>((*rotated)[digit],
                                                     hmux_input.ct, amount);
                switch_inputs[2 * digit] = &(*rotated)[digit][0];
                switch_inputs[2 * digit + 1] = &(*rotated)[digit][1];
                switch_keys[2 * digit] = &stage.branches[digit].body_key;
                switch_keys[2 * digit + 1] = &stage.branches[digit].mask_key;
            }
            if (!TFHEpp::gl_detail::accumulateSmallKeySwitchSumNTT<GLP,
                                                                   SwitchKey>(
                    *accumulated, switch_inputs, switch_keys, &workspace))
                std::terminate();
            TFHEpp::gl_ship_detail::reduce<GLP, Schedule::half_bootstrap_log_q +
                                                    GLP::auxiliary_log_q>(
                *accumulated);
            TFHEpp::gl_ship_detail::rescaleBase<
                GLP, Schedule::half_bootstrap_log_q + GLP::auxiliary_log_q,
                GLP::auxiliary_log_q>(legacy_hmux_batch[batch].ct,
                                      *accumulated);
        }
    }
    const double legacy_hmux_batch_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      legacy_hmux_batch_start)
            .count();
    for (const auto &batch : legacy_hmux_batch)
        if (batch.ct != hmux_output.ct) return false;
    std::vector<Ciphertext> hmux_batch(batch_count);
    const auto hmux_batch_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::gl_ship_detail::HMuxNTTWorkspace<Schedule> workspace;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            TFHEpp::gl_ship_detail::hmux<Schedule>(
                hmux_batch[batch], hmux_input, hmux_key, &workspace);
    }
    const double hmux_batch_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      hmux_batch_start)
            .count();
    for (const auto &batch : hmux_batch)
        if (batch.ct != hmux_output.ct) return false;
#ifdef USE_HEXL
    std::vector<Ciphertext> batched_hmux_batch(batch_count);
    const auto batched_hmux_batch_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::gl_ship_detail::HMuxNTTWorkspace<Schedule> workspace;
        workspace.switch_workspace.use_batched_vector_mac = true;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            TFHEpp::gl_ship_detail::hmux<Schedule>(
                batched_hmux_batch[batch], hmux_input, hmux_key, &workspace);
    }
    const double batched_hmux_batch_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      batched_hmux_batch_start)
            .count();
    for (const auto &batch : batched_hmux_batch)
        if (batch.ct != hmux_output.ct) return false;
    std::vector<Ciphertext> blocked_hmux_batch(batch_count);
    const auto blocked_hmux_batch_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::gl_ship_detail::HMuxNTTWorkspace<Schedule> workspace;
        workspace.switch_workspace.use_batched_vector_mac = true;
        workspace.switch_workspace.use_coefficient_blocked_key_layout = true;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            TFHEpp::gl_ship_detail::hmux<Schedule>(
                blocked_hmux_batch[batch], hmux_input, hmux_key, &workspace);
    }
    const double blocked_hmux_batch_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      blocked_hmux_batch_start)
            .count();
    for (const auto &batch : blocked_hmux_batch)
        if (batch.ct != hmux_output.ct) return false;
#endif
    std::cout << "n512 complete warm HMux stage: " << hmux_seconds << " s; "
              << batch_count << " calls in " << hmux_batch_seconds << " s";
#ifdef USE_HEXL
    std::cout << " vs " << batched_hmux_batch_seconds << " s batched-MAC";
    std::cout << " vs " << blocked_hmux_batch_seconds << " s blocked-MAC";
#endif
    std::cout << " vs " << legacy_hmux_batch_seconds
              << " s with unhoisted transforms" << std::endl;
    return hmux_output[0][0] != typename GLP::T{} ||
           hmux_output[1][0] != typename GLP::T{};
}

bool benchmarkProductionXTrace()
{
    if (std::getenv("TFHEPP_GL_N512_X_TRACE_BENCH") == nullptr) return true;
    using Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;
    using GLP = typename Schedule::Parameter;
    std::mt19937_64 rng(0x4e35313258545243ULL);
    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    for (std::size_t component = 0; component < 2; component++)
        for (auto &slice : (*input)[component])
            for (auto &coefficient : slice) {
                coefficient.limb[0] = rng();
                coefficient.limb[1] = rng() & ((std::uint64_t{1} << 21) - 1);
            }
    using Product =
        TFHEpp::GLRawProductCiphertext<GLP, Schedule::input_log_q,
                                       Schedule::input_log_delta,
                                       Schedule::x_transform_log_scale>;
    auto output = std::make_unique<Product>();
    const auto start = std::chrono::steady_clock::now();
    TFHEpp::GLXTransformMatrixMultiplyRaw<Schedule>(*output, *input);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cout << "n512 compact exact X-transform trace: " << seconds << " s"
              << std::endl;
    return (*output)[0][0][0] != typename GLP::T{} ||
           (*output)[1][0][0] != typename GLP::T{};
}

bool benchmarkProductionWTransform()
{
    if (std::getenv("TFHEPP_GL_N512_W_BENCH") == nullptr) return true;
    using Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;
    using GLP = typename Schedule::Parameter;
    using AfterX =
        TFHEpp::GLRawProductCiphertext<GLP, Schedule::input_log_q,
                                       Schedule::input_log_delta,
                                       Schedule::x_transform_log_scale>;
    using BeforeRescale =
        TFHEpp::GLRawProductCiphertext<GLP, Schedule::input_log_q,
                                       Schedule::input_log_delta +
                                           Schedule::x_transform_log_scale,
                                       Schedule::w_transform_log_scale>;
    std::mt19937_64 rng(0x4e3531325754524eULL);
    std::vector<AfterX> baby_rotations(Schedule::w_baby_step);
    for (auto &ciphertext : baby_rotations)
        for (std::size_t component = 0; component < 2; component++)
            for (auto &slice : ciphertext[component])
                for (auto &coefficient : slice) {
                    coefficient.limb[0] = rng();
                    coefficient.limb[1] =
                        rng() & ((std::uint64_t{1} << 21) - 1);
                }

    TFHEpp::gl_ship_detail::FusedWTransformNTTState<Schedule> state;
    const auto prepare_start = std::chrono::steady_clock::now();
    if (!TFHEpp::gl_ship_detail::prepareFusedWTransformNTT<Schedule>(
            state, baby_rotations))
        return false;
    const double prepare_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      prepare_start)
            .count();
    auto group = std::make_unique<BeforeRescale>();
    const auto groups_start = std::chrono::steady_clock::now();
    for (std::uint32_t b = 0; b < Schedule::w_giant_steps; b++)
        TFHEpp::gl_ship_detail::fusedWTransformGroup<Schedule>(*group, state,
                                                               b);
    const double groups_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      groups_start)
            .count();
    std::cout << "n512 fused W transform: " << groups_seconds << " s for "
              << Schedule::w_giant_steps << " groups (" << prepare_seconds
              << " s spectrum preparation)" << std::endl;
    return (*group)[0][0][0] != typename GLP::T{} ||
           (*group)[1][0][0] != typename GLP::T{};
}

bool benchmarkProductionBatchRotation()
{
    if (std::getenv("TFHEPP_GL_N512_ROTATION_BENCH") == nullptr) return true;
    using Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;
    using GLP = typename Schedule::Parameter;
    using RotationKey =
        TFHEpp::GLBatchRotationKey<GLP, Schedule::input_log_q,
                                   Schedule::primary_bit, Schedule::bbar_bit>;
    using SwitchKey = TFHEpp::GLDDSmallKeySwitchKey<
        GLP, Schedule::input_log_q, Schedule::primary_bit, Schedule::bbar_bit>;
    std::mt19937_64 rng(0x4e353132524f544eULL);
    auto rotation_key = std::make_unique<RotationKey>();
    rotation_key->amount = 1;
    rotation_key->multiplier = TFHEpp::gl_detail::powMod(GLP::primitive_root, 1,
                                                         GLP::cyclotomic_order);
    rotation_key->switch_key.allocate();
    for (auto &ciphertext : rotation_key->switch_key.data)
        for (auto &component : ciphertext)
            for (auto &coefficient : component)
                coefficient = static_cast<std::int16_t>(rng());

    using AfterX =
        TFHEpp::GLRawProductCiphertext<GLP, Schedule::input_log_q,
                                       Schedule::input_log_delta,
                                       Schedule::x_transform_log_scale>;
    auto input = std::make_unique<AfterX>();
    for (std::size_t component = 0; component < 2; component++)
        for (auto &slice : (*input)[component])
            for (auto &coefficient : slice) {
                coefficient.limb[0] = rng();
                coefficient.limb[1] = rng() & ((std::uint64_t{1} << 21) - 1);
            }
    const auto cache_start = std::chrono::steady_clock::now();
    TFHEpp::gl_detail::prepareSmallKeySwitchNTTCache<GLP, SwitchKey>(
        rotation_key->switch_key);
    const double cache_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      cache_start)
            .count();
    auto output = std::make_unique<AfterX>();
    auto body = std::make_unique<TFHEpp::GLPolynomial<GLP>>();
    auto mask = std::make_unique<TFHEpp::GLPolynomial<GLP>>();
    const auto start = std::chrono::steady_clock::now();
    TFHEpp::gl_detail::polynomialAutomorphism<GLP>(*body, (*input)[0], 1, 1, 1,
                                                   rotation_key->multiplier);
    TFHEpp::gl_detail::polynomialAutomorphism<GLP>(*mask, (*input)[1], 1, 1, 1,
                                                   rotation_key->multiplier);
    const auto switch_start = std::chrono::steady_clock::now();
    TFHEpp::GLDDSmallKeySwitch(output->ct, *mask, rotation_key->switch_key);
    const auto add_start = std::chrono::steady_clock::now();
    TFHEpp::gl_detail::addInPlace<GLP>((*output)[0], *body);
    TFHEpp::gl_detail::reduce<GLP, Schedule::input_log_q>((*output)[0]);
    TFHEpp::gl_detail::reduce<GLP, Schedule::input_log_q>((*output)[1]);
    const auto finish = std::chrono::steady_clock::now();
    const double seconds =
        std::chrono::duration<double>(finish - start).count();
    const double automorphism_seconds =
        std::chrono::duration<double>(switch_start - start).count();
    const double switch_seconds =
        std::chrono::duration<double>(add_start - switch_start).count();
    const double add_seconds =
        std::chrono::duration<double>(finish - add_start).count();
    std::cout << "n512 warm batch rotation: " << seconds << " s ("
              << automorphism_seconds << " s automorphisms, " << switch_seconds
              << " s DD switch, " << add_seconds << " s add/reduce, "
              << cache_seconds << " s cache preparation)" << std::endl;
    return (*output)[0][0][0] != typename GLP::T{} ||
           (*output)[1][0][0] != typename GLP::T{};
}

bool benchmarkProductionStC()
{
    if (std::getenv("TFHEPP_GL_N512_STC_BENCH") == nullptr) return true;
    using Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;
    using GLP = typename Schedule::Parameter;
    using StCKey = TFHEpp::GLSHIPSlotsToCoefficientsKey<Schedule>;
    std::mt19937_64 rng(0x4e35313253544342ULL);
    auto key = std::make_unique<StCKey>();
    key->conjugate_transpose_key.allocate();
    for (auto &ciphertext : key->conjugate_transpose_key.data)
        for (std::size_t component = 0; component < 2; component++)
            for (std::size_t y = 0; y < GLP::matrix_dimension; y++)
                for (auto &coefficient : ciphertext[component][y])
                    coefficient = static_cast<std::int16_t>(rng());

    constexpr std::array<std::uint32_t, 6> rotation_amounts{1, 2, 3, 4, 8, 12};
    key->w_rotation_keys.reserve(rotation_amounts.size());
    for (const std::uint32_t amount : rotation_amounts) {
        typename StCKey::RotationEntry entry;
        entry.amount = amount;
        entry.key.amount = amount;
        entry.key.multiplier = TFHEpp::gl_detail::powMod(
            GLP::primitive_root, amount, GLP::cyclotomic_order);
        entry.key.switch_key.allocate();
        for (auto &ciphertext : entry.key.switch_key.data)
            for (auto &component : ciphertext)
                for (auto &coefficient : component)
                    coefficient = static_cast<std::int16_t>(rng());
        key->w_rotation_keys.push_back(std::move(entry));
    }

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    for (std::size_t component = 0; component < 2; component++)
        for (auto &slice : (*input)[component])
            for (auto &coefficient : slice) {
                coefficient.limb[0] = rng();
                coefficient.limb[1] = rng() & ((std::uint64_t{1} << 21) - 1);
            }
    auto output = std::make_unique<typename Schedule::CoefficientCiphertext>();
    const auto start = std::chrono::steady_clock::now();
    TFHEpp::GLSHIPSlotsToCoefficients<Schedule>(*output, *input, *key);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cout << "n512 end-to-end SlotsToCoefficients: " << seconds << " s"
              << std::endl;
    return (*output)[0][0][0] != typename GLP::T{} ||
           (*output)[1][0][0] != typename GLP::T{};
}

bool benchmarkProductionMaskedColumn()
{
    if (std::getenv("TFHEPP_GL_N512_MASKED_BENCH") == nullptr) return true;
    using Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;
    using GLP = typename Schedule::Parameter;
    using P = typename GLP::baseP;
    constexpr std::size_t candidate_count = 48;
    std::mt19937_64 rng(0x4e3531324d41534bULL);
    auto key = std::make_unique<TFHEpp::GLSHIPMaskedColumnKey<Schedule>>();
    key->candidates.reserve(candidate_count);
    key->encrypted_masks.resize(candidate_count);
    for (std::size_t candidate = 0; candidate < candidate_count; candidate++) {
        key->candidates.push_back(
            {static_cast<std::uint32_t>(candidate % Schedule::theta),
             static_cast<std::uint32_t>((candidate / Schedule::theta) %
                                        GLP::phi),
             static_cast<std::uint32_t>(candidate % 4)});
        for (auto &component : key->encrypted_masks[candidate])
            for (auto &coefficient : component) {
                for (auto &limb : coefficient.limb) limb = rng();
                coefficient = TFHEpp::ckks_detail::reduceToLevel<
                    P, TFHEpp::GLSHIPMaskedColumnKey<Schedule>::key_log_q>(
                    coefficient);
            }
    }
    auto mask = std::make_unique<TFHEpp::GLBasePolynomial<GLP>>();
    for (auto &coefficient : *mask) {
        coefficient.limb[0] = rng();
        coefficient = TFHEpp::ckks_detail::reduceToLevel<P, Schedule::q0_log_q>(
            coefficient);
    }
    using Output = TFHEpp::GLBaseCiphertext<GLP, Schedule::half_bootstrap_log_q,
                                            Schedule::tree_log_delta>;
    auto output = std::make_unique<Output>();
    TFHEpp::GLBasePlaintext<GLP,
                            TFHEpp::GLSHIPMaskedColumnKey<Schedule>::key_log_q,
                            Schedule::tree_log_delta>
        sample_plaintext;
    TFHEpp::gl_ship_detail::buildCandidatePlaintext<Schedule>(
        sample_plaintext, *mask, key->candidates.front(), 0);
    std::uint32_t sample_bits = 0;
    for (const auto &coefficient : sample_plaintext.poly) {
        const auto [negative, magnitude] =
            TFHEpp::ckks_detail::smallSignedMagnitude<
                P, TFHEpp::GLSHIPMaskedColumnKey<Schedule>::key_log_q>(
                coefficient);
        const std::int64_t signed_coefficient =
            negative ? -static_cast<std::int64_t>(magnitude)
                     : static_cast<std::int64_t>(magnitude);
        sample_bits =
            std::max(sample_bits,
                     TFHEpp::gl_detail::signedTorusBitWidth(
                         TFHEpp::gl_detail::signedI128ToTorus<typename GLP::T>(
                             signed_coefficient)));
    }
    const auto cache_start = std::chrono::steady_clock::now();
    const auto cache =
        TFHEpp::gl_ship_detail::prepareMaskedColumnNTTCache<Schedule>(*key);
    const double cache_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      cache_start)
            .count();
    const double cache_mib =
        (cache->spectra[0].size() + cache->spectra[1].size()) *
        sizeof(std::uint64_t) / (1024.0 * 1024.0);
    const auto start = std::chrono::steady_clock::now();
    TFHEpp::gl_ship_detail::maskedColumn<Schedule>(*output, *mask, 0, *key);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();

    std::size_t batch_count = 1;
#ifdef _OPENMP
    batch_count = 4 * static_cast<std::size_t>(omp_get_max_threads());
#endif
    if (std::getenv("TFHEPP_GL_N512_FACTOR_BENCH") != nullptr)
        batch_count *= 2;
    std::vector<Output> output_batch(batch_count);
    const auto batch_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::gl_ship_detail::MaskedColumnNTTWorkspace<Schedule> workspace;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            TFHEpp::gl_ship_detail::maskedColumn<Schedule>(
                output_batch[batch], *mask, 0, *key, &workspace);
    }
    const double batch_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      batch_start)
            .count();
    for (const auto &batch : output_batch)
        if (batch.ct != output->ct) return false;
    std::cout << "n512 masked column (" << candidate_count << " candidates, "
              << sample_bits << "-bit encoded coefficients): " << seconds
              << " s" << " (" << cache_seconds << " s cache preparation, "
              << cache_mib << " MiB cache); " << batch_count
              << " warm calls with reused per-thread scratch in "
              << batch_seconds << " s" << std::endl;
    if (std::getenv("TFHEPP_GL_N512_FACTOR_BENCH") != nullptr) {
        using SwitchKey = TFHEpp::GLDDSmallKeySwitchKey<
            GLP, Schedule::half_bootstrap_log_q, Schedule::primary_bit,
            Schedule::bbar_bit>;
        TFHEpp::GLSHIPHMuxKey<Schedule> hmux_key;
        hmux_key.stages.resize(1);
        auto &stage = hmux_key.stages.front();
        stage.step = 1;
        stage.branches.resize(Schedule::hmux_radix);
        std::vector<SwitchKey *> hmux_switches;
        for (auto &branch : stage.branches) {
            for (SwitchKey *switch_key :
                 {&branch.body_key, &branch.mask_key}) {
                switch_key->allocate();
                for (auto &ciphertext : switch_key->data)
                    for (auto &component : ciphertext)
                        for (auto &coefficient : component)
                            coefficient = static_cast<std::int16_t>(rng());
                hmux_switches.push_back(switch_key);
            }
        }
#pragma omp parallel for schedule(dynamic, 1)
        for (std::size_t i = 0; i < hmux_switches.size(); i++)
            TFHEpp::gl_detail::prepareSmallKeySwitchNTTCache<GLP, SwitchKey>(
                *hmux_switches[i]);
#ifdef USE_HEXL
#pragma omp parallel for schedule(dynamic, 1)
        for (std::size_t i = 0; i < hmux_switches.size(); i++)
            TFHEpp::gl_detail::prepareSmallKeySwitchBlockedNTTCache<
                GLP, SwitchKey>(*hmux_switches[i]);
#endif

        std::vector<Output> fused_outputs(batch_count);
        const auto fused_start = std::chrono::steady_clock::now();
#pragma omp parallel
        {
            TFHEpp::gl_ship_detail::MaskedColumnNTTWorkspace<Schedule>
                masked_workspace;
            TFHEpp::gl_ship_detail::HMuxNTTWorkspace<Schedule> hmux_workspace;
            auto selected = std::make_unique<Output>();
#pragma omp for schedule(dynamic)
            for (std::size_t batch = 0; batch < batch_count; batch++) {
                TFHEpp::gl_ship_detail::maskedColumn<Schedule>(
                    *selected, *mask, static_cast<std::uint32_t>(batch & 1U),
                    *key, &masked_workspace);
                TFHEpp::gl_ship_detail::hmux<Schedule>(
                    fused_outputs[batch], *selected, hmux_key,
                    &hmux_workspace);
            }
        }
        const double fused_factor_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          fused_start)
                .count();

        int maximum_workers = 1;
#ifdef _OPENMP
        maximum_workers = omp_get_max_threads();
#endif
        const int hmux_workers = std::max(1, maximum_workers / 2);
        for (const std::size_t requested_tile : {256U, 128U, 64U}) {
#ifdef USE_HEXL
            constexpr std::array blocked_layouts{false, true};
#else
            constexpr std::array blocked_layouts{false};
#endif
            for (const bool blocked_layout : blocked_layouts) {
                const std::size_t tile_capacity =
                    std::min(requested_tile, batch_count);
                std::vector<Output> selected_tile(tile_capacity);
                std::vector<Output> staged_outputs(batch_count);
                const auto staged_start = std::chrono::steady_clock::now();
                for (std::size_t tile_begin = 0; tile_begin < batch_count;
                     tile_begin += tile_capacity) {
                    const std::size_t tile_count =
                        std::min(tile_capacity, batch_count - tile_begin);
#pragma omp parallel
                    {
                        TFHEpp::gl_ship_detail::MaskedColumnNTTWorkspace<
                            Schedule>
                            workspace;
#pragma omp for schedule(dynamic)
                        for (std::size_t local = 0; local < tile_count;
                             local++)
                            TFHEpp::gl_ship_detail::maskedColumn<Schedule>(
                                selected_tile[local], *mask,
                                static_cast<std::uint32_t>(
                                    (tile_begin + local) & 1U),
                                *key, &workspace);
                    }
#pragma omp parallel num_threads(hmux_workers)
                    {
                        TFHEpp::gl_ship_detail::HMuxNTTWorkspace<Schedule>
                            workspace;
#ifdef USE_HEXL
                        workspace.switch_workspace.use_batched_vector_mac =
                            true;
                        workspace.switch_workspace
                            .use_coefficient_blocked_key_layout =
                            blocked_layout;
#endif
#pragma omp for schedule(dynamic)
                        for (std::size_t local = 0; local < tile_count;
                             local++)
                            TFHEpp::gl_ship_detail::hmux<Schedule>(
                                staged_outputs[tile_begin + local],
                                selected_tile[local], hmux_key, &workspace);
                    }
                }
                const double staged_seconds =
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - staged_start)
                        .count();
                for (std::size_t batch = 0; batch < batch_count; batch++)
                    if (staged_outputs[batch].ct != fused_outputs[batch].ct)
                        return false;
                std::cout
                    << "n512 one-stage sparse factors, tile " << requested_tile
                    << (blocked_layout ? " blocked" : " row-major") << ": "
                    << staged_seconds << " s (" << maximum_workers
                    << " masked / " << hmux_workers << " HMux workers) vs "
                    << fused_factor_seconds << " s fused" << std::endl;
            }
        }
    }
    return (*output)[0][0] != typename GLP::T{} ||
           (*output)[1][0] != typename GLP::T{};
}

template <std::size_t Level>
bool benchmarkProductionProductLevel(double &projected_total_seconds)
{
    using Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;
    using GLP = typename Schedule::Parameter;
    using P = typename GLP::baseP;
    using RelinKey = std::tuple_element_t<
        Level, typename TFHEpp::GLSHIPProductRelinKeyChain<Schedule>::Tuple>;
    constexpr std::uint32_t input_log_q =
        Schedule::half_bootstrap_log_q -
        static_cast<std::uint32_t>(Level) * Schedule::tree_log_delta;
    constexpr std::uint32_t output_log_q =
        input_log_q - Schedule::tree_log_delta;
    static_assert(RelinKey::log_q == input_log_q);
    using Input =
        TFHEpp::GLBaseCiphertext<GLP, input_log_q, Schedule::tree_log_delta>;
    using Output =
        TFHEpp::GLBaseCiphertext<GLP, output_log_q, Schedule::tree_log_delta>;

    std::mt19937_64 rng(0x50524f4454524545ULL + Level);
    auto key = std::make_unique<RelinKey>();
    key->allocate();
    for (auto &ciphertext : key->data)
        for (auto &component : ciphertext)
            for (auto &coefficient : component)
                coefficient = static_cast<std::int16_t>(rng());
    TFHEpp::gl_detail::prepareSmallKeySwitchNTTCache<GLP, RelinKey>(*key);

    auto lhs = std::make_unique<Input>();
    auto rhs = std::make_unique<Input>();
    for (std::size_t component = 0; component < 2; component++) {
        for (auto &coefficient : (*lhs)[component]) {
            for (auto &limb : coefficient.limb) limb = rng();
            coefficient =
                TFHEpp::ckks_detail::reduceToLevel<P, input_log_q>(coefficient);
        }
        for (auto &coefficient : (*rhs)[component]) {
            for (auto &limb : coefficient.limb) limb = rng();
            coefficient =
                TFHEpp::ckks_detail::reduceToLevel<P, input_log_q>(coefficient);
        }
    }

    auto expected = std::make_unique<Output>();
    TFHEpp::GLBaseHadamardWorkspace<GLP> single_workspace;
    TFHEpp::GLBaseHadamardMultiply<GLP, input_log_q, Schedule::tree_log_delta,
                                   input_log_q, Schedule::tree_log_delta,
                                   Schedule::tree_log_delta,
                                   Schedule::primary_bit, Schedule::bbar_bit>(
        *expected, *lhs, *rhs, *key, &single_workspace);

    std::size_t batch_count = 1;
#ifdef _OPENMP
    batch_count = 4 * static_cast<std::size_t>(omp_get_max_threads());
#endif
    std::vector<Output> products(batch_count);
    const auto start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::GLBaseHadamardWorkspace<GLP> workspace;
        TFHEpp::gl_detail::SmallKeySwitchSumNTTWorkspace<GLP, RelinKey, 1>
            switch_workspace;
#ifdef USE_HEXL
        switch_workspace.use_batched_vector_mac = true;
#endif
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            TFHEpp::GLBaseHadamardMultiply<
                GLP, input_log_q, Schedule::tree_log_delta, input_log_q,
                Schedule::tree_log_delta, Schedule::tree_log_delta,
                Schedule::primary_bit, Schedule::bbar_bit>(
                products[batch], *lhs, *rhs, *key, &workspace,
                &switch_workspace);
    }
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    for (const auto &product : products)
        if (product.ct != expected->ct) return false;

    std::vector<Output> scalar_relin_products(batch_count);
    const auto scalar_relin_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        TFHEpp::GLBaseHadamardWorkspace<GLP> workspace;
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++)
            TFHEpp::GLBaseHadamardMultiply<
                GLP, input_log_q, Schedule::tree_log_delta, input_log_q,
                Schedule::tree_log_delta, Schedule::tree_log_delta,
                Schedule::primary_bit, Schedule::bbar_bit>(
                scalar_relin_products[batch], *lhs, *rhs, *key, &workspace);
    }
    const double scalar_relin_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      scalar_relin_start)
            .count();
    for (const auto &product : scalar_relin_products)
        if (product.ct != expected->ct) return false;

    std::vector<Output> legacy_products(batch_count);
    const auto legacy_start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        auto tensor =
            std::make_unique<std::array<TFHEpp::GLBasePolynomial<GLP>, 4>>();
        auto square_term =
            std::make_unique<TFHEpp::GLBaseCiphertextData<GLP>>();
        auto relinearized =
            std::make_unique<TFHEpp::GLBaseCiphertextData<GLP>>();
#pragma omp for schedule(static)
        for (std::size_t batch = 0; batch < batch_count; batch++) {
            TFHEpp::gl_detail::baseMultiplyAtLevel<GLP, input_log_q>(
                (*tensor)[0], (*lhs)[0], (*rhs)[0]);
            TFHEpp::gl_detail::baseMultiplyAtLevel<GLP, input_log_q>(
                (*tensor)[1], (*lhs)[0], (*rhs)[1]);
            TFHEpp::gl_detail::baseMultiplyAtLevel<GLP, input_log_q>(
                (*tensor)[2], (*lhs)[1], (*rhs)[0]);
            TFHEpp::gl_detail::baseMultiplyAtLevel<GLP, input_log_q>(
                (*tensor)[3], (*lhs)[1], (*rhs)[1]);
            TFHEpp::gl_detail::addInPlace<GLP>((*tensor)[1], (*tensor)[2]);
            TFHEpp::gl_detail::reduce<GLP, input_log_q>((*tensor)[1]);
            TFHEpp::GLDDSmallKeySwitchBase(*square_term, (*tensor)[3], *key);
            (*relinearized)[0] = (*tensor)[0];
            (*relinearized)[1] = (*tensor)[1];
            TFHEpp::gl_detail::addInPlace<GLP>((*relinearized)[0],
                                               (*square_term)[0]);
            TFHEpp::gl_detail::addInPlace<GLP>((*relinearized)[1],
                                               (*square_term)[1]);
            TFHEpp::gl_ship_detail::reduce<GLP, input_log_q>(*relinearized);
            TFHEpp::gl_ship_detail::rescaleBase<GLP, input_log_q,
                                                Schedule::tree_log_delta>(
                legacy_products[batch].ct, *relinearized);
            TFHEpp::gl_ship_detail::reduce<GLP, output_log_q>(
                legacy_products[batch].ct);
        }
    }
    const double legacy_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      legacy_start)
            .count();
    for (const auto &product : legacy_products)
        if (product.ct != expected->ct) return false;

    constexpr std::uint64_t multiplication_count =
        static_cast<std::uint64_t>(2) * GLP::matrix_dimension *
        (Schedule::factor_count >> (Level + 1));
    const double projected_seconds =
        seconds * multiplication_count / batch_count;
    projected_total_seconds += projected_seconds;
    std::cout << "n512 product level " << Level << " (logQ=" << input_log_q
              << "): " << batch_count << " calls in " << seconds << " s vs "
              << scalar_relin_seconds << " s with scalar relinearization and "
              << legacy_seconds
              << " s with separate tensor transforms; projected "
              << projected_seconds << " s" << std::endl;
    return true;
}

template <std::size_t... Levels>
bool benchmarkProductionProductTreeImpl(std::index_sequence<Levels...>)
{
    double projected_total_seconds = 0;
    const bool valid =
        (benchmarkProductionProductLevel<Levels>(projected_total_seconds) &&
         ...);
    std::cout << "n512 projected complete product tree: "
              << projected_total_seconds << " s" << std::endl;
    return valid;
}

bool benchmarkProductionProductTree()
{
    if (std::getenv("TFHEPP_GL_N512_PRODUCT_BENCH") == nullptr) return true;
    using Schedule = TFHEpp::GLSHIP512p17FusedDDSchedule;
    return benchmarkProductionProductTreeImpl(
        std::make_index_sequence<Schedule::tree_depth>{});
}

bool benchmarkProductionBigSwitch()
{
    if (std::getenv("TFHEPP_GL_N512_BIG_BENCH") == nullptr) return true;
    using GLP = TFHEpp::GL512p17Parameter;
    using P = typename GLP::baseP;
    using SwitchKey = TFHEpp::GLDDBigKeySwitchKey<GLP, 85, 85, 16>;
    std::mt19937_64 rng(0x4e35313242494744ULL);
    auto switch_key = std::make_unique<SwitchKey>();
    switch_key->allocate();
    for (auto &ciphertext : switch_key->data)
        for (std::size_t component = 0; component < 2; component++)
            for (std::size_t y = 0; y < ciphertext[component].size(); y++)
                for (auto &coefficient : ciphertext[component][y])
                    coefficient = static_cast<std::int16_t>(rng());

    auto input = std::make_unique<TFHEpp::GLPolynomial<GLP>>();
    for (auto &slice : *input)
        for (auto &coefficient : slice) {
            coefficient.limb[0] = rng();
            coefficient.limb[1] = rng();
            coefficient =
                TFHEpp::ckks_detail::reduceToLevel<P, 85>(coefficient);
        }
    auto output = std::make_unique<TFHEpp::GLCiphertextData<GLP>>();
    const auto start = std::chrono::steady_clock::now();
    TFHEpp::GLDDBigKeySwitch(*output, *input, *switch_key);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cout << "n512 dense 85x16 full DD big switch: " << seconds << " s"
              << std::endl;

    using Cache = TFHEpp::gl_detail::BigKeySwitchNTTCache<GLP, SwitchKey>;
    auto cache = std::make_unique<Cache>();
    const auto cache_start = std::chrono::steady_clock::now();
    TFHEpp::gl_detail::prepareBigKeySwitchNTTCache<GLP, SwitchKey>(*cache,
                                                                   *switch_key);
    const double cache_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      cache_start)
            .count();
    auto cached_output = std::make_unique<TFHEpp::GLCiphertextData<GLP>>();
    const auto cached_start = std::chrono::steady_clock::now();
    TFHEpp::GLDDBigKeySwitch(*cached_output, *input, *switch_key, cache.get());
    const double cached_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      cached_start)
            .count();
    std::cout << "n512 cached full DD big switch: " << cached_seconds << " s ("
              << cache_seconds << " s cache preparation, "
              << TFHEpp::gl_detail::bigKeySwitchNTTCacheBytes<GLP, SwitchKey> /
                     (1024.0 * 1024.0 * 1024.0)
              << " GiB cache)" << std::endl;
    for (std::size_t component = 0; component < 2; component++)
        for (std::size_t y = 0; y < GLP::matrix_dimension; y++)
            if ((*cached_output)[component][y] != (*output)[component][y])
                return false;
    return (*output)[0][0][0] != typename GLP::T{} ||
           (*output)[1][0][0] != typename GLP::T{};
}

bool checkProductionShape()
{
    using GLP = TFHEpp::GL512p17Parameter;
    using P = typename GLP::baseP;
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

    for (auto &coefficient : *lhs) {
        for (auto &limb : coefficient.limb) limb = rng();
        coefficient = TFHEpp::ckks_detail::reduceToLevel<P, 338>(coefficient);
    }
    for (auto &coefficient : *rhs) {
        for (auto &limb : coefficient.limb) limb = rng();
        coefficient = TFHEpp::ckks_detail::reduceToLevel<P, 338>(coefficient);
    }
    const auto double_chunk_start = std::chrono::steady_clock::now();
    if (!TFHEpp::gl_detail::baseMultiplyNTT<GLP>(*got, *lhs, *rhs))
        return false;
    const double double_chunk_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      double_chunk_start)
            .count();
    std::cout << "n512 p17 dense 338x338 base multiply: "
              << double_chunk_seconds << " s" << std::endl;
    return true;
}

}  // namespace

int main()
{
    if (!checkDirectBaseAutomorphism()) {
        std::cerr << "direct GL base automorphism regression failed"
                  << std::endl;
        return 1;
    }
    if (!checkSymmetricDecompositionSpectrumHoist()) {
        std::cerr << "GL symmetric decomposition/spectrum hoist failed"
                  << std::endl;
        return 1;
    }
    if (!checkFusedCiphertextTensorMultiply()) {
        std::cerr << "fused GL ciphertext tensor multiply failed" << std::endl;
        return 1;
    }
    if (!checkSmallDense()) {
        std::cerr << "small GL NTT multiplication regression failed"
                  << std::endl;
        return 1;
    }
    if (!checkFullPolynomialDense()) {
        std::cerr << "full GL polynomial NTT regression failed" << std::endl;
        return 1;
    }
    if (!checkDDSwitchAccumulation()) {
        std::cerr << "GL DD transform accumulation regression failed"
                  << std::endl;
        return 1;
    }
    if (!checkBigDDSwitchAccumulation()) {
        std::cerr << "GL big DD transform accumulation regression failed"
                  << std::endl;
        return 1;
    }
    if (!checkOnePrimeDDSwitchAccumulation()) {
        std::cerr << "one-prime GL DD transform accumulation regression failed"
                  << std::endl;
        return 1;
    }
    if (!checkFusedWTransform()) {
        std::cerr << "fused GL W-transform regression failed" << std::endl;
        return 1;
    }
    if (!checkFusedMaskedColumn()) {
        std::cerr << "fused GL masked-column regression failed" << std::endl;
        return 1;
    }
    if (!benchmarkProductionDDBaseSwitch()) {
        std::cerr << "n512 GL DD base-switch benchmark produced zero output"
                  << std::endl;
        return 1;
    }
    if (!benchmarkProductionFullTransform()) {
        std::cerr << "n512 full-ring NTT roundtrip failed" << std::endl;
        return 1;
    }
    if (!benchmarkProductionHMuxSwitchFusion()) {
        std::cerr << "n512 HMux DD fusion regression failed" << std::endl;
        return 1;
    }
    if (!benchmarkProductionXTrace()) {
        std::cerr << "n512 compact X trace produced zero output" << std::endl;
        return 1;
    }
    if (!benchmarkProductionWTransform()) {
        std::cerr << "n512 fused W transform produced zero output" << std::endl;
        return 1;
    }
    if (!benchmarkProductionBatchRotation()) {
        std::cerr << "n512 batch rotation produced zero output" << std::endl;
        return 1;
    }
    if (!benchmarkProductionStC()) {
        std::cerr << "n512 StC produced zero output" << std::endl;
        return 1;
    }
    if (!benchmarkProductionMaskedColumn()) {
        std::cerr << "n512 masked column produced zero output" << std::endl;
        return 1;
    }
    if (!benchmarkProductionProductTree()) {
        std::cerr << "n512 product-tree regression failed" << std::endl;
        return 1;
    }
    if (!benchmarkProductionBigSwitch()) {
        std::cerr << "n512 big DD benchmark produced zero output" << std::endl;
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
