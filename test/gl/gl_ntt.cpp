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
constexpr std::uint32_t dd_log_q = 96;
constexpr std::uint32_t dd_primary_bit = 48;
constexpr std::uint32_t dd_bbar_bit = 16;
using SmallDDSwitchKey =
    TFHEpp::GLDDSmallKeySwitchKey<SmallGLP, dd_log_q, dd_primary_bit,
                                  dd_bbar_bit>;
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
    const double per_slice_seconds =
        std::max(0.0, slices_seconds - cache_seconds) / measured_slices;
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
    return (*output)[0][0] != typename GLP::T{} ||
           (*output)[1][0] != typename GLP::T{};
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
    if (!checkDDSwitchAccumulation()) {
        std::cerr << "GL DD transform accumulation regression failed"
                  << std::endl;
        return 1;
    }
    if (!checkOnePrimeDDSwitchAccumulation()) {
        std::cerr << "one-prime GL DD transform accumulation regression failed"
                  << std::endl;
        return 1;
    }
    if (!benchmarkProductionDDBaseSwitch()) {
        std::cerr << "n512 GL DD base-switch benchmark produced zero output"
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
