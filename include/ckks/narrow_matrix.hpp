#pragma once

#include <modular_ntt.hpp>
#include <algorithm>
#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace TFHEpp::narrow_detail {

// Lemma 5.1: multiplying odd NTT exponents permutes natural-order points.
// Ring degrees are powers of two, so only the low bits of the product matter.
// Keep the selector exponent instead of a degree-sized index table per row.
template <std::size_t N>
inline std::size_t narrowNTTPermutationIndex(std::size_t exponent, std::size_t point)
{
    static_assert(N > 0 && (N & (N - 1)) == 0);
    return ((exponent * (2 * point + 1)) & (2 * N - 1)) >> 1;
}

// Generate Shoup constants for many canonical selectors with one invariant
// divisor. With R=floor(2^128/p), floor(s*R/2^64) underestimates
// floor(s*2^64/p) by at most one for s<p<2^63. The remainder corrects it.
class SelectorMultiplierFactory {
    std::uint64_t modulus_;
    unsigned __int128 reciprocal_ = 0;
public:
    explicit SelectorMultiplierFactory(std::uint64_t modulus) : modulus_(modulus)
    {
        if (modulus == 0 || modulus >= (std::uint64_t{1} << 63))
            throw std::invalid_argument("invalid selector multiplier modulus");
        if (modulus > 1) {
            constexpr auto maximum = ~static_cast<unsigned __int128>(0);
            reciprocal_ = maximum / modulus + (maximum % modulus == modulus - 1);
        }
    }
    modular_ntt::ShoupMultiplier make(std::uint64_t selector) const
    {
        if (selector >= modulus_) throw std::invalid_argument("noncanonical selector multiplier");
        modular_ntt::ShoupMultiplier result;
        result.value = selector;
        const auto estimate = static_cast<std::uint64_t>((selector * reciprocal_) >> 64);
        const auto remainder = (static_cast<unsigned __int128>(selector) << 64)
            - static_cast<unsigned __int128>(estimate) * modulus_;
        result.quotient = estimate + (remainder >= modulus_);
        return result;
    }
};

// Both signed tables use the same selector spectrum. Share the quotient
// calculation while retaining a separate canonical modular sum for each table.
// All inputs must be canonical residues; the two destination rows are distinct.
inline void addSelectorTablePair(std::uint64_t *first, std::uint64_t *second,
                                 const std::uint64_t *first_table,
                                 const std::uint64_t *second_table,
                                 std::size_t columns, std::uint64_t selector,
                                 std::uint64_t modulus)
{
    const modular_ntt::ShoupMultiplier fixed(selector, modulus);
    for (std::size_t column = 0; column < columns; column++) {
        first[column] = modular_ntt::add(first[column], fixed.apply(first_table[column], modulus), modulus);
        second[column] = modular_ntt::add(second[column], fixed.apply(second_table[column], modulus), modulus);
    }
}

#if defined(__AVX2__)
namespace matrix_vector_detail {

// AVX2 has unsigned 32x32 multiplication but no 64x64 high product.
inline __m256i highProduct(__m256i a, __m256i b)
{
    const auto mask = _mm256_set1_epi64x(0xffffffffULL);
    const auto ah = _mm256_srli_epi64(a, 32), bh = _mm256_srli_epi64(b, 32);
    const auto ll = _mm256_mul_epu32(a, b);
    const auto lh = _mm256_mul_epu32(a, bh);
    const auto hl = _mm256_mul_epu32(ah, b);
    const auto hh = _mm256_mul_epu32(ah, bh);
    const auto middle = _mm256_add_epi64(_mm256_srli_epi64(ll, 32),
        _mm256_add_epi64(_mm256_and_si256(lh, mask), _mm256_and_si256(hl, mask)));
    return _mm256_add_epi64(hh, _mm256_add_epi64(_mm256_srli_epi64(lh, 32),
        _mm256_add_epi64(_mm256_srli_epi64(hl, 32), _mm256_srli_epi64(middle, 32))));
}

inline __m256i lowProduct(__m256i a, __m256i b)
{
    const auto cross = _mm256_add_epi64(_mm256_mul_epu32(a, _mm256_srli_epi64(b, 32)),
                                       _mm256_mul_epu32(_mm256_srli_epi64(a, 32), b));
    return _mm256_add_epi64(_mm256_mul_epu32(a, b), _mm256_slli_epi64(cross, 32));
}

inline __m256i canonicalize(__m256i a, __m256i modulus)
{
    const auto reduce = _mm256_cmpgt_epi64(a, _mm256_sub_epi64(modulus, _mm256_set1_epi64x(1)));
    return _mm256_sub_epi64(a, _mm256_and_si256(reduce, modulus));
}
}  // namespace matrix_vector_detail
#endif

// Exact vector counterpart; keep the scalar entry point for paired controls.
// Signed lane comparisons are valid because both the lazy product and the
// sum of canonical residues stay below 2*modulus < 2^63.
inline void addSelectorTablePairPrepared(std::uint64_t *first, std::uint64_t *second,
                                      const std::uint64_t *first_table,
                                      const std::uint64_t *second_table,
                                      std::size_t columns, const modular_ntt::ShoupMultiplier &fixed,
                                      std::uint64_t modulus)
{
#if defined(__AVX2__)
    if (modulus < (std::uint64_t{1} << 62) && columns >= 4) {
        const auto value = _mm256_set1_epi64x(fixed.value);
        const auto quotient = _mm256_set1_epi64x(fixed.quotient);
        const auto prime = _mm256_set1_epi64x(modulus);
        const auto accumulate = [&](std::uint64_t *dest, const std::uint64_t *table) {
            using namespace matrix_vector_detail;
            const auto input = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(table));
            const auto estimate = highProduct(input, quotient);
            const auto lazy = _mm256_sub_epi64(lowProduct(input, value), lowProduct(estimate, prime));
            const auto product = canonicalize(lazy, prime);
            const auto sum = _mm256_add_epi64(product,
                _mm256_loadu_si256(reinterpret_cast<const __m256i *>(dest)));
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(dest), canonicalize(sum, prime));
        };
        std::size_t column = 0;
        for (; column + 4 <= columns; column += 4) {
            accumulate(first + column, first_table + column);
            accumulate(second + column, second_table + column);
        }
        for (; column < columns; ++column) {
            first[column] = modular_ntt::add(first[column], fixed.apply(first_table[column], modulus), modulus);
            second[column] = modular_ntt::add(second[column], fixed.apply(second_table[column], modulus), modulus);
        }
        return;
    }
#endif
    for (std::size_t column = 0; column < columns; column++) {
        first[column] = modular_ntt::add(first[column], fixed.apply(first_table[column], modulus), modulus);
        second[column] = modular_ntt::add(second[column], fixed.apply(second_table[column], modulus), modulus);
    }
}

inline void addSelectorTablePairVector(std::uint64_t *first, std::uint64_t *second,
                                      const std::uint64_t *first_table,
                                      const std::uint64_t *second_table,
                                      std::size_t columns, std::uint64_t selector,
                                      std::uint64_t modulus)
{
    const modular_ntt::ShoupMultiplier fixed(selector, modulus);
    addSelectorTablePairPrepared(first, second, first_table, second_table, columns, fixed, modulus);
}

struct SelectorMatrixRow {
    const std::uint64_t *spectrum;
    std::size_t exponent;
};

// For p < 2^61, retain each destination in [0,2p). A Shoup lazy
// product is also below 2p, so their sum is below 4p < 2^63.
// One reduction by 2p suffices, including signed AVX2 lane comparisons.
// Canonicalize once after all selector rows, before publishing the sums.
inline void addSelectorTablePairLazy(std::uint64_t *first, std::uint64_t *second,
    const std::uint64_t *first_table, const std::uint64_t *second_table,
    std::size_t columns, const modular_ntt::ShoupMultiplier &fixed, std::uint64_t modulus)
{
    const auto twice = 2 * modulus;
    std::size_t column = 0;
#if defined(__AVX2__)
    const auto value = _mm256_set1_epi64x(fixed.value);
    const auto quotient = _mm256_set1_epi64x(fixed.quotient);
    const auto prime = _mm256_set1_epi64x(modulus);
    const auto double_prime = _mm256_set1_epi64x(twice);
    const auto accumulate = [&](std::uint64_t *dest, const std::uint64_t *table) {
        using namespace matrix_vector_detail;
        const auto input = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(table));
        const auto estimate = highProduct(input, quotient);
        const auto product = _mm256_sub_epi64(lowProduct(input, value), lowProduct(estimate, prime));
        const auto sum = _mm256_add_epi64(product,
            _mm256_loadu_si256(reinterpret_cast<const __m256i *>(dest)));
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dest), canonicalize(sum, double_prime));
    };
    for (; column + 4 <= columns; column += 4) {
        accumulate(first + column, first_table + column);
        accumulate(second + column, second_table + column);
    }
#endif
    const auto accumulate_scalar = [&](std::uint64_t destination, std::uint64_t table) {
        const auto estimate = static_cast<std::uint64_t>(
            (static_cast<unsigned __int128>(table) * fixed.quotient) >> 64);
        const auto product = table * fixed.value - estimate * modulus;
        const auto sum = destination + product;
        return sum >= twice ? sum - twice : sum;
    };
    for (; column < columns; ++column) {
        first[column] = accumulate_scalar(first[column], first_table[column]);
        second[column] = accumulate_scalar(second[column], second_table[column]);
    }
}

// Retain a small destination block in cache while processing selector rows.
// Each point still accumulates selectors in their original order; there is no
// rounding or reconstruction at a block boundary.
template <std::size_t N, std::size_t PointTile = 256, bool LazyAccumulation = true, class GetRow>
inline void addSelectorMatrixRowsPair(std::uint64_t *first, std::uint64_t *second,
    const std::uint64_t *first_table, const std::uint64_t *second_table,
    std::size_t columns, std::size_t count, GetRow get_row,
    const SelectorMultiplierFactory &factory, std::uint64_t modulus)
{
    static_assert(PointTile > 0);
    const bool lazy = LazyAccumulation && modulus < (std::uint64_t{1} << 61);
    for (std::size_t begin = 0; begin < N; begin += PointTile) {
        const auto end = begin + std::min(PointTile, N - begin);
        for (std::size_t local = 0; local < count; ++local) {
            const auto row = get_row(local);
            for (std::size_t i = begin; i < end; ++i) {
                const auto permuted = narrowNTTPermutationIndex<N>(row.exponent, i);
                const auto fixed = factory.make(row.spectrum[i]);
                if (lazy)
                    addSelectorTablePairLazy(first + i * columns, second + i * columns,
                        first_table + permuted * columns, second_table + permuted * columns,
                        columns, fixed, modulus);
                else
                    addSelectorTablePairPrepared(first + i * columns, second + i * columns,
                        first_table + permuted * columns, second_table + permuted * columns,
                        columns, fixed, modulus);
            }
        }
        if (lazy)
            for (std::size_t i = begin * columns; i < end * columns; ++i) {
                if (first[i] >= modulus) first[i] -= modulus;
                if (second[i] >= modulus) second[i] -= modulus;
            }
    }
}

}  // namespace TFHEpp::narrow_detail
