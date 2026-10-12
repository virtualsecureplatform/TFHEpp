#pragma once

#include "ckks.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Public signed plaintext digits, reusable across ciphertexts and components.
// A negative coefficient contributes negative magnitude digits at every radix
// position. This folds the two unsigned products into one exact convolution.
template <class P, std::uint32_t LogQ>
struct CKKSPreparedSmallPlainFFT {
    static_assert(use_multilimb_digit_fft_v<P>);
    static constexpr int digit_bits = torus_unsigned_digit_bits_v<P>;
    static constexpr int digits = (64 + digit_bits - 1) / digit_bits;
    static_assert(2 * digit_bits + P::nbit + 3 < std::numeric_limits<double>::digits);
    std::array<PolynomialInFD<P>, digits> spectra{};
    std::array<bool, digits> nonzero{};
};

template <class P, std::uint32_t LogQ>
inline std::unique_ptr<CKKSPreparedSmallPlainFFT<P, LogQ>>
makeCKKSPreparedSmallPlainFFT(const Polynomial<P> &plain)
{
    using Cache = CKKSPreparedSmallPlainFFT<P, LogQ>;
    using T = typename P::T;
    constexpr std::uint64_t mask = (std::uint64_t{1} << Cache::digit_bits) - 1;
    auto result = std::make_unique<Cache>();
    auto digit = std::make_unique<Polynomial<P>>();
    for (int d = 0; d < Cache::digits; ++d) {
        for (std::size_t i = 0; i < P::n; ++i) {
            const auto [negative, magnitude] = ckks_detail::smallSignedMagnitude<P, LogQ>(plain[i]);
            const auto value = (magnitude >> (d * Cache::digit_bits)) & mask;
            (*digit)[i] = negative ? T{} - T{value} : T{value};
            result->nonzero[d] = result->nonzero[d] || value != 0;
        }
        if (result->nonzero[d]) TwistIFFTDigit<P>(result->spectra[d], *digit);
    }
    return result;
}

template <class P, std::uint32_t LogQ>
inline void CKKSPreparedSmallPlainMul(Polynomial<P> &result,
    const Polynomial<P> &torus, const CKKSPreparedSmallPlainFFT<P, LogQ> &plain)
{
    using Cache = CKKSPreparedSmallPlainFFT<P, LogQ>;
    using T = typename P::T;
    constexpr int width = std::numeric_limits<T>::digits;
    constexpr int bits = Cache::digit_bits;
    constexpr T half = T{1} << (bits - 1), mask = (T{1} << bits) - T{1};
    constexpr T offset = [] {
        T value{};
        for (int j = 0; j < width / bits; ++j)
            value += (T{1} << (bits - 1)) << (width - (j + 1) * bits);
        return value;
    }();
    if constexpr (P::n < 8192) {
        auto digit = std::make_unique<Polynomial<P>>();
        auto product = std::make_unique<Polynomial<P>>();
        auto spectrum = std::make_unique<PolynomialInFD<P>>();
        auto multiplied = std::make_unique<PolynomialInFD<P>>();
        result.fill(T{});
        for (int j = (width - LogQ) / bits; j < width / bits; ++j) {
            const int torus_shift = width - (j + 1) * bits;
            for (std::size_t i = 0; i < P::n; ++i)
                (*digit)[i] = (((torus[i] + offset) >> torus_shift) & mask) - half;
            TwistIFFTDigit<P>(*spectrum, *digit);
            for (int d = 0; d < Cache::digits; ++d) {
                const int shift = torus_shift + d * bits;
                if (shift >= LogQ || !plain.nonzero[d]) continue;
                MulInFD<P::n>(*multiplied, *spectrum, plain.spectra[d]);
                TwistFFTDigitProduct<P>(*product, *multiplied);
                for (std::size_t i = 0; i < P::n; ++i) result[i] += (*product)[i] << shift;
            }
        }
    }
    else {
        constexpr int first = (width - LogQ) / bits;
        constexpr int rows = width / bits - first;
        // Rows are independent until their exact torus sums are combined.
        std::vector<Polynomial<P>> row_products(rows);
    #pragma omp parallel if(P::n >= 8192)
        {
            auto digit = std::make_unique<Polynomial<P>>();
            auto product = std::make_unique<Polynomial<P>>();
            auto spectrum = std::make_unique<PolynomialInFD<P>>();
            auto multiplied = std::make_unique<PolynomialInFD<P>>();
    #pragma omp for schedule(static)
            for (int row = 0; row < rows; ++row) {
                const int j = first + row;
                const int torus_shift = width - (j + 1) * bits;
                auto &sum = row_products[row];
                for (std::size_t i = 0; i < P::n; ++i)
                    (*digit)[i] = (((torus[i] + offset) >> torus_shift) & mask) - half;
                TwistIFFTDigit<P>(*spectrum, *digit);
                for (int d = 0; d < Cache::digits; ++d) {
                    const int shift = torus_shift + d * bits;
                    if (shift >= LogQ || !plain.nonzero[d]) continue;
                    MulInFD<P::n>(*multiplied, *spectrum, plain.spectra[d]);
                    TwistFFTDigitProduct<P>(*product, *multiplied);
                    for (std::size_t i = 0; i < P::n; ++i) sum[i] += (*product)[i] << shift;
                }
            }
        }
        result.fill(T{});
        for (const auto &sum : row_products)
            for (std::size_t i = 0; i < P::n; ++i) result[i] += sum[i];
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSPreparedPlainMulRescale(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &result,
    const CKKSCiphertext<P, LogQ, LogDelta> &input,
    const CKKSPreparedSmallPlainFFT<P, LogQ> &plain)
{
    auto product = std::make_unique<Polynomial<P>>();
    for (std::size_t c = 0; c <= P::k; ++c) {
        CKKSPreparedSmallPlainMul<P, LogQ>(*product, input.ct[c], plain);
        for (std::size_t i = 0; i < P::n; ++i)
            result.ct[c][i] = ckks_detail::roundedLevelRightShift<P, LogQ, PlainLogDelta>((*product)[i]);
    }
    ckks_detail::reduceTRLWEToLevel<P, LogQ - PlainLogDelta>(result.ct);
}

} // namespace TFHEpp
#endif
