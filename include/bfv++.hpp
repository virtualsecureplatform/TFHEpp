#pragma once
#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>

#include "keyswitch.hpp"
#include "mulfft.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"

// #include "hexl/hexl.hpp"

namespace TFHEpp {

// Standard TRLWE multiplication without relinearization
// Uses PolyMulRescaleUnsigned for rescaling by Δ
template <class P>
void TRLWEMultWithoutRelinerization(TRLWE3<P> &res, const TRLWE<P> &a,
                                    const TRLWE<P> &b)
{
    alignas(64) PolynomialInFD<P> ffta, fftb, fftc;
    TwistIFFTUInt<P>(ffta, a[0]);
    TwistIFFTUInt<P>(fftb, b[1]);
    MulInFD<P::n>(fftc, ffta, fftb);
    TwistIFFTUInt<P>(ffta, a[1]);
    TwistIFFTUInt<P>(fftb, b[0]);
    FMAInFD<P::n>(fftc, ffta, fftb);
    TwistFFTrescale<P>(res[0], fftc);

    PolyMulRescaleUnsigned<P>(res[1], a[1], b[1]);
    PolyMulRescaleUnsigned<P>(res[2], a[0], b[0]);
}

// Relinearization key switch - automatically handles DD when l̅ > 1
template <class P>
inline void relinKeySwitch(TRLWE<P> &res, const Polynomial<P> &poly,
                           const relinKeyFFT<P> &relinkeyfft)
{
    if constexpr (is_multilimb_uint_v<typename P::T>) {
        static_assert(P::l̅ > 1,
                      "multi-limb relinKeySwitch currently expects DD rows");

        auto decvec = std::make_unique<DecomposedPolynomial<P>>();
        auto decvecfft = std::make_unique<PolynomialInFD<P>>();
        Decomposition<P>(*decvec, poly);

        auto resfft_dd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
        for (int j = 0; j < P::l̅; j++)
            for (int m = 0; m <= P::k; m++)
                for (int n = 0; n < static_cast<int>(P::n); n++)
                    (*resfft_dd)[j][m][n] = 0.0;

        for (int i = 0; i < static_cast<int>(P::l); i++) {
            TwistIFFTDigit<P>(*decvecfft, (*decvec)[i]);
            for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                const int row_idx = i * P::l̅ + j;
                for (int m = 0; m <= static_cast<int>(P::k); m++)
                    FMAInFD<P::n>((*resfft_dd)[j][m], *decvecfft,
                                  relinkeyfft[row_idx][m]);
            }
        }

        auto results_dd = std::make_unique<std::array<TRLWE<P>, P::l̅>>();
        for (int j = 0; j < static_cast<int>(P::l̅); j++)
            for (int k = 0; k <= static_cast<int>(P::k); k++)
                TwistFFTDigitProduct<P>((*results_dd)[j][k],
                                        (*resfft_dd)[j][k]);

        RecombineTRLWEFromDD<P, false>(res, *results_dd);
    }
    else if constexpr (P::l̅ > 1) {
        alignas(64) DecomposedPolynomial<P> decvec;
        Decomposition<P>(decvec, poly);
        alignas(64) PolynomialInFD<P> decvecfft;

        // Double Decomposition path: l̅ separate accumulators
        alignas(64) std::array<TRLWEInFD<P>, P::l̅> resfft_dd;

        // Initialize all accumulators to zero
        for (int j = 0; j < P::l̅; j++)
            for (int m = 0; m <= P::k; m++)
                for (int n = 0; n < P::n; n++)
                    resfft_dd[j][m][n] = 0.0;

        // Process with standard decomposition (l levels), accumulate into l̅ results
        for (int i = 0; i < P::l; i++) {
            TwistIFFT<P>(decvecfft, decvec[i]);
            // Each decomposition level i multiplies with l̅ relinkey rows
            for (int j = 0; j < P::l̅; j++) {
                const int row_idx = i * P::l̅ + j;
                for (int m = 0; m <= P::k; m++) {
                    FMAInFD<P::n>(resfft_dd[j][m], decvecfft,
                                  relinkeyfft[row_idx][m]);
                }
            }
        }

        // FFT back to coefficient domain for each accumulator and recombine
        std::array<TRLWE<P>, P::l̅> results_dd;
        for (int j = 0; j < P::l̅; j++)
            for (int k = 0; k <= P::k; k++)
                TwistFFT<P>(results_dd[j][k], resfft_dd[j][k]);

        // Recombine the l̅ TRLWEs back to single TRLWE
        RecombineTRLWEFromDD<P, false>(res, results_dd);
    }
    else {
        alignas(64) DecomposedPolynomial<P> decvec;
        Decomposition<P>(decvec, poly);
        alignas(64) PolynomialInFD<P> decvecfft;

        // Standard path
        TRLWEInFD<P> resfft;
        TwistIFFT<P>(decvecfft, decvec[0]);
        MulInFD<P::n>(resfft[0], decvecfft, relinkeyfft[0][0]);
        MulInFD<P::n>(resfft[1], decvecfft, relinkeyfft[0][1]);
        for (int i = 1; i < P::l; i++) {
            TwistIFFT<P>(decvecfft, decvec[i]);
            FMAInFD<P::n>(resfft[0], decvecfft, relinkeyfft[i][0]);
            FMAInFD<P::n>(resfft[1], decvecfft, relinkeyfft[i][1]);
        }
        TwistFFT<P>(res[0], resfft[0]);
        TwistFFT<P>(res[1], resfft[1]);
    }
}

// Relinearization - automatically handles DD when l̅ > 1
template <class P>
inline void Relinearization(TRLWE<P> &res, const TRLWE3<P> &mult,
                            const relinKeyFFT<P> &relinkeyfft)
{
    TRLWE<P> squareterm;
    relinKeySwitch<P>(squareterm, mult[2], relinkeyfft);
    for (int i = 0; i < P::n; i++) res[0][i] = mult[0][i] + squareterm[0][i];
    for (int i = 0; i < P::n; i++) res[1][i] = mult[1][i] + squareterm[1][i];
}

// TRLWE multiplication with relinearization - automatically handles DD when l̅ > 1
template <class P>
inline void TRLWEMult(TRLWE<P> &res, const TRLWE<P> &a, const TRLWE<P> &b,
                      const relinKeyFFT<P> &relinkeyfft)
{
    TRLWE3<P> resmult;
    TRLWEMultWithoutRelinerization<P>(resmult, a, b);
    Relinearization<P>(res, resmult, relinkeyfft);
}

// ---------------------------------------------------------------------------
// 384-bit unsigned accumulator for BFV FullDD multiplication.
// Stores a 6-limb (6 × 64-bit) unsigned value and supports wide shift-add.
// The raw product c_a·c_b (256-bit) is accumulated using the positional weights
// of the B̅g digits (up to 2^224), then divided by Δ = floor(Q/t) at the end.
// ---------------------------------------------------------------------------
struct Wide384 {
    uint64_t w[6] = {};  // w[0] = LSB, w[5] = MSB

    // Add a SIGNED 128-bit value left-shifted by `shift` bits.
    // Handles negative values correctly (sign-extends before shifting).
    void add_shifted(__int128_t svalue, int shift)
    {
        if (shift <= -128 || shift >= 384) return;
        if (shift < 0) {
            svalue >>= (-shift);
            shift = 0;
        }
        // For negative values: decompose into magnitude + sign, add/subtract
        const bool negative = (svalue < 0);
        __uint128_t value = negative ? static_cast<__uint128_t>(-svalue)
                                     : static_cast<__uint128_t>(svalue);

        const int limb_offset = shift / 64;
        const int bit_offset  = shift % 64;

        uint64_t v0 = static_cast<uint64_t>(value);
        uint64_t v1 = static_cast<uint64_t>(value >> 64);

        uint64_t s0, s1, s2;
        if (bit_offset == 0) {
            s0 = v0; s1 = v1; s2 = 0;
        } else {
            s0 = v0 << bit_offset;
            s1 = (v1 << bit_offset) | (v0 >> (64 - bit_offset));
            s2 = v1 >> (64 - bit_offset);
        }

        if (!negative) {
            // Add s0:s1:s2 with carry
            uint64_t carry = 0;
            for (int k = 0; k < 3 && limb_offset + k < 6; k++) {
                uint64_t sv = (k == 0) ? s0 : (k == 1) ? s1 : s2;
                __uint128_t sum = static_cast<__uint128_t>(w[limb_offset + k])
                                + sv + carry;
                w[limb_offset + k] = static_cast<uint64_t>(sum);
                carry = static_cast<uint64_t>(sum >> 64);
            }
            for (int k = 3; carry && limb_offset + k < 6; k++) {
                __uint128_t sum = static_cast<__uint128_t>(w[limb_offset + k]) + carry;
                w[limb_offset + k] = static_cast<uint64_t>(sum);
                carry = static_cast<uint64_t>(sum >> 64);
            }
        } else {
            // Subtract s0:s1:s2 with borrow
            uint64_t borrow = 0;
            for (int k = 0; k < 3 && limb_offset + k < 6; k++) {
                uint64_t sv = (k == 0) ? s0 : (k == 1) ? s1 : s2;
                uint64_t wk = w[limb_offset + k];
                // borrow detection: wk < sv + borrow (with care for overflow of sv+borrow)
                uint64_t new_borrow = (wk < sv) || (wk - sv < borrow) ? 1 : 0;
                w[limb_offset + k] = wk - sv - borrow;
                borrow = new_borrow;
            }
            for (int k = 3; borrow && limb_offset + k < 6; k++) {
                uint64_t wk = w[limb_offset + k];
                uint64_t new_borrow = (wk < borrow) ? 1 : 0;
                w[limb_offset + k] = wk - borrow;
                borrow = new_borrow;
            }
        }
    }

    __uint128_t div128_unsigned(__uint128_t divisor) const
    {
        // Schoolbook long division with 128-bit "digits":
        // Treat the 384-bit value as three 128-bit limbs: [L2:L1:L0]
        // where L0 = w[0]:w[1], L1 = w[2]:w[3], L2 = w[4]:w[5].
        __uint128_t L0 = static_cast<__uint128_t>(w[1]) << 64 | w[0];
        __uint128_t L1 = static_cast<__uint128_t>(w[3]) << 64 | w[2];
        __uint128_t L2 = static_cast<__uint128_t>(w[5]) << 64 | w[4];

        // Long division from MSB: quotient digits q2, q1, q0
        // Remainder fits in 128 bits at each step.
        __uint128_t rem = L2 % divisor;

        // Now divide (rem:L1) by divisor. rem < divisor, so this fits.
        // Use: (rem * 2^128 + L1) / divisor
        // Since rem < divisor < 2^128, rem*2^128 + L1 < 2*2^256 — needs care.
        // Use iterative approach: process 64 bits at a time.
        auto div256by128 = [](__uint128_t hi, __uint128_t lo,
                              __uint128_t d) -> std::pair<__uint128_t, __uint128_t> {
            // Divide (hi:lo) by d where hi < d. Returns (quotient, remainder).
            // Process bit by bit for the high part, then use hardware div for low.
            if (hi == 0) return {lo / d, lo % d};
            __uint128_t q = 0;
            __uint128_t r = hi;
            // Process 128 bits of lo from MSB to LSB
            for (int bit = 127; bit >= 0; bit--) {
                // r = r * 2 + next_bit
                bool overflow = (r >> 127) != 0;
                r = (r << 1) | (static_cast<__uint128_t>((lo >> bit) & 1));
                if (overflow || r >= d) {
                    r -= d;
                    q |= static_cast<__uint128_t>(1) << bit;
                }
            }
            return {q, r};
        };

        auto [q1, rem1] = div256by128(rem, L1, divisor);
        auto [q0, rem0] = div256by128(rem1, L0, divisor);

        // The full quotient is (q2 * 2^256 + q1 * 2^128 + q0).
        // We only need the low 128 bits for torus arithmetic.
        return q0;
    }

    // Divide this signed 384-bit two's-complement value by a 128-bit divisor,
    // return the quotient modulo 2^128.
    __uint128_t div128(__uint128_t divisor) const
    {
        const bool negative = (w[5] >> 63) != 0;
        if (!negative) return div128_unsigned(divisor);

        Wide384 magnitude;
        uint64_t carry = 1;
        for (int i = 0; i < 6; i++) {
            const uint64_t limb = ~w[i];
            const uint64_t next = limb + carry;
            magnitude.w[i] = next;
            carry = (carry != 0 && next == 0) ? 1 : 0;
        }

        const __uint128_t q = magnitude.div128_unsigned(divisor);
        return static_cast<__uint128_t>(-q);
    }
};

// Full Double Decomposition TRLWE Multiplication (without relinearization)
//
// For __uint128_t params (lvl3simdparam): uses 384-bit accumulators and
// divides by Δ=floor(Q/t) — required because the raw 256-bit product
// overflows 128-bit torus arithmetic.
//
// For smaller types (uint64_t): uses the original shift-based approach
// with power-of-2 Δ cancellation folded into the shift.
template <class P>
void TRLWEMultWithoutRelinearizationFullDD(TRLWE3<P> &res, const TRLWE<P> &a,
                                            const TRLWE<P> &b)
{
    constexpr int width = std::numeric_limits<typename P::T>::digits;

    if constexpr (std::is_same_v<typename P::T, __uint128_t>) {
        auto a_dec = std::make_unique<std::array<TRLWE<P>, P::l̅>>();
        auto b_dec = std::make_unique<std::array<TRLWE<P>, P::l̅>>();
        TRLWEBaseBbarDecompose<P>(*a_dec, a);
        TRLWEBaseBbarDecompose<P>(*b_dec, b);

        // 128-bit torus: 384-bit accumulation + Δ-division.  Digit
        // convolutions fit exactly in double for lvl3simd:
        // n * (Bbar/2)^2 with a small safety margin below the 53-bit mantissa.
        constexpr bool use_fft_digits =
            2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3 <
            std::numeric_limits<double>::digits;

        auto acc = std::make_unique<std::array<std::array<Wide384, P::n>, 3>>();

        if constexpr (use_fft_digits) {
            auto a_fd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
            auto b_fd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    TwistIFFT<P>((*a_fd)[i][c], (*a_dec)[i][c]);
                    TwistIFFT<P>((*b_fd)[i][c], (*b_dec)[i][c]);
                }
            }

            alignas(64) PolynomialInFD<P> prod_fd;
            Polynomial<P> prod;
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                    const int k = i + j;
                    const int shift =
                        2 * width - (k + 2) * static_cast<int>(P::B̅gbit);
                    if (shift <= -width) continue;

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][0], (*b_fd)[j][0]);
                    TwistFFT<P>(prod, prod_fd);
                    for (uint32_t n = 0; n < P::n; n++)
                        (*acc)[2][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][1], (*b_fd)[j][1]);
                    TwistFFT<P>(prod, prod_fd);
                    for (uint32_t n = 0; n < P::n; n++)
                        (*acc)[1][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][0], (*b_fd)[j][1]);
                    FMAInFD<P::n>(prod_fd, (*a_fd)[i][1], (*b_fd)[j][0]);
                    TwistFFT<P>(prod, prod_fd);
                    for (uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                }
            }
        }
        else {
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                    const int k = i + j;
                    const int shift =
                        2 * width - (k + 2) * static_cast<int>(P::B̅gbit);
                    if (shift <= -width) continue;

                    Polynomial<P> prod;
                    PolyMulNaive<P>(prod, (*a_dec)[i][0], (*b_dec)[j][0]);
                    for (uint32_t n = 0; n < P::n; n++)
                        (*acc)[2][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][1], (*b_dec)[j][1]);
                    for (uint32_t n = 0; n < P::n; n++)
                        (*acc)[1][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][0], (*b_dec)[j][1]);
                    for (uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][1], (*b_dec)[j][0]);
                    for (uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                }
            }
        }

        constexpr typename P::T delta = P::delta_int;
        for (int c = 0; c < 3; c++)
            for (uint32_t n = 0; n < P::n; n++)
                res[c][n] = (*acc)[c][n].div128(delta);
    }
    else if constexpr (is_multilimb_uint_v<typename P::T>) {
        using Torus = typename P::T;
        static_assert(P::l̅ * P::B̅gbit == width,
                      "FullDD multi-limb multiply requires Bbar digits to cover T");

        auto a_dec = std::make_unique<std::array<TRLWE<P>, P::l̅>>();
        auto b_dec = std::make_unique<std::array<TRLWE<P>, P::l̅>>();
        TRLWEBaseBbarDecompose<P>(*a_dec, a);
        TRLWEBaseBbarDecompose<P>(*b_dec, b);

        using Acc = WideSignedLimbAccumulator<2 * Torus::limbs + 2>;
        auto acc = std::make_unique<std::array<std::array<Acc, P::n>, 3>>();

        auto accumulate = [&](int comp, const Polynomial<P> &prod, int shift) {
            for (uint32_t n = 0; n < P::n; n++)
                (*acc)[comp][n].add_shifted_i64(
                    multilimb_to_signed_i64(prod[n]), shift);
        };

        constexpr bool use_fft_digits =
            use_multilimb_digit_fft_v<P> &&
            2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 4 <
                std::numeric_limits<double>::digits;

        if constexpr (use_fft_digits) {
            constexpr int digit_product_bits =
                2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3;
            constexpr int cross_term_product_bits = digit_product_bits + 1;
            constexpr int fd_batch_slack =
                std::numeric_limits<double>::digits - cross_term_product_bits;
            static_assert(fd_batch_slack > 0);
            constexpr int fd_batch_size =
                std::min(static_cast<int>(P::l̅), 1 << fd_batch_slack);

            auto a_fd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
            auto b_fd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    TwistIFFTDigit<P>((*a_fd)[i][c], (*a_dec)[i][c]);
                    TwistIFFTDigit<P>((*b_fd)[i][c], (*b_dec)[i][c]);
                }
            }

            auto sum_fd = std::make_unique<std::array<PolynomialInFD<P>, 3>>();
            auto prod = std::make_unique<Polynomial<P>>();
            for (int digit_sum = 0;
                 digit_sum <= 2 * static_cast<int>(P::l̅) - 2; digit_sum++) {
                const int i_begin =
                    std::max(0, digit_sum - static_cast<int>(P::l̅) + 1);
                const int i_end =
                    std::min(static_cast<int>(P::l̅) - 1, digit_sum);
                const int shift =
                    2 * width -
                    (digit_sum + 2) * static_cast<int>(P::B̅gbit);

                for (int batch_begin = i_begin; batch_begin <= i_end;
                     batch_begin += fd_batch_size) {
                    for (int c = 0; c < 3; c++) (*sum_fd)[c].fill(0.0);
                    const int batch_end =
                        std::min(i_end, batch_begin + fd_batch_size - 1);
                    for (int i = batch_begin; i <= batch_end; i++) {
                        const int j = digit_sum - i;
                        FMAInFD<P::n>((*sum_fd)[2], (*a_fd)[i][0],
                                      (*b_fd)[j][0]);
                        FMAInFD<P::n>((*sum_fd)[1], (*a_fd)[i][1],
                                      (*b_fd)[j][1]);
                        FMAInFD<P::n>((*sum_fd)[0], (*a_fd)[i][0],
                                      (*b_fd)[j][1]);
                        FMAInFD<P::n>((*sum_fd)[0], (*a_fd)[i][1],
                                      (*b_fd)[j][0]);
                    }
                    for (int c = 0; c < 3; c++) {
                        TwistFFTDigitProduct<P>(*prod, (*sum_fd)[c]);
                        accumulate(c, *prod, shift);
                    }
                }
            }
        }
        else {
            auto prod = std::make_unique<Polynomial<P>>();
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                    const int k = i + j;
                    const int shift =
                        2 * width - (k + 2) * static_cast<int>(P::B̅gbit);

                    PolyMulNaive<P>(*prod, (*a_dec)[i][0], (*b_dec)[j][0]);
                    accumulate(2, *prod, shift);
                    PolyMulNaive<P>(*prod, (*a_dec)[i][1], (*b_dec)[j][1]);
                    accumulate(1, *prod, shift);
                    PolyMulNaive<P>(*prod, (*a_dec)[i][0], (*b_dec)[j][1]);
                    accumulate(0, *prod, shift);
                    PolyMulNaive<P>(*prod, (*a_dec)[i][1], (*b_dec)[j][0]);
                    accumulate(0, *prod, shift);
                }
            }
        }

        constexpr Torus delta = P::delta_int;
        for (int c = 0; c < 3; c++)
            for (uint32_t n = 0; n < P::n; n++)
                res[c][n] = (*acc)[c][n].template div_to_torus<Torus::limbs>(
                    delta);
    }
    else {
        std::array<TRLWE<P>, P::l̅> a_dec, b_dec;
        TRLWEBaseBbarDecompose<P>(a_dec, a);
        TRLWEBaseBbarDecompose<P>(b_dec, b);

        // 32/64-bit torus: shift-based with power-of-2 Δ cancellation
        for (int n = 0; n < static_cast<int>(P::n); n++) {
            res[0][n] = 0; res[1][n] = 0; res[2][n] = 0;
        }
        for (int i = 0; i < static_cast<int>(P::l̅); i++) {
            for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                const int k = i + j;
                const int shift = width - (k + 2) * static_cast<int>(P::B̅gbit)
                                + static_cast<int>(P::plain_modulusbit);
                if (shift >= width || (shift < 0 && -shift >= width)) continue;

                Polynomial<P> prod;
                auto apply = [&](int comp) {
                    if (shift >= 0) {
                        for (int n = 0; n < static_cast<int>(P::n); n++)
                            res[comp][n] += prod[n] << shift;
                    } else {
                        for (int n = 0; n < static_cast<int>(P::n); n++)
                            res[comp][n] += static_cast<typename P::T>(
                                static_cast<std::make_signed_t<typename P::T>>(prod[n]) >> (-shift));
                    }
                };
                PolyMulNaive<P>(prod, a_dec[i][0], b_dec[j][0]); apply(2);
                PolyMulNaive<P>(prod, a_dec[i][1], b_dec[j][1]); apply(1);
                PolyMulNaive<P>(prod, a_dec[i][0], b_dec[j][1]); apply(0);
                PolyMulNaive<P>(prod, a_dec[i][1], b_dec[j][0]); apply(0);
            }
        }
    }
}

// Full DD TRLWE multiplication with relinearization
// Decomposes both input TRLWEs for the multiplication step (not just relinearization)
template <class P>
inline void TRLWEMultFullDD(TRLWE<P> &res, const TRLWE<P> &a, const TRLWE<P> &b,
                             const relinKeyFFT<P> &relinkeyfft)
{
    TRLWE3<P> resmult;
    TRLWEMultWithoutRelinearizationFullDD<P>(resmult, a, b);
    Relinearization<P>(res, resmult, relinkeyfft);
}

// TLWE multiplication - automatically handles DD when l̅ > 1
template <class P>
inline void TLWEMult(TLWE<typename P::targetP> &res,
                     const TLWE<typename P::domainP> &a,
                     const TLWE<typename P::domainP> &b,
                     const relinKeyFFT<typename P::targetP> &relinkeyfft,
                     const PrivateKeySwitchingKey<P> &privksk)
{
    TRLWE<typename P::targetP> trlweres, trlwea, trlweb;
    PrivKeySwitch<P>(trlwea, a, privksk);
    PrivKeySwitch<P>(trlweb, b, privksk);
    TRLWEMult<typename P::targetP>(trlweres, trlwea, trlweb, relinkeyfft);
    SampleExtractIndex<typename P::targetP>(res, trlweres, 0);
}

}  // namespace TFHEpp
