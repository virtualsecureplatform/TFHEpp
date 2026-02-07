#pragma once
#include <cstdint>

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
    alignas(64) DecomposedPolynomial<P> decvec;
    Decomposition<P>(decvec, poly);
    alignas(64) PolynomialInFD<P> decvecfft;

    if constexpr (P::l̅ > 1) {
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

// Full Double Decomposition TRLWE Multiplication (without relinearization)
// Both TRLWEs are decomposed by l̅, multiplication is polynomial-like in decomposition indices
// Algorithm:
//   1. Decompose a[k] and b[k] into l̅ components each using base B̅g
//   2. For polynomial product, compute convolution in decomposition index space:
//      (Σᵢ aᵢ·B̅g^i) × (Σⱼ bⱼ·B̅g^j) = Σₖ (Σᵢ₊ⱼ₌ₖ aᵢ·bⱼ)·B̅g^k
//   3. Each aᵢ·bⱼ is computed via FFT polynomial multiplication
//   4. IFFT to recover coefficients, rescale by Δ
//   5. Recombine the 2l̅-1 terms back to proper scaling
template <class P>
void TRLWEMultWithoutRelinearizationFullDD(TRLWE3<P> &res, const TRLWE<P> &a,
                                            const TRLWE<P> &b)
{
    constexpr int width = std::numeric_limits<typename P::T>::digits;

    // Decompose all input polynomials into l̅ components
    // TRLWEBaseBbarDecompose gives: a = Σⱼ a_dec[j] * 2^(width - (j+1)*B̅gbit)
    // where a_dec[j] has coefficients in [-B̅g/2, B̅g/2)
    std::array<TRLWE<P>, P::l̅> a_dec, b_dec;
    TRLWEBaseBbarDecompose<P>(a_dec, a);
    TRLWEBaseBbarDecompose<P>(b_dec, b);

    // FFT all decomposed components
    std::array<std::array<PolynomialInFD<P>, P::l̅>, P::k + 1> a_fft, b_fft;
    for (int poly_idx = 0; poly_idx <= P::k; poly_idx++) {
        for (int j = 0; j < P::l̅; j++) {
            TwistIFFT<P>(a_fft[poly_idx][j], a_dec[j][poly_idx]);
            TwistIFFT<P>(b_fft[poly_idx][j], b_dec[j][poly_idx]);
        }
    }

    // Compute c[0] = a[0]*b[1] + a[1]*b[0] using DD polynomial multiplication
    std::array<PolynomialInFD<P>, 2 * P::l̅ - 1> c0_fft;
    for (int k = 0; k < 2 * P::l̅ - 1; k++) {
        for (int n = 0; n < P::n; n++) c0_fft[k][n] = 0.0;
    }
    for (int i = 0; i < P::l̅; i++) {
        for (int j = 0; j < P::l̅; j++) {
            FMAInFD<P::n>(c0_fft[i + j], a_fft[0][i], b_fft[1][j]);
            FMAInFD<P::n>(c0_fft[i + j], a_fft[1][i], b_fft[0][j]);
        }
    }

    // Compute c[1] = a[1]*b[1]
    std::array<PolynomialInFD<P>, 2 * P::l̅ - 1> c1_fft;
    for (int k = 0; k < 2 * P::l̅ - 1; k++) {
        for (int n = 0; n < P::n; n++) c1_fft[k][n] = 0.0;
    }
    for (int i = 0; i < P::l̅; i++) {
        for (int j = 0; j < P::l̅; j++) {
            FMAInFD<P::n>(c1_fft[i + j], a_fft[1][i], b_fft[1][j]);
        }
    }

    // Compute c[2] = a[0]*b[0]
    std::array<PolynomialInFD<P>, 2 * P::l̅ - 1> c2_fft;
    for (int k = 0; k < 2 * P::l̅ - 1; k++) {
        for (int n = 0; n < P::n; n++) c2_fft[k][n] = 0.0;
    }
    for (int i = 0; i < P::l̅; i++) {
        for (int j = 0; j < P::l̅; j++) {
            FMAInFD<P::n>(c2_fft[i + j], a_fft[0][i], b_fft[0][j]);
        }
    }

    // Initialize results to zero
    for (int n = 0; n < P::n; n++) {
        res[0][n] = 0;
        res[1][n] = 0;
        res[2][n] = 0;
    }

    // IFFT and recombine with Δ rescaling
    // The decomposition was: x = Σⱼ x_j * h̅[j] where h̅[j] = 2^(width - (j+1)*B̅gbit)
    // Product of positions i and j: scale = h̅[i] * h̅[j] = 2^(2*width - (i+j+2)*B̅gbit)
    // At convolution position k = i+j: scale = 2^(2*width - (k+2)*B̅gbit)
    //
    // The original ciphertext coefficients are scaled by Δ = 2^(width - plain_modulusbit).
    // After multiplying a × b = Δ² × plaintext_product.
    // We need to divide by Δ to get back to Δ scaling.
    //
    // For convolution position k:
    //   raw_scale = 2^(2*width - (k+2)*B̅gbit)
    //   after /Δ: scale = raw_scale / Δ = 2^(2*width - (k+2)*B̅gbit) / 2^(width - plain_modulusbit)
    //                   = 2^(width - (k+2)*B̅gbit + plain_modulusbit)
    //
    // So the shift for position k is: width - (k+2)*B̅gbit + plain_modulusbit

    for (int k = 0; k < 2 * P::l̅ - 1; k++) {
        Polynomial<P> temp0, temp1, temp2;
        TwistFFT<P>(temp0, c0_fft[k]);
        TwistFFT<P>(temp1, c1_fft[k]);
        TwistFFT<P>(temp2, c2_fft[k]);

        // Shift includes the Δ division: width - (k+2)*B̅gbit + plain_modulusbit
        // Positive shift means left shift, negative means right shift
        const int shift = width - (k + 2) * P::B̅gbit + P::plain_modulusbit;

        if (shift >= 0 && shift < width) {
            // Left shift
            for (int n = 0; n < P::n; n++) {
                res[0][n] += temp0[n] << shift;
                res[1][n] += temp1[n] << shift;
                res[2][n] += temp2[n] << shift;
            }
        } else if (shift < 0 && -shift < width) {
            // Right shift
            const int right_shift = -shift;
            for (int n = 0; n < P::n; n++) {
                res[0][n] += static_cast<typename P::T>(
                    static_cast<std::make_signed_t<typename P::T>>(temp0[n]) >> right_shift);
                res[1][n] += static_cast<typename P::T>(
                    static_cast<std::make_signed_t<typename P::T>>(temp1[n]) >> right_shift);
                res[2][n] += static_cast<typename P::T>(
                    static_cast<std::make_signed_t<typename P::T>>(temp2[n]) >> right_shift);
            }
        }
        // If |shift| >= width, contribution is zero
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
