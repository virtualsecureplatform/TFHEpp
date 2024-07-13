#pragma once

// Poted from here 
// https://github.com/zama-ai/tfhe-rs/blob/main/tfhe/src/core_crypto/fft_impl/fft64/math/fft/x86.rs

//! For documentation on the various intrinsics used here, refer to Intel's intrinsics guide.
//! <https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html>
//!
//! currently we dispatch based on the availability of
//!  - avx+avx2(advanced vector extensions for 256 intrinsics)+fma(fused multiply add for complex
//!    multiplication, usually comes with avx+avx2),
//!  - or the availability of avx512f[+avx512dq(doubleword/quadword intrinsics for conversion of f64
//!    to/from i64. usually comes with avx512f on modern cpus)]
//!
//! more dispatch options may be added in the future


#include <immintrin.h>

namespace SPQLIOS{

// Convert a vector of f64 values to a vector of i64 values.
// See `f64_to_i64_bit_twiddles` in `fft/tests.rs` for the scalar version.
inline __m256i mm256_cvtpd_epi64(const __m256d x) {
    // reinterpret the bits as u64 values
    const __m256i bits = _mm256_castpd_si256(x);

    // mask that covers the first 52 bits
    const __m256i mantissa_mask = _mm256_set1_epi64x(0xFFFFFFFFFFFFF);

    // mask that covers the 52nd bit
    const __m256i explicit_mantissa_bit = _mm256_set1_epi64x(0x10000000000000);

    // mask that covers the first 11 bits
    const __m256i exp_mask = _mm256_set1_epi64x(0x7FF);

    // extract the first 52 bits and add the implicit bit
    const __m256i mantissa = _mm256_or_si256(
        _mm256_and_si256(bits, mantissa_mask),
        explicit_mantissa_bit
    );

    // extract the 52nd to 63rd (excluded) bits for the biased exponent
    const __m256i biased_exp = _mm256_and_si256(_mm256_srli_epi64(bits, 52), exp_mask);

    // extract the 63rd sign bit
    const __m256i sign_is_negative_mask = _mm256_sub_epi64(
        _mm256_setzero_si256(),
        _mm256_srli_epi64(bits, 63)
    );

    // we need to shift the mantissa by some value that may be negative, so we first shift it to
    // the left by the maximum amount, then shift it to the right by our value plus the offset we
    // just shifted by
    //
    // the 52nd bit is set to 1, so we shift to the left by 11 so the 63rd (last) bit is set.
    const __m256i mantissa_lshift = _mm256_slli_epi64(mantissa, 11);

    // shift to the right and apply the exponent bias
    // If biased_exp == 0 then we have 0 or a subnormal value which should return 0, here we will
    // shift to the right by 1086 which will return 0 as we are shifting in 0s from the left, so
    // subnormals are already covered
    const __m256i mantissa_shift = _mm256_srlv_epi64(
        mantissa_lshift,
        _mm256_sub_epi64(_mm256_set1_epi64x(1086), biased_exp)
    );

    // if the sign bit is unset, we keep our result
    const __m256i value_if_positive = mantissa_shift;
    // otherwise, we negate it
    const __m256i value_if_negative = _mm256_sub_epi64(_mm256_setzero_si256(), value_if_positive);

    // if the biased exponent is all zeros, we have a subnormal value (or zero)

    // Select the value based on the sign mask
    return _mm256_blendv_epi8(value_if_positive, value_if_negative, sign_is_negative_mask);
}

// Convert a vector of i64 values to a vector of f64 values. Not sure how it works.
// Ported from <https://stackoverflow.com/a/41148578>.

inline __m256d mm256_cvtepi64_pd(const __m256i x) {
    // Shift right arithmetic by 16 bits
    __m256i x_hi = _mm256_srai_epi32(x, 16);

    // Blend lower 16 bits with zero
    x_hi = _mm256_blend_epi16(x_hi, _mm256_setzero_si256(), 0x33);

    // Add a large constant to x_hi
    x_hi = _mm256_add_epi64(
        x_hi,
        _mm256_castpd_si256(_mm256_set1_pd(442721857769029238784.0)) // 3*2^67
    );

    // Blend specific bits of x with a constant
    __m256i x_lo = _mm256_blend_epi16(
        x,
        _mm256_castpd_si256(_mm256_set1_pd(4503599627370496.0)), // 2^52
        0x88
    );

    // Subtract a large constant from the floating-point representation of x_hi
    __m256d f = _mm256_sub_pd(
        _mm256_castsi256_pd(x_hi),
        _mm256_set1_pd(442726361368656609280.0) // 3*2^67 + 2^52
    );

    // Add the floating-point representation of x_lo to the result
    return _mm256_add_pd(f, _mm256_castsi256_pd(x_lo));
}

// Convert a vector of i64 values to a vector of i32 values without AVX512.
// https://stackoverflow.com/questions/69408063/how-to-convert-int-64-to-int-32-with-avx-but-without-avx-512
// Slower
// __m128 mm256_cvtepi64_epi32_avx(const __m256i v)
// {
//    const __m256 vf = _mm256_castsi256_ps( v );      // free
//    const __m128 hi = _mm256_extractf128_ps(vf, 1);  // vextractf128
//    const __m128 lo = _mm256_castps256_ps128( vf );  // also free
//    // take the bottom 32 bits of each 64-bit chunk in lo and hi
//    const __m128 packed = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(2, 0, 2, 0));  // shufps
//    //return _mm_castps_si128(packed);  // if you want
//    return packed;
// }

// 2x 256 -> 1x 256-bit result
__m256i pack64to32(__m256i a, __m256i b)
{
    // grab the 32-bit low halves of 64-bit elements into one vector
   __m256 combined = _mm256_shuffle_ps(_mm256_castsi256_ps(a),
                                       _mm256_castsi256_ps(b), _MM_SHUFFLE(2,0,2,0));
    // {b3,b2, a3,a2 | b1,b0, a1,a0}  from high to low

    // re-arrange pairs of 32-bit elements with vpermpd (or vpermq if you want)
    __m256d ordered = _mm256_permute4x64_pd(_mm256_castps_pd(combined), _MM_SHUFFLE(3,1,2,0));
    return _mm256_castpd_si256(ordered);
}


inline void convert_f64_to_u32(uint32_t* const res, const double* const real_inout_direct, const int32_t N) {
    for (int32_t i = 0; i < N; i += 8) {
        // Load 4 double values
        const __m256d real_vals1 = _mm256_loadu_pd(&real_inout_direct[i]);
        const __m256d real_vals2 = _mm256_loadu_pd(&real_inout_direct[i + 4]);

        // Convert double to int64
        const __m256i int64_vals1 = mm256_cvtpd_epi64(real_vals1);
        const __m256i int64_vals2 = mm256_cvtpd_epi64(real_vals2);

        const __m256i packed32 = pack64to32(int64_vals1,int64_vals2);

        // Store the result
        _mm256_storeu_si256((__m256i*)&res[i], packed32);
    }
}
}