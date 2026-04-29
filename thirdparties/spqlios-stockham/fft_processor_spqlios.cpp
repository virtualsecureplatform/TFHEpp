// Stockham radix-4 FFT with split real/imaginary layout.
//
// Data layout: re[0..ns2-1], im[0..ns2-1] (same as original SPQLIOS)
// Algorithm: Stockham auto-sort radix-4 DIF (out-of-place)
//
// The split-format radix-4 butterfly:
//   apc_re = a_re+c_re,  apc_im = a_im+c_im
//   amc_re = a_re-c_re,  amc_im = a_im-c_im
//   bpd_re = b_re+d_re,  bpd_im = b_im+d_im
//   bmd_re = b_re-d_re,  bmd_im = b_im-d_im
//   jbmd_re = -bmd_im,   jbmd_im = bmd_re   (mul by j for fwd)
//
//   out0 = apc + bpd
//   out1 = (amc - jbmd) * w1
//   out2 = (apc - bpd) * w2
//   out3 = (amc + jbmd) * w3
//
// The j-multiply is free: just swap re/im and negate the real part.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <immintrin.h>

#include <params.hpp>
#include "fft_processor_spqlios.h"

// ── Trig helpers ─────────────────────────────────────────────────────────────
static double accurate_cos(int32_t i, int32_t n) {
    i = ((i % n) + n) % n;
    if (i >= 3*n/4) return cos(2.*M_PI*(n-i)/double(n));
    if (i >= 2*n/4) return -cos(2.*M_PI*(i-n/2)/double(n));
    if (i >= 1*n/4) return -cos(2.*M_PI*(n/2-i)/double(n));
    return cos(2.*M_PI*i/double(n));
}
static double accurate_sin(int32_t i, int32_t n) {
    i = ((i % n) + n) % n;
    if (i >= 3*n/4) return -sin(2.*M_PI*(n-i)/double(n));
    if (i >= 2*n/4) return -sin(2.*M_PI*(i-n/2)/double(n));
    if (i >= 1*n/4) return sin(2.*M_PI*(n/2-i)/double(n));
    return sin(2.*M_PI*i/double(n));
}

// ── Split-format complex multiply with AVX2 ──────────────────────────────────
// res_re = a_re*w_re - a_im*w_sin
// res_im = a_im*w_re + a_re*w_sin (actually a_re*w_sin + a_im*w_cos = a_im*w_re + a_re*w_sin)
// Wait, standard complex multiply:
//   (a_re + i*a_im) * (w_re + i*w_im) = (a_re*w_re - a_im*w_im) + i*(a_re*w_im + a_im*w_re)

static inline void cmul_split(__m256d a_re, __m256d a_im,
                               __m256d w_re, __m256d w_im,
                               __m256d &r_re, __m256d &r_im) {
    r_re = _mm256_mul_pd(a_re, w_re);
    r_re = _mm256_fnmadd_pd(a_im, w_im, r_re);  // a_re*w_re - a_im*w_im
    r_im = _mm256_mul_pd(a_im, w_re);
    r_im = _mm256_fmadd_pd(a_re, w_im, r_im);   // a_im*w_re + a_re*w_im
}

// ── Table structure ──────────────────────────────────────────────────────────
struct STOCKHAM_PRECOMP {
    int32_t n;
    int32_t ns2;
    double *trig_fwd;    // butterfly twiddles for forward FFT
    double *trig_inv;    // butterfly twiddles for inverse FFT
    double *twist_fwd;   // forward twist twiddles
    double *twist_inv;   // inverse twist twiddles
    double *scratch_re;  // scratch buffer (re part)
    double *scratch_im;  // scratch buffer (im part)
    double *data;        // data buffer for execute_direct
    void *buf;
};

// ── Stockham radix-4 with split layout ──────────────────────────────────────

static void stockham_r4_split(
    int32_t ns2, bool fwd,
    const double *trig,
    double *re, double *im,
    double *sre, double *sim)
{
    double *src_re = re, *src_im = im;
    double *dst_re = sre, *dst_im = sim;
    const double *tw = trig;
    int32_t q = ns2 / 4;
    int32_t s = 1;
    const __m256d neg = _mm256_set1_pd(-0.0);

    while (q >= 4) {
        const int32_t stride = q * s;  // = ns2/4

        // j=0 pass: no twiddle (w1=w2=w3=1)
        for (int32_t p = 0; p < q; p += 4) {
            const int32_t idx = s * p;  // j=0, so idx = s*p
            const double *sr = src_re + idx, *si = src_im + idx;
            __m256d a_re = _mm256_loadu_pd(sr);
            __m256d a_im = _mm256_loadu_pd(si);
            __m256d b_re = _mm256_loadu_pd(sr + stride);
            __m256d b_im = _mm256_loadu_pd(si + stride);
            __m256d c_re = _mm256_loadu_pd(sr + 2*stride);
            __m256d c_im = _mm256_loadu_pd(si + 2*stride);
            __m256d d_re = _mm256_loadu_pd(sr + 3*stride);
            __m256d d_im = _mm256_loadu_pd(si + 3*stride);

            __m256d apc_re = _mm256_add_pd(a_re, c_re);
            __m256d apc_im = _mm256_add_pd(a_im, c_im);
            __m256d amc_re = _mm256_sub_pd(a_re, c_re);
            __m256d amc_im = _mm256_sub_pd(a_im, c_im);
            __m256d bpd_re = _mm256_add_pd(b_re, d_re);
            __m256d bpd_im = _mm256_add_pd(b_im, d_im);
            __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
            __m256d bmd_im = _mm256_sub_pd(b_im, d_im);

            __m256d jbmd_re, jbmd_im;
            if (fwd) { jbmd_re = _mm256_xor_pd(bmd_im, neg); jbmd_im = bmd_re; }
            else     { jbmd_re = bmd_im; jbmd_im = _mm256_xor_pd(bmd_re, neg); }

            _mm256_storeu_pd(dst_re + p, _mm256_add_pd(apc_re, bpd_re));
            _mm256_storeu_pd(dst_im + p, _mm256_add_pd(apc_im, bpd_im));
            _mm256_storeu_pd(dst_re + q + p, _mm256_sub_pd(amc_re, jbmd_re));
            _mm256_storeu_pd(dst_im + q + p, _mm256_sub_pd(amc_im, jbmd_im));
            _mm256_storeu_pd(dst_re + 2*q + p, _mm256_sub_pd(apc_re, bpd_re));
            _mm256_storeu_pd(dst_im + 2*q + p, _mm256_sub_pd(apc_im, bpd_im));
            _mm256_storeu_pd(dst_re + 3*q + p, _mm256_add_pd(amc_re, jbmd_re));
            _mm256_storeu_pd(dst_im + 3*q + p, _mm256_add_pd(amc_im, jbmd_im));
        }

        // j=1..s-1 passes: with twiddle multiply
        for (int32_t j = 1; j < s; j++) {
            const double *twj = tw + j * 6;
            __m256d w1_re = _mm256_set1_pd(twj[0]), w1_im = _mm256_set1_pd(twj[1]);
            __m256d w2_re = _mm256_set1_pd(twj[2]), w2_im = _mm256_set1_pd(twj[3]);
            __m256d w3_re = _mm256_set1_pd(twj[4]), w3_im = _mm256_set1_pd(twj[5]);
            const int32_t base = q * 4 * j;

            for (int32_t p = 0; p < q; p += 4) {
                const int32_t idx = j + s * p;
                const double *sr = src_re + idx, *si = src_im + idx;
                __m256d a_re = _mm256_loadu_pd(sr);
                __m256d a_im = _mm256_loadu_pd(si);
                __m256d b_re = _mm256_loadu_pd(sr + stride);
                __m256d b_im = _mm256_loadu_pd(si + stride);
                __m256d c_re = _mm256_loadu_pd(sr + 2*stride);
                __m256d c_im = _mm256_loadu_pd(si + 2*stride);
                __m256d d_re = _mm256_loadu_pd(sr + 3*stride);
                __m256d d_im = _mm256_loadu_pd(si + 3*stride);

                __m256d apc_re = _mm256_add_pd(a_re, c_re);
                __m256d apc_im = _mm256_add_pd(a_im, c_im);
                __m256d amc_re = _mm256_sub_pd(a_re, c_re);
                __m256d amc_im = _mm256_sub_pd(a_im, c_im);
                __m256d bpd_re = _mm256_add_pd(b_re, d_re);
                __m256d bpd_im = _mm256_add_pd(b_im, d_im);
                __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
                __m256d bmd_im = _mm256_sub_pd(b_im, d_im);

                __m256d jbmd_re, jbmd_im;
                if (fwd) { jbmd_re = _mm256_xor_pd(bmd_im, neg); jbmd_im = bmd_re; }
                else     { jbmd_re = bmd_im; jbmd_im = _mm256_xor_pd(bmd_re, neg); }

                // out0 = apc + bpd (no twiddle)
                _mm256_storeu_pd(dst_re + base + p, _mm256_add_pd(apc_re, bpd_re));
                _mm256_storeu_pd(dst_im + base + p, _mm256_add_pd(apc_im, bpd_im));

                // out1 = (amc - jbmd) * w1
                __m256d t_re = _mm256_sub_pd(amc_re, jbmd_re);
                __m256d t_im = _mm256_sub_pd(amc_im, jbmd_im);
                __m256d r_re = _mm256_mul_pd(t_re, w1_re);
                r_re = _mm256_fnmadd_pd(t_im, w1_im, r_re);
                __m256d r_im = _mm256_mul_pd(t_im, w1_re);
                r_im = _mm256_fmadd_pd(t_re, w1_im, r_im);
                _mm256_storeu_pd(dst_re + base + q + p, r_re);
                _mm256_storeu_pd(dst_im + base + q + p, r_im);

                // out2 = (apc - bpd) * w2
                t_re = _mm256_sub_pd(apc_re, bpd_re);
                t_im = _mm256_sub_pd(apc_im, bpd_im);
                r_re = _mm256_mul_pd(t_re, w2_re);
                r_re = _mm256_fnmadd_pd(t_im, w2_im, r_re);
                r_im = _mm256_mul_pd(t_im, w2_re);
                r_im = _mm256_fmadd_pd(t_re, w2_im, r_im);
                _mm256_storeu_pd(dst_re + base + 2*q + p, r_re);
                _mm256_storeu_pd(dst_im + base + 2*q + p, r_im);

                // out3 = (amc + jbmd) * w3
                t_re = _mm256_add_pd(amc_re, jbmd_re);
                t_im = _mm256_add_pd(amc_im, jbmd_im);
                r_re = _mm256_mul_pd(t_re, w3_re);
                r_re = _mm256_fnmadd_pd(t_im, w3_im, r_re);
                r_im = _mm256_mul_pd(t_im, w3_re);
                r_im = _mm256_fmadd_pd(t_re, w3_im, r_im);
                _mm256_storeu_pd(dst_re + base + 3*q + p, r_re);
                _mm256_storeu_pd(dst_im + base + 3*q + p, r_im);
            }
        }

        tw += s * 6;
        double *t;
        t = src_re; src_re = dst_re; dst_re = t;
        t = src_im; src_im = dst_im; dst_im = t;
        s *= 4; q /= 4;
    }

    // Final pass
    if (q == 1) {
        const int32_t stride = s;
        for (int32_t j = 0; j < s; j += 4) {
            __m256d a_re = _mm256_loadu_pd(src_re + j);
            __m256d a_im = _mm256_loadu_pd(src_im + j);
            __m256d c_re = _mm256_loadu_pd(src_re + j + 2*stride);
            __m256d c_im = _mm256_loadu_pd(src_im + j + 2*stride);
            __m256d apc_re = _mm256_add_pd(a_re, c_re);
            __m256d apc_im = _mm256_add_pd(a_im, c_im);
            __m256d amc_re = _mm256_sub_pd(a_re, c_re);
            __m256d amc_im = _mm256_sub_pd(a_im, c_im);
            __m256d b_re = _mm256_loadu_pd(src_re + j + stride);
            __m256d b_im = _mm256_loadu_pd(src_im + j + stride);
            __m256d d_re = _mm256_loadu_pd(src_re + j + 3*stride);
            __m256d d_im = _mm256_loadu_pd(src_im + j + 3*stride);
            __m256d bpd_re = _mm256_add_pd(b_re, d_re);
            __m256d bpd_im = _mm256_add_pd(b_im, d_im);
            __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
            __m256d bmd_im = _mm256_sub_pd(b_im, d_im);
            __m256d jbmd_re, jbmd_im;
            if (fwd) { jbmd_re = _mm256_xor_pd(bmd_im, neg); jbmd_im = bmd_re; }
            else     { jbmd_re = bmd_im; jbmd_im = _mm256_xor_pd(bmd_re, neg); }
            // Scatter to Stockham positions
            double t_re[16], t_im[16];
            _mm256_storeu_pd(t_re+0, _mm256_add_pd(apc_re, bpd_re));
            _mm256_storeu_pd(t_im+0, _mm256_add_pd(apc_im, bpd_im));
            _mm256_storeu_pd(t_re+4, _mm256_sub_pd(amc_re, jbmd_re));
            _mm256_storeu_pd(t_im+4, _mm256_sub_pd(amc_im, jbmd_im));
            _mm256_storeu_pd(t_re+8, _mm256_sub_pd(apc_re, bpd_re));
            _mm256_storeu_pd(t_im+8, _mm256_sub_pd(apc_im, bpd_im));
            _mm256_storeu_pd(t_re+12, _mm256_add_pd(amc_re, jbmd_re));
            _mm256_storeu_pd(t_im+12, _mm256_add_pd(amc_im, jbmd_im));
            for (int k = 0; k < 4; k++)
                for (int r = 0; r < 4; r++) {
                    dst_re[4*(j+k)+r] = t_re[r*4+k];
                    dst_im[4*(j+k)+r] = t_im[r*4+k];
                }
        }
    } else if (q == 2) {
        const int32_t half = ns2 / 2;
        for (int32_t j = 0; j < half; j += 4) {
            __m256d a_re = _mm256_loadu_pd(src_re + j);
            __m256d a_im = _mm256_loadu_pd(src_im + j);
            __m256d b_re = _mm256_loadu_pd(src_re + j + half);
            __m256d b_im = _mm256_loadu_pd(src_im + j + half);
            double ar[4], ai[4], br[4], bi[4];
            _mm256_storeu_pd(ar, _mm256_add_pd(a_re, b_re));
            _mm256_storeu_pd(ai, _mm256_add_pd(a_im, b_im));
            _mm256_storeu_pd(br, _mm256_sub_pd(a_re, b_re));
            _mm256_storeu_pd(bi, _mm256_sub_pd(a_im, b_im));
            for (int k = 0; k < 4; k++) {
                dst_re[2*(j+k)]   = ar[k]; dst_im[2*(j+k)]   = ai[k];
                dst_re[2*(j+k)+1] = br[k]; dst_im[2*(j+k)+1] = bi[k];
            }
        }
    }

    if (dst_re != re) {
        memcpy(re, dst_re, ns2 * sizeof(double));
        memcpy(im, dst_im, ns2 * sizeof(double));
    }
}


// ── Table construction ──────────────────────────────────────────────────────

static STOCKHAM_PRECOMP *build_stockham_tables(int32_t nn) {
    int32_t n = 2 * nn;
    int32_t ns2 = nn / 2;
    auto *reps = new STOCKHAM_PRECOMP;
    reps->n = n;
    reps->ns2 = ns2;

    // Count twiddle doubles for butterfly passes
    int32_t bf_doubles = 0;
    for (int32_t s = 1, q = ns2/4; q >= 4; s *= 4, q /= 4)
        bf_doubles += s * 6;  // s twiddle entries × 6 doubles (3 complex)

    int32_t twist_doubles = ns2 * 2;  // cos+sin for ns2/4 values... actually ns2 cos + ns2 sin
    // Twist: for each j=0..ns2-1: cos(-j/n) and sin(-j/n)
    // Split: cos[0..ns2-1], sin[0..ns2-1]
    twist_doubles = ns2 * 2;

    int32_t total = 2 * (bf_doubles + twist_doubles) + 2 * ns2 + nn;
    reps->buf = aligned_alloc(64, total * sizeof(double));
    double *ptr = (double *)reps->buf;

    reps->trig_fwd = ptr; ptr += bf_doubles + twist_doubles;
    reps->trig_inv = ptr; ptr += bf_doubles + twist_doubles;
    reps->scratch_re = ptr; ptr += ns2;
    reps->scratch_im = ptr; ptr += ns2;
    reps->data = ptr;  // nn doubles for execute_direct scratch

    // Build forward twiddles
    double *fwd = reps->trig_fwd;
    for (int32_t s = 1, q = ns2/4; q >= 4; s *= 4, q /= 4) {
        for (int32_t j = 0; j < s; j++) {
            int32_t denom = 4 * s;
            // w1 = exp(-2πi*j/(4s))
            *fwd++ = accurate_cos(-j, denom);
            *fwd++ = accurate_sin(-j, denom);
            // w2 = exp(-2πi*2j/(4s))
            *fwd++ = accurate_cos(-2*j, denom);
            *fwd++ = accurate_sin(-2*j, denom);
            // w3 = exp(-2πi*3j/(4s))
            *fwd++ = accurate_cos(-3*j, denom);
            *fwd++ = accurate_sin(-3*j, denom);
        }
    }
    reps->twist_fwd = fwd;
    for (int32_t j = 0; j < ns2; j++) {
        *fwd++ = accurate_cos(-j, n);
        *fwd++ = accurate_sin(-j, n);
    }

    // Build inverse twiddles
    double *inv = reps->trig_inv;
    for (int32_t s = 1, q = ns2/4; q >= 4; s *= 4, q /= 4) {
        for (int32_t j = 0; j < s; j++) {
            int32_t denom = 4 * s;
            *inv++ = accurate_cos(j, denom);
            *inv++ = accurate_sin(j, denom);
            *inv++ = accurate_cos(2*j, denom);
            *inv++ = accurate_sin(2*j, denom);
            *inv++ = accurate_cos(3*j, denom);
            *inv++ = accurate_sin(3*j, denom);
        }
    }
    reps->twist_inv = inv;
    for (int32_t j = 0; j < ns2; j++) {
        *inv++ = accurate_cos(j, n);
        *inv++ = accurate_sin(j, n);
    }

    return reps;
}

// ── FFT/IFFT wrappers ──────────────────────────────────────────────────────

static void split_fft(const STOCKHAM_PRECOMP *tables, double *data) {
    int32_t ns2 = tables->ns2;
    double *re = data;
    double *im = data + ns2;

    stockham_r4_split(ns2, true, tables->trig_fwd,
                      re, im, tables->scratch_re, tables->scratch_im);

    // Apply twist
    const double *tw = tables->twist_fwd;
    for (int32_t j = 0; j < ns2; j += 4) {
        __m256d r = _mm256_loadu_pd(re + j);
        __m256d i = _mm256_loadu_pd(im + j);
        // Load cos and sin (stored as [cos0,sin0,cos1,sin1,...])
        // Need to deinterleave
        __m256d t0 = _mm256_loadu_pd(tw + j*2);      // [c0,s0,c1,s1]
        __m256d t1 = _mm256_loadu_pd(tw + j*2 + 4);  // [c2,s2,c3,s3]
        __m256d tc = _mm256_shuffle_pd(t0, t1, 0b0000); // [c0,c1,c2,c3]
        __m256d ts = _mm256_shuffle_pd(t0, t1, 0b1111); // [s0,s1,s2,s3]
        // Hmm, shuffle_pd with 0b0000 gives [t0[0],t1[0],t0[2],t1[2]] = [c0,c2,c1,c3]
        // That's wrong. Let me use unpacklo/unpackhi:
        tc = _mm256_unpacklo_pd(t0, t1);  // [c0,c2,s0,s2]... also wrong for cross-lane
        // Actually for 256-bit, unpacklo operates per 128-bit lane:
        // lane0: [t0[0],t1[0]] = [c0,c2], lane1: [t0[2],t1[2]] = [c1,c3]
        // = [c0,c2,c1,c3]... still not what we want.
        // Just do scalar for now:
        double cos_vals[4], sin_vals[4];
        for (int k = 0; k < 4; k++) {
            cos_vals[k] = tw[(j+k)*2];
            sin_vals[k] = tw[(j+k)*2 + 1];
        }
        tc = _mm256_loadu_pd(cos_vals);
        ts = _mm256_loadu_pd(sin_vals);

        __m256d nr = _mm256_mul_pd(r, tc);
        nr = _mm256_fnmadd_pd(i, ts, nr);
        __m256d ni = _mm256_mul_pd(i, tc);
        ni = _mm256_fmadd_pd(r, ts, ni);
        _mm256_storeu_pd(re + j, nr);
        _mm256_storeu_pd(im + j, ni);
    }
}

static void split_ifft(const STOCKHAM_PRECOMP *tables, double *data) {
    int32_t ns2 = tables->ns2;
    double *re = data;
    double *im = data + ns2;

    // Apply twist first
    const double *tw = tables->twist_inv;
    for (int32_t j = 0; j < ns2; j += 4) {
        __m256d r = _mm256_loadu_pd(re + j);
        __m256d i = _mm256_loadu_pd(im + j);
        double cos_vals[4], sin_vals[4];
        for (int k = 0; k < 4; k++) {
            cos_vals[k] = tw[(j+k)*2];
            sin_vals[k] = tw[(j+k)*2 + 1];
        }
        __m256d tc = _mm256_loadu_pd(cos_vals);
        __m256d ts = _mm256_loadu_pd(sin_vals);
        __m256d nr = _mm256_mul_pd(r, tc);
        nr = _mm256_fnmadd_pd(i, ts, nr);
        __m256d ni = _mm256_mul_pd(i, tc);
        ni = _mm256_fmadd_pd(r, ts, ni);
        _mm256_storeu_pd(re + j, nr);
        _mm256_storeu_pd(im + j, ni);
    }

    stockham_r4_split(ns2, false, tables->trig_inv,
                      re, im, tables->scratch_re, tables->scratch_im);
}

// ── Conversion helpers ──────────────────────────────────────────────────────

static inline __m256i magic_cvtpd_epi64(__m256d x) {
    const __m256d m = _mm256_set1_pd(6755399441055744.0);
    const __m256i mi = _mm256_set1_epi64x(0x4338000000000000LL);
    return _mm256_sub_epi64(_mm256_castpd_si256(_mm256_add_pd(x, m)), mi);
}
static inline __m256d magic_cvtepi64_pd(__m256i x) {
    const __m256i mi = _mm256_set1_epi64x(0x4338000000000000LL);
    const __m256d m = _mm256_set1_pd(6755399441055744.0);
    return _mm256_sub_pd(_mm256_castsi256_pd(_mm256_add_epi64(x, mi)), m);
}

// ── FFT_Processor implementation ────────────────────────────────────────────

static int32_t rev(int32_t x, int32_t M) {
    int32_t r = 0;
    for (int32_t j = M; j > 1; j /= 2) { r = 2*r+(x%2); x /= 2; }
    return r;
}

FFT_Processor_Spqlios::FFT_Processor_Spqlios(const int32_t N)
    : _2N(2*N), N(N), Ns2(N/2) {
    auto *tables = build_stockham_tables(N);
    tables_direct = tables;
    tables_reverse = tables;
    real_inout_direct = tables->data;
    reva = new int32_t[Ns2];
    cosomegaxminus1 = new double[2*_2N];
    sinomegaxminus1 = cosomegaxminus1 + _2N;
    int32_t r1 = rev(1, _2N), r3 = rev(3, _2N);
    for (int32_t ri = r1; ri < r3; ri++) reva[ri-r1] = rev(ri, _2N);
    for (int32_t j = 0; j < _2N; j++) {
        cosomegaxminus1[j] = cos(2*M_PI*j/_2N) - 1.;
        sinomegaxminus1[j] = sin(2*M_PI*j/_2N);
    }
}

FFT_Processor_Spqlios::~FFT_Processor_Spqlios() {
    auto *tables = (STOCKHAM_PRECOMP *)tables_direct;
    free(tables->buf);
    delete tables;
    delete[] cosomegaxminus1;
    delete[] reva;
}

void FFT_Processor_Spqlios::execute_reverse_torus32(double *res, const uint32_t *a) {
    const int32_t *aa = (const int32_t*)a;
    for (int32_t i = 0; i < N; i++) res[i] = (double)aa[i];
    split_ifft((const STOCKHAM_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios::execute_reverse_int(double *res, const int32_t *a) {
    for (int32_t i = 0; i < N; i++) res[i] = (double)a[i];
    split_ifft((const STOCKHAM_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios::execute_reverse_uint(double *res, const uint32_t *a) {
    for (int32_t i = 0; i < N; i++) res[i] = (double)a[i];
    split_ifft((const STOCKHAM_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios::execute_reverse_torus64(double *res, const uint64_t *a) {
    const int64_t *aa = (const int64_t*)a;
    for (int32_t i = 0; i < N; i++) res[i] = (double)aa[i];
    split_ifft((const STOCKHAM_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios::execute_reverse_torus64_uint(double *res, const uint64_t *a) {
    for (int32_t i = 0; i < N; i++) res[i] = (double)a[i];
    split_ifft((const STOCKHAM_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios::execute_direct_torus32(uint32_t *res, const double *a) {
    const double s = 2.0 / N;
    for (int32_t i = 0; i < N; i++) real_inout_direct[i] = a[i] * s;
    split_fft((const STOCKHAM_PRECOMP*)tables_direct, real_inout_direct);
    for (int32_t i = 0; i < N; i++)
        res[i] = (uint32_t)(int64_t)real_inout_direct[i];
}

void FFT_Processor_Spqlios::execute_direct_torus32_add(uint32_t *res, const double *a) {
    const double s = 2.0 / N;
    for (int32_t i = 0; i < N; i++) real_inout_direct[i] = a[i] * s;
    split_fft((const STOCKHAM_PRECOMP*)tables_direct, real_inout_direct);
    for (int32_t i = 0; i < N; i++)
        res[i] += (uint32_t)(int64_t)real_inout_direct[i];
}

void FFT_Processor_Spqlios::execute_direct_torus64(uint64_t *res, double *a) {
    const double s = 2.0 / N;
    for (int32_t i = 0; i < N; i++) a[i] *= s;
    split_fft((const STOCKHAM_PRECOMP*)tables_direct, a);
    const double magic = 6755399441055744.0;
    const int64_t magic_i = 0x4338000000000000LL;
    for (int32_t i = 0; i < N; i++) {
        union { double d; int64_t l; } u;
        u.d = a[i] + magic; res[i] = (uint64_t)(u.l - magic_i);
    }
}

void FFT_Processor_Spqlios::execute_direct_torus64_add(uint64_t *res, double *a) {
    const double s = 2.0 / N;
    for (int32_t i = 0; i < N; i++) a[i] *= s;
    split_fft((const STOCKHAM_PRECOMP*)tables_direct, a);
    const double magic = 6755399441055744.0;
    const int64_t magic_i = 0x4338000000000000LL;
    for (int32_t i = 0; i < N; i++) {
        union { double d; int64_t l; } u;
        u.d = a[i] + magic; res[i] += (uint64_t)(u.l - magic_i);
    }
}

void FFT_Processor_Spqlios::execute_direct_torus32_q(uint32_t *res, const double *a, const uint32_t q) {
    const double s = 2.0 / N;
    for (int32_t i = 0; i < N; i++) real_inout_direct[i] = a[i] * s;
    split_fft((const STOCKHAM_PRECOMP*)tables_direct, real_inout_direct);
    for (int32_t i = 0; i < N; i++)
        res[i] = uint32_t((int64_t(real_inout_direct[i])%q+q)%q);
}

void FFT_Processor_Spqlios::execute_direct_torus32_rescale(uint32_t *res, const double *a, const double D) {
    const double s = 2.0 / N;
    for (int32_t i = 0; i < N; i++) real_inout_direct[i] = a[i] * s;
    split_fft((const STOCKHAM_PRECOMP*)tables_direct, real_inout_direct);
    for (int32_t i = 0; i < N; i++)
        res[i] = (uint32_t)(int64_t)(real_inout_direct[i]/D);
}

void FFT_Processor_Spqlios::execute_direct_torus32_rescale_clpx(
    uint32_t *res, const double *a, const double q, const uint32_t plain_modulus) {
    execute_direct_torus32(res, a);
}

void FFT_Processor_Spqlios::execute_direct_torus64_rescale(uint64_t *res, const double *a, const double D) {
    const double s = 2.0 / N;
    alignas(64) double tmp[N];
    for (int32_t i = 0; i < N; i++) tmp[i] = a[i] * s;
    split_fft((const STOCKHAM_PRECOMP*)tables_direct, tmp);
    for (int32_t i = 0; i < N; i++)
        res[i] = (uint64_t)std::round(tmp[i]/D);
}

void FFT_Processor_Spqlios::execute_direct_torus64_rescale_clpx(
    uint64_t *res, const double *a, const uint32_t plain_modulus) {
    execute_direct_torus64(res, const_cast<double*>(a));
}

thread_local FFT_Processor_Spqlios fftplvl1(TFHEpp::lvl1param::n);
thread_local FFT_Processor_Spqlios fftplvl2(TFHEpp::lvl2param::n);
thread_local FFT_Processor_Spqlios fftplvl3(TFHEpp::lvl3param::n);
thread_local FFT_Processor_Spqlios fftplvl5(1 << 14);
