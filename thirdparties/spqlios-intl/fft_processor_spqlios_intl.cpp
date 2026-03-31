// Interleaved-format FFT using Stockham radix-4 algorithm.
//
// Data layout: [re0, im0, re1, im1, ...] (N doubles = N/2 complex values)
// Each YMM register holds 2 complex values (c64x2).
//
// Based on the Stockham auto-sort algorithm (same approach as OTFFT/tfhe-rs):
// - Out-of-place: alternates between data and scratch buffers
// - Radix-4: processes 4 inputs per butterfly, halving the number of passes
// - Mixed radix: uses radix-2 final pass when n is not a power of 4
// - mul_j optimization: multiplication by ±j is free (swap+negate)

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <immintrin.h>

#include <params.hpp>
#include "fft_processor_spqlios_intl.h"

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

#ifdef USE_AVX512

// ── Complex arithmetic with AVX512 ──────────────────────────────────────────
// Each ZMM holds 4 complex values: [re0,im0, re1,im1, re2,im2, re3,im3]

static inline __m512d cmul512(__m512d a, __m512d w) {
    __m512d w_swap = _mm512_permute_pd(w, 0x55);        // swap re/im pairs
    __m512d a_re = _mm512_unpacklo_pd(a, a);             // broadcast re
    __m512d a_im = _mm512_unpackhi_pd(a, a);             // broadcast im
    return _mm512_fmaddsub_pd(a_re, w, _mm512_mul_pd(a_im, w_swap));
}

// Multiply by j: [re,im] → [-im,re]
static inline __m512d mul_j_fwd512(__m512d x) {
    __m512d swapped = _mm512_permute_pd(x, 0x55);
    return _mm512_xor_pd(swapped,
        _mm512_set_pd(0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0));
}

// Multiply by -j: [re,im] → [im,-re]
static inline __m512d mul_j_inv512(__m512d x) {
    __m512d swapped = _mm512_permute_pd(x, 0x55);
    return _mm512_xor_pd(swapped,
        _mm512_set_pd(-0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0));
}

// Radix-4 butterfly with twiddle (4 complex per ZMM) — templated on direction
template<bool Fwd>
static inline void radix4_dit_butterfly512(
    __m512d a, __m512d b, __m512d c, __m512d d,
    __m512d w1, __m512d w2, __m512d w3,
    __m512d &out0, __m512d &out1, __m512d &out2, __m512d &out3)
{
    __m512d apc = _mm512_add_pd(a, c);
    __m512d amc = _mm512_sub_pd(a, c);
    __m512d bpd = _mm512_add_pd(b, d);
    __m512d bmd = _mm512_sub_pd(b, d);
    __m512d jbmd = Fwd ? mul_j_fwd512(bmd) : mul_j_inv512(bmd);

    out0 = _mm512_add_pd(apc, bpd);
    out1 = cmul512(_mm512_sub_pd(amc, jbmd), w1);
    out2 = cmul512(_mm512_sub_pd(apc, bpd), w2);
    out3 = cmul512(_mm512_add_pd(amc, jbmd), w3);
}

// Last radix-4 butterfly (no twiddle) — templated on direction
template<bool Fwd>
static inline void radix4_last_butterfly512(
    __m512d a, __m512d b, __m512d c, __m512d d,
    __m512d &out0, __m512d &out1, __m512d &out2, __m512d &out3)
{
    __m512d apc = _mm512_add_pd(a, c);
    __m512d amc = _mm512_sub_pd(a, c);
    __m512d bpd = _mm512_add_pd(b, d);
    __m512d bmd = _mm512_sub_pd(b, d);
    __m512d jbmd = Fwd ? mul_j_fwd512(bmd) : mul_j_inv512(bmd);

    out0 = _mm512_add_pd(apc, bpd);
    out1 = _mm512_sub_pd(amc, jbmd);
    out2 = _mm512_sub_pd(apc, bpd);
    out3 = _mm512_add_pd(amc, jbmd);
}

#else  // AVX2

// ── Complex arithmetic with AVX2 ────────────────────────────────────────────

// Complex multiply: a * w
static inline __m256d cmul(__m256d a, __m256d w) {
    __m256d w_swap = _mm256_permute_pd(w, 0b0101);
    __m256d a_re = _mm256_unpacklo_pd(a, a);
    __m256d a_im = _mm256_unpackhi_pd(a, a);
    return _mm256_fmaddsub_pd(a_re, w, _mm256_mul_pd(a_im, w_swap));
}

// Multiply by j (imaginary unit): [re,im] → [-im,re]
static inline __m256d mul_j_fwd(__m256d x) {
    __m256d swapped = _mm256_permute_pd(x, 0b0101);
    return _mm256_xor_pd(swapped, _mm256_set_pd(0.0, -0.0, 0.0, -0.0));
}

// Multiply by -j: [re,im] → [im,-re]
static inline __m256d mul_j_inv(__m256d x) {
    __m256d swapped = _mm256_permute_pd(x, 0b0101);
    return _mm256_xor_pd(swapped, _mm256_set_pd(-0.0, 0.0, -0.0, 0.0));
}

// ── Stockham radix-4 DIF butterfly ──────────────────────────────────────────
// Processes 2 complex values per YMM.
// Input: a,b,c,d from 4 quarters of the source array
// Output: 4 results with twiddle factors applied
// w1,w2,w3 are broadcast twiddles (same for both complex values in YMM)

static inline void radix4_dit_butterfly(
    __m256d a, __m256d b, __m256d c, __m256d d,
    __m256d w1, __m256d w2, __m256d w3, bool fwd,
    __m256d &out0, __m256d &out1, __m256d &out2, __m256d &out3)
{
    __m256d apc = _mm256_add_pd(a, c);
    __m256d amc = _mm256_sub_pd(a, c);
    __m256d bpd = _mm256_add_pd(b, d);
    __m256d bmd = _mm256_sub_pd(b, d);
    __m256d jbmd = fwd ? mul_j_fwd(bmd) : mul_j_inv(bmd);

    out0 = _mm256_add_pd(apc, bpd);
    out1 = cmul(_mm256_sub_pd(amc, jbmd), w1);
    out2 = cmul(_mm256_sub_pd(apc, bpd), w2);
    out3 = cmul(_mm256_add_pd(amc, jbmd), w3);
}

// Last radix-4 butterfly (no twiddle) — separate fwd/inv to eliminate branch
static inline void radix4_last_butterfly_fwd(
    __m256d a, __m256d b, __m256d c, __m256d d,
    __m256d &out0, __m256d &out1, __m256d &out2, __m256d &out3)
{
    __m256d apc = _mm256_add_pd(a, c);
    __m256d amc = _mm256_sub_pd(a, c);
    __m256d bpd = _mm256_add_pd(b, d);
    __m256d jbmd = mul_j_fwd(_mm256_sub_pd(b, d));
    out0 = _mm256_add_pd(apc, bpd);
    out1 = _mm256_sub_pd(amc, jbmd);
    out2 = _mm256_sub_pd(apc, bpd);
    out3 = _mm256_add_pd(amc, jbmd);
}

static inline void radix4_last_butterfly_inv(
    __m256d a, __m256d b, __m256d c, __m256d d,
    __m256d &out0, __m256d &out1, __m256d &out2, __m256d &out3)
{
    __m256d apc = _mm256_add_pd(a, c);
    __m256d amc = _mm256_sub_pd(a, c);
    __m256d bpd = _mm256_add_pd(b, d);
    __m256d jbmd = mul_j_inv(_mm256_sub_pd(b, d));
    out0 = _mm256_add_pd(apc, bpd);
    out1 = _mm256_sub_pd(amc, jbmd);
    out2 = _mm256_sub_pd(apc, bpd);
    out3 = _mm256_add_pd(amc, jbmd);
}

// Separate fwd/inv radix-4 DIF butterfly with twiddle
static inline void radix4_butterfly_fwd(
    __m256d a, __m256d b, __m256d c, __m256d d,
    __m256d w1, __m256d w2, __m256d w3,
    __m256d &out0, __m256d &out1, __m256d &out2, __m256d &out3)
{
    __m256d apc = _mm256_add_pd(a, c);
    __m256d amc = _mm256_sub_pd(a, c);
    __m256d bpd = _mm256_add_pd(b, d);
    __m256d jbmd = mul_j_fwd(_mm256_sub_pd(b, d));
    out0 = _mm256_add_pd(apc, bpd);
    out1 = cmul(_mm256_sub_pd(amc, jbmd), w1);
    out2 = cmul(_mm256_sub_pd(apc, bpd), w2);
    out3 = cmul(_mm256_add_pd(amc, jbmd), w3);
}

static inline void radix4_butterfly_inv(
    __m256d a, __m256d b, __m256d c, __m256d d,
    __m256d w1, __m256d w2, __m256d w3,
    __m256d &out0, __m256d &out1, __m256d &out2, __m256d &out3)
{
    __m256d apc = _mm256_add_pd(a, c);
    __m256d amc = _mm256_sub_pd(a, c);
    __m256d bpd = _mm256_add_pd(b, d);
    __m256d jbmd = mul_j_inv(_mm256_sub_pd(b, d));
    out0 = _mm256_add_pd(apc, bpd);
    out1 = cmul(_mm256_sub_pd(amc, jbmd), w1);
    out2 = cmul(_mm256_sub_pd(apc, bpd), w2);
    out3 = cmul(_mm256_add_pd(amc, jbmd), w3);
}

#endif  // USE_AVX512

// ── Table structures ────────────────────────────────────────────────────────

struct INTL_FFT_PRECOMP {
    int32_t n;             // = 2*N
    int32_t ns2;           // = N/2 = number of complex points
    double *trig_fwd;      // forward twiddles + twist
    double *trig_inv;      // inverse twiddles + twist
    double *scratch;       // scratch buffer for out-of-place FFT
    double *data;          // data buffer for execute_direct
    void *buf;             // single allocation
};

#ifdef USE_AVX512

// ── AVX512 interleave/de-interleave index vectors ──────────────────────────
// Used by conversion loops: 8 re + 8 im → [re0,im0,...,re7,im7] (2 ZMMs)
static const __m512i idx_intl_lo = _mm512_set_epi64(11, 3, 10, 2,  9, 1,  8, 0);
static const __m512i idx_intl_hi = _mm512_set_epi64(15, 7, 14, 6, 13, 5, 12, 4);
// Reverse: [re0,im0,...] (2 ZMMs) → 8 re, 8 im
static const __m512i idx_deinl_re = _mm512_set_epi64(14, 12, 10, 8, 6, 4, 2, 0);
static const __m512i idx_deinl_im = _mm512_set_epi64(15, 13, 11, 9, 7, 5, 3, 1);

// ── Stockham radix-4 FFT (AVX512) ──────────────────────────────────────────
// Templated on direction to eliminate the runtime fwd/inv branch.
// Hand-tuned asm inner loop: w1/w2/w3 preloaded in zmm0-2, jmask in zmm3.
// zmm4-zmm11 used as temporaries; zmm12-zmm31 remain free for OOO scheduling.

template<bool Fwd>
static void stockham_r4_from(int32_t ns2, const double *trig,
                             const double *src_ro, double *x, double *y) {
    const double *src = src_ro;
    double *dst = y;
    const double *tw = trig;
    int32_t q = ns2 / 4;
    int32_t s = 1;

    // j-multiply sign mask, compile-time selected:
    //   Fwd: swap re/im then negate even (re) positions → multiply by +j
    //   Inv: swap re/im then negate odd  (im) positions → multiply by -j
    const __m512d jmask = Fwd ?
        _mm512_set_pd(0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0) :
        _mm512_set_pd(-0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0);

    while (q >= 4) {
        int32_t stride = q * s;
        if (s == 1) {
            // First pass: no twiddle multiply; process 4 consecutive p values per ZMM.
            for (int32_t p = 0; p < q; p += 4) {
                __m512d a = _mm512_loadu_pd(src + p * 2);
                __m512d b = _mm512_loadu_pd(src + (p + stride) * 2);
                __m512d c = _mm512_loadu_pd(src + (p + 2*stride) * 2);
                __m512d d = _mm512_loadu_pd(src + (p + 3*stride) * 2);
                __m512d r0, r1, r2, r3;
                radix4_last_butterfly512<Fwd>(a, b, c, d, r0, r1, r2, r3);
                _mm512_storeu_pd(dst + p * 2, r0);
                _mm512_storeu_pd(dst + (p + q) * 2, r1);
                _mm512_storeu_pd(dst + (p + 2*q) * 2, r2);
                _mm512_storeu_pd(dst + (p + 3*q) * 2, r3);
            }
        } else {
            // Subsequent passes: preload w1/w2/w3 into zmm0-2 and jmask into zmm3,
            // then run hand-tuned asm p-loop to guarantee no ZMM register spills.
            int64_t stride_bytes = (int64_t)stride * 16;  // stride × 2 doubles × 8 bytes
            for (int32_t j = 0; j < s; j += 4) {
                int64_t q_bytes = (int64_t)q * 16;  // q × 2 doubles × 8 bytes

                // Preload twiddles and jmask before the p-loop.
                // zmm0-3 clobbered; GCC will not assign C variables to them
                // between this block and the inner asm, leaving our values intact.
                __asm__ __volatile__ (
                    "vmovupd    (%[tw]),   %%zmm0\n\t"  // w1: 4 complex twiddles
                    "vmovupd  64(%[tw]),   %%zmm1\n\t"  // w2
                    "vmovupd 128(%[tw]),   %%zmm2\n\t"  // w3
                    "vmovapd    %[jm],     %%zmm3\n\t"  // jmask (compile-time constant)
                    : : [tw] "r"(tw), [jm] "x"(jmask)
                    : "zmm0","zmm1","zmm2","zmm3"
                );

                for (int32_t p = 0; p < q; p++) {
                    const double *sptr = src + (int64_t)(j + s * p) * 2;
                    double *dptr = dst + (int64_t)(p + q * (4 * j)) * 2;

                    // zmm0=w1, zmm1=w2, zmm2=w3, zmm3=jmask (preserved across p-loop).
                    // zmm4-zmm11: temporaries; zmm12-zmm31: untouched.
                    __asm__ __volatile__ (
                        // ── Load a, b, c, d ──────────────────────────────────
                        "vmovupd    (%[src]),           %%zmm4\n\t"  // a
                        "vmovupd    (%[src],%[st]),     %%zmm5\n\t"  // b
                        "vmovupd    (%[src],%[st],2),   %%zmm6\n\t"  // c
                        "vmovupd    (%[src],%[st3]),    %%zmm7\n\t"  // d

                        // ── Radix-4 butterfly sums/diffs ─────────────────────
                        "vaddpd     %%zmm6, %%zmm4, %%zmm8\n\t"   // apc = a+c
                        "vsubpd     %%zmm6, %%zmm4, %%zmm9\n\t"   // amc = a-c
                        "vaddpd     %%zmm7, %%zmm5, %%zmm10\n\t"  // bpd = b+d
                        "vsubpd     %%zmm7, %%zmm5, %%zmm4\n\t"   // bmd = b-d  (reuse zmm4)
                        // j-multiply: swap re/im pairs, then apply sign mask
                        "vpermilpd  $0x55, %%zmm4, %%zmm4\n\t"
                        "vxorpd     %%zmm3, %%zmm4, %%zmm11\n\t"  // jbmd  (zmm3=jmask)
                        // live: apc(8), amc(9), bpd(10), jbmd(11)

                        // ── out0 = apc + bpd  (store immediately) ────────────
                        "vaddpd     %%zmm10, %%zmm8, %%zmm4\n\t"
                        "vmovupd    %%zmm4, (%[dst])\n\t"

                        // ── out2 = (apc − bpd) × w2 ──────────────────────────
                        "vsubpd     %%zmm10, %%zmm8, %%zmm4\n\t"   // t = apc-bpd
                        "vpermilpd  $0x55, %%zmm1, %%zmm5\n\t"     // w2_swap
                        "vunpckhpd  %%zmm4, %%zmm4, %%zmm6\n\t"    // t_im broadcast
                        "vunpcklpd  %%zmm4, %%zmm4, %%zmm4\n\t"    // t_re broadcast
                        "vmulpd     %%zmm5, %%zmm6, %%zmm6\n\t"    // t_im × w2_swap
                        "vfmaddsub231pd %%zmm1, %%zmm4, %%zmm6\n\t"// t_re×w2 ± t_im×w2_swap
                        "vmovupd    %%zmm6, (%[dst],%[q2])\n\t"

                        // ── out1 = (amc − jbmd) × w1 ─────────────────────────
                        "vsubpd     %%zmm11, %%zmm9, %%zmm4\n\t"   // t = amc-jbmd
                        "vpermilpd  $0x55, %%zmm0, %%zmm5\n\t"     // w1_swap
                        "vunpckhpd  %%zmm4, %%zmm4, %%zmm6\n\t"
                        "vunpcklpd  %%zmm4, %%zmm4, %%zmm4\n\t"
                        "vmulpd     %%zmm5, %%zmm6, %%zmm6\n\t"
                        "vfmaddsub231pd %%zmm0, %%zmm4, %%zmm6\n\t"
                        "vmovupd    %%zmm6, (%[dst],%[q1])\n\t"

                        // ── out3 = (amc + jbmd) × w3 ─────────────────────────
                        "vaddpd     %%zmm11, %%zmm9, %%zmm4\n\t"
                        "vpermilpd  $0x55, %%zmm2, %%zmm5\n\t"     // w3_swap
                        "vunpckhpd  %%zmm4, %%zmm4, %%zmm6\n\t"
                        "vunpcklpd  %%zmm4, %%zmm4, %%zmm4\n\t"
                        "vmulpd     %%zmm5, %%zmm6, %%zmm6\n\t"
                        "vfmaddsub231pd %%zmm2, %%zmm4, %%zmm6\n\t"
                        "vmovupd    %%zmm6, (%[dst],%[q3])\n\t"

                        : : [src] "r"(sptr), [dst] "r"(dptr),
                            [st]  "r"(stride_bytes),
                            [st3] "r"(stride_bytes * 3),
                            [q1]  "r"(q_bytes),
                            [q2]  "r"(2 * q_bytes),
                            [q3]  "r"(3 * q_bytes)
                        : "zmm4","zmm5","zmm6","zmm7",
                          "zmm8","zmm9","zmm10","zmm11","memory"
                    );
                }
                tw += 24;  // advance past this j-group (3 ZMMs × 8 doubles)
            }
        }
        // Switch ping-pong buffers
        if (s == 1) { src = dst; dst = x; }
        else { const double *tmp = src; src = dst; dst = const_cast<double*>(tmp); }
        s *= 4; q /= 4;
    }

    // Final q=1 or q=2 pass (no twiddle multiply)
    if (q == 1) {
        int32_t stride = s;
        for (int32_t j = 0; j < s; j += 4) {
            __m512d a = _mm512_loadu_pd(src + j * 2);
            __m512d b = _mm512_loadu_pd(src + (j + stride) * 2);
            __m512d c = _mm512_loadu_pd(src + (j + 2*stride) * 2);
            __m512d d = _mm512_loadu_pd(src + (j + 3*stride) * 2);
            __m512d r0, r1, r2, r3;
            radix4_last_butterfly512<Fwd>(a, b, c, d, r0, r1, r2, r3);
            _mm512_storeu_pd(dst + (4*j) * 2, r0);
            _mm512_storeu_pd(dst + (4*j + 4) * 2, r1);
            _mm512_storeu_pd(dst + (4*j + 8) * 2, r2);
            _mm512_storeu_pd(dst + (4*j + 12) * 2, r3);
        }
    } else if (q == 2) {
        int32_t half = ns2 / 2;
        for (int32_t j = 0; j < half; j += 4) {
            __m512d a = _mm512_loadu_pd(src + j * 2);
            __m512d b = _mm512_loadu_pd(src + (j + half) * 2);
            __m512d sum = _mm512_add_pd(a, b);
            __m512d diff = _mm512_sub_pd(a, b);
            __m512d out0 = _mm512_shuffle_f64x2(sum, diff, 0x44);
            __m512d out1 = _mm512_shuffle_f64x2(sum, diff, 0xEE);
            _mm512_storeu_pd(dst + (2*j) * 2, out0);
            _mm512_storeu_pd(dst + (2*j + 4) * 2, out1);
        }
    }
    // Ensure result ends up in x
    if (dst != x) memcpy(x, dst, ns2 * 2 * sizeof(double));
}

// Convenience wrapper: in-place (src == x)
template<bool Fwd>
static void stockham_r4(int32_t ns2, const double *trig, double *x, double *y) {
    stockham_r4_from<Fwd>(ns2, trig, x, x, y);
}

#else  // AVX2

// ── Stockham radix-4 DIF FFT (AVX2) ────────────────────────────────────────
// Verified correct algorithm:
//   Read:  src[p + s*(j + m*q)]  for m=0..3
//   BF:    standard radix-4 DIF
//   Twid:  exp(sign*2*pi*i * m * j * s / ns2)  applied AFTER butterfly to outputs m=1,2,3
//   Write: dst[p + s*(m + 4*j)]  (auto-sort)
//   Update: s *= 4, q /= 4
//
// Twiddle table layout (paired):
//   For each pass, for pair k (j=2k, j=2k+1):
//     w1: [cos(2k*angle), sin(2k*angle), cos((2k+1)*angle), sin((2k+1)*angle)]
//     w2: same with 2*angle
//     w3: same with 3*angle
//   Total per pass: ceil(q/2) * 12 doubles
//
// For s=1: vectorize over j (step 2), 128-bit output stores.
// For s>=4: vectorize over p (step 2), broadcast twiddle from paired table.

template<bool Fwd>
static void stockham_r4_avx2(int32_t ns2, const double *bf_twiddles,
                              double *x, double *scratch) {
    double *src = x, *dst = scratch;
    const double *tw = bf_twiddles;
    int32_t q = ns2 / 4, s = 1;

    // Radix-4 passes
    while (q >= 1) {
        int64_t arm_stride = (int64_t)s * q;  // = ns2/4

        if (s == 1) {
            // ── First pass (s=1): vectorize over j, step 2 ──────────────
            // Each YMM holds complex values for (j, j+1).
            // Source: src[j + m*q] adjacent → YMM load.
            // Output: dst[m + 4*j] and dst[m + 4*(j+1)] are 4 apart → 128-bit stores.
            for (int32_t j = 0; j < q; j += 2) {
                __m256d a = _mm256_loadu_pd(src + (int64_t)j * 2);
                __m256d b = _mm256_loadu_pd(src + (int64_t)(j + q) * 2);
                __m256d c = _mm256_loadu_pd(src + (int64_t)(j + 2*q) * 2);
                __m256d d = _mm256_loadu_pd(src + (int64_t)(j + 3*q) * 2);

                __m256d apc = _mm256_add_pd(a, c);
                __m256d amc = _mm256_sub_pd(a, c);
                __m256d bpd = _mm256_add_pd(b, d);
                __m256d bmd = _mm256_sub_pd(b, d);
                __m256d jbmd = Fwd ? mul_j_fwd(bmd) : mul_j_inv(bmd);

                __m256d r0 = _mm256_add_pd(apc, bpd);
                __m256d r1 = _mm256_sub_pd(amc, jbmd);
                __m256d r2 = _mm256_sub_pd(apc, bpd);
                __m256d r3 = _mm256_add_pd(amc, jbmd);

                // Load paired twiddle for (j, j+1) and apply.
                // j=0 twiddle is (1,0) stored in the table, so multiply is harmless.
                int32_t pair_k = j / 2;
                __m256d tw1 = _mm256_loadu_pd(tw + pair_k * 12);
                __m256d tw2 = _mm256_loadu_pd(tw + pair_k * 12 + 4);
                __m256d tw3 = _mm256_loadu_pd(tw + pair_k * 12 + 8);
                r1 = cmul(r1, tw1);
                r2 = cmul(r2, tw2);
                r3 = cmul(r3, tw3);

                // Output: dst[m + 4*j] for j, dst[m + 4*(j+1)] for j+1
                int64_t out_j0 = (int64_t)(4 * j) * 2;
                int64_t out_j1 = (int64_t)(4 * (j+1)) * 2;
                _mm_storeu_pd(dst + out_j0,     _mm256_castpd256_pd128(r0));
                _mm_storeu_pd(dst + out_j1,     _mm256_extractf128_pd(r0, 1));
                _mm_storeu_pd(dst + out_j0 + 2, _mm256_castpd256_pd128(r1));
                _mm_storeu_pd(dst + out_j1 + 2, _mm256_extractf128_pd(r1, 1));
                _mm_storeu_pd(dst + out_j0 + 4, _mm256_castpd256_pd128(r2));
                _mm_storeu_pd(dst + out_j1 + 4, _mm256_extractf128_pd(r2, 1));
                _mm_storeu_pd(dst + out_j0 + 6, _mm256_castpd256_pd128(r3));
                _mm_storeu_pd(dst + out_j1 + 6, _mm256_extractf128_pd(r3, 1));
            }
            tw += (q / 2) * 12;  // q/2 pairs
        } else {
            // ── General passes (s>=4): vectorize over p, step 2 ─────────
            // Each YMM holds (p, p+1). Twiddle broadcast from paired table.
            for (int32_t j = 0; j < q; j++) {
                __m256d w1, w2, w3;
                bool has_twiddle = (j != 0);
                if (has_twiddle) {
                    // Twiddles in paired format. Extract the lane for this j.
                    int32_t pair_k = j / 2;
                    const double *twp = tw + pair_k * 12;
                    __m256d p1 = _mm256_loadu_pd(twp);
                    __m256d p2 = _mm256_loadu_pd(twp + 4);
                    __m256d p3 = _mm256_loadu_pd(twp + 8);
                    if (j % 2 == 0) {
                        w1 = _mm256_permute2f128_pd(p1, p1, 0x00);
                        w2 = _mm256_permute2f128_pd(p2, p2, 0x00);
                        w3 = _mm256_permute2f128_pd(p3, p3, 0x00);
                    } else {
                        w1 = _mm256_permute2f128_pd(p1, p1, 0x11);
                        w2 = _mm256_permute2f128_pd(p2, p2, 0x11);
                        w3 = _mm256_permute2f128_pd(p3, p3, 0x11);
                    }
                }

                int64_t base_in = (int64_t)s * j;
                int64_t base_out = (int64_t)s * 4 * j;
                int64_t stride_s = (int64_t)s;

                for (int32_t p = 0; p < s; p += 2) {
                    const double *sp = src + (base_in + p) * 2;
                    __m256d a = _mm256_loadu_pd(sp);
                    __m256d b = _mm256_loadu_pd(sp + arm_stride * 2);
                    __m256d c = _mm256_loadu_pd(sp + arm_stride * 4);
                    __m256d d = _mm256_loadu_pd(sp + arm_stride * 6);

                    __m256d apc = _mm256_add_pd(a, c);
                    __m256d amc = _mm256_sub_pd(a, c);
                    __m256d bpd = _mm256_add_pd(b, d);
                    __m256d bmd = _mm256_sub_pd(b, d);
                    __m256d jbmd = Fwd ? mul_j_fwd(bmd) : mul_j_inv(bmd);

                    __m256d r0 = _mm256_add_pd(apc, bpd);
                    __m256d r1 = _mm256_sub_pd(amc, jbmd);
                    __m256d r2 = _mm256_sub_pd(apc, bpd);
                    __m256d r3 = _mm256_add_pd(amc, jbmd);

                    if (has_twiddle) {
                        r1 = cmul(r1, w1);
                        r2 = cmul(r2, w2);
                        r3 = cmul(r3, w3);
                    }

                    double *dp = dst + (base_out + p) * 2;
                    _mm256_storeu_pd(dp,                    r0);
                    _mm256_storeu_pd(dp + stride_s * 2,     r1);
                    _mm256_storeu_pd(dp + stride_s * 4,     r2);
                    _mm256_storeu_pd(dp + stride_s * 6,     r3);
                }
            }
            tw += ((q + 1) / 2) * 12;  // ceil(q/2) pairs
        }

        double *tmp = src; src = dst; dst = tmp;
        s *= 4; q /= 4;
    }

    // Radix-2 final pass if ns2 is 2*4^k (e.g., 512 = 2*256)
    if (s < ns2) {
        for (int32_t p = 0; p < s; p += 2) {
            __m256d a = _mm256_loadu_pd(src + p * 2);
            __m256d b = _mm256_loadu_pd(src + (p + s) * 2);
            _mm256_storeu_pd(dst + p * 2, _mm256_add_pd(a, b));
            _mm256_storeu_pd(dst + (p + s) * 2, _mm256_sub_pd(a, b));
        }
        double *tmp = src; src = dst; dst = tmp;
    }

    if (src != x) memcpy(x, src, ns2 * 2 * sizeof(double));
}

#endif  // USE_AVX512

// ── Forward/Inverse FFT wrappers ────────────────────────────────────────────

// Forward FFT: Stockham + twist.  Reads from src_in (may be != c).
// Scale factor (2/N) is pre-baked into the forward twist twiddles.
static void intl_fft_from(const INTL_FFT_PRECOMP *tables, const double *src_in,
                          double *c) {
    const int32_t ns2 = tables->ns2;
    const double *trig = tables->trig_fwd;

#ifdef USE_AVX512
    stockham_r4_from<true>(ns2, trig, src_in, c, tables->scratch);

    const double *tw_twist = trig;
    for (int32_t s = 4; 4*s < ns2; s *= 4)
        tw_twist += (s / 4) * 24;

    for (int32_t j = 0; j < ns2; j += 4) {
        double *p = c + j * 2;
        __m512d a = _mm512_load_pd(p);
        __m512d w = _mm512_load_pd(tw_twist + j * 2);
        _mm512_store_pd(p, cmul512(a, w));
    }
#else
    // AVX2: copy input, run Stockham, then twist
    if (src_in != c) memcpy(c, src_in, ns2 * 2 * sizeof(double));
    stockham_r4_avx2<true>(ns2, trig, c, tables->scratch);

    // Advance past butterfly twiddles to reach the twist table
    const double *tw_twist = trig;
    {
        int32_t ss = 1, qq = ns2 / 4;
        while (qq >= 1) {
            if (ss == 1)
                tw_twist += (qq / 2) * 12;         // q/2 pairs
            else
                tw_twist += ((qq + 1) / 2) * 12;   // ceil(q/2) pairs
            ss *= 4; qq /= 4;
        }
    }

    for (int32_t j = 0; j < ns2; j += 2) {
        double *p = c + j * 2;
        __m256d a = _mm256_load_pd(p);
        __m256d w = _mm256_load_pd(tw_twist + j * 2);
        _mm256_store_pd(p, cmul(a, w));
    }
#endif
}

static void intl_fft(const INTL_FFT_PRECOMP *tables, double *c) {
    intl_fft_from(tables, c, c);
}

static void intl_ifft(const INTL_FFT_PRECOMP *tables, double *c) {
    const int32_t ns2 = tables->ns2;
    const double *trig = tables->trig_inv;
    const double *tw_twist = trig;

#ifdef USE_AVX512
    // Apply twist first
    for (int32_t j = 0; j < ns2; j += 4) {
        double *p = c + j * 2;
        __m512d a = _mm512_load_pd(p);
        __m512d w = _mm512_load_pd(tw_twist + j * 2);
        _mm512_store_pd(p, cmul512(a, w));
    }
    const double *tw_bf = trig + ns2 * 2;
    stockham_r4<false>(ns2, tw_bf, c, tables->scratch);
#else
    // AVX2: Separate twist + Stockham inverse
    for (int32_t j = 0; j < ns2; j += 2) {
        double *p = c + j * 2;
        __m256d a = _mm256_load_pd(p);
        __m256d w = _mm256_load_pd(tw_twist + j * 2);
        _mm256_store_pd(p, cmul(a, w));
    }

    const double *tw_bf = trig + ns2 * 2;
    stockham_r4_avx2<false>(ns2, tw_bf, c, tables->scratch);
#endif
}

// ── Table construction ──────────────────────────────────────────────────────

static void build_tables(int32_t nn, INTL_FFT_PRECOMP *reps) {
    int32_t n = 2 * nn;
    int32_t ns2 = nn / 2;
    reps->n = n;
    reps->ns2 = ns2;

#ifdef USE_AVX512
    // AVX512: twiddles packed as 4 complex (8 doubles) per ZMM.
    // First pass (s=1) uses no twiddles.
    // Subsequent passes (s=4,16,...): j steps by 4, s/4 groups per pass, 3 ZMMs each = 24 doubles.
    int32_t bf_twiddle_doubles = 0;
    for (int32_t s = 4; 4*s < ns2; s *= 4) {
        bf_twiddle_doubles += (s / 4) * 24;
    }
#else
    // AVX2: paired twiddle table. For each radix-4 pass:
    //   s=1:    q/2 pairs * 12 doubles
    //   s>=4:   ceil(q/2) pairs * 12 doubles
    int32_t bf_twiddle_doubles = 0;
    {
        int32_t ss = 1, qq = ns2 / 4;
        while (qq >= 1) {
            if (ss == 1)
                bf_twiddle_doubles += (qq / 2) * 12;
            else
                bf_twiddle_doubles += ((qq + 1) / 2) * 12;
            ss *= 4; qq /= 4;
        }
    }
#endif

    int32_t twist_doubles = ns2 * 2;
    int32_t total_per_dir = bf_twiddle_doubles + twist_doubles;

    int32_t total_doubles = 2 * total_per_dir + nn + nn;
    reps->buf = aligned_alloc(64, total_doubles * sizeof(double));
    double *ptr = (double *)reps->buf;

    reps->trig_fwd = ptr; ptr += total_per_dir;
    reps->trig_inv = ptr; ptr += total_per_dir;
    reps->scratch = ptr; ptr += nn;
    reps->data = ptr;

#ifdef USE_AVX512
    // Build forward butterfly twiddles (AVX512: 4 complex per ZMM)
    double *fwd = reps->trig_fwd;
    for (int32_t s = 4; 4*s < ns2; s *= 4) {
        int32_t denom = 4 * s;
        for (int32_t j = 0; j < s; j += 4) {
            // w1: exp(-2πi*(j+jj)/(4s)) for jj=0..3
            for (int32_t jj = 0; jj < 4; jj++) {
                fwd[jj*2]   = accurate_cos(-(j+jj), denom);
                fwd[jj*2+1] = accurate_sin(-(j+jj), denom);
            }
            fwd += 8;
            // w2: exp(-2πi*2*(j+jj)/(4s))
            for (int32_t jj = 0; jj < 4; jj++) {
                fwd[jj*2]   = accurate_cos(-2*(j+jj), denom);
                fwd[jj*2+1] = accurate_sin(-2*(j+jj), denom);
            }
            fwd += 8;
            // w3: exp(-2πi*3*(j+jj)/(4s))
            for (int32_t jj = 0; jj < 4; jj++) {
                fwd[jj*2]   = accurate_cos(-3*(j+jj), denom);
                fwd[jj*2+1] = accurate_sin(-3*(j+jj), denom);
            }
            fwd += 8;
        }
    }
    // Forward twist with 2/N scale baked in
    {
        const double scale = 2.0 / nn;
        for (int32_t k = 0; k < ns2; k++) {
            *fwd++ = scale * accurate_cos(-k, n);
            *fwd++ = scale * accurate_sin(-k, n);
        }
    }

    // Build inverse twiddles
    double *inv = reps->trig_inv;
    // Inverse twist first
    for (int32_t k = 0; k < ns2; k++) {
        *inv++ = accurate_cos(k, n);
        *inv++ = accurate_sin(k, n);
    }
    // Inverse butterfly twiddles (positive angles)
    for (int32_t s = 4; 4*s < ns2; s *= 4) {
        int32_t denom = 4 * s;
        for (int32_t j = 0; j < s; j += 4) {
            for (int32_t jj = 0; jj < 4; jj++) {
                inv[jj*2]   = accurate_cos(j+jj, denom);
                inv[jj*2+1] = accurate_sin(j+jj, denom);
            }
            inv += 8;
            for (int32_t jj = 0; jj < 4; jj++) {
                inv[jj*2]   = accurate_cos(2*(j+jj), denom);
                inv[jj*2+1] = accurate_sin(2*(j+jj), denom);
            }
            inv += 8;
            for (int32_t jj = 0; jj < 4; jj++) {
                inv[jj*2]   = accurate_cos(3*(j+jj), denom);
                inv[jj*2+1] = accurate_sin(3*(j+jj), denom);
            }
            inv += 8;
        }
    }

#else  // AVX2
    // Build forward butterfly twiddles.
    // Correct twiddle angle: m * j * s / ns2  (in units of full turn)
    // For each pass (s, q=ns2/(4*s)), store paired twiddles for (j=2k, j=2k+1).
    // Forward: negative angle (exp(-2*pi*i * m * j * s / ns2)).
    double *fwd = reps->trig_fwd;
    {
        int32_t ss = 1, qq = ns2 / 4;
        while (qq >= 1) {
            int32_t npairs = (ss == 1) ? (qq / 2) : ((qq + 1) / 2);
            for (int32_t k = 0; k < npairs; k++) {
                int32_t j0 = 2 * k;
                int32_t j1 = 2 * k + 1;
                // angle_j = j * ss (numerator), denom = ns2
                // w_m = exp(-2*pi*i * m * j * ss / ns2) = accurate_cos(-m*j*ss, ns2)
                for (int32_t m = 1; m <= 3; m++) {
                    fwd[0] = accurate_cos(-m * j0 * ss, ns2);
                    fwd[1] = accurate_sin(-m * j0 * ss, ns2);
                    if (j1 < qq) {
                        fwd[2] = accurate_cos(-m * j1 * ss, ns2);
                        fwd[3] = accurate_sin(-m * j1 * ss, ns2);
                    } else {
                        fwd[2] = fwd[0];  // pad with copy of j0
                        fwd[3] = fwd[1];
                    }
                    fwd += 4;
                }
            }
            ss *= 4; qq /= 4;
        }
    }
    // Forward twist with 2/N scale baked in
    {
        const double scale = 2.0 / nn;
        for (int32_t k = 0; k < ns2; k++) {
            *fwd++ = scale * accurate_cos(-k, n);
            *fwd++ = scale * accurate_sin(-k, n);
        }
    }

    // Build inverse twiddles: twist first, then butterfly twiddles.
    double *inv = reps->trig_inv;
    for (int32_t k = 0; k < ns2; k++) {
        *inv++ = accurate_cos(k, n);
        *inv++ = accurate_sin(k, n);
    }
    // Inverse butterfly twiddles: positive angle.
    {
        int32_t ss = 1, qq = ns2 / 4;
        while (qq >= 1) {
            int32_t npairs = (ss == 1) ? (qq / 2) : ((qq + 1) / 2);
            for (int32_t k = 0; k < npairs; k++) {
                int32_t j0 = 2 * k;
                int32_t j1 = 2 * k + 1;
                for (int32_t m = 1; m <= 3; m++) {
                    inv[0] = accurate_cos(m * j0 * ss, ns2);
                    inv[1] = accurate_sin(m * j0 * ss, ns2);
                    if (j1 < qq) {
                        inv[2] = accurate_cos(m * j1 * ss, ns2);
                        inv[3] = accurate_sin(m * j1 * ss, ns2);
                    } else {
                        inv[2] = inv[0];
                        inv[3] = inv[1];
                    }
                    inv += 4;
                }
            }
            ss *= 4; qq /= 4;
        }
    }
#endif  // USE_AVX512
}

// ── Conversion helpers ──────────────────────────────────────────────────────

static inline __m256i magic_cvtpd_epi64(__m256d x) {
    const __m256d m = _mm256_set1_pd(6755399441055744.0);
    const __m256i mi = _mm256_set1_epi64x(0x4338000000000000LL);
    return _mm256_sub_epi64(_mm256_castpd_si256(_mm256_add_pd(x, m)), mi);
}

// ── FFT_Processor implementation ────────────────────────────────────────────

static int32_t rev(int32_t x, int32_t M) {
    int32_t r = 0;
    for (int32_t j = M; j > 1; j /= 2) { r = 2*r+(x%2); x /= 2; }
    return r;
}

FFT_Processor_Spqlios_Intl::FFT_Processor_Spqlios_Intl(const int32_t N)
    : _2N(2*N), N(N), Ns2(N/2) {
    auto *tables = new INTL_FFT_PRECOMP;
    build_tables(N, tables);
    tables_direct = tables;
    tables_reverse = tables;  // same struct, different trig pointers
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

FFT_Processor_Spqlios_Intl::~FFT_Processor_Spqlios_Intl() {
    auto *tables = (INTL_FFT_PRECOMP *)tables_direct;
    free(tables->buf);
    delete tables;
    delete[] cosomegaxminus1;
    delete[] reva;
}

void FFT_Processor_Spqlios_Intl::execute_reverse_torus32(double *res, const uint32_t *a) {
    const int32_t *aa = (const int32_t*)a;
#ifdef USE_AVX512
    // Vectorized int32→double with interleave
    for (int32_t i = 0; i < Ns2; i += 8) {
        __m256i re_i32 = _mm256_loadu_si256((const __m256i*)(aa + i));
        __m256i im_i32 = _mm256_loadu_si256((const __m256i*)(aa + i + Ns2));
        __m512d re = _mm512_cvtepi32_pd(re_i32);
        __m512d im = _mm512_cvtepi32_pd(im_i32);
        _mm512_storeu_pd(res + 2*i,     _mm512_permutex2var_pd(re, idx_intl_lo, im));
        _mm512_storeu_pd(res + 2*i + 8, _mm512_permutex2var_pd(re, idx_intl_hi, im));
    }
#else
    // AVX2: convert 4 int32 → 4 double, interleave re/im with unpack
    for (int32_t i = 0; i < Ns2; i += 4) {
        __m128i re_i32 = _mm_loadu_si128((const __m128i*)(aa + i));
        __m128i im_i32 = _mm_loadu_si128((const __m128i*)(aa + i + Ns2));
        __m256d re = _mm256_cvtepi32_pd(re_i32);
        __m256d im = _mm256_cvtepi32_pd(im_i32);
        // Interleave: [re0,re1,re2,re3] + [im0,im1,im2,im3] → [re0,im0,re1,im1], [re2,im2,re3,im3]
        __m256d lo = _mm256_unpacklo_pd(re, im);  // [re0,im0,re2,im2]
        __m256d hi = _mm256_unpackhi_pd(re, im);  // [re1,im1,re3,im3]
        _mm256_storeu_pd(res + 2*i,     _mm256_permute2f128_pd(lo, hi, 0x20)); // [re0,im0,re1,im1]
        _mm256_storeu_pd(res + 2*i + 4, _mm256_permute2f128_pd(lo, hi, 0x31)); // [re2,im2,re3,im3]
    }
#endif
    intl_ifft((const INTL_FFT_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios_Intl::execute_reverse_int(double *res, const int32_t *a) {
#ifdef USE_AVX512
    for (int32_t i = 0; i < Ns2; i += 8) {
        __m256i re_i32 = _mm256_loadu_si256((const __m256i*)(a + i));
        __m256i im_i32 = _mm256_loadu_si256((const __m256i*)(a + i + Ns2));
        __m512d re = _mm512_cvtepi32_pd(re_i32);
        __m512d im = _mm512_cvtepi32_pd(im_i32);
        _mm512_storeu_pd(res + 2*i,     _mm512_permutex2var_pd(re, idx_intl_lo, im));
        _mm512_storeu_pd(res + 2*i + 8, _mm512_permutex2var_pd(re, idx_intl_hi, im));
    }
#else
    for (int32_t i = 0; i < Ns2; i += 4) {
        __m128i re_i32 = _mm_loadu_si128((const __m128i*)(a + i));
        __m128i im_i32 = _mm_loadu_si128((const __m128i*)(a + i + Ns2));
        __m256d re = _mm256_cvtepi32_pd(re_i32);
        __m256d im = _mm256_cvtepi32_pd(im_i32);
        __m256d lo = _mm256_unpacklo_pd(re, im);
        __m256d hi = _mm256_unpackhi_pd(re, im);
        _mm256_storeu_pd(res + 2*i,     _mm256_permute2f128_pd(lo, hi, 0x20));
        _mm256_storeu_pd(res + 2*i + 4, _mm256_permute2f128_pd(lo, hi, 0x31));
    }
#endif
    intl_ifft((const INTL_FFT_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios_Intl::execute_reverse_uint(double *res, const uint32_t *a) {
#ifdef USE_AVX512
    for (int32_t i = 0; i < Ns2; i += 8) {
        __m256i re_i32 = _mm256_loadu_si256((const __m256i*)(a + i));
        __m256i im_i32 = _mm256_loadu_si256((const __m256i*)(a + i + Ns2));
        __m512d re = _mm512_cvtepu32_pd(re_i32);
        __m512d im = _mm512_cvtepu32_pd(im_i32);
        _mm512_storeu_pd(res + 2*i,     _mm512_permutex2var_pd(re, idx_intl_lo, im));
        _mm512_storeu_pd(res + 2*i + 8, _mm512_permutex2var_pd(re, idx_intl_hi, im));
    }
#else
    for (int32_t i = 0; i < Ns2; i += 4) {
        __m128i re_i32 = _mm_loadu_si128((const __m128i*)(a + i));
        __m128i im_i32 = _mm_loadu_si128((const __m128i*)(a + i + Ns2));
        // AVX2 has no _mm256_cvtepu32_pd. Emulate: convert signed, then fix
        // negative results by adding 2^32 (for values where bit 31 was set).
        __m256d re = _mm256_cvtepi32_pd(re_i32);
        __m256d im = _mm256_cvtepi32_pd(im_i32);
        const __m256d fix = _mm256_set1_pd(4294967296.0);  // 2^32
        re = _mm256_add_pd(re, _mm256_and_pd(fix, _mm256_cmp_pd(re, _mm256_setzero_pd(), _CMP_LT_OQ)));
        im = _mm256_add_pd(im, _mm256_and_pd(fix, _mm256_cmp_pd(im, _mm256_setzero_pd(), _CMP_LT_OQ)));
        __m256d lo = _mm256_unpacklo_pd(re, im);
        __m256d hi = _mm256_unpackhi_pd(re, im);
        _mm256_storeu_pd(res + 2*i,     _mm256_permute2f128_pd(lo, hi, 0x20));
        _mm256_storeu_pd(res + 2*i + 4, _mm256_permute2f128_pd(lo, hi, 0x31));
    }
#endif
    intl_ifft((const INTL_FFT_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios_Intl::execute_reverse_torus64(double *res, const uint64_t *a) {
    const int64_t *aa = (const int64_t*)a;
    for (int32_t i = 0; i < Ns2; i++) {
        res[2*i]     = (double)aa[i];
        res[2*i + 1] = (double)aa[i + Ns2];
    }
    intl_ifft((const INTL_FFT_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios_Intl::execute_reverse_torus64_uint(double *res, const uint64_t *a) {
    for (int32_t i = 0; i < Ns2; i++) {
        res[2*i]     = (double)a[i];
        res[2*i + 1] = (double)a[i + Ns2];
    }
    intl_ifft((const INTL_FFT_PRECOMP*)tables_reverse, res);
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus32(uint32_t *res, const double *a) {
    auto *tables = (const INTL_FFT_PRECOMP*)tables_direct;
    // Scale is baked into forward twist — read input directly, no copy needed
    intl_fft_from(tables, a, real_inout_direct);
#ifdef USE_AVX512
    // Vectorized de-interleave + double→int32 extraction
    for (int32_t i = 0; i < Ns2; i += 8) {
        __m512d in0 = _mm512_load_pd(real_inout_direct + 2*i);
        __m512d in1 = _mm512_load_pd(real_inout_direct + 2*i + 8);
        __m512d re = _mm512_permutex2var_pd(in0, idx_deinl_re, in1);
        __m512d im = _mm512_permutex2var_pd(in0, idx_deinl_im, in1);
        __m256i re_i32 = _mm512_cvtepi64_epi32(_mm512_cvttpd_epi64(re));
        __m256i im_i32 = _mm512_cvtepi64_epi32(_mm512_cvttpd_epi64(im));
        _mm256_storeu_si256((__m256i*)(res + i), re_i32);
        _mm256_storeu_si256((__m256i*)(res + i + Ns2), im_i32);
    }
#else
    // AVX2: de-interleave [re0,im0,re1,im1,...] → split re[]+im[] with magic f64→i32
    {
        const __m256d magic_d = _mm256_set1_pd(6755399441055744.0);
        const __m256i magic_i = _mm256_set1_epi64x(0x4338000000000000LL);
        for (int32_t i = 0; i < Ns2; i += 4) {
            __m256d v0 = _mm256_loadu_pd(real_inout_direct + 2*i);       // [re0,im0,re1,im1]
            __m256d v1 = _mm256_loadu_pd(real_inout_direct + 2*i + 4);   // [re2,im2,re3,im3]
            // De-interleave: use shuffle_pd for within-lane, permute4x64 for cross-lane
            __m256d re_raw = _mm256_shuffle_pd(v0, v1, 0b0000); // [re0,re2,re1,re3]
            __m256d im_raw = _mm256_shuffle_pd(v0, v1, 0b1111); // [im0,im2,im1,im3]
            re_raw = _mm256_permute4x64_pd(re_raw, 0b11011000);  // [re0,re1,re2,re3]
            im_raw = _mm256_permute4x64_pd(im_raw, 0b11011000);  // [im0,im1,im2,im3]
            // Magic f64→i64→i32
            __m256i re_i64 = _mm256_sub_epi64(_mm256_castpd_si256(_mm256_add_pd(re_raw, magic_d)), magic_i);
            __m256i im_i64 = _mm256_sub_epi64(_mm256_castpd_si256(_mm256_add_pd(im_raw, magic_d)), magic_i);
            // Pack i64→i32: shuffle + permute
            __m256 combined = _mm256_shuffle_ps(_mm256_castsi256_ps(re_i64),
                                                _mm256_castsi256_ps(im_i64), _MM_SHUFFLE(2,0,2,0));
            __m256d ordered = _mm256_permute4x64_pd(_mm256_castps_pd(combined), _MM_SHUFFLE(3,1,2,0));
            __m128i re_i32 = _mm256_castsi256_si128(_mm256_castpd_si256(ordered));
            __m128i im_i32 = _mm256_extracti128_si256(_mm256_castpd_si256(ordered), 1);
            _mm_storeu_si128((__m128i*)(res + i), re_i32);
            _mm_storeu_si128((__m128i*)(res + i + Ns2), im_i32);
        }
    }
#endif
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus32_add(uint32_t *res, const double *a) {
    auto *tables = (const INTL_FFT_PRECOMP*)tables_direct;
    intl_fft_from(tables, a, real_inout_direct);
#ifdef USE_AVX512
    for (int32_t i = 0; i < Ns2; i += 8) {
        __m512d in0 = _mm512_load_pd(real_inout_direct + 2*i);
        __m512d in1 = _mm512_load_pd(real_inout_direct + 2*i + 8);
        __m512d re = _mm512_permutex2var_pd(in0, idx_deinl_re, in1);
        __m512d im = _mm512_permutex2var_pd(in0, idx_deinl_im, in1);
        __m256i re_i32 = _mm512_cvtepi64_epi32(_mm512_cvttpd_epi64(re));
        __m256i im_i32 = _mm512_cvtepi64_epi32(_mm512_cvttpd_epi64(im));
        __m256i old_re = _mm256_loadu_si256((__m256i*)(res + i));
        __m256i old_im = _mm256_loadu_si256((__m256i*)(res + i + Ns2));
        _mm256_storeu_si256((__m256i*)(res + i), _mm256_add_epi32(old_re, re_i32));
        _mm256_storeu_si256((__m256i*)(res + i + Ns2), _mm256_add_epi32(old_im, im_i32));
    }
#else
    {
        const __m256d magic_d = _mm256_set1_pd(6755399441055744.0);
        const __m256i magic_i = _mm256_set1_epi64x(0x4338000000000000LL);
        for (int32_t i = 0; i < Ns2; i += 4) {
            __m256d v0 = _mm256_loadu_pd(real_inout_direct + 2*i);
            __m256d v1 = _mm256_loadu_pd(real_inout_direct + 2*i + 4);
            __m256d re_raw = _mm256_permute4x64_pd(_mm256_shuffle_pd(v0, v1, 0b0000), 0b11011000);
            __m256d im_raw = _mm256_permute4x64_pd(_mm256_shuffle_pd(v0, v1, 0b1111), 0b11011000);
            __m256i re_i64 = _mm256_sub_epi64(_mm256_castpd_si256(_mm256_add_pd(re_raw, magic_d)), magic_i);
            __m256i im_i64 = _mm256_sub_epi64(_mm256_castpd_si256(_mm256_add_pd(im_raw, magic_d)), magic_i);
            __m256 combined = _mm256_shuffle_ps(_mm256_castsi256_ps(re_i64),
                                                _mm256_castsi256_ps(im_i64), _MM_SHUFFLE(2,0,2,0));
            __m256d ordered = _mm256_permute4x64_pd(_mm256_castps_pd(combined), _MM_SHUFFLE(3,1,2,0));
            __m128i re_i32 = _mm256_castsi256_si128(_mm256_castpd_si256(ordered));
            __m128i im_i32 = _mm256_extracti128_si256(_mm256_castpd_si256(ordered), 1);
            __m128i old_re = _mm_loadu_si128((__m128i*)(res + i));
            __m128i old_im = _mm_loadu_si128((__m128i*)(res + i + Ns2));
            _mm_storeu_si128((__m128i*)(res + i), _mm_add_epi32(old_re, re_i32));
            _mm_storeu_si128((__m128i*)(res + i + Ns2), _mm_add_epi32(old_im, im_i32));
        }
    }
#endif
}

// Full-range f64→u64 via IEEE754 bit extraction (from non-interleaved SPQLIOS).
// Uses offset 1075 (= bias + mantissa_bits) so that for values > 2^52, the left
// shift naturally wraps modulo 2^64, matching torus arithmetic. This is critical
// because polynomial multiplication results can exceed 2^63 before modular reduction.
static inline __m256i f64_to_i64(__m256d x) {
    const __m256i bits = _mm256_castpd_si256(x);
    const __m256i mantissa = _mm256_or_si256(
        _mm256_and_si256(bits, _mm256_set1_epi64x(0xFFFFFFFFFFFFF)),
        _mm256_set1_epi64x(0x10000000000000));
    const __m256i biased_exp = _mm256_and_si256(
        _mm256_srli_epi64(bits, 52), _mm256_set1_epi64x(0x7FF));
    const __m256i sign_mask = _mm256_sub_epi64(
        _mm256_setzero_si256(), _mm256_srli_epi64(bits, 63));
    const __m256i offset = _mm256_set1_epi64x(1075);
    const __m256i shift = _mm256_sub_epi64(biased_exp, offset);
    const __m256i neg_shift = _mm256_sub_epi64(offset, biased_exp);
    const __m256i val = _mm256_or_si256(
        _mm256_sllv_epi64(mantissa, shift),
        _mm256_srlv_epi64(mantissa, neg_shift));
    return _mm256_sub_epi64(_mm256_xor_si256(val, sign_mask), sign_mask);
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus64(uint64_t *res, double *a) {
    auto *tables = (const INTL_FFT_PRECOMP*)tables_direct;
    intl_fft_from(tables, a, real_inout_direct);
    // De-interleave + convert f64→u64 with full-range IEEE754 bit extraction.
    for (int32_t i = 0; i < Ns2; i += 4) {
        __m256d v0 = _mm256_loadu_pd(real_inout_direct + 2*i);
        __m256d v1 = _mm256_loadu_pd(real_inout_direct + 2*i + 4);
        __m256d re = _mm256_permute4x64_pd(
            _mm256_shuffle_pd(v0, v1, 0b0000), 0b11011000);
        __m256d im = _mm256_permute4x64_pd(
            _mm256_shuffle_pd(v0, v1, 0b1111), 0b11011000);
        _mm256_storeu_si256((__m256i*)(res + i), f64_to_i64(re));
        _mm256_storeu_si256((__m256i*)(res + i + Ns2), f64_to_i64(im));
    }
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus64_add(uint64_t *res, double *a) {
    auto *tables = (const INTL_FFT_PRECOMP*)tables_direct;
    intl_fft_from(tables, a, real_inout_direct);
    for (int32_t i = 0; i < Ns2; i += 4) {
        __m256d v0 = _mm256_loadu_pd(real_inout_direct + 2*i);
        __m256d v1 = _mm256_loadu_pd(real_inout_direct + 2*i + 4);
        __m256d re = _mm256_permute4x64_pd(
            _mm256_shuffle_pd(v0, v1, 0b0000), 0b11011000);
        __m256d im = _mm256_permute4x64_pd(
            _mm256_shuffle_pd(v0, v1, 0b1111), 0b11011000);
        __m256i old_re = _mm256_loadu_si256((const __m256i*)(res + i));
        __m256i old_im = _mm256_loadu_si256((const __m256i*)(res + i + Ns2));
        _mm256_storeu_si256((__m256i*)(res + i),
                            _mm256_add_epi64(old_re, f64_to_i64(re)));
        _mm256_storeu_si256((__m256i*)(res + i + Ns2),
                            _mm256_add_epi64(old_im, f64_to_i64(im)));
    }
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus32_q(uint32_t *res, const double *a, const uint32_t q) {
    auto *tables = (const INTL_FFT_PRECOMP*)tables_direct;
    intl_fft_from(tables, a, real_inout_direct);
    for (int32_t i = 0; i < Ns2; i++) {
        res[i] = uint32_t((int64_t(real_inout_direct[2*i])%q+q)%q);
        res[i+Ns2] = uint32_t((int64_t(real_inout_direct[2*i+1])%q+q)%q);
    }
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus32_rescale(uint32_t *res, const double *a, const double D) {
    auto *tables = (const INTL_FFT_PRECOMP*)tables_direct;
    intl_fft_from(tables, a, real_inout_direct);
    for (int32_t i = 0; i < Ns2; i++) {
        res[i] = (uint32_t)(int64_t)(real_inout_direct[2*i]/D);
        res[i+Ns2] = (uint32_t)(int64_t)(real_inout_direct[2*i+1]/D);
    }
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus32_rescale_clpx(
    uint32_t *res, const double *a, const double q, const uint32_t plain_modulus) {
    execute_direct_torus32(res, a);
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus64_rescale(uint64_t *res, const double *a, const double D) {
    auto *tables = (const INTL_FFT_PRECOMP*)tables_direct;
    alignas(64) double tmp[N];
    intl_fft_from(tables, a, tmp);
    for (int32_t i = 0; i < Ns2; i++) {
        // Cast through int64_t to avoid UB on negative doubles → uint64_t.
        // std::llround gives the correctly rounded signed integer.
        res[i] = static_cast<uint64_t>(static_cast<int64_t>(std::llround(tmp[2*i]/D)));
        res[i+Ns2] = static_cast<uint64_t>(static_cast<int64_t>(std::llround(tmp[2*i+1]/D)));
    }
}

void FFT_Processor_Spqlios_Intl::execute_direct_torus64_rescale_clpx(
    uint64_t *res, const double *a, const uint32_t plain_modulus) {
    execute_direct_torus64(res, const_cast<double*>(a));
}

thread_local FFT_Processor_Spqlios_Intl fftplvl1(TFHEpp::lvl1param::n);
thread_local FFT_Processor_Spqlios_Intl fftplvl2(TFHEpp::lvl2param::n);
thread_local FFT_Processor_Spqlios_Intl fftplvl3(TFHEpp::lvl3param::n);
