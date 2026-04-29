// Stockham radix-4 FFT with AVX2 FMA — split real/imaginary layout.
//
// Drop-in replacement for spqlios-fft-fma.s / spqlios-ifft-fma.s:
//   void fft(const void *tables, double *data);
//   void ifft(const void *tables, double *data);
//   void *new_fft_table(int32_t nn);
//   void *new_ifft_table(int32_t nn);
//   double *fft_table_get_buffer(const void *);
//   double *ifft_table_get_buffer(const void *);
//
// Data layout: re[0..ns4-1], im[0..ns4-1]  where ns4 = nn/2 = n/4
// Twiddle layout (per stage): [cos0..3|sin0..3|cos4..7|sin4..7|...] × 3 sets
//
// Key optimizations over the assembly radix-2 version:
//   - Radix-4: 4 passes instead of 8+1 for ns4=256 (N=512)
//   - Forward twist bakes in 2/N scaling (eliminates separate scaling loop)
//   - Fused last butterfly + twist (saves one full data sweep)
//   - Fused first butterfly + inverse twist (saves one full data sweep)

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <immintrin.h>

#include <params.hpp>
#include "fft_processor_spqlios.h"

// ── Trig helpers ────────────────────────────────────────────────────────────
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

// ── Complex multiply helper (split format, AVX2) ───────────────────────────
// res_re = a_re*w_re - a_im*w_im
// res_im = a_im*w_re + a_re*w_im
static inline void cmul_split_avx2(
    __m256d a_re, __m256d a_im,
    __m256d w_re, __m256d w_im,
    __m256d &r_re, __m256d &r_im)
{
    r_re = _mm256_mul_pd(a_re, w_re);
    r_re = _mm256_fnmadd_pd(a_im, w_im, r_re);
    r_im = _mm256_mul_pd(a_im, w_re);
    r_im = _mm256_fmadd_pd(a_re, w_im, r_im);
}

// ── Vectorized q=1 final pass for inverse FFT ──────────────────────────────
// Processes 4 j-values at once: loads from stride-s, inverse butterfly + twiddle,
// 4×4 transpose, stores contiguous.
static void inv_final_pass_vec(
    int32_t s, const double *tw,
    const double *src_re, const double *src_im,
    double *dst_re, double *dst_im)
{
    const __m256d neg = _mm256_set1_pd(-0.0);
    for (int32_t j = 0; j < s; j += 4) {
        __m256d a_re = _mm256_loadu_pd(src_re + j);
        __m256d a_im = _mm256_loadu_pd(src_im + j);
        __m256d b_re = _mm256_loadu_pd(src_re + j + s);
        __m256d b_im = _mm256_loadu_pd(src_im + j + s);
        __m256d c_re = _mm256_loadu_pd(src_re + j + 2*s);
        __m256d c_im = _mm256_loadu_pd(src_im + j + 2*s);
        __m256d d_re = _mm256_loadu_pd(src_re + j + 3*s);
        __m256d d_im = _mm256_loadu_pd(src_im + j + 3*s);

        __m256d apc_re = _mm256_add_pd(a_re, c_re);
        __m256d apc_im = _mm256_add_pd(a_im, c_im);
        __m256d amc_re = _mm256_sub_pd(a_re, c_re);
        __m256d amc_im = _mm256_sub_pd(a_im, c_im);
        __m256d bpd_re = _mm256_add_pd(b_re, d_re);
        __m256d bpd_im = _mm256_add_pd(b_im, d_im);
        __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
        __m256d bmd_im = _mm256_sub_pd(b_im, d_im);
        // Inverse j: [bmd_im, -bmd_re]
        __m256d jbmd_re = bmd_im;
        __m256d jbmd_im = _mm256_xor_pd(bmd_re, neg);

        __m256d o0_re = _mm256_add_pd(apc_re, bpd_re);
        __m256d o0_im = _mm256_add_pd(apc_im, bpd_im);

        __m256d t1_re = _mm256_sub_pd(amc_re, jbmd_re);
        __m256d t1_im = _mm256_sub_pd(amc_im, jbmd_im);
        __m256d t2_re = _mm256_sub_pd(apc_re, bpd_re);
        __m256d t2_im = _mm256_sub_pd(apc_im, bpd_im);
        __m256d t3_re = _mm256_add_pd(amc_re, jbmd_re);
        __m256d t3_im = _mm256_add_pd(amc_im, jbmd_im);

        // Gather per-j twiddles
        __m256d w1r = _mm256_set_pd(tw[(j+3)*6+0], tw[(j+2)*6+0], tw[(j+1)*6+0], (j==0)?1.0:tw[j*6+0]);
        __m256d w1i = _mm256_set_pd(tw[(j+3)*6+1], tw[(j+2)*6+1], tw[(j+1)*6+1], (j==0)?0.0:tw[j*6+1]);
        __m256d o1_re, o1_im;
        cmul_split_avx2(t1_re, t1_im, w1r, w1i, o1_re, o1_im);

        __m256d w2r = _mm256_set_pd(tw[(j+3)*6+2], tw[(j+2)*6+2], tw[(j+1)*6+2], (j==0)?1.0:tw[j*6+2]);
        __m256d w2i = _mm256_set_pd(tw[(j+3)*6+3], tw[(j+2)*6+3], tw[(j+1)*6+3], (j==0)?0.0:tw[j*6+3]);
        __m256d o2_re, o2_im;
        cmul_split_avx2(t2_re, t2_im, w2r, w2i, o2_re, o2_im);

        __m256d w3r = _mm256_set_pd(tw[(j+3)*6+4], tw[(j+2)*6+4], tw[(j+1)*6+4], (j==0)?1.0:tw[j*6+4]);
        __m256d w3i = _mm256_set_pd(tw[(j+3)*6+5], tw[(j+2)*6+5], tw[(j+1)*6+5], (j==0)?0.0:tw[j*6+5]);
        __m256d o3_re, o3_im;
        cmul_split_avx2(t3_re, t3_im, w3r, w3i, o3_re, o3_im);

        // 4×4 transpose
        __m256d t01_lo_re = _mm256_unpacklo_pd(o0_re, o1_re);
        __m256d t01_hi_re = _mm256_unpackhi_pd(o0_re, o1_re);
        __m256d t23_lo_re = _mm256_unpacklo_pd(o2_re, o3_re);
        __m256d t23_hi_re = _mm256_unpackhi_pd(o2_re, o3_re);
        __m256d row0_re = _mm256_permute2f128_pd(t01_lo_re, t23_lo_re, 0x20);
        __m256d row1_re = _mm256_permute2f128_pd(t01_hi_re, t23_hi_re, 0x20);
        __m256d row2_re = _mm256_permute2f128_pd(t01_lo_re, t23_lo_re, 0x31);
        __m256d row3_re = _mm256_permute2f128_pd(t01_hi_re, t23_hi_re, 0x31);

        __m256d t01_lo_im = _mm256_unpacklo_pd(o0_im, o1_im);
        __m256d t01_hi_im = _mm256_unpackhi_pd(o0_im, o1_im);
        __m256d t23_lo_im = _mm256_unpacklo_pd(o2_im, o3_im);
        __m256d t23_hi_im = _mm256_unpackhi_pd(o2_im, o3_im);
        __m256d row0_im = _mm256_permute2f128_pd(t01_lo_im, t23_lo_im, 0x20);
        __m256d row1_im = _mm256_permute2f128_pd(t01_hi_im, t23_hi_im, 0x20);
        __m256d row2_im = _mm256_permute2f128_pd(t01_lo_im, t23_lo_im, 0x31);
        __m256d row3_im = _mm256_permute2f128_pd(t01_hi_im, t23_hi_im, 0x31);

        int32_t ob = 4 * j;
        _mm256_storeu_pd(dst_re + ob, row0_re);      _mm256_storeu_pd(dst_im + ob, row0_im);
        _mm256_storeu_pd(dst_re + ob + 4, row1_re);  _mm256_storeu_pd(dst_im + ob + 4, row1_im);
        _mm256_storeu_pd(dst_re + ob + 8, row2_re);  _mm256_storeu_pd(dst_im + ob + 8, row2_im);
        _mm256_storeu_pd(dst_re + ob + 12, row3_re); _mm256_storeu_pd(dst_im + ob + 12, row3_im);
    }
}

// ── Table structure ─────────────────────────────────────────────────────────
// For Stockham radix-4 with ns4 complex values:
//   - nstages = log4(ns4) butterfly stages
//   - Per stage (except last): 3*(ns4/4) twiddle entries (w1,w2,w3 per group)
//     stored as blocks of 8 doubles: [cos0..3|sin0..3] for each w
//   - Forward twist: ns4 entries with 2/N baked in
//   - Inverse twist: ns4 entries (no scaling)
struct STOCKHAM_R4_PRECOMP {
    int32_t n;        // = 2*nn
    int32_t ns4;      // = nn/2 = number of complex values
    double *trig_fwd; // butterfly twiddles + fused twist for forward
    double *trig_inv; // fused twist + butterfly twiddles for inverse
    double *data;     // working buffer
    double *scratch;  // scratch buffer for out-of-place
    void *buf;        // raw allocation
};

// ── Stockham radix-4 core (forward DIF) ─────────────────────────────────────
// Operates on ns4 complex values in split re/im layout.
// src_re/im → dst_re/im (out of place, caller ping-pongs buffers)
// Restructured butterfly: compute a+c, a-c, b+d, b-d first, then produce each
// output and store immediately. This reduces live values from ~30 to ~14, fitting
// in 16 YMM registers without spills. The j=0 case (no twiddle) is split into a
// separate loop to avoid branching in the hot path.
static void stockham_r4_fwd_pass(
    int32_t ns4, int32_t q, int32_t s,
    const double *tw,
    const double *src_re, const double *src_im,
    double *dst_re, double *dst_im)
{
    const __m256d neg = _mm256_set1_pd(-0.0);
    const int32_t stride = q * s;

    // j=0: no twiddle multiply
    for (int32_t p = 0; p < q; p += 4) {
        const int32_t idx = s * p;
        __m256d a_re = _mm256_loadu_pd(src_re + idx);
        __m256d a_im = _mm256_loadu_pd(src_im + idx);
        __m256d c_re = _mm256_loadu_pd(src_re + idx + 2*stride);
        __m256d c_im = _mm256_loadu_pd(src_im + idx + 2*stride);
        __m256d apc_re = _mm256_add_pd(a_re, c_re);
        __m256d apc_im = _mm256_add_pd(a_im, c_im);
        __m256d amc_re = _mm256_sub_pd(a_re, c_re);
        __m256d amc_im = _mm256_sub_pd(a_im, c_im);

        __m256d b_re = _mm256_loadu_pd(src_re + idx + stride);
        __m256d b_im = _mm256_loadu_pd(src_im + idx + stride);
        __m256d d_re = _mm256_loadu_pd(src_re + idx + 3*stride);
        __m256d d_im = _mm256_loadu_pd(src_im + idx + 3*stride);
        __m256d bpd_re = _mm256_add_pd(b_re, d_re);
        __m256d bpd_im = _mm256_add_pd(b_im, d_im);
        __m256d bmd_im = _mm256_sub_pd(b_im, d_im);
        __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
        __m256d jbmd_re = _mm256_xor_pd(bmd_im, neg);

        _mm256_storeu_pd(dst_re + p, _mm256_add_pd(apc_re, bpd_re));
        _mm256_storeu_pd(dst_im + p, _mm256_add_pd(apc_im, bpd_im));
        _mm256_storeu_pd(dst_re + q + p, _mm256_sub_pd(amc_re, jbmd_re));
        _mm256_storeu_pd(dst_im + q + p, _mm256_sub_pd(amc_im, bmd_re));
        _mm256_storeu_pd(dst_re + 2*q + p, _mm256_sub_pd(apc_re, bpd_re));
        _mm256_storeu_pd(dst_im + 2*q + p, _mm256_sub_pd(apc_im, bpd_im));
        _mm256_storeu_pd(dst_re + 3*q + p, _mm256_add_pd(amc_re, jbmd_re));
        _mm256_storeu_pd(dst_im + 3*q + p, _mm256_add_pd(amc_im, bmd_re));
    }

    // j=1..s-1: with twiddle multiply.
    // Hand-written inline asm for zero YMM spills. Key trick: compute out0 and out2
    // first (uses apc,bpd), freeing 4 regs before computing out1 and out3 (uses amc,jbmd).
    // ymm0=w1_re, ymm1=w1_im, ymm2=w2_re, ymm3=w2_im, ymm4=w3_re, ymm5=w3_im, ymm6=neg
    for (int32_t j = 1; j < s; j++) {
        const double *twj = tw + j * 6;
        const int64_t stride_bytes = (int64_t)stride * 8;
        const int64_t q_bytes = (int64_t)q * 8;
        const int64_t out_off_bytes = (int64_t)(q * 4 * j) * 8;

        __asm__ __volatile__ (
            "vbroadcastsd   (%[tw]),    %%ymm0\n\t"
            "vbroadcastsd  8(%[tw]),    %%ymm1\n\t"
            "vbroadcastsd 16(%[tw]),    %%ymm2\n\t"
            "vbroadcastsd 24(%[tw]),    %%ymm3\n\t"
            "vbroadcastsd 32(%[tw]),    %%ymm4\n\t"
            "vbroadcastsd 40(%[tw]),    %%ymm5\n\t"
            : : [tw] "r"(twj)
            : "ymm0","ymm1","ymm2","ymm3","ymm4","ymm5"
        );

        for (int32_t p = 0; p < q; p += 4) {
            // Compute pointers for this iteration
            const double *sre = src_re + j + (int64_t)s * p;
            const double *sim = src_im + j + (int64_t)s * p;
            double *dre = dst_re + q * 4 * j + p;
            double *dim = dst_im + q * 4 * j + p;
            // All destination offsets relative to dre/dim in bytes
            __asm__ __volatile__ (
                // --- Load a,c → apc, amc ---
                "vmovupd    (%[sre]),         %%ymm15\n\t"
                "vmovupd    (%[sre],%[st2]),  %%ymm7\n\t"
                "vaddpd     %%ymm7,  %%ymm15, %%ymm8\n\t"
                "vsubpd     %%ymm7,  %%ymm15, %%ymm9\n\t"
                "vmovupd    (%[sim]),         %%ymm15\n\t"
                "vmovupd    (%[sim],%[st2]),  %%ymm7\n\t"
                "vaddpd     %%ymm7,  %%ymm15, %%ymm10\n\t"
                "vsubpd     %%ymm7,  %%ymm15, %%ymm14\n\t"
                // --- Load b,d → bpd, bmd ---
                "vmovupd    (%[sre],%[st1]),  %%ymm15\n\t"
                "vmovupd    (%[sre],%[st3]),  %%ymm7\n\t"
                "vaddpd     %%ymm7,  %%ymm15, %%ymm11\n\t"
                "vsubpd     %%ymm7,  %%ymm15, %%ymm13\n\t"
                "vmovupd    (%[sim],%[st1]),  %%ymm15\n\t"
                "vmovupd    (%[sim],%[st3]),  %%ymm7\n\t"
                "vaddpd     %%ymm7,  %%ymm15, %%ymm12\n\t"
                "vsubpd     %%ymm7,  %%ymm15, %%ymm7\n\t"
                // --- out0 = apc + bpd ---
                "vaddpd     %%ymm11, %%ymm8,  %%ymm15\n\t"
                "vmovupd    %%ymm15, (%[dre])\n\t"
                "vaddpd     %%ymm12, %%ymm10, %%ymm15\n\t"
                "vmovupd    %%ymm15, (%[dim])\n\t"
                // --- out2 = (apc - bpd) * w2 ---
                "vsubpd     %%ymm11, %%ymm8,  %%ymm8\n\t"
                "vsubpd     %%ymm12, %%ymm10, %%ymm10\n\t"
                "vmulpd     %%ymm2,  %%ymm8,  %%ymm11\n\t"
                "vfnmadd231pd %%ymm3,%%ymm10, %%ymm11\n\t"
                "vmovupd    %%ymm11, (%[dre],%[q2])\n\t"
                "vmulpd     %%ymm2,  %%ymm10, %%ymm12\n\t"
                "vfmadd231pd %%ymm3, %%ymm8,  %%ymm12\n\t"
                "vmovupd    %%ymm12, (%[dim],%[q2])\n\t"
                // --- jbmd ---
                "vxorpd     %%ymm6,  %%ymm7,  %%ymm8\n\t"
                // --- out1 = (amc - jbmd) * w1 ---
                "vsubpd     %%ymm8,  %%ymm9,  %%ymm10\n\t"
                "vsubpd     %%ymm13, %%ymm14, %%ymm11\n\t"
                "vmulpd     %%ymm0,  %%ymm10, %%ymm12\n\t"
                "vfnmadd231pd %%ymm1,%%ymm11, %%ymm12\n\t"
                "vmovupd    %%ymm12, (%[dre],%[q1])\n\t"
                "vmulpd     %%ymm0,  %%ymm11, %%ymm15\n\t"
                "vfmadd231pd %%ymm1, %%ymm10, %%ymm15\n\t"
                "vmovupd    %%ymm15, (%[dim],%[q1])\n\t"
                // --- out3 = (amc + jbmd) * w3 ---
                "vaddpd     %%ymm8,  %%ymm9,  %%ymm10\n\t"
                "vaddpd     %%ymm13, %%ymm14, %%ymm11\n\t"
                "vmulpd     %%ymm4,  %%ymm10, %%ymm12\n\t"
                "vfnmadd231pd %%ymm5,%%ymm11, %%ymm12\n\t"
                "vmovupd    %%ymm12, (%[dre],%[q3])\n\t"
                "vmulpd     %%ymm4,  %%ymm11, %%ymm15\n\t"
                "vfmadd231pd %%ymm5, %%ymm10, %%ymm15\n\t"
                "vmovupd    %%ymm15, (%[dim],%[q3])\n\t"
                : : [sre] "r"(sre), [sim] "r"(sim),
                    [st1] "r"(stride_bytes), [st2] "r"(stride_bytes*2), [st3] "r"(stride_bytes*3),
                    [dre] "r"(dre), [dim] "r"(dim),
                    [q1] "r"(q_bytes), [q2] "r"(2*q_bytes), [q3] "r"(3*q_bytes)
                : "ymm7","ymm8","ymm9","ymm10","ymm11","ymm12","ymm13","ymm14","ymm15","memory"
            );
        }
    }
}

// Inverse pass with same register-pressure optimizations as forward
static void stockham_r4_inv_pass(
    int32_t ns4, int32_t q, int32_t s,
    const double *tw,
    const double *src_re, const double *src_im,
    double *dst_re, double *dst_im)
{
    const __m256d neg = _mm256_set1_pd(-0.0);
    const int32_t stride = q * s;

    // j=0: no twiddle
    for (int32_t p = 0; p < q; p += 4) {
        const int32_t idx = s * p;
        __m256d a_re = _mm256_loadu_pd(src_re + idx);
        __m256d a_im = _mm256_loadu_pd(src_im + idx);
        __m256d c_re = _mm256_loadu_pd(src_re + idx + 2*stride);
        __m256d c_im = _mm256_loadu_pd(src_im + idx + 2*stride);
        __m256d apc_re = _mm256_add_pd(a_re, c_re);
        __m256d apc_im = _mm256_add_pd(a_im, c_im);
        __m256d amc_re = _mm256_sub_pd(a_re, c_re);
        __m256d amc_im = _mm256_sub_pd(a_im, c_im);
        __m256d b_re = _mm256_loadu_pd(src_re + idx + stride);
        __m256d b_im = _mm256_loadu_pd(src_im + idx + stride);
        __m256d d_re = _mm256_loadu_pd(src_re + idx + 3*stride);
        __m256d d_im = _mm256_loadu_pd(src_im + idx + 3*stride);
        __m256d bpd_re = _mm256_add_pd(b_re, d_re);
        __m256d bpd_im = _mm256_add_pd(b_im, d_im);
        __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
        __m256d bmd_im = _mm256_sub_pd(b_im, d_im);
        // Inverse j: jbmd = [bmd_im, -bmd_re]
        __m256d jbmd_re = bmd_im;
        __m256d neg_bmd_re = _mm256_xor_pd(bmd_re, neg);

        _mm256_storeu_pd(dst_re + p, _mm256_add_pd(apc_re, bpd_re));
        _mm256_storeu_pd(dst_im + p, _mm256_add_pd(apc_im, bpd_im));
        _mm256_storeu_pd(dst_re + q + p, _mm256_sub_pd(amc_re, jbmd_re));
        _mm256_storeu_pd(dst_im + q + p, _mm256_sub_pd(amc_im, neg_bmd_re));
        _mm256_storeu_pd(dst_re + 2*q + p, _mm256_sub_pd(apc_re, bpd_re));
        _mm256_storeu_pd(dst_im + 2*q + p, _mm256_sub_pd(apc_im, bpd_im));
        _mm256_storeu_pd(dst_re + 3*q + p, _mm256_add_pd(amc_re, jbmd_re));
        _mm256_storeu_pd(dst_im + 3*q + p, _mm256_add_pd(amc_im, neg_bmd_re));
    }

    for (int32_t j = 1; j < s; j++) {
        const double *twj = tw + j * 6;
        __m256d w1_re = _mm256_broadcast_sd(twj);
        __m256d w1_im = _mm256_broadcast_sd(twj + 1);
        __m256d w2_re = _mm256_broadcast_sd(twj + 2);
        __m256d w2_im = _mm256_broadcast_sd(twj + 3);
        __m256d w3_re = _mm256_broadcast_sd(twj + 4);
        __m256d w3_im = _mm256_broadcast_sd(twj + 5);
        const int32_t out_off = q * 4 * j;

        for (int32_t p = 0; p < q; p += 4) {
            const int32_t idx = j + s * p;
            __m256d a_re = _mm256_loadu_pd(src_re + idx);
            __m256d a_im = _mm256_loadu_pd(src_im + idx);
            __m256d c_re = _mm256_loadu_pd(src_re + idx + 2*stride);
            __m256d c_im = _mm256_loadu_pd(src_im + idx + 2*stride);
            __m256d apc_re = _mm256_add_pd(a_re, c_re);
            __m256d apc_im = _mm256_add_pd(a_im, c_im);
            __m256d amc_re = _mm256_sub_pd(a_re, c_re);
            __m256d amc_im = _mm256_sub_pd(a_im, c_im);
            __m256d b_re = _mm256_loadu_pd(src_re + idx + stride);
            __m256d b_im = _mm256_loadu_pd(src_im + idx + stride);
            __m256d d_re = _mm256_loadu_pd(src_re + idx + 3*stride);
            __m256d d_im = _mm256_loadu_pd(src_im + idx + 3*stride);
            __m256d bpd_re = _mm256_add_pd(b_re, d_re);
            __m256d bpd_im = _mm256_add_pd(b_im, d_im);
            __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
            __m256d bmd_im = _mm256_sub_pd(b_im, d_im);
            __m256d jbmd_re = bmd_im;
            __m256d neg_bmd_re = _mm256_xor_pd(bmd_re, neg);

            _mm256_storeu_pd(dst_re + out_off + p, _mm256_add_pd(apc_re, bpd_re));
            _mm256_storeu_pd(dst_im + out_off + p, _mm256_add_pd(apc_im, bpd_im));

            __m256d t_re = _mm256_sub_pd(amc_re, jbmd_re);
            __m256d t_im = _mm256_sub_pd(amc_im, neg_bmd_re);
            __m256d r_re = _mm256_mul_pd(t_re, w1_re);
            r_re = _mm256_fnmadd_pd(t_im, w1_im, r_re);
            __m256d r_im = _mm256_mul_pd(t_im, w1_re);
            r_im = _mm256_fmadd_pd(t_re, w1_im, r_im);
            _mm256_storeu_pd(dst_re + out_off + q + p, r_re);
            _mm256_storeu_pd(dst_im + out_off + q + p, r_im);

            t_re = _mm256_sub_pd(apc_re, bpd_re);
            t_im = _mm256_sub_pd(apc_im, bpd_im);
            r_re = _mm256_mul_pd(t_re, w2_re);
            r_re = _mm256_fnmadd_pd(t_im, w2_im, r_re);
            r_im = _mm256_mul_pd(t_im, w2_re);
            r_im = _mm256_fmadd_pd(t_re, w2_im, r_im);
            _mm256_storeu_pd(dst_re + out_off + 2*q + p, r_re);
            _mm256_storeu_pd(dst_im + out_off + 2*q + p, r_im);

            t_re = _mm256_add_pd(amc_re, jbmd_re);
            t_im = _mm256_add_pd(amc_im, neg_bmd_re);
            r_re = _mm256_mul_pd(t_re, w3_re);
            r_re = _mm256_fnmadd_pd(t_im, w3_im, r_re);
            r_im = _mm256_mul_pd(t_im, w3_re);
            r_im = _mm256_fmadd_pd(t_re, w3_im, r_im);
            _mm256_storeu_pd(dst_re + out_off + 3*q + p, r_re);
            _mm256_storeu_pd(dst_im + out_off + 3*q + p, r_im);
        }
    }
}

// ── Full forward FFT (multi-pass + fused twist) ─────────────────────────────
static void stockham_fft_forward(
    int32_t ns4,
    const double *butterfly_tw,  // butterfly twiddles for all passes
    const double *twist_tw,      // twist twiddles (with 2/N baked in), split [cos|sin]
    double *re, double *im,
    double *sre, double *sim)
{
    double *cur_re = re, *cur_im = im;
    double *dst_re = sre, *dst_im = sim;
    const double *tw = butterfly_tw;

    int32_t q = ns4 / 4;
    int32_t s = 1;

    // All but the last pass: standard Stockham radix-4
    while (q >= 4) {
        stockham_r4_fwd_pass(ns4, q, s, tw, cur_re, cur_im, dst_re, dst_im);
        tw += s * 6;
        // Swap buffers
        double *t;
        t = cur_re; cur_re = dst_re; dst_re = t;
        t = cur_im; cur_im = dst_im; dst_im = t;
        s *= 4; q /= 4;
    }

    // Last pass: q=1, fused with twist.
    // Process 4 j-values at once: load from stride-s positions (contiguous per group),
    // do 4 butterflies in parallel, 4×4 transpose, apply twist, store contiguous.
    {
        const __m256d neg = _mm256_set1_pd(-0.0);
        const double *tw_cos = twist_tw;
        const double *tw_sin = twist_tw + ns4;

        for (int32_t j = 0; j < s; j += 4) {
            // Load 4 consecutive values from each quarter: a[j..j+3], b[j+s..j+s+3], etc.
            __m256d a_re = _mm256_loadu_pd(cur_re + j);
            __m256d a_im = _mm256_loadu_pd(cur_im + j);
            __m256d b_re = _mm256_loadu_pd(cur_re + j + s);
            __m256d b_im = _mm256_loadu_pd(cur_im + j + s);
            __m256d c_re = _mm256_loadu_pd(cur_re + j + 2*s);
            __m256d c_im = _mm256_loadu_pd(cur_im + j + 2*s);
            __m256d d_re = _mm256_loadu_pd(cur_re + j + 3*s);
            __m256d d_im = _mm256_loadu_pd(cur_im + j + 3*s);

            // Radix-4 butterfly (4 butterflies in parallel)
            __m256d apc_re = _mm256_add_pd(a_re, c_re);
            __m256d apc_im = _mm256_add_pd(a_im, c_im);
            __m256d amc_re = _mm256_sub_pd(a_re, c_re);
            __m256d amc_im = _mm256_sub_pd(a_im, c_im);
            __m256d bpd_re = _mm256_add_pd(b_re, d_re);
            __m256d bpd_im = _mm256_add_pd(b_im, d_im);
            __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
            __m256d bmd_im = _mm256_sub_pd(b_im, d_im);
            __m256d jbmd_re = _mm256_xor_pd(bmd_im, neg); // -bmd_im
            __m256d jbmd_im = bmd_re;

            // out0 = (apc + bpd) — no twiddle
            __m256d o0_re = _mm256_add_pd(apc_re, bpd_re);
            __m256d o0_im = _mm256_add_pd(apc_im, bpd_im);

            // out1 = (amc - jbmd) * w1[j] — each lane has different twiddle
            __m256d t1_re = _mm256_sub_pd(amc_re, jbmd_re);
            __m256d t1_im = _mm256_sub_pd(amc_im, jbmd_im);
            // Load per-j twiddles: tw[j*6+0..3] for w1_re, tw[j*6+1..3] for w1_im
            // But twiddles are stored per-j (6 doubles each), need to gather
            __m256d w1r = _mm256_set_pd(tw[(j+3)*6+0], tw[(j+2)*6+0], tw[(j+1)*6+0], (j==0)?1.0:tw[j*6+0]);
            __m256d w1i = _mm256_set_pd(tw[(j+3)*6+1], tw[(j+2)*6+1], tw[(j+1)*6+1], (j==0)?0.0:tw[j*6+1]);
            __m256d o1_re, o1_im;
            cmul_split_avx2(t1_re, t1_im, w1r, w1i, o1_re, o1_im);

            // out2 = (apc - bpd) * w2[j]
            __m256d t2_re = _mm256_sub_pd(apc_re, bpd_re);
            __m256d t2_im = _mm256_sub_pd(apc_im, bpd_im);
            __m256d w2r = _mm256_set_pd(tw[(j+3)*6+2], tw[(j+2)*6+2], tw[(j+1)*6+2], (j==0)?1.0:tw[j*6+2]);
            __m256d w2i = _mm256_set_pd(tw[(j+3)*6+3], tw[(j+2)*6+3], tw[(j+1)*6+3], (j==0)?0.0:tw[j*6+3]);
            __m256d o2_re, o2_im;
            cmul_split_avx2(t2_re, t2_im, w2r, w2i, o2_re, o2_im);

            // out3 = (amc + jbmd) * w3[j]
            __m256d t3_re = _mm256_add_pd(amc_re, jbmd_re);
            __m256d t3_im = _mm256_add_pd(amc_im, jbmd_im);
            __m256d w3r = _mm256_set_pd(tw[(j+3)*6+4], tw[(j+2)*6+4], tw[(j+1)*6+4], (j==0)?1.0:tw[j*6+4]);
            __m256d w3i = _mm256_set_pd(tw[(j+3)*6+5], tw[(j+2)*6+5], tw[(j+1)*6+5], (j==0)?0.0:tw[j*6+5]);
            __m256d o3_re, o3_im;
            cmul_split_avx2(t3_re, t3_im, w3r, w3i, o3_re, o3_im);

            // 4×4 transpose: o0[j0,j1,j2,j3], o1[j0,j1,j2,j3], o2[...], o3[...]
            // → row0[o0_j0, o1_j0, o2_j0, o3_j0], row1[o0_j1, ...], ...
            // Using unpack + permute2f128
            // Step 1: interleave pairs
            __m256d t01_lo_re = _mm256_unpacklo_pd(o0_re, o1_re); // [o0_j0, o1_j0, o0_j2, o1_j2]
            __m256d t01_hi_re = _mm256_unpackhi_pd(o0_re, o1_re); // [o0_j1, o1_j1, o0_j3, o1_j3]
            __m256d t23_lo_re = _mm256_unpacklo_pd(o2_re, o3_re); // [o2_j0, o3_j0, o2_j2, o3_j2]
            __m256d t23_hi_re = _mm256_unpackhi_pd(o2_re, o3_re); // [o2_j1, o3_j1, o2_j3, o3_j3]
            // Step 2: combine 128-bit halves
            __m256d row0_re = _mm256_permute2f128_pd(t01_lo_re, t23_lo_re, 0x20); // [o0_j0, o1_j0, o2_j0, o3_j0]
            __m256d row1_re = _mm256_permute2f128_pd(t01_hi_re, t23_hi_re, 0x20); // [o0_j1, o1_j1, o2_j1, o3_j1]
            __m256d row2_re = _mm256_permute2f128_pd(t01_lo_re, t23_lo_re, 0x31); // [o0_j2, o1_j2, o2_j2, o3_j2]
            __m256d row3_re = _mm256_permute2f128_pd(t01_hi_re, t23_hi_re, 0x31); // [o0_j3, o1_j3, o2_j3, o3_j3]

            // Same for imaginary
            __m256d t01_lo_im = _mm256_unpacklo_pd(o0_im, o1_im);
            __m256d t01_hi_im = _mm256_unpackhi_pd(o0_im, o1_im);
            __m256d t23_lo_im = _mm256_unpacklo_pd(o2_im, o3_im);
            __m256d t23_hi_im = _mm256_unpackhi_pd(o2_im, o3_im);
            __m256d row0_im = _mm256_permute2f128_pd(t01_lo_im, t23_lo_im, 0x20);
            __m256d row1_im = _mm256_permute2f128_pd(t01_hi_im, t23_hi_im, 0x20);
            __m256d row2_im = _mm256_permute2f128_pd(t01_lo_im, t23_lo_im, 0x31);
            __m256d row3_im = _mm256_permute2f128_pd(t01_hi_im, t23_hi_im, 0x31);

            // Apply fused twist (with 2/N scaling) — 4 contiguous outputs per row
            int32_t ob = 4 * j;
            __m256d tc, ts, tr, ti;

            tc = _mm256_loadu_pd(tw_cos + ob); ts = _mm256_loadu_pd(tw_sin + ob);
            cmul_split_avx2(row0_re, row0_im, tc, ts, tr, ti);
            _mm256_storeu_pd(dst_re + ob, tr); _mm256_storeu_pd(dst_im + ob, ti);

            tc = _mm256_loadu_pd(tw_cos + ob + 4); ts = _mm256_loadu_pd(tw_sin + ob + 4);
            cmul_split_avx2(row1_re, row1_im, tc, ts, tr, ti);
            _mm256_storeu_pd(dst_re + ob + 4, tr); _mm256_storeu_pd(dst_im + ob + 4, ti);

            tc = _mm256_loadu_pd(tw_cos + ob + 8); ts = _mm256_loadu_pd(tw_sin + ob + 8);
            cmul_split_avx2(row2_re, row2_im, tc, ts, tr, ti);
            _mm256_storeu_pd(dst_re + ob + 8, tr); _mm256_storeu_pd(dst_im + ob + 8, ti);

            tc = _mm256_loadu_pd(tw_cos + ob + 12); ts = _mm256_loadu_pd(tw_sin + ob + 12);
            cmul_split_avx2(row3_re, row3_im, tc, ts, tr, ti);
            _mm256_storeu_pd(dst_re + ob + 12, tr); _mm256_storeu_pd(dst_im + ob + 12, ti);
        }
    }

    // If result ended up in scratch buffer, copy back
    if (dst_re != re) {
        memcpy(re, dst_re, ns4 * sizeof(double));
        memcpy(im, dst_im, ns4 * sizeof(double));
    }
}

// ── Full inverse FFT (fused twist + multi-pass) ─────────────────────────────
static void stockham_fft_inverse(
    int32_t ns4,
    const double *twist_tw,       // inverse twist [cos|sin], split format
    const double *butterfly_tw,   // butterfly twiddles
    double *re, double *im,
    double *sre, double *sim)
{
    double *cur_re = re, *cur_im = im;
    double *dst_re = sre, *dst_im = sim;
    int32_t q = ns4 / 4;
    int32_t s = 1;

    // First pass: q=ns4/4, s=1. Fuse inverse twist with input load.
    {
        const double *tw_cos = twist_tw;
        const double *tw_sin = twist_tw + ns4;
        int32_t stride = q; // q * s = q

        // j=0 pass: twist + butterfly (no butterfly twiddle for j=0)
        for (int32_t p = 0; p < q; p += 4) {
            // Apply inverse twist to inputs
            __m256d tc0 = _mm256_loadu_pd(tw_cos + p);
            __m256d ts0 = _mm256_loadu_pd(tw_sin + p);
            __m256d tc1 = _mm256_loadu_pd(tw_cos + p + stride);
            __m256d ts1 = _mm256_loadu_pd(tw_sin + p + stride);
            __m256d tc2 = _mm256_loadu_pd(tw_cos + p + 2*stride);
            __m256d ts2 = _mm256_loadu_pd(tw_sin + p + 2*stride);
            __m256d tc3 = _mm256_loadu_pd(tw_cos + p + 3*stride);
            __m256d ts3 = _mm256_loadu_pd(tw_sin + p + 3*stride);

            __m256d raw_a_re = _mm256_loadu_pd(cur_re + p);
            __m256d raw_a_im = _mm256_loadu_pd(cur_im + p);
            __m256d a_re, a_im;
            cmul_split_avx2(raw_a_re, raw_a_im, tc0, ts0, a_re, a_im);

            __m256d raw_b_re = _mm256_loadu_pd(cur_re + p + stride);
            __m256d raw_b_im = _mm256_loadu_pd(cur_im + p + stride);
            __m256d b_re, b_im;
            cmul_split_avx2(raw_b_re, raw_b_im, tc1, ts1, b_re, b_im);

            __m256d raw_c_re = _mm256_loadu_pd(cur_re + p + 2*stride);
            __m256d raw_c_im = _mm256_loadu_pd(cur_im + p + 2*stride);
            __m256d c_re, c_im;
            cmul_split_avx2(raw_c_re, raw_c_im, tc2, ts2, c_re, c_im);

            __m256d raw_d_re = _mm256_loadu_pd(cur_re + p + 3*stride);
            __m256d raw_d_im = _mm256_loadu_pd(cur_im + p + 3*stride);
            __m256d d_re, d_im;
            cmul_split_avx2(raw_d_re, raw_d_im, tc3, ts3, d_re, d_im);

            // Radix-4 inverse butterfly
            __m256d neg = _mm256_set1_pd(-0.0);
            __m256d apc_re = _mm256_add_pd(a_re, c_re);
            __m256d apc_im = _mm256_add_pd(a_im, c_im);
            __m256d amc_re = _mm256_sub_pd(a_re, c_re);
            __m256d amc_im = _mm256_sub_pd(a_im, c_im);
            __m256d bpd_re = _mm256_add_pd(b_re, d_re);
            __m256d bpd_im = _mm256_add_pd(b_im, d_im);
            __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
            __m256d bmd_im = _mm256_sub_pd(b_im, d_im);
            // Inverse j-multiply: [bmd_im, -bmd_re]
            __m256d jbmd_re = bmd_im;
            __m256d jbmd_im = _mm256_xor_pd(bmd_re, neg);

            _mm256_storeu_pd(dst_re + p, _mm256_add_pd(apc_re, bpd_re));
            _mm256_storeu_pd(dst_im + p, _mm256_add_pd(apc_im, bpd_im));
            _mm256_storeu_pd(dst_re + q + p, _mm256_sub_pd(amc_re, jbmd_re));
            _mm256_storeu_pd(dst_im + q + p, _mm256_sub_pd(amc_im, jbmd_im));
            _mm256_storeu_pd(dst_re + 2*q + p, _mm256_sub_pd(apc_re, bpd_re));
            _mm256_storeu_pd(dst_im + 2*q + p, _mm256_sub_pd(apc_im, bpd_im));
            _mm256_storeu_pd(dst_re + 3*q + p, _mm256_add_pd(amc_re, jbmd_re));
            _mm256_storeu_pd(dst_im + 3*q + p, _mm256_add_pd(amc_im, jbmd_im));
        }

        // Swap buffers
        double *t;
        t = cur_re; cur_re = dst_re; dst_re = t;
        t = cur_im; cur_im = dst_im; dst_im = t;

        q /= 4; // now q = ns4/16
        s = 4;
    }

    // Remaining passes: standard Stockham radix-4 inverse
    // Skip first stage's butterfly twiddles (s=1 → 1*6=6 doubles)
    const double *tw = butterfly_tw + 1 * 6;
    while (q >= 1) {
        if (q >= 4) {
            stockham_r4_inv_pass(ns4, q, s, tw, cur_re, cur_im, dst_re, dst_im);
        } else {
            inv_final_pass_vec(s, tw, cur_re, cur_im, dst_re, dst_im);
        }
        tw += s * 6;
        double *t;
        t = cur_re; cur_re = dst_re; dst_re = t;
        t = cur_im; cur_im = dst_im; dst_im = t;
        s *= 4; q /= 4;
    }

    // Copy back if needed
    if (cur_re != re) {
        memcpy(re, cur_re, ns4 * sizeof(double));
        memcpy(im, cur_im, ns4 * sizeof(double));
    }
}

// ── Table construction ──────────────────────────────────────────────────────
static STOCKHAM_R4_PRECOMP *build_tables(int32_t nn, bool is_forward) {
    int32_t n = 2 * nn;
    int32_t ns4 = nn / 2;

    auto *reps = new STOCKHAM_R4_PRECOMP;
    reps->n = n;
    reps->ns4 = ns4;

    // Count twiddle storage needed
    // Butterfly twiddles: for each stage, s entries × 6 doubles
    int32_t bf_doubles = 0;
    for (int32_t s = 1, q = ns4/4; q >= 1; s *= 4, q /= 4)
        bf_doubles += s * 6;

    // Twist twiddles: ns4 cos + ns4 sin = 2*ns4 doubles
    int32_t twist_doubles = 2 * ns4;

    int32_t total = bf_doubles + twist_doubles + 2 * ns4 + nn;
    reps->buf = aligned_alloc(64, total * sizeof(double));
    double *ptr = (double *)reps->buf;

    if (is_forward) {
        reps->trig_fwd = ptr;
        // Build butterfly twiddles
        double *bfp = ptr;
        for (int32_t s = 1, q = ns4/4; q >= 1; s *= 4, q /= 4) {
            for (int32_t j = 0; j < s; j++) {
                int32_t denom = 4 * s;
                // w1 = exp(-2πij/(4s)), w2 = exp(-2πi·2j/(4s)), w3 = exp(-2πi·3j/(4s))
                *bfp++ = accurate_cos(-j, denom);
                *bfp++ = accurate_sin(-j, denom);
                *bfp++ = accurate_cos(-2*j, denom);
                *bfp++ = accurate_sin(-2*j, denom);
                *bfp++ = accurate_cos(-3*j, denom);
                *bfp++ = accurate_sin(-3*j, denom);
            }
        }
        // Build twist twiddles with 2/N scaling baked in
        double scale = 2.0 / nn;
        double *twist = bfp;
        // cos part
        for (int32_t j = 0; j < ns4; j++)
            *bfp++ = scale * accurate_cos(-j, n);
        // sin part
        for (int32_t j = 0; j < ns4; j++)
            *bfp++ = scale * accurate_sin(-j, n);

        ptr = bfp;
    } else {
        reps->trig_inv = ptr;
        // Inverse twist twiddles (no scaling)
        double *twist = ptr;
        // cos part
        for (int32_t j = 0; j < ns4; j++)
            *ptr++ = accurate_cos(j, n);
        // sin part
        for (int32_t j = 0; j < ns4; j++)
            *ptr++ = accurate_sin(j, n);
        // Build butterfly twiddles (conjugate direction)
        for (int32_t s = 1, q = ns4/4; q >= 1; s *= 4, q /= 4) {
            for (int32_t j = 0; j < s; j++) {
                int32_t denom = 4 * s;
                *ptr++ = accurate_cos(j, denom);
                *ptr++ = accurate_sin(j, denom);
                *ptr++ = accurate_cos(2*j, denom);
                *ptr++ = accurate_sin(2*j, denom);
                *ptr++ = accurate_cos(3*j, denom);
                *ptr++ = accurate_sin(3*j, denom);
            }
        }
    }

    reps->scratch = ptr;
    ptr += ns4;         // scratch re
    ptr += ns4;         // scratch im (contiguous)
    reps->data = ptr;   // nn doubles for data buffer

    return reps;
}

// ── C interface ─────────────────────────────────────────────────────────────
extern "C" {

void *new_fft_table(int32_t nn) {
    return build_tables(nn, true);
}

void *new_ifft_table(int32_t nn) {
    return build_tables(nn, false);
}

double *fft_table_get_buffer(const void *tables) {
    return ((STOCKHAM_R4_PRECOMP *)tables)->data;
}

double *ifft_table_get_buffer(const void *tables) {
    return ((STOCKHAM_R4_PRECOMP *)tables)->data;
}

void fft(const void *tables, double *c) {
    auto *t = (STOCKHAM_R4_PRECOMP *)tables;
    int32_t ns4 = t->ns4;
    double *re = c;
    double *im = c + ns4;

    // Butterfly twiddles start at trig_fwd
    const double *bf_tw = t->trig_fwd;
    // Count butterfly doubles to find twist start
    int32_t bf_doubles = 0;
    for (int32_t s = 1, q = ns4/4; q >= 1; s *= 4, q /= 4)
        bf_doubles += s * 6;
    const double *twist_tw = t->trig_fwd + bf_doubles;

    stockham_fft_forward(ns4, bf_tw, twist_tw, re, im, t->scratch, t->scratch + ns4);
}

void ifft(const void *tables, double *c) {
    auto *t = (STOCKHAM_R4_PRECOMP *)tables;
    int32_t ns4 = t->ns4;
    double *re = c;
    double *im = c + ns4;

    const double *twist_tw = t->trig_inv;
    const double *bf_tw = t->trig_inv + 2 * ns4;

    stockham_fft_inverse(ns4, twist_tw, bf_tw, re, im, t->scratch, t->scratch + ns4);
}

// Fused IFFT from int32_t input: converts int32→double and does IFFT in
// one fewer pass by fusing the conversion into the first butterfly+twist pass.
void ifft_from_i32(const void *tables, double *out, const int32_t *in_re, const int32_t *in_im) {
    auto *t = (STOCKHAM_R4_PRECOMP *)tables;
    int32_t ns4 = t->ns4;
    double *dst_re = t->scratch;
    double *dst_im = t->scratch + ns4;
    const double *twist_tw = t->trig_inv;
    const double *tw_cos = twist_tw;
    const double *tw_sin = twist_tw + ns4;

    int32_t q = ns4 / 4;
    int32_t stride = q;

    // Fused first pass: convert int32→double + inverse twist + radix-4 butterfly
    for (int32_t p = 0; p < q; p += 4) {
        // Load twist twiddles
        __m256d tc0 = _mm256_loadu_pd(tw_cos + p);
        __m256d ts0 = _mm256_loadu_pd(tw_sin + p);
        __m256d tc1 = _mm256_loadu_pd(tw_cos + p + stride);
        __m256d ts1 = _mm256_loadu_pd(tw_sin + p + stride);
        __m256d tc2 = _mm256_loadu_pd(tw_cos + p + 2*stride);
        __m256d ts2 = _mm256_loadu_pd(tw_sin + p + 2*stride);
        __m256d tc3 = _mm256_loadu_pd(tw_cos + p + 3*stride);
        __m256d ts3 = _mm256_loadu_pd(tw_sin + p + 3*stride);

        // Convert int32→double and apply twist in one step
        __m128i ia = _mm_loadu_si128((const __m128i*)(in_re + p));
        __m256d raw_a_re = _mm256_cvtepi32_pd(ia);
        ia = _mm_loadu_si128((const __m128i*)(in_im + p));
        __m256d raw_a_im = _mm256_cvtepi32_pd(ia);
        __m256d a_re, a_im;
        cmul_split_avx2(raw_a_re, raw_a_im, tc0, ts0, a_re, a_im);

        __m128i ib = _mm_loadu_si128((const __m128i*)(in_re + p + stride));
        __m256d raw_b_re = _mm256_cvtepi32_pd(ib);
        ib = _mm_loadu_si128((const __m128i*)(in_im + p + stride));
        __m256d raw_b_im = _mm256_cvtepi32_pd(ib);
        __m256d b_re, b_im;
        cmul_split_avx2(raw_b_re, raw_b_im, tc1, ts1, b_re, b_im);

        __m128i ic = _mm_loadu_si128((const __m128i*)(in_re + p + 2*stride));
        __m256d raw_c_re = _mm256_cvtepi32_pd(ic);
        ic = _mm_loadu_si128((const __m128i*)(in_im + p + 2*stride));
        __m256d raw_c_im = _mm256_cvtepi32_pd(ic);
        __m256d c_re, c_im;
        cmul_split_avx2(raw_c_re, raw_c_im, tc2, ts2, c_re, c_im);

        __m128i id = _mm_loadu_si128((const __m128i*)(in_re + p + 3*stride));
        __m256d raw_d_re = _mm256_cvtepi32_pd(id);
        id = _mm_loadu_si128((const __m128i*)(in_im + p + 3*stride));
        __m256d raw_d_im = _mm256_cvtepi32_pd(id);
        __m256d d_re, d_im;
        cmul_split_avx2(raw_d_re, raw_d_im, tc3, ts3, d_re, d_im);

        // Radix-4 inverse butterfly
        __m256d neg = _mm256_set1_pd(-0.0);
        __m256d apc_re = _mm256_add_pd(a_re, c_re);
        __m256d apc_im = _mm256_add_pd(a_im, c_im);
        __m256d amc_re = _mm256_sub_pd(a_re, c_re);
        __m256d amc_im = _mm256_sub_pd(a_im, c_im);
        __m256d bpd_re = _mm256_add_pd(b_re, d_re);
        __m256d bpd_im = _mm256_add_pd(b_im, d_im);
        __m256d bmd_re = _mm256_sub_pd(b_re, d_re);
        __m256d bmd_im = _mm256_sub_pd(b_im, d_im);
        __m256d jbmd_re = bmd_im;
        __m256d jbmd_im = _mm256_xor_pd(bmd_re, neg);

        _mm256_storeu_pd(dst_re + p, _mm256_add_pd(apc_re, bpd_re));
        _mm256_storeu_pd(dst_im + p, _mm256_add_pd(apc_im, bpd_im));
        _mm256_storeu_pd(dst_re + q + p, _mm256_sub_pd(amc_re, jbmd_re));
        _mm256_storeu_pd(dst_im + q + p, _mm256_sub_pd(amc_im, jbmd_im));
        _mm256_storeu_pd(dst_re + 2*q + p, _mm256_sub_pd(apc_re, bpd_re));
        _mm256_storeu_pd(dst_im + 2*q + p, _mm256_sub_pd(apc_im, bpd_im));
        _mm256_storeu_pd(dst_re + 3*q + p, _mm256_add_pd(amc_re, jbmd_re));
        _mm256_storeu_pd(dst_im + 3*q + p, _mm256_add_pd(amc_im, jbmd_im));
    }

    // Remaining passes
    double *cur_re = dst_re, *cur_im = dst_im;
    dst_re = out; dst_im = out + ns4;
    q /= 4;
    int32_t s = 4;
    const double *bf_tw = t->trig_inv + 2 * ns4 + 1 * 6; // skip first stage twiddles

    while (q >= 1) {
        if (q >= 4) {
            stockham_r4_inv_pass(ns4, q, s, bf_tw, cur_re, cur_im, dst_re, dst_im);
        } else {
            // q=1 final pass (scalar gather, contiguous store)
            for (int32_t j = 0; j < s; j++) {
                double a_r = cur_re[j], a_i = cur_im[j];
                double b_r = cur_re[j+s], b_i = cur_im[j+s];
                double c_r = cur_re[j+2*s], c_i = cur_im[j+2*s];
                double d_r = cur_re[j+3*s], d_i = cur_im[j+3*s];
                double apc_r = a_r+c_r, apc_i = a_i+c_i;
                double amc_r = a_r-c_r, amc_i = a_i-c_i;
                double bpd_r = b_r+d_r, bpd_i = b_i+d_i;
                double bmd_r = b_r-d_r, bmd_i = b_i-d_i;
                double jbmd_r = bmd_i, jbmd_i = -bmd_r;
                int32_t ob = 4*j;
                dst_re[ob] = apc_r + bpd_r; dst_im[ob] = apc_i + bpd_i;
                double o1_r = amc_r-jbmd_r, o1_i = amc_i-jbmd_i;
                double o2_r = apc_r-bpd_r, o2_i = apc_i-bpd_i;
                double o3_r = amc_r+jbmd_r, o3_i = amc_i+jbmd_i;
                if (j != 0) {
                    const double *twj = bf_tw + j*6;
                    double t;
                    t = o1_r*twj[0]-o1_i*twj[1]; o1_i = o1_i*twj[0]+o1_r*twj[1]; o1_r = t;
                    t = o2_r*twj[2]-o2_i*twj[3]; o2_i = o2_i*twj[2]+o2_r*twj[3]; o2_r = t;
                    t = o3_r*twj[4]-o3_i*twj[5]; o3_i = o3_i*twj[4]+o3_r*twj[5]; o3_r = t;
                }
                dst_re[ob+1] = o1_r; dst_im[ob+1] = o1_i;
                dst_re[ob+2] = o2_r; dst_im[ob+2] = o2_i;
                dst_re[ob+3] = o3_r; dst_im[ob+3] = o3_i;
            }
        }
        bf_tw += s * 6;
        double *tmp;
        tmp = cur_re; cur_re = dst_re; dst_re = tmp;
        tmp = cur_im; cur_im = dst_im; dst_im = tmp;
        s *= 4; q /= 4;
    }

    // Copy to output if needed
    if (cur_re != out) {
        memcpy(out, cur_re, ns4 * sizeof(double));
        memcpy(out + ns4, cur_im, ns4 * sizeof(double));
    }
}

// Model functions (for debugging, same interface)
void fft_model(const void *tables) {
    auto *t = (STOCKHAM_R4_PRECOMP *)tables;
    fft(tables, t->data);
}

void ifft_model(void *tables) {
    auto *t = (STOCKHAM_R4_PRECOMP *)tables;
    ifft(tables, t->data);
}

} // extern "C"

// ── FFT processor ───────────────────────────────────────────────────────────
static int32_t rev(int32_t x, int32_t M) {
    int32_t r = 0;
    for (int32_t j = M; j > 1; j /= 2) { r = 2*r+(x%2); x /= 2; }
    return r;
}

FFT_Processor_Spqlios::FFT_Processor_Spqlios(const int32_t N)
    : _2N(2*N), N(N), Ns2(N/2) {
    tables_direct = new_fft_table(N);
    tables_reverse = new_ifft_table(N);
    real_inout_direct = fft_table_get_buffer(tables_direct);
    imag_inout_direct = real_inout_direct + Ns2;
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
    auto *ft = (STOCKHAM_R4_PRECOMP *)tables_direct;
    auto *it = (STOCKHAM_R4_PRECOMP *)tables_reverse;
    free(ft->buf); delete ft;
    if (it != ft) { free(it->buf); delete it; }
    delete[] cosomegaxminus1;
    delete[] reva;
}

thread_local FFT_Processor_Spqlios fftplvl1(TFHEpp::lvl1param::n);
thread_local FFT_Processor_Spqlios fftplvl2(TFHEpp::lvl2param::n);
thread_local FFT_Processor_Spqlios fftplvl3(TFHEpp::lvl3param::n);
thread_local FFT_Processor_Spqlios fftplvl5(1 << 14);
