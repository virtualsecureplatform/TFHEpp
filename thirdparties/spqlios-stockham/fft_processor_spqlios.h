#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>

// Stockham radix-4 FFT processor with split real/imaginary layout.
// Same data layout as original SPQLIOS: re[0..ns2-1], im[0..ns2-1]
// Same external interface — drop-in replacement.

class FFT_Processor_Spqlios {
public:
    const int32_t _2N;
    const int32_t N;
    const int32_t Ns2;

private:
    double *real_inout_direct;
    void *tables_direct;
    void *tables_reverse;
public:
    double *cosomegaxminus1;
    double *sinomegaxminus1;
    int32_t *reva;

    FFT_Processor_Spqlios(const int32_t N);

    void execute_reverse_int(double *res, const int32_t *a);
    void execute_reverse_uint(double *res, const uint32_t *a);
    void execute_reverse_torus32(double *res, const uint32_t *a);
    void execute_direct_torus32(uint32_t *res, const double *a);
    void execute_direct_torus32_q(uint32_t *res, const double *a, const uint32_t q);
    void execute_direct_torus32_rescale(uint32_t *res, const double *a, const double Δ);
    void execute_direct_torus32_rescale_clpx(
        uint32_t *res, const double *a, const double q,
        const uint32_t plain_modulus);
    void execute_reverse_torus64(double *res, const uint64_t *a);
    void execute_reverse_torus64_uint(double *res, const uint64_t *a);
    void execute_direct_torus64(uint64_t *res, double *a);
    void execute_direct_torus32_add(uint32_t *res, const double *a);
    void execute_direct_torus64_add(uint64_t *res, double *a);
    void execute_direct_torus64_rescale(uint64_t *res, const double *a, const double Δ);
    void execute_direct_torus64_rescale_clpx(
        uint64_t *res, const double *a, const uint32_t plain_modulus);

    ~FFT_Processor_Spqlios();
};

extern thread_local FFT_Processor_Spqlios fftplvl1;
extern thread_local FFT_Processor_Spqlios fftplvl2;
extern thread_local FFT_Processor_Spqlios fftplvl3;
extern thread_local FFT_Processor_Spqlios fftplvl5;
extern thread_local FFT_Processor_Spqlios fftplvl6;
