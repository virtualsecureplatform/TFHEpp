#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>

// Interleaved-format FFT processor.
// Data layout: [re0, im0, re1, im1, ...] — N doubles total for N/2 complex values.
// This matches the USE_INTERLEAVED_FORMAT convention in TFHEpp.

class FFT_Processor_Spqlios_Intl {
public:
    const int32_t _2N;
    const int32_t N;
    const int32_t Ns2;

private:
    double *real_inout_direct;  // scratch buffer for forward FFT (N doubles)
    void *tables_direct;        // FFT trig tables (interleaved format)
    void *tables_reverse;       // IFFT trig tables (interleaved format)
public:
    double *cosomegaxminus1;
    double *sinomegaxminus1;
    int32_t *reva;

    FFT_Processor_Spqlios_Intl(const int32_t N);

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

    ~FFT_Processor_Spqlios_Intl();
};

extern thread_local FFT_Processor_Spqlios_Intl fftplvl1;
extern thread_local FFT_Processor_Spqlios_Intl fftplvl2;
extern thread_local FFT_Processor_Spqlios_Intl fftplvl3;
extern thread_local FFT_Processor_Spqlios_Intl fftplvl5;
