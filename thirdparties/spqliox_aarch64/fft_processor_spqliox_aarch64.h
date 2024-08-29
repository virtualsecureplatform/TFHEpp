#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>

#include "spqliox_aarch64_impl.h"

struct alignas(32) FFT_PRECOMP {
    uint64_t n;
    double *trig_tables;
};

class FFT_Processor_Spqliox_AArch64 {
public:
    const int32_t _2N;
    const int32_t N;
    const int32_t Ns2;

private:
    double *real_inout;
    FFT_PRECOMP tables_direct_, tables_reverse_, table_negation_forward_,
        table_negation_reverse_;
    unsigned long (*fft_)(double *, const double *, const double *,
                          const double *, const FFT_PRECOMP *, double *,
                          double *);
    unsigned long (*ifft_)(double *, const double *, const double *,
                           const FFT_PRECOMP *, double *, double *);
    spqliox_aarch64::fft_code fft_code_;
    spqliox_aarch64::ifft_code ifft_code_;

public:
    FFT_Processor_Spqliox_AArch64(const int32_t N);

    void execute_reverse_int(double *res, const int32_t *a);

    void execute_reverse_uint(double *res, const uint32_t *a);

    void execute_reverse_torus32(double *res, const uint32_t *a);

    void execute_direct_torus32(uint32_t *res, const double *a);

    void execute_direct_torus32_rescale(uint32_t *res, const double *a,
                                        const double Δ);

    void execute_direct_torus64_rescale(uint64_t *res, const double *a,
                                        const double Δ);

    void execute_reverse_torus64(double *res, const uint64_t *a);

    void execute_direct_torus64(uint64_t *res, const double *a);

    ~FFT_Processor_Spqliox_AArch64();
};

extern thread_local FFT_Processor_Spqliox_AArch64 fftplvl1;
extern thread_local FFT_Processor_Spqliox_AArch64 fftplvl2;
