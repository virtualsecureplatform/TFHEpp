#pragma once

#include <fftw3.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

class FFT_Processor_FFTW {
public:
    const int32_t _2N;
    const int32_t N;
    const int32_t Ns2;

private:
    std::vector<std::complex<double>> twist;
    fftw_plan plan_forward;
    fftw_plan plan_backward;
    fftw_complex *inbuf;
    fftw_complex *outbuf;

public:
    FFT_Processor_FFTW(const int32_t N);

    void execute_reverse_int(double *res, const int32_t *a);

    void execute_reverse_uint(double *res, const uint32_t *a);

    void execute_reverse_torus32(double *res, const uint32_t *a);

    void execute_direct_torus32(uint32_t *res, const double *a);

    void execute_direct_torus32_rescale(uint32_t *res, const double *a,
                                        const double Δ);

    void execute_reverse_torus64(double *res, const uint64_t *a);

    void execute_direct_torus64(uint64_t *res, const double *a);

    void execute_direct_torus64_rescale(uint64_t *res, const double *a,
                                        const double Δ);

    ~FFT_Processor_FFTW();
};

// FFT_Processor_FFTW is thread-safe
extern FFT_Processor_FFTW fftplvl1;
extern FFT_Processor_FFTW fftplvl2;