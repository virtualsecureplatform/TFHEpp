#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>

extern "C" {
#include "../thirdparties/spqlios-arithmetic/spqlios/reim/reim_fft.h"
}

class FFT_Processor_Spqlios_Arithmetic {
public:
    const int32_t _2N;
    const int32_t N;
    const int32_t Ns2;

private:
    REIM_FFT_PRECOMP* tables_direct;
    REIM_IFFT_PRECOMP* tables_reverse;
    double* buffer_direct;
    double* buffer_reverse;

public:
    double* cosomegaxminus1;
    double* sinomegaxminus1;
    int32_t* reva; // rev(2i+1,_2N)

    FFT_Processor_Spqlios_Arithmetic(const int32_t N);

    void execute_reverse_int(double* res, const int32_t* a);

    void execute_reverse_uint(double* res, const uint32_t* a);

    void execute_reverse_torus32(double* res, const uint32_t* a);

    void execute_direct_torus32(uint32_t* res, const double* a);

    void execute_direct_torus32_q(uint32_t* res, const double* a, const uint32_t q);

    void execute_direct_torus32_rescale(uint32_t* res, const double* a, const double Δ);

    void execute_reverse_torus64(double* res, const uint64_t* a);

    void execute_direct_torus64(uint64_t* res, const double* a);

    void execute_direct_torus64_rescale(uint64_t* res, const double* a, const double Δ);

    ~FFT_Processor_Spqlios_Arithmetic();
};

extern thread_local FFT_Processor_Spqlios_Arithmetic fftplvl1;
extern thread_local FFT_Processor_Spqlios_Arithmetic fftplvl2;