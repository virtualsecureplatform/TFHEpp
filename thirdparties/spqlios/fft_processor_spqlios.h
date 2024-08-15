#pragma once

#include <cassert>
#include <cmath>
#include<cstdint>
extern "C"{
    void *new_fft_table(int32_t nn);
    double *fft_table_get_buffer(const void *tables);
    void *new_ifft_table(int32_t nn);
    double *ifft_table_get_buffer(const void *tables);
    void fft_model(const void *tables);
    void ifft_model(void *tables);
    void fft(const void *tables, double *data);
    void ifft(const void *tables, double *data);
}

class FFT_Processor_Spqlios {
public:
    const int32_t _2N;
    const int32_t N;
    const int32_t Ns2;

private:
    double *real_inout_direct;
    double *imag_inout_direct;
    void *tables_direct;
    void *tables_reverse;
public:
    double *cosomegaxminus1;
    double *sinomegaxminus1;
    int32_t *reva; //rev(2i+1,_2N)

    FFT_Processor_Spqlios(const int32_t N);

    void execute_reverse_int(double *res, const int32_t *a);

    void execute_reverse_uint(double *res, const uint32_t *a);

    void execute_reverse_torus32(double *res, const uint32_t *a);

    void execute_direct_torus32(uint32_t *res, const double *a);
    
    void execute_direct_torus32_q(uint32_t *res, const double *a, const uint32_t q);

    void execute_direct_torus32_rescale(uint32_t *res, const double *a, const double Δ);

    void execute_reverse_torus64(double* res, const uint64_t* a);
    
    void execute_direct_torus64(uint64_t* res, const double* a);

    void execute_direct_torus64_rescale(uint64_t* res, const double* a, const double Δ);

    ~FFT_Processor_Spqlios();
};

extern thread_local FFT_Processor_Spqlios fftplvl1;
extern thread_local FFT_Processor_Spqlios fftplvl2;