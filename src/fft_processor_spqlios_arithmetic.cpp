#include <fft_processor_spqlios_arithmetic.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <params.hpp>

extern "C" {
#include "../thirdparties/spqlios-arithmetic/spqlios/reim/reim_fft.h"
}

using namespace std;

static int32_t rev(int32_t x, int32_t M)
{
    int32_t reps = 0;
    for (int32_t j = M; j > 1; j /= 2) {
        reps = 2 * reps + (x % 2);
        x /= 2;
    }
    return reps;
}

FFT_Processor_Spqlios_Arithmetic::FFT_Processor_Spqlios_Arithmetic(
    const int32_t N)
    : _2N(2 * N), N(N), Ns2(N / 2)
{
    // Initialize FFT and IFFT tables (may not be needed for simple API)
    tables_direct = nullptr;
    tables_reverse = nullptr;

    // Allocate buffers
    buffer_direct = new double[N];
    buffer_reverse = new double[N];

    // Initialize auxiliary arrays
    reva = new int32_t[Ns2];
    cosomegaxminus1 = new double[2 * _2N];
    sinomegaxminus1 = cosomegaxminus1 + _2N;

    // Compute reversal indices
    int32_t rev1 = rev(1, _2N);
    int32_t rev3 = rev(3, _2N);
    for (int32_t revi = rev1; revi < rev3; revi++)
        reva[revi - rev1] = rev(revi, _2N);

    // Compute twiddle factors
    for (int32_t j = 0; j < _2N; j++) {
        cosomegaxminus1[j] = cos(2 * M_PI * j / _2N) - 1.;
        sinomegaxminus1[j] = sin(2 * M_PI * j / _2N);
    }
}

void FFT_Processor_Spqlios_Arithmetic::execute_reverse_uint(double* res,
                                                            const uint32_t* a)
{
    // Convert uint32_t to double
    for (int32_t i = 0; i < N; i++) {
        res[i] = (double)a[i];
    }

    // Execute IFFT using spqlios-arithmetic simple API
    reim_ifft_simple(Ns2, res);
}

void FFT_Processor_Spqlios_Arithmetic::execute_reverse_int(double* res,
                                                           const int32_t* a)
{
    // Convert int32_t to double
    for (int32_t i = 0; i < N; i++) {
        res[i] = (double)a[i];
    }

    // Execute IFFT using simple API
    reim_ifft_simple(Ns2, res);
}

void FFT_Processor_Spqlios_Arithmetic::execute_reverse_torus32(
    double* res, const uint32_t* a)
{
    int32_t* aa = (int32_t*)a;
    execute_reverse_int(res, aa);
}

void FFT_Processor_Spqlios_Arithmetic::execute_reverse_torus64(
    double* res, const uint64_t* a)
{
    // Convert int64_t to double
    int64_t* aa = (int64_t*)a;
    for (int i = 0; i < N; i++) {
        res[i] = (double)aa[i];
    }

    // Execute IFFT using simple API
    reim_ifft_simple(Ns2, res);
}

void FFT_Processor_Spqlios_Arithmetic::execute_direct_torus32(uint32_t* res,
                                                              const double* a)
{
    static const double _2sN = double(2) / double(N);

    // Copy and scale input to buffer
    for (int32_t i = 0; i < N; i++) {
        buffer_direct[i] = a[i] * _2sN;
    }

    // Execute FFT using simple API
    reim_fft_simple(Ns2, buffer_direct);

    // Convert back to torus32
    for (int32_t i = 0; i < Ns2; i++) {
        double re = buffer_direct[i];
        double im = buffer_direct[i + Ns2];

        // Twist
        double temp0 = re * cosomegaxminus1[i] - im * sinomegaxminus1[i];
        double temp1 = re * sinomegaxminus1[i] + im * cosomegaxminus1[i];
        temp0 += re;

        int32_t t0 = int32_t(int64_t(temp0 * (1L << 32)));
        int32_t t1 = int32_t(int64_t(temp1 * (1L << 32)));
        res[i] = t0;
        res[i + Ns2] = t1;
    }
}

void FFT_Processor_Spqlios_Arithmetic::execute_direct_torus32_q(
    uint32_t* res, const double* a, const uint32_t q)
{
    static const double _2sN = double(2) / double(N);

    // Copy and scale input to buffer
    for (int32_t i = 0; i < N; i++) {
        buffer_direct[i] = a[i] * _2sN;
    }

    // Execute FFT using simple API
    reim_fft_simple(Ns2, buffer_direct);

    // Convert back to torus32 with modulus q
    for (int32_t i = 0; i < Ns2; i++) {
        double re = buffer_direct[i];
        double im = buffer_direct[i + Ns2];

        // Twist
        double temp0 = re * cosomegaxminus1[i] - im * sinomegaxminus1[i];
        double temp1 = re * sinomegaxminus1[i] + im * cosomegaxminus1[i];
        temp0 += re;

        int64_t t0 = int64_t(temp0 * double(q));
        int64_t t1 = int64_t(temp1 * double(q));
        res[i] = uint32_t(t0 % q);
        res[i + Ns2] = uint32_t(t1 % q);
    }
}

void FFT_Processor_Spqlios_Arithmetic::execute_direct_torus32_rescale(
    uint32_t* res, const double* a, const double Δ)
{
    static const double _2sN = double(2) / double(N);

    // Copy and scale input to buffer
    for (int32_t i = 0; i < N; i++) {
        buffer_direct[i] = a[i] * _2sN;
    }

    // Execute FFT using simple API
    reim_fft_simple(Ns2, buffer_direct);

    // Convert back to torus32 with rescaling
    for (int32_t i = 0; i < Ns2; i++) {
        double re = buffer_direct[i];
        double im = buffer_direct[i + Ns2];

        // Twist
        double temp0 = re * cosomegaxminus1[i] - im * sinomegaxminus1[i];
        double temp1 = re * sinomegaxminus1[i] + im * cosomegaxminus1[i];
        temp0 += re;

        res[i] = uint32_t(int32_t(temp0 / Δ));
        res[i + Ns2] = uint32_t(int32_t(temp1 / Δ));
    }
}

void FFT_Processor_Spqlios_Arithmetic::execute_direct_torus64(uint64_t* res,
                                                              const double* a)
{
    static const double _2sN = double(2) / double(N);

    // Copy and scale input to buffer
    for (int32_t i = 0; i < N; i++) {
        buffer_direct[i] = a[i] * _2sN;
    }

    // Execute FFT using simple API
    reim_fft_simple(Ns2, buffer_direct);

    // Convert back to torus64
    for (int32_t i = 0; i < Ns2; i++) {
        double re = buffer_direct[i];
        double im = buffer_direct[i + Ns2];

        // Twist
        double temp0 = re * cosomegaxminus1[i] - im * sinomegaxminus1[i];
        double temp1 = re * sinomegaxminus1[i] + im * cosomegaxminus1[i];
        temp0 += re;

        // Scale to 64-bit torus
        int64_t t0 = int64_t(temp0 * double(1ULL << 63)) * 2;
        int64_t t1 = int64_t(temp1 * double(1ULL << 63)) * 2;
        res[i] = uint64_t(t0);
        res[i + Ns2] = uint64_t(t1);
    }
}

void FFT_Processor_Spqlios_Arithmetic::execute_direct_torus64_rescale(
    uint64_t* res, const double* a, const double Δ)
{
    static const double _2sN = double(2) / double(N);

    // Copy and scale input to buffer
    for (int32_t i = 0; i < N; i++) {
        buffer_direct[i] = a[i] * _2sN;
    }

    // Execute FFT using simple API
    reim_fft_simple(Ns2, buffer_direct);

    // Convert back to torus64 with rescaling
    for (int32_t i = 0; i < Ns2; i++) {
        double re = buffer_direct[i];
        double im = buffer_direct[i + Ns2];

        // Twist
        double temp0 = re * cosomegaxminus1[i] - im * sinomegaxminus1[i];
        double temp1 = re * sinomegaxminus1[i] + im * cosomegaxminus1[i];
        temp0 += re;

        res[i] = uint64_t(int64_t(temp0 / Δ));
        res[i + Ns2] = uint64_t(int64_t(temp1 / Δ));
    }
}

FFT_Processor_Spqlios_Arithmetic::~FFT_Processor_Spqlios_Arithmetic()
{
    delete[] buffer_direct;
    delete[] buffer_reverse;
    delete[] reva;
    delete[] cosomegaxminus1;
}

// Thread-local FFT processors for level 1 and level 2
thread_local FFT_Processor_Spqlios_Arithmetic fftplvl1(TFHEpp::lvl1param::n);
thread_local FFT_Processor_Spqlios_Arithmetic fftplvl2(TFHEpp::lvl2param::n);