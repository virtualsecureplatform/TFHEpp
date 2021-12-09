#include "fft_processor_fftw.h"

#include <fftw3.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <params.hpp>
#include <vector>

#define CAST_DOUBLE_TO_UINT32(d) ((uint32_t)((int64_t)(d)))

FFT_Processor_FFTW::FFT_Processor_FFTW(const int32_t N)
    : _2N(2 * N), N(N), Ns2(N / 2)
{
    auto in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);
    auto out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);
    plan_forward = fftw_plan_dft_1d(Ns2, in, out, FFTW_FORWARD, FFTW_MEASURE);
    plan_backward = fftw_plan_dft_1d(Ns2, in, out, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_free(in);
    fftw_free(out);

    for (int i = 0; i < Ns2; i++) {
        double value = (double)i * M_PI / (double)N;
        twist.push_back(std::complex<double>(std::cos(value), std::sin(value)));
    }
}

void FFT_Processor_FFTW::execute_reverse_int(double *res, const int32_t *a)
{
    auto in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);
    auto out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);

    for (int i = 0; i < Ns2; i++) {
        auto tmp = twist[i] * std::complex((double)a[i], (double)a[Ns2 + i]);
        in[i][0] = tmp.real();
        in[i][1] = tmp.imag();
    }
    fftw_execute_dft(plan_forward, in, out);
    for (int i = 0; i < Ns2; i++) {
        res[i] = out[i][0];
        res[i + Ns2] = out[i][1];
    }

    fftw_free(in);
    fftw_free(out);
}

void FFT_Processor_FFTW::execute_reverse_torus32(double *res, const uint32_t *a)
{
    int32_t *aa = (int32_t *)a;
    execute_reverse_int(res, aa);
}

void FFT_Processor_FFTW::execute_reverse_torus64(double *res, const uint64_t *a)
{
    auto in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);
    auto out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);

    for (int i = 0; i < Ns2; i++) {
        auto tmp = twist[i] * std::complex((double)((int64_t)a[i]),
                                           (double)((int64_t)a[Ns2 + i]));
        in[i][0] = tmp.real();
        in[i][1] = tmp.imag();
    }
    fftw_execute_dft(plan_forward, in, out);
    for (int i = 0; i < Ns2; i++) {
        res[i] = out[i][0];
        res[i + Ns2] = out[i][1];
    }

    fftw_free(in);
    fftw_free(out);
}

void FFT_Processor_FFTW::execute_direct_torus32(uint32_t *res, const double *a)
{
    auto in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);
    auto out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);

    for (int i = 0; i < Ns2; i++) {
        in[i][0] = a[i] / Ns2;
        in[i][1] = a[Ns2 + i] / Ns2;
    }
    fftw_execute_dft(plan_backward, in, out);
    for (int i = 0; i < Ns2; i++) {
        auto res_tmp =
            std::complex<double>(out[i][0], out[i][1]) * std::conj(twist[i]);
        res[i] = CAST_DOUBLE_TO_UINT32(res_tmp.real());
        res[i + Ns2] = CAST_DOUBLE_TO_UINT32(res_tmp.imag());
    }

    fftw_free(in);
    fftw_free(out);
}

void FFT_Processor_FFTW::execute_direct_torus32_rescale(uint32_t *res,
                                                        const double *a,
                                                        const double Δ)
{
    auto in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);
    auto out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);

    for (int i = 0; i < Ns2; i++) {
        in[i][0] = a[i] / Ns2;
        in[i][1] = a[Ns2 + i] / Ns2;
    }
    fftw_execute_dft(plan_backward, in, out);
    for (int i = 0; i < Ns2; i++) {
        auto res_tmp =
            std::complex<double>(out[i][0], out[i][1]) * std::conj(twist[i]);
        res[i] = CAST_DOUBLE_TO_UINT32(res_tmp.real() / (Δ / 4));
        res[i + Ns2] = CAST_DOUBLE_TO_UINT32(res_tmp.imag() / (Δ / 4));
    }

    fftw_free(in);
    fftw_free(out);
}

void FFT_Processor_FFTW::execute_direct_torus64(uint64_t *res, const double *a)
{
    auto in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);
    auto out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Ns2);

    for (int i = 0; i < Ns2; i++) {
        in[i][0] = a[i] / Ns2;
        in[i][1] = a[Ns2 + i] / Ns2;
    }
    fftw_execute_dft(plan_backward, in, out);
    static const uint64_t valmask0 = 0x000FFFFFFFFFFFFFul;
    static const uint64_t valmask1 = 0x0010000000000000ul;
    static const uint16_t expmask0 = 0x07FFu;
    double tmp[N];
    for (int i = 0; i < Ns2; i++) {
        auto res_tmp =
            std::complex<double>(out[i][0], out[i][1]) * std::conj(twist[i]);
        tmp[i] = res_tmp.real();
        tmp[i + Ns2] = res_tmp.imag();
    }
    const uint64_t *const vals = (const uint64_t *)tmp;
    for (int i = 0; i < N; i++) {
        uint64_t val = (vals[i] & valmask0) | valmask1;  // mantissa on 53 bits
        uint16_t expo = (vals[i] >> 52) & expmask0;      // exponent 11 bits
        // 1023 -> 52th pos -> 0th pos
        // 1075 -> 52th pos -> 52th pos
        int16_t trans = expo - 1075;
        uint64_t val2 = trans > 0 ? (val << trans) : (val >> -trans);
        res[i] = (vals[i] >> 63) ? -val2 : val2;
    }

    fftw_free(in);
    fftw_free(out);
}

FFT_Processor_FFTW::~FFT_Processor_FFTW()
{
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_cleanup();
}

// FFT_Processor_FFTW is thread-safe
FFT_Processor_FFTW fftplvl1(TFHEpp::lvl1param::n);
FFT_Processor_FFTW fftplvl2(TFHEpp::lvl2param::n);