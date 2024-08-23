#include "fft_processor_spqliox_aarch64.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <params.hpp>

#include "spqliox_aarch64_impl.h"

using namespace std;

double accurate_cos(int32_t i, int32_t n)  // cos(2pi*i/n)
{
    constexpr double pi = M_PI;
    i = ((i % n) + n) % n;
    if (i >= 3 * n / 4)
        return std::cos(2. * pi * (n - i) / static_cast<double>(n));
    if (i >= 2 * n / 4)
        return -std::cos(2. * pi * (i - n / 2) / static_cast<double>(n));
    if (i >= 1 * n / 4)
        return -std::cos(2. * pi * (n / 2 - i) / static_cast<double>(n));
    return std::cos(2. * pi * (i) / static_cast<double>(n));
}

double accurate_sin(int32_t i, int32_t n)  // sin(2pi*i/n)
{
    constexpr double pi = M_PI;
    i = ((i % n) + n) % n;
    if (i >= 3 * n / 4)
        return -std::sin(2. * pi * (n - i) / static_cast<double>(n));
    if (i >= 2 * n / 4)
        return -std::sin(2. * pi * (i - n / 2) / static_cast<double>(n));
    if (i >= 1 * n / 4)
        return std::sin(2. * pi * (n / 2 - i) / static_cast<double>(n));
    return std::sin(2. * pi * (i) / static_cast<double>(n));
}

double *direct_trig_tables(int N)
{
    int32_t n = 2 * N, ns4 = n / 4;
    double *head = new (std::align_val_t(32)) double[n];
    double *ptr = head;

    // subsequent iterations
    for (int32_t halfnn = 4; halfnn < ns4; halfnn *= 2) {
        int32_t nn = 2 * halfnn;
        int32_t j = n / nn;
        // cerr << "- b: " << halfnn  << "(offset: " <<
        // (ptr-reps->trig_tables) << ", mult: " << j << ")" << endl;
        for (int32_t i = 0; i < halfnn; i += 4) {
            // cerr << "--- i: " << i << endl;
            for (int32_t k = 0; k < 4; k++)
                *(ptr++) = accurate_cos(-j * (i + k), n);
            for (int32_t k = 0; k < 4; k++)
                *(ptr++) = accurate_sin(-j * (i + k), n);
        }
    }
    // last iteration
    for (int32_t i = 0; i < ns4; i += 4) {
        for (int32_t k = 0; k < 4; k++) *(ptr++) = accurate_cos(-(i + k), n);
        for (int32_t k = 0; k < 4; k++) *(ptr++) = accurate_sin(-(i + k), n);
    }

    return head;
}

double *reverse_trig_tables(int N)
{
    int32_t n = 2 * N, ns4 = n / 4;
    double *head = new (std::align_val_t(32)) double[n];
    double *ptr = head;

    // first iteration
    for (int32_t j = 0; j < ns4; j += 4) {
        for (int32_t k = 0; k < 4; k++) *(ptr++) = accurate_cos(j + k, n);
        for (int32_t k = 0; k < 4; k++) *(ptr++) = accurate_sin(j + k, n);
    }
    // subsequent iterations
    for (int32_t nn = ns4; nn >= 8; nn /= 2) {
        int32_t halfnn = nn / 2;
        int32_t j = n / nn;
        // cerr << "- b: " << nn  << "(offset: " <<
        // (ptr-reps->trig_tables) << ", mult: " << j << ")" << endl;
        for (int32_t i = 0; i < halfnn; i += 4) {
            // cerr << "--- i: " << i << endl;
            for (int32_t k = 0; k < 4; k++)
                *(ptr++) = accurate_cos(j * (i + k), n);
            for (int32_t k = 0; k < 4; k++)
                *(ptr++) = accurate_sin(j * (i + k), n);
        }
    }

    return head;
}

double *negation_table_forward()
{
    double *head = new (std::align_val_t(32)) double[16];
    // for v24.d2
    head[0] = 1.0;
    head[1] = -1.0;
    // for v25.d2
    head[2] = 1.0;
    head[3] = -1.0;
    // for v26.d2
    head[4] = 1.0;
    head[5] = 1.0;
    // for v27.d2
    head[6] = -1.0;
    head[7] = -1.0;
    // for v28.d2
    head[8] = 1.0;
    head[9] = -1.0;
    // for v29.d2
    head[10] = -1.0;
    head[11] = 1.0;
    // for v30.d2
    head[12] = 1.0;
    head[13] = 1.0;
    // for v31.d2
    head[14] = 1.0;
    head[15] = -1.0;

    return head;
}

double *negation_table_reverse()
{
    double *head = new (std::align_val_t(32)) double[16];
    // for v24.d2
    head[0] = 1.0;
    head[1] = 1.0;
    // for v25.d2
    head[2] = 1.0;
    head[3] = -1.0;
    // for v26.d2
    head[4] = 1.0;
    head[5] = 1.0;
    // for v27.d2
    head[6] = -1.0;
    head[7] = 1.0;
    // for v28.d2
    head[8] = 1.0;
    head[9] = 1.0;
    // for v29.d2
    head[10] = -1.0;
    head[11] = -1.0;
    // for v30.d2
    head[12] = 1.0;
    head[13] = -1.0;
    // for v31.d2
    head[14] = 1.0;
    head[15] = -1.0;

    return head;
}

FFT_Processor_Spqliox_AArch64::FFT_Processor_Spqliox_AArch64(const int32_t N)
    : _2N(2 * N), N(N), Ns2(N / 2)
{
    assert(N >= 16);
    assert((N & (N - 1)) == 0);

    tables_direct_.n = 2 * N;
    tables_direct_.trig_tables = direct_trig_tables(N);
    tables_reverse_.n = 2 * N;
    tables_reverse_.trig_tables = reverse_trig_tables(N);
    table_negation_forward_.n = 16;
    table_negation_forward_.trig_tables = negation_table_forward();
    table_negation_reverse_.n = 16;
    table_negation_reverse_.trig_tables = negation_table_reverse();

    real_inout = new (std::align_val_t(32)) double[N];

    fft_code_.ready();
    fft_ = fft_code_.getCode<unsigned long (*)(
        double *, const double *, const double *, const double *,
        const FFT_PRECOMP *, double *, double *)>();
    ifft_code_.ready();
    ifft_ = ifft_code_.getCode<unsigned long (*)(
        double *, const double *, const double *, const FFT_PRECOMP *, double *,
        double *)>();
}

void FFT_Processor_Spqliox_AArch64::execute_reverse_int(double *res,
                                                        const int32_t *a)
{
    for (size_t i = 0; i < N; i++) real_inout[i] = (double)a[i];

    ifft_(real_inout, NULL, NULL, &tables_reverse_, real_inout,
          table_negation_reverse_.trig_tables);

    for (size_t i = 0; i < N; i++) res[i] = real_inout[i];
}

void FFT_Processor_Spqliox_AArch64::execute_reverse_uint(double *res,
                                                         const uint32_t *a)
{
    for (size_t i = 0; i < N; i++) real_inout[i] = (double)a[i];

    ifft_(real_inout, NULL, NULL, &tables_reverse_, real_inout,
          table_negation_reverse_.trig_tables);

    for (size_t i = 0; i < N; i++) res[i] = real_inout[i];
}

void FFT_Processor_Spqliox_AArch64::execute_reverse_torus32(double *res,
                                                            const uint32_t *a)
{
    int32_t *aa = (int32_t *)a;
    execute_reverse_int(res, aa);
}

void FFT_Processor_Spqliox_AArch64::execute_reverse_torus64(double *res,
                                                            const uint64_t *a)
{
    int64_t *aa = (int64_t *)a;
    for (int i = 0; i < N; i++) real_inout[i] = (double)aa[i];
    ifft_(real_inout, NULL, NULL, &tables_reverse_, real_inout,
          table_negation_reverse_.trig_tables);
    for (int i = 0; i < N; i++) res[i] = real_inout[i];
}

void FFT_Processor_Spqliox_AArch64::execute_direct_torus32(uint32_t *res,
                                                           const double *a)
{
    const double *sit = a;
    const double *send = a + N;
    static const double _2sN = double(2) / double(N);
    const double *bla = &_2sN;
    double *dst = real_inout;

    fft_(dst, sit, send, bla, &tables_direct_, real_inout,
         table_negation_forward_.trig_tables);

    for (int32_t i = 0; i < N; i++) res[i] = uint32_t(int64_t(real_inout[i]));
}

void FFT_Processor_Spqliox_AArch64::execute_direct_torus32_rescale(
    uint32_t *res, const double *a, const double Δ)
{
    const double *sit = a;
    const double *send = a + N;
    static const double _2sN = double(2) / double(N);
    const double *bla = &_2sN;
    double *dst = real_inout;
    fft_(dst, sit, send, bla, &tables_direct_, real_inout,
         table_negation_forward_.trig_tables);
    for (int32_t i = 0; i < N; i++)
        res[i] = uint32_t(int64_t(real_inout[i] / Δ));
}

void FFT_Processor_Spqliox_AArch64::execute_direct_torus64_rescale(
    uint64_t *res, const double *a, const double Δ)
{
    const double *sit = a;
    const double *send = a + N;
    static const double _2sN = double(2) / double(N);
    const double *bla = &_2sN;
    double *dst = real_inout;
    fft_(dst, sit, send, bla, &tables_direct_, real_inout,
         table_negation_forward_.trig_tables);
    for (int32_t i = 0; i < N; i++)
        res[i] = uint64_t(int64_t(real_inout[i] / Δ));
}

void FFT_Processor_Spqliox_AArch64::execute_direct_torus64(uint64_t *res,
                                                           const double *a)
{
    const double *sit = a;
    const double *send = a + N;
    static const double _2sN = double(2) / double(N);
    const double *bla = &_2sN;
    double *dst = real_inout;
    fft_(dst, sit, send, bla, &tables_direct_, real_inout,
         table_negation_forward_.trig_tables);
    const uint64_t *const vals = (const uint64_t *)real_inout;
    static const uint64_t valmask0 = 0x000FFFFFFFFFFFFFul;
    static const uint64_t valmask1 = 0x0010000000000000ul;
    static const uint16_t expmask0 = 0x07FFu;
    for (int i = 0; i < N; i++) {
        uint64_t val = (vals[i] & valmask0) | valmask1;  // mantissa on 53 bits
        uint16_t expo = (vals[i] >> 52) & expmask0;      // exponent 11 bits
        // 1023 -> 52th pos -> 0th pos
        // 1075 -> 52th pos -> 52th pos
        int16_t trans = expo - 1075;
        uint64_t val2 = trans > 0 ? (val << trans) : (val >> -trans);
        res[i] = (vals[i] >> 63) ? -val2 : val2;
    }
}

FFT_Processor_Spqliox_AArch64::~FFT_Processor_Spqliox_AArch64() {}

thread_local FFT_Processor_Spqliox_AArch64 fftplvl1(TFHEpp::lvl1param::n);
thread_local FFT_Processor_Spqliox_AArch64 fftplvl2(TFHEpp::lvl2param::n);
