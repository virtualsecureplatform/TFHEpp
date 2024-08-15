#include <cassert>
#include <cmath>
#include<cstdint>

#include<params.hpp>

#include "x86.h"
#include "fft_processor_spqlios.h"

using namespace std;

int32_t rev(int32_t x, int32_t M) {
    int32_t reps = 0;
    for (int32_t j = M; j > 1; j /= 2) {
        reps = 2 * reps + (x % 2);
        x /= 2;
    }
    return reps;
}

FFT_Processor_Spqlios::FFT_Processor_Spqlios(const int32_t N) : _2N(2 * N), N(N), Ns2(N / 2) {
    tables_direct = new_fft_table(N);
    tables_reverse = new_ifft_table(N);
    real_inout_direct = fft_table_get_buffer(tables_direct);
    imag_inout_direct = real_inout_direct + Ns2;
    reva = new int32_t[Ns2];
    cosomegaxminus1 = new double[2 * _2N];
    sinomegaxminus1 = cosomegaxminus1 + _2N;
    int32_t rev1 = rev(1, _2N);
    int32_t rev3 = rev(3, _2N);
    //printf("rev-interval: %d, %d\n",rev1,rev3);
    for (int32_t revi = rev1; revi < rev3; revi++)
        reva[revi - rev1] = rev(revi, _2N);
    for (int32_t j = 0; j < _2N; j++) {
        cosomegaxminus1[j] = cos(2 * M_PI * j / _2N) - 1.;
        sinomegaxminus1[j] = sin(2 * M_PI * j / _2N);
    }
}

void FFT_Processor_Spqlios::execute_reverse_uint(double *res, const uint32_t *a) {
    //for (int32_t i=0; i<N; i++) real_inout_rev[i]=(double)a[i];
    {
        double *dst = res;
        // double *dst = real_inout_rev;
        const uint32_t *ait = a;
        const uint32_t *aend = a + N;
        __asm__ __volatile__ (
        "0:\n"
                "vmovupd (%1),%%xmm0\n"
                "vcvtudq2pd %%xmm0,%%ymm1\n"
                "vmovapd %%ymm1,(%0)\n"
                "addq $16,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 0b\n"
        : "=r"(dst), "=r"(ait), "=r"(aend)
        : "0"(dst), "1"(ait), "2"(aend)
        : "%xmm0", "%ymm1", "memory"
        );
    }
    ifft(tables_reverse, res);
}

void FFT_Processor_Spqlios::execute_reverse_int(double *res, const int32_t *a) {
    //for (int32_t i=0; i<N; i++) real_inout_rev[i]=(double)a[i];
    {
        double *dst = res;
        // double *dst = real_inout_rev;
        const int32_t *ait = a;
        const int32_t *aend = a + N;
        __asm__ __volatile__ (
        "0:\n"
                "vmovupd (%1),%%xmm0\n"
                "vcvtdq2pd %%xmm0,%%ymm1\n"
                "vmovapd %%ymm1,(%0)\n"
                "addq $16,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 0b\n"
        : "=r"(dst), "=r"(ait), "=r"(aend)
        : "0"(dst), "1"(ait), "2"(aend)
        : "%xmm0", "%ymm1", "memory"
        );
    }
    ifft(tables_reverse, res);
}

void FFT_Processor_Spqlios::execute_reverse_torus32(double *res, const uint32_t *a) {
    int32_t *aa = (int32_t *) a;
    execute_reverse_int(res, aa);
}

void FFT_Processor_Spqlios::execute_reverse_torus64(double* res, const uint64_t* a) {
    #ifdef USE_AVX512
    __m512d * ri512 = (__m512d *) res;
    __m512i * aa = (__m512i *) a;
    for (size_t i = 0; i < proc->N/8; i++) ri512[i] = _mm512_cvtepi64_pd (aa[i]);
    #else
    int64_t *aa = (int64_t *)a;
    for (int i=0; i<N; i++) res[i]=(double)aa[i];
    #endif
    ifft(tables_reverse,res);
}

void FFT_Processor_Spqlios::execute_direct_torus32(uint32_t *res, const double *a) {
    //TODO: parallelization
    static const double _2sN = double(2) / double(N);
    //for (int32_t i=0; i<N; i++) real_inout_direct[i]=a[i]*_2sn;
    {
        double *dst = real_inout_direct;
        const double *sit = a;
        const double *send = a + N;
        //double __2sN = 2./N;
        const double *bla = &_2sN;
        __asm__ __volatile__ (
        "vbroadcastsd (%3),%%ymm2\n"
                "1:\n"
                "vmovupd (%1),%%ymm0\n"
                "vmulpd	%%ymm2,%%ymm0,%%ymm0\n"
                "vmovapd %%ymm0,(%0)\n"
                "addq $32,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 1b\n"
        : "=r"(dst), "=r"(sit), "=r"(send), "=r"(bla)
        : "0"(dst), "1"(sit), "2"(send), "3"(bla)
        : "%ymm0", "%ymm2", "memory"
        );
    }
    fft(tables_direct, real_inout_direct);
    // for (int32_t i = 0; i < N; i++) res[i] = uint32_t(int64_t(real_inout_direct[i]));
    SPQLIOS::convert_f64_to_u32(res,real_inout_direct,N);
}

void FFT_Processor_Spqlios::execute_direct_torus32_q(uint32_t *res, const double *a, const uint32_t q) {
    //TODO: parallelization
    static const double _2sN = double(2) / double(N);
    //for (int32_t i=0; i<N; i++) real_inout_direct[i]=a[i]*_2sn;
    {
        double *dst = real_inout_direct;
        const double *sit = a;
        const double *send = a + N;
        //double __2sN = 2./N;
        const double *bla = &_2sN;
        __asm__ __volatile__ (
        "vbroadcastsd (%3),%%ymm2\n"
                "1:\n"
                "vmovupd (%1),%%ymm0\n"
                "vmulpd	%%ymm2,%%ymm0,%%ymm0\n"
                "vmovapd %%ymm0,(%0)\n"
                "addq $32,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 1b\n"
        : "=r"(dst), "=r"(sit), "=r"(send), "=r"(bla)
        : "0"(dst), "1"(sit), "2"(send), "3"(bla)
        : "%ymm0", "%ymm2", "memory"
        );
    }
    fft(tables_direct, real_inout_direct);
    for (int32_t i = 0; i < N; i++) res[i] = uint32_t((int64_t(real_inout_direct[i])%q+q)%q);
}

void FFT_Processor_Spqlios::execute_direct_torus32_rescale(uint32_t *res, const double *a, const double Δ) {
    //TODO: parallelization
    static const double _2sN = double(2) / double(N);
    //for (int32_t i=0; i<N; i++) real_inout_direct[i]=a[i]*_2sn;
    {
        double *dst = real_inout_direct;
        const double *sit = a;
        const double *send = a + N;
        //double __2sN = 2./N;
        const double *bla = &_2sN;
        __asm__ __volatile__ (
        "vbroadcastsd (%3),%%ymm2\n"
                "1:\n"
                "vmovupd (%1),%%ymm0\n"
                "vmulpd	%%ymm2,%%ymm0,%%ymm0\n"
                "vmovapd %%ymm0,(%0)\n"
                "addq $32,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 1b\n"
        : "=r"(dst), "=r"(sit), "=r"(send), "=r"(bla)
        : "0"(dst), "1"(sit), "2"(send), "3"(bla)
        : "%ymm0", "%ymm2", "memory"
        );
    }
    fft(tables_direct, real_inout_direct);
    for (int32_t i = 0; i < N; i++) res[i] = static_cast<uint32_t>(int64_t(real_inout_direct[i]/Δ));
}

void FFT_Processor_Spqlios::execute_direct_torus64(uint64_t* res, const double* a) {
    static const double _2sN = double(2)/double(N);
    //static const double _2p64 = pow(2.,64);
    //for (int i=0; i<N; i++) real_inout_direct[i]=a[i]*_2sn;
    {
    	double* dst = real_inout_direct;
	const double* sit = a;
	const double* send = a+N;
	//double __2sN = 2./N;
	const double* bla = &_2sN;
	__asm__ __volatile__ (
		"vbroadcastsd (%3),%%ymm2\n"
		"1:\n"
		"vmovupd (%1),%%ymm0\n"
		"vmulpd	%%ymm2,%%ymm0,%%ymm0\n"
		"vmovapd %%ymm0,(%0)\n"
		"addq $32,%1\n"
		"addq $32,%0\n"
		"cmpq %2,%1\n"
		"jb 1b\n"
		: "=r"(dst),"=r"(sit),"=r"(send),"=r"(bla)
		: "0"(dst),"1"(sit),"2"(send),"3"(bla)
		: "%ymm0","%ymm2","memory"
		);
    }
    fft(tables_direct,real_inout_direct); 
    #ifdef USE_AVX512
    __m512d * ri512 = (__m512d *) proc->real_inout_direct;
    __m512i * res512 = (__m512i *) res;
    const __m512d modc = {64, 64, 64, 64, 64, 64, 64, 64};
    for (size_t i = 0; i < proc->N/8; i++) {
        const __m512d _1 = _mm512_scalef_pd (ri512[i], -modc);
        const __m512d _2 = _mm512_reduce_pd (_1, 0);
        const __m512d _3 = _mm512_scalef_pd (_2, modc);
        res512[i] = _mm512_cvtpd_epi64 (_3);
    }
    #else
    const uint64_t* const vals = (const uint64_t*) real_inout_direct;
    static const uint64_t valmask0 = 0x000FFFFFFFFFFFFFul;
    static const uint64_t valmask1 = 0x0010000000000000ul;
    static const uint16_t expmask0 = 0x07FFu;
    for (int i=0; i<N; i++) {
        uint64_t val = (vals[i]&valmask0)|valmask1; //mantissa on 53 bits
        uint16_t expo = (vals[i]>>52)&expmask0; //exponent 11 bits
        // 1023 -> 52th pos -> 0th pos
        // 1075 -> 52th pos -> 52th pos
        int16_t trans = expo-1075;
        uint64_t val2 = trans>0?(val<<trans):(val>>-trans);
        res[i]=(vals[i]>>63)?-val2:val2;
    }
    #endif
}

void FFT_Processor_Spqlios::execute_direct_torus64_rescale(uint64_t* res, const double* a, const double Δ) {
    static const double _2sN = double(2)/double(N);
    //static const double _2p64 = pow(2.,64);
    //for (int i=0; i<N; i++) real_inout_direct[i]=a[i]*_2sn;
    {
    	double* dst = real_inout_direct;
	const double* sit = a;
	const double* send = a+N;
	//double __2sN = 2./N;
	const double* bla = &_2sN;
	__asm__ __volatile__ (
		"vbroadcastsd (%3),%%ymm2\n"
		"1:\n"
		"vmovupd (%1),%%ymm0\n"
		"vmulpd	%%ymm2,%%ymm0,%%ymm0\n"
		"vmovapd %%ymm0,(%0)\n"
		"addq $32,%1\n"
		"addq $32,%0\n"
		"cmpq %2,%1\n"
		"jb 1b\n"
		: "=r"(dst),"=r"(sit),"=r"(send),"=r"(bla)
		: "0"(dst),"1"(sit),"2"(send),"3"(bla)
		: "%ymm0","%ymm2","memory"
		);
    }
    fft(tables_direct,real_inout_direct); 
    for (int i=0; i<N; i++) res[i] = uint64_t(std::round(real_inout_direct[i]/Δ));
}

FFT_Processor_Spqlios::~FFT_Processor_Spqlios() {
    delete (tables_direct);
    delete (tables_reverse);
    delete[] cosomegaxminus1;
}

thread_local FFT_Processor_Spqlios fftplvl1(TFHEpp::lvl1param::n);
thread_local FFT_Processor_Spqlios fftplvl2(TFHEpp::lvl2param::n);
