#include <cassert>
#include <cmath>
#include<cstdint>

#include<params.hpp>

#include "spqlios-fft.h"

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
    real_inout_rev = fft_table_get_buffer(tables_reverse);
    imag_inout_rev = real_inout_rev + Ns2;
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

void FFT_Processor_Spqlios::execute_reverse_int(double *res, const int32_t *a) {
    //for (int32_t i=0; i<N; i++) real_inout_rev[i]=(double)a[i];
    {
        double *dst = real_inout_rev;
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
    ifft(tables_reverse, real_inout_rev);
    //for (int32_t i=0; i<N; i++) res[i]=real_inout_rev[i];
    {
        double *dst = res;
        double *sit = real_inout_rev;
        double *send = real_inout_rev + N;
        __asm__ __volatile__ (
        "1:\n"
                "vmovapd (%1),%%ymm0\n"
                "vmovupd %%ymm0,(%0)\n"
                "addq $32,%1\n"
                "addq $32,%0\n"
                "cmpq %2,%1\n"
                "jb 1b\n"
                "vzeroall\n"
        : "=r"(dst), "=r"(sit), "=r"(send)
        : "0"(dst), "1"(sit), "2"(send)
        : "%ymm0", "memory"
        );
    }
}

void FFT_Processor_Spqlios::execute_reverse_torus32(double *res, const uint32_t *a) {
    int32_t *aa = (int32_t *) a;
    execute_reverse_int(res, aa);
}

void FFT_Processor_Spqlios::execute_reverse_torus64(double* res, const uint64_t* a) {
    int64_t *aa = (int64_t *)a;
    for (int i=0; i<N; i++) real_inout_rev[i]=(double)aa[i];
    ifft(tables_reverse,real_inout_rev);
    for (int i=0; i<N; i++) res[i]=real_inout_rev[i];
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
    for (int32_t i = 0; i < N; i++) res[i] = uint32_t(int64_t(real_inout_direct[i]));
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
}

FFT_Processor_Spqlios::~FFT_Processor_Spqlios() {
    //delete (tables_direct);
    //delete (tables_reverse);
    delete[] cosomegaxminus1;
}

thread_local FFT_Processor_Spqlios fftplvl1(TFHEpp::DEF_N);
thread_local FFT_Processor_Spqlios fftplvl2(TFHEpp::DEF_nbar);