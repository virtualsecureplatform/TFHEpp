#include <cstdint>
#include <complex>
#include <array>
#include <memory>

#include<iostream>

#include<spqlios-fft.h>

#include <params.hpp>
#include <utils.hpp>
#include<cassert>

using namespace std;

namespace TFHEpp{

    template <uint32_t N = DEF_N>
    array<uint32_t, N> TwistFFT(const array<double,N> &a, FFT_Processor_Spqlios &fftp){
        array<uint32_t,N> res;
        fftp.execute_direct_torus32(&res[0],&a[0]);
        return res;
    }

    template <uint32_t N = DEF_N>
    array<double, N> TwistIFFT(const array<uint32_t,N> &a, FFT_Processor_Spqlios &fftp){
        array<double,N> res;
        fftp.execute_reverse_torus32(&res[0],&a[0]);
        return res;
    }

    template <typename T = uint32_t, uint32_t N = DEF_N>
    inline array<T,N> PolyMul(const array<T,N> &a, const array<T,N> &b, FFT_Processor_Spqlios &fftp){
        array<double,N> ffta = TwistIFFT<N>(a,fftp);
        array<double,N> fftb = TwistIFFT<N>(b,fftp);
        array<double,N> fftres;
        for (int32_t i = 0;i < N/2;i++){
            fftres[i] = ffta[i] * fftb[i] - ffta[i+N/2] * fftb[i+N/2];
            fftres[i+N/2] = ffta[i] * fftb[i+N/2] + ffta[i+N/2] * fftb[i];
        }
        array<T,N> res = TwistFFT<N>(fftres,fftp);
        return res;
    }
    
    array<uint32_t,DEF_N> PolyMullvl1(const array<uint32_t,DEF_N> &a, const array<uint32_t,DEF_N> &b){
        array<uint32_t,DEF_N> res = PolyMul<uint32_t,DEF_N>(a,b,fftplvl1);
        return res;
    }
}