#include <cstdint>
#include <complex>
#include <array>
#include <memory>

#include<ipp.h>

#include <params.hpp>
#include <utils.hpp>

using namespace std;

namespace TFHEpp{

    template <uint32_t N = DEF_N>
    array<Ipp64f, N> TwistGen(){
        array<Ipp64f, N> twist;
        for(int32_t i = 0; i<N/2; i++){
            twist[i] = cos(i * M_PI / N);
            twist[i+DEF_N/2] = sin(i * M_PI / N);
        }
        return twist;
    }

    static const array<Ipp64f,DEF_N> twistlvl1 = TwistGen<DEF_N>();
    static const array<Ipp64f,DEF_nbar> twistlvl2 = TwistGen<DEF_nbar>();
    static IppsFFTSpec_C_64f* speclvl1 = 0;
    static Ipp8u* bufflvl1 = 0;
    static Ipp8u* pMemSpec = 0;
 
    void FFTInit(){
        int32_t sizeSpec = 0;
        int32_t sizeInit = 0;
        int32_t sizeBuffer = 0;
        ippsFFTGetSize_C_64f(DEF_Nbit, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuffer);
        Ipp8u* pMemInit = 0;
        if (sizeInit > 0){
            pMemInit = (Ipp8u*)ippMalloc(sizeInit);
        }

        if (sizeBuffer > 0){
            bufflvl1 = (Ipp8u*) ippMalloc(sizeBuffer);
        }
        ippsFFTInit_C_64f(&speclvl1, DEF_Nbit, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, pMemSpec, pMemInit);

        /// free initialization buffer
        if (sizeInit > 0){
            ippFree(pMemInit);
        }
    }

    template <uint32_t N = DEF_N>
    array<double, N> TwistFFT(array<double,N> a, array<Ipp64f,N> twist){
        array<Ipp64f,N> twisted;
        for(int32_t i = 0; i < N/2; i++){
            twisted[i] = a[i] * twist[i] - a[i + N/2] * twist[i + N/2];
            twisted[i + N/2] = a[i] * twist[i + N/2] + a[i + N/2] * twist[i];
        }
        array<Ipp64f,N> res;
        ippsFFTFwd_CToC_64f(&twisted[0],&twisted[N/2], &res[0], &res[N/2], speclvl1, bufflvl1);
        return res;
    }

    template <uint32_t N = DEF_N>
    array<double, N> TwistIFFT(array<double,N> a, array<Ipp64f,N> twist){
        array<Ipp64f,N> dst;
        ippsFFTInv_CToC_64f(&dst[0],&dst[N/2], &a[0], &a[N/2], speclvl1, bufflvl1);
        array<Ipp64f,N> untwisted;
        for(int32_t i = 0; i < N/2; i++){
            untwisted[i] = dst[i] * twist[i] + dst[i + N/2] * twist[i + N/2];
            untwisted[i + N/2] = dst[i] * -twist[i + N/2] + dst[i + N/2] * twist[i];
        }
        return untwisted;
    }

    template <typename T = uint32_t, uint32_t N = DEF_N>
    inline array<T,N> PolyMul(array<T,N> a, array<T,N>b, const array<Ipp64f,N> twist){
        array<Ipp64f,N> ffta = TwistFFT<N>(ttod32_array<N>(a),twist); 
        array<Ipp64f,N> fftb = TwistFFT<N>(ttod32_array<N>(b),twist); 
        array<Ipp64f,N> fftres;
        for (int32_t i = 0;i < N/2;i++){
            fftres[i] = ffta[i] * fftb[i] - ffta[i+N/2] * fftb[i+N/2];
            fftres[i+N/2] = ffta[i] * fftb[i+N/2] + ffta[i+N/2] * fftb[i];
        }
        array<double,N> doubleres = TwistIFFT<N>(fftres,twist);
        array<T,N> res;
        for(int32_t i = 0; i<N; i++){
            res[i] = static_cast<T>(fmod(round(doubleres[i]),numeric_limits<T>::max()/2));
        }
        return res;
    }
    
    array<uint32_t,DEF_N> PolyMullvl1(array<uint32_t,DEF_N> a, array<uint32_t,DEF_N> b){
        return PolyMul<uint32_t,DEF_N>(a,b,twistlvl1);
    }
}