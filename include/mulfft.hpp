#pragma once

#include<spqlios-fft.h>

#include<params.hpp>
#include<utils.hpp>

namespace TFHEpp{
    using namespace std;

    inline void TwistFFTlvl1(array<uint32_t,DEF_N> &res,const array<double,DEF_N> &a){
        fftplvl1.execute_direct_torus32(res.data(),a.data());
    }

    inline void TwistIFFTlvl1(array<double,DEF_N> &res, const array<uint32_t,DEF_N> &a){
        fftplvl1.execute_reverse_torus32(res.data(),a.data());
    }

    inline void PolyMullvl1(array<uint32_t,DEF_N> &res, const array<uint32_t,DEF_N> &a, const array<uint32_t,DEF_N> &b){
        array<double,DEF_N> ffta;
        TwistIFFTlvl1(ffta,a);
        array<double,DEF_N> fftb;
        TwistIFFTlvl1(fftb,b);
        MulInFD<DEF_N>(ffta,ffta,fftb);
        TwistFFTlvl1(res,ffta);
    }

    inline void TwistFFTlvl2(array<uint64_t,DEF_nbar> &res,const array<double,DEF_nbar> &a){
        fftplvl2.execute_direct_torus64(res.data(),a.data());
    }

    inline void TwistIFFTlvl2(array<double,DEF_nbar> &res, const array<uint64_t,DEF_nbar> &a){
        fftplvl2.execute_reverse_torus64(res.data(),a.data());
    }

    inline void PolyMullvl2(array<uint64_t,DEF_nbar> &res, const array<uint64_t,DEF_nbar> &a, const array<uint64_t,DEF_nbar> &b){
        array<double,DEF_nbar> ffta;
        TwistIFFTlvl2(ffta,a);
        array<double,DEF_nbar> fftb;
        TwistIFFTlvl2(fftb,b);
        array<double,DEF_nbar> fftres;
        for (int32_t i = 0;i < DEF_nbar/2;i++){
            fftres[i] = ffta[i] * fftb[i] - ffta[i+DEF_nbar/2] * fftb[i+DEF_nbar/2];
            fftres[i+DEF_nbar/2] = ffta[i] * fftb[i+DEF_nbar/2] + ffta[i+DEF_nbar/2] * fftb[i];
        }
        TwistFFTlvl2(res,fftres);
    }
}