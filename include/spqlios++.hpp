#pragma once

#include<utility>

#include "params.hpp"   
#include "utils.hpp"

namespace SPQLIOSpp{
    using namespace std;
    using namespace TFHEpp;

    template<uint32_t N>
    array<double,N> TwistGen(){
        array<double,N> twist;
        for(uint32_t i = 0;i<N/2;i++){
            twist[i] = cos(i*M_PI/N);
            twist[i+N/2] = sin(i*M_PI/N);
        }
        return twist;
    }

    template<uint32_t N>
    array<double,2*N> TableGen(){
        array<double,2*N> table;
        for(uint32_t i = 0;i<N;i++){
            table[i] = cos(2*i*M_PI/N);
            table[i+N] = sin(2*i*M_PI/N);
        }
        return table;
    }

    static const array<double,DEF_N> twistlvl1 = TwistGen<DEF_N>();
    static const array<double,DEF_N> tablelvl1 = TableGen<DEF_N/2>();

    template<uint32_t size, uint32_t N>
    inline void ButterflyAdd(double* const a){
        double* const b = a+size/2;
        double* const are = a;
        double* const aim = a+N;
        double* const bre = b;
        double* const bim = b+N;
        
        for (int i = 0;i<size/2;i++){
            const double tempre = are[i];
            are[i] += bre[i];
            bre[i] = tempre - bre[i];
            const double tempim = aim[i];
            aim[i] += bim[i];
            bim[i] = tempim - bim[i];
        }
    }

    template <uint32_t Nbit = DEF_Nbit, uint32_t step, uint32_t size, uint32_t stride, bool isinvert =true>
    inline void TwiddleMul(double* const a, const array<double,1<<(Nbit+1)> &table){
        constexpr uint32_t N = 1<<Nbit;

        double* const are = a;
        double* const aim = a+N;
        for(int i = 1; i<size; i++){
            const double bre = table[stride*(1<<step)*i];
            const double bim = isinvert?table[stride*(1<<step)*i + N]:-table[stride*(1<<step)*i + N];

            const double aimbim = aim[i] * bim;
            const double arebim = are[i] * bim;
            are[i] = are[i] * bre - aimbim;
            aim[i] = aim[i] * bre + arebim;
        }
    }
    template <uint32_t Nbit = DEF_Nbit-1, uint32_t size, bool isinvert =true>
    inline void Radix4TwiddleMul(double* const a){
        constexpr uint32_t N = 1<<Nbit;

        double* const are = a;
        double* const aim = a+N;
        for(int i = 0;i<size;i++){
            swap(are[i],aim[i]);
            if constexpr(isinvert){
                are[i]*=-1;
            }
            else{
                aim[i]*=-1;
            }
        }
    }

    template<uint32_t Nbit = DEF_Nbit-1, int step = 0>
    void IFFT(double* const res, const array<double,1<<(Nbit+1)> &table){
        constexpr uint32_t N = 1<<Nbit;
        constexpr uint32_t size = 1<<(Nbit-step);
        
        if constexpr(size == 2){
            ButterflyAdd<size,N>(res);
        }
        else if constexpr(size==4){
            ButterflyAdd<size,N>(res);
            Radix4TwiddleMul<Nbit,size/4,true>(res+size*3/4);
            ButterflyAdd<size/2,N>(res);
            ButterflyAdd<size/2,N>(res+size/2);
        }
        else{
            ButterflyAdd<size,N>(res);

            Radix4TwiddleMul<Nbit,size/4,true>(res+size*3/4);

            ButterflyAdd<size/2,N>(res+size/2);

            TwiddleMul<Nbit,step,size/4,1,true>(res+size/2,table);
            TwiddleMul<Nbit,step,size/4,3,true>(res+size*3/4,table);

            IFFT<Nbit,step+1>(res, table);
            IFFT<Nbit,step+2>(res+size/2, table);
            IFFT<Nbit,step+2>(res+size*3/4, table);
        }
    }

    void TwistIFFTlvl1(array<double ,DEF_N> &res, const array<uint32_t,DEF_N> &a){
        for (int i = 0; i < DEF_N / 2; i++) {
            const double are = static_cast<double>(static_cast<int32_t>(a[i]));
            const double aim = static_cast<double>(static_cast<int32_t>(a[i+DEF_N/2]));
            const double aimbim = aim * twistlvl1[i + DEF_N / 2];
            const double arebim = are * twistlvl1[i + DEF_N / 2];
            res[i] = are * twistlvl1[i] - aimbim;
            res[i + DEF_N / 2] = aim * twistlvl1[i] + arebim;
        }
        IFFT<DEF_Nbit-1,0>(res.data(),tablelvl1);
    }

    template<uint32_t Nbit = DEF_Nbit-1, int step = 0>
    void FFT(double* const res, const array<double,1<<(Nbit+1)> &table){
        constexpr uint32_t N = 1<<Nbit;
        constexpr uint32_t size = 1<<(Nbit-step);

        if constexpr(size == 2){
            ButterflyAdd<size,N>(res);
        }
        else if constexpr(size == 4){
            ButterflyAdd<size/2,N>(res);

            ButterflyAdd<size/2,N>(res+size/2);
            Radix4TwiddleMul<Nbit,size/4,false>(res+size*3/4);
            ButterflyAdd<size,N>(res);
        }
        else{
            FFT<Nbit,step+1>(res, table);
            FFT<Nbit,step+2>(res+size/2, table);
            FFT<Nbit,step+2>(res+size*3/4, table);

            TwiddleMul<Nbit,step,size/4,1,false>(res+size/2,table);
            TwiddleMul<Nbit,step,size/4,3,false>(res+size*3/4,table);
            ButterflyAdd<size/2,N>(res+size/2);
            Radix4TwiddleMul<Nbit,size/4,false>(res+size*3/4);
            ButterflyAdd<size,N>(res);
        }
    }

    void TwistFFTlvl1(array<uint32_t,DEF_N> &res, array<double,DEF_N> &a){
        FFT<DEF_Nbit-1,0>(a.data(),tablelvl1);
        for (int i = 0; i < DEF_N / 2; i++) {
            const double aimbim = a[i + DEF_N / 2] * -twistlvl1[i + DEF_N / 2];
            const double arebim = a[i] * -twistlvl1[i + DEF_N / 2];
            res[i] = static_cast<int64_t>((a[i] * twistlvl1[i] - aimbim)*(2.0/DEF_N));
            res[i + DEF_N / 2] = static_cast<int64_t>((a[i + DEF_N / 2] * twistlvl1[i] + arebim)*(2.0/DEF_N));
        }
    }

    inline void PolyMullvl1(Polynomiallvl1 &res, const Polynomiallvl1 &a,
                        const Polynomiallvl1 &b)
    {
        PolynomialInFDlvl1 ffta;
        TwistIFFTlvl1(ffta, a);
        PolynomialInFDlvl1 fftb;
        TwistIFFTlvl1(fftb, b);
        MulInFD<DEF_N>(ffta, ffta, fftb);
        TwistFFTlvl1(res, ffta);
    }
}