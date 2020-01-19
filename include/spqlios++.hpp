#pragma once

#include<utility>
#include <type_traits>

#include "params.hpp"   
#include "utils.hpp"

namespace SPQLIOSpp{
    using namespace std;
    using namespace TFHEpp;

    template<uint32_t N>
    inline array<double,N> TwistGen(){
        array<double,N> twist;
        for(uint32_t i = 0;i<N/2;i++){
            twist[i] = cos(i*M_PI/N);
            twist[i+N/2] = sin(i*M_PI/N);
        }
        return twist;
    }

    template<uint32_t N>
    inline array<double,2*N> TableGen(){
        array<double,2*N> table;
        for(uint32_t i = 0;i<N;i++){
            table[i] = cos(2*i*M_PI/N);
            table[i+N] = sin(2*i*M_PI/N);
        }
        return table;
    }

    static const array<double,DEF_N> twistlvl1 = TwistGen<DEF_N>();
    static const array<double,DEF_N> tablelvl1 = TableGen<DEF_N/2>();
    static const array<double,DEF_nbar> twistlvl2 = TwistGen<DEF_nbar>();
    static const array<double,DEF_nbar> tablelvl2 = TableGen<DEF_nbar/2>();

    template<typename T = uint32_t,uint32_t N = DEF_N>
    inline void TwistMulInvert(array<double,N> &res, const array<T,N> &a, const array<double,N> &twist){
        for (int i = 0; i < N / 2; i++) {
            const double are = static_cast<double>(static_cast<typename make_signed<T>::type>(a[i]));
            const double aim = static_cast<double>(static_cast<typename make_signed<T>::type>(a[i+N/2]));
            const double aimbim = aim * twist[i + N / 2];
            const double arebim = are * twist[i + N / 2];
            res[i] = are * twist[i] - aimbim;
            res[i + N / 2] = aim * twist[i] + arebim;
        }
    }

    template<uint32_t N = DEF_N>
    inline void TwistMulDirectlvl1(array<uint32_t,N> &res, const array<double,N> &a, const array<double,N> &twist){
        for (int i = 0; i < N / 2; i++) {
            const double aimbim = a[i + N / 2] * -twist[i + N / 2];
            const double arebim = a[i] * -twist[i + N / 2];
            res[i] = static_cast<int64_t>((a[i] * twist[i] - aimbim)*(2.0/N));
            res[i + N / 2] = static_cast<int64_t>((a[i + N / 2] * twist[i] + arebim)*(2.0/N));
        }
    }

    template<uint32_t N = DEF_nbar>
    inline void TwistMulDirectlvl2(array<uint64_t,N> &res, const array<double,N> &a, const array<double,N> &twist){
        constexpr uint64_t valmask0 = 0x000FFFFFFFFFFFFFul;
        constexpr uint64_t valmask1 = 0x0010000000000000ul;
        constexpr uint16_t expmask0 = 0x07FFu;
        for (int i = 0; i < N / 2; i++) {
            const double aimbim = a[i + N / 2] * -twist[i + N / 2];
            const double arebim = a[i] * -twist[i + N / 2];
            const double resdoublere = (a[i] * twist[i] - aimbim)*(2.0/N);
            const double resdoubleim = (a[i + N / 2] * twist[i] + arebim)*(2.0/N);
            const uint64_t resre = reinterpret_cast<const uint64_t&>(resdoublere);
            const uint64_t resim = reinterpret_cast<const uint64_t&>(resdoubleim);

            uint64_t val = (resre&valmask0)|valmask1; //mantissa on 53 bits
            uint16_t expo = (resre>>52)&expmask0; //exponent 11 bits
            // 1023 -> 52th pos -> 0th pos
            // 1075 -> 52th pos -> 52th pos
            int16_t trans = expo-1075;
            uint64_t val2 = trans>0?(val<<trans):(val>>-trans);
            res[i] = (resre>>63)?-val2:val2;

            val = (resim&valmask0)|valmask1; //mantissa on 53 bits
            expo = (resim>>52)&expmask0; //exponent 11 bits
            // 1023 -> 52th pos -> 0th pos
            // 1075 -> 52th pos -> 52th pos
            trans = expo-1075;
            val2 = trans>0?(val<<trans):(val>>-trans);
            res[i + N / 2] = (resim>>63)?-val2:val2;
            
        }
    }

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

    template <uint32_t Nbit = DEF_Nbit-1, uint32_t size, bool isinvert =true>
    inline void Radix8TwiddleMulStrideOne(double* const a){
        constexpr uint32_t N = 1<<Nbit;

        double* const are = a;
        double* const aim = a+N;
        const double _1sroot2 = 1/sqrt(2);
        for(int i = 0;i<size;i++){
            const double aimbim = isinvert?aim[i]:-aim[i];
            const double arebim = isinvert?are[i]:-are[i];
            are[i] = _1sroot2*(are[i] - aimbim);
            aim[i] = _1sroot2*(aim[i] + arebim);
        }
    }

    template <uint32_t Nbit = DEF_Nbit-1, uint32_t size, bool isinvert =true>
    inline void Radix8TwiddleMulStrideThree(double* const a){
        constexpr uint32_t N = 1<<Nbit;

        double* const are = a;
        double* const aim = a+N;
        const double _1sroot2 = 1/sqrt(2);
        for(int i = 0;i<size;i++){
            const double aimbim = isinvert?aim[i]:-aim[i];
            const double arebim = isinvert?are[i]:-are[i];
            are[i] = _1sroot2*(-are[i] - aimbim);
            aim[i] = _1sroot2*(-aim[i] + arebim);
        }
    }

    template<uint32_t Nbit = DEF_Nbit-1, int step = 0>
    inline void IFFT(double* const res, const array<double,1<<(Nbit+1)> &table){
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
        else if constexpr(size == 8){
            ButterflyAdd<size,N>(res);

            IFFT<Nbit,step+1>(res, table);

            Radix8TwiddleMulStrideOne<Nbit,size/8,true>(res+size*5/8);
            Radix4TwiddleMul<Nbit,size/8,true>(res+size*3/4);
            Radix8TwiddleMulStrideThree<Nbit,size/8,true>(res+size*7/8);

            ButterflyAdd<size/2,N>(res+size/2);
            
            Radix4TwiddleMul<Nbit,size/8,true>(res+size*7/8);

            ButterflyAdd<size/4,N>(res+size/2);
            ButterflyAdd<size/4,N>(res+size*3/4);

        }
        else{
            ButterflyAdd<size,N>(res);

            IFFT<Nbit,step+1>(res, table);

            Radix8TwiddleMulStrideOne<Nbit,size/8,true>(res+size*5/8);
            Radix4TwiddleMul<Nbit,size/8,true>(res+size*3/4);
            Radix8TwiddleMulStrideThree<Nbit,size/8,true>(res+size*7/8);

            ButterflyAdd<size/2,N>(res+size/2);

            Radix4TwiddleMul<Nbit,size/8,true>(res+size*7/8);

            ButterflyAdd<size/4,N>(res+size/2);
            ButterflyAdd<size/4,N>(res+size*3/4);

            TwiddleMul<Nbit,step,size/8,1,true>(res+size/2,table);
            TwiddleMul<Nbit,step,size/8,5,true>(res+size*5/8,table);
            TwiddleMul<Nbit,step,size/8,3,true>(res+size*3/4,table);
            TwiddleMul<Nbit,step,size/8,7,true>(res+size*7/8,table);

            IFFT<Nbit,step+3>(res+size/2, table);
            IFFT<Nbit,step+3>(res+size*5/8, table);
            IFFT<Nbit,step+3>(res+size*3/4, table);
            IFFT<Nbit,step+3>(res+size*7/8, table);
        }
    }

    inline void TwistIFFTlvl1(array<double ,DEF_N> &res, const array<uint32_t,DEF_N> &a){
        TwistMulInvert<uint32_t,DEF_N>(res,a,twistlvl1);
        IFFT<DEF_Nbit-1,0>(res.data(),tablelvl1);
    }

    inline void TwistIFFTlvl2(array<double ,DEF_nbar> &res, const array<uint64_t,DEF_nbar> &a){
        TwistMulInvert<uint64_t,DEF_nbar>(res,a,twistlvl2);
        IFFT<DEF_nbarbit-1,0>(res.data(),tablelvl2);
    }

    template<uint32_t Nbit = DEF_Nbit-1, int step = 0>
    inline void FFT(double* const res, const array<double,1<<(Nbit+1)> &table){
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
        else if constexpr(size==8){
            FFT<Nbit,step+1>(res, table);

            ButterflyAdd<size/4,N>(res+size/2);
            ButterflyAdd<size/4,N>(res+size*3/4);

            Radix4TwiddleMul<Nbit,size/8,false>(res+size*7/8);

            ButterflyAdd<size/2,N>(res+size/2);

            Radix8TwiddleMulStrideOne<Nbit,size/8,false>(res+size*5/8);
            Radix4TwiddleMul<Nbit,size/8,false>(res+size*3/4);
            Radix8TwiddleMulStrideThree<Nbit,size/8,false>(res+size*7/8);

            ButterflyAdd<size,N>(res);
        }
        else{
            FFT<Nbit,step+1>(res, table);
            FFT<Nbit,step+3>(res+size/2, table);
            FFT<Nbit,step+3>(res+size*5/8, table);
            FFT<Nbit,step+3>(res+size*3/4, table);
            FFT<Nbit,step+3>(res+size*7/8, table);
            
            TwiddleMul<Nbit,step,size/8,1,false>(res+size/2,table);
            TwiddleMul<Nbit,step,size/8,5,false>(res+size*5/8,table);
            TwiddleMul<Nbit,step,size/8,3,false>(res+size*3/4,table);
            TwiddleMul<Nbit,step,size/8,7,false>(res+size*7/8,table);

            ButterflyAdd<size/4,N>(res+size/2);
            ButterflyAdd<size/4,N>(res+size*3/4);

            Radix4TwiddleMul<Nbit,size/8,false>(res+size*7/8);

            ButterflyAdd<size/2,N>(res+size/2);

            Radix8TwiddleMulStrideOne<Nbit,size/8,false>(res+size*5/8);
            Radix4TwiddleMul<Nbit,size/8,false>(res+size*3/4);
            Radix8TwiddleMulStrideThree<Nbit,size/8,false>(res+size*7/8);

            ButterflyAdd<size,N>(res);
        }
    }

    inline void TwistFFTlvl1(array<uint32_t,DEF_N> &res, array<double,DEF_N> &a){
        FFT<DEF_Nbit-1,0>(a.data(),tablelvl1);
        TwistMulDirectlvl1<DEF_N>(res,a,twistlvl1);
    }

    inline void TwistFFTlvl2(array<uint64_t,DEF_nbar> &res, array<double,DEF_nbar> &a){
        FFT<DEF_nbarbit-1,0>(a.data(),tablelvl2);
        TwistMulDirectlvl2<DEF_nbar>(res,a,twistlvl2);
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

    inline void PolyMulNaievelvl2(Polynomiallvl2 &res, const Polynomiallvl2 &a,
                              const Polynomiallvl2 &b)
    {
        for (int i = 0; i < DEF_nbar; i++) {
            uint64_t ri = 0;
            for (int j = 0; j <= i; j++)
                ri += static_cast<int64_t>(a[j]) * b[i - j];
            for (int j = i + 1; j < DEF_nbar; j++)
                ri -= static_cast<int64_t>(a[j]) * b[DEF_nbar + i - j];
            res[i] = ri;
        }
    }
}