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

    static const array<double,lvl1param::n> twistlvl1 = TwistGen<lvl1param::n>();
    static const array<double,lvl1param::n> tablelvl1 = TableGen<lvl1param::n/2>();
    static const array<double,lvl2param::n> twistlvl2 = TwistGen<lvl2param::n>();
    static const array<double,lvl2param::n> tablelvl2 = TableGen<lvl2param::n/2>();

    template<typename T = uint32_t,uint32_t N>
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

    template<uint32_t N>
    inline void TwistMulDirectlvl1(array<uint32_t,N> &res, const array<double,N> &a, const array<double,N> &twist){
        for (int i = 0; i < N / 2; i++) {
            const double aimbim = a[i + N / 2] * -twist[i + N / 2];
            const double arebim = a[i] * -twist[i + N / 2];
            res[i] = static_cast<int64_t>((a[i] * twist[i] - aimbim)*(2.0/N));
            res[i + N / 2] = static_cast<int64_t>((a[i + N / 2] * twist[i] + arebim)*(2.0/N));
        }
    }

    template<uint32_t N>
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

    template<uint32_t N>
    inline void ButterflyAdd(double* const a, double* const b, int i){
        double* const are = a;
        double* const aim = a+N;
        double* const bre = b;
        double* const bim = b+N;
        
        const double tempre = are[i];
        are[i] += bre[i];
        bre[i] = tempre - bre[i];
        const double tempim = aim[i];
        aim[i] += bim[i];
        bim[i] = tempim - bim[i];
    }

    template <uint32_t Nbit, uint32_t step, uint32_t stride, bool isinvert =true>
    inline void TwiddleMul(double* const a, const array<double,1<<(Nbit+1)> &table,int i){
        constexpr uint32_t N = 1<<Nbit;

        double* const are = a;
        double* const aim = a+N;
        const double bre = table[stride*(1<<step)*i];
        const double bim = isinvert?table[stride*(1<<step)*i + N]:-table[stride*(1<<step)*i + N];

        const double aimbim = aim[i] * bim;
        const double arebim = are[i] * bim;
        are[i] = are[i] * bre - aimbim;
        aim[i] = aim[i] * bre + arebim;
    }

    template <uint32_t Nbit, bool isinvert =true>
    inline void Radix4TwiddleMul(double* const a, int i){
        constexpr uint32_t N = 1<<Nbit;

        double* const are = a;
        double* const aim = a+N;
        swap(are[i],aim[i]);
        if constexpr(isinvert){
            are[i]*=-1;
        }
        else{
            aim[i]*=-1;
        }
    }

    template <uint32_t Nbit, bool isinvert =true>
    inline void Radix8TwiddleMulStrideOne(double* const a,int i){
        constexpr uint32_t N = 1<<Nbit;

        double* const are = a;
        double* const aim = a+N;
        const double _1sroot2 = 1/sqrt(2);
        const double aimbim = isinvert?aim[i]:-aim[i];
        const double arebim = isinvert?are[i]:-are[i];
        are[i] = _1sroot2*(are[i] - aimbim);
        aim[i] = _1sroot2*(aim[i] + arebim);
    }

    template <uint32_t Nbit, bool isinvert =true>
    inline void Radix8TwiddleMulStrideThree(double* const a, int i){
        constexpr uint32_t N = 1<<Nbit;

        double* const are = a;
        double* const aim = a+N;
        const double _1sroot2 = 1/sqrt(2);
        const double aimbim = isinvert?aim[i]:-aim[i];
        const double arebim = isinvert?are[i]:-are[i];
        are[i] = _1sroot2*(-are[i] - aimbim);
        aim[i] = _1sroot2*(-aim[i] + arebim);
    }

    template<uint32_t Nbit, int step = 0>
    inline void IFFT(double* const res, const array<double,1<<(Nbit+1)> &table){
        constexpr uint32_t N = 1<<Nbit;
        constexpr uint32_t size = 1<<(Nbit-step);
        
        if constexpr(size == 2){
            double* const res0 = &res[0];
            double* const res1 = &res[size/2];
            ButterflyAdd<N>(res0,res1,0);
        }
        else if constexpr(size==4){
            double* const res0 = &res[0];
            double* const res1 = &res[size/4];
            double* const res2 = &res[size*2/4];
            double* const res3 = &res[size*3/4];
            ButterflyAdd<N>(res0,res2,0);
            ButterflyAdd<N>(res1,res3,0);
            Radix4TwiddleMul<Nbit,true>(res3,0);
            ButterflyAdd<N>(res0,res1,0);
            ButterflyAdd<N>(res2,res3,0);
        }
        else if constexpr(size == 8){
            double* const res0 = &res[0];
            double* const res1 = &res[size/8];
            double* const res2 = &res[size*2/8];
            double* const res3 = &res[size*3/8];
            double* const res4 = &res[size*4/8];
            double* const res5 = &res[size*5/8];
            double* const res6 = &res[size*6/8];
            double* const res7 = &res[size*7/8];

            ButterflyAdd<N>(res0,res4,0);
            ButterflyAdd<N>(res1,res5,0);
            ButterflyAdd<N>(res2,res6,0);
            ButterflyAdd<N>(res3,res7,0);

            Radix8TwiddleMulStrideOne<Nbit,true>(res5,0);
            Radix4TwiddleMul<Nbit,true>(res6,0);
            Radix8TwiddleMulStrideThree<Nbit,true>(res7,0);

            ButterflyAdd<N>(res0,res2,0);
            ButterflyAdd<N>(res1,res3,0);
            ButterflyAdd<N>(res4,res6,0);
            ButterflyAdd<N>(res5,res7,0);
            
            Radix4TwiddleMul<Nbit,true>(res3,0);
            Radix4TwiddleMul<Nbit,true>(res7,0);

            ButterflyAdd<N>(res0,res1,0);
            ButterflyAdd<N>(res2,res3,0);
            ButterflyAdd<N>(res4,res5,0);
            ButterflyAdd<N>(res6,res7,0);

        }
        else{
            double* const res0 = &res[0];
            double* const res1 = &res[size/8];
            double* const res2 = &res[size*2/8];
            double* const res3 = &res[size*3/8];
            double* const res4 = &res[size*4/8];
            double* const res5 = &res[size*5/8];
            double* const res6 = &res[size*6/8];
            double* const res7 = &res[size*7/8];
            for(int i = 0;i<size/8;i++){

                ButterflyAdd<N>(res0,res4,i);
                ButterflyAdd<N>(res1,res5,i);
                ButterflyAdd<N>(res2,res6,i);
                ButterflyAdd<N>(res3,res7,i);

                Radix8TwiddleMulStrideOne<Nbit,true>(res5,i);
                Radix4TwiddleMul<Nbit,true>(res6,i);
                Radix8TwiddleMulStrideThree<Nbit,true>(res7,i);

                ButterflyAdd<N>(res0,res2,i);
                ButterflyAdd<N>(res1,res3,i);
                ButterflyAdd<N>(res4,res6,i);
                ButterflyAdd<N>(res5,res7,i);
                
                Radix4TwiddleMul<Nbit,true>(res3,i);
                Radix4TwiddleMul<Nbit,true>(res7,i);

                ButterflyAdd<N>(res0,res1,i);
                ButterflyAdd<N>(res2,res3,i);
                ButterflyAdd<N>(res4,res5,i);
                ButterflyAdd<N>(res6,res7,i);

                if(i!=0){
                    TwiddleMul<Nbit,step,4,true>(res1,table,i);
                    TwiddleMul<Nbit,step,2,true>(res2,table,i);
                    TwiddleMul<Nbit,step,6,true>(res3,table,i);
                    TwiddleMul<Nbit,step,1,true>(res4,table,i);
                    TwiddleMul<Nbit,step,5,true>(res5,table,i);
                    TwiddleMul<Nbit,step,3,true>(res6,table,i);
                    TwiddleMul<Nbit,step,7,true>(res7,table,i);
                }
            }

            IFFT<Nbit,step+3>(res0, table);
            IFFT<Nbit,step+3>(res1, table);
            IFFT<Nbit,step+3>(res2, table);
            IFFT<Nbit,step+3>(res3, table);
            IFFT<Nbit,step+3>(res4, table);
            IFFT<Nbit,step+3>(res5, table);
            IFFT<Nbit,step+3>(res6, table);
            IFFT<Nbit,step+3>(res7, table);
        }
    }

    inline void TwistIFFTlvl1(array<double ,lvl1param::n> &res, const array<uint32_t,lvl1param::n> &a){
        TwistMulInvert<uint32_t,lvl1param::n>(res,a,twistlvl1);
        IFFT<lvl1param::nbit-1,0>(res.data(),tablelvl1);
    }

    inline void TwistIFFTlvl2(array<double ,lvl2param::n> &res, const array<uint64_t,lvl2param::n> &a){
        TwistMulInvert<uint64_t,lvl2param::n>(res,a,twistlvl2);
        IFFT<lvl2param::nbit-1,0>(res.data(),tablelvl2);
    }

    template<uint32_t Nbit, int step = 0>
    inline void FFT(double* const res, const array<double,1<<(Nbit+1)> &table){
        constexpr uint32_t N = 1<<Nbit;
        constexpr uint32_t size = 1<<(Nbit-step);

        if constexpr(size == 2){
            double* const res0 = &res[0];
            double* const res1 = &res[size/2];
            ButterflyAdd<N>(res0,res1,0);
        }
        else if constexpr(size == 4){
            double* const res0 = &res[0];
            double* const res1 = &res[size/4];
            double* const res2 = &res[size*2/4];
            double* const res3 = &res[size*3/4];

            ButterflyAdd<N>(res0,res1,0);

            ButterflyAdd<N>(res2,res3,0);

            Radix4TwiddleMul<Nbit,false>(res3,0);
            ButterflyAdd<N>(res0,res2,0);
            ButterflyAdd<N>(res1,res3,0);
        }
        else if constexpr(size==8){
            double* const res0 = &res[0];
            double* const res1 = &res[size/8];
            double* const res2 = &res[size*2/8];
            double* const res3 = &res[size*3/8];
            double* const res4 = &res[size*4/8];
            double* const res5 = &res[size*5/8];
            double* const res6 = &res[size*6/8];
            double* const res7 = &res[size*7/8];

            ButterflyAdd<N>(res0,res1,0);
            ButterflyAdd<N>(res2,res3,0);
            ButterflyAdd<N>(res4,res5,0);
            ButterflyAdd<N>(res6,res7,0);

            Radix4TwiddleMul<Nbit,false>(res3,0);
            Radix4TwiddleMul<Nbit,false>(res7,0);

            ButterflyAdd<N>(res0,res2,0);
            ButterflyAdd<N>(res1,res3,0);
            ButterflyAdd<N>(res4,res6,0);
            ButterflyAdd<N>(res5,res7,0);

            Radix8TwiddleMulStrideOne<Nbit,false>(res5,0);
            Radix4TwiddleMul<Nbit,false>(res6,0);
            Radix8TwiddleMulStrideThree<Nbit,false>(res7,0);

            ButterflyAdd<N>(res0,res4,0);
            ButterflyAdd<N>(res1,res5,0);
            ButterflyAdd<N>(res2,res6,0);
            ButterflyAdd<N>(res3,res7,0);
        }
        else{
            
            double* const res0 = &res[0];
            double* const res1 = &res[size/8];
            double* const res2 = &res[size*2/8];
            double* const res3 = &res[size*3/8];
            double* const res4 = &res[size*4/8];
            double* const res5 = &res[size*5/8];
            double* const res6 = &res[size*6/8];
            double* const res7 = &res[size*7/8];

            FFT<Nbit,step+3>(res0, table);
            FFT<Nbit,step+3>(res1, table);
            FFT<Nbit,step+3>(res2, table);
            FFT<Nbit,step+3>(res3, table);
            FFT<Nbit,step+3>(res4, table);
            FFT<Nbit,step+3>(res5, table);
            FFT<Nbit,step+3>(res6, table);
            FFT<Nbit,step+3>(res7, table);

            for(int i = 0;i<size/8;i++){
            
                if(i!=0){
                    TwiddleMul<Nbit,step,4,false>(res1,table,i);
                    TwiddleMul<Nbit,step,2,false>(res2,table,i);
                    TwiddleMul<Nbit,step,6,false>(res3,table,i);
                    TwiddleMul<Nbit,step,1,false>(res4,table,i);
                    TwiddleMul<Nbit,step,5,false>(res5,table,i);
                    TwiddleMul<Nbit,step,3,false>(res6,table,i);
                    TwiddleMul<Nbit,step,7,false>(res7,table,i);
                }

                ButterflyAdd<N>(res0,res1,i);
                ButterflyAdd<N>(res2,res3,i);
                ButterflyAdd<N>(res4,res5,i);
                ButterflyAdd<N>(res6,res7,i);

                Radix4TwiddleMul<Nbit,false>(res3,i);
                Radix4TwiddleMul<Nbit,false>(res7,i);

                ButterflyAdd<N>(res0,res2,i);
                ButterflyAdd<N>(res1,res3,i);
                ButterflyAdd<N>(res4,res6,i);
                ButterflyAdd<N>(res5,res7,i);

                Radix8TwiddleMulStrideOne<Nbit,false>(res5,i);
                Radix4TwiddleMul<Nbit,false>(res6,i);
                Radix8TwiddleMulStrideThree<Nbit,false>(res7,i);

                ButterflyAdd<N>(res0,res4,i);
                ButterflyAdd<N>(res1,res5,i);
                ButterflyAdd<N>(res2,res6,i);
                ButterflyAdd<N>(res3,res7,i);
            }
        }
    }

    inline void TwistFFTlvl1(array<uint32_t,lvl1param::n> &res, array<double,lvl1param::n> &a){
        FFT<lvl1param::nbit-1,0>(a.data(),tablelvl1);
        TwistMulDirectlvl1<lvl1param::n>(res,a,twistlvl1);
    }

    inline void TwistFFTlvl2(array<uint64_t,lvl2param::n> &res, array<double,lvl2param::n> &a){
        FFT<lvl2param::nbit-1,0>(a.data(),tablelvl2);
        TwistMulDirectlvl2<lvl2param::n>(res,a,twistlvl2);
    }

    inline void PolyMullvl1(Polynomial<lvl1param> &res, const Polynomial<lvl1param> &a,
                        const Polynomial<lvl1param> &b)
    {
        PolynomialInFD<lvl1param> ffta;
        TwistIFFTlvl1(ffta, a);
        PolynomialInFD<lvl1param> fftb;
        TwistIFFTlvl1(fftb, b);
        MulInFD<lvl1param::n>(ffta, ffta, fftb);
        TwistFFTlvl1(res, ffta);
    }
}