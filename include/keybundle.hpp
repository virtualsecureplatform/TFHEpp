# pragma once 

#include "params.hpp"
#include "mulfft.hpp"

namespace TFHEpp{

// template <typename T = uint32_t, uint32_t N = DEF_N>
// static void PolynomialMulByXaiSub(array<T, N> &res,
//                                        const array<T, N> &poly, const array<T, N> &mpoly,const T a){
//     if (a < N) {
//         for (int i = 0; i < a; i++) res[i] = -poly[i - a + N] - mpoly[i];
//         for (int i = a; i < N; i++) res[i] = poly[i - a] - mpoly[i];
//     }
//     else {
//         const T aa = a - N;
//         for (int i = 0; i < aa; i++) res[i] = poly[i - aa + N] - mpoly[i];
//         for (int i = aa; i < N; i++) res[i] = -poly[i - aa] - mpoly[i];
//     }
// }

// template <class P>
// void KeyBundle(TRGSWFFT<P>& kbfft, const array<array<array<array<T, N>, 2>, 2 * l>,2*DEF_Addends>&bk, const array<T,DEF_Addends> &bara){
//     array<T,N> temp,mtemp;
//     for(int i = 0;i<2*l;i++){
//        for(int j = 0;j<2;j++){
//            PolynomialMulByXaiSub<T,N>(temp,bk[0][i][j],bk[1][i][j],bara[1]);
//            PolynomialMulByXaiSub<T,N>(mtemp,bk[2][i][j],bk[3][i][j],bara[1]);
//            PolynomialMulByXaiSub<T,N>(mtemp,temp,mtemp,bara[0]);
//            if constexpr (std::is_same_v<P, lvl1param>){
//                TwistIFFTlvl1(kbfft[i][j],mtemp);
//            }else{
//                TwistIFFTlvl2(kbfft[i][j],mtemp);
//            }
//        } 
//     }
// }

template <class P>
void KeyBundleFFT(TRGSWFFT<typename P::targetP>& kbfft, const BootstrappingKeyElementFFT<P> &bkfft, const std::array<typename P::domainP::T,P::Addends> &bara){
    for(int i = 0;i<2*P::targetP::l;i++){
       for(int j = 0;j<P::targetP::k+1;j++){
           constexpr uint32_t indexmask = 2*P::targetP::n-1;
           kbfft[i][j] = bkfft[3][i][j];
            if constexpr (std::is_same_v<typename P::targetP, lvl1param>){
                FMAInFD<P::targetP::n>(kbfft[i][j],bkfft[2][i][j],xaittlvl1[bara[1] & indexmask]);
                FMAInFD<P::targetP::n>(kbfft[i][j],bkfft[1][i][j],xaittlvl1[bara[0] & indexmask]);
                FMAInFD<P::targetP::n>(kbfft[i][j],bkfft[0][i][j],xaittlvl1[(bara[0]+bara[1]) & indexmask]);
           }else{
                FMAInFD<P::targetP::n>(kbfft[i][j],bkfft[2][i][j],xaittlvl2[bara[1] & indexmask]);
                FMAInFD<P::targetP::n>(kbfft[i][j],bkfft[1][i][j],xaittlvl2[bara[0] & indexmask]);
                FMAInFD<P::targetP::n>(kbfft[i][j],bkfft[0][i][j],xaittlvl2[(bara[0]+bara[1]) & indexmask]);
           }
       } 
    }
}
}