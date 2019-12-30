#include<params.hpp>
#include<trlwe.hpp>
#include<mulfft.hpp>
#include<utils.hpp>

#include<array>

namespace TFHEpp{

    inline array<array<array<uint32_t,DEF_N>,2>,2*DEF_l> trgswSymEncryptlvl1(int32_t p, double α, array<uint32_t,DEF_N> &key){
        array<uint32_t,DEF_l> h;
        for(int  i = 0;i<DEF_l;i++) h[i] = 1U<<(32-(i+1)*DEF_Bgbit);
        array<array<array<uint32_t,DEF_N>,2>,2*DEF_l> trgsw;
        for(array<array<uint32_t,DEF_N>,2> &trlwe:trgsw) trlwe = trlweSymEncryptZerolvl1(α,key);
        for(int i = 0; i<DEF_l; i++){
            trgsw[i][0][0] += static_cast<uint32_t>(p)*h[i];
            trgsw[i+DEF_l][1][0] += static_cast<uint32_t>(p)*h[i];
        }
        return trgsw;
    }

    array<array<array<double,DEF_N>,2>,2*DEF_l> trgswfftSymEncryptlvl1(int32_t p, double α, array<uint32_t,DEF_N> &key){
        array<array<array<uint32_t,DEF_N>,2>,2*DEF_l> trgsw = trgswSymEncryptlvl1(p,α,key);
        array<array<array<double,DEF_N>,2>,2*DEF_l> trgswfft;
        for(int i = 0;i<2*DEF_l;i++) for(int j = 0; j<2; j++) TwistIFFTlvl1(trgswfft[i][j],trgsw[i][j]);
        return trgswfft;
    }

    template<typename T = uint32_t, uint32_t N = DEF_N, uint32_t l = DEF_l, uint32_t Bgbit = DEF_Bgbit>
    inline void Decomposition(array<array<T,N>,2*l> &decvec, const array<array<T,N>,2> &trlwe, const T offset){
        const T mask = static_cast<T>((1<<Bgbit) - 1);
        array<array<T,N>,2> siftedtrlwe;
        for(int i = 0;i<N;i++) siftedtrlwe[0][i] = trlwe[0][i] + offset;
        for(int i = 0;i<N;i++) siftedtrlwe[1][i] = trlwe[1][i] + offset;

        for(int i = 0; i<l; i++){
            for(int j = 0; j<N; j++) decvec[i][j] = ((siftedtrlwe[0][j]>>(32-(i+1)*Bgbit)) & mask) - (1U<<(Bgbit-1));
            for(int j = 0; j<N; j++) decvec[i+l][j] = ((siftedtrlwe[1][j]>>(32-(i+1)*Bgbit)) & mask) - (1U<<(Bgbit-1));
        }
    }

    constexpr uint32_t offsetgenlvl1(){
        uint32_t offset = 0;
        for(int i = 1;i <= DEF_l; i++) offset += DEF_Bg/2 * (1U << (32-i*DEF_Bgbit));
        return offset;
    }

    inline void Decompositionlvl1(array<array<uint32_t,DEF_N>,2*DEF_l> &decvec, const array<array<uint32_t,DEF_N>,2> &trlwe){
        constexpr uint32_t offset = offsetgenlvl1();
        Decomposition<uint32_t,DEF_N,DEF_l,DEF_Bgbit>(decvec,trlwe,offset);
    }

    inline void DecompositionFFTlvl1(array<array<double,DEF_N>,2*DEF_l> &decvecfft, const array<array<uint32_t,DEF_N>,2> &trlwe){
        array<array<uint32_t,DEF_N>,2*DEF_l> decvec;
        Decompositionlvl1(decvec,trlwe);
        for(int i =0;i<2*DEF_l;i++) TwistIFFTlvl1(decvecfft[i],decvec[i]);
    }

    void trgswfftExternalProductlvl1(array<array<uint32_t,DEF_N>,2> &restrlwe, const array<array<array<double,DEF_N>,2>,2*DEF_l> trgswfft, const array<array<uint32_t,DEF_N>,2> &trlwe){
        array<array<double,DEF_N>,2*DEF_l> decvecfft;
        DecompositionFFTlvl1(decvecfft,trlwe);
        array<array<double,DEF_N>,2> restrlwefft;
        MulInFD<DEF_N>(restrlwefft[0],decvecfft[0],trgswfft[0][0]);
        MulInFD<DEF_N>(restrlwefft[1],decvecfft[0],trgswfft[0][1]);
        for(int i = 1; i < 2*DEF_l; i++){
            FMAInFD<DEF_N>(restrlwefft[0],decvecfft[i],trgswfft[i][0]);
            FMAInFD<DEF_N>(restrlwefft[1],decvecfft[i],trgswfft[i][1]);
        }
        TwistFFTlvl1(restrlwe[0],restrlwefft[0]);
        TwistFFTlvl1(restrlwe[1],restrlwefft[1]);
    }
}