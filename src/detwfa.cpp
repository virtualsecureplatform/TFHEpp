#include "detwfa.hpp"


namespace TFHEpp {
template <class P>
void CMUXFFT(TRLWE<P> &res, const TRGSWFFT<P> &cs, const TRLWE<P> &c1,
             const TRLWE<P> &c0)
{
    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) res[k][i] = c1[k][i] - c0[k][i];
    trgswfftExternalProduct<P>(res, res, cs);
    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) res[k][i] += c0[k][i];
}
#define INST(P)                                                     \
    template void CMUXFFT<P>(TRLWE<P> & res, const TRGSWFFT<P> &cs, \
                             const TRLWE<P> &c1, const TRLWE<P> &c0)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class bkP>
void CMUXFFTwithPolynomialMulByXaiMinusOne(TRLWE<typename bkP::targetP> &acc, const BootstrappingKeyElementFFT<bkP> &cs,
                                           const int a)
{
    TRLWE<typename bkP::targetP> temp;
    if constexpr(bkP::domainP::key_value_diff == 1){
        for (int k = 0; k < bkP::targetP::k + 1; k++)
            PolynomialMulByXaiMinusOne<typename bkP::targetP>(temp[k], acc[k], a);
        trgswfftExternalProduct<typename bkP::targetP>(temp, temp, cs[0]);
        for (int k = 0; k < bkP::targetP::k + 1; k++)
            for (int i = 0; i < bkP::targetP::n; i++) acc[k][i] += temp[k][i];
    }else{
        TRGSWFFT<typename bkP::targetP> trgsw = {};
        int count = 0;
        PolynomialInFD<typename bkP::targetP> poly;
        for(int i = bkP::domainP::key_value_min; i <= bkP::domainP::key_value_max; i++){
            if(i!=0){
                for(int j = 0; j < (bkP::targetP::k+1)*bkP::targetP::l; j++){
                    for (int k = 0; k < bkP::targetP::k + 1; k++){
                        PolynomialMulByXaiMinusOneInFD<typename bkP::targetP>(poly, cs[count][j][k], a*i);
                        for(int l = 0; l < bkP::targetP::n; l++) trgsw[j][k][l] += poly[l];
                    }
                }
                count++;
            }
        }
        trgswfftExternalProduct<typename bkP::targetP>(temp, acc, trgsw);
        for (int k = 0; k < bkP::targetP::k + 1; k++)
            for (int i = 0; i < bkP::targetP::n; i++) acc[k][i] += temp[k][i];
    }
}
#define INST(bkP)                                             \
    template void CMUXFFTwithPolynomialMulByXaiMinusOne<bkP>( \
        TRLWE<typename bkP::targetP> & acc, const BootstrappingKeyElementFFT<bkP> &cs, const int a)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void CMUXNTTwithPolynomialMulByXaiMinusOne(TRLWE<P> &acc, const TRGSWNTT<P> &cs,
                                           const typename P::T a)
{
    TRLWE<P> temp;
    for (int k = 0; k < P::k + 1; k++)
        PolynomialMulByXaiMinusOne<P>(temp[k], acc[k], a);
    trgswnttExternalProduct<P>(temp, temp, cs);
    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) acc[k][i] += temp[k][i];
}
#define INST(P)                                             \
    template void CMUXNTTwithPolynomialMulByXaiMinusOne<P>( \
        TRLWE<P> & acc, const TRGSWNTT<P> &cs, const typename P::T a)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

}  // namespace TFHEpp
