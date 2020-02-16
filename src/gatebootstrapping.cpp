#include <cloudkey.hpp>
#include <keyswitch.hpp>
#include <mulfft.hpp>
#include <params.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;

template <typename T = uint32_t, uint32_t N = DEF_N>
array<array<double,N>,2*N> XaittGen(){
    array<array<double,N>, 2*N> xaitt;
    for(int i = 0;i<2*N;i++){
        array<T,N> xai = {};
        if(i<N) xai[i] = 1;
        else xai[i-N] = -1;
        if constexpr (N==DEF_N) TwistIFFTlvl1(xaitt[i],xai);
        else TwistIFFTlvl2(xaitt[i],xai);
    }
    return xaitt;
}

static const array<PolynomialInFDlvl1,2*DEF_N> xaittlvl1 = XaittGen<uint32_t,DEF_N>();
static const array<PolynomialInFDlvl2,2*DEF_nbar> xaittlvl2 = XaittGen<uint64_t,DEF_nbar>();

template <typename T = uint32_t, uint32_t N = DEF_N>
inline void PolynomialMulByXai(array<T, N> &res, const array<T, N> &poly,
                               const T a)
{
    if (a == 0)
        return;
    else if (a < N) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + N];
        for (int i = a; i < N; i++) res[i] = poly[i - a];
    }
    else {
        const T aa = a - N;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + N];
        for (int i = aa; i < N; i++) res[i] = -poly[i - aa];
    }
}

void PolynomialMulByXailvl1(Polynomiallvl1 &res, const Polynomiallvl1 &poly,
                            const uint32_t a)
{
    PolynomialMulByXai<uint32_t, DEF_N>(res, poly, a);
}

void PolynomialMulByXailvl2(Polynomiallvl2 &res, const Polynomiallvl2 &poly,
                            const uint64_t a)
{
    PolynomialMulByXai<uint64_t, DEF_nbar>(res, poly, a);
}

template <typename T = uint32_t, uint32_t N = DEF_N>
inline void PolynomialMulByXaiMinusOne(array<T, N> &res,
                                       const array<T, N> &poly, const T a)
{
    if (a == 0)
        return;
    else if (a < N) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + N] - poly[i];
        for (int i = a; i < N; i++) res[i] = poly[i - a] - poly[i];
    }
    else {
        const T aa = a - N;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + N] - poly[i];
        for (int i = aa; i < N; i++) res[i] = -poly[i - aa] - poly[i];
    }
}

void PolynomialMulByXaiMinusOnelvl1(Polynomiallvl1 &res,
                                    const Polynomiallvl1 &poly,
                                    const uint32_t a)
{
    PolynomialMulByXaiMinusOne<uint32_t, DEF_N>(res, poly, a);
}

inline void PolynomialMulByXaiMinusOnelvl2(Polynomiallvl2 &res,
                                           const Polynomiallvl2 &poly,
                                           const uint64_t a)
{
    PolynomialMulByXaiMinusOne<uint64_t, DEF_nbar>(res, poly, a);
}

template <typename T = uint32_t, uint32_t N = DEF_N>
static void PolynomialMulByXaiSub(array<T, N> &res,
                                       const array<T, N> &poly, const array<T, N> &mpoly,const T a){
    if (a < N) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + N] - mpoly[i];
        for (int i = a; i < N; i++) res[i] = poly[i - a] - mpoly[i];
    }
    else {
        const T aa = a - N;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + N] - mpoly[i];
        for (int i = aa; i < N; i++) res[i] = -poly[i - aa] - mpoly[i];
    }
}

template <typename T = uint32_t, uint32_t N = DEF_N>
inline void RotatedTestVector(array<array<T, N>, 2> &testvector, uint32_t bara,
                              const T μ)
{
    testvector[0] = {};
    if (bara < N) {
        for (int i = 0; i < bara; i++) testvector[1][i] = -μ;
        for (int i = bara; i < N; i++) testvector[1][i] = μ;
    }
    else {
        const T baraa = bara - N;
        for (int i = 0; i < baraa; i++) testvector[1][i] = μ;
        for (int i = baraa; i < N; i++) testvector[1][i] = -μ;
    }
}

template <typename T = uint32_t, uint32_t N = DEF_N, uint32_t l = DEF_l>
void KeyBundle(array<array<array<double, N>, 2>, 2 * l>& kbfft, const array<array<array<array<T, N>, 2>, 2 * l>,2*DEF_Addends>&bk, const array<T,DEF_Addends> &bara){
    array<T,N> temp,mtemp;
    for(int i = 0;i<2*l;i++){
       for(int j = 0;j<2;j++){
           PolynomialMulByXaiSub<T,N>(temp,bk[0][i][j],bk[1][i][j],bara[1]);
           PolynomialMulByXaiSub<T,N>(mtemp,bk[2][i][j],bk[3][i][j],bara[1]);
           PolynomialMulByXaiSub<T,N>(mtemp,temp,mtemp,bara[0]);
           if constexpr (N==DEF_N){
               TwistIFFTlvl1(kbfft[i][j],mtemp);
           }else if constexpr (N==DEF_nbar){
               TwistIFFTlvl2(kbfft[i][j],mtemp);
           }
       } 
    }
}

template <typename T = uint32_t, uint32_t N = DEF_N, uint32_t l = DEF_l>
void KeyBundleFFT(array<array<array<double, N>, 2>, 2 * l>& kbfft, const array<array<array<array<double, N>, 2>, 2 * l>,2*DEF_Addends>&bkfft, const array<T,DEF_Addends> &bara){
    for(int i = 0;i<2*l;i++){
       for(int j = 0;j<2;j++){
           constexpr uint32_t indexmask = 2*N-1;
           kbfft[i][j] = bkfft[3][i][j];
            if constexpr (N==DEF_N){
                FMAInFD<N>(kbfft[i][j],bkfft[2][i][j],xaittlvl1[bara[1] & indexmask]);
                FMAInFD<N>(kbfft[i][j],bkfft[1][i][j],xaittlvl1[bara[0] & indexmask]);
                FMAInFD<N>(kbfft[i][j],bkfft[0][i][j],xaittlvl1[(bara[0]+bara[1]) & indexmask]);
           }else if constexpr(N==DEF_nbar){
                FMAInFD<N>(kbfft[i][j],bkfft[2][i][j],xaittlvl2[bara[1] & indexmask]);
                FMAInFD<N>(kbfft[i][j],bkfft[1][i][j],xaittlvl2[bara[0] & indexmask]);
                FMAInFD<N>(kbfft[i][j],bkfft[0][i][j],xaittlvl2[(bara[0]+bara[1]) & indexmask]);
           }
       } 
    }
}

void GateBootstrappingTLWE2TRLWEFFTlvl01(TRLWElvl1 &acc, const TLWElvl0 &tlwe,
                                         const GateKey &gk)
{
    TRGSWFFTlvl1 kbfft;
    const uint32_t bara = 2 * DEF_N - modSwitchFromTorus32<2 * DEF_N>(tlwe[DEF_n]);
    RotatedTestVector<uint32_t, DEF_N>(acc, bara, DEF_μ);
    for (int i = 0; i < DEF_n/DEF_Addends; i++) {
        array<uint32_t,DEF_Addends> bara = {modSwitchFromTorus32<2 * DEF_N>(tlwe[2*i]),modSwitchFromTorus32<2 * DEF_N>(tlwe[2*i+1])};
        KeyBundleFFT<uint32_t,DEF_N,DEF_l>(kbfft, gk.bkfftlvl01[i], bara);
        trgswfftExternalProductlvl1(acc, acc, kbfft);
    }
}

void GateBootstrappingTLWE2TLWEFFTlvl01(TLWElvl1 &res, const TLWElvl0 &tlwe,
                                        const GateKey &gk)
{
    TRLWElvl1 acc;
    GateBootstrappingTLWE2TRLWEFFTlvl01(acc, tlwe, gk);
    SampleExtractIndexlvl1(res, acc, 0);
}

inline void GateBootstrappingTLWE2TLWEFFTlvl02(TLWElvl2 &res,
                                               const TLWElvl0 &tlwe,
                                               const CircuitKey &ck,
                                               const uint64_t μs2)
{
    TRGSWFFTlvl2 kbfft;
    TRLWElvl2 acc;
    uint32_t bara =
        2 * DEF_nbar - modSwitchFromTorus64<2 * DEF_nbar>(tlwe[DEF_n]);
    RotatedTestVector<uint64_t, DEF_nbar>(acc, bara, μs2);
    for (int i = 0; i < DEF_n/DEF_Addends; i++) {
        array<uint64_t,DEF_Addends> bara = {modSwitchFromTorus64<2 * DEF_nbar>(tlwe[2*i]),modSwitchFromTorus64<2 * DEF_nbar>(tlwe[2*i+1])};
        KeyBundleFFT<uint64_t,DEF_nbar,DEF_lbar>(kbfft, ck.bkfftlvl02[i], bara);
        trgswfftExternalProductlvl2(acc, acc, kbfft);
    }
    SampleExtractIndexlvl2(res, acc, 0);
    res[DEF_nbar] += μs2;
}

void GateBootstrapping(TLWElvl0 &res, const TLWElvl0 &tlwe, const GateKey &gk)
{
    TLWElvl1 tlwelvl1;
    GateBootstrappingTLWE2TLWEFFTlvl01(tlwelvl1, tlwe, gk);
    IdentityKeySwitchlvl10(res, tlwelvl1, gk.ksk);
}

inline void CircuitBootstrapping(TRGSWlvl1 &trgsw, const TLWElvl0 &tlwe,
                                 const CircuitKey &ck)
{
    TLWElvl2 tlwelvl2;
    for (int i = 0; i < DEF_l; i++) {
        GateBootstrappingTLWE2TLWEFFTlvl02(
            tlwelvl2, tlwe, ck, 1UL << (64 - (i + 1) * DEF_Bgbit - 1));
        PrivKeySwitchlvl21(trgsw[i], tlwelvl2, 0, ck.privksk);
        PrivKeySwitchlvl21(trgsw[i + DEF_l], tlwelvl2, 1, ck.privksk);
    }
}

void CircuitBootstrappingFFT(TRGSWFFTlvl1 &trgswfft, const TLWElvl0 &tlwe,
                             const CircuitKey &ck)
{
    TRGSWlvl1 trgsw;
    CircuitBootstrapping(trgsw, tlwe, ck);
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++) TwistIFFTlvl1(trgswfft[i][j], trgsw[i][j]);
}

void CircuitBootstrappingFFTwithInv(TRGSWFFTlvl1 &trgswfft,
                                    TRGSWFFTlvl1 &invtrgswfft,
                                    const TLWElvl0 &tlwe, const CircuitKey &ck)
{
    TRGSWlvl1 trgsw;
    array<uint32_t, DEF_l> h;
    for (int i = 0; i < DEF_l; i++) h[i] = 1U << (32 - (i + 1) * DEF_Bgbit);
    CircuitBootstrapping(trgsw, tlwe, ck);
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++) TwistIFFTlvl1(trgswfft[i][j], trgsw[i][j]);
    for (int i = 0; i < DEF_l; i++) {
        for (int j = 0; j < DEF_N; j++) {
            trgsw[i][0][j] *= -1;
            trgsw[i][1][j] *= -1;
            trgsw[i + DEF_l][0][j] *= -1;
            trgsw[i + DEF_l][1][j] *= -1;
        }
        trgsw[i][0][0] += h[i];
        trgsw[i + DEF_l][1][0] += h[i];
    }
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++)
            TwistIFFTlvl1(invtrgswfft[i][j], trgsw[i][j]);
}
}  // namespace TFHEpp