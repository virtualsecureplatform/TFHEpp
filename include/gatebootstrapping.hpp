#pragma once

#include <cmath>
#include <limits>

#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

#ifdef USE_KEY_BUNDLE
#include "keybundle.hpp"
#endif

namespace TFHEpp {

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyFFT<P> &bkfft,
                 const Polynomial<typename P::targetP> &testvector)
{
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    const uint32_t b̄ = 2 * P::targetP::n -
                       ((tlwe[P::domainP::k * P::domainP::n] >>
                         (std::numeric_limits<typename P::domainP::T>::digits -
                          1 - P::targetP::nbit + bitwidth))
                        << bitwidth);
    res = {};
    PolynomialMulByXai<typename P::targetP>(res[P::targetP::k], testvector, b̄);
    #ifdef USE_KEY_BUNDLE
    alignas(64) std::array<TRGSWFFT<typename P::targetP>,  P::domainP::k * P::domainP::n/P::Addends> BKadded;
    #pragma omp parallel for num_threads(8)
    for (int i = 0; i < P::domainP::k * P::domainP::n/P::Addends; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        std::array<typename P::domainP::T,P::Addends> bara;
        bara[0] = (tlwe[2*i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        bara[1] = (tlwe[2*i+1] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        KeyBundleFFT<P>(BKadded[i], bkfft[i], bara);
    }
    for (int i = 0; i < P::domainP::k * P::domainP::n/P::Addends; i++){
        trgswfftExternalProduct<typename P::targetP>(res, res, BKadded[i]);
    }
    #else
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        const uint32_t ā =
            (tlwe[i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        if (ā == 0) continue;
        // Do not use CMUXFFT to avoid unnecessary copy.
            {
            alignas(64) TRLWE<typename P::targetP> temp;
            int count = 0;
            std::array<TRLWEInFD<typename P::targetP>,P::targetP::k + 1> temptrlwefft = {};
            #pragma omp parallel num_threads(2) 
            {
            for(int val = P::domainP::key_value_min; val <= P::domainP::key_value_max; val++){
                if(val!=0){
                    #pragma omp single
                    {
                    const int mod = (ā*val) % (2*P::targetP::n);
                    const int index = mod>0?mod:mod+(2*P::targetP::n);
                    for (int k = 0; k < P::targetP::k + 1; k++)
                        PolynomialMulByXaiMinusOne<typename P::targetP>(temp[k], res[k], index);
                    temptrlwefft = {};
                    }

                    const TRGSWFFT<typename P::targetP>& trgswfft = bkfft[i][count];
                    #pragma omp for
                    for(int j = 0; j < P::targetP::k + 1; j++){
                        for (int k = 0; k < P::targetP::l; k++) {
                            DecomposedPolynomialInFD<typename P::targetP> decpolyfft;
                            __builtin_prefetch(trgswfft[k + j * P::targetP::l].data());
                            DecompositionPolynomialFFT<typename P::targetP>(decpolyfft, temp[j], k);
                            for (int m = 0; m < P::targetP::k + 1; m++)
                                FMAInFD<P::targetP::n>(temptrlwefft[j][m], decpolyfft, trgswfft[k + j * P::targetP::l][m]);
                        }
                    }
                    #pragma omp single
                    {
                    for(int j = 1; j < P::targetP::k + 1; j++) for(int i = 0; i < P::targetP::k + 1; i++)  for(int k = 0; k < P::targetP::n; k++) temptrlwefft[0][i][k] += temptrlwefft[j][i][k];
                    for (int k = 0; k < P::targetP::k + 1; k++) TwistFFT<typename P::targetP>(temp[k], temptrlwefft[0][k]);

                    for (int k = 0; k < P::targetP::k + 1; k++)
                        for (int j = 0; j < P::targetP::n; j++) res[k][j] += temp[k][j];
                    count++;
                    }
                }
            }
            }
        }
    }
    #endif
}

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyFFT<P> &bkfft,
                 const TRLWE<typename P::targetP> &testvector)
{
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    const uint32_t b̄ = 2 * P::targetP::n -
                       ((tlwe[P::domainP::k * P::domainP::n] >>
                         (std::numeric_limits<typename P::domainP::T>::digits -
                          1 - P::targetP::nbit + bitwidth))
                        << bitwidth);
    for (int k = 0; k < P::targetP::k + 1; k++)
        PolynomialMulByXai<typename P::targetP>(res[k], testvector[k], b̄);
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        const uint32_t ā =
            (tlwe[i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        if (ā == 0) continue;
        // Do not use CMUXFFT to avoid unnecessary copy.
        CMUXFFTwithPolynomialMulByXaiMinusOne<typename P::targetP>(res,
                                                                   bkfft[i], ā);
    }
}

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyNTT<P> &bkntt,
                 const Polynomial<typename P::targetP> &testvector)
{
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    const uint32_t b̄ = 2 * P::targetP::n -
                       ((tlwe[P::domainP::k * P::domainP::n] >>
                         (std::numeric_limits<typename P::domainP::T>::digits -
                          1 - P::targetP::nbit + bitwidth))
                        << bitwidth);
    res = {};
    PolynomialMulByXai<typename P::targetP>(res[P::targetP::k], testvector, b̄);
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        const uint32_t ā =
            (tlwe[i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        if (ā == 0) continue;
        // Do not use CMUXNTT to avoid unnecessary copy.
        CMUXNTTwithPolynomialMulByXaiMinusOne<typename P::targetP>(res,
                                                                   bkntt[i], ā);
    }
}

template <class P>
void GateBootstrappingTLWE2TLWEFFT(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector);

template <class P>
void GateBootstrappingTLWE2TLWENTT(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyNTT<P> &bkntt,
    const Polynomial<typename P::targetP> &testvector){
        TRLWE<typename P::targetP> acc;
        BlindRotate<P>(acc, tlwe, bkntt, testvector);
        SampleExtractIndex<typename P::targetP>(res, acc, 0);
    }

template <class P, uint32_t num_out>
void GateBootstrappingManyLUT(
    std::array<TLWE<typename P::targetP>, num_out> &res,
    const TLWE<typename P::domainP> &tlwe, const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    TRLWE<typename P::targetP> acc;
    BlindRotate<P, num_out>(acc, tlwe, bkfft, testvector);
    for (int i = 0; i < num_out; i++)
        SampleExtractIndex<typename P::targetP>(res[i], acc, i);
}

template <class P, typename P::T μ>
constexpr Polynomial<P> μpolygen()
{
    Polynomial<P> poly;
    for (typename P::T &p : poly) p = μ;
    return poly;
}

template <typename lvl1param::T μ = lvl1param::μ>
void GateBootstrapping(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe,
                       const EvalKey &ek)
{
    TLWE<lvl1param> tlwelvl1;
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(tlwelvl1, tlwe, *ek.bkfftlvl01,
                                              μpolygen<lvl1param, μ>());
    IdentityKeySwitch<lvl10param>(res, tlwelvl1, *ek.iksklvl10);
}

template <typename lvl1param::T μ = lvl1param::μ>
void GateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                       const EvalKey &ek)
{
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlwe, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01,
                                              μpolygen<lvl1param, μ>());
}

template <typename lvl1param::T μ = lvl1param::μ>
void GateBootstrappingNTT(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe,
                       const EvalKey &ek)
{
    TLWE<lvl1param> tlwelvl1;
    GateBootstrappingTLWE2TLWENTT<lvl01param>(tlwelvl1, tlwe, *ek.bknttlvl01,
                                              μpolygen<lvl1param, μ>());
    IdentityKeySwitch<lvl10param>(res, tlwelvl1, *ek.iksklvl10);
}

template <typename lvl1param::T μ = lvl1param::μ>
void GateBootstrappingNTT(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                       const EvalKey &ek)
{
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlwe, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWENTT<lvl01param>(res, tlwelvl0, *ek.bknttlvl01,
                                              μpolygen<lvl1param, μ>());
}
}  // namespace TFHEpp