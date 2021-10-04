#pragma once

#include <array>

#include "params.hpp"

namespace TFHEpp {

template <class P>
void IdentityKeySwitch(TLWE<typename P::targetP> &res,
                       const TLWE<typename P::domainP> &tlwe,
                       const KeySwitchingKey<P> &ksk);

template <class P>
void TLWE2TRLWEIKS(TRLWE<typename P::targetP> &res,
                   const TLWE<typename P::domainP> &tlwe,
                   const TLWE2TRLWEIKSKey<P> &iksk);

template <class P>
void EvalAuto(TRLWE<P> &res, const TRLWE<P> &trlwe, const int d,
              const TRGSWFFT<P> &autokey);

template <class P>
void AnnihilateKeySwitching(TRLWE<P> &res, const TRLWE<P> &trlwe,
                            const AnnihilateKey<P> &ahk);

template <class P, uint num_func>
void AnnihilatePrivateKeySwitching(
    std::array<TRLWE<P>, num_func> &res, const TRLWE<P> &trlwe,
    const AnnihilateKey<P> &ahk,
    const std::array<TRGSWFFT<P>, num_func> &privks)
{
    static_assert(num_func > 0, "num_func must be bigger than 0");
    res[num_func - 1] = trlwe;
    for (int i = 0; i < P::nbit - 1; i++) {
        TRLWE<P> evaledauto;
        EvalAuto<P>(evaledauto, res[num_func - 1], (1 << (P::nbit - i)) + 1,
                    ahk[i]);
        for (int j = 0; j < 2 * P::n; j++)
            res[num_func - 1][0][j] += evaledauto[0][j];
    }
    for (int i = 0; i < num_func; i++) {
        TRLWE<P> evaledauto;
        EvalAuto<P>(evaledauto, res[num_func - 1], (1 << (P::nbit - i)) + 1,
                    privks[i]);
        for (int j = 0; j < 2 * P::n; j++)
            res[i][0][j] += res[num_func - 1][0][j] + evaledauto[0][j];
    }
}

template <class P>
void PrivKeySwitch(TRLWE<typename P::targetP> &res,
                   const TLWE<typename P::domainP> &tlwe,
                   const PrivateKeySwitchingKey<P> &privksk);
}  // namespace TFHEpp