#pragma once

#include <array>

#include "params.hpp"

namespace TFHEpp {
using namespace std;

template <class P>
void IdentityKeySwitch(TLWE<typename P::targetP> &res,
                       const TLWE<typename P::domainP> &tlwe,
                       const KeySwitchingKey<P> &ksk);

template <class P>
void TLWE2TRLWEIKS(TRLWE<typename P::targetP> &res,
                       const TLWE<typename P::domainP> &tlwe,
                       const TLWE2TRLWEIKSKey<P> &iksk);

template <class P>
void AnnihilateKeySwitching(TRLWE<P> &res, const TRLWE<P> &trlwe, const AnnihilateKey<P> &ahk);

template <class P>
void PrivKeySwitch(TRLWE<typename P::targetP> &res,
                   const TLWE<typename P::domainP> &tlwe,
                   const PrivKeySwitchKey<P> &privksk);
}  // namespace TFHEpp