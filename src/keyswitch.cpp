#include <cstdint>
#include <keyswitch.hpp>
#include <limits>
#include <trgsw.hpp>

#include "utils.hpp"

namespace TFHEpp {
using namespace std;

template <class P>
void IdentityKeySwitch(TLWE<typename P::targetP> &res,
                       const TLWE<typename P::domainP> &tlwe,
                       const KeySwitchingKey<P> &ksk)
{
    constexpr typename P::domainP::T prec_offset =
        1ULL << (numeric_limits<typename P::domainP::T>::digits -
                 (1 + P::basebit * P::t));
    constexpr uint32_t mask = (1U << P::basebit) - 1;
    res = {};
    constexpr uint32_t domain_digit =
        std::numeric_limits<typename P::domainP::T>::digits;
    constexpr uint32_t target_digit =
        std::numeric_limits<typename P::targetP::T>::digits;
    if constexpr (domain_digit == target_digit)
        res[P::targetP::k * P::targetP::n] =
            tlwe[P::domainP::k * P::domainP::n];
    else if constexpr (domain_digit > target_digit)
        res[P::targetP::k * P::targetP::n] =
            (tlwe[P::domainP::k * P::domainP::n] +
             (1ULL << (domain_digit - target_digit - 1))) >>
            (domain_digit - target_digit);
    else if constexpr (domain_digit < target_digit)
        res[P::targetP::k * P::targetP::n] = tlwe[P::domainP::k * P::domainP::n]
                                             << (target_digit - domain_digit);
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        const typename P::domainP::T aibar = tlwe[i] + prec_offset;
        for (int j = 0; j < P::t; j++) {
            const uint32_t aij =
                (aibar >> (numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;
            if (aij != 0)
                for (int k = 0; k <= P::targetP::k * P::targetP::n; k++)
                    res[k] -= ksk[i][j][aij - 1][k];
        }
    }
}
#define INST(P)                                                               \
    template void IdentityKeySwitch<P>(TLWE<typename P::targetP> & res,       \
                                       const TLWE<typename P::domainP> &tlwe, \
                                       const KeySwitchingKey<P> &ksk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(INST)
#undef INST

template <class P>
void TLWE2TRLWEIKS(TRLWE<typename P::targetP> &res,
                   const TLWE<typename P::domainP> &tlwe,
                   const TLWE2TRLWEIKSKey<P> &iksk)
{
    constexpr typename P::domainP::T prec_offset =
        1ULL << (numeric_limits<typename P::domainP::T>::digits -
                 (1 + P::basebit * P::t));
    constexpr uint32_t mask = (1U << P::basebit) - 1;
    res = {};
    constexpr uint32_t domain_digit =
        std::numeric_limits<typename P::domainP::T>::digits;
    constexpr uint32_t target_digit =
        std::numeric_limits<typename P::targetP::T>::digits;
    if constexpr (domain_digit == target_digit)
        res[1][0] = tlwe[P::domainP::n];
    else if constexpr (domain_digit > target_digit)
        res[1][0] = (tlwe[P::domainP::n] +
                     (1ULL << (domain_digit - target_digit - 1))) >>
                    (domain_digit - target_digit);
    else if constexpr (domain_digit < target_digit)
        res[1][0] = tlwe[P::domainP::n] << (target_digit - domain_digit);
    for (int i = 0; i < P::domainP::n; i++) {
        const typename P::domainP::T aibar = tlwe[i] + prec_offset;
        for (int j = 0; j < P::t; j++) {
            const uint32_t aij =
                (aibar >> (numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;
            if (aij != 0)
                for (int k = 0; k < P::targetP::n; k++) {
                    res[0][k] -= iksk[i][j][aij - 1][0][k];
                    res[1][k] -= iksk[i][j][aij - 1][1][k];
                }
        }
    }
}
#define INST(P)                                                           \
    template void TLWE2TRLWEIKS<P>(TRLWE<typename P::targetP> & res,      \
                                   const TLWE<typename P::domainP> &tlwe, \
                                   const TLWE2TRLWEIKSKey<P> &iksk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

template <class P>
void EvalAuto(TRLWE<P> &res, const TRLWE<P> &trlwe, const int d,
              const TRGSWFFT<P> &autokey)
{
    Polynomial<P> polyb;
    Automorphism<P>(polyb, trlwe[1], d);
    res = {};
    Automorphism<P>(res[1], trlwe[0], d);
    trgswfftExternalProduct<P>(res, res, autokey);
    for (int i = 0; i < P::n; i++) {
        res[0][i] = -res[0][i];
        res[1][i] = polyb[i] - res[1][i];
    }
}
#define INST(P)                                                      \
    template void EvalAuto<P>(TRLWE<P> & res, const TRLWE<P> &trlwe, \
                              const int d, const TRGSWFFT<P> &autokey)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
void AnnihilateKeySwitching(TRLWE<P> &res, const TRLWE<P> &trlwe,
                            const AnnihilateKey<P> &ahk)
{
    res = trlwe;
    for (int i = 0; i < P::nbit; i++) {
        TRLWE<P> evaledauto;
        EvalAuto<P>(evaledauto, res, (1 << (P::nbit - i)) + 1, ahk[i]);
        for (int j = 0; j < 2 * P::n; j++) res[0][j] += evaledauto[0][j];
    }
}
#define INST(P)                              \
    template void AnnihilateKeySwitching<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const AnnihilateKey<P> &ahk)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
void PrivKeySwitch(TRLWE<typename P::targetP> &res,
                   const TLWE<typename P::domainP> &tlwe,
                   const PrivateKeySwitchingKey<P> &privksk)
{
    constexpr uint32_t mask = (1 << P::basebit) - 1;
    constexpr uint64_t prec_offset =
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                 (1 + P::basebit * P::t));

    res = {};
    for (int i = 0; i <= P::domainP::k * P::domainP::n; i++) {
        const typename P::domainP::T aibar = tlwe[i] + prec_offset;

        for (int j = 0; j < P::t; j++) {
            const typename P::domainP::T aij =
                (aibar >> (std::numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;

            if (aij != 0) {
                for (int k = 0; k < P::targetP::k + 1; k++)
                    for (int p = 0; p < P::targetP::n; p++)
                        res[k][p] -= privksk[i][j][aij - 1][k][p];
            }
        }
    }
}
#define INST(P)                                                           \
    template void PrivKeySwitch<P>(TRLWE<typename P::targetP> & res,      \
                                   const TLWE<typename P::domainP> &tlwe, \
                                   const PrivateKeySwitchingKey<P> &privksk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(INST)
#undef INST

}  // namespace TFHEpp
