#pragma once

#include <array>

#include "params.hpp"
#include "trgsw.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class P>
constexpr typename P::domainP::T iksoffsetgen()
{
    typename P::domainP::T offset = 0;
    for (int i = 1; i <= P::t; i++)
        offset +=
            (1ULL << P::basebit) / 2 *
            (1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                      i * P::basebit));
    return offset;
}

template <class P>
void IdentityKeySwitch(TLWE<typename P::targetP> &res,
                       const TLWE<typename P::domainP> &tlwe,
                       const KeySwitchingKey<P> &ksk)
{
    res = {};
    constexpr uint domain_digit =
        std::numeric_limits<typename P::domainP::T>::digits;
    constexpr uint target_digit =
        std::numeric_limits<typename P::targetP::T>::digits;
    constexpr typename P::domainP::T roundoffset =
        (P::basebit * P::t) < domain_digit
            ? 1ULL << (domain_digit - (1 + P::basebit * P::t))
            : 0;
    if constexpr (domain_digit == target_digit)
        res[P::targetP::k * P::targetP::n] =
            tlwe[P::domainP::k * P::domainP::n];
    else if constexpr (domain_digit > target_digit)
        res[P::targetP::k * P::targetP::n] =
            (tlwe[P::domainP::k * P::domainP::n] +
             (1ULL << (domain_digit - target_digit - 1))) >>
            (domain_digit - target_digit);
    else if constexpr (domain_digit < target_digit)
        res[P::targetP::k * P::targetP::n] =
            static_cast<typename P::targetP::T>(
                tlwe[P::domainP::k * P::domainP::n])
            << (target_digit - domain_digit);

    // Koga's Optimization
    constexpr typename P::domainP::T offset = iksoffsetgen<P>();
    constexpr typename P::domainP::T mask = (1ULL << P::basebit) - 1;
    constexpr typename P::domainP::T halfbase = 1ULL << (P::basebit - 1);

    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        const typename P::domainP::T aibar = tlwe[i] + offset + roundoffset;
        for (int j = 0; j < P::t; j++) {
            const int32_t aij =
                ((aibar >>
                  (std::numeric_limits<typename P::domainP::T>::digits -
                   (j + 1) * P::basebit)) &
                 mask) -
                halfbase;
            if (aij > 0)
                for (int k = 0; k <= P::targetP::k * P::targetP::n; k++)
                    res[k] -= ksk[i][j][aij - 1][k];
            else if (aij < 0)
                for (int k = 0; k <= P::targetP::k * P::targetP::n; k++)
                    res[k] += ksk[i][j][-aij - 1][k];
        }
    }
}

template <class P, uint numcat>
void CatIdentityKeySwitch(
    std::array<TLWE<typename P::targetP>, numcat> &res,
    const std::array<TLWE<typename P::domainP>, numcat> &tlwe,
    const KeySwitchingKey<P> &ksk)
{
    res = {};
    constexpr uint domain_digit =
        std::numeric_limits<typename P::domainP::T>::digits;
    constexpr uint target_digit =
        std::numeric_limits<typename P::targetP::T>::digits;
    constexpr typename P::domainP::T roundoffset =
        (P::basebit * P::t) < domain_digit
            ? 1ULL << (domain_digit - (1 + P::basebit * P::t))
            : 0;

    for (int cat = 0; cat < numcat; cat++) {
        if constexpr (domain_digit == target_digit)
            res[cat][P::targetP::k * P::targetP::n] =
                tlwe[cat][P::domainP::k * P::domainP::n];
        else if constexpr (domain_digit > target_digit)
            res[cat][P::targetP::k * P::targetP::n] =
                (tlwe[cat][P::domainP::k * P::domainP::n] +
                 (1ULL << (domain_digit - target_digit - 1))) >>
                (domain_digit - target_digit);
        else if constexpr (domain_digit < target_digit)
            res[cat][P::targetP::k * P::targetP::n] =
                static_cast<typename P::targetP::T>(
                    tlwe[cat][P::domainP::k * P::domainP::n])
                << (target_digit - domain_digit);
    }

    // Koga's Optimization
    constexpr typename P::domainP::T offset = iksoffsetgen<P>();
    constexpr typename P::domainP::T mask = (1ULL << P::basebit) - 1;
    constexpr typename P::domainP::T halfbase = 1ULL << (P::basebit - 1);
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        std::array<typename P::domainP::T, numcat> aibarcat;
        for (int cat = 0; cat < numcat; cat++)
            aibarcat[cat] = tlwe[cat][i] + offset + roundoffset;
        for (int j = 0; j < P::t; j++) {
            for (int cat = 0; cat < numcat; cat++) {
                const int32_t aij =
                    ((aibarcat[cat] >>
                      (std::numeric_limits<typename P::domainP::T>::digits -
                       (j + 1) * P::basebit)) &
                     mask) -
                    halfbase;
                if (aij > 0)
                    for (int k = 0; k <= P::targetP::k * P::targetP::n; k++)
                        res[cat][k] -= ksk[i][j][aij - 1][k];
                else if (aij < 0)
                    for (int k = 0; k <= P::targetP::k * P::targetP::n; k++)
                        res[cat][k] += ksk[i][j][-aij - 1][k];
            }
        }
    }
}

template <class P>
void SubsetIdentityKeySwitch(TLWE<typename P::targetP> &res,
                             const TLWE<typename P::domainP> &tlwe,
                             const SubsetKeySwitchingKey<P> &ksk)
{
    constexpr typename P::domainP::T prec_offset =
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                 (1 + P::basebit * P::t));
    constexpr uint32_t mask = (1U << P::basebit) - 1;
    res = {};
    constexpr uint32_t domain_digit =
        std::numeric_limits<typename P::domainP::T>::digits;
    constexpr uint32_t target_digit =
        std::numeric_limits<typename P::targetP::T>::digits;
    if constexpr (domain_digit == target_digit) {
        for (int i = 0; i < P::targetP::k * P::targetP::n; i++)
            res[i] = tlwe[i];
        res[P::targetP::k * P::targetP::n] =
            tlwe[P::domainP::k * P::domainP::n];
    }
    else if constexpr (domain_digit > target_digit) {
        for (int i = 0; i < P::targetP::k * P::targetP::n; i++)
            res[i] = (tlwe[i] + (1ULL << (domain_digit - target_digit - 1))) >>
                     (domain_digit - target_digit);
        res[P::targetP::k * P::targetP::n] =
            (tlwe[P::domainP::k * P::domainP::n] +
             (1ULL << (domain_digit - target_digit - 1))) >>
            (domain_digit - target_digit);
    }
    else if constexpr (domain_digit < target_digit) {
        for (int i = 0; i < P::targetP::k * P::targetP::n; i++)
            res[i] = tlwe[i] << (target_digit - domain_digit);
        res[P::targetP::k * P::targetP::n] = tlwe[P::domainP::k * P::domainP::n]
                                             << (target_digit - domain_digit);
    }
    for (int i = 0;
         i < P::domainP::k * P::domainP::n - P::targetP::k * P::targetP::n;
         i++) {
        const typename P::domainP::T aibar =
            tlwe[i + P::targetP::n] + prec_offset;
        for (int j = 0; j < P::t; j++) {
            const uint32_t aij =
                (aibar >> (std::numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;
            if (aij != 0)
                for (int k = 0; k <= P::targetP::k * P::targetP::n; k++)
                    res[k] -= ksk[i][j][aij - 1][k];
        }
    }
}

template <class P>
void TLWE2TRLWEIKS(TRLWE<typename P::targetP> &res,
                   const TLWE<typename P::domainP> &tlwe,
                   const TLWE2TRLWEIKSKey<P> &iksk)
{
    constexpr typename P::domainP::T prec_offset =
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                 (1 + P::basebit * P::t));
    constexpr uint32_t mask = (1U << P::basebit) - 1;
    res = {};
    constexpr uint32_t domain_digit =
        std::numeric_limits<typename P::domainP::T>::digits;
    constexpr uint32_t target_digit =
        std::numeric_limits<typename P::targetP::T>::digits;
    if constexpr (domain_digit == target_digit)
        res[P::targetP::k][0] = tlwe[P::domainP::n];
    else if constexpr (domain_digit > target_digit)
        res[P::targetP::k][0] = (tlwe[P::domainP::n] +
                                 (1ULL << (domain_digit - target_digit - 1))) >>
                                (domain_digit - target_digit);
    else if constexpr (domain_digit < target_digit)
        res[P::targetP::k][0] = tlwe[P::domainP::n]
                                << (target_digit - domain_digit);
    for (int i = 0; i < P::domainP::n; i++) {
        const typename P::domainP::T aibar = tlwe[i] + prec_offset;
        for (int j = 0; j < P::t; j++) {
            const uint32_t aij =
                (aibar >> (std::numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;
            if (aij != 0)
                for (int l = 0; l < P::targetP::k + 1; l++)
                    for (int k = 0; k < P::targetP::n; k++)
                        res[l][k] -= iksk[i][j][aij - 1][l][k];
        }
    }
}

template <class P>
void EvalAuto(TRLWE<P> &res, const TRLWE<P> &trlwe, const int d,
              const EvalAutoKey<P> &autokey)
{
    res = {};
    Automorphism<P>(res[P::k], trlwe[P::k], d);

    for(int i = 0; i < P::k; i++){
        Polynomial<P> temppoly;
        TRLWE<P> temptrlwe;
        Automorphism<P>(temppoly, trlwe[i], d);
        halftrgswfftExternalProduct<P>(temptrlwe, temppoly, autokey[i]);
        for(int j = 0; j < P::k+1; j++)
            for (int k = 0; k < P::n; k++)
                res[j][k] -= temptrlwe[j][k];
    }
}

template <class P>
void AnnihilateKeySwitching(TRLWE<P> &res, const TRLWE<P> &trlwe,
                            const AnnihilateKey<P> &ahk)
{
    res = trlwe;
    for (int i = 0; i < P::nbit; i++) {
        TRLWE<P> evaledauto;
        for (int j = 0; j < (P::k+1) * P::n; j++) res[0][j] /= 2;
        EvalAuto<P>(evaledauto, res, (1 << (P::nbit - i)) + 1, ahk[i]);
        for (int j = 0; j < (P::k+1) * P::n; j++) res[0][j] += evaledauto[0][j];
    }
}

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
        for (int j = 0; j < (P::k+1) * P::n; j++)
            res[num_func - 1][0][j] += evaledauto[0][j];
    }
    for (int i = 0; i < num_func; i++) {
        TRLWE<P> evaledauto;
        EvalAuto<P>(evaledauto, res[num_func - 1], (1 << (P::nbit - i)) + 1,
                    privks[i]);
        for (int j = 0; j < (P::k+1) * P::n; j++)
            res[i][0][j] += res[num_func - 1][0][j] + evaledauto[0][j];
    }
}

template <class P>
void PrivKeySwitch(TRLWE<typename P::targetP> &res,
                   const TLWE<typename P::domainP> &tlwe,
                   const PrivateKeySwitchingKey<P> &privksk)
{
    constexpr typename P::domainP::T roundoffset =
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                 (1 + P::basebit * P::t));

    // Koga's Optimization
    constexpr typename P::domainP::T offset = iksoffsetgen<P>();
    constexpr typename P::domainP::T mask = (1ULL << P::basebit) - 1;
    constexpr typename P::domainP::T halfbase = 1ULL << (P::basebit - 1);
    res = {};
    for (int i = 0; i <= P::domainP::k * P::domainP::n; i++) {
        const typename P::domainP::T aibar = tlwe[i] + offset + roundoffset;

        for (int j = 0; j < P::t; j++) {
            const int32_t aij =
                ((aibar >>
                  (std::numeric_limits<typename P::domainP::T>::digits -
                   (j + 1) * P::basebit)) &
                 mask) -
                halfbase;

            if (aij > 0)
                for (int k = 0; k < P::targetP::k + 1; k++)
                    for (int p = 0; p < P::targetP::n; p++)
                        res[k][p] -= privksk[i][j][aij - 1][k][p];
            else if (aij < 0)
                for (int k = 0; k < P::targetP::k + 1; k++)
                    for (int p = 0; p < P::targetP::n; p++)
                        res[k][p] += privksk[i][j][abs(aij) - 1][k][p];
        }
    }
}

template <class P>
void SubsetPrivKeySwitch(TRLWE<typename P::targetP> &res,
                         const TLWE<typename P::targetP> &tlwe,
                         const SubsetPrivateKeySwitchingKey<P> &privksk)
{
    constexpr uint32_t mask = (1 << P::basebit) - 1;
    constexpr uint64_t prec_offset =
        1ULL << (std::numeric_limits<typename P::targetP::T>::digits -
                 (1 + P::basebit * P::t));

    res = {};
    for (int i = 0; i <= P::targetP::k * P::targetP::n; i++) {
        const typename P::targetP::T aibar = tlwe[i] + prec_offset;

        for (int j = 0; j < P::t; j++) {
            const typename P::domainP::T aij =
                (aibar >> (std::numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;

            if (aij != 0) {
                for (int p = 0; p < P::targetP::n; p++)
                    for (int k = 0; k < P::targetP::k + 1; k++)
                        res[k][p] -= privksk[i][j][aij - 1][k][p];
            }
        }
    }
}

}  // namespace TFHEpp