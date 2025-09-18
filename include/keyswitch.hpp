#pragma once

#include <array>
#include <bit>
#include <span>

#include "params.hpp"
#include "tlwe.hpp"
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
                TLWESub<typename P::targetP>(res, res, ksk[i][j][aij - 1]);
            else if (aij < 0)
                TLWEAdd<typename P::targetP>(res, res, ksk[i][j][-aij - 1]);
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
                    TLWESub<typename P::targetP>(res[cat], res[cat],
                                                 ksk[i][j][aij - 1]);
                else if (aij < 0)
                    TLWEAdd<typename P::targetP>(res[cat], res[cat],
                                                 ksk[i][j][-aij - 1]);
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
            tlwe[i + P::targetP::k * P::targetP::n] + prec_offset;
        for (int j = 0; j < P::t; j++) {
            const uint32_t aij =
                (aibar >> (std::numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;
            if (aij != 0)
                TLWESub<typename P::targetP>(res, res, ksk[i][j][aij - 1]);
        }
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
                TRLWESub<typename P::targetP>(res, res, privksk[i][j][aij - 1]);
            // for (int k = 0; k < P::targetP::k + 1; k++)
            //     for (int p = 0; p < P::targetP::n; p++)
            //         res[k][p] -= privksk[i][j][aij - 1][k][p];
            else if (aij < 0)
                TRLWEAdd<typename P::targetP>(res, res,
                                              privksk[i][j][-aij - 1]);
            // for (int k = 0; k < P::targetP::k + 1; k++)
            //     for (int p = 0; p < P::targetP::n; p++)
            //         res[k][p] += privksk[i][j][abs(aij) - 1][k][p];
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
                (aibar >> (std::numeric_limits<typename P::targetP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;

            if (aij != 0)
                TRLWEAdd<typename P::targetP>(res, res, privksk[i][j][aij - 1]);
            // for (int p = 0; p < P::targetP::n; p++)
            // for (int k = 0; k < P::targetP::k + 1; k++)
            // res[k][p] -= privksk[i][j][aij - 1][k][p];
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
                TRLWESub<typename P::targetP>(res, res, iksk[i][j][aij - 1]);
            // for (int l = 0; l < P::targetP::k + 1; l++)
            // for (int k = 0; k < P::targetP::n; k++)
            // res[l][k] -= iksk[i][j][aij - 1][l][k];
        }
    }
}

template <class P>
void EvalAuto(TRLWE<P> &res, const TRLWE<P> &trlwe, const int d,
              const EvalAutoKey<P> &autokey)
{
    res = {};
    Automorphism<P>(res[P::k], trlwe[P::k], d);

    for (int i = 0; i < P::k; i++) {
        Polynomial<P> temppoly;
        TRLWE<P> temptrlwe;
        Automorphism<P>(temppoly, trlwe[i], d);
        halftrgswfftExternalProduct<P>(temptrlwe, temppoly, autokey[i]);
        TRLWESub<P>(res, res, temptrlwe);
    }
}

// https://eprint.iacr.org/2024/1318
// Reversed order but this is easily proved by packing trivial all 0 TRLWE.
// TODO: They says we should divide by N first, not by 2 for each step. Why?
template <class P>
void AnnihilateKeySwitching(TRLWE<P> &res, const TRLWE<P> &trlwe,
                            const AnnihilateKey<P> &ahk)
{
    res = trlwe;
    // for (int j = 0; j < (P::k + 1) * P::n; j++) res[0][j] /= P::n;
    for (int i = 0; i < P::nbit; i++) {
        for (int j = 0; j < (P::k + 1) * P::n; j++) res[0][j] /= 2;
        TRLWE<P> evaledauto;
        EvalAuto<P>(evaledauto, res, (1 << (i + 1)) + 1, ahk[i]);
        TRLWEAdd<P>(res, res, evaledauto);
    }
}

// template <class P, uint num_func>
// void AnnihilatePrivateKeySwitching(
//     std::array<TRLWE<P>, num_func> &res, const TRLWE<P> &trlwe,
//     const AnnihilateKey<P> &ahk,
//     const std::array<TRGSWFFT<P>, num_func> &privks)
// {
//     static_assert(num_func > 0, "num_func must be bigger than 0");
//     res[num_func - 1] = trlwe;
//     for (int i = 0; i < P::nbit - 1; i++) {
//         TRLWE<P> evaledauto;
//         EvalAuto<P>(evaledauto, res[num_func - 1], (1 << (P::nbit - i)) + 1,
//                     ahk[i]);
//         for (int j = 0; j < (P::k + 1) * P::n; j++)
//             res[num_func - 1][0][j] += evaledauto[0][j];
//     }
//     for (int i = 0; i < num_func; i++) {
//         TRLWE<P> evaledauto;
//         EvalAuto<P>(evaledauto, res[num_func - 1], (1 << (P::nbit - i)) + 1,
//                     privks[i]);
//         for (int j = 0; j < (P::k + 1) * P::n; j++)
//             res[i][0][j] += res[num_func - 1][0][j] + evaledauto[0][j];
//     }
// }

// template <class P, uint num_tlwe>
// void AnnihilatePacking(TRLWE<P> &res, const std::array<TLWE<P>, num_tlwe>
// &tlwes,
//                             const AnnihilateKey<P> &ahk)
// {
//     static_assert(std::has_single_bit(num_tlwe), "Currently, num_tlwe must be
//     power of 2"); std::array<TRLWE<P>, num_tlwe> trlwes; constexpr uint l =
//     std::count_zero(num_tlwe); for (int i = 0; i < num_tlwe; i++) {
//         InvSampleExtractIndex<P>(trlwes[i], tlwes[i], 0);
//         for (int j = 0; j <= P::k * P::n; j++)//rest are known to be 0
//             trlwes[i][0][j] /= P::n;
//     }
//     // Using res as a temporary variable
//     for (int i = 0; i < l; i++){
//         constexpr uint stride = 1 << (l - i - 1);
//         for(int j = 0; j < stride; j++){
//             PolynomialMulByXai<P>(res, trlwes[stride+j], P::n >> i);
//             for(int k = 0; i < (P::k+1) * P::n; k++)
//                 trlwes[stride+j][k] = trlwes[j][k] - res[k];
//             for(int k = 0; i < (P::k+1) * P::n; k++)
//                 trlwes[j][k] += res[k];
//             EvalAuto<P>(res, trlwes[stride+j], (1 << (P::nbit - i)) + 1,
//             ahk[i]); for(int k = 0; i < (P::k+1) * P::n; k++)
//                 trlwes[j][k] += res[k];
//         }
//     }
//     res = trlwes[0];
//     // using trlews[0] and trlwes[1] as temporary variables
//     for (int i = l; i < P::nbit; i++) {
//         PolynomialMulByXai<P>(res, trlwes[(1<<i)+j], P::n >> i);
//         EvalAuto<P>(evaledauto, res, (1 << (P::nbit - i)) + 1, ahk[i]);
//         for (int j = 0; j < (P::k + 1) * P::n; j++)
//             res[0][j] += evaledauto[0][j];
//     }
// }

template <class P, class Container>
void PackLWEs(TRLWE<P> &res, const Container &tlwe, const AnnihilateKey<P> &ahk,
              const uint l, const uint offset, const uint interval)
{
    if (l == 0)
        InvSampleExtractIndex<P>(res, tlwe[offset], 0);
    else {
        TRLWE<P> tempeven;
        PackLWEs<P>(tempeven, tlwe, ahk, l - 1, offset, interval * 2);
        TRLWE<P> tempodd;
        PackLWEs<P>(tempodd, tlwe, ahk, l - 1, offset + interval, interval * 2);
        TRLWE<P> tempoddmul;
        for (int i = 0; i < P::k + 1; i++) {
            PolynomialMulByXai<P>(tempoddmul[i], tempodd[i], P::n >> l);
            for (int j = 0; j < P::n; j++) {
                tempeven[i][j] /= 2;
                tempoddmul[i][j] /= 2;
                tempodd[i][j] = tempeven[i][j] - tempoddmul[i][j];
            }
        }
        EvalAuto<P>(res, tempodd, (1 << l) + 1, ahk[l - 1]);
        TRLWEAdd<P>(res, res, tempeven, tempoddmul);
        // for (int i = 0; i < P::k + 1; i++)
        //     for (int j = 0; j < P::n; j++)
        //         res[i][j] += tempeven[i][j] + tempoddmul[i][j];
    }
}

template <class P>
void TLWE2TRLWEChensPacking(TRLWE<P> &res, std::vector<TLWE<P>> &tlwe,
                            const AnnihilateKey<P> &ahk)
{
    uint l = std::bit_width(tlwe.size()) - 1;
    if (!std::has_single_bit(tlwe.size())) {
        l++;
        tlwe.resize(1 << l);
    }
    PackLWEs<P>(res, tlwe, ahk, l, 0, 1);
    for (int i = l; i < P::nbit; i++) {
        TRLWE<P> evaledauto;
        for (int j = 0; j < (P::k + 1) * P::n; j++) res[0][j] /= 2;
        EvalAuto<P>(evaledauto, res, (1 << (i + 1)) + 1, ahk[i]);
        TRLWEAdd<P>(res, res, evaledauto);
        // for (int j = 0; j < (P::k + 1) * P::n; j++)
        // res[0][j] += evaledauto[0][j];
    }
}

template <class P, uint num_tlwe>
void TLWE2TablePacking(TRLWE<P> &res, std::array<TLWE<P>, num_tlwe> &tlwe,
                       const AnnihilateKey<P> &ahk)
{
    static_assert(std::has_single_bit(num_tlwe),
                  "Currently, num_tlwe must be power of 2");
    constexpr uint l = std::countr_zero(num_tlwe);
    PackLWEs<P>(res, tlwe, ahk, l, 0, 1);
    for (int i = l; i < P::nbit; i++) {
        TRLWE<P> tempmul;
        for (int j = 0; j < P::k + 1; j++)
            PolynomialMulByXai<P>(tempmul[j], res[j], P::n >> (i + 1));
        TRLWE<P> tempsub;
        for (int j = 0; j < (P::k + 1) * P::n; j++) {
            res[0][j] /= 2;
            tempmul[0][j] /= 2;
            tempsub[0][j] = res[0][j] - tempmul[0][j];
            res[0][j] += tempmul[0][j];
        }
        // reuse tempmul
        EvalAuto<P>(tempmul, tempsub, (1 << (i + 1)) + 1, ahk[i]);
        TRLWEAdd<P>(res, res, tempmul);
        // for (int j = 0; j < (P::k + 1) * P::n; j++) res[0][j] +=
        // tempmul[0][j];
    }
}

template <class P, uint num_tlwe, uint num_func>
void TLWE2TablePackingManyLUT(
    TRLWE<P> &res, std::array<std::array<TLWE<P>, num_tlwe>, num_func> &tlwe,
    const AnnihilateKey<P> &ahk)
{
    static_assert(std::has_single_bit(num_tlwe),
                  "Currently, num_tlwe must be power of 2");
    constexpr uint l = std::countr_zero(num_tlwe);
    static_assert(num_func == 2, "Currently, num_func must be 2");
    constexpr uint f = std::countr_zero(num_func);
    std::array<TRLWE<P>, num_func> temptrlwe;
    for (int index = 0; index < num_func; index++) {
        PackLWEs<P>(temptrlwe[index], tlwe[index], ahk, l, 0, 1);
        for (int i = l; i < P::nbit - f; i++) {
            TRLWE<P> tempmul;
            for (int j = 0; j < P::k + 1; j++)
                PolynomialMulByXai<P>(tempmul[j], temptrlwe[index][j],
                                      P::n >> (i + 1));
            TRLWE<P> tempsub;
            for (int j = 0; j < (P::k + 1) * P::n; j++) {
                temptrlwe[index][0][j] /= 2;
                tempmul[0][j] /= 2;
                tempsub[0][j] = temptrlwe[index][0][j] - tempmul[0][j];
                temptrlwe[index][0][j] += tempmul[0][j];
            }
            // reuse tempmul
            EvalAuto<P>(tempmul, tempsub, (1 << (i + 1)) + 1, ahk[i]);
            TRLWEAdd<P>(temptrlwe[index], temptrlwe[index], tempmul);
            // for (int j = 0; j < (P::k + 1) * P::n; j++)
            // temptrlwe[index][0][j] += tempmul[0][j];
        }
    }
    {
        TRLWE<P> tempoddmul;
        for (int i = 0; i < P::k + 1; i++) {
            PolynomialMulByXai<P>(tempoddmul[i], temptrlwe[1][i], 1);
            for (int j = 0; j < P::n; j++) {
                temptrlwe[0][i][j] /= 2;
                tempoddmul[i][j] /= 2;
                temptrlwe[1][i][j] = temptrlwe[0][i][j] - tempoddmul[i][j];
            }
        }
        EvalAuto<P>(res, temptrlwe[1], (1 << P::nbit) + 1, ahk[P::nbit - 1]);
        TRLWEAdd<P>(res, res, temptrlwe[0], tempoddmul);
        // for (int i = 0; i < P::k + 1; i++)
        // for (int j = 0; j < P::n; j++)
        // res[i][j] += temptrlwe[0][i][j] + tempoddmul[i][j];
    }
}

template <class P>
void PackLWEsLSB(TRLWE<P> &res, const std::vector<TLWE<P>> &tlwe,
                 const AnnihilateKey<P> &ahk, const uint l, const uint offset,
                 const uint interval)
{
    if (offset >= tlwe.size())
        res = {};
    else if (l == 0)
        InvSampleExtractIndex<P>(res, tlwe[offset], 0);
    else {
        TRLWE<P> tempeven;
        PackLWEsLSB<P>(tempeven, tlwe, ahk, l - 1, offset, interval * 2);
        TRLWE<P> tempodd;
        PackLWEsLSB<P>(tempodd, tlwe, ahk, l - 1, offset + interval,
                       interval * 2);
        TRLWE<P> tempoddmul;
        for (int i = 0; i < P::k + 1; i++) {
            PolynomialMulByXai<P>(tempoddmul[i], tempodd[i], P::n >> l);
            for (int j = 0; j < P::n; j++) {
                tempeven[i][j] /= 2;
                tempoddmul[i][j] /= 2;
                tempodd[i][j] = tempeven[i][j] - tempoddmul[i][j];
            }
        }
        EvalAuto<P>(res, tempodd, (1 << l) + 1, ahk[l - 1]);
        TRLWEAdd<P>(res, res, tempeven, tempoddmul);
        // for (int i = 0; i < P::k + 1; i++)
        // for (int j = 0; j < P::n; j++)
        // res[i][j] += tempeven[i][j] + tempoddmul[i][j];
    }
}

template <class P>
void TLWE2TRLWEPacking(TRLWE<P> &res, std::vector<TLWE<P>> &tlwe,
                       const AnnihilateKey<P> &ahk)
{
    PackLWEsLSB<P>(res, tlwe, ahk, P::nbit, 0, 1);
}
}  // namespace TFHEpp