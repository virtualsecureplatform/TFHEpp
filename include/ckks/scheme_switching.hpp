#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "ckks.hpp"
#include "tfhe/evalkeygens.hpp"
#include "tfhe/keyswitch.hpp"

namespace TFHEpp {

namespace ckks_scheme_switch_detail {

template <class P>
inline std::int64_t ternarySecretToInt(const typename P::T value)
{
    if (value == typename P::T{0}) return 0;
    if (value == typename P::T{1}) return 1;
    if (value == typename P::T{0} - typename P::T{1}) return -1;
    throw std::invalid_argument(
        "CKKS-to-TFHE key switching requires a ternary source secret");
}

template <class P, std::uint32_t LogQ>
inline void exactActiveTorusMulBySmall(Polynomial<P> &res,
                                       const Polynomial<P> &torus,
                                       const Polynomial<P> &small)
{
    ckks_detail::exactActiveTorusMulBySmall<P, LogQ>(res, torus, small);
}

template <class P, std::uint32_t LogQ>
inline void exactEncryptPolynomialAtLevel(CKKSSeededKeySwitchRow<P, LogQ> &row,
                                          const Polynomial<P> &poly,
                                          const Key<P> &key,
                                          CKKSNoise noise)
{
    for (std::uint32_t i = 0; i < P::n; i++)
        row.body[i] = ckks_detail::reduceToLevel<P, LogQ>(
            poly[i] + ckks_detail::sampleNoise<P, LogQ>(noise));

    auto partkey = std::make_unique<Polynomial<P>>();
    auto mask = std::make_unique<Polynomial<P>>();
    auto mask_phase = std::make_unique<Polynomial<P>>();
    for (std::uint32_t component = 0; component < P::k; component++) {
        row.mask_seeds[component] = ckks_detail::randomSeed();
        ckks_detail::fillSeededUniformPolynomialAtLevel<P, LogQ>(
            *mask, row.mask_seeds[component]);
        for (std::uint32_t i = 0; i < P::n; i++)
            (*partkey)[i] = key[component * P::n + i];
        exactActiveTorusMulBySmall<P, LogQ>(*mask_phase, *mask, *partkey);
        for (std::uint32_t i = 0; i < P::n; i++)
            row.body[i] = ckks_detail::reduceToLevel<P, LogQ>(
                row.body[i] + (*mask_phase)[i]);
    }
}

template <class P, std::uint32_t LogQ>
inline void exactSeededSecretKeySwitchKeyGen(
    CKKSSeededSecretKeySwitchKey<P, LogQ> &switch_key,
    const Key<P> &source_key, const Key<P> &target_key, CKKSNoise noise)
{
    constexpr int width = std::numeric_limits<typename P::T>::digits;
    constexpr std::uint32_t row_count = CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row = CKKSKeySwitchFirstRowForLevel<P, LogQ>();
    auto partkey = std::make_unique<Polynomial<P>>();
    auto gadget = std::make_unique<Polynomial<P>>();
    for (std::uint32_t component = 0; component < P::k; component++) {
        for (std::uint32_t i = 0; i < P::n; i++)
            (*partkey)[i] = source_key[component * P::n + i];
        for (std::uint32_t row = 0; row < row_count; row++) {
            const std::uint32_t full_row = first_row + row;
            const int shift = width -
                (full_row + 1) * static_cast<int>(P::B̅gbit);
            for (std::uint32_t i = 0; i < P::n; i++)
                (*gadget)[i] = ckks_detail::reduceToLevel<P, LogQ>(
                    (*partkey)[i] << shift);
            exactEncryptPolynomialAtLevel<P, LogQ>(
                switch_key[component][row], *gadget, target_key, noise);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void exactSeededKeySwitchRows(
    TRLWE<P> &res, const Polynomial<P> &poly,
    const typename CKKSSeededSecretKeySwitchKey<P, LogQ>::Rows &rows)
{
    constexpr std::uint32_t row_count = CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row = CKKSKeySwitchFirstRowForLevel<P, LogQ>();
    auto centered = std::make_unique<Polynomial<P>>();
    auto decomposed = std::make_unique<std::array<Polynomial<P>, row_count>>();
    ckks_detail::centeredPolynomialAtLevel<P, LogQ>(*centered, poly);
    ckks_detail::baseBbarDecomposePolynomialRows<P, first_row>(
        *decomposed, *centered);

    for (std::uint32_t component = 0; component <= P::k; component++)
        for (std::uint32_t i = 0; i < P::n; i++) res[component][i] = 0;
    auto mask = std::make_unique<Polynomial<P>>();
    auto product = std::make_unique<Polynomial<P>>();
    for (std::uint32_t row = 0; row < row_count; row++) {
        const auto &key_row = rows[row];
        for (std::uint32_t component = 0; component < P::k; component++) {
            ckks_detail::fillSeededUniformPolynomialAtLevel<P, LogQ>(
                *mask, key_row.mask_seeds[component]);
            exactActiveTorusMulBySmall<P, LogQ>(
                *product, *mask, (*decomposed)[row]);
            for (std::uint32_t i = 0; i < P::n; i++)
                res[component][i] = ckks_detail::reduceToLevel<P, LogQ>(
                    res[component][i] + (*product)[i]);
        }
        exactActiveTorusMulBySmall<P, LogQ>(
            *product, key_row.body, (*decomposed)[row]);
        for (std::uint32_t i = 0; i < P::n; i++)
            res[P::k][i] = ckks_detail::reduceToLevel<P, LogQ>(
                res[P::k][i] + (*product)[i]);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void exactSeededSecretKeySwitch(
    CKKSCiphertext<P, LogQ, LogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSSeededSecretKeySwitchKey<P, LogQ> &switch_key)
{
    auto switched = std::make_unique<TRLWE<P>>();
    auto out = std::make_unique<TRLWE<P>>();
    for (std::uint32_t i = 0; i < P::n; i++)
        (*out)[P::k][i] = ct.ct[P::k][i];
    for (std::uint32_t component = 0; component < P::k; component++) {
        exactSeededKeySwitchRows<P, LogQ>(
            *switched, ct.ct[component], switch_key[component]);
        CKKSSubTRLWEInPlace<P, LogQ>(*out, *switched);
    }
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(*out);
    res.ct = *out;
}

}  // namespace ckks_scheme_switch_detail

// Parameters for the extraction/key-switch portion of CKKS -> TFHE scheme
// switching.  The CKKS ciphertext is first decoded into coefficient form (as
// in OpenFHE's EvalCKKStoFHEW), then a coefficient is sample-extracted and
// key-switched to a TFHE LWE secret key.
//
// This deliberately keeps the CKKS and TFHE parameter families independent.
// In particular, CKKSParam's torus can be a multi-limb type while TFHEParam
// can use the conventional 32- or 64-bit torus.
template <class CKKSParam, class TFHEParam, std::uint32_t BaseBit,
          std::uint32_t Levels>
struct CKKSToTFHEKeySwitchParams {
    using domainP = CKKSParam;
    using targetP = TFHEParam;

    static_assert(BaseBit > 0);
    static_assert(Levels > 0);
    static_assert(BaseBit < 31,
                  "Identity-key-switch digit tables use 32-bit indices");
    static_assert(BaseBit * Levels <=
                      std::numeric_limits<typename CKKSParam::T>::digits,
                  "key-switch decomposition exceeds the CKKS torus width");
    static_assert(BaseBit * Levels <=
                      std::numeric_limits<typename TFHEParam::T>::digits,
                  "key-switch decomposition exceeds the TFHE torus width");
    static_assert(std::numeric_limits<typename TFHEParam::T>::digits <= 64,
                  "CKKS-to-TFHE currently supports up to 64-bit TFHE torus");

    static constexpr std::uint32_t basebit = BaseBit;
    static constexpr std::uint32_t t = Levels;
};

template <class CKKSParam, class TFHEParam, std::uint32_t BaseBit,
          std::uint32_t Levels>
using CKKSToTFHEKeySwitchKey = KeySwitchingKey<
    CKKSToTFHEKeySwitchParams<CKKSParam, TFHEParam, BaseBit, Levels>>;

// Generates the RLWE-coefficient-secret -> TFHE-LWE-secret switching key.
// The caller owns both secrets, exactly as for the OpenFHE setup/keygen API.
template <class CKKSParam, class TFHEParam, std::uint32_t BaseBit,
          std::uint32_t Levels>
inline void CKKSToTFHEKeyGen(
    CKKSToTFHEKeySwitchKey<CKKSParam, TFHEParam, BaseBit, Levels> &switch_key,
    const Key<CKKSParam> &ckks_key, const Key<TFHEParam> &tfhe_key)
{
    using Params =
        CKKSToTFHEKeySwitchParams<CKKSParam, TFHEParam, BaseBit, Levels>;
    using TargetT = typename TFHEParam::T;
    constexpr std::uint32_t target_bits =
        std::numeric_limits<TargetT>::digits;
    for (std::uint32_t ring = 0; ring < CKKSParam::k; ring++)
        for (std::uint32_t coefficient = 0; coefficient < CKKSParam::n;
             coefficient++) {
            const std::int64_t secret =
                ckks_scheme_switch_detail::ternarySecretToInt<CKKSParam>(
                    ckks_key[ring * CKKSParam::n + coefficient]);
            for (std::uint32_t level = 0; level < Params::t; level++)
                for (std::uint32_t digit = 0;
                     digit < (std::uint32_t{1} << (Params::basebit - 1));
                     digit++) {
                    const TargetT weight =
                        TargetT{digit + 1}
                        << (target_bits - (level + 1) * Params::basebit);
                    const TargetT message =
                        secret < 0 ? TargetT{0} - weight :
                                     (secret == 0 ? TargetT{0} : weight);
                    tlweSymEncrypt<TFHEParam>(
                        switch_key[ring * CKKSParam::n + coefficient][level]
                                  [digit],
                        message, tfhe_key);
                }
        }
}

// Compact CKKS -> TFHE switching material following OpenFHE's intermediate
// RLWE construction.  Instead of keeping one TLWE key-switch ciphertext for
// every CKKS secret coefficient, it key-switches the CKKS ciphertext once to
// a CKKS-ring secret whose first LWE entries are the TFHE secret and whose
// remaining entries are zero.  Sample extraction can then discard those zero
// secret coordinates and directly yield a TFHE TLWE.
//
// The underlying key is seeded, so its public representation is O(N) per
// CKKS gadget row rather than O(N * LWE_dimension) TLWE coefficients.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ>
struct CKKSToTFHERingSwitchKey {
    static constexpr std::uint32_t lwe_dimension =
        TFHEParam::k * TFHEParam::n;
    static_assert(is_multilimb_uint_v<typename CKKSParam::T>,
                  "the compact CKKS-to-TFHE ring switch requires a multi-limb "
                  "CKKS torus for exact low-level key switching");
    static_assert(CKKSParam::n == (std::uint32_t{1} << CKKSParam::nbit) &&
                      CKKSParam::nbit <= 15,
                  "the exact compact CKKS-to-TFHE switch supports rings through 2^15");
    static_assert(CKKSParam::B̅gbit <= 20,
                  "the exact compact CKKS-to-TFHE switch requires small Bbar digits");
    static_assert(lwe_dimension <= CKKSParam::k * CKKSParam::n,
                  "TFHE LWE secret does not fit in the CKKS ring secret");

    CKKSSeededSecretKeySwitchKey<CKKSParam, LogQ> ring_switch{};
};

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ>
inline void CKKSToTFHERingSwitchKeyGen(
    CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, LogQ> &switch_key,
    const Key<CKKSParam> &ckks_key, const Key<TFHEParam> &tfhe_key,
    CKKSNoise noise = {CKKSParam::α, 0})
{
    Key<CKKSParam> padded_tfhe_key{};
    padded_tfhe_key.fill(typename CKKSParam::T{0});
    constexpr std::uint32_t lwe_dimension =
        CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, LogQ>::lwe_dimension;
    for (std::uint32_t i = 0; i < lwe_dimension; i++) {
        const std::int64_t secret =
            ckks_scheme_switch_detail::ternarySecretToInt<TFHEParam>(
                tfhe_key[i]);
        padded_tfhe_key[i] =
            secret < 0 ? typename CKKSParam::T{0} - typename CKKSParam::T{1}
                       : static_cast<typename CKKSParam::T>(secret);
    }
    ckks_scheme_switch_detail::exactSeededSecretKeySwitchKeyGen<CKKSParam,
                                                                  LogQ>(
        switch_key.ring_switch, ckks_key, padded_tfhe_key, noise);
}

namespace ckks_scheme_switch_detail {

// Modulus-switch an *active-level* CKKS torus value to a conventional TFHE
// torus.  Work directly at InputLogQ: going through the enclosing CKKS torus
// width is needlessly error-prone for multi-limb torus types.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ>
inline typename TFHEParam::T ckksActiveTorusToTFHE(
    const typename CKKSParam::T value)
{
    using SourceT = typename CKKSParam::T;
    using TargetT = typename TFHEParam::T;
    constexpr std::uint32_t target_bits =
        std::numeric_limits<TargetT>::digits;
    const SourceT active = ckks_detail::reduceToLevel<CKKSParam, InputLogQ>(
        value);
    if constexpr (InputLogQ > target_bits) {
        constexpr SourceT rounding =
            SourceT{1} << (InputLogQ - target_bits - 1);
        return static_cast<TargetT>(static_cast<std::uint64_t>(
            (active + rounding) >> (InputLogQ - target_bits)));
    }
    else if constexpr (InputLogQ < target_bits) {
        return static_cast<TargetT>(static_cast<std::uint64_t>(active))
               << (target_bits - InputLogQ);
    }
    else {
        return static_cast<TargetT>(static_cast<std::uint64_t>(active));
    }
}

template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta>
inline void extractRingSwitchedCKKSCoeffToTFHE(
    TLWE<TFHEParam> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &switched,
    std::uint32_t coefficient)
{
    if (coefficient >= CKKSParam::n)
        throw std::out_of_range("CKKS coefficient index is out of range");
    TLWE<CKKSParam> extracted{};
    SampleExtractIndex<CKKSParam>(extracted, switched.ct,
                                  static_cast<int>(coefficient));
    out = {};
    constexpr std::uint32_t lwe_dimension =
        CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam,
                                InputLogQ>::lwe_dimension;
    for (std::uint32_t i = 0; i < lwe_dimension; i++)
        out[i] = ckksActiveTorusToTFHE<CKKSParam, TFHEParam, InputLogQ>(
            extracted[i]);
    out[lwe_dimension] =
        ckksActiveTorusToTFHE<CKKSParam, TFHEParam, InputLogQ>(
            extracted[CKKSParam::k * CKKSParam::n]);
}

}  // namespace ckks_scheme_switch_detail

// Switch a coefficient-domain CKKS ciphertext through the compact ring key.
// The target TFHE phase has the same scale interpretation as
// CKKSCoeffToTFHE, but the ciphertext is key-switched only once regardless of
// how many coefficients will subsequently be extracted.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta>
inline void CKKSCoeffToTFHEViaRingSwitch(
    TLWE<TFHEParam> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t coefficient,
    const CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, InputLogQ>
        &switch_key)
{
    static_assert(InputLogQ <= ckks_detail::torus_width_v<CKKSParam>);
    static_assert(InputLogDelta < InputLogQ);
    static_assert(static_cast<std::int32_t>(
                      std::numeric_limits<typename TFHEParam::T>::digits) -
                      static_cast<std::int32_t>(InputLogQ) +
                      static_cast<std::int32_t>(InputLogDelta) > 0,
                  "TFHE torus is too small for the requested CKKS precision");
    auto switched = std::make_unique<
        CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta>>();
    ckks_scheme_switch_detail::exactSeededSecretKeySwitch<
        CKKSParam, InputLogQ, InputLogDelta>(*switched, in,
                                              switch_key.ring_switch);
    ckks_scheme_switch_detail::extractRingSwitchedCKKSCoeffToTFHE<
        CKKSParam, TFHEParam, InputLogQ, InputLogDelta>(out, *switched,
                                                         coefficient);
}

// Extract the first num_coefficients coefficient-domain CKKS values using one
// compact ring key switch.  This is the intended batched entry point: all
// sample extractions reuse the same switched RLWE ciphertext.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta>
inline void CKKSCoeffsToTFHEViaRingSwitch(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_coefficients,
    const CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, InputLogQ>
        &switch_key)
{
    static_assert(InputLogQ <= ckks_detail::torus_width_v<CKKSParam>);
    static_assert(InputLogDelta < InputLogQ);
    static_assert(static_cast<std::int32_t>(
                      std::numeric_limits<typename TFHEParam::T>::digits) -
                      static_cast<std::int32_t>(InputLogQ) +
                      static_cast<std::int32_t>(InputLogDelta) > 0,
                  "TFHE torus is too small for the requested CKKS precision");
    if (num_coefficients == 0 || num_coefficients > CKKSParam::n)
        throw std::out_of_range("invalid number of CKKS coefficients");

    auto switched = std::make_unique<
        CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta>>();
    ckks_scheme_switch_detail::exactSeededSecretKeySwitch<
        CKKSParam, InputLogQ, InputLogDelta>(*switched, in,
                                              switch_key.ring_switch);
    out.resize(num_coefficients);
    for (std::uint32_t i = 0; i < num_coefficients; i++)
        ckks_scheme_switch_detail::extractRingSwitchedCKKSCoeffToTFHE<
            CKKSParam, TFHEParam, InputLogQ, InputLogDelta>(out[i], *switched,
                                                             i);
}

// Batched compact extraction at caller-selected coefficient indices.  The
// output order is exactly coefficient_indices' order, which accommodates the
// spaced extraction layout used when a CKKS application employs sparse slots.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta>
inline void CKKSCoeffIndicesToTFHEViaRingSwitch(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    const std::vector<std::uint32_t> &coefficient_indices,
    const CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, InputLogQ>
        &switch_key)
{
    static_assert(InputLogQ <= ckks_detail::torus_width_v<CKKSParam>);
    static_assert(InputLogDelta < InputLogQ);
    static_assert(static_cast<std::int32_t>(
                      std::numeric_limits<typename TFHEParam::T>::digits) -
                      static_cast<std::int32_t>(InputLogQ) +
                      static_cast<std::int32_t>(InputLogDelta) > 0,
                  "TFHE torus is too small for the requested CKKS precision");
    if (coefficient_indices.empty())
        throw std::out_of_range("no CKKS coefficients were requested");
    for (const std::uint32_t coefficient : coefficient_indices)
        if (coefficient >= CKKSParam::n)
            throw std::out_of_range("CKKS coefficient index is out of range");

    auto switched = std::make_unique<
        CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta>>();
    ckks_scheme_switch_detail::exactSeededSecretKeySwitch<
        CKKSParam, InputLogQ, InputLogDelta>(*switched, in,
                                              switch_key.ring_switch);
    out.resize(coefficient_indices.size());
    for (std::size_t i = 0; i < coefficient_indices.size(); i++)
        ckks_scheme_switch_detail::extractRingSwitchedCKKSCoeffToTFHE<
            CKKSParam, TFHEParam, InputLogQ, InputLogDelta>(
            out[i], *switched, coefficient_indices[i]);
}

// Precision-preserving compact extraction for large CKKS rings.  Some CKKS
// backends multiply key-switch rows through a floating FFT; switching at a
// small active modulus can leave its rounding error inside the TFHE word.
// Lift the ciphertext scale to SwitchLogQ first, key-switch there, and retain
// only the target TFHE bits during extraction.  The output fixed-point scale
// is unchanged from the SourceLogQ representation.
template <class CKKSParam, class TFHEParam, std::uint32_t SourceLogQ,
          std::uint32_t SourceLogDelta, std::uint32_t SwitchLogQ>
inline void CKKSCoeffsToTFHEViaRingSwitchLifted(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, SourceLogQ, SourceLogDelta> &in,
    std::uint32_t num_coefficients,
    const CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, SwitchLogQ>
        &switch_key)
{
    static_assert(SwitchLogQ >= SourceLogQ,
                  "ring-switch level must not be below the source level");
    static_assert(SwitchLogQ <= ckks_detail::torus_width_v<CKKSParam>);
    constexpr std::uint32_t lift_bits = SwitchLogQ - SourceLogQ;
    constexpr std::uint32_t lifted_log_delta = SourceLogDelta + lift_bits;
    static_assert(lifted_log_delta < SwitchLogQ,
                  "lifted CKKS scale must remain below the switching level");
    static_assert(static_cast<std::int32_t>(
                      std::numeric_limits<typename TFHEParam::T>::digits) -
                      static_cast<std::int32_t>(SourceLogQ) +
                      static_cast<std::int32_t>(SourceLogDelta) > 0,
                  "TFHE torus is too small for the requested CKKS precision");
    if (num_coefficients == 0 || num_coefficients > CKKSParam::n)
        throw std::out_of_range("invalid number of CKKS coefficients");

    auto raised = std::make_unique<
        CKKSCiphertext<CKKSParam, SwitchLogQ, SourceLogDelta>>();
    if constexpr (lift_bits == 0) {
        *raised = in;
    }
    else {
        CKKSModRaise<CKKSParam, SourceLogQ, SwitchLogQ, SourceLogDelta>(
            *raised, in);
    }
    auto lifted = std::make_unique<
        CKKSCiphertext<CKKSParam, SwitchLogQ, lifted_log_delta>>();
    if constexpr (lift_bits == 0) {
        *lifted = *raised;
    }
    else {
        CKKSScaleUpAtSameLevel<CKKSParam, SwitchLogQ, SourceLogDelta,
                                lifted_log_delta>(*lifted, *raised);
    }

    auto switched = std::make_unique<
        CKKSCiphertext<CKKSParam, SwitchLogQ, lifted_log_delta>>();
    ckks_scheme_switch_detail::exactSeededSecretKeySwitch<
        CKKSParam, SwitchLogQ, lifted_log_delta>(*switched, *lifted,
                                                  switch_key.ring_switch);
    out.resize(num_coefficients);
    for (std::uint32_t i = 0; i < num_coefficients; i++)
        ckks_scheme_switch_detail::extractRingSwitchedCKKSCoeffToTFHE<
            CKKSParam, TFHEParam, SwitchLogQ, lifted_log_delta>(out[i],
                                                                  *switched, i);
}

// Number of fractional bits in the returned TFHE torus phase.  If a CKKS
// coefficient encrypts x at scale 2^InputLogDelta, the output phase encrypts
// x * 2^output_scale_bits in the TFHE torus.  A positive value is required to
// retain a usable message/noise gap after modulus switching.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta>
inline constexpr std::int32_t CKKSToTFHEOutputScaleBits =
    static_cast<std::int32_t>(std::numeric_limits<typename TFHEParam::T>::digits) -
    static_cast<std::int32_t>(InputLogQ) +
    static_cast<std::int32_t>(InputLogDelta);

// Extract one *coefficient-domain* CKKS value and key-switch it into a TFHE
// TLWE.  Packed CKKS slots must first be homomorphically decoded to
// coefficients; that linear transform is intentionally separate because its
// rotation-key schedule depends on the selected CKKS bootstrap schedule.
//
// The key-switch decomposition operates on the most-significant torus bits.
// TFHEpp CKKS ciphertexts store their active modulus in the low InputLogQ
// bits, so the extracted LWE is centred and lifted before key switching.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t BaseBit,
          std::uint32_t Levels>
inline void CKKSCoeffToTFHE(
    TLWE<TFHEParam> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t coefficient,
    const CKKSToTFHEKeySwitchKey<CKKSParam, TFHEParam, BaseBit, Levels>
        &switch_key)
{
    static_assert(InputLogQ <= ckks_detail::torus_width_v<CKKSParam>);
    static_assert(InputLogDelta < InputLogQ);
    static_assert(BaseBit * Levels >= InputLogQ,
                  "key-switch decomposition must cover the active CKKS modulus");
    static_assert(CKKSToTFHEOutputScaleBits<CKKSParam, TFHEParam, InputLogQ,
                                             InputLogDelta> > 0,
                  "TFHE torus is too small for the requested CKKS precision");

    if (coefficient >= CKKSParam::n)
        throw std::out_of_range("CKKS coefficient index is out of range");

    TLWE<CKKSParam> extracted{};
    SampleExtractIndex<CKKSParam>(extracted, in.ct,
                                  static_cast<int>(coefficient));

    // Convert Z/(2^InputLogQ) to the full CKKS torus so IdentityKeySwitch's
    // MSB-first digit decomposition performs the intended modulus switch.
    // centeredLevelToTorus sign-extends the active value; the following shift
    // also moves its modulus from the low InputLogQ bits to the torus MSBs.
    constexpr std::uint32_t lift_bits =
        ckks_detail::torus_width_v<CKKSParam> - InputLogQ;
    for (auto &value : extracted) {
        value = ckks_detail::centeredLevelToTorus<CKKSParam, InputLogQ>(value);
        if constexpr (lift_bits != 0) value <<= lift_bits;
    }

    using SourceT = typename CKKSParam::T;
    using TargetT = typename TFHEParam::T;
    constexpr std::uint32_t source_bits =
        std::numeric_limits<SourceT>::digits;
    constexpr std::uint32_t target_bits =
        std::numeric_limits<TargetT>::digits;
    constexpr SourceT digit_mask = (SourceT{1} << BaseBit) - SourceT{1};
    constexpr std::int32_t half_base =
        std::int32_t{1} << (BaseBit - 1);
    constexpr SourceT decomposition_offset = [] {
        SourceT value = 0;
        for (std::uint32_t level = 1; level <= Levels; level++)
            value += (SourceT{1} << (BaseBit - 1))
                     << (source_bits - level * BaseBit);
        return value;
    }();
    constexpr SourceT rounding_offset =
        BaseBit * Levels < source_bits
            ? SourceT{1} << (source_bits - 1 - BaseBit * Levels)
            : SourceT{0};

    out = {};
    const SourceT body = extracted[CKKSParam::k * CKKSParam::n];
    if constexpr (source_bits == target_bits) {
        out[TFHEParam::k * TFHEParam::n] =
            static_cast<TargetT>(static_cast<std::uint64_t>(body));
    }
    else if constexpr (source_bits > target_bits) {
        const SourceT rounded =
            body + (SourceT{1} << (source_bits - target_bits - 1));
        out[TFHEParam::k * TFHEParam::n] = static_cast<TargetT>(
            static_cast<std::uint64_t>(rounded >> (source_bits - target_bits)));
    }
    else {
        out[TFHEParam::k * TFHEParam::n] =
            static_cast<TargetT>(static_cast<std::uint64_t>(body))
            << (target_bits - source_bits);
    }

    for (std::uint32_t i = 0; i < CKKSParam::k * CKKSParam::n; i++) {
        const SourceT aibar = extracted[i] + decomposition_offset +
                               rounding_offset;
        for (std::uint32_t level = 0; level < Levels; level++) {
            const SourceT raw =
                (aibar >> (source_bits - (level + 1) * BaseBit)) & digit_mask;
            const std::int32_t digit = static_cast<std::int32_t>(
                                           static_cast<std::uint64_t>(raw)) -
                                       half_base;
            if (digit > 0)
                TLWESub<TFHEParam>(out, out, switch_key[i][level][digit - 1]);
            else if (digit < 0)
                TLWEAdd<TFHEParam>(out, out,
                                    switch_key[i][level][-digit - 1]);
        }
    }
}

// Extract an arbitrary ordered set of coefficient-domain CKKS values through
// the direct RLWE-to-TLWE key switch.  This is useful for OpenFHE's sparse
// packed layout, where decoded values are separated by a deterministic gap.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t BaseBit,
          std::uint32_t Levels>
inline void CKKSCoeffIndicesToTFHE(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    const std::vector<std::uint32_t> &coefficient_indices,
    const CKKSToTFHEKeySwitchKey<CKKSParam, TFHEParam, BaseBit, Levels>
        &switch_key)
{
    if (coefficient_indices.empty())
        throw std::out_of_range("empty CKKS coefficient-index set");
    out.resize(coefficient_indices.size());
    for (std::size_t i = 0; i < coefficient_indices.size(); i++)
        CKKSCoeffToTFHE<CKKSParam, TFHEParam, InputLogQ, InputLogDelta,
                        BaseBit, Levels>(out[i], in, coefficient_indices[i],
                                         switch_key);
}

// OpenFHE extracts every RingDim/(2*num_slots) coefficient after SlotToCoeff
// when a CKKS ciphertext uses fewer than N/2 logical slots.  Derive those
// indices explicitly so the direct and compact paths share the same sparse
// packing convention.
template <class CKKSParam>
inline void CKKSOpenFHEExtractionIndices(
    std::vector<std::uint32_t> &coefficient_indices,
    std::uint32_t num_ciphertexts, std::uint32_t num_slots)
{
    constexpr std::uint32_t slot_count = CKKSParam::n / 2;
    if (num_slots == 0 || num_slots > slot_count ||
        (slot_count % num_slots) != 0 || num_ciphertexts == 0 ||
        num_ciphertexts > num_slots)
        throw std::out_of_range("invalid OpenFHE CKKS extraction layout");
    const std::uint32_t gap = slot_count / num_slots;
    coefficient_indices.resize(num_ciphertexts);
    for (std::uint32_t i = 0; i < num_ciphertexts; i++)
        coefficient_indices[i] = i * gap;
}

// Homomorphically converts the usual CKKS slot representation to the ring
// coefficient representation used by CKKSCoeffToTFHE.  This is the same
// mathematical SlotToCoeff transform used by OpenFHE before RLWE extraction;
// it is exposed separately so applications can select their own BSGS step and
// rotation-key layout.
template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSSlotsToCoeffs(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &out,
    const CKKSCiphertext<P, LogQ, LogDelta> &in, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key)
{
    static_assert(PlainLogDelta > 0);
    static_assert(PlainLogDelta < LogQ);
    if (bsgs_step <= 0 || bsgs_step > static_cast<int>(P::n / 2))
        throw std::invalid_argument("invalid CKKS SlotToCoeff BSGS step");

    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> plan;
    CKKSBuildPackedSlotToCoeffPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, bsgs_step);
    CKKSRealLinearTransform<P, LogQ, LogDelta, PlainLogDelta>(
        out, in, plan, input_galois_key, output_galois_key);
}

// Collect exactly the binary rotation keys consumed by the reference
// SlotToCoeff plan.  This is preferable to a full CKKSGaloisKey for scheme
// switching: its footprint follows the chosen BSGS step rather than every
// rotation supported by the ring.  The input set includes conjugation when
// the real SlotToCoeff transform needs it.
template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSSlotsToCoeffsRotationKeyIndices(
    CKKSRotationKeyIndexSet<P> &input_keys,
    CKKSRotationKeyIndexSet<P> &output_keys, int bsgs_step)
{
    if (bsgs_step <= 0 || bsgs_step > static_cast<int>(P::n / 2))
        throw std::invalid_argument("invalid CKKS SlotToCoeff BSGS step");
    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> plan;
    CKKSBuildPackedSlotToCoeffPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, bsgs_step);
    CKKSClearRotationKeyIndexSet<P>(input_keys);
    CKKSClearRotationKeyIndexSet<P>(output_keys);
    if (plan.has_direct)
        CKKSCollectLinearTransformPlanRotationKeyIndices<
            P, LogQ, LogDelta, PlainLogDelta>(input_keys, output_keys,
                                               plan.direct);
    if (plan.has_conjugate) {
        CKKSMarkConjugationKeyIndex<P>(input_keys);
        CKKSCollectLinearTransformPlanRotationKeyIndices<
            P, LogQ, LogDelta, PlainLogDelta>(input_keys, output_keys,
                                               plan.conjugate);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSSlotsToCoeffsSparseGaloisKeyGen(
    CKKSSparseGaloisKey<P, LogQ> &input_galois_key,
    CKKSSparseGaloisKey<P, LogQ - PlainLogDelta> &output_galois_key,
    const Key<P> &key, int bsgs_step, CKKSNoise noise = {P::α, 0})
{
    CKKSRotationKeyIndexSet<P> input_indices{};
    CKKSRotationKeyIndexSet<P> output_indices{};
    CKKSSlotsToCoeffsRotationKeyIndices<P, LogQ, LogDelta, PlainLogDelta>(
        input_indices, output_indices, bsgs_step);
    CKKSSparseGaloisKeyGen<P, LogQ>(input_galois_key, key, input_indices,
                                    noise);
    CKKSSparseGaloisKeyGen<P, LogQ - PlainLogDelta>(
        output_galois_key, key, output_indices, noise);
}

// Hybrid variant of the sparse schedule.  Baby rotations can be composed from
// binary keys while BSGS giant rotations can use a direct automorphism key.
// The collector retains whichever representation has no more key entries at
// each level (and always retains binary conjugation when it is required).
template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSSlotsToCoeffsHybridRotationKeyIndices(
    CKKSRotationKeyIndexSet<P> &input_binary,
    CKKSDirectRotationKeyIndexSet<P> &input_direct,
    CKKSRotationKeyIndexSet<P> &output_binary,
    CKKSDirectRotationKeyIndexSet<P> &output_direct, int bsgs_step)
{
    if (bsgs_step <= 0 || bsgs_step > static_cast<int>(P::n / 2))
        throw std::invalid_argument("invalid CKKS SlotToCoeff BSGS step");
    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> plan;
    CKKSBuildPackedSlotToCoeffPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, bsgs_step);
    CKKSClearRotationKeyIndexSet<P>(input_binary);
    CKKSClearDirectRotationKeyIndexSet<P>(input_direct);
    CKKSClearRotationKeyIndexSet<P>(output_binary);
    CKKSClearDirectRotationKeyIndexSet<P>(output_direct);

    const auto collect = [&](const auto &linear_plan) {
        CKKSCollectLinearTransformPlanRotationKeyIndices<
            P, LogQ, LogDelta, PlainLogDelta>(input_binary, output_binary,
                                               linear_plan);
        std::vector<bool> baby_used;
        CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta,
                                              PlainLogDelta>(baby_used,
                                                              linear_plan);
        for (int baby = 1; baby < linear_plan.k_step; baby++)
            if (baby_used[static_cast<std::size_t>(baby)])
                CKKSMarkDirectRotationKeyIndex<P>(input_direct, baby);
        for (const auto &group : linear_plan.groups)
            if (group.giant_step != 0)
                CKKSMarkDirectRotationKeyIndex<P>(
                    output_direct, linear_plan.k_step * group.giant_step);
    };
    if (plan.has_direct) collect(plan.direct);
    if (plan.has_conjugate) {
        CKKSMarkConjugationKeyIndex<P>(input_binary);
        collect(plan.conjugate);
    }

    const bool keep_input_direct =
        CKKSDirectRotationKeyIndexSetCount<P>(input_direct) <=
        CKKSRotationKeyIndexSetCount<P>(input_binary) -
            (plan.has_conjugate ? 1U : 0U);
    if (keep_input_direct) {
        const bool conjugate = input_binary[P::nbit];
        CKKSClearRotationKeyIndexSet<P>(input_binary);
        input_binary[P::nbit] = conjugate;
    }
    else {
        CKKSClearDirectRotationKeyIndexSet<P>(input_direct);
    }

    const bool keep_output_direct =
        CKKSDirectRotationKeyIndexSetCount<P>(output_direct) <=
        CKKSRotationKeyIndexSetCount<P>(output_binary);
    if (keep_output_direct) {
        CKKSClearRotationKeyIndexSet<P>(output_binary);
    }
    else {
        CKKSClearDirectRotationKeyIndexSet<P>(output_direct);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSSlotsToCoeffsHybridSparseGaloisKeyGen(
    CKKSHybridSparseGaloisKey<P, LogQ> &input_galois_key,
    CKKSHybridSparseGaloisKey<P, LogQ - PlainLogDelta> &output_galois_key,
    const Key<P> &key, int bsgs_step, CKKSNoise noise = {P::α, 0})
{
    CKKSRotationKeyIndexSet<P> input_binary{};
    CKKSDirectRotationKeyIndexSet<P> input_direct{};
    CKKSRotationKeyIndexSet<P> output_binary{};
    CKKSDirectRotationKeyIndexSet<P> output_direct{};
    CKKSSlotsToCoeffsHybridRotationKeyIndices<P, LogQ, LogDelta,
                                               PlainLogDelta>(
        input_binary, input_direct, output_binary, output_direct, bsgs_step);
    CKKSHybridSparseGaloisKeyGen<P, LogQ>(
        input_galois_key, key, input_binary, input_direct, noise);
    CKKSHybridSparseGaloisKeyGen<P, LogQ - PlainLogDelta>(
        output_galois_key, key, output_binary, output_direct, noise);
}

// OpenFHE-style CKKS-slot -> TFHE extraction.  The output contains the first
// num_coefficients decoded ring coefficients.  It intentionally exposes the
// SlotToCoeff result level: that level is the intermediate CKKS modulus which
// must fit the target TFHE torus, matching OpenFHE's Qswitch constraint.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t BaseBit, std::uint32_t Levels, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSSlotsToTFHE(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_coefficients, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHEKeySwitchKey<CKKSParam, TFHEParam, BaseBit, Levels>
        &switch_key)
{
    using Decoded =
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>;
    static_assert(BaseBit * Levels >= Decoded::log_q,
                  "key-switch decomposition must cover the decoded CKKS modulus");
    static_assert(CKKSToTFHEOutputScaleBits<CKKSParam, TFHEParam,
                                             Decoded::log_q,
                                             InputLogDelta> > 0,
                  "TFHE torus is too small for the decoded CKKS precision");

    if (num_coefficients == 0 || num_coefficients > CKKSParam::n)
        throw std::out_of_range("invalid number of decoded CKKS coefficients");

    auto decoded = std::make_unique<Decoded>();
    CKKSSlotsToCoeffs<CKKSParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        *decoded, in, bsgs_step, input_galois_key, output_galois_key);

    out.resize(num_coefficients);
    for (std::uint32_t i = 0; i < num_coefficients; i++)
        CKKSCoeffToTFHE<CKKSParam, TFHEParam, Decoded::log_q,
                        InputLogDelta, BaseBit, Levels>(out[i], *decoded, i,
                                                        switch_key);
}

// Direct SlotToCoeff extraction with OpenFHE's sparse packed-slot gap.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t BaseBit, std::uint32_t Levels, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSSlotsToTFHEOpenFHEPacked(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_ciphertexts, std::uint32_t num_slots, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHEKeySwitchKey<CKKSParam, TFHEParam, BaseBit, Levels>
        &switch_key)
{
    using Decoded =
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>;
    std::vector<std::uint32_t> indices;
    CKKSOpenFHEExtractionIndices<CKKSParam>(indices, num_ciphertexts,
                                             num_slots);
    auto decoded = std::make_unique<Decoded>();
    CKKSSlotsToCoeffs<CKKSParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        *decoded, in, bsgs_step, input_galois_key, output_galois_key);
    CKKSCoeffIndicesToTFHE<CKKSParam, TFHEParam, Decoded::log_q,
                           InputLogDelta, BaseBit, Levels>(out, *decoded,
                                                            indices, switch_key);
}

// Compact OpenFHE-style slot extraction.  SlotToCoeff is evaluated once,
// followed by one RLWE key switch and num_coefficients sample extractions.
// Generate ring_switch at Decoded::log_q, not the input CKKS level.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          class InputGaloisKey, class OutputGaloisKey>
inline void CKKSSlotsToTFHEViaRingSwitch(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_coefficients, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHERingSwitchKey<
        CKKSParam, TFHEParam,
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>::log_q> &ring_switch)
{
    using Decoded =
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>;
    auto decoded = std::make_unique<Decoded>();
    CKKSSlotsToCoeffs<CKKSParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        *decoded, in, bsgs_step, input_galois_key, output_galois_key);
    CKKSCoeffsToTFHEViaRingSwitch<CKKSParam, TFHEParam, Decoded::log_q,
                                  InputLogDelta>(out, *decoded,
                                                 num_coefficients, ring_switch);
}

// Slot-to-coefficient decode followed by compact selected-index extraction.
// This is the counterpart to OpenFHE's configurable sample-extraction gap.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          class InputGaloisKey, class OutputGaloisKey>
inline void CKKSSlotsToTFHEViaRingSwitchAtIndices(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    const std::vector<std::uint32_t> &coefficient_indices, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHERingSwitchKey<
        CKKSParam, TFHEParam,
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>::log_q> &ring_switch)
{
    using Decoded =
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>;
    auto decoded = std::make_unique<Decoded>();
    CKKSSlotsToCoeffs<CKKSParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        *decoded, in, bsgs_step, input_galois_key, output_galois_key);
    CKKSCoeffIndicesToTFHEViaRingSwitch<CKKSParam, TFHEParam,
                                        Decoded::log_q, InputLogDelta>(
        out, *decoded, coefficient_indices, ring_switch);
}

// Compact SlotToCoeff extraction with the same RingDim/(2*num_slots) spacing
// used by OpenFHE's EvalCKKStoFHEW.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          class InputGaloisKey, class OutputGaloisKey>
inline void CKKSSlotsToTFHEViaRingSwitchOpenFHEPacked(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_ciphertexts, std::uint32_t num_slots, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHERingSwitchKey<
        CKKSParam, TFHEParam,
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>::log_q> &ring_switch)
{
    std::vector<std::uint32_t> indices;
    CKKSOpenFHEExtractionIndices<CKKSParam>(indices, num_ciphertexts,
                                             num_slots);
    CKKSSlotsToTFHEViaRingSwitchAtIndices<
        CKKSParam, TFHEParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        out, in, indices, bsgs_step, input_galois_key, output_galois_key,
        ring_switch);
}

// OpenFHE's EvalCKKStoFHEW performs a modulus switch after SlotToCoeff and
// before its cross-key switch.  In TFHEpp's power-of-two active-level model,
// this is a level reduction to ExtractionLogQ.  Expose it explicitly so the
// extraction modulus can be chosen for the target TFHE precision/security
// trade-off rather than being tied to the SlotToCoeff output level.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t ExtractionLogQ, std::uint32_t BaseBit,
          std::uint32_t Levels, class InputGaloisKey, class OutputGaloisKey>
inline void CKKSSlotsToTFHEAtLevel(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_coefficients, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHEKeySwitchKey<CKKSParam, TFHEParam, BaseBit, Levels>
        &switch_key)
{
    using Decoded =
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>;
    static_assert(ExtractionLogQ <= Decoded::log_q,
                  "extraction level cannot exceed SlotToCoeff output level");
    static_assert(InputLogDelta < ExtractionLogQ,
                  "extraction level must retain the CKKS scale");
    static_assert(BaseBit * Levels >= ExtractionLogQ,
                  "key-switch decomposition must cover the extraction level");
    if (num_coefficients == 0 || num_coefficients > CKKSParam::n)
        throw std::out_of_range("invalid number of decoded CKKS coefficients");

    auto decoded = std::make_unique<Decoded>();
    CKKSSlotsToCoeffs<CKKSParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        *decoded, in, bsgs_step, input_galois_key, output_galois_key);
    auto reduced = std::make_unique<
        CKKSCiphertext<CKKSParam, ExtractionLogQ, InputLogDelta>>();
    CKKSLevelReduce<CKKSParam, Decoded::log_q, ExtractionLogQ,
                    InputLogDelta>(*reduced, *decoded);
    out.resize(num_coefficients);
    for (std::uint32_t i = 0; i < num_coefficients; i++)
        CKKSCoeffToTFHE<CKKSParam, TFHEParam, ExtractionLogQ,
                        InputLogDelta, BaseBit, Levels>(out[i], *reduced, i,
                                                        switch_key);
}

// Explicit-Q' direct extraction at arbitrary decoded coefficient indices.
// This composes both non-consecutive OpenFHE sample extraction and its
// post-SlotToCoeff modulus switch.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t ExtractionLogQ, std::uint32_t BaseBit,
          std::uint32_t Levels, class InputGaloisKey, class OutputGaloisKey>
inline void CKKSSlotsToTFHEAtLevelAtIndices(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    const std::vector<std::uint32_t> &coefficient_indices, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHEKeySwitchKey<CKKSParam, TFHEParam, BaseBit, Levels>
        &switch_key)
{
    using Decoded =
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>;
    static_assert(ExtractionLogQ <= Decoded::log_q,
                  "extraction level cannot exceed SlotToCoeff output level");
    static_assert(InputLogDelta < ExtractionLogQ,
                  "extraction level must retain the CKKS scale");
    static_assert(BaseBit * Levels >= ExtractionLogQ,
                  "key-switch decomposition must cover the extraction level");
    auto decoded = std::make_unique<Decoded>();
    CKKSSlotsToCoeffs<CKKSParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        *decoded, in, bsgs_step, input_galois_key, output_galois_key);
    auto reduced = std::make_unique<
        CKKSCiphertext<CKKSParam, ExtractionLogQ, InputLogDelta>>();
    CKKSLevelReduce<CKKSParam, Decoded::log_q, ExtractionLogQ,
                    InputLogDelta>(*reduced, *decoded);
    CKKSCoeffIndicesToTFHE<CKKSParam, TFHEParam, ExtractionLogQ,
                           InputLogDelta, BaseBit, Levels>(
        out, *reduced, coefficient_indices, switch_key);
}

// OpenFHE's Q'-level direct extraction with its standard sparse-slot gap.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t ExtractionLogQ, std::uint32_t BaseBit,
          std::uint32_t Levels, class InputGaloisKey, class OutputGaloisKey>
inline void CKKSSlotsToTFHEAtLevelOpenFHEPacked(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_ciphertexts, std::uint32_t num_slots, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHEKeySwitchKey<CKKSParam, TFHEParam, BaseBit, Levels>
        &switch_key)
{
    std::vector<std::uint32_t> indices;
    CKKSOpenFHEExtractionIndices<CKKSParam>(indices, num_ciphertexts,
                                             num_slots);
    CKKSSlotsToTFHEAtLevelAtIndices<
        CKKSParam, TFHEParam, InputLogQ, InputLogDelta, PlainLogDelta,
        ExtractionLogQ, BaseBit, Levels>(
        out, in, indices, bsgs_step, input_galois_key, output_galois_key,
        switch_key);
}

// Compact counterpart of CKKSSlotsToTFHEAtLevel.  The seeded ring switch is
// generated at ExtractionLogQ and is amortized across all extractions.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t ExtractionLogQ, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSSlotsToTFHEViaRingSwitchAtLevel(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_coefficients, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, ExtractionLogQ>
        &ring_switch)
{
    using Decoded =
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>;
    static_assert(ExtractionLogQ <= Decoded::log_q,
                  "extraction level cannot exceed SlotToCoeff output level");
    static_assert(InputLogDelta < ExtractionLogQ,
                  "extraction level must retain the CKKS scale");
    auto decoded = std::make_unique<Decoded>();
    CKKSSlotsToCoeffs<CKKSParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        *decoded, in, bsgs_step, input_galois_key, output_galois_key);
    auto reduced = std::make_unique<
        CKKSCiphertext<CKKSParam, ExtractionLogQ, InputLogDelta>>();
    CKKSLevelReduce<CKKSParam, Decoded::log_q, ExtractionLogQ,
                    InputLogDelta>(*reduced, *decoded);
    CKKSCoeffsToTFHEViaRingSwitch<CKKSParam, TFHEParam, ExtractionLogQ,
                                  InputLogDelta>(out, *reduced,
                                                 num_coefficients, ring_switch);
}

// Explicit-Q' compact extraction at arbitrary decoded coefficient indices.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t ExtractionLogQ, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSSlotsToTFHEViaRingSwitchAtLevelAtIndices(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    const std::vector<std::uint32_t> &coefficient_indices, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, ExtractionLogQ>
        &ring_switch)
{
    using Decoded =
        CKKSPlainMulResult<CKKSParam, InputLogQ, InputLogDelta,
                           PlainLogDelta>;
    static_assert(ExtractionLogQ <= Decoded::log_q,
                  "extraction level cannot exceed SlotToCoeff output level");
    static_assert(InputLogDelta < ExtractionLogQ,
                  "extraction level must retain the CKKS scale");
    auto decoded = std::make_unique<Decoded>();
    CKKSSlotsToCoeffs<CKKSParam, InputLogQ, InputLogDelta, PlainLogDelta>(
        *decoded, in, bsgs_step, input_galois_key, output_galois_key);
    auto reduced = std::make_unique<
        CKKSCiphertext<CKKSParam, ExtractionLogQ, InputLogDelta>>();
    CKKSLevelReduce<CKKSParam, Decoded::log_q, ExtractionLogQ,
                    InputLogDelta>(*reduced, *decoded);
    CKKSCoeffIndicesToTFHEViaRingSwitch<
        CKKSParam, TFHEParam, ExtractionLogQ, InputLogDelta>(
        out, *reduced, coefficient_indices, ring_switch);
}

// Compact counterpart of CKKSSlotsToTFHEAtLevelOpenFHEPacked.
template <class CKKSParam, class TFHEParam, std::uint32_t InputLogQ,
          std::uint32_t InputLogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t ExtractionLogQ, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSSlotsToTFHEViaRingSwitchAtLevelOpenFHEPacked(
    std::vector<TLWE<TFHEParam>> &out,
    const CKKSCiphertext<CKKSParam, InputLogQ, InputLogDelta> &in,
    std::uint32_t num_ciphertexts, std::uint32_t num_slots, int bsgs_step,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const CKKSToTFHERingSwitchKey<CKKSParam, TFHEParam, ExtractionLogQ>
        &ring_switch)
{
    std::vector<std::uint32_t> indices;
    CKKSOpenFHEExtractionIndices<CKKSParam>(indices, num_ciphertexts,
                                             num_slots);
    CKKSSlotsToTFHEViaRingSwitchAtLevelAtIndices<
        CKKSParam, TFHEParam, InputLogQ, InputLogDelta, PlainLogDelta,
        ExtractionLogQ>(out, in, indices, bsgs_step, input_galois_key,
                         output_galois_key, ring_switch);
}

// CKKS encryption of a TFHE LWE secret.  This is the switching key used by
// the reverse direction: multiplying it by a public LWE mask homomorphically
// evaluates a·s without exposing s to the evaluator.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta>
struct TFHEToCKKSSwitchKey {
    static constexpr std::uint32_t lwe_dimension =
        TFHEParam::k * TFHEParam::n;
    static_assert(lwe_dimension <= CKKSParam::n,
                  "TFHE LWE secret does not fit in the CKKS ring");

    CKKSCiphertext<CKKSParam, LogQ, LogDelta> encrypted_secret{};
};

namespace ckks_scheme_switch_detail {

template <class P>
inline double tfheSecretCoefficientToDouble(const typename P::T value)
{
    if (value == typename P::T{0}) return 0.0;
    if (value == typename P::T{1}) return 1.0;
    if (value == typename P::T{0} - typename P::T{1}) return -1.0;
    // TFHEpp's built-in TFHE secrets are ternary.  Keeping this fallback makes
    // experimental small-integer secret distributions usable as well.
    return static_cast<double>(value);
}

template <class P>
inline double tfheTorusToUnitInterval(const typename P::T value)
{
    return std::ldexp(static_cast<double>(value),
                      -static_cast<int>(std::numeric_limits<typename P::T>::digits));
}

}  // namespace ckks_scheme_switch_detail

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta>
inline void TFHEToCKKSKeyGen(
    TFHEToCKKSSwitchKey<CKKSParam, TFHEParam, LogQ, LogDelta> &switch_key,
    const Key<CKKSParam> &ckks_key, const Key<TFHEParam> &tfhe_key,
    CKKSNoise noise = {CKKSParam::α, 0})
{
    std::array<double, CKKSParam::n> secret{};
    secret.fill(0.0);
    for (std::uint32_t i = 0;
         i < TFHEToCKKSSwitchKey<CKKSParam, TFHEParam, LogQ,
                                  LogDelta>::lwe_dimension;
         i++)
        secret[i] =
            ckks_scheme_switch_detail::tfheSecretCoefficientToDouble<TFHEParam>(
                tfhe_key[i]);
    ckksEncrypt<CKKSParam, LogQ, LogDelta>(switch_key.encrypted_secret,
                                           secret, ckks_key, noise);
}

// Slot-encoded counterpart of TFHEToCKKSSwitchKey.  A complex CKKS slot can
// carry two real values, but this reference implementation deliberately uses
// only the real part so that every linear-transform diagonal is real.  This
// preserves the real phase required by periodic reduction and supports up to
// N/2 TLWE-secret coefficients in one packed CKKS ciphertext.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta>
struct TFHEToCKKSPackedSwitchKey {
    static constexpr std::uint32_t lwe_dimension =
        TFHEParam::k * TFHEParam::n;
    static_assert(lwe_dimension <= CKKSParam::n / 2,
                  "TFHE LWE secret does not fit in the real CKKS slots");
    static constexpr std::uint32_t max_ciphertexts =
        (CKKSParam::n / 2) / lwe_dimension;

    CKKSCiphertext<CKKSParam, LogQ, LogDelta> encrypted_secret_slots{};
};

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta>
inline void TFHEToCKKSPackedKeyGen(
    TFHEToCKKSPackedSwitchKey<CKKSParam, TFHEParam, LogQ, LogDelta>
        &switch_key,
    const Key<CKKSParam> &ckks_key, const Key<TFHEParam> &tfhe_key,
    CKKSNoise noise = {CKKSParam::α, 0})
{
    CKKSSlotVector<CKKSParam> secret{};
    secret.fill({0.0, 0.0});
    constexpr std::uint32_t lwe_dimension =
        TFHEToCKKSPackedSwitchKey<CKKSParam, TFHEParam, LogQ,
                                  LogDelta>::lwe_dimension;
    for (std::uint32_t block = 0;
         block < TFHEToCKKSPackedSwitchKey<CKKSParam, TFHEParam, LogQ,
                                           LogDelta>::max_ciphertexts;
         block++)
        for (std::uint32_t i = 0; i < lwe_dimension; i++)
            secret[block * lwe_dimension + i] = {
                ckks_scheme_switch_detail::
                    tfheSecretCoefficientToDouble<TFHEParam>(tfhe_key[i]),
                0.0};
    ckksSlotEncrypt<CKKSParam, LogQ, LogDelta>(
        switch_key.encrypted_secret_slots, secret, ckks_key, noise);
}

// Dense-slot counterpart of TFHEToCKKSPackedSwitchKey.  The secret occurs
// once in the first lwe_dimension slots; a rectangular CKKS linear transform
// then evaluates every requested public LWE mask against that same encrypted
// vector.  Unlike the block-replicated reference layout below, this supports
// up to N/2 output phases in one ciphertext, as in OpenFHE's
// EvalPartialHomDecryption.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta>
struct TFHEToCKKSSlotPackedSwitchKey {
    static constexpr std::uint32_t lwe_dimension =
        TFHEParam::k * TFHEParam::n;
    static_assert(lwe_dimension <= CKKSParam::n / 2,
                  "TFHE LWE secret does not fit in the real CKKS slots");

    CKKSCiphertext<CKKSParam, LogQ, LogDelta> encrypted_secret_slots{};
};

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta>
inline void TFHEToCKKSSlotPackedKeyGen(
    TFHEToCKKSSlotPackedSwitchKey<CKKSParam, TFHEParam, LogQ, LogDelta>
        &switch_key,
    const Key<CKKSParam> &ckks_key, const Key<TFHEParam> &tfhe_key,
    CKKSNoise noise = {CKKSParam::α, 0})
{
    CKKSSlotVector<CKKSParam> secret{};
    secret.fill({0.0, 0.0});
    for (std::uint32_t i = 0;
         i < TFHEToCKKSSlotPackedSwitchKey<CKKSParam, TFHEParam, LogQ,
                                           LogDelta>::lwe_dimension;
         i++)
        secret[i] = {
            ckks_scheme_switch_detail::tfheSecretCoefficientToDouble<TFHEParam>(
                tfhe_key[i]),
            0.0};
    ckksSlotEncrypt<CKKSParam, LogQ, LogDelta>(
        switch_key.encrypted_secret_slots, secret, ckks_key, noise);
}

// Homomorphically evaluates b-a*s into consecutive real CKKS slots.  This is
// the dense rectangular partial-decryption form used by OpenFHE.  The phase is
// unwrapped; apply CKKSEvalModBoundedCosNormalized (or a dense bootstrap) to
// reduce it modulo one.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          class InputGaloisKey, class OutputGaloisKey>
inline void TFHEToCKKSSlotPackedPhase(
    CKKSPlainMulResult<CKKSParam, LogQ, LogDelta, PlainLogDelta> &out,
    const std::vector<TLWE<TFHEParam>> &in,
    const TFHEToCKKSSlotPackedSwitchKey<CKKSParam, TFHEParam, LogQ, LogDelta>
        &switch_key,
    int bsgs_step, const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key, double phase_scale = 1.0)
{
    static_assert(PlainLogDelta > 0);
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t lwe_dimension =
        TFHEToCKKSSlotPackedSwitchKey<CKKSParam, TFHEParam, LogQ,
                                      LogDelta>::lwe_dimension;
    constexpr int slot_count = static_cast<int>(CKKSParam::n / 2);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    if (in.empty() || in.size() > static_cast<std::size_t>(slot_count))
        throw std::out_of_range("invalid number of dense packed TFHE ciphertexts");
    if (bsgs_step <= 0 || bsgs_step > slot_count)
        throw std::invalid_argument("invalid dense TFHE-to-CKKS BSGS step");
    if (!std::isfinite(phase_scale) || phase_scale == 0.0)
        throw std::invalid_argument("invalid dense packed TFHE phase scale");

    std::vector<int> offsets;
    std::vector<CKKSSlotVector<CKKSParam>> diagonals;
    for (std::uint32_t slot = 0; slot < in.size(); slot++) {
        for (std::uint32_t secret_index = 0; secret_index < lwe_dimension;
             secret_index++) {
            const int offset =
                (static_cast<int>(secret_index) - static_cast<int>(slot) +
                 slot_count) %
                slot_count;
            std::size_t diagonal = 0;
            while (diagonal < offsets.size() && offsets[diagonal] != offset)
                diagonal++;
            if (diagonal == offsets.size()) {
                offsets.push_back(offset);
                diagonals.emplace_back();
                diagonals.back().fill({0.0, 0.0});
            }
            diagonals[diagonal][slot] = {
                -phase_scale *
                    ckks_scheme_switch_detail::
                        tfheTorusToUnitInterval<TFHEParam>(
                            in[slot][secret_index]),
                0.0};
        }
    }

    CKKSLinearTransformPlan<CKKSParam, LogQ, LogDelta, PlainLogDelta> plan;
    CKKSBuildLinearTransformBSGSPlan<CKKSParam, LogQ, LogDelta,
                                     PlainLogDelta>(
        plan, diagonals, offsets, bsgs_step);
    CKKSLinearTransformBSGS<CKKSParam, LogQ, LogDelta, PlainLogDelta>(
        out, switch_key.encrypted_secret_slots, plan, input_galois_key,
        output_galois_key);

    auto body = std::make_unique<CKKSSlotVector<CKKSParam>>();
    body->fill({0.0, 0.0});
    for (std::uint32_t slot = 0; slot < in.size(); slot++)
        (*body)[slot] = {
            phase_scale *
                ckks_scheme_switch_detail::tfheTorusToUnitInterval<TFHEParam>(
                    in[slot][lwe_dimension]),
            0.0};
    auto encoded_body = std::make_unique<Polynomial<CKKSParam>>();
    ckksSlotEncode<CKKSParam, out_log_q, LogDelta>(*encoded_body, *body);
    for (std::uint32_t i = 0; i < CKKSParam::n; i++)
        out.ct[CKKSParam::k][i] =
            ckks_detail::reduceToLevel<CKKSParam, out_log_q>(
                out.ct[CKKSParam::k][i] + (*encoded_body)[i]);
}

// Rotation-key schedule for TFHEToCKKSSlotPackedPhase.  It derives the BSGS
// baby and giant rotations from the dense rectangular mask layout without
// materializing the public-mask plaintexts.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta>
inline void TFHEToCKKSSlotPackedPhaseRotationKeyIndices(
    CKKSRotationKeyIndexSet<CKKSParam> &input_keys,
    CKKSRotationKeyIndexSet<CKKSParam> &output_keys,
    std::uint32_t num_ciphertexts, int bsgs_step)
{
    static_assert(PlainLogDelta > 0);
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t lwe_dimension =
        TFHEToCKKSSlotPackedSwitchKey<CKKSParam, TFHEParam, LogQ,
                                      LogDelta>::lwe_dimension;
    constexpr int slot_count = static_cast<int>(CKKSParam::n / 2);
    if (num_ciphertexts == 0 || num_ciphertexts > static_cast<std::uint32_t>(slot_count))
        throw std::out_of_range("invalid number of dense packed TFHE ciphertexts");
    if (bsgs_step <= 0 || bsgs_step > slot_count)
        throw std::invalid_argument("invalid dense TFHE-to-CKKS BSGS step");

    CKKSClearRotationKeyIndexSet<CKKSParam>(input_keys);
    CKKSClearRotationKeyIndexSet<CKKSParam>(output_keys);
    for (std::uint32_t slot = 0; slot < num_ciphertexts; slot++)
        for (std::uint32_t secret_index = 0; secret_index < lwe_dimension;
             secret_index++) {
            const int rotation =
                (static_cast<int>(secret_index) - static_cast<int>(slot) +
                 slot_count) %
                slot_count;
            const int baby = rotation % bsgs_step;
            const int giant = rotation - baby;
            if (baby != 0)
                CKKSMarkRotationPowerKeyIndices<CKKSParam>(input_keys, baby);
            if (giant != 0)
                CKKSMarkRotationPowerKeyIndices<CKKSParam>(output_keys,
                                                            giant);
        }
}

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta>
inline void TFHEToCKKSSlotPackedPhaseSparseGaloisKeyGen(
    CKKSSparseGaloisKey<CKKSParam, LogQ> &input_galois_key,
    CKKSSparseGaloisKey<CKKSParam, LogQ - PlainLogDelta> &output_galois_key,
    const Key<CKKSParam> &key, std::uint32_t num_ciphertexts, int bsgs_step,
    CKKSNoise noise = {CKKSParam::α, 0})
{
    CKKSRotationKeyIndexSet<CKKSParam> input_indices{};
    CKKSRotationKeyIndexSet<CKKSParam> output_indices{};
    TFHEToCKKSSlotPackedPhaseRotationKeyIndices<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta>(
        input_indices, output_indices, num_ciphertexts, bsgs_step);
    CKKSSparseGaloisKeyGen<CKKSParam, LogQ>(input_galois_key, key,
                                            input_indices, noise);
    CKKSSparseGaloisKeyGen<CKKSParam, LogQ - PlainLogDelta>(
        output_galois_key, key, output_indices, noise);
}

// Homomorphically evaluates b-a*s for a batch of TFHE TLWEs and returns one
// unwrapped normalized phase per secret-sized CKKS slot block.  The i-th
// phase is at slot i*lwe_dimension; the other slots are unspecified.  As in
// TFHEToCKKSPhase, the caller must follow it with a periodic
// EvalMod/bootstrapping stage to obtain canonical messages.  The secret is
// replicated by block at key generation, so all rotations are in
// [0,lwe_dimension), avoiding a near-full-ring rotation for every packed LWE.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          class InputGaloisKey, class OutputGaloisKey>
inline void TFHEToCKKSPackedPhase(
    CKKSPlainMulResult<CKKSParam, LogQ, LogDelta, PlainLogDelta> &out,
    const std::vector<TLWE<TFHEParam>> &in,
    const TFHEToCKKSPackedSwitchKey<CKKSParam, TFHEParam, LogQ, LogDelta>
        &switch_key,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key, double phase_scale = 1.0)
{
    static_assert(PlainLogDelta > 0);
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t lwe_dimension =
        TFHEToCKKSPackedSwitchKey<CKKSParam, TFHEParam, LogQ,
                                  LogDelta>::lwe_dimension;
    constexpr std::uint32_t max_ciphertexts =
        TFHEToCKKSPackedSwitchKey<CKKSParam, TFHEParam, LogQ,
                                  LogDelta>::max_ciphertexts;
    constexpr int slot_count = static_cast<int>(CKKSParam::n / 2);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;

    if (in.empty() || in.size() > max_ciphertexts)
        throw std::out_of_range("invalid number of packed TFHE ciphertexts");
    if (!std::isfinite(phase_scale) || phase_scale == 0.0)
        throw std::invalid_argument("invalid packed TFHE phase scale");

    std::vector<int> offsets;
    std::vector<CKKSSlotVector<CKKSParam>> diagonals;
    for (std::uint32_t ciphertext_index = 0;
         ciphertext_index < in.size(); ciphertext_index++) {
        const std::uint32_t slot = ciphertext_index * lwe_dimension;
        for (std::uint32_t secret_index = 0; secret_index < lwe_dimension;
             secret_index++) {
            const int offset = static_cast<int>(secret_index);
            std::size_t diagonal = 0;
            while (diagonal < offsets.size() && offsets[diagonal] != offset)
                diagonal++;
            if (diagonal == offsets.size()) {
                offsets.push_back(offset);
                diagonals.emplace_back();
                diagonals.back().fill({0.0, 0.0});
            }
            diagonals[diagonal][slot] = {
                -phase_scale *
                    ckks_scheme_switch_detail::
                        tfheTorusToUnitInterval<TFHEParam>(
                            in[ciphertext_index][secret_index]),
                0.0};
        }
    }

    CKKSLinearTransformPlan<CKKSParam, LogQ, LogDelta, PlainLogDelta> plan;
    CKKSBuildLinearTransformBSGSPlan<CKKSParam, LogQ, LogDelta,
                                     PlainLogDelta>(
        plan, diagonals, offsets, slot_count);
    CKKSLinearTransformBSGS<CKKSParam, LogQ, LogDelta, PlainLogDelta>(
        out, switch_key.encrypted_secret_slots, plan, input_galois_key,
        output_galois_key);

    auto body = std::make_unique<CKKSSlotVector<CKKSParam>>();
    body->fill({0.0, 0.0});
    for (std::uint32_t ciphertext_index = 0;
         ciphertext_index < in.size(); ciphertext_index++) {
        const std::uint32_t slot = ciphertext_index * lwe_dimension;
        (*body)[slot] = {
            phase_scale *
                ckks_scheme_switch_detail::tfheTorusToUnitInterval<TFHEParam>(
                    in[ciphertext_index][lwe_dimension]),
            0.0};
    }
    auto encoded_body = std::make_unique<Polynomial<CKKSParam>>();
    ckksSlotEncode<CKKSParam, out_log_q, LogDelta>(*encoded_body, *body);
    for (std::uint32_t i = 0; i < CKKSParam::n; i++)
        out.ct[CKKSParam::k][i] =
            ckks_detail::reduceToLevel<CKKSParam, out_log_q>(
                out.ct[CKKSParam::k][i] + (*encoded_body)[i]);
}

// The block-replicated packed reverse transform uses only rotations by TFHE
// secret-coordinate offsets.  Expose their sparse binary-key index sets so a
// switching context never needs a complete Galois key at either level.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta>
inline void TFHEToCKKSPackedPhaseRotationKeyIndices(
    CKKSRotationKeyIndexSet<CKKSParam> &input_keys,
    CKKSRotationKeyIndexSet<CKKSParam> &output_keys)
{
    static_assert(PlainLogDelta > 0);
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t lwe_dimension =
        TFHEToCKKSPackedSwitchKey<CKKSParam, TFHEParam, LogQ,
                                  LogDelta>::lwe_dimension;
    CKKSClearRotationKeyIndexSet<CKKSParam>(input_keys);
    CKKSClearRotationKeyIndexSet<CKKSParam>(output_keys);
    for (std::uint32_t offset = 1; offset < lwe_dimension; offset++)
        CKKSMarkRotationPowerKeyIndices<CKKSParam>(
            input_keys, static_cast<int>(offset));
}

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta>
inline void TFHEToCKKSPackedPhaseSparseGaloisKeyGen(
    CKKSSparseGaloisKey<CKKSParam, LogQ> &input_galois_key,
    CKKSSparseGaloisKey<CKKSParam, LogQ - PlainLogDelta> &output_galois_key,
    const Key<CKKSParam> &key, CKKSNoise noise = {CKKSParam::α, 0})
{
    CKKSRotationKeyIndexSet<CKKSParam> input_indices{};
    CKKSRotationKeyIndexSet<CKKSParam> output_indices{};
    TFHEToCKKSPackedPhaseRotationKeyIndices<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta>(input_indices,
                                                              output_indices);
    CKKSSparseGaloisKeyGen<CKKSParam, LogQ>(input_galois_key, key,
                                            input_indices, noise);
    CKKSSparseGaloisKeyGen<CKKSParam, LogQ - PlainLogDelta>(
        output_galois_key, key, output_indices, noise);
}

// Evaluation material for the complete packed TFHE -> CKKS conversion. EvalMod
// sees message_ratio*phase/(phase_bound*message_ratio), so phase_bound must
// bound |b/Q-a*s/Q| for every input ciphertext. For ternary secrets a
// conservative choice is lwe_dimension+1.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
struct TFHEToCKKSPackedEvalModKey {
    static_assert(LogQ > PlainLogDelta);
    static constexpr std::uint32_t phase_log_q = LogQ - PlainLogDelta;
    using SwitchKey =
        TFHEToCKKSPackedSwitchKey<CKKSParam, TFHEParam, LogQ, LogDelta>;
    using EvalModKey = CKKSEvalModBoundedCosRelinKeys<
        CKKSParam, phase_log_q, LogDelta, CoeffLogDelta, Degree, DoubleAngle,
        InvDegree>;

    SwitchKey switch_key{};
    EvalModKey evalmod_key{};
    CKKSBoundedCosEvalModPolynomial polynomial{};
    double phase_normalizer = 0.0;
    double phase_scale = 0.0;
};

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
inline void TFHEToCKKSPackedEvalModKeyGen(
    TFHEToCKKSPackedEvalModKey<CKKSParam, TFHEParam, LogQ, LogDelta,
                                 PlainLogDelta, CoeffLogDelta, Degree,
                                 DoubleAngle, InvDegree> &key,
    const Key<CKKSParam> &ckks_key, const Key<TFHEParam> &tfhe_key,
    std::uint32_t phase_bound, std::uint32_t log_message_ratio = 0,
    CKKSNoise noise = {CKKSParam::α, 0})
{
    if (phase_bound == 0)
        throw std::invalid_argument("TFHE phase bound must be positive");
    TFHEToCKKSPackedKeyGen<CKKSParam, TFHEParam, LogQ, LogDelta>(
        key.switch_key, ckks_key, tfhe_key, noise);
    CKKSEvalModBoundedCosKeyGen<CKKSParam, LogQ - PlainLogDelta, LogDelta,
                                CoeffLogDelta, Degree, DoubleAngle, InvDegree>(
        key.evalmod_key, ckks_key, noise);
    key.polynomial = CKKSBuildBoundedCosEvalModPolynomial(
        phase_bound, static_cast<std::uint32_t>(Degree), log_message_ratio,
        DoubleAngle, 1.0, InvDegree != 0);
    key.phase_normalizer = static_cast<double>(phase_bound) *
                           key.polynomial.message_ratio *
                           key.polynomial.q_diff;
    // EvalMod takes (message_ratio * phase) / phase_normalizer.  Encoding the
    // phase divided by phase_bound therefore gives exactly that normalized
    // input while retaining the original message scale after reduction.
    key.phase_scale = 1.0 / static_cast<double>(phase_bound);
}

// OpenFHE-style packed reverse switch: encrypted partial decryption followed
// by bounded cosine EvalMod.  The result contains canonical approximate CKKS
// messages at the block-leading slots described by TFHEToCKKSPackedPhase.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0,
          class InputGaloisKey, class OutputGaloisKey>
inline void TFHEToCKKSPackedEvalMod(
    CKKSEvalModBoundedCosResult<CKKSParam, LogQ - PlainLogDelta, LogDelta,
                                CoeffLogDelta, Degree, DoubleAngle,
                                InvDegree> &out,
    const std::vector<TLWE<TFHEParam>> &in,
    const TFHEToCKKSPackedEvalModKey<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta, CoeffLogDelta,
        Degree, DoubleAngle, InvDegree> &key,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key)
{
    if (!(key.phase_normalizer > 0.0) ||
        !std::isfinite(key.phase_normalizer) || !(key.phase_scale > 0.0) ||
        !std::isfinite(key.phase_scale))
        throw std::invalid_argument("invalid TFHE-to-CKKS EvalMod key");
    auto phase = std::make_unique<
        CKKSPlainMulResult<CKKSParam, LogQ, LogDelta, PlainLogDelta>>();
    TFHEToCKKSPackedPhase<CKKSParam, TFHEParam, LogQ, LogDelta,
                          PlainLogDelta>(
        *phase, in, key.switch_key, input_galois_key, output_galois_key,
        key.phase_scale);
    CKKSEvalModBoundedCosNormalized<
        CKKSParam, LogQ - PlainLogDelta, LogDelta, CoeffLogDelta, Degree,
        DoubleAngle, InvDegree>(out, *phase, key.polynomial, key.evalmod_key);
}

// Complete dense-slot counterpart of TFHEToCKKSPackedEvalModKey.  It keeps a
// single encrypted TFHE secret and therefore supports OpenFHE-style
// consecutive-slot packing through TFHEToCKKSSlotPackedPhase.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
struct TFHEToCKKSSlotPackedEvalModKey {
    static_assert(LogQ > PlainLogDelta);
    static constexpr std::uint32_t phase_log_q = LogQ - PlainLogDelta;
    using SwitchKey =
        TFHEToCKKSSlotPackedSwitchKey<CKKSParam, TFHEParam, LogQ, LogDelta>;
    using EvalModKey = CKKSEvalModBoundedCosRelinKeys<
        CKKSParam, phase_log_q, LogDelta, CoeffLogDelta, Degree, DoubleAngle,
        InvDegree>;

    SwitchKey switch_key{};
    EvalModKey evalmod_key{};
    CKKSBoundedCosEvalModPolynomial polynomial{};
    double phase_normalizer = 0.0;
    double phase_scale = 0.0;
};

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
inline void TFHEToCKKSSlotPackedEvalModKeyGen(
    TFHEToCKKSSlotPackedEvalModKey<CKKSParam, TFHEParam, LogQ, LogDelta,
                                     PlainLogDelta, CoeffLogDelta, Degree,
                                     DoubleAngle, InvDegree> &key,
    const Key<CKKSParam> &ckks_key, const Key<TFHEParam> &tfhe_key,
    std::uint32_t phase_bound, std::uint32_t log_message_ratio = 0,
    CKKSNoise noise = {CKKSParam::α, 0})
{
    if (phase_bound == 0)
        throw std::invalid_argument("TFHE phase bound must be positive");
    TFHEToCKKSSlotPackedKeyGen<CKKSParam, TFHEParam, LogQ, LogDelta>(
        key.switch_key, ckks_key, tfhe_key, noise);
    CKKSEvalModBoundedCosKeyGen<CKKSParam, LogQ - PlainLogDelta, LogDelta,
                                CoeffLogDelta, Degree, DoubleAngle, InvDegree>(
        key.evalmod_key, ckks_key, noise);
    key.polynomial = CKKSBuildBoundedCosEvalModPolynomial(
        phase_bound, static_cast<std::uint32_t>(Degree), log_message_ratio,
        DoubleAngle, 1.0, InvDegree != 0);
    key.phase_normalizer = static_cast<double>(phase_bound) *
                           key.polynomial.message_ratio *
                           key.polynomial.q_diff;
    key.phase_scale = 1.0 / static_cast<double>(phase_bound);
}

template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0,
          class InputGaloisKey, class OutputGaloisKey>
inline void TFHEToCKKSSlotPackedEvalMod(
    CKKSEvalModBoundedCosResult<CKKSParam, LogQ - PlainLogDelta, LogDelta,
                                CoeffLogDelta, Degree, DoubleAngle,
                                InvDegree> &out,
    const std::vector<TLWE<TFHEParam>> &in,
    const TFHEToCKKSSlotPackedEvalModKey<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta, CoeffLogDelta,
        Degree, DoubleAngle, InvDegree> &key,
    int bsgs_step, const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key)
{
    if (!(key.phase_normalizer > 0.0) ||
        !std::isfinite(key.phase_normalizer) || !(key.phase_scale > 0.0) ||
        !std::isfinite(key.phase_scale))
        throw std::invalid_argument("invalid dense TFHE-to-CKKS EvalMod key");
    auto phase = std::make_unique<
        CKKSPlainMulResult<CKKSParam, LogQ, LogDelta, PlainLogDelta>>();
    TFHEToCKKSSlotPackedPhase<CKKSParam, TFHEParam, LogQ, LogDelta,
                              PlainLogDelta>(
        *phase, in, key.switch_key, bsgs_step, input_galois_key,
        output_galois_key, key.phase_scale);
    CKKSEvalModBoundedCosNormalized<
        CKKSParam, LogQ - PlainLogDelta, LogDelta, CoeffLogDelta, Degree,
        DoubleAngle, InvDegree>(out, *phase, key.polynomial, key.evalmod_key);
}

// Complete dense reverse switch with OpenFHE-style post-processing.  EvalMod
// yields the canonical fractional phase; this final plaintext multiplication
// and addition maps it to output_scale * phase + output_bias.  It consumes
// PostPlainLogDelta CKKS level bits, just as OpenFHE's final post-scale does.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::uint32_t PostPlainLogDelta,
          std::size_t InvDegree = 0, class InputGaloisKey,
          class OutputGaloisKey>
inline void TFHEToCKKSSlotPackedEvalModAffine(
    CKKSPlainMulResult<
        CKKSParam,
        CKKSEvalModBoundedCosTraits<CKKSParam, LogQ - PlainLogDelta,
                                    LogDelta, CoeffLogDelta, Degree,
                                    DoubleAngle, InvDegree>::log_q,
        LogDelta, PostPlainLogDelta> &out,
    const std::vector<TLWE<TFHEParam>> &in,
    const TFHEToCKKSSlotPackedEvalModKey<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta, CoeffLogDelta,
        Degree, DoubleAngle, InvDegree> &key,
    int bsgs_step, const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key, double output_scale,
    double output_bias)
{
    using EvalOut = CKKSEvalModBoundedCosResult<
        CKKSParam, LogQ - PlainLogDelta, LogDelta, CoeffLogDelta, Degree,
        DoubleAngle, InvDegree>;
    static_assert(PostPlainLogDelta > 0);
    static_assert(PostPlainLogDelta < EvalOut::log_q,
                  "post-processing level drop exceeds EvalMod output level");
    if (!std::isfinite(output_scale) || !std::isfinite(output_bias))
        throw std::invalid_argument("invalid TFHE-to-CKKS output affine map");

    auto periodic = std::make_unique<EvalOut>();
    TFHEToCKKSSlotPackedEvalMod<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta, CoeffLogDelta,
        Degree, DoubleAngle, InvDegree>(*periodic, in, key, bsgs_step,
                                        input_galois_key, output_galois_key);
    CKKSPlainMulByReal<CKKSParam, EvalOut::log_q, LogDelta,
                       PostPlainLogDelta>(out, *periodic, output_scale);
    CKKSAddPlainRealInPlace<
        CKKSParam, CKKSPlainMulResult<CKKSParam, EvalOut::log_q, LogDelta,
                                      PostPlainLogDelta>::log_q,
        LogDelta>(out, output_bias);
}

// OpenFHE's EvalFHEWtoCKKS post-processing convention.  Its periodic stage
// returns a fractional phase, then maps it with a modulus-dependent scale and
// (for a requested [pmin, pmax] interval) a public bias.  Keep this wrapper
// beside the general affine API so callers can reproduce OpenFHE's `p`,
// `pmin`, and `pmax` arguments without re-deriving that convention.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::uint32_t PostPlainLogDelta,
          std::size_t InvDegree = 0, class InputGaloisKey,
          class OutputGaloisKey>
inline void TFHEToCKKSSlotPackedEvalModOpenFHE(
    CKKSPlainMulResult<
        CKKSParam,
        CKKSEvalModBoundedCosTraits<CKKSParam, LogQ - PlainLogDelta,
                                    LogDelta, CoeffLogDelta, Degree,
                                    DoubleAngle, InvDegree>::log_q,
        LogDelta, PostPlainLogDelta> &out,
    const std::vector<TLWE<TFHEParam>> &in,
    const TFHEToCKKSSlotPackedEvalModKey<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta, CoeffLogDelta,
        Degree, DoubleAngle, InvDegree> &key,
    int bsgs_step, const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key, std::uint32_t plaintext_modulus,
    double plaintext_min = 0.0, double plaintext_max = 0.0)
{
    if (plaintext_modulus == 0)
        throw std::invalid_argument("OpenFHE plaintext modulus must be positive");
    if (!std::isfinite(plaintext_min) || !std::isfinite(plaintext_max))
        throw std::invalid_argument("invalid OpenFHE plaintext range");

    constexpr double two_pi = 6.283185307179586476925286766559005768;
    double output_scale = plaintext_modulus <= 4
                              ? two_pi
                              : static_cast<double>(plaintext_modulus);
    double output_bias = 0.0;
    if (plaintext_min != 0.0) {
        if (!(plaintext_max > plaintext_min))
            throw std::invalid_argument("invalid OpenFHE plaintext range");
        const double half_range = (plaintext_max - plaintext_min) / 4.0;
        output_scale *= half_range;
        output_bias = half_range;
    }
    TFHEToCKKSSlotPackedEvalModAffine<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta, CoeffLogDelta,
        Degree, DoubleAngle, PostPlainLogDelta, InvDegree>(
        out, in, key, bsgs_step, input_galois_key, output_galois_key,
        output_scale, output_bias);
}

// OpenFHE's EvalFHEWtoCKKS can return a sparsely packed ciphertext: a batch
// occupying the first `filled_slots` slots is repeated until `output_slots`.
// Keep that operation separate from partial decryption so callers that want
// full packing do not pay for its rotations.  The ratios must be powers of two
// because the implementation uses the same doubling tree as OpenFHE.
template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          class GaloisKey>
inline void CKKSSchemeSwitchReplicateFirstSlotsInPlace(
    CKKSCiphertext<P, LogQ, LogDelta> &ct, std::uint32_t filled_slots,
    std::uint32_t output_slots, const GaloisKey &galois_key)
{
    constexpr std::uint32_t slot_count = P::n / 2;
    if (filled_slots == 0 || filled_slots > output_slots ||
        output_slots > slot_count ||
        (output_slots % filled_slots) != 0)
        throw std::out_of_range("invalid scheme-switch output slot count");
    const std::uint32_t repetitions = output_slots / filled_slots;
    if ((repetitions & (repetitions - 1)) != 0)
        throw std::invalid_argument(
            "scheme-switch output replication ratio must be a power of two");

    for (std::uint32_t rotation = filled_slots; rotation < output_slots;
         rotation <<= 1) {
        auto rotated = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
        CKKSRotateSlots<P, LogQ>(rotated->ct, ct.ct,
                                 static_cast<int>(rotation), galois_key);
        CKKSAddInPlace<P, LogQ, LogDelta>(ct, *rotated);
    }
}

// Sparse binary rotation-key schedule for
// CKKSSchemeSwitchReplicateFirstSlotsInPlace.  It is intentionally generated
// at the final EvalMod level, rather than at the partial-decryption level.
template <class P>
inline void CKKSSchemeSwitchReplicationRotationKeyIndices(
    CKKSRotationKeyIndexSet<P> &indices, std::uint32_t filled_slots,
    std::uint32_t output_slots)
{
    constexpr std::uint32_t slot_count = P::n / 2;
    if (filled_slots == 0 || filled_slots > output_slots ||
        output_slots > slot_count ||
        (output_slots % filled_slots) != 0)
        throw std::out_of_range("invalid scheme-switch output slot count");
    const std::uint32_t repetitions = output_slots / filled_slots;
    if ((repetitions & (repetitions - 1)) != 0)
        throw std::invalid_argument(
            "scheme-switch output replication ratio must be a power of two");

    CKKSClearRotationKeyIndexSet<P>(indices);
    for (std::uint32_t rotation = filled_slots; rotation < output_slots;
         rotation <<= 1)
        CKKSMarkRotationPowerKeyIndices<P>(indices,
                                           static_cast<int>(rotation));
}

template <class P, std::uint32_t LogQ>
inline void CKKSSchemeSwitchReplicationSparseGaloisKeyGen(
    CKKSSparseGaloisKey<P, LogQ> &galois_key, const Key<P> &key,
    std::uint32_t filled_slots, std::uint32_t output_slots,
    CKKSNoise noise = {P::α, 0})
{
    CKKSRotationKeyIndexSet<P> indices{};
    CKKSSchemeSwitchReplicationRotationKeyIndices<P>(
        indices, filled_slots, output_slots);
    CKKSSparseGaloisKeyGen<P, LogQ>(galois_key, key, indices, noise);
}

// Complete dense TFHE -> CKKS conversion followed by OpenFHE-style sparse
// output packing.  `in.size()` values are repeated over `output_slots`; the
// caller may then use the ciphertext with the same logical slot count as the
// originating CKKS application.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0,
          class InputGaloisKey, class OutputGaloisKey,
          class ReplicationGaloisKey>
inline void TFHEToCKKSSlotPackedEvalModToSlots(
    CKKSEvalModBoundedCosResult<CKKSParam, LogQ - PlainLogDelta, LogDelta,
                                CoeffLogDelta, Degree, DoubleAngle,
                                InvDegree> &out,
    const std::vector<TLWE<TFHEParam>> &in,
    const TFHEToCKKSSlotPackedEvalModKey<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta, CoeffLogDelta,
        Degree, DoubleAngle, InvDegree> &key,
    int bsgs_step, std::uint32_t output_slots,
    const InputGaloisKey &input_galois_key,
    const OutputGaloisKey &output_galois_key,
    const ReplicationGaloisKey &replication_galois_key)
{
    TFHEToCKKSSlotPackedEvalMod<
        CKKSParam, TFHEParam, LogQ, LogDelta, PlainLogDelta, CoeffLogDelta,
        Degree, DoubleAngle, InvDegree>(out, in, key, bsgs_step,
                                        input_galois_key, output_galois_key);
    CKKSSchemeSwitchReplicateFirstSlotsInPlace<
        CKKSParam,
        CKKSEvalModBoundedCosTraits<CKKSParam, LogQ - PlainLogDelta,
                                    LogDelta, CoeffLogDelta, Degree,
                                    DoubleAngle, InvDegree>::log_q,
        LogDelta>(
        out, static_cast<std::uint32_t>(in.size()), output_slots,
        replication_galois_key);
}

// Homomorphically constructs an encryption of the *unwrapped* normalized LWE
// phase b/Q - a·s/Q in coefficient zero.  Values differing by an integer are
// equivalent after the periodic modular reduction used in the next stage.
// Other CKKS coefficients are intentionally unspecified: this primitive is a
// building block for the packed linear transform, not the final packed API.
template <class CKKSParam, class TFHEParam, std::uint32_t LogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta>
inline void TFHEToCKKSPhase(
    CKKSPlainMulResult<CKKSParam, LogQ, LogDelta, PlainLogDelta> &out,
    const TLWE<TFHEParam> &in,
    const TFHEToCKKSSwitchKey<CKKSParam, TFHEParam, LogQ, LogDelta>
        &switch_key)
{
    static_assert(PlainLogDelta > 0);
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t lwe_dimension =
        TFHEToCKKSSwitchKey<CKKSParam, TFHEParam, LogQ,
                             LogDelta>::lwe_dimension;

    // In Z[X]/(X^N+1), coefficient zero of this product is -a·s/Q:
    // p[0] = -a[0]/Q and p[N-i] = a[i]/Q for i > 0.
    auto reverse_mask = std::make_unique<std::array<double, CKKSParam::n>>();
    reverse_mask->fill(0.0);
    (*reverse_mask)[0] =
        -ckks_scheme_switch_detail::tfheTorusToUnitInterval<TFHEParam>(in[0]);
    for (std::uint32_t i = 1; i < lwe_dimension; i++)
        (*reverse_mask)[CKKSParam::n - i] =
            ckks_scheme_switch_detail::tfheTorusToUnitInterval<TFHEParam>(
                in[i]);

    auto plain = std::make_unique<
        CKKSPlaintext<CKKSParam, LogQ, PlainLogDelta>>();
    ckksEncodePolynomial<CKKSParam, LogQ, PlainLogDelta>(plain->poly,
                                                           *reverse_mask);
    CKKSPlainMulRescale<CKKSParam, LogQ, LogDelta, PlainLogDelta>(
        out, switch_key.encrypted_secret, *plain);

    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    out.ct[CKKSParam::k][0] = ckks_detail::reduceToLevel<CKKSParam, out_log_q>(
        out.ct[CKKSParam::k][0] +
        ckksEncodeCoeff<CKKSParam, out_log_q, LogDelta>(
            ckks_scheme_switch_detail::tfheTorusToUnitInterval<TFHEParam>(
                in[lwe_dimension])));
}

}  // namespace TFHEpp
