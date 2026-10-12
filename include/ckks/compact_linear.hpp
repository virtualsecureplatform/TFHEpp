#pragma once

#include "ckks.hpp"
#include <bit>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// An explicit, lossless format for low-scale plaintext matrices. Evaluation
// expands one giant group and uses the existing BSGS kernel, including its
// fused accumulation and rounding. Legacy plans and archives are unchanged.
template <class P, std::uint32_t Q, std::uint32_t Delta, std::uint32_t PlainDelta>
struct CKKSCompactLinearTransformPlan {
    static_assert(is_multilimb_uint_v<typename P::T>);
    struct Term {
        int baby_step = 0;
        std::array<std::int64_t, P::n> coefficients{};
        template <class Archive> void serialize(Archive &ar) { ar(baby_step, coefficients); }
    };
    struct Group {
        int giant_step = 0;
        std::vector<Term> terms;
        template <class Archive> void serialize(Archive &ar) { ar(giant_step, terms); }
    };
    int k_step = 0;
    std::vector<Group> groups;

    static std::int64_t pack(typename P::T value)
    {
        const auto centered = ckks_detail::centeredLevelToTorus<P, Q>(value);
        const auto digit = std::bit_cast<std::int64_t>(centered.limb[0]);
        if (P::T::from_signed_i64(digit) != centered)
            throw std::invalid_argument("CKKS plaintext coefficient exceeds compact int64 range");
        return digit;
    }
    template <class Archive> void serialize(Archive &ar)
    {
        const std::array<std::uint64_t, 7> expected{
            0x434b4b53434c5031ULL, P::n, P::k, Q, Delta, PlainDelta, 64};
        auto metadata = expected;
        ar(metadata);
        if (metadata != expected) throw std::runtime_error("compact CKKS plan metadata mismatch");
        ar(k_step, groups);
    }
};

template <class P, std::uint32_t Q, std::uint32_t Delta, std::uint32_t PlainDelta>
inline void CKKSLinearTransformBSGS(
    CKKSPlainMulResult<P, Q, Delta, PlainDelta> &res,
    const CKKSCiphertext<P, Q, Delta> &ct,
    const CKKSCompactLinearTransformPlan<P, Q, Delta, PlainDelta> &plan,
    const auto &input_gk, const auto &output_gk)
{
    if (plan.k_step <= 0 || plan.k_step > static_cast<int>(P::n / 2) || plan.groups.empty())
        throw std::invalid_argument("invalid compact CKKS transform schedule");
    std::vector<bool> used(plan.k_step, false);
    for (const auto &group : plan.groups) {
        if (group.giant_step < 0 || group.giant_step * std::int64_t{plan.k_step} >= P::n / 2 || group.terms.empty())
            throw std::invalid_argument("invalid compact CKKS giant step");
        for (const auto &term : group.terms) {
            if (term.baby_step < 0 || term.baby_step >= plan.k_step)
                throw std::invalid_argument("invalid compact CKKS baby step");
            used[term.baby_step] = true;
        }
    }
    auto expanded = std::make_unique<CKKSLinearTransformGiantStep<P, Q, PlainDelta>>();
    CKKSLinearTransformBSGSWithGroupProvider<P, Q, Delta, PlainDelta>(
        res, ct, plan.k_step, used, plan.groups.size(),
        [&](std::size_t index) -> const auto & {
            const auto &group = plan.groups[index];
            expanded->giant_step = group.giant_step;
            expanded->terms.resize(group.terms.size());
            for (std::size_t i = 0; i < group.terms.size(); i++) {
                expanded->terms[i].baby_step = group.terms[i].baby_step;
                for (std::size_t n = 0; n < P::n; n++)
                    expanded->terms[i].plain.poly[n] = ckks_detail::reduceToLevel<P, Q>(
                        P::T::from_signed_i64(group.terms[i].coefficients[n]));
            }
            return *expanded;
        }, input_gk, output_gk);
}

template <class P, std::uint32_t Q, std::uint32_t Delta, std::uint32_t PlainDelta>
struct CKKSCompactRealLinearTransformPlan {
    bool has_direct = false, has_conjugate = false;
    CKKSCompactLinearTransformPlan<P, Q, Delta, PlainDelta> direct, conjugate;
};

// Plaintext diagonals depend on PlainDelta, not the input ciphertext's scale.
// A different modulus also preserves their signed int64 coefficients when
// both levels have at least 64 bits. Reuse the encoded coefficients directly.
// The destination archive records its new input scale; source metadata and
// payload remain unchanged. This does not change ciphertext or key levels.
template <class P, std::uint32_t NewQ, std::uint32_t NewDelta,
          std::uint32_t OldQ, std::uint32_t OldDelta, std::uint32_t PlainDelta>
inline void CKKSRebindCompactLinearInputScale(
    CKKSCompactLinearTransformPlan<P, NewQ, NewDelta, PlainDelta> &result,
    const CKKSCompactLinearTransformPlan<P, OldQ, OldDelta, PlainDelta> &source)
{
    static_assert(NewDelta < NewQ);
    static_assert(NewQ == OldQ || (NewQ >= 64 && OldQ >= 64));
    if constexpr (NewDelta == OldDelta && NewQ == OldQ) {
        result = source;
    }
    else {
        result.k_step = source.k_step;
        result.groups.resize(source.groups.size());
        for (std::size_t g = 0; g < source.groups.size(); ++g) {
            auto &out = result.groups[g];
            const auto &in = source.groups[g];
            out.giant_step = in.giant_step;
            out.terms.resize(in.terms.size());
            for (std::size_t t = 0; t < in.terms.size(); ++t) {
                out.terms[t].baby_step = in.terms[t].baby_step;
                out.terms[t].coefficients = in.terms[t].coefficients;
            }
        }
    }
}

template <class P, std::uint32_t NewQ, std::uint32_t NewDelta,
          std::uint32_t OldQ, std::uint32_t OldDelta, std::uint32_t PlainDelta>
inline void CKKSRebindCompactLinearInputScale(
    CKKSCompactRealLinearTransformPlan<P, NewQ, NewDelta, PlainDelta> &result,
    const CKKSCompactRealLinearTransformPlan<P, OldQ, OldDelta, PlainDelta> &source)
{
    result.has_direct = source.has_direct;
    result.has_conjugate = source.has_conjugate;
    CKKSRebindCompactLinearInputScale(result.direct, source.direct);
    CKKSRebindCompactLinearInputScale(result.conjugate, source.conjugate);
}

template <class P, std::uint32_t Q, std::uint32_t Delta, std::uint32_t PlainDelta>
inline void CKKSRealLinearTransform(
    CKKSPlainMulResult<P, Q, Delta, PlainDelta> &res,
    const CKKSCiphertext<P, Q, Delta> &ct,
    const CKKSCompactRealLinearTransformPlan<P, Q, Delta, PlainDelta> &plan,
    const auto &input_gk, const auto &output_gk)
{
    if (!plan.has_direct && !plan.has_conjugate)
        throw std::invalid_argument("empty compact CKKS real transform");
    auto term = std::make_unique<CKKSPlainMulResult<P, Q, Delta, PlainDelta>>();
    if (plan.has_direct) {
        CKKSLinearTransformBSGS<P, Q, Delta, PlainDelta>(*term, ct, plan.direct, input_gk, output_gk);
        res = *term;
    }
    if (plan.has_conjugate) {
        auto conjugated = std::make_unique<CKKSCiphertext<P, Q, Delta>>();
        CKKSConjugateSlots<P, Q>(conjugated->ct, ct.ct, input_gk);
        CKKSLinearTransformBSGS<P, Q, Delta, PlainDelta>(*term, *conjugated, plan.conjugate, input_gk, output_gk);
        if (!plan.has_direct) res = *term;
        else CKKSAddTRLWEInPlace<P, Q - PlainDelta>(res.ct, term->ct);
    }
    ckks_detail::reduceTRLWEToLevel<P, Q - PlainDelta>(res.ct);
}

// Packed S2C has exactly one diagonal at each offset. Generate and encode
// those diagonals directly; neither the complex N-by-N matrix nor a legacy
// full-width plan is ever allocated. Formula and group ordering match the
// ordinary builder coefficient for coefficient.
template <class P, std::uint32_t Q, std::uint32_t Delta, std::uint32_t PlainDelta>
inline void CKKSBuildPackedSlotToCoeffPlan(
    CKKSCompactRealLinearTransformPlan<P, Q, Delta, PlainDelta> &plan, int k_step)
{
    constexpr int half = P::n / 2;
    if (k_step <= 0 || k_step > half)
        throw std::invalid_argument("invalid compact CKKS BSGS step");
    const auto &slot_to_eval = ckks_detail::ckksSlotToEvalIndex<P>();
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    auto diagonal = std::make_unique<CKKSSlotVector<P>>();
    auto rotated = std::make_unique<CKKSSlotVector<P>>();
    auto plain = std::make_unique<Polynomial<P>>();
    auto build = [&](auto &out, bool conjugate) {
        out.k_step = k_step;
        out.groups.clear();
        out.groups.reserve((half + k_step - 1) / k_step);
        for (int d = 0; d < half; d++) {
            const int j2 = d / k_step;
            if (d % k_step == 0) {
                out.groups.emplace_back();
                out.groups.back().giant_step = j2;
                out.groups.back().terms.reserve(std::min(k_step, half - d));
            }
            for (int i = 0; i < half; i++) {
                const long double h = static_cast<long double>(2 * slot_to_eval[i] + 1);
                const long double theta = pi * h / static_cast<long double>(P::n);
                const int c = (i + d) % half;
                const auto e0 = std::complex<long double>(std::cos(theta * c), std::sin(theta * c));
                const auto e1 = std::complex<long double>(std::cos(theta * (c + half)), std::sin(theta * (c + half)));
                const auto value = static_cast<std::complex<double>>(0.5L *
                    (conjugate ? e0 + std::complex<long double>(0, 1) * e1 :
                                 e0 - std::complex<long double>(0, 1) * e1));
                (*diagonal)[i] = std::complex<double>{0, 0} + value;
            }
            rotateCKKSSlotVector<P>(*rotated, *diagonal, -k_step * j2);
            ckksSlotEncode<P, Q, PlainDelta>(*plain, *rotated);
            out.groups.back().terms.emplace_back();
            auto &term = out.groups.back().terms.back();
            term.baby_step = d % k_step;
            for (std::size_t n = 0; n < P::n; n++) term.coefficients[n] = out.pack((*plain)[n]);
        }
    };
    build(plan.direct, false);
    build(plan.conjugate, true);
    plan.has_direct = plan.has_conjugate = true;
}

}  // namespace TFHEpp
#endif
