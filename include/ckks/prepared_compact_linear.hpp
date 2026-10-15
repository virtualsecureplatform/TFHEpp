#pragma once

#include "compact_linear.hpp"
#include "prepared_auto.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Prepared key-independent plaintext spectra. Ciphertext rotations and
// products are still computed for every evaluation.
template <class P, std::uint32_t Q, std::uint32_t Delta,
          std::uint32_t PlainDelta,
          std::uint32_t Bits = torus_unsigned_digit_bits_v<P>>
struct CKKSPreparedCompactRealLinearTransformFFT {
    static_assert(use_multilimb_digit_fft_v<P>);
    static_assert(Bits > 0 && Bits <= 32 && 64 % Bits == 0);
    static constexpr std::size_t plain_digits = 64 / Bits;
    using PlainFFT = std::array<PolynomialInFD<P>,plain_digits>;
    struct Term { int baby_step; };
    struct Group {
        int giant_step;
        std::vector<Term> terms;
        std::size_t prepared_plain_digits;
        std::vector<std::vector<PolynomialInFD<P>>> prepared_plain_fft;
    };
    struct Plan { int k_step; std::vector<Group> groups; };
    bool has_direct, has_conjugate;
    Plan direct, conjugate;
    std::uint64_t spectrum_bytes = 0;

    CKKSPreparedCompactRealLinearTransformFFT(
        const CKKSCompactRealLinearTransformPlan<P,Q,Delta,PlainDelta> &source,
        std::uint64_t budget_bytes)
        : has_direct(source.has_direct), has_conjugate(source.has_conjugate)
    {
        if (!has_direct && !has_conjugate)
            throw std::invalid_argument("empty prepared compact CKKS plan");
        const auto validate = [&](const auto &plan) {
            if (plan.k_step <= 0 || plan.k_step > P::n/2 || plan.groups.empty())
                throw std::invalid_argument("invalid prepared compact CKKS schedule");
            for (const auto &group : plan.groups) {
                if (group.giant_step < 0 || std::int64_t(group.giant_step)*plan.k_step >= P::n/2 || group.terms.empty())
                    throw std::invalid_argument("invalid prepared compact CKKS giant step");
                std::uint64_t maximum=0;
                for (const auto &term:group.terms)
                    for (const auto value:term.coefficients)
                        maximum=std::max(maximum,value<0?std::uint64_t{0}-std::uint64_t(value):std::uint64_t(value));
                const auto digits=std::max<std::size_t>(1,(std::bit_width(maximum)+Bits-1)/Bits);
                const auto term_bytes=digits*sizeof(PolynomialInFD<P>);
                for (const auto &term : group.terms) {
                    if (term.baby_step < 0 || term.baby_step >= plan.k_step)
                        throw std::invalid_argument("invalid prepared compact CKKS baby step");
                    if (spectrum_bytes > budget_bytes || term_bytes > budget_bytes-spectrum_bytes)
                        throw std::invalid_argument("prepared compact CKKS spectrum budget exceeded");
                    spectrum_bytes += term_bytes;
                }
            }
        };
        if (has_direct) validate(source.direct);
        if (has_conjugate) validate(source.conjugate);
        const auto prepare = [&](auto &out,const auto &plan) {
            out.k_step=plan.k_step; out.groups.resize(plan.groups.size());
            for (std::size_t g=0;g<plan.groups.size();++g) {
                const auto &raw=plan.groups[g]; auto &group=out.groups[g];
                group.giant_step=raw.giant_step;
                group.terms.resize(raw.terms.size());
                std::uint64_t maximum=0;
                for (const auto &term:raw.terms)
                    for (const auto value:term.coefficients)
                        maximum=std::max(maximum,value<0?std::uint64_t{0}-std::uint64_t(value):std::uint64_t(value));
                group.prepared_plain_digits=std::max<std::size_t>(1,(std::bit_width(maximum)+Bits-1)/Bits);
                group.prepared_plain_fft.resize(raw.terms.size());
                for (auto &term:group.prepared_plain_fft) term.resize(group.prepared_plain_digits);
                for (std::size_t t=0;t<raw.terms.size();++t)
                    group.terms[t].baby_step=raw.terms[t].baby_step;
#pragma omp parallel
                {
                    auto digit=std::make_unique<Polynomial<P>>();
#pragma omp for collapse(2) schedule(static)
                    for (std::size_t t=0;t<raw.terms.size();++t)
                        for (std::size_t d=0;d<group.prepared_plain_digits;++d) {
                            for (std::size_t i=0;i<P::n;++i) {
                                const auto value=ckks_detail::reduceToLevel<P,Q>(
                                    P::T::from_signed_i64(raw.terms[t].coefficients[i]));
                                const auto [negative,magnitude]=ckks_detail::smallSignedMagnitude<P,Q>(value);
                                const typename P::T part{(magnitude>>(d*Bits)) & ((std::uint64_t{1}<<Bits)-1)};
                                (*digit)[i]=negative ? typename P::T{}-part : part;
                            }
                            TwistIFFTDigit<P>(group.prepared_plain_fft[t][d],*digit);
                        }
                }
            }
        };
        if (has_direct) prepare(direct,source.direct);
        if (has_conjugate) prepare(conjugate,source.conjugate);
    }
};

template <class P,std::uint32_t Q,std::uint32_t Delta,std::uint32_t PlainDelta,
          std::uint32_t Bits,class InputKey,class OutputKey>
inline void CKKSRealLinearTransform(
    CKKSPlainMulResult<P,Q,Delta,PlainDelta> &result,
    const CKKSCiphertext<P,Q,Delta> &input,
    const CKKSPreparedCompactRealLinearTransformFFT<P,Q,Delta,PlainDelta,Bits> &plan,
    const InputKey &input_key,const OutputKey &output_key)
{
    const auto evaluate=[&](auto &out,const auto &ct,const auto &part) {
        std::vector<bool> used(part.k_step,false);
        for (const auto &group:part.groups)
            for (const auto &term:group.terms) used[term.baby_step]=true;
        const auto provider=[&](std::size_t i)->const auto& { return part.groups[i]; };
        CKKSLinearTransformBSGSWithGroupProvider<P,Q,Delta,PlainDelta,
            InputKey,OutputKey,decltype(provider)&,Bits>(
            out,ct,part.k_step,used,part.groups.size(),provider,input_key,output_key);
    };
    auto term=std::make_unique<CKKSPlainMulResult<P,Q,Delta,PlainDelta>>();
    if (plan.has_direct) { evaluate(*term,input,plan.direct); result=*term; }
    if (plan.has_conjugate) {
        auto conjugated=std::make_unique<CKKSCiphertext<P,Q,Delta>>();
        CKKSConjugateSlots<P,Q>(conjugated->ct,input.ct,input_key);
        evaluate(*term,*conjugated,plan.conjugate);
        if (!plan.has_direct) result=*term;
        else CKKSAddTRLWEInPlace<P,Q-PlainDelta>(result.ct,term->ct);
    }
    ckks_detail::reduceTRLWEToLevel<P,Q-PlainDelta>(result.ct);
}

// Explicit retention policy: default evaluation-local caches still clear at
// BSGS boundaries. This adapter is prepared before use and retains all public
// rotation/conjugation spectra across repeated serial evaluations.
template <class P,std::uint32_t Q>
struct CKKSRetainedPreparedGaloisKeyFFT : CKKSLazyPreparedGaloisKeyFFT<P,Q> {
    using Base=CKKSLazyPreparedGaloisKeyFFT<P,Q>;
    explicit CKKSRetainedPreparedGaloisKeyFFT(const CKKSGaloisKey<P,Q> &raw)
        : Base(raw)
    {
        for (std::size_t bit=0;bit<P::nbit-1;++bit) this->get(bit);
        this->get(P::nbit);
    }
    void clear_cached_keys() const {}
};
} // namespace TFHEpp
#endif
