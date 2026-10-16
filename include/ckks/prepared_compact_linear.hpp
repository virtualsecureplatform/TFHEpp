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
        std::vector<std::uint64_t> prepared_plain_l1;
        bool prepared_adaptive_batches=true;
        std::vector<std::vector<PolynomialInFD<P>>> prepared_plain_fft;
    };
    struct Plan { int k_step; std::vector<Group> groups; };
    bool has_direct, has_conjugate;
    Plan direct, conjugate;
    std::uint64_t spectrum_bytes = 0;

    CKKSPreparedCompactRealLinearTransformFFT(
        const CKKSCompactRealLinearTransformPlan<P,Q,Delta,PlainDelta> &source,
        std::uint64_t budget_bytes,bool adaptive_batches=true)
        : has_direct(source.has_direct), has_conjugate(source.has_conjugate)
    {
        if (!has_direct && !has_conjugate)
            throw std::invalid_argument("empty prepared compact CKKS plan");
        struct Work { std::size_t group; std::vector<std::size_t> terms; std::size_t digits,bits; };
        const auto shape = [&](const auto &plan) {
            std::vector<Work> work;
            if (plan.k_step <= 0 || plan.k_step > P::n/2 || plan.groups.empty())
                throw std::invalid_argument("invalid prepared compact CKKS schedule");
            for (std::size_t g=0;g<plan.groups.size();++g) {
                const auto &group=plan.groups[g];
                if (group.giant_step < 0 || std::int64_t(group.giant_step)*plan.k_step >= P::n/2 || group.terms.empty())
                    throw std::invalid_argument("invalid prepared compact CKKS giant step");
                Work item{g,{},0,0}; std::uint64_t maximum=0;
                for (std::size_t t=0;t<group.terms.size();++t) {
                    const auto &term=group.terms[t];
                    if (term.baby_step < 0 || term.baby_step >= plan.k_step)
                        throw std::invalid_argument("invalid prepared compact CKKS baby step");
                    std::uint64_t term_maximum=0;
                    for (const auto value:term.coefficients)
                        term_maximum=std::max(term_maximum,value<0?std::uint64_t{0}-std::uint64_t(value):std::uint64_t(value));
                    if (term_maximum) item.terms.push_back(t);
                    maximum=std::max(maximum,term_maximum);
                }
                if (!maximum) continue;
                item.bits=std::bit_width(maximum);
                item.digits=(item.bits+Bits-1)/Bits;
                const auto term_bytes=item.digits*sizeof(PolynomialInFD<P>);
                for (std::size_t i=0;i<item.terms.size();++i) {
                    if (spectrum_bytes>budget_bytes || term_bytes>budget_bytes-spectrum_bytes)
                        throw std::invalid_argument("prepared compact CKKS spectrum budget exceeded");
                    spectrum_bytes+=term_bytes;
                }
                work.push_back(std::move(item));
            }
            return work;
        };
        std::vector<Work> direct_work,conjugate_work;
        if (has_direct) direct_work=shape(source.direct);
        if (has_conjugate) conjugate_work=shape(source.conjugate);
        has_direct=!direct_work.empty();has_conjugate=!conjugate_work.empty();
        const auto prepare = [&](auto &out,const auto &plan,const auto &work) {
            out.k_step=plan.k_step; out.groups.resize(work.size());
            for (std::size_t g=0;g<work.size();++g) {
                const auto &item=work[g];const auto &raw=plan.groups[item.group];auto &group=out.groups[g];
                group.giant_step=raw.giant_step;
                group.terms.resize(item.terms.size());
                group.prepared_plain_digits=item.digits;
                group.prepared_adaptive_batches=adaptive_batches;
                group.prepared_plain_fft.resize(item.terms.size());
                group.prepared_plain_l1.resize(item.digits);
                std::vector<std::vector<std::uint64_t>> norms(item.terms.size(),std::vector<std::uint64_t>(item.digits));
                for (auto &term:group.prepared_plain_fft) term.resize(group.prepared_plain_digits);
                for (std::size_t t=0;t<item.terms.size();++t)
                    group.terms[t].baby_step=raw.terms[item.terms[t]].baby_step;
#pragma omp parallel
                {
                    auto digit=std::make_unique<Polynomial<P>>();
#pragma omp for collapse(2) schedule(static)
                    for (std::size_t t=0;t<item.terms.size();++t)
                        for (std::size_t d=0;d<group.prepared_plain_digits;++d) {
                            std::uint64_t norm=0;
                            for (std::size_t i=0;i<P::n;++i) {
                                const auto value=ckks_detail::reduceToLevel<P,Q>(
                                    P::T::from_signed_i64(raw.terms[item.terms[t]].coefficients[i]));
                                const auto [negative,magnitude]=ckks_detail::smallSignedMagnitude<P,Q>(value);
                                const typename P::T part{(magnitude>>(d*Bits)) & ((std::uint64_t{1}<<Bits)-1)};
                                norm+=std::uint64_t(part);
                                (*digit)[i]=negative ? typename P::T{}-part : part;
                            }
                            norms[t][d]=norm;
                            TwistIFFTDigit<P>(group.prepared_plain_fft[t][d],*digit);
                        }
                }
                for(const auto &term:norms)for(std::size_t d=0;d<item.digits;++d)
                    group.prepared_plain_l1[d]=std::max(group.prepared_plain_l1[d],term[d]);
            }
        };
        if (has_direct) prepare(direct,source.direct,direct_work);
        if (has_conjugate) prepare(conjugate,source.conjugate,conjugate_work);
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
    if (!plan.has_direct && !plan.has_conjugate) {
        for (auto &component:result.ct) component.fill({});
        return;
    }
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
// rotation/conjugation spectra across evaluations. The explicit coordinator
// hook makes subsequent concurrent rotation lookups read-only.
template <class P,std::uint32_t Q,std::uint32_t PrimaryBits=P::B̅gbit>
struct CKKSRetainedPreparedGaloisKeyFFT : CKKSLazyPreparedGaloisKeyFFT<P,Q,PrimaryBits> {
    using Base=CKKSLazyPreparedGaloisKeyFFT<P,Q,PrimaryBits>;
    explicit CKKSRetainedPreparedGaloisKeyFFT(const CKKSGaloisKey<P,Q> &raw)
        : Base(raw)
    {
        for (std::size_t bit=0;bit<P::nbit-1;++bit) this->get(bit);
        this->get(P::nbit);
    }
    void clear_cached_keys() const {}
    void prepare_parallel_rotations(std::uint32_t mask) const
    {
        constexpr auto allowed=(std::uint32_t{1}<<(P::nbit-1))-1;
        if(mask&~allowed)throw std::invalid_argument("invalid retained rotation bit mask");
        for(std::size_t bit=0;bit<P::nbit-1;++bit)if((mask>>bit)&1)this->get(bit);
    }
};
} // namespace TFHEpp
#endif
