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
namespace ckks_detail {
// res = a*b or res += a*b on complex indices [begin,end), in MulInFD's layout.
template <std::uint32_t N>
inline void mulInFDRange(double *__restrict__ res,const double *__restrict__ a,
    const double *__restrict__ b,std::size_t begin,std::size_t end,bool accumulate)
{
#if defined(USE_INTERLEAVED_FORMAT) && defined(__AVX2__)
    constexpr std::size_t imaginary=1,stride=2;
#else
    constexpr std::size_t imaginary=N/2,stride=1;
#endif
    const double *are=a+stride*begin,*aim=are+imaginary,*bre=b+stride*begin,*bim=bre+imaginary;
    double *rre=res+stride*begin,*rim=rre+imaginary;
    const std::size_t count=end-begin;
    if(accumulate){
#pragma omp simd
        for(std::size_t i=0;i<count;++i){
            const double x=are[stride*i],y=aim[stride*i],u=bre[stride*i],v=bim[stride*i];
            rre[stride*i]+=x*u-y*v;
            rim[stride*i]+=x*v+y*u;
        }
    }else{
#pragma omp simd
        for(std::size_t i=0;i<count;++i){
            const double x=are[stride*i],y=aim[stride*i],u=bre[stride*i],v=bim[stride*i];
            rre[stride*i]=x*u-y*v;
            rim[stride*i]=x*v+y*u;
        }
    }
}
} // namespace ckks_detail

// Same map as CKKSRealLinearTransform on a prepared plan, arranged for many
// cores. The serial-per-group kernel streams every plaintext spectrum once
// per ciphertext digit and every baby spectrum once per group, which is
// bounded by memory bandwidth. Here
//   * baby rotations of equal depth run concurrently,
//   * products are blocked on spectrum index, so the baby spectra of a block
//     stay in cache across all groups and each plaintext block is read once,
//   * giant rotations are merged by a binary tree: giant-1 single-key
//     rotations, each level concurrent, in place of a chain per group.
// Both keys must be retained (read-only lookups) and k_step a power of two.
template <class P,std::uint32_t Q,std::uint32_t Delta,std::uint32_t PlainDelta,
          std::uint32_t Bits,class InputKey,class OutputKey>
inline void CKKSRealLinearTransformBlocked(
    CKKSPlainMulResult<P,Q,Delta,PlainDelta> &result,
    const CKKSCiphertext<P,Q,Delta> &input,
    const CKKSPreparedCompactRealLinearTransformFFT<P,Q,Delta,PlainDelta,Bits> &plan,
    const InputKey &input_key,const OutputKey &output_key)
{
    using T=typename P::T;
    static_assert(requires { input_key.prepare_parallel_rotations(std::uint32_t{0}); });
    static_assert(requires { output_key.prepare_parallel_rotations(std::uint32_t{0}); });
    constexpr std::uint32_t out_q=Q-PlainDelta;
    constexpr int width=std::numeric_limits<T>::digits;
    constexpr int torus_digits=(Q+Bits-1)/Bits;
    constexpr int torus_first=width/Bits-torus_digits;
    constexpr T torus_mask=(T{1}<<Bits)-T{1};
    constexpr T offset=[]{
        T value=0;
        for(int j=0;j<std::numeric_limits<T>::digits/int(Bits);j++)
            value+=(T{1}<<(Bits-1))<<(std::numeric_limits<T>::digits-(j+1)*int(Bits));
        return value;
    }();
    constexpr int fd_margin=std::numeric_limits<double>::digits-(2*int(Bits)+int(P::nbit)+3);
    static_assert(fd_margin>0);
    constexpr std::size_t max_fused_terms=std::size_t{1}<<std::min(fd_margin-1,30);
    constexpr std::size_t components=P::k+1,spectrum=P::n/2;
    // 32 complex values: one block of all baby spectra fits a private cache.
    constexpr std::size_t block=std::min<std::size_t>(32,spectrum);
    if(!plan.has_direct&&!plan.has_conjugate){
        for(auto &component:result.ct)component.fill({});
        return;
    }
#ifdef _OPENMP
    const std::size_t threads=static_cast<std::size_t>(omp_get_max_threads());
#else
    const std::size_t threads=1;
#endif
    const auto evaluate=[&](TRLWE<P> &out,const TRLWE<P> &ct,const auto &part){
        const int k_step=part.k_step;
        const bool trace=ckks_detail::ckks_linear_trace_enabled();
        auto mark=std::chrono::steady_clock::now();
        const auto stage=[&](const char *name){
            if(!trace)return;
            ckks_detail::ckks_trace_linear_event(name,ckks_detail::ckks_trace_elapsed_ms(mark));
            mark=std::chrono::steady_clock::now();
        };
        if(k_step<=0||(k_step&(k_step-1)))
            throw std::invalid_argument("blocked CKKS linear transform needs a power-of-two baby step");
        std::vector<bool> used(k_step,false);
        int giant_limit=1;
        for(const auto &group:part.groups){
            for(const auto &term:group.terms)used[term.baby_step]=true;
            giant_limit=std::max(giant_limit,group.giant_step+1);
        }
        auto baby=std::make_unique<ckks_detail::CKKSSparseBabyStepTable<P>>(k_step);
        ckks_detail::CKKSBuildSparseBabyStepRotationTable<P,Q>(*baby,ct,used,input_key,threads);
        stage("blocked_baby_rotation_table");

        using BabyFFT=std::array<std::array<PolynomialInFD<P>,torus_digits>,components>;
        std::vector<std::unique_ptr<BabyFFT>> baby_fft(k_step);
        for(int j1=0;j1<k_step;++j1)if(used[j1])baby_fft[j1]=std::make_unique<BabyFFT>();
#pragma omp parallel
        {
            auto digit=std::make_unique<std::array<std::int64_t,P::n>>();
            auto adjusted=std::make_unique<Polynomial<P>>();
#pragma omp for collapse(2) schedule(dynamic)
            for(int j1=0;j1<k_step;++j1)
                for(std::size_t c=0;c<components;++c){
                    if(!used[j1])continue;
                    const auto &source=(*(*baby)[j1])[c];
                    for(std::size_t n=0;n<P::n;++n)
                        (*adjusted)[n]=ckks_detail::centeredLevelToTorus<P,Q>(source[n])+offset;
                    for(int j=0;j<torus_digits;++j){
                        const int shift=width-(torus_first+j+1)*int(Bits);
                        for(std::size_t n=0;n<P::n;++n)
                            (*digit)[n]=static_cast<std::int64_t>(static_cast<std::uint64_t>(
                                ((*adjusted)[n]>>shift)&torus_mask))-(std::int64_t{1}<<(Bits-1));
                        TwistIFFTDigit<P>((*baby_fft[j1])[c][j],*digit);
                    }
                }
        }
        baby.reset();
        stage("blocked_baby_digit_fft");

        // One accumulator per exact batch of fused digit products.
        struct Slot { std::uint32_t group,component,row,digit,begin,end; int shift; };
        std::vector<Slot> slots;
        std::vector<std::size_t> first(part.groups.size()+1,0);
        for(std::size_t g=0;g<part.groups.size();++g){
            const auto &group=part.groups[g];
            const std::size_t digits=group.prepared_plain_digits,terms=group.terms.size();
            if(group.prepared_plain_fft.size()!=terms||group.prepared_plain_l1.size()!=digits)
                throw std::invalid_argument("prepared CKKS plaintext spectrum count mismatch");
            for(std::size_t c=0;c<components;++c)
                for(int j=0;j<torus_digits;++j)
                    for(std::size_t d=0;d<digits;++d){
                        const int shift=width-(torus_first+j+1)*int(Bits)+int(d*Bits);
                        if(shift>=int(Q))continue;
                        std::size_t fused=max_fused_terms;
                        if(group.prepared_adaptive_batches){
                            const int input_bits=j==0?int(Q)-(torus_digits-1)*int(Bits):int(Bits);
                            const auto norm=group.prepared_plain_l1[d];
                            if(!norm)continue;
                            constexpr int budget_bits=std::max(47,int(2*Bits+P::nbit)-1);
                            fused=std::min<std::uint64_t>(std::uint64_t{1}<<30,
                                (std::uint64_t{1}<<budget_bits)/((std::uint64_t{1}<<(input_bits-1))*norm));
                            if(!fused)throw std::invalid_argument("prepared CKKS plaintext norm exceeds digit range");
                        }
                        for(std::size_t begin=0;begin<terms;begin+=fused)
                            slots.push_back({std::uint32_t(g),std::uint32_t(c),std::uint32_t(j),
                                std::uint32_t(d),std::uint32_t(begin),
                                std::uint32_t(std::min(terms,begin+fused)),shift});
                    }
            first[g+1]=slots.size();
        }
        const auto sums=std::make_unique_for_overwrite<PolynomialInFD<P>[]>(slots.size());
#pragma omp parallel for schedule(dynamic)
        for(std::size_t begin=0;begin<spectrum;begin+=block){
            const std::size_t end=begin+block;
            for(std::size_t s=0;s<slots.size();++s){
                const auto &slot=slots[s];
                const auto &group=part.groups[slot.group];
                for(std::size_t t=slot.begin;t<slot.end;++t)
                    ckks_detail::mulInFDRange<P::n>(sums[s].data(),
                        (*baby_fft[group.terms[t].baby_step])[slot.component][slot.row].data(),
                        group.prepared_plain_fft[t][slot.digit].data(),begin,end,t!=slot.begin);
            }
        }
        baby_fft.clear();
        stage("blocked_products");

        // Recombine each group, then merge the giant rotations bit by bit:
        // sum_g rot_{g k}(x_g) = sum_lo rot_{lo k}(x_lo + rot_{2^b k}(x_{lo+2^b})).
        const std::size_t giant_bits=std::bit_width(static_cast<unsigned>(giant_limit-1));
        std::vector<std::unique_ptr<TRLWE<P>>> giant(std::size_t{1}<<giant_bits);
        for(const auto &group:part.groups){
            if(giant[group.giant_step])throw std::invalid_argument("repeated CKKS giant step");
            giant[group.giant_step]=std::make_unique<TRLWE<P>>();
        }
#pragma omp parallel
        {
            auto product=std::make_unique<Polynomial<P>>();
#pragma omp for schedule(dynamic)
            for(std::size_t g=0;g<part.groups.size();++g){
                auto &inner=*giant[part.groups[g].giant_step];
                for(auto &component:inner)component.fill({});
                for(std::size_t s=first[g];s<first[g+1];++s){
                    TwistFFTDigitProduct<P>(*product,sums[s]);
                    auto &component=inner[slots[s].component];
                    for(std::size_t n=0;n<P::n;++n)component[n]+=(*product)[n]<<slots[s].shift;
                }
                for(auto &component:inner)
                    for(auto &value:component)
                        value=ckks_detail::roundedLevelRightShift<P,Q,PlainDelta>(value);
                ckks_detail::reduceTRLWEToLevel<P,out_q>(inner);
            }
        }
        stage("blocked_recombine");
        for(std::size_t bit=giant_bits;bit-->0;){
            const std::size_t half=std::size_t{1}<<bit;
            const int steps=k_step<<bit;
#pragma omp parallel
            {
                auto rotated=std::make_unique<TRLWE<P>>();
#pragma omp for schedule(dynamic)
                for(std::size_t low=0;low<half;++low){
                    if(!giant[low+half])continue;
                    CKKSRotateSlots<P,out_q>(*rotated,*giant[low+half],steps,output_key);
                    if(giant[low])CKKSAddTRLWEInPlace<P,out_q>(*giant[low],*rotated);
                    else giant[low]=std::make_unique<TRLWE<P>>(*rotated);
                    giant[low+half].reset();
                }
            }
        }
        stage("blocked_giant_tree");
        out=*giant[0];
        ckks_detail::reduceTRLWEToLevel<P,out_q>(out);
    };
    if(plan.has_direct)evaluate(result.ct,input.ct,plan.direct);
    if(plan.has_conjugate){
        auto conjugated=std::make_unique<CKKSCiphertext<P,Q,Delta>>();
        CKKSConjugateSlots<P,Q>(conjugated->ct,input.ct,input_key);
        if(!plan.has_direct)evaluate(result.ct,conjugated->ct,plan.conjugate);
        else{
            auto term=std::make_unique<TRLWE<P>>();
            evaluate(*term,conjugated->ct,plan.conjugate);
            CKKSAddTRLWEInPlace<P,Q-PlainDelta>(result.ct,*term);
        }
    }
    ckks_detail::reduceTRLWEToLevel<P,Q-PlainDelta>(result.ct);
}
} // namespace TFHEpp
#endif
