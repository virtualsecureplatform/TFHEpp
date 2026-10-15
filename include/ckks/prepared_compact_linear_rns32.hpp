#pragma once
#include "prepared_compact_linear.hpp"
#include <modular_rns32.hpp>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
struct CKKSLinearRNSTimings {
    double baby_rotations=0,baby_ntt=0,matrix=0,inverse_crt=0,giant_rotations=0;
    std::size_t giant_rotation_calls=0,giant_rotation_bits=0;
};

// The public diagonals use a complete residue representation. Each giant
// group is reconstructed modulo 2^Q and rounded once, exactly as in BSGS.
template<class P,std::uint32_t Q,std::uint32_t Delta,std::uint32_t PlainDelta,
         std::size_t Count=6>
class CKKSPreparedCompactRealLinearTransformRNS32 {
    static_assert(is_multilimb_uint_v<typename P::T> && Q>PlainDelta && Q<512);
    static_assert(Count>0 && Count<=modular_rns32::prime_count);
    using Input=CKKSCiphertext<P,Q,Delta>;
    using Result=CKKSPlainMulResult<P,Q,Delta,PlainDelta>;
    using Small=MultiLimbUInt<(Q+63)/64>;
    static constexpr auto output_q=Q-PlainDelta;
    static const auto &plans(){
        static const auto value=[] {
            std::vector<modular_ntt32::NegacyclicNTTPlan> result;result.reserve(Count);
            for(std::size_t p=0;p<Count;++p)result.emplace_back(P::n,modular_rns32::primes[p]);
            return result;
        }();return value;
    }
public:
    static constexpr std::uint64_t term_spectrum_bytes=std::uint64_t{Count}*P::n*sizeof(std::uint32_t);
    struct Term{int baby_step;};
    struct Group{int giant_step;std::vector<Term> terms;std::array<std::vector<std::uint32_t>,Count> spectra;};
    struct Plan{int k_step=0;std::vector<Group> groups;};
    bool has_direct,has_conjugate;
    Plan direct,conjugate;
    std::uint64_t spectrum_bytes=0;
    CKKSPreparedCompactRealLinearTransformRNS32(
        const CKKSCompactRealLinearTransformPlan<P,Q,Delta,PlainDelta>&source,std::uint64_t budget)
        :has_direct(source.has_direct),has_conjugate(source.has_conjugate)
    {
        if(!has_direct&&!has_conjugate)throw std::invalid_argument("empty prepared RNS linear plan");
        struct Work{std::size_t group;std::vector<std::size_t> terms;};
        const auto shape=[&](const auto &plan){
            std::vector<Work> work;
            if(plan.k_step<=0 || plan.k_step>P::n/2 || plan.groups.empty())throw std::invalid_argument("invalid RNS linear schedule");
            for(std::size_t g=0;g<plan.groups.size();++g){
                const auto &raw=plan.groups[g];
                if(raw.giant_step<0 || std::int64_t(raw.giant_step)*plan.k_step>=P::n/2 || raw.terms.empty())throw std::invalid_argument("invalid RNS giant step");
                Work item{g,{}};unsigned __int128 norm=0;
                for(std::size_t t=0;t<raw.terms.size();++t){
                    const auto &term=raw.terms[t];
                    if(term.baby_step<0 || term.baby_step>=plan.k_step)throw std::invalid_argument("invalid RNS baby step");
                    unsigned __int128 term_norm=0;
                    for(const auto v:term.coefficients)term_norm+=v<0?std::uint64_t{0}-std::uint64_t(v):std::uint64_t(v);
                    if(term_norm){item.terms.push_back(t);norm+=term_norm;}
                }
                if(item.terms.empty())continue;
                const auto norm_bits=std::uint64_t(norm>>64)?64+std::bit_width(std::uint64_t(norm>>64)):std::bit_width(std::uint64_t(norm));
                if(Q+norm_bits+1>512)throw std::invalid_argument("RNS linear convolution bound too large");
                MultiLimbUInt<8> bound;bound.limb[0]=std::uint64_t(norm);bound.limb[1]=std::uint64_t(norm>>64);bound<<=Q+1;
                if(bound>=modular_rns32::basisProduct<Count>())throw std::invalid_argument("RNS linear convolution exceeds centered CRT range");
                for(std::size_t t=0;t<item.terms.size();++t){
                    if(spectrum_bytes>budget || term_spectrum_bytes>budget-spectrum_bytes)throw std::invalid_argument("RNS linear spectrum budget exceeded");
                    spectrum_bytes+=term_spectrum_bytes;
                }
                work.push_back(std::move(item));
            }
            return work;
        };
        std::vector<Work> dw,cw;if(has_direct)dw=shape(source.direct);if(has_conjugate)cw=shape(source.conjugate);
        has_direct=!dw.empty();has_conjugate=!cw.empty();
        const auto prepare=[&](auto&out,const auto&plan,const auto&work){
            out.k_step=plan.k_step;out.groups.resize(work.size());
            for(std::size_t g=0;g<work.size();++g){
                const auto &item=work[g];const auto &raw=plan.groups[item.group];auto &group=out.groups[g];
                group.giant_step=raw.giant_step;group.terms.resize(item.terms.size());
                for(std::size_t t=0;t<item.terms.size();++t)group.terms[t].baby_step=raw.terms[item.terms[t]].baby_step;
                for(auto &spectrum:group.spectra)spectrum.resize(item.terms.size()*P::n);
#pragma omp parallel for collapse(2) schedule(static)
                for(std::size_t p=0;p<Count;++p)
                    for(std::size_t t=0;t<item.terms.size();++t){
                        const auto &arithmetic=plans()[p].arithmetic();const auto modulus=arithmetic.modulus();
                        auto values=std::span(group.spectra[p]).subspan(t*P::n,P::n);
                        for(std::size_t i=0;i<P::n;++i){
                            const auto v=raw.terms[item.terms[t]].coefficients[i];const auto magnitude=v<0?std::uint64_t{0}-std::uint64_t(v):std::uint64_t(v);
                            const auto reduced=magnitude%modulus;values[i]=arithmetic.encode(v<0&&reduced?modulus-reduced:reduced);
                        }
                        plans()[p].forward(values);
                    }
            }
        };
        // Initialize the transform plans before entering parallel preparation.
        if(has_direct||has_conjugate)(void)plans();
        if(has_direct)prepare(direct,source.direct,dw);if(has_conjugate)prepare(conjugate,source.conjugate,cw);
    }
private:
    template<class InputKey,class OutputKey>
    void applyPart(Result&result,const Input&input,const Plan&plan,const InputKey&ik,const OutputKey&ok,CKKSLinearRNSTimings*timings,bool fused_terms,bool horner_giants)const {
        using Clock=std::chrono::steady_clock;
        auto start=Clock::now();std::vector<bool> used(plan.k_step,false);
        for(const auto &g:plan.groups)for(const auto &t:g.terms)used[t.baby_step]=true;
        auto baby=std::make_unique<ckks_detail::CKKSSparseBabyStepTable<P>>(plan.k_step);
        ckks_detail::CKKSBuildSparseBabyStepRotationTable<P,Q>(*baby,input.ct,used,ik);
        ckks_detail::maybe_clear_cached_keys(ik);
        if(timings)timings->baby_rotations+=std::chrono::duration<double>(Clock::now()-start).count();start=Clock::now();
        using Rows=std::array<std::vector<std::uint32_t>,Count>;
        std::vector<std::array<Rows,P::k+1>> transformed(plan.k_step);
        static const modular_rns32::CanonicalEncoder<P::T::limbs,Q,Count> encoder;
#pragma omp parallel for schedule(static)
        for(std::size_t b=0;b<used.size();++b)if(used[b])
            for(std::size_t c=0;c<=P::k;++c)for(std::size_t p=0;p<Count;++p){
                auto &values=transformed[b][c][p];values.resize(P::n);
                for(std::size_t i=0;i<P::n;++i)values[i]=encoder.residue((*(*baby)[b])[c][i],p);
                plans()[p].forward(values);
            }
        baby.reset();if(timings)timings->baby_ntt+=std::chrono::duration<double>(Clock::now()-start).count();
        std::array<Rows,P::k+1> sums;
        auto inner=std::make_unique<TRLWE<P>>(),rotated=std::make_unique<TRLWE<P>>();
        static const modular_rns32::CenteredCRT<Small::limbs,Count> crt;
        std::vector<std::size_t> order;
        if(horner_giants){
            order.resize(plan.groups.size());for(std::size_t i=0;i<order.size();++i)order[i]=i;
            std::stable_sort(order.begin(),order.end(),[&](auto a,auto b){return plan.groups[a].giant_step>plan.groups[b].giant_step;});
        }
        const auto rotate=[&](auto&out,const auto&in,int steps){
            if(timings){++timings->giant_rotation_calls;timings->giant_rotation_bits+=std::popcount(std::uint32_t(steps));}
            CKKSRotateSlots<P,output_q>(out,in,steps,ok);
        };
        bool initialized=false;int previous_giant=0;
        for(std::size_t g=0;g<plan.groups.size();++g){
            const auto &group=plan.groups[horner_giants?order[g]:g];
            for(auto &component:sums)for(auto &row:component)row.assign(P::n,0);
            start=Clock::now();
#pragma omp parallel for collapse(3) schedule(static)
            for(std::size_t c=0;c<=P::k;++c)for(std::size_t p=0;p<Count;++p)for(std::size_t begin=0;begin<P::n;begin+=1024){
                const auto &arithmetic=plans()[p].arithmetic();auto *out=sums[c][p].data();
                const auto end=std::min(begin+1024,std::size_t(P::n));
                std::size_t first=0;
                if(fused_terms)for(;first+4<=group.terms.size();first+=4){
                    std::array<const std::uint32_t*,4> left,right;
                    for(std::size_t t=0;t<4;++t){
                        left[t]=transformed[group.terms[first+t].baby_step][c][p].data();
                        right[t]=group.spectra[p].data()+(first+t)*P::n;
                    }
                    std::size_t i=begin;
#if defined(__AVX2__)
                    for(;i+8<=end;i+=8){
                        const auto sum=arithmetic.accumulateDot4Vector(
                            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(out+i)),
                            [&](std::size_t t){return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(left[t]+i));},
                            [&](std::size_t t){return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(right[t]+i));});
                        _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+i),sum);
                    }
#endif
                    for(;i<end;++i){
                        std::uint64_t dot=0;for(std::size_t t=0;t<4;++t)dot+=std::uint64_t(left[t][i])*right[t][i];
                        const auto sum=out[i]+arithmetic.reduceLazy(dot);
                        out[i]=sum>=2*arithmetic.modulus()?sum-2*arithmetic.modulus():sum;
                    }
                }
                for(std::size_t t=first;t<group.terms.size();++t){
                    const auto *left=transformed[group.terms[t].baby_step][c][p].data();const auto *right=group.spectra[p].data()+t*P::n;std::size_t i=begin;
#if defined(__AVX2__)
                    for(;i+8<=end;i+=8){
                        const auto sum=arithmetic.accumulateVector(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(left+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(right+i)));
                        _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+i),sum);
                    }
#endif
                    for(;i<end;++i){const auto v=out[i]+arithmetic.multiplyLazy(left[i],right[i]);out[i]=v>=2*arithmetic.modulus()?v-2*arithmetic.modulus():v;}
                }
            }
            if(timings)timings->matrix+=std::chrono::duration<double>(Clock::now()-start).count();start=Clock::now();
#pragma omp parallel for collapse(2) schedule(static)
            for(std::size_t c=0;c<=P::k;++c)for(std::size_t p=0;p<Count;++p){
                auto &values=sums[c][p];const auto &a=plans()[p].arithmetic();for(auto &v:values)v=a.canonical(v);plans()[p].inverse(values);
            }
#pragma omp parallel
            {
                std::vector<Small> recovered(1024);
#pragma omp for collapse(2) schedule(static)
                for(std::size_t c=0;c<=P::k;++c)for(std::size_t begin=0;begin<P::n;begin+=1024){
                    const auto count=std::min(std::size_t{1024},std::size_t(P::n)-begin);
                    std::array<std::span<std::uint32_t>,Count> residues;for(std::size_t p=0;p<Count;++p)residues[p]=std::span(sums[c][p]).subspan(begin,count);
                    crt.reconstruct(residues,std::span(recovered).first(count));
                    for(std::size_t i=0;i<count;++i){typename P::T value{};for(std::size_t l=0;l<Small::limbs;++l)value.limb[l]=recovered[i].limb[l];(*inner)[c][begin+i]=ckks_detail::roundedLevelRightShift<P,Q,PlainDelta>(value);}
                }
            }
            if(timings)timings->inverse_crt+=std::chrono::duration<double>(Clock::now()-start).count();start=Clock::now();
            if(horner_giants){
                // Maintain sum_h Rot((h-g)*k,inner_h) at the current g.
                // Adjacent dense groups need just one fixed binary key.
                if(initialized && previous_giant!=group.giant_step){
                    rotate(*rotated,result.ct,plan.k_step*(previous_giant-group.giant_step));result.ct=*rotated;
                }
                if(!initialized){result.ct=*inner;initialized=true;}else CKKSAddTRLWEInPlace<P,output_q>(result.ct,*inner);
                previous_giant=group.giant_step;
            }else{
                const auto *term=inner.get();if(group.giant_step){rotate(*rotated,*inner,plan.k_step*group.giant_step);term=rotated.get();}
                if(!initialized){result.ct=*term;initialized=true;}else CKKSAddTRLWEInPlace<P,output_q>(result.ct,*term);
            }
            if(timings)timings->giant_rotations+=std::chrono::duration<double>(Clock::now()-start).count();
        }
        if(horner_giants && initialized && previous_giant){
            start=Clock::now();rotate(*rotated,result.ct,plan.k_step*previous_giant);result.ct=*rotated;
            if(timings)timings->giant_rotations+=std::chrono::duration<double>(Clock::now()-start).count();
        }
    }
public:
    // Horner preserves per-group rounding but changes switching noise.
    // Dense giant steps reuse one binary key; sparse plans remain opt-in.
    template<class InputKey,class OutputKey>
    void apply(Result&result,const Input&input,const InputKey&ik,const OutputKey&ok,CKKSLinearRNSTimings*timings=nullptr,bool fused_terms=true,bool horner_giants=false)const {
        if(!has_direct&&!has_conjugate){for(auto &component:result.ct)component.fill({});return;}
        auto term=std::make_unique<Result>();if(has_direct){applyPart(*term,input,direct,ik,ok,timings,fused_terms,horner_giants);result=*term;}
        if(has_conjugate){auto ct=std::make_unique<Input>();CKKSConjugateSlots<P,Q>(ct->ct,input.ct,ik);applyPart(*term,*ct,conjugate,ik,ok,timings,fused_terms,horner_giants);if(!has_direct)result=*term;else CKKSAddTRLWEInPlace<P,output_q>(result.ct,term->ct);}
    }
};
template<class P,std::uint32_t Q,std::uint32_t Delta,std::uint32_t PlainDelta,std::size_t Count,class InputKey,class OutputKey>
inline void CKKSRealLinearTransform(CKKSPlainMulResult<P,Q,Delta,PlainDelta>&result,const CKKSCiphertext<P,Q,Delta>&input,
    const CKKSPreparedCompactRealLinearTransformRNS32<P,Q,Delta,PlainDelta,Count>&plan,const InputKey&ik,const OutputKey&ok,CKKSLinearRNSTimings*timings=nullptr,bool horner_giants=false)
{plan.apply(result,input,ik,ok,timings,true,horner_giants);}
} // namespace TFHEpp
#endif
