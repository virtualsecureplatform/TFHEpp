#pragma once
#include "narrow_frequency_basis.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
// Section 7 reconstruction with explicitly supplied prepared public backends.
// No legacy key placeholders, encryption, or cache construction during use.
template <class Schedule, class ProductKey, class Conjugate, class RootConjugate=std::nullptr_t>
inline void CKKSNarrowEvalSeriesFromBasisRootsPrepared(
    typename Schedule::Output &result,
    std::span<const typename Schedule::Root> roots,
    const CKKSNarrowFrequencyBasis &plan,
    const std::vector<std::complex<double>> &coefficients,
    const ProductKey &product_key, Conjugate conjugate, RootConjugate root_conjugate=nullptr)
{
    using P = typename Schedule::Param;
    static_assert(ProductKey::log_q == Schedule::root_log_q);
    if (roots.size() != plan.frequencies().size() || coefficients.size() != plan.degree()+1 ||
        coefficients[0].imag() != 0)
        throw std::invalid_argument("invalid prepared narrow basis roots or coefficients");
    if constexpr(std::is_same_v<RootConjugate,std::nullptr_t>)
        if(plan.usesDifferences())throw std::invalid_argument("frequency differences require a root conjugation backend");
    plan.validateCoefficients(coefficients);
    std::vector<std::unique_ptr<typename Schedule::Root>> negative_roots(roots.size());
    std::vector<std::array<std::int64_t,2>> weights;
    for (const auto coefficient : coefficients)
        weights.push_back({narrow_detail::frequencyCoefficient<Schedule::coefficient_bits>(coefficient.real()),
                          narrow_detail::frequencyCoefficient<Schedule::coefficient_bits>(coefficient.imag())});
    for (auto &component : result.ct) component.fill(0);
    auto product = std::make_unique<typename Schedule::Product>();
    for (std::size_t k = 1; k <= plan.degree(); ++k) {
        if (!weights[k][0] && !weights[k][1]) continue;
        const auto expression = plan.expression(k);
        if (expression[1] == roots.size()) narrow_detail::liftFrequencyRoot<Schedule>(*product,roots[expression[0]]);
        else {
            const auto *rhs=&roots[expression[1]];
            if(plan.conjugatesSecond(k)){
                if constexpr(!std::is_same_v<RootConjugate,std::nullptr_t>){
                    auto &negative=negative_roots[expression[1]];
                    if(!negative){negative=std::make_unique<typename Schedule::Root>();root_conjugate(negative->ct,rhs->ct);}
                    rhs=negative.get();
                }
            }
            CKKSMultRelinThenRescaleBy<P, Schedule::product_rescale_bits>(*product,roots[expression[0]],*rhs,product_key);
        }
        narrow_detail::addFrequencyWeight<Schedule>(result,*product,weights[k][0],weights[k][1]);
    }
    auto negative = std::make_unique<typename Schedule::Output>();
    conjugate(negative->ct,result.ct);
    CKKSAddInPlace<P>(result,*negative);
    result.ct[P::k][0] = ckks_detail::reduceToLevel<P,Schedule::output_log_q>(
        result.ct[P::k][0]+(P::T::from_signed_i64(weights[0][0]) << Schedule::product_log_delta));
}
// Weight exact modular tensors first, then relinearize and rescale once.
// The larger intermediate scale must fit at the root modulus. This changes
// rounding/key-switch noise relative to termwise reconstruction, so callers
// must qualify encrypted precision rather than expect identical ciphertexts.
template <class Schedule, class ProductKey, class Conjugate, class RootConjugate=std::nullptr_t>
inline void CKKSNarrowEvalSeriesFromBasisRootsAccumulatedPrepared(
    typename Schedule::Output &result,
    std::span<const typename Schedule::Root> roots,
    const CKKSNarrowFrequencyBasis &plan,
    const std::vector<std::complex<double>> &coefficients,
    const ProductKey &product_key, Conjugate conjugate, RootConjugate root_conjugate=nullptr,
    bool group_tensors=false)
{
    using P=typename Schedule::Param;
    constexpr auto Q=Schedule::root_log_q;
    constexpr auto delta=2*Schedule::root_log_delta+Schedule::coefficient_bits;
    static_assert(P::k==1 && ProductKey::log_q==Q);
    static_assert(delta<Q,"weighted unscaled tensor needs modulus headroom");
    if(roots.size()!=plan.frequencies().size() || coefficients.size()!=plan.degree()+1 || coefficients[0].imag()!=0)
        throw std::invalid_argument("invalid accumulated narrow basis roots or coefficients");
    if constexpr(std::is_same_v<RootConjugate,std::nullptr_t>)
        if(plan.usesDifferences())throw std::invalid_argument("frequency differences require a root conjugation backend");
    plan.validateCoefficients(coefficients);
    std::vector<std::array<std::int64_t,2>> weights;
    for(const auto coefficient:coefficients)
        weights.push_back({narrow_detail::frequencyCoefficient<Schedule::coefficient_bits>(coefficient.real()),
                          narrow_detail::frequencyCoefficient<Schedule::coefficient_bits>(coefficient.imag())});
    auto accumulated=std::make_unique<TRLWE3<P>>();
    auto tensor=std::make_unique<TRLWE3<P>>();
    std::vector<std::unique_ptr<typename Schedule::Root>> negative_roots(roots.size());
    // Tensor multiplication is bilinear over Z_Q[X]/(X^N+1). The integer
    // complex weights are ring elements real + imag*X^(N/2), so grouping
    // weighted right operands preserves the entire tensor before switching.
    std::vector<std::unique_ptr<typename Schedule::Root>> grouped_rhs(group_tensors?roots.size():0);
    for(std::size_t k=1;k<=plan.degree();++k){
        if(!weights[k][0] && !weights[k][1])continue;
        const auto expression=plan.expression(k);
        if(expression[1]==roots.size()){
            for(auto&component:*tensor)component.fill(0);
            for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)
                (*tensor)[c][i]=ckks_detail::reduceToLevel<P,Q>(roots[expression[0]].ct[c][i]<<Schedule::root_log_delta);
        }else{
            const auto*rhs=&roots[expression[1]];
            if(plan.conjugatesSecond(k)){
                if constexpr(!std::is_same_v<RootConjugate,std::nullptr_t>){
                    auto&negative=negative_roots[expression[1]];
                    if(!negative){negative=std::make_unique<typename Schedule::Root>();root_conjugate(negative->ct,rhs->ct);}
                    rhs=negative.get();
                }
            }
            if(group_tensors){
                auto&right=grouped_rhs[expression[0]];
                if(!right)right=std::make_unique<typename Schedule::Root>();
                for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i){
                    const auto imaginary=i<P::n/2?-rhs->ct[c][i+P::n/2]:rhs->ct[c][i-P::n/2];
                    right->ct[c][i]=ckks_detail::reduceToLevel<P,Q>(right->ct[c][i]+rhs->ct[c][i]*weights[k][0]+imaginary*weights[k][1]);
                }
                continue;
            }
            if constexpr(requires{product_key.tensorProduct(*tensor,roots[expression[0]].ct,rhs->ct);})
                product_key.tensorProduct(*tensor,roots[expression[0]].ct,rhs->ct);
            else CKKSTensorProductRescale<P,Q,Q,0>(*tensor,roots[expression[0]].ct,rhs->ct);
        }
        for(std::size_t c=0;c<3;++c)for(std::size_t i=0;i<P::n;++i){
            const auto imaginary=i<P::n/2?-(*tensor)[c][i+P::n/2]:(*tensor)[c][i-P::n/2];
            (*accumulated)[c][i]=ckks_detail::reduceToLevel<P,Q>((*accumulated)[c][i]+(*tensor)[c][i]*weights[k][0]+imaginary*weights[k][1]);
        }
    }
    for(std::size_t j=0;j<grouped_rhs.size();++j)if(grouped_rhs[j]){
        if constexpr(requires{product_key.tensorProduct(*tensor,roots[j].ct,grouped_rhs[j]->ct);})
            product_key.tensorProduct(*tensor,roots[j].ct,grouped_rhs[j]->ct);
        else CKKSTensorProductRescale<P,Q,Q,0>(*tensor,roots[j].ct,grouped_rhs[j]->ct);
        for(std::size_t c=0;c<3;++c)for(std::size_t i=0;i<P::n;++i)
            (*accumulated)[c][i]=ckks_detail::reduceToLevel<P,Q>((*accumulated)[c][i]+(*tensor)[c][i]);
    }
    auto unscaled=std::make_unique<CKKSCiphertext<P,Q,delta>>();
    CKKSRelinearization<P,Q>(unscaled->ct,*accumulated,product_key);
    CKKSScaleDown<P,Q,delta,Schedule::product_rescale_bits>(result,*unscaled);
    auto negative=std::make_unique<typename Schedule::Output>();
    conjugate(negative->ct,result.ct);CKKSAddInPlace<P>(result,*negative);
    result.ct[P::k][0]=ckks_detail::reduceToLevel<P,Schedule::output_log_q>(
        result.ct[P::k][0]+(P::T::from_signed_i64(weights[0][0])<<Schedule::product_log_delta));
}
} // namespace TFHEpp
#endif
