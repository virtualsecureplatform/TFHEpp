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
        if (expression[1] == roots.size()) CKKSLevelReduce<P>(*product,roots[expression[0]]);
        else {
            const auto *rhs=&roots[expression[1]];
            if(plan.conjugatesSecond(k)){
                if constexpr(!std::is_same_v<RootConjugate,std::nullptr_t>){
                    auto &negative=negative_roots[expression[1]];
                    if(!negative){negative=std::make_unique<typename Schedule::Root>();root_conjugate(negative->ct,rhs->ct);}
                    rhs=negative.get();
                }
            }
            CKKSMultRelinThenRescale<P>(*product,roots[expression[0]],*rhs,product_key);
        }
        narrow_detail::addFrequencyWeight<Schedule>(result,*product,weights[k][0],weights[k][1]);
    }
    auto negative = std::make_unique<typename Schedule::Output>();
    conjugate(negative->ct,result.ct);
    CKKSAddInPlace<P>(result,*negative);
    result.ct[P::k][0] = ckks_detail::reduceToLevel<P,Schedule::output_log_q>(
        result.ct[P::k][0]+(P::T::from_signed_i64(weights[0][0]) << Schedule::root_log_delta));
}
} // namespace TFHEpp
#endif
