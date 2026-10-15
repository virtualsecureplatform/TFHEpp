#pragma once
#include "narrow_frequency_basis.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
// Section 7 reconstruction with explicitly supplied prepared public backends.
// No legacy key placeholders, encryption, or cache construction during use.
template <class Schedule, class ProductKey, class Conjugate>
inline void CKKSNarrowEvalSeriesFromBasisRootsPrepared(
    typename Schedule::Output &result,
    std::span<const typename Schedule::Root> roots,
    const CKKSNarrowFrequencyBasis &plan,
    const std::vector<std::complex<double>> &coefficients,
    const ProductKey &product_key, Conjugate conjugate)
{
    using P = typename Schedule::Param;
    static_assert(ProductKey::log_q == Schedule::root_log_q);
    if (roots.size() != plan.frequencies().size() || coefficients.size() != plan.degree()+1 ||
        coefficients[0].imag() != 0)
        throw std::invalid_argument("invalid prepared narrow basis roots or coefficients");
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
        else CKKSMultRelinThenRescale<P>(*product,roots[expression[0]],roots[expression[1]],product_key);
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
