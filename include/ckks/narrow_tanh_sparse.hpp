#pragma once
#include "narrow_approximation.hpp"
#include <cstdint>

namespace TFHEpp {
// Offline bounded-L1 fit on [-4,4], period16. Integer imaginary weights
// already match the evaluator's 28-bit coefficient precision, so there is no
// additional approximation-to-evaluator coefficient rounding.
// Reproduce the continuous-domain certificate with sparse_tanh_certificate.cpp.
inline CKKSNarrowTrigonometricApproximation CKKSBuildNarrowSparseTanhApproximation()
{
    constexpr std::array<unsigned,20> frequencies{1,2,3,5,6,7,8,9,10,14,15,16,17,18,19,20,22,24,25,30};
    constexpr std::array<std::int64_t,20> weights{-159877184,-761407,-33621499,-9306132,-1415660,-381343,-2655022,1173886,-977945,-148264,85245,3444,-83323,81208,-47670,13393,-1773,406,-221,-4};
    CKKSNarrowTrigonometricApproximation result;
    result.coefficients.resize(31);
    for(std::size_t i=0;i<frequencies.size();++i)
        result.coefficients[frequencies[i]]={0,std::ldexp(static_cast<double>(weights[i]),-28)};
    result.message_bound=4;result.gap=16;
    // Conservative rounded envelope, including the certificate's between-grid
    // curvature bound and floating-point allowance; not an encrypted guarantee.
    result.max_sampled_error=1.041e-8;
    result.uniform_error_bound=1.06e-8;
    return result;
}
} // namespace TFHEpp
