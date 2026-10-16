// No keys or ciphertexts. Build as in tools/frequency_budget.cpp.
#include <ckks/narrow_tanh_sparse.hpp>
#include <iostream>
#include <iomanip>
int main()
{
    const auto a=TFHEpp::CKKSBuildNarrowSparseTanhApproximation();
    constexpr std::size_t intervals=1U<<20;
    const long double omega=2*std::acos(-1.L)/a.gap;
    long double l1=0,curvature=2,sampled=0;
    for(std::size_t k=1;k<a.coefficients.size();++k){
        const long double amplitude=2*std::abs(a.coefficients[k]);
        l1+=amplitude;curvature+=amplitude*omega*omega*k*k;
    }
    // Odd target and series: [0,4] suffices for the declared [-4,4] domain.
    for(std::size_t i=0;i<=intervals;++i){
        const long double x=4.L*i/intervals;long double value=0;
        for(std::size_t k=1;k<a.coefficients.size();++k)
            value-=2*static_cast<long double>(a.coefficients[k].imag())*std::sin(omega*k*x);
        sampled=std::max(sampled,std::abs(value-std::tanh(x)));
    }
    // Linear interpolation remainder |r''|h^2/8, |tanh''| <=2.
    const long double h=4.L/intervals;
    const long double rounding=1e-12L; // Deliberately exceeds floating evaluation error.
    const long double uniform=sampled+curvature*h*h/8+rounding;
    std::cout<<std::setprecision(17)<<"{\"intervals\":"<<intervals<<",\"sampled\":"<<sampled<<",\"curvature_bound\":"<<curvature<<",\"between_grid_bound\":"<<curvature*h*h/8<<",\"floating_allowance\":"<<rounding<<",\"uniform_bound\":"<<uniform<<",\"published_bound\":"<<*a.uniform_error_bound<<",\"l1\":"<<l1<<"}\n";
    return uniform<*a.uniform_error_bound?0:1;
}
