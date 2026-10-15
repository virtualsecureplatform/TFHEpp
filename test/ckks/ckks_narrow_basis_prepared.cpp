#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_basis_prepared.hpp>
#include <ckks/prepared_dd_ntt.hpp>
#include <iostream>

using namespace TFHEpp;
using P = narrow_benchmark::BenchmarkParam<9,7>;
using Narrow = CKKSNarrowBootstrapSchedule<P,48,357,35,35,31,true,13,true>;
using Basis = CKKSNarrowFrequencyBasisSchedule<Narrow,32>;

int main()
{
    CKKSNarrowFrequencyBasis plan(12,{1,3,5,6});
    std::vector<Basis::Root> roots(plan.frequencies().size());
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    for (std::size_t j = 0; j < roots.size(); ++j) {
        for (std::size_t i = 0; i < P::n/2; ++i)
            (*slots)[i] = std::polar(1.0,plan.frequencies()[j]*(0.3+0.7*i/(P::n/2)));
        ckksSlotEncode<P,Basis::root_log_q,Basis::root_log_delta>(roots[j].ct[1],*slots);
    }
    // Transparent ciphertexts with a zero secret need no encrypted evaluation
    // keys. This isolates reconstruction and rounding without generating keys.
    auto raw = std::make_unique<CKKSNarrowFrequencyBasisKey<Basis>>();
    auto fft = makeCKKSDDRelinKeyFFT(raw->product);
    auto auto_fft = makeCKKSPreparedAutoKeyFFT<P,Basis::output_log_q>(raw->conjugation);
    auto compact = std::make_unique<CKKSCompactDDRelinKey<P,247,64,32>>();
    auto ntt = makeCKKSDDRelinKeyNTT(*compact);
    CKKSAuxiliaryDDRelinKeyNTT<P,147,100> auxiliary(ntt);
    std::vector<std::complex<double>> coefficients(13);
    coefficients[0] = {0.5,0};
    for (std::size_t k = 1; k < coefficients.size(); ++k)
        coefficients[k] = {0.015*k/(k+1),-0.01*k/(k+1)};
    auto reference = std::make_unique<Basis::Output>();
    auto prepared = std::make_unique<Basis::Output>();
    auto integer = std::make_unique<Basis::Output>();
    CKKSNarrowEvalSeriesFromBasisRoots<Basis>(*reference,roots,*raw,plan,coefficients);
    const auto conjugate = [&](auto &result,const auto &input) {
        CKKSEvalAuto<P,Basis::output_log_q>(result,input,2*P::n-1,*auto_fft);
    };
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*prepared,roots,plan,coefficients,*fft,conjugate);
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*integer,roots,plan,coefficients,auxiliary,conjugate);
    if (reference->ct != prepared->ct || reference->ct != integer->ct)
        throw std::runtime_error("prepared basis backends differ from existing reconstruction");
    ckksSlotDecode<P,Basis::output_log_q,Basis::log_delta>(*slots,integer->ct[1]);
    double maximum = 0;
    for (std::size_t i = 0; i < P::n/2; ++i) {
        double expected = coefficients[0].real();
        for (std::size_t k = 1; k < coefficients.size(); ++k)
            expected += 2*(coefficients[k]*std::polar(1.0,k*(0.3+0.7*i/(P::n/2)))).real();
        maximum = std::max(maximum,std::abs((*slots)[i]-expected));
    }
    if (maximum >= std::ldexp(1.0,-23)) throw std::runtime_error("independent basis series gate failed");
    std::cout << "prepared FFT and auxiliary NTT basis byte-identical; series error=" << maximum << '\n';
}
