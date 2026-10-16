#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_basis_prepared.hpp>
#include <ckks/prepared_dd_ntt.hpp>
#include <ckks/tensor_rns32.hpp>
#include <ckks/narrow_two_layer_basis.hpp>
#include <ckks/key_projection.hpp>
#include <ckks/narrow_tanh_sparse.hpp>
#include <iostream>

using namespace TFHEpp;
using P = narrow_benchmark::BenchmarkParam<9,7>;
using Narrow = CKKSNarrowBootstrapSchedule<P,48,357,35,35,31,true,13,true>;
using Basis = CKKSNarrowFrequencyBasisSchedule<Narrow,32>;

template<class BP>struct TwoLayerMiddle {using Param=BP;static constexpr std::uint32_t output_log_q=156,log_delta=37;using Output=CKKSCiphertext<BP,156,37>;};
template<std::uint32_t Nbit>void checkTwoLayerBasis()
{
    using BP=narrow_benchmark::BenchmarkParam<Nbit,7>;
    using Seed=CKKSCiphertext<BP,193,37>;
    using Middle=TwoLayerMiddle<BP>;
    using Series=CKKSNarrowFrequencyBasisSchedule<Middle,28>;
    static_assert(Series::output_log_q==119 && Series::log_delta==65);
    static_assert(119-28+(44-37)>=93); // Final93/44 and subsequent49/44 fit.
    CKKSNarrowTwoLayerFrequencyBasis plan(40,{6,9,10},{3,6,9,12,15,18,19,20});
    std::vector<Seed> seeds(3);auto slots=std::make_unique<CKKSSlotVector<BP>>();
    const auto angle=[&](std::size_t i){return -std::acos(-1.0)/2+std::acos(-1.0)*i/(BP::n/2-1);};
    for(std::size_t j=0;j<seeds.size();++j){
        for(std::size_t i=0;i<BP::n/2;++i)(*slots)[i]=std::polar(1.0,plan.seeds()[j]*angle(i));
        ckksSlotEncode<BP,193,37>(seeds[j].ct[1],*slots);
    }
    // Transparent rows exercise arithmetic without generating encrypted keys.
    auto first_raw=std::make_unique<CKKSCompactDDRelinKey<BP,241,20,20>>();
    auto second_raw=std::make_unique<CKKSCompactDDRelinKey<BP,204,20,20>>();
    CKKSAuxiliaryDDRelinKeyNTT<BP,193,48,20,20> first(makeCKKSDDRelinKeyNTT<20>(*first_raw));
    CKKSAuxiliaryDDRelinKeyNTT<BP,156,48,20,20> second(makeCKKSDDRelinKeyNTT<20>(*second_raw));
    CKKSRelinKeyTensorRNS32<BP,decltype(first)> first_rns{first};
    CKKSRelinKeyTensorRNS32<BP,decltype(second)> second_rns{second};
    std::vector<typename Series::Root> intermediate;
    std::size_t conjugations=0;
    CKKSNarrowDeriveFrequencyRoots<BP,193,37>(intermediate,seeds,plan,first_rns,
        [&](auto&out,const auto&in){++conjugations;for(std::size_t c=0;c<2;++c)Automorphism<BP>(out[c],in[c],2*BP::n-1);ckks_detail::reduceTRLWEToLevel<BP,193>(out);});
    if(conjugations!=1 || intermediate.size()!=8)throw std::runtime_error("two-layer root construction shape differs");
    std::vector<std::complex<double>> coefficients(41);coefficients[0]={0.5,0};
    for(std::size_t k=1;k<coefficients.size();++k)coefficients[k]={0.015*k/(k+1),-0.01*k/(k+1)};
    auto weighted=std::make_unique<typename Series::Output>();
    const auto conjugate=[&](auto&out,const auto&in){for(std::size_t c=0;c<2;++c)Automorphism<BP>(out[c],in[c],2*BP::n-1);ckks_detail::reduceTRLWEToLevel<BP,119>(out);};
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Series>(*weighted,intermediate,plan.finalPlan(),coefficients,second_rns,conjugate,
        [&](auto&out,const auto&in){for(std::size_t c=0;c<2;++c)Automorphism<BP>(out[c],in[c],2*BP::n-1);ckks_detail::reduceTRLWEToLevel<BP,156>(out);});
    auto reduced=std::make_unique<CKKSCiphertext<BP,91,37>>();CKKSScaleDown<BP,119,65,28>(*reduced,*weighted);
    auto output=std::make_unique<CKKSCiphertext<BP,93,44>>();
    for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<BP::n;++i)output->ct[c][i]=ckks_detail::reduceToLevel<BP,93>(reduced->ct[c][i]<<7);
    ckksSlotDecode<BP,93,44>(*slots,output->ct[1]);double maximum=0;
    for(std::size_t i=0;i<BP::n/2;++i){double expected=coefficients[0].real();for(std::size_t k=1;k<=40;++k)expected+=2*(coefficients[k]*std::polar(1.0,k*angle(i))).real();maximum=std::max(maximum,std::abs((*slots)[i]-expected));}
    if(maximum>=std::ldexp(1.0,-23))throw std::runtime_error("two-layer37-bit independent series gate failed");
    auto post_raw=makeCKKSCompactDDRelinKeyAtLevel<141>(*first_raw);
    CKKSAuxiliaryDDRelinKeyNTT<BP,93,48,20,20> post(makeCKKSDDRelinKeyNTT<20>(*post_raw));
    CKKSRelinKeyTensorRNS32<BP,decltype(post)> post_rns{post};
    auto square=std::make_unique<CKKSMultResult<BP,93,44,93,44>>();CKKSMultRelinThenRescale<BP>(*square,*output,*output,post_rns);
    ckksSlotDecode<BP,49,44>(*slots,square->ct[1]);double square_error=0;
    for(std::size_t i=0;i<BP::n/2;++i){double expected=coefficients[0].real();for(std::size_t k=1;k<=40;++k)expected+=2*(coefficients[k]*std::polar(1.0,k*angle(i))).real();square_error=std::max(square_error,std::abs((*slots)[i]-expected*expected));}
    if(square_error>=std::ldexp(1.0,-23))throw std::runtime_error("two-layer remaining-level square gate failed");
    bool rejected=false;try{CKKSNarrowTwoLayerFrequencyBasis invalid(40,{6,9,10},{3,6,9,12,15,18,19,21});}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("uncovered two-layer plan accepted");
    std::cout<<"two-layer basis N="<<BP::n<<" seeds=3 derived=8 delta=37 output=93/44 independent error="<<maximum<<" square49/44 error="<<square_error<<'\n';
}

namespace TFHEpp {
template<class BP,class Key>struct CountedBasisProduct {
    static constexpr auto log_q=Key::log_q;
    const Key&key;std::size_t&switches;std::size_t&tensors;
    void tensorProduct(TRLWE3<BP>&out,const TRLWE<BP>&a,const TRLWE<BP>&b)const
    {++tensors;key.tensorProduct(out,a,b);}
};
namespace ckks_detail {
template<class BP,class Key>struct is_ckks_dd_relin_key<CountedBasisProduct<BP,Key>> : std::true_type {};
}
template<class BP,class Key>void CKKSRelinKeySwitch(TRLWE<BP>&out,const Polynomial<BP>&in,const CountedBasisProduct<BP,Key>&key)
{++key.switches;CKKSRelinKeySwitch<BP>(out,in,key.key);}
}

template<class BP>struct SparseMiddle {using Param=BP;static constexpr std::uint32_t output_log_q=165,log_delta=44;using Output=CKKSCiphertext<BP,165,44>;};
template<std::uint32_t Nbit>void checkSparseTanh()
{
    using BP=narrow_benchmark::BenchmarkParam<Nbit,7>;
    using Series=CKKSNarrowFrequencyBasisSchedule<SparseMiddle<BP>,28>;
    auto plan=CKKSNarrowFrequencyBasis::forSparseSeries(30,{7,9,10,15});
    auto approximation=CKKSBuildNarrowSparseTanhApproximation();
    plan.validateCoefficients(approximation.coefficients);
    if(plan.covers(4) || !plan.covers(30))throw std::runtime_error("sparse frequency coverage differs");
    std::vector<typename Series::Root> roots(4);auto slots=std::make_unique<CKKSSlotVector<BP>>();
    const auto x=[&](std::size_t i){return -4+8.0*i/(BP::n/2-1);};
    for(std::size_t j=0;j<roots.size();++j){
        for(std::size_t i=0;i<BP::n/2;++i)(*slots)[i]=std::polar(1.0,2*std::acos(-1.0)*x(i)*plan.frequencies()[j]/16);
        ckksSlotEncode<BP,165,44>(roots[j].ct[1],*slots);
    }
    auto raw=std::make_unique<CKKSCompactDDRelinKey<BP,213,20,20>>();
    CKKSAuxiliaryDDRelinKeyNTT<BP,165,48,20,20> product(makeCKKSDDRelinKeyNTT<20>(*raw));
    CKKSRelinKeyTensorRNS32<BP,decltype(product)> rns{product};
    auto weighted=std::make_unique<typename Series::Output>();
    const auto conjugate=[](auto&out,const auto&in){for(std::size_t c=0;c<2;++c)Automorphism<BP>(out[c],in[c],2*BP::n-1);ckks_detail::reduceTRLWEToLevel<BP,121>(out);};
    const auto root_conjugate=[](auto&out,const auto&in){for(std::size_t c=0;c<2;++c)Automorphism<BP>(out[c],in[c],2*BP::n-1);ckks_detail::reduceTRLWEToLevel<BP,165>(out);};
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Series>(*weighted,roots,plan,approximation.coefficients,rns,conjugate,root_conjugate);
    auto fft=std::make_unique<typename Series::Output>();
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Series>(*fft,roots,plan,approximation.coefficients,product,conjugate,root_conjugate);
    if(fft->ct!=weighted->ct)throw std::runtime_error("sparse RNS and FFT basis tensors differ");
    auto bad=approximation.coefficients;bad[4]={1e-300,0};bool rejected=false;
    try{CKKSNarrowEvalSeriesFromBasisRootsPrepared<Series>(*weighted,roots,plan,bad,rns,conjugate,root_conjugate);}catch(const std::invalid_argument&){rejected=true;}
    if(!rejected || fft->ct!=weighted->ct)throw std::runtime_error("unsupported sparse coefficient was not rejected before output mutation");
    std::size_t switches=0,tensors=0;
    CountedBasisProduct<BP,decltype(rns)> counted{rns,switches,tensors};
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Series>(*weighted,roots,plan,approximation.coefficients,counted,conjugate,root_conjugate);
    if(switches!=16 || tensors!=16)throw std::runtime_error("termwise sparse reconstruction product count differs");
    switches=tensors=0;
    auto accumulated=std::make_unique<typename Series::Output>();
    CKKSNarrowEvalSeriesFromBasisRootsAccumulatedPrepared<Series>(*accumulated,roots,plan,approximation.coefficients,counted,conjugate,root_conjugate);
    if(switches!=1 || tensors!=16)throw std::runtime_error("accumulated basis did not use exactly one key switch");
    CKKSNarrowEvalSeriesFromBasisRootsAccumulatedPrepared<Series>(*fft,roots,plan,approximation.coefficients,product,conjugate,root_conjugate);
    if(fft->ct!=accumulated->ct)throw std::runtime_error("accumulated sparse RNS and FFT tensors differ");
    CKKSNarrowEvalSeriesFromBasisRootsAccumulatedPrepared<Series>(*weighted,roots,plan,approximation.coefficients,rns,conjugate,root_conjugate);
    auto output=std::make_unique<CKKSCiphertext<BP,93,44>>();CKKSScaleDown<BP,121,72,28>(*output,*weighted);
    ckksSlotDecode<BP,93,44>(*slots,output->ct[1]);double error=0;
    for(std::size_t i=0;i<BP::n/2;++i)error=std::max(error,std::abs((*slots)[i]-std::tanh(x(i))));
    auto postraw=makeCKKSCompactDDRelinKeyAtLevel<141>(*raw);
    CKKSAuxiliaryDDRelinKeyNTT<BP,93,48,20,20> post(makeCKKSDDRelinKeyNTT<20>(*postraw));
    auto square=std::make_unique<CKKSMultResult<BP,93,44,93,44>>();CKKSMultRelinThenRescale<BP>(*square,*output,*output,post);
    ckksSlotDecode<BP,49,44>(*slots,square->ct[1]);double square_error=0;
    for(std::size_t i=0;i<BP::n/2;++i)square_error=std::max(square_error,std::abs((*slots)[i]-std::pow(std::tanh(x(i)),2)));
    if(error>=std::ldexp(1.0,-23) || square_error>=3*std::ldexp(1.0,-23))throw std::runtime_error("sparse Tanh independent function or square failed");
    std::cout<<"sparse Tanh N="<<BP::n<<" roots=4 delta=44 output=93/44 error="<<error<<" square49/44 error="<<square_error<<'\n';
}

int main(int argc,char**argv)
{
    if(argc>1){if(argc!=2)throw std::invalid_argument("expected one full-ring option");if(std::string(argv[1])=="--two-layer-full")checkTwoLayerBasis<14>();else if(std::string(argv[1])=="--sparse-full")checkSparseTanh<14>();else throw std::invalid_argument("unknown full-ring option");return 0;}
    checkSparseTanh<9>();
    checkTwoLayerBasis<9>();
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
    CKKSRelinKeyTensorRNS32<P,decltype(auxiliary)> rns_product{auxiliary};
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
    auto rns=std::make_unique<Basis::Output>();
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*rns,roots,plan,coefficients,rns_product,conjugate);
    if(rns->ct!=integer->ct)throw std::runtime_error("RNS basis tensor differs from FFT tensor");
    if (reference->ct != prepared->ct || reference->ct != integer->ct)
        throw std::runtime_error("prepared basis backends differ from existing reconstruction");
    auto accumulated_dense=std::make_unique<Basis::Output>();
    CKKSNarrowEvalSeriesFromBasisRootsAccumulatedPrepared<Basis>(*accumulated_dense,roots,plan,coefficients,rns_product,conjugate);
    auto accumulated_slots=std::make_unique<CKKSSlotVector<P>>();
    ckksSlotDecode<P,Basis::output_log_q,Basis::log_delta>(*accumulated_slots,accumulated_dense->ct[1]);
    ckksSlotDecode<P,Basis::output_log_q,Basis::log_delta>(*slots,integer->ct[1]);
    double maximum = 0;
    for (std::size_t i = 0; i < P::n/2; ++i) {
        double expected = coefficients[0].real();
        for (std::size_t k = 1; k < coefficients.size(); ++k)
            expected += 2*(coefficients[k]*std::polar(1.0,k*(0.3+0.7*i/(P::n/2)))).real();
        maximum = std::max({maximum,std::abs((*slots)[i]-expected),std::abs((*accumulated_slots)[i]-expected)});
    }
    if (maximum >= std::ldexp(1.0,-23)) throw std::runtime_error("independent basis series gate failed");
    {
        CKKSNarrowFrequencyBasis three(10,{2,4,5},true),four(10,{1,3,5,6});
        std::vector<Basis::Root> small_roots(3),control_roots(4);
        const auto angle=[&](std::size_t i){return -std::acos(-1.0)/2+std::acos(-1.0)*i/(P::n/2-1);};
        for(auto pair:{std::pair{&three,&small_roots},std::pair{&four,&control_roots}})
            for(std::size_t j=0;j<pair.second->size();++j){
                for(std::size_t i=0;i<P::n/2;++i)(*slots)[i]=std::polar(1.0,pair.first->frequencies()[j]*angle(i));
                ckksSlotEncode<P,Basis::root_log_q,Basis::root_log_delta>((*pair.second)[j].ct[1],*slots);
            }
        auto small_coefficients=coefficients;small_coefficients.resize(11);
        std::size_t conjugations=0;
        CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*integer,small_roots,three,small_coefficients,auxiliary,conjugate,
            [&](auto&output,const auto&input){
                ++conjugations;
                for(std::size_t c=0;c<2;++c)Automorphism<P>(output[c],input[c],2*P::n-1);
                ckks_detail::reduceTRLWEToLevel<P,Basis::root_log_q>(output);
            });
        if(conjugations!=2)throw std::runtime_error("three-root conjugations were not shared");
        CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*prepared,control_roots,four,small_coefficients,auxiliary,conjugate);
        auto control=std::make_unique<CKKSSlotVector<P>>();
        ckksSlotDecode<P,Basis::output_log_q,Basis::log_delta>(*slots,integer->ct[1]);
        ckksSlotDecode<P,Basis::output_log_q,Basis::log_delta>(*control,prepared->ct[1]);
        double error=0;
        for(std::size_t i=0;i<P::n/2;++i){
            double expected=small_coefficients[0].real();
            for(std::size_t k=1;k<=10;++k)expected+=2*(small_coefficients[k]*std::polar(1.0,k*angle(i))).real();
            error=std::max({error,std::abs((*slots)[i]-expected),std::abs((*control)[i]-expected),std::abs((*slots)[i]-(*control)[i])});
        }
        if(error>=std::ldexp(1.0,-23))throw std::runtime_error("three-root reconstruction differs from independent series or four-root control");
        std::cout<<"three-root degree10 independent series error="<<error<<'\n';
    }
    CKKSNarrowFrequencyBasis signed_plan(40,{6,12,15,16,17,19,20},true);
    if(!signed_plan.usesDifferences())throw std::runtime_error("difference plan was not recognized");
    bool rejected=false;
    try{CKKSNarrowFrequencyBasis invalid(40,{6,12,15,16,17,19,20});}
    catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("unsigned gap accepted");
    std::vector<Basis::Root> signed_roots(signed_plan.frequencies().size());
    for(std::size_t j=0;j<signed_roots.size();++j){
        for(std::size_t i=0;i<P::n/2;++i)(*slots)[i]=std::polar(1.0,signed_plan.frequencies()[j]*(0.3+0.7*i/(P::n/2)));
        ckksSlotEncode<P,Basis::root_log_q,Basis::root_log_delta>(signed_roots[j].ct[1],*slots);
    }
    std::vector<std::complex<double>> signed_coefficients(41);signed_coefficients[0]={0.5,0};
    for(std::size_t k=1;k<=40;++k){
        signed_coefficients[k]={0.015*k/(k+1),-0.01*k/(k+1)};
        const auto expression=signed_plan.expression(k);
        const auto a=signed_plan.frequencies()[expression[0]];
        const auto value=expression[1]==signed_roots.size()?a:signed_plan.conjugatesSecond(k)?a-signed_plan.frequencies()[expression[1]]:a+signed_plan.frequencies()[expression[1]];
        if(value!=k)throw std::runtime_error("signed frequency expression is incorrect");
    }
    rejected=false;
    try{CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*integer,signed_roots,signed_plan,signed_coefficients,auxiliary,conjugate);}
    catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("missing root conjugation backend accepted");
    std::size_t root_conjugations=0;
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*integer,signed_roots,signed_plan,signed_coefficients,auxiliary,conjugate,
        [&](auto&output,const auto&input){
            ++root_conjugations;
            for(std::size_t c=0;c<2;++c)Automorphism<P>(output[c],input[c],2*P::n-1);
            ckks_detail::reduceTRLWEToLevel<P,Basis::root_log_q>(output);
        });
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*rns,signed_roots,signed_plan,signed_coefficients,rns_product,conjugate,
        [&](auto&output,const auto&input){
            for(std::size_t c=0;c<2;++c)Automorphism<P>(output[c],input[c],2*P::n-1);
            ckks_detail::reduceTRLWEToLevel<P,Basis::root_log_q>(output);
        });
    if(rns->ct!=integer->ct)throw std::runtime_error("signed RNS basis tensors differ from FFT tensors");
    if(root_conjugations!=3)throw std::runtime_error("conjugated roots were not reused");
    ckksSlotDecode<P,Basis::output_log_q,Basis::log_delta>(*slots,integer->ct[1]);
    double signed_error=0;
    for(std::size_t i=0;i<P::n/2;++i){
        double expected=signed_coefficients[0].real();
        for(std::size_t k=1;k<=40;++k)expected+=2*(signed_coefficients[k]*std::polar(1.0,k*(0.3+0.7*i/(P::n/2)))).real();
        signed_error=std::max(signed_error,std::abs((*slots)[i]-expected));
    }
    if(signed_error>=std::ldexp(1.0,-23))throw std::runtime_error("signed basis independent series gate failed");
    // A six-bit smaller product drop retains six extra scale bits. Exercise
    // direct roots, sums, differences and the nonzero constant independently.
    using Mixed = CKKSNarrowFrequencyBasisSchedule<Narrow,32,29>;
    static_assert(Mixed::product_log_delta == Basis::root_log_delta + 6);
    auto mixed = std::make_unique<Mixed::Output>();
    const auto mixed_conjugate = [&](auto& output,const auto& input) {
        for(std::size_t c=0;c<2;++c)Automorphism<P>(output[c],input[c],2*P::n-1);
        ckks_detail::reduceTRLWEToLevel<P,Mixed::output_log_q>(output);
    };
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Mixed>(*mixed,signed_roots,signed_plan,signed_coefficients,auxiliary,mixed_conjugate,
        [&](auto& output,const auto& input) {
            for(std::size_t c=0;c<2;++c)Automorphism<P>(output[c],input[c],2*P::n-1);
            ckks_detail::reduceTRLWEToLevel<P,Mixed::root_log_q>(output);
        });
    auto mixed_raw=std::make_unique<CKKSNarrowFrequencyBasisKey<Mixed>>();
    auto mixed_legacy=std::make_unique<Mixed::Output>();
    auto mixed_prepared=std::make_unique<Mixed::Output>();
    CKKSNarrowEvalSeriesFromBasisRoots<Mixed>(*mixed_legacy,roots,*mixed_raw,plan,coefficients);
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Mixed>(*mixed_prepared,roots,plan,coefficients,auxiliary,mixed_conjugate);
    auto mixed_rns=std::make_unique<Mixed::Output>();
    CKKSNarrowEvalSeriesFromBasisRootsPrepared<Mixed>(*mixed_rns,roots,plan,coefficients,rns_product,mixed_conjugate);
    if(mixed_rns->ct!=mixed_prepared->ct)throw std::runtime_error("mixed-scale RNS basis tensor differs from FFT tensor");
    if(mixed_legacy->ct!=mixed_prepared->ct)throw std::runtime_error("mixed-scale legacy and prepared backends differ");
    ckksSlotDecode<P,Mixed::output_log_q,Mixed::log_delta>(*slots,mixed->ct[1]);
    double mixed_error=0;
    for(std::size_t i=0;i<P::n/2;++i){
        double expected=signed_coefficients[0].real();
        for(std::size_t k=1;k<=40;++k)expected+=2*(signed_coefficients[k]*std::polar(1.0,k*(0.3+0.7*i/(P::n/2)))).real();
        mixed_error=std::max(mixed_error,std::abs((*slots)[i]-expected));
    }
    if(mixed_error>=std::ldexp(1.0,-23))throw std::runtime_error("mixed-scale independent series gate failed");
    std::cout<<"mixed-scale signed basis series error="<<mixed_error<<'\n';
    std::cout<<"seven roots cover degree40 with one product level; root conjugations="<<root_conjugations<<" series error="<<signed_error<<'\n';
    std::cout << "prepared FFT and auxiliary NTT basis byte-identical; series error=" << maximum << '\n';
}
