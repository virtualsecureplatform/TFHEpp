// Experimental packed-selector real-input functional replay using retained LL14
// component keys. Never reconstruct selectors or generate private/evaluation
// keys here. The supplied public key artifacts must match the checkpoint.
#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_dd_ntt.hpp>
#include <ckks/narrow_auto_ntt.hpp>
#include <ckks/key_projection.hpp>
#include <ckks/narrow_selector_projection.hpp>
#include <ckks/narrow_factor_ntt32.hpp>
#include <iomanip>
using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<14,7>;
struct PackedP:P { static constexpr uint32_t B̅gbit=8,l̅=56; };
using PackedSource=CKKSNarrowBootstrapSchedule<PackedP,48,387,57,34,31>;
struct Source : CKKSNarrowBootstrapSchedule<P,48,379,49,34,31,false,0,true> {
 static constexpr uint32_t selector_matrix_digit_bits=60;
};
using SourceR=CKKSNarrowRefreshSchedule<Source,41,37,0,true>;
using PackedSchedule=CKKSNarrowBootstrapSchedule<P,48,387,57,34,31>;
using SelectorArchive=CKKSNarrowPackedReconstructionSchedule<PackedSchedule,42>;
// High-scale factor construction precedes rescaling to the tree schedule.
struct TableSource : CKKSNarrowBootstrapSchedule<P,48,379,57,34,31,false,0,true> {
 static constexpr uint32_t plain_log_delta=42,factor_log_q=337,factor_key_log_q=337;
 static constexpr uint32_t selector_matrix_digit_bits=54;
 using Factor=CKKSCiphertext<P,337,57>;
};
struct Narrow : CKKSNarrowBootstrapSchedule<P,48,366,44,42,31,true,13,true> {
 static constexpr uint32_t selector_matrix_digit_bits=54;
};
using Refresh=CKKSNarrowRefreshSchedule<Narrow,44,37,0,true,8>;
using Output=CKKSCiphertext<P,101,44>;
using Square=CKKSMultResult<P,101,44,101,44>;
struct SourceLinear {
 SourceR::S2CPlan plan;
 CKKSGaloisKey<P,85> input;
 CKKSGaloisKey<P,48> output;
 CKKSSecretKeySwitchKey<P,48> sparse;
};
int main(int argc,char**argv)try {
 if(argc<9 || argc>11)throw std::runtime_error("usage: checkpoint-directory dd-chain-file switch56-file sigmoid|gelu|tanh output-directory reconstructed-selector-directory saved-dd433 saved-packed435 [--parallel-pairs] [--small-prime-factors]");
 bool parallel_pairs=false,small_prime_factors=false;
 for(int i=9;i<argc;++i) {
  const std::string option=argv[i];
  if(option=="--parallel-pairs" && !parallel_pairs)parallel_pairs=true;
  else if(option=="--small-prime-factors" && !small_prime_factors)small_prime_factors=true;
  else throw std::invalid_argument("unknown or duplicate benchmark option");
 }
 const auto checkpoint=std::filesystem::path(argv[1]),root=std::filesystem::path(argv[5]);const std::string function=argv[4];
 if(function!="sigmoid" && function!="gelu" && function!="tanh")throw std::invalid_argument("unsupported function");
 if(!std::filesystem::is_directory(root))throw std::runtime_error("output directory missing");
 const std::string origin="f70c9f110231eac0b4fffe300cccd9ba89a519809e1cc89c45e8814f8f7e53af";
 const auto start=std::chrono::steady_clock::now();
 auto dense=std::make_unique<Key<P>>(),sparse=std::make_unique<Key<P>>();auto savedslots=std::make_unique<CKKSSlotVector<P>>();auto savedinput=std::make_unique<SourceR::Input>();
 auto old=std::make_unique<SourceLinear>();auto key=std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
 {std::ifstream f(checkpoint/"checkpoint.bin",std::ios::binary);if(!f)throw std::runtime_error("checkpoint missing");cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,9> meta;std::string savedfunction;ar(meta,savedfunction);
  if(meta!=std::array<uint64_t,9>{0x4e415252504b3031ULL,P::n,85,41,379,49,34,100,448} || savedfunction!="tanh")throw std::runtime_error("checkpoint profile mismatch");
  ar(*dense,*sparse,*savedslots,*savedinput,old->plan.has_direct,old->plan.has_conjugate,old->plan.direct,old->plan.conjugate,old->input,old->output,old->sparse);
  uint64_t head,tail;ar(head,tail);if(head || tail)throw std::runtime_error("expected file-backed selectors");auto original=std::make_unique<CKKSGaloisKey<P,345>>();ar(*original);CKKSReduceGaloisKey<P,337>(key->narrow.factor_galois,*original);}
 CKKSRebindCompactLinearInputScale(key->s2c_plan,old->plan);
 CKKSReduceGaloisKey<P,93>(key->s2c_input_galois,key->narrow.factor_galois);
 CKKSReduceGaloisKey<P,56>(key->s2c_output_galois,key->narrow.factor_galois);
 CKKSReduceAutoKey<P,117>(key->narrow.output_conjugation,key->narrow.factor_galois[P::nbit]);
 old.reset();savedslots.reset();savedinput.reset();
 {std::ifstream f(argv[3],std::ios::binary);if(!f)throw std::runtime_error("saved switch key missing");cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,4> meta;ar(meta,key->to_sparse);if(meta!=std::array<uint64_t,4>{0x5332435357303031ULL,P::n,448,56} || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("switch key metadata/length mismatch");}
 auto post=std::make_unique<CKKSCompactDDRelinKey<P,101>>();
 {std::ifstream f(argv[2],std::ios::binary);if(!f)throw std::runtime_error("saved DD chain missing");cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,6> meta;std::string stored;ar(meta,stored,key->narrow.relin,*post);if(meta!=std::array<uint64_t,6>{0x4444434841494e31ULL,P::n,448,368,44,5} || stored!=origin || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("DD chain origin/metadata/length mismatch");}
 using File=CKKSNarrowSelectorFile<SelectorArchive>;auto raw=std::make_shared<File>(argv[6],File::Existing{});
 std::shared_ptr<CKKSNarrowFactorSource<Narrow>> source;
 if(small_prime_factors)source=std::make_shared<CKKSNarrowDeferredRescaleFactorSource32<Narrow,TableSource,SelectorArchive>>(raw);
 else source=std::make_shared<CKKSNarrowDeferredRescaleFactorSource<Narrow,TableSource,SelectorArchive>>(raw);
 key->narrow.factor_source=std::make_shared<CKKSNarrowPublicMaskNormalizedFactorSource<Narrow,SelectorArchive>>(source);
 CKKSNarrowUseStagedFFT(key->narrow,40);
 auto master=std::make_unique<CKKSCompactDDRelinKey<P,433,64,16>>();
 if(std::filesystem::file_size(argv[7])!=1+8+64+56+sizeof(master->data))throw std::runtime_error("DD433 length mismatch");
 {std::ifstream f(argv[7],std::ios::binary);cereal::PortableBinaryInputArchive ar(f);uint64_t length;ar(length);if(length!=64)throw std::runtime_error("DD433 origin length mismatch");std::string stored(64,' ');ar(cereal::binary_data(stored.data(),stored.size()),*master);if(stored!=origin||f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("DD433 origin mismatch");}
 if(parallel_pairs)CKKSNarrowUseAuxiliaryNTTDD<100,32,true>(key->narrow,*master);
 else CKKSNarrowUseAuxiliaryNTTDD(key->narrow,*master);
 auto autoStart=std::chrono::steady_clock::now();
 auto packed=CKKSNarrowLoadPackedKey<PackedSource>(argv[8]);
 CKKSNarrowUsePublicNTTAutomorphisms<337>(key->narrow,*packed);packed.reset();
 const auto auto_setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-autoStart).count();

 auto input=std::make_unique<Refresh::Input>();auto slots=std::make_unique<CKKSSlotVector<P>>();const auto inputfile=root/"candidate-input.bin";
 const std::array<uint64_t,5> inputmeta{0x4c4c3134494e3031ULL,P::n,448,85,44};bool reused=std::filesystem::exists(inputfile);
 if(reused){std::ifstream f(inputfile,std::ios::binary);cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,5> meta;std::string stored;ar(meta,stored,*slots,*input);if(meta!=inputmeta || stored!=origin || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("retained input mismatch");}
 else throw std::runtime_error("retained candidate-input.bin required; no input encryption in this replay");
 const auto setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 auto target=[&](double x){if(function=="sigmoid")return 1/(1+std::exp(-x));if(function=="gelu")return 0.5*x*(1+std::erf(x/std::sqrt(2.0)));return std::tanh(x);};
 auto curvature=[&](double x){double y=target(x);if(function=="sigmoid")return y*(1-y)*(1-2*y);if(function=="gelu")return (2-x*x)*std::exp(-x*x/2)/std::sqrt(2*std::acos(-1.0));return -2*y*(1-y*y);};
 const uint32_t degree=function=="tanh"?40:12;
 auto approx=function=="tanh"?CKKSBuildNarrowWindowedTrigonometricApproximation(target,[](double x){double y=std::tanh(x);return 1-y*y;},curvature,{1,1,2,4,40},4,16,degree):CKKSBuildNarrowTrigonometricApproximation(target,4,16,degree,function=="gelu"?2.0:1.0,65536,40,2.0,curvature,160.0);
 size_t groups=0;key->narrow.factor_observer=[&](bool tail,size_t group,const auto& timing){++groups;std::cout<<function<<" completed "<<(tail?"tail":"head")<<" group="<<group<<" factor_groups="<<groups
  <<" load="<<timing.selector_load_seconds<<" table_ntt="<<timing.table_ntt_seconds
  <<" selector_ntt="<<timing.selector_ntt_seconds<<" matrix="<<timing.matrix_seconds
  <<" inverse_crt="<<timing.inverse_crt_seconds<<" conjugation="<<timing.conjugation_seconds<<std::endl;};
 auto prepared=std::make_unique<Refresh::Coefficients>();auto t=std::chrono::steady_clock::now();CKKSNarrowPrepareCoefficients<Refresh>(*prepared,*input,*key);const auto s2c_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();std::cout<<function<<" S2C completed seconds="<<s2c_seconds<<std::endl;
 key->narrow.root_stage_observer=[&](const auto&t){std::cout<<function<<" stages factor="<<t.leaf_factors_and_inline_alignment_seconds<<" align="<<t.bounded_tail_alignment_seconds<<" dd_cache="<<t.dd_cache_preparation_seconds<<" trees="<<t.dd_trees_and_beta_seconds<<std::endl;};
 auto originalout=std::make_unique<Refresh::Output>();t=std::chrono::steady_clock::now();CKKSNarrowFunctionalRefreshPrepared<Refresh>(*originalout,*prepared,*key,approx.coefficients,false,degree);const auto blind_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
 auto out=std::make_unique<Output>();CKKSLevelReduce<P>(*out,*originalout);auto square=std::make_unique<Square>();auto postkey=makeCKKSCompactDDRelinKeyAtLevel<201>(*master);auto postntt=makeCKKSDDRelinKeyNTT(*postkey);CKKSAuxiliaryDDRelinKeyNTT<P,101,100> auxiliary(postntt);CKKSMultRelinThenRescale<P>(*square,*out,*out,auxiliary);
 auto decoded=std::make_unique<CKKSSlotVector<P>>();ckksSlotDecrypt<P>(*decoded,*out,*dense);double error=0,amplitude=0;
 for(size_t i=0;i<P::n/2;++i){const double v=target((*slots)[i].real());error=narrow_benchmark::finiteMaximum(error,std::abs((*decoded)[i]-std::complex<double>{v,0}));amplitude=narrow_benchmark::finiteMaximum(amplitude,std::abs(v));}
 ckksSlotDecrypt<P>(*decoded,*square,*dense);double product_error=0;for(size_t i=0;i<P::n/2;++i){double v=target((*slots)[i].real());product_error=narrow_benchmark::finiteMaximum(product_error,std::abs((*decoded)[i]-std::complex<double>{v*v,0}));}
 const auto tolerance=std::ldexp(1.0,-23),product_tolerance=(2*amplitude+1)*tolerance+tolerance*tolerance;bool passed=error<tolerance && product_error<product_tolerance;
 {std::ofstream f(root/(function+"-outputs.bin"),std::ios::binary);cereal::PortableBinaryOutputArchive ar(f);const std::array<uint64_t,9> meta{0x4c4c313446425331ULL,P::n,448,85,44,48,101,44,57};ar(meta,function,*slots,*input,*prepared,*out,*square);}
 {std::ofstream f(root/(function+"-report.json"));f<<std::setprecision(17)<<"{\"function\":\""<<function<<"\",\"n\":"<<P::n<<",\"complex_input\":false,\"ntt_dd\":true,\"parallel_product_pairs\":"<<(parallel_pairs?"true":"false")<<",\"small_prime_factors\":"<<(small_prime_factors?"true":"false")<<",\"dd_primary_bits\":64,\"dd_secondary_bits\":32,\"dd_auxiliary_bits\":100,\"public_derived_ntt_automorphisms\":true,\"conjugation_primary_bits\":64,\"automorphism_setup_seconds\":"<<auto_setup_seconds<<",\"degree\":"<<degree<<",\"selector_matrix_digit_bits\":"<<(small_prime_factors?52:54)<<",\"selector_chunk_size\":40,\"fresh_private_keys_generated\":false,\"new_evaluation_keys_generated\":0,\"input_reused\":"<<(reused?"true":"false")<<",\"key_load_and_input_setup_seconds\":"<<setup_seconds<<",\"s2c_seconds\":"<<s2c_seconds<<",\"blind_rotation_and_series_seconds\":"<<blind_seconds<<",\"functional_evaluation_seconds\":"<<s2c_seconds+blind_seconds<<",\"output_log_q\":101,\"output_log_delta\":44,\"post_product_log_q\":57,\"max_error\":"<<error<<",\"post_product_max_error\":"<<product_error<<",\"product_tolerance\":"<<product_tolerance<<",\"accuracy_passed\":"<<(passed?"true":"false")<<",\"security_validated\":false}";}
 std::cout<<function<<" error="<<error<<" square error="<<product_error<<" seconds="<<s2c_seconds+blind_seconds<<std::endl;if(!passed)throw std::runtime_error("functional accuracy gate failed");
} catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
