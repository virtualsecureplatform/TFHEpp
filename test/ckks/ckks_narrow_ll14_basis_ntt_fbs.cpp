// Experimental direct-selector real-input frequency-basis replay using retained LL14
// component keys. Never reconstruct selectors or generate private/evaluation
// keys here. The supplied public key artifacts must match the checkpoint.
#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_basis_prepared.hpp>
#include <ckks/prepared_compact_linear.hpp>
#include <ckks/narrow_dd_ntt.hpp>
#include <ckks/narrow_auto_ntt.hpp>
#include <ckks/key_projection.hpp>
#include <ckks/narrow_selector_projection.hpp>
#include <ckks/narrow_factor_ntt32.hpp>
#include <ckks/narrow_factor_rns32.hpp>
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
struct TableSource : CKKSNarrowBootstrapSchedule<P,48,440,57,42,31,false,0,true> {
 static constexpr uint32_t plain_log_delta=42,factor_log_q=398,factor_key_log_q=398;
 static constexpr uint32_t selector_matrix_digit_bits=54;
 using Factor=CKKSCiphertext<P,398,57>;
};
struct Narrow : CKKSNarrowBootstrapSchedule<P,48,427,44,42,31,true,13,true> {
 static constexpr uint32_t selector_matrix_digit_bits=54;
};
using Refresh=CKKSNarrowRefreshSchedule<Narrow,44,37,0,true,8>;
using Basis=CKKSNarrowFrequencyBasisSchedule<Narrow,28>;
using Output=CKKSCiphertext<P,93,44>;
using Square=CKKSMultResult<P,93,44,93,44>;
struct SourceLinear {
 SourceR::S2CPlan plan;
 CKKSGaloisKey<P,85> input;
 CKKSGaloisKey<P,48> output;
 CKKSSecretKeySwitchKey<P,48> sparse;
};
int main(int argc,char**argv)try {
 if(argc<8 || (argc-8)%2)throw std::runtime_error("usage: checkpoint-directory switch56-file sigmoid|gelu|tanh|all output-directory direct440-selectors saved-dd433-primary20 saved-packed435 [--cache-gib 90] [--s2c-cache-gib 9] [--direct-rns 1]");
 uint64_t cache_gib=90,s2c_cache_gib=9;bool direct_rns=true;
 for(int i=8;i<argc;i+=2){
  const std::string option=argv[i];const auto value=std::stoull(argv[i+1]);
  if(value>1024)throw std::invalid_argument("cache budget exceeds supported range");
  if(option=="--cache-gib")cache_gib=value;
  else if(option=="--s2c-cache-gib")s2c_cache_gib=value;
  else if(option=="--direct-rns"){if(value>1)throw std::invalid_argument("direct RNS flag must be 0 or 1");direct_rns=value;}
  else throw std::invalid_argument("unknown option");
 }
 const bool parallel_pairs=true,small_prime_factors=true;
 const auto checkpoint=std::filesystem::path(argv[1]),root=std::filesystem::path(argv[4]);const std::string requested=argv[3];
 if(requested!="all" && requested!="sigmoid" && requested!="gelu" && requested!="tanh")throw std::invalid_argument("unsupported function");
 if(!std::filesystem::is_directory(root))throw std::runtime_error("output directory missing");
 const std::string origin="f70c9f110231eac0b4fffe300cccd9ba89a519809e1cc89c45e8814f8f7e53af";
 const auto start=std::chrono::steady_clock::now();
 auto dense=std::make_unique<Key<P>>(),sparse=std::make_unique<Key<P>>();auto savedslots=std::make_unique<CKKSSlotVector<P>>();auto savedinput=std::make_unique<SourceR::Input>();
 auto old=std::make_unique<SourceLinear>();auto key=std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
 {std::ifstream f(checkpoint/"checkpoint.bin",std::ios::binary);if(!f)throw std::runtime_error("checkpoint missing");cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,9> meta;std::string savedfunction;ar(meta,savedfunction);
  if(meta!=std::array<uint64_t,9>{0x4e415252504b3031ULL,P::n,85,41,379,49,34,100,448} || savedfunction!="tanh")throw std::runtime_error("checkpoint profile mismatch");
  ar(*dense,*sparse,*savedslots,*savedinput,old->plan.has_direct,old->plan.has_conjugate,old->plan.direct,old->plan.conjugate,old->input,old->output,old->sparse);
  uint64_t head,tail;ar(head,tail);if(head || tail)throw std::runtime_error("expected file-backed selectors");auto original=std::make_unique<CKKSGaloisKey<P,345>>();ar(*original);CKKSReduceGaloisKey<P,93>(key->s2c_input_galois,*original);CKKSReduceGaloisKey<P,56>(key->s2c_output_galois,*original);CKKSReduceAutoKey<P,178>(key->narrow.output_conjugation,(*original)[P::nbit]);}
 CKKSRebindCompactLinearInputScale(key->s2c_plan,old->plan);
 old.reset();savedslots.reset();savedinput.reset();
 {std::ifstream f(argv[2],std::ios::binary);if(!f)throw std::runtime_error("saved switch key missing");cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,4> meta;ar(meta,key->to_sparse);if(meta!=std::array<uint64_t,4>{0x5332435357303031ULL,P::n,448,56} || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("switch key metadata/length mismatch");}
 using PreparedS2C=CKKSPreparedCompactRealLinearTransformFFT<P,93,52,37>;
 std::unique_ptr<PreparedS2C> prepared_s2c;
 std::unique_ptr<CKKSRetainedPreparedGaloisKeyFFT<P,93>> s2c_input_cache;
 std::unique_ptr<CKKSRetainedPreparedGaloisKeyFFT<P,56>> s2c_output_cache;
 auto s2c_setup_start=std::chrono::steady_clock::now();
 if(s2c_cache_gib){
  constexpr uint64_t galois_bytes=uint64_t{P::nbit}*sizeof(PolynomialInFD<P>)*
   (CKKSPreparedAutoKeyFFT<P,93>::spectrum_count+CKKSPreparedAutoKeyFFT<P,56>::spectrum_count);
  const uint64_t budget=s2c_cache_gib<<30;
  if(budget<galois_bytes)throw std::invalid_argument("S2C spectrum budget too small");
  prepared_s2c=std::make_unique<PreparedS2C>(key->s2c_plan,budget-galois_bytes);
  s2c_input_cache=std::make_unique<CKKSRetainedPreparedGaloisKeyFFT<P,93>>(key->s2c_input_galois);
  s2c_output_cache=std::make_unique<CKKSRetainedPreparedGaloisKeyFFT<P,56>>(key->s2c_output_galois);
 }
 const auto s2c_setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-s2c_setup_start).count();
 const uint64_t s2c_retained_bytes=prepared_s2c?prepared_s2c->spectrum_bytes+s2c_input_cache->spectrum_bytes()+s2c_output_cache->spectrum_bytes():0;
 std::cout<<"public S2C setup seconds="<<s2c_setup_seconds<<" spectrum bytes="<<s2c_retained_bytes<<std::endl;
 using File=CKKSNarrowSelectorFile<TableSource>;auto raw=std::make_shared<File>(argv[5],File::Existing{});
 uint64_t selector_retained_bytes=0;
 if(direct_rns){
  using Cache=CKKSNarrowFactorSelectorCacheRNS32<TableSource>;
  std::shared_ptr<const Cache> retained;
  if(cache_gib){retained=std::make_shared<Cache>(raw,cache_gib<<30,40,
   [&](bool tail,size_t group){std::cout<<"prepared direct RNS cache "<<(tail?"tail":"head")<<" group="<<group<<std::endl;});selector_retained_bytes=Cache::retained_bytes;}
  key->narrow.factor_source=std::make_shared<CKKSNarrowDeferredRescaleFactorSourceRNS32<Narrow,TableSource,TableSource>>(raw,retained);
 }
 else{
  using Cache=CKKSNarrowFactorSelectorCache32<TableSource>;
  std::shared_ptr<const Cache> retained;
  if(cache_gib){retained=std::make_shared<Cache>(raw,cache_gib<<30,40,
   [&](bool tail,size_t group){std::cout<<"prepared selector cache "<<(tail?"tail":"head")<<" group="<<group<<std::endl;});selector_retained_bytes=Cache::retained_bytes;}
  key->narrow.factor_source=std::make_shared<CKKSNarrowDeferredRescaleFactorSource32<Narrow,TableSource,TableSource>>(raw,retained);
 }
 CKKSNarrowUseStagedFFT(key->narrow,40);
 auto master=std::make_unique<CKKSCompactDDRelinKey<P,433,20,20>>();
 if(std::filesystem::file_size(argv[6])!=1+8+64+56+sizeof(master->data))throw std::runtime_error("DD433 length mismatch");
 {std::ifstream f(argv[6],std::ios::binary);cereal::PortableBinaryInputArchive ar(f);uint64_t length;ar(length);if(length!=64)throw std::runtime_error("DD433 origin length mismatch");std::string stored(64,' ');ar(cereal::binary_data(stored.data(),stored.size()),*master);if(stored!=origin||f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("DD433 origin mismatch");}
 CKKSNarrowUseAuxiliaryNTTDD<48,20,true>(key->narrow,*master);
 auto autoStart=std::chrono::steady_clock::now();
 auto packed=CKKSNarrowLoadPackedKey<PackedSource>(argv[7]);
 CKKSNarrowUsePublicNTTAutomorphisms<398,20,20>(key->narrow,*packed);
 PublicPackedAutoDeriver<P,PackedSource,20> derive(*packed);auto [rows,steps]=derive.derive(2*P::n-1);
 auto basis_auto=std::make_unique<PreparedPublicAutoNTT<P,435,20,20>>(*rows);packed.reset();
 auto basis_raw=makeCKKSCompactDDRelinKeyAtLevel<213>(*master);auto basis_ntt=makeCKKSDDRelinKeyNTT<20>(*basis_raw);
 CKKSAuxiliaryDDRelinKeyNTT<P,165,48,20,20> basis_product(basis_ntt);
 const auto auto_setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-autoStart).count();

 auto input=std::make_unique<Refresh::Input>();auto slots=std::make_unique<CKKSSlotVector<P>>();const auto inputfile=root/"candidate-input.bin";
 const std::array<uint64_t,5> inputmeta{0x4c4c3134494e3031ULL,P::n,448,85,44};bool reused=std::filesystem::exists(inputfile);
 if(reused){std::ifstream f(inputfile,std::ios::binary);cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,5> meta;std::string stored;ar(meta,stored,*slots,*input);if(meta!=inputmeta || stored!=origin || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("retained input mismatch");}
 else throw std::runtime_error("retained candidate-input.bin required; no input encryption in this replay");
 const auto setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 for(const auto& function:requested=="all"?std::vector<std::string>{"gelu","tanh","sigmoid"}:std::vector<std::string>{requested}) {
 auto target=[&](double x){if(function=="sigmoid")return 1/(1+std::exp(-x));if(function=="gelu")return 0.5*x*(1+std::erf(x/std::sqrt(2.0)));return std::tanh(x);};
 auto curvature=[&](double x){double y=target(x);if(function=="sigmoid")return y*(1-y)*(1-2*y);if(function=="gelu")return (2-x*x)*std::exp(-x*x/2)/std::sqrt(2*std::acos(-1.0));return -2*y*(1-y*y);};
 const uint32_t degree=function=="tanh"?40:12;
 auto approx=function=="tanh"?CKKSBuildNarrowWindowedTrigonometricApproximation(target,[](double x){double y=std::tanh(x);return 1-y*y;},curvature,{1,1,2,4,40},4,16,degree):CKKSBuildNarrowTrigonometricApproximation(target,4,16,degree,function=="gelu"?2.0:1.0,65536,40,2.0,curvature,160.0);
 size_t groups=0;key->narrow.factor_observer=[&](bool tail,size_t group,const auto& timing){++groups;std::cout<<function<<" completed "<<(tail?"tail":"head")<<" group="<<group<<" factor_groups="<<groups
  <<" load="<<timing.selector_load_seconds<<" table_ntt="<<timing.table_ntt_seconds
  <<" selector_ntt="<<timing.selector_ntt_seconds<<" matrix="<<timing.matrix_seconds
  <<" inverse_crt="<<timing.inverse_crt_seconds<<" conjugation="<<timing.conjugation_seconds<<std::endl;};
 auto prepared=std::make_unique<Refresh::Coefficients>();auto t=std::chrono::steady_clock::now();if(prepared_s2c)
 CKKSNarrowPrepareCoefficients<Refresh>(*prepared,*input,*key,[&](auto &result,const auto &ct){
  CKKSRealLinearTransform<P,93,52,37>(result,ct,*prepared_s2c,*s2c_input_cache,*s2c_output_cache);
 });
 else CKKSNarrowPrepareCoefficients<Refresh>(*prepared,*input,*key);
 const auto s2c_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();std::cout<<function<<" S2C completed seconds="<<s2c_seconds<<std::endl;
 key->narrow.root_stage_observer=[&](const auto&t){std::cout<<function<<" stages factor="<<t.leaf_factors_and_inline_alignment_seconds<<" align="<<t.bounded_tail_alignment_seconds<<" dd_cache="<<t.dd_cache_preparation_seconds<<" trees="<<t.dd_trees_and_beta_seconds<<std::endl;};
 t=std::chrono::steady_clock::now();CKKSNarrowFrequencyBasis plan(degree,degree==12?std::vector<uint32_t>{1,3,5,6}:std::vector<uint32_t>{1,3,4,9,11,16,17,19,20});
 std::vector<CKKSNarrowBlindRotationRequest<Narrow,44>> requests;for(auto freq:plan.frequencies())requests.push_back({prepared.get(),freq,{1,0}});
 std::vector<Narrow::Output> roots;CKKSNarrowBlindRotateBatch<Narrow,44>(roots,requests,key->narrow,plan.frequencies().size());
 auto weighted=std::make_unique<Basis::Output>();CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*weighted,roots,plan,approx.coefficients,basis_product,
 [&](auto& result,const auto& input){basis_auto->template apply<121>(result,input);});const auto blind_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
 auto out=std::make_unique<Output>();CKKSScaleDown<P,121,72,28>(*out,*weighted);auto square=std::make_unique<Square>();auto postkey=makeCKKSCompactDDRelinKeyAtLevel<141>(*master);auto postntt=makeCKKSDDRelinKeyNTT<20>(*postkey);CKKSAuxiliaryDDRelinKeyNTT<P,93,48,20,20> auxiliary(postntt);CKKSMultRelinThenRescale<P>(*square,*out,*out,auxiliary);
 auto decoded=std::make_unique<CKKSSlotVector<P>>();ckksSlotDecrypt<P>(*decoded,*out,*dense);double error=0,amplitude=0;
 for(size_t i=0;i<P::n/2;++i){const double v=target((*slots)[i].real());error=narrow_benchmark::finiteMaximum(error,std::abs((*decoded)[i]-std::complex<double>{v,0}));amplitude=narrow_benchmark::finiteMaximum(amplitude,std::abs(v));}
 ckksSlotDecrypt<P>(*decoded,*square,*dense);double product_error=0;for(size_t i=0;i<P::n/2;++i){double v=target((*slots)[i].real());product_error=narrow_benchmark::finiteMaximum(product_error,std::abs((*decoded)[i]-std::complex<double>{v*v,0}));}
 const auto tolerance=std::ldexp(1.0,-23),product_tolerance=(2*amplitude+1)*tolerance+tolerance*tolerance;bool passed=error<tolerance && product_error<product_tolerance;
 {std::ofstream f(root/(function+"-outputs.bin"),std::ios::binary);cereal::PortableBinaryOutputArchive ar(f);const std::array<uint64_t,9> meta{0x4c4c313446425331ULL,P::n,448,85,44,48,93,44,49};ar(meta,function,*slots,*input,*prepared,*out,*square);}
 {std::ofstream f(root/(function+"-report.json"));f<<std::setprecision(17)<<"{\"function\":\""<<function<<"\",\"n\":"<<P::n<<",\"complex_input\":false,\"ntt_dd\":true,\"parallel_product_pairs\":"<<(parallel_pairs?"true":"false")<<",\"small_prime_factors\":"<<(small_prime_factors?"true":"false")<<",\"dd_primary_bits\":20,\"dd_secondary_bits\":20,\"dd_auxiliary_bits\":48,\"public_derived_ntt_automorphisms\":true,\"conjugation_primary_bits\":20,\"automorphism_setup_seconds\":"<<auto_setup_seconds<<",\"direct_rns_factors\":"<<(direct_rns?"true":"false")<<",\"selector_residue_rows\":"<<(direct_rns?17:36)<<",\"frequency_basis\":true,\"selector_cache_retained_bytes\":"<<selector_retained_bytes<<",\"degree\":"<<degree<<",\"selector_matrix_digit_bits\":"<<(direct_rns?0:52)<<",\"selector_chunk_size\":40,\"fresh_private_keys_generated\":false,\"new_evaluation_keys_generated\":0,\"input_reused\":"<<(reused?"true":"false")<<",\"key_load_and_input_setup_seconds\":"<<setup_seconds<<",\"s2c_public_setup_seconds\":"<<s2c_setup_seconds<<",\"s2c_retained_spectrum_bytes\":"<<s2c_retained_bytes<<",\"s2c_public_keys_prepared\":"<<(prepared_s2c?s2c_input_cache->prepared_keys+s2c_output_cache->prepared_keys:0)<<",\"s2c_seconds\":"<<s2c_seconds<<",\"blind_rotation_and_series_seconds\":"<<blind_seconds<<",\"functional_evaluation_seconds\":"<<s2c_seconds+blind_seconds<<",\"output_log_q\":93,\"output_log_delta\":44,\"post_product_log_q\":49,\"max_error\":"<<error<<",\"post_product_max_error\":"<<product_error<<",\"product_tolerance\":"<<product_tolerance<<",\"accuracy_passed\":"<<(passed?"true":"false")<<",\"security_validated\":false}";}
 std::cout<<function<<" error="<<error<<" square error="<<product_error<<" seconds="<<s2c_seconds+blind_seconds<<std::endl;if(!passed)throw std::runtime_error("functional accuracy gate failed");
}
} catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
