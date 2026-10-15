// Experimental direct-selector real-input frequency-basis replay using retained LL14
// component keys. Mixed tree43/product37 profile; requires fresh qualification.
// Never reconstruct selectors or generate private/evaluation
// keys here. The supplied public key artifacts must match the checkpoint.
#include "ckks_narrow_basis_benchmark.hpp"
#include "ckks_benchmark_affinity.hpp"
#include <ckks/narrow_basis_prepared.hpp>
#include <ckks/prepared_compact_linear_rns32.hpp>
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
// Experimental mixed-scale profile: selectors are publicly reduced, not rescaled.
// Its lower data levels and wider gadgets require fresh accuracy qualification.
struct RawTableSource : CKKSNarrowBootstrapSchedule<P,48,440,57,42,31,false,0,true> {
 static constexpr uint32_t plain_log_delta=42,factor_log_q=398,factor_key_log_q=398;
 static constexpr uint32_t selector_matrix_digit_bits=54;
 using Factor=CKKSCiphertext<P,398,57>;
};
struct TableSource : CKKSNarrowBootstrapSchedule<P,48,433,57,42,31,false,0,true> {
 static constexpr uint32_t plain_log_delta=42,factor_log_q=391,factor_key_log_q=391;
 static constexpr uint32_t selector_matrix_digit_bits=52,selector_rns_prime_count=17;
 using Factor=CKKSCiphertext<P,391,57>;
};
struct Narrow : CKKSNarrowBootstrapSchedule<P,48,419,43,42,31,true,14,true> {
 static constexpr uint32_t selector_matrix_digit_bits=54;
};
using Refresh=CKKSNarrowRefreshSchedule<Narrow,44,37,0,true,8>;
using Basis=CKKSNarrowFrequencyBasisSchedule<Narrow,28,37>;
using Output=CKKSCiphertext<P,92,44>;
using Square=CKKSMultResult<P,92,44,92,44>;
struct SourceLinear {
 SourceR::S2CPlan plan;
 CKKSGaloisKey<P,85> input;
 CKKSGaloisKey<P,48> output;
 CKKSSecretKeySwitchKey<P,48> sparse;
};
int main(int argc,char**argv)try {
 if(argc<8 || (argc-8)%2)throw std::runtime_error("usage: checkpoint-directory switch56-file sigmoid|gelu|tanh|all output-directory direct440-selectors saved-dd433-primary20 saved-packed435 [--cache-gib 90] [--s2c-cache-gib 9] [--direct-rns 1] [--s2c-rns 1] [--tensor-rns 1] [--auto-rns 1] [--signed-basis 1] [--basis-pair 0] [--dd-wide-rns 1] [--dd-pair 0] [--thread-sweep 0] [--matrix-tiled 1] [--matrix-pair 0]");
 uint64_t cache_gib=90,s2c_cache_gib=9;bool direct_rns=true,s2c_rns=true,tensor_rns=true,auto_rns=true,signed_basis=true,basis_pair=false,dd_wide_rns=true,dd_pair=false,thread_sweep=false,matrix_tiled=true,matrix_pair=false;
 for(int i=8;i<argc;i+=2){
  const std::string option=argv[i];const auto value=std::stoull(argv[i+1]);
  if(value>1024)throw std::invalid_argument("cache budget exceeds supported range");
  if(option=="--cache-gib")cache_gib=value;
  else if(option=="--s2c-cache-gib")s2c_cache_gib=value;
  else if(option=="--matrix-tiled"){if(value>1)throw std::invalid_argument("matrix tiled flag must be 0 or 1");matrix_tiled=value;}
  else if(option=="--matrix-pair"){if(value>1)throw std::invalid_argument("matrix pair flag must be 0 or 1");matrix_pair=value;}
  else if(option=="--thread-sweep"){if(value>1)throw std::invalid_argument("thread sweep flag must be 0 or 1");thread_sweep=value;}
  else if(option=="--dd-wide-rns"){if(value>1)throw std::invalid_argument("DD wide RNS flag must be 0 or 1");dd_wide_rns=value;}
  else if(option=="--dd-pair"){if(value>1)throw std::invalid_argument("DD pair flag must be 0 or 1");dd_pair=value;}
  else if(option=="--basis-pair"){if(value>1)throw std::invalid_argument("basis pair flag must be 0 or 1");basis_pair=value;}
  else if(option=="--signed-basis"){if(value>1)throw std::invalid_argument("signed basis flag must be 0 or 1");signed_basis=value;}
  else if(option=="--auto-rns"){if(value>1)throw std::invalid_argument("automorphism RNS flag must be 0 or 1");auto_rns=value;}
  else if(option=="--tensor-rns"){if(value>1)throw std::invalid_argument("tensor RNS flag must be 0 or 1");tensor_rns=value;}
  else if(option=="--s2c-rns"){if(value>1)throw std::invalid_argument("S2C RNS flag must be 0 or 1");s2c_rns=value;}
  else if(option=="--direct-rns"){if(value>1)throw std::invalid_argument("direct RNS flag must be 0 or 1");direct_rns=value;}
  else throw std::invalid_argument("unknown option");
 }
 const bool parallel_pairs=true,small_prime_factors=true;
 const auto checkpoint=std::filesystem::path(argv[1]),root=std::filesystem::path(argv[4]);const std::string requested=argv[3];
 if(requested!="all" && requested!="sigmoid" && requested!="gelu" && requested!="tanh")throw std::invalid_argument("unsupported function");
 if((basis_pair||dd_pair||thread_sweep||matrix_pair)&&requested!="tanh")throw std::invalid_argument("paired replay requires tanh");
 if(int(basis_pair)+int(dd_pair)+int(thread_sweep)+int(matrix_pair)>1)throw std::invalid_argument("paired replay controls must be separate");
 if(matrix_pair&&!direct_rns)throw std::invalid_argument("matrix pair requires direct RNS factors");
 const auto cpu_profiles=(thread_sweep||matrix_pair)?narrow_benchmark::coreSweepProfiles():std::vector<narrow_benchmark::CPUProfile>{};
 if(thread_sweep||matrix_pair)narrow_benchmark::applyCPUProfile(cpu_profiles.front());
 if(!std::filesystem::is_directory(root))throw std::runtime_error("output directory missing");
 const std::string origin="f70c9f110231eac0b4fffe300cccd9ba89a519809e1cc89c45e8814f8f7e53af";
 const auto start=std::chrono::steady_clock::now();
 auto dense=std::make_unique<Key<P>>(),sparse=std::make_unique<Key<P>>();auto savedslots=std::make_unique<CKKSSlotVector<P>>();auto savedinput=std::make_unique<SourceR::Input>();
 auto old=std::make_unique<SourceLinear>();auto key=std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
 {std::ifstream f(checkpoint/"checkpoint.bin",std::ios::binary);if(!f)throw std::runtime_error("checkpoint missing");cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,9> meta;std::string savedfunction;ar(meta,savedfunction);
  if(meta!=std::array<uint64_t,9>{0x4e415252504b3031ULL,P::n,85,41,379,49,34,100,448} || savedfunction!="tanh")throw std::runtime_error("checkpoint profile mismatch");
  ar(*dense,*sparse,*savedslots,*savedinput,old->plan.has_direct,old->plan.has_conjugate,old->plan.direct,old->plan.conjugate,old->input,old->output,old->sparse);
  uint64_t head,tail;ar(head,tail);if(head || tail)throw std::runtime_error("expected file-backed selectors");auto original=std::make_unique<CKKSGaloisKey<P,345>>();ar(*original);CKKSReduceGaloisKey<P,93>(key->s2c_input_galois,*original);CKKSReduceGaloisKey<P,56>(key->s2c_output_galois,*original);CKKSReduceAutoKey<P,Narrow::output_key_log_q>(key->narrow.output_conjugation,(*original)[P::nbit]);}
 CKKSRebindCompactLinearInputScale(key->s2c_plan,old->plan);
 old.reset();savedslots.reset();savedinput.reset();
 {std::ifstream f(argv[2],std::ios::binary);if(!f)throw std::runtime_error("saved switch key missing");cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,4> meta;ar(meta,key->to_sparse);if(meta!=std::array<uint64_t,4>{0x5332435357303031ULL,P::n,448,56} || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("switch key metadata/length mismatch");}
 using PreparedS2C=CKKSPreparedCompactRealLinearTransformFFT<P,93,52,37>;
 std::unique_ptr<PreparedS2C> prepared_s2c;
 using PreparedS2CRNS=CKKSPreparedCompactRealLinearTransformRNS32<P,93,52,37,5>;
 std::unique_ptr<PreparedS2CRNS> prepared_s2c_rns;
 std::unique_ptr<CKKSRetainedPreparedGaloisKeyFFT<P,93>> s2c_input_cache;
 std::unique_ptr<CKKSRetainedPreparedGaloisKeyFFT<P,56>> s2c_output_cache;
 auto s2c_setup_start=std::chrono::steady_clock::now();
 if(s2c_cache_gib){
  constexpr uint64_t galois_bytes=uint64_t{P::nbit}*sizeof(PolynomialInFD<P>)*
   (CKKSPreparedAutoKeyFFT<P,93>::spectrum_count+CKKSPreparedAutoKeyFFT<P,56>::spectrum_count);
  const uint64_t budget=s2c_cache_gib<<30;
  if(budget<galois_bytes)throw std::invalid_argument("S2C spectrum budget too small");
  if(s2c_rns)prepared_s2c_rns=std::make_unique<PreparedS2CRNS>(key->s2c_plan,budget-galois_bytes);
  else prepared_s2c=std::make_unique<PreparedS2C>(key->s2c_plan,budget-galois_bytes);
  s2c_input_cache=std::make_unique<CKKSRetainedPreparedGaloisKeyFFT<P,93>>(key->s2c_input_galois);
  s2c_output_cache=std::make_unique<CKKSRetainedPreparedGaloisKeyFFT<P,56>>(key->s2c_output_galois);
 }
 const bool have_prepared_s2c=bool(prepared_s2c)||bool(prepared_s2c_rns);
 if(have_prepared_s2c)key->s2c_plan={};
 const auto s2c_setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-s2c_setup_start).count();
 const uint64_t s2c_retained_bytes=have_prepared_s2c?(prepared_s2c_rns?prepared_s2c_rns->spectrum_bytes:prepared_s2c->spectrum_bytes)+s2c_input_cache->spectrum_bytes()+s2c_output_cache->spectrum_bytes():0;
 std::cout<<"public S2C setup seconds="<<s2c_setup_seconds<<" spectrum bytes="<<s2c_retained_bytes<<std::endl;
 using File=CKKSNarrowSelectorFile<RawTableSource>;auto archived=std::make_shared<File>(argv[5],File::Existing{});
 auto raw=std::make_shared<CKKSNarrowReducedSelectorSource<TableSource,RawTableSource>>(archived);
 uint64_t selector_retained_bytes=0;
 decltype(key->narrow.factor_source) coarse_matrix,tiled_matrix_source;
 if(direct_rns){
  using Cache=CKKSNarrowFactorSelectorCacheRNS32<TableSource>;
  std::shared_ptr<const Cache> retained;
  if(cache_gib){retained=std::make_shared<Cache>(raw,cache_gib<<30,40,
   [&](bool tail,size_t group){std::cout<<"prepared direct RNS cache "<<(tail?"tail":"head")<<" group="<<group<<std::endl;});selector_retained_bytes=Cache::retained_bytes;}
  using Source=CKKSNarrowDeferredRescaleFactorSourceRNS32<Narrow,TableSource,TableSource>;
  if(matrix_pair){coarse_matrix=std::make_shared<Source>(raw,retained,false);tiled_matrix_source=std::make_shared<Source>(raw,retained,true);key->narrow.factor_source=coarse_matrix;}
  else key->narrow.factor_source=std::make_shared<Source>(raw,retained,matrix_tiled);
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
 auto grouped_master=makeCKKSCompactDDRelinKeyWithPrimaryBase<40>(*master);
 decltype(key->narrow.product_tree_evaluator) digit_tree,wide_tree;
 const auto prepare_tree=[&](bool wide){
  if(wide){
   if(tensor_rns)CKKSNarrowUseAuxiliaryNTTDD<56,20,true,true,true>(key->narrow,*grouped_master);
   else CKKSNarrowUseAuxiliaryNTTDD<56,20,true,false,true>(key->narrow,*grouped_master);
  }else{
   if(tensor_rns)CKKSNarrowUseAuxiliaryNTTDD<56,20,true,true>(key->narrow,*grouped_master);
   else CKKSNarrowUseAuxiliaryNTTDD<56,20,true>(key->narrow,*grouped_master);
  }
 };
 if(dd_pair){prepare_tree(false);digit_tree=key->narrow.product_tree_evaluator;prepare_tree(true);wide_tree=key->narrow.product_tree_evaluator;}
 else prepare_tree(dd_wide_rns);
 auto autoStart=std::chrono::steady_clock::now();
 auto packed=CKKSNarrowLoadPackedKey<PackedSource>(argv[7]);
 if(auto_rns)CKKSNarrowUsePublicNTTAutomorphisms<391,20,20,true,125>(key->narrow,*packed);
 else CKKSNarrowUsePublicNTTAutomorphisms<391,20,20,false,125>(key->narrow,*packed);
 packed.reset();
 auto basis_raw=makeCKKSCompactDDRelinKeyAtLevel<218>(*grouped_master);
 using BasisKey=CKKSAuxiliaryDDRelinKeyNTT<P,162,56,40,20,true>;
 auto basis_cache=std::make_shared<BasisKey::Cache>(*basis_raw);BasisKey basis_product(basis_cache);
 CKKSRelinKeyTensorRNS32<P,BasisKey> basis_tensor{basis_product};
 const auto auto_setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-autoStart).count();

 auto input=std::make_unique<Refresh::Input>();auto slots=std::make_unique<CKKSSlotVector<P>>();const auto inputfile=root/"candidate-input.bin";
 const std::array<uint64_t,5> inputmeta{0x4c4c3134494e3031ULL,P::n,448,85,44};bool reused=std::filesystem::exists(inputfile);
 if(reused){std::ifstream f(inputfile,std::ios::binary);cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,5> meta;std::string stored;ar(meta,stored,*slots,*input);if(meta!=inputmeta || stored!=origin || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("retained input mismatch");}
 else throw std::runtime_error("retained candidate-input.bin required; no input encryption in this replay");
 const auto setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 auto functions=requested=="all"?std::vector<std::string>{"gelu","tanh","sigmoid"}:std::vector<std::string>{requested};
 if(basis_pair||dd_pair||matrix_pair)functions.assign(4,"tanh");
 if(thread_sweep)functions.assign(cpu_profiles.size(),"tanh");
 std::size_t trial_index=0;
 for(const auto& function:functions) {
 const bool use_matrix_tiled=matrix_pair?(trial_index==1||trial_index==2):matrix_tiled;
 const auto *active_profile=thread_sweep?&cpu_profiles[trial_index]:matrix_pair?&(use_matrix_tiled?cpu_profiles[2]:cpu_profiles.front()):nullptr;
 if(active_profile)narrow_benchmark::applyCPUProfile(*active_profile);
 if(matrix_pair)key->narrow.factor_source=use_matrix_tiled?tiled_matrix_source:coarse_matrix;
 const bool use_signed_basis=basis_pair?(trial_index==1||trial_index==2):signed_basis;
 const bool use_dd_wide=dd_pair?(trial_index==1||trial_index==2):dd_wide_rns;
 if(dd_pair)key->narrow.product_tree_evaluator=use_dd_wide?wide_tree:digit_tree;
 const auto artifact=function+((basis_pair||dd_pair||thread_sweep||matrix_pair)?"-trial"+std::to_string(trial_index):"");
 auto target=[&](double x){if(function=="sigmoid")return 1/(1+std::exp(-x));if(function=="gelu")return 0.5*x*(1+std::erf(x/std::sqrt(2.0)));return std::tanh(x);};
 auto curvature=[&](double x){double y=target(x);if(function=="sigmoid")return y*(1-y)*(1-2*y);if(function=="gelu")return (2-x*x)*std::exp(-x*x/2)/std::sqrt(2*std::acos(-1.0));return -2*y*(1-y*y);};
 const uint32_t degree=function=="tanh"?40:12;
 auto approx=function=="tanh"?CKKSBuildNarrowWindowedTrigonometricApproximation(target,[](double x){double y=std::tanh(x);return 1-y*y;},curvature,{1,1,2,4,40},4,16,degree):CKKSBuildNarrowTrigonometricApproximation(target,4,16,degree,function=="gelu"?2.0:1.0,65536,40,2.0,curvature,160.0);
 size_t groups=0;key->narrow.factor_observer=[&](bool tail,size_t group,const auto& timing){++groups;std::cout<<artifact<<" completed "<<(tail?"tail":"head")<<" group="<<group<<" factor_groups="<<groups
  <<" load="<<timing.selector_load_seconds<<" table_ntt="<<timing.table_ntt_seconds
  <<" selector_ntt="<<timing.selector_ntt_seconds<<" matrix="<<timing.matrix_seconds
  <<" inverse_crt="<<timing.inverse_crt_seconds<<" conjugation="<<timing.conjugation_seconds<<std::endl;};
 auto prepared=std::make_unique<Refresh::Coefficients>();auto t=std::chrono::steady_clock::now();if(prepared_s2c_rns)
 CKKSNarrowPrepareCoefficients<Refresh>(*prepared,*input,*key,[&](auto &result,const auto &ct){
  CKKSRealLinearTransform<P,93,52,37>(result,ct,*prepared_s2c_rns,*s2c_input_cache,*s2c_output_cache);
 });
 else if(prepared_s2c)
 CKKSNarrowPrepareCoefficients<Refresh>(*prepared,*input,*key,[&](auto &result,const auto &ct){
  CKKSRealLinearTransform<P,93,52,37>(result,ct,*prepared_s2c,*s2c_input_cache,*s2c_output_cache);
 });
 else CKKSNarrowPrepareCoefficients<Refresh>(*prepared,*input,*key);
 const auto s2c_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();std::cout<<artifact<<" S2C completed seconds="<<s2c_seconds<<std::endl;
 key->narrow.root_stage_observer=[&](const auto&t){std::cout<<artifact<<" stages factor="<<t.leaf_factors_and_inline_alignment_seconds<<" align="<<t.bounded_tail_alignment_seconds<<" dd_cache="<<t.dd_cache_preparation_seconds<<" trees="<<t.dd_trees_and_beta_seconds<<std::endl;};
 t=std::chrono::steady_clock::now();CKKSNarrowFrequencyBasis plan(degree,degree==12?std::vector<uint32_t>{1,3,5,6}:use_signed_basis?std::vector<uint32_t>{6,12,15,16,17,19,20}:std::vector<uint32_t>{1,3,4,9,11,16,17,19,20},use_signed_basis);
 std::vector<CKKSNarrowBlindRotationRequest<Narrow,44>> requests;for(auto freq:plan.frequencies())requests.push_back({prepared.get(),freq,{1,0}});
 std::vector<Narrow::Output> roots;CKKSNarrowBlindRotateBatch<Narrow,44>(roots,requests,key->narrow,plan.frequencies().size());
 auto weighted=std::make_unique<Basis::Output>();CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*weighted,roots,plan,approx.coefficients,basis_tensor,
 [&](auto& result,const auto& input){key->narrow.automorphism_evaluator(result,input,125,2*P::n-1);},
 [&](auto& result,const auto& input){key->narrow.automorphism_evaluator(result,input,162,2*P::n-1);});const auto blind_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
 auto out=std::make_unique<Output>();CKKSScaleDown<P,125,77,33>(*out,*weighted);auto square=std::make_unique<Square>();auto postkey=makeCKKSCompactDDRelinKeyAtLevel<148>(*grouped_master);using PostKey=CKKSAuxiliaryDDRelinKeyNTT<P,92,56,40,20,true>;auto postcache=std::make_shared<PostKey::Cache>(*postkey);PostKey auxiliary(postcache);CKKSRelinKeyTensorRNS32<P,PostKey> post_tensor{auxiliary};CKKSMultRelinThenRescale<P>(*square,*out,*out,post_tensor);
 auto decoded=std::make_unique<CKKSSlotVector<P>>();ckksSlotDecrypt<P>(*decoded,*out,*dense);double error=0,amplitude=0;
 for(size_t i=0;i<P::n/2;++i){const double v=target((*slots)[i].real());error=narrow_benchmark::finiteMaximum(error,std::abs((*decoded)[i]-std::complex<double>{v,0}));amplitude=narrow_benchmark::finiteMaximum(amplitude,std::abs(v));}
 ckksSlotDecrypt<P>(*decoded,*square,*dense);double product_error=0;for(size_t i=0;i<P::n/2;++i){double v=target((*slots)[i].real());product_error=narrow_benchmark::finiteMaximum(product_error,std::abs((*decoded)[i]-std::complex<double>{v*v,0}));}
 const auto tolerance=std::ldexp(1.0,-23),product_tolerance=(2*amplitude+1)*tolerance+tolerance*tolerance;bool passed=error<tolerance && product_error<product_tolerance;
 {std::ofstream f(root/(artifact+"-outputs.bin"),std::ios::binary);cereal::PortableBinaryOutputArchive ar(f);const std::array<uint64_t,9> meta{0x4c4c313446425331ULL,P::n,448,85,44,48,92,44,48};ar(meta,function,*slots,*input,*prepared,*out,*square);}
 {std::ofstream f(root/(artifact+"-report.json"));f<<std::setprecision(17)<<"{\"function\":\""<<function<<"\",\"n\":"<<P::n<<",\"complex_input\":false,\"ntt_dd\":true,\"parallel_product_pairs\":"<<(parallel_pairs?"true":"false")<<",\"small_prime_factors\":"<<(small_prime_factors?"true":"false")<<",\"dd_small_prime_backend\":"<<((use_dd_wide||CKKSPreparedDDRelinKeyNTT<P,433,40,20>::uses_small_primes)?"true":"false")<<",\"dd_convolution_primes\":"<<(use_dd_wide?CKKSPreparedGadgetRNS32<P,433,40>::prime_count:CKKSPreparedDDRelinKeyNTT<P,433,40,20>::prime_count)<<",\"dd_cache_chunk_bits\":"<<(use_dd_wide?80:20)<<",\"profile\":\"mixed43-product37-public-reuse\",\"selector_log_q\":433,\"selector_log_delta\":57,\"factor_log_q\":377,\"factor_log_delta\":43,\"dd_primary_bits\":40,\"dd_secondary_bits\":20,\"dd_auxiliary_bits\":56,\"public_derived_ntt_automorphisms\":true,\"conjugation_primary_bits\":20,\"basis_shares_public_conjugation_cache\":true,\"automorphism_setup_seconds\":"<<auto_setup_seconds<<",\"auto_rns\":"<<(auto_rns?"true":"false")<<",\"tensor_rns\":"<<(tensor_rns?"true":"false")<<",\"direct_rns_factors\":"<<(direct_rns?"true":"false")<<",\"selector_residue_rows\":"<<(direct_rns?17:36)<<",\"omp_threads\":"<<narrow_benchmark::benchmarkThreads()<<",\"thread_sweep_trial\":"<<(thread_sweep?std::to_string(trial_index):"null")<<",\"pinned_cpu_ids\":"<<(active_profile?narrow_benchmark::CPUProfileJSON(*active_profile):"null")<<",\"matrix_tiled\":"<<((direct_rns&&use_matrix_tiled)?"true":"false")<<",\"matrix_pair_trial\":"<<(matrix_pair?std::to_string(trial_index):"null")<<",\"dd_wide_rns\":"<<(use_dd_wide?"true":"false")<<",\"dd_pair_trial\":"<<(dd_pair?std::to_string(trial_index):"null")<<",\"basis_pair_trial\":"<<(basis_pair?std::to_string(trial_index):"null")<<",\"signed_frequency_basis\":"<<(plan.usesDifferences()?"true":"false")<<",\"basis_root_count\":"<<plan.frequencies().size()<<",\"frequency_basis\":true,\"selector_cache_retained_bytes\":"<<selector_retained_bytes<<",\"degree\":"<<degree<<",\"selector_matrix_digit_bits\":"<<(direct_rns?0:52)<<",\"selector_chunk_size\":40,\"fresh_private_keys_generated\":false,\"new_evaluation_keys_generated\":0,\"input_reused\":"<<(reused?"true":"false")<<",\"key_load_and_input_setup_seconds\":"<<setup_seconds<<",\"s2c_rns\":"<<(prepared_s2c_rns?"true":"false")<<",\"s2c_rns_primes\":"<<(prepared_s2c_rns?5:0)<<",\"s2c_public_setup_seconds\":"<<s2c_setup_seconds<<",\"s2c_prepared_direct_groups\":"<<(prepared_s2c_rns?prepared_s2c_rns->direct.groups.size():prepared_s2c?prepared_s2c->direct.groups.size():0)<<",\"s2c_prepared_conjugate_groups\":"<<(prepared_s2c_rns?prepared_s2c_rns->conjugate.groups.size():prepared_s2c?prepared_s2c->conjugate.groups.size():0)<<",\"s2c_retained_spectrum_bytes\":"<<s2c_retained_bytes<<",\"s2c_public_keys_prepared\":"<<(have_prepared_s2c?s2c_input_cache->prepared_keys+s2c_output_cache->prepared_keys:0)<<",\"s2c_seconds\":"<<s2c_seconds<<",\"blind_rotation_and_series_seconds\":"<<blind_seconds<<",\"functional_evaluation_seconds\":"<<s2c_seconds+blind_seconds<<",\"output_log_q\":92,\"output_log_delta\":44,\"post_product_log_q\":48,\"max_error\":"<<error<<",\"post_product_max_error\":"<<product_error<<",\"product_tolerance\":"<<product_tolerance<<",\"accuracy_passed\":"<<(passed?"true":"false")<<",\"security_validated\":false}";}
 std::cout<<artifact<<" error="<<error<<" square error="<<product_error<<" seconds="<<s2c_seconds+blind_seconds<<std::endl;if(!passed)throw std::runtime_error("functional accuracy gate failed");
 ++trial_index;
}
} catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
