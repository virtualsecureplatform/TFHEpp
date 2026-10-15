// Experimental direct-selector real-input frequency-basis replay using retained LL14
// component keys. Never reconstruct selectors or generate private/evaluation
// keys here. The supplied public key artifacts must match the checkpoint.
#include "ckks_narrow_basis_benchmark.hpp"
#include "ckks_benchmark_affinity.hpp"
#include <ckks/narrow_basis_prepared.hpp>
#include <ckks/prepared_compact_linear_rns32.hpp>
#include <ckks/narrow_dd_ntt.hpp>
#include <ckks/narrow_auto_ntt.hpp>
#include <ckks/native_auto_file.hpp>
#include <set>
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
 if(argc<8 || (argc-8)%2)throw std::runtime_error("usage: checkpoint-directory switch56-file sigmoid|gelu|tanh|all output-directory direct440-selectors saved-dd433-primary20 saved-packed435 [--cache-gib 90] [--s2c-cache-gib 9] [--direct-rns 1] [--s2c-rns 1] [--tensor-rns 1] [--auto-rns 1] [--signed-basis 1] [--basis-pair 0] [--dd-wide-rns 1] [--dd-pair 0] [--dd-primary-pair 0] [--parallel-conjugation 0] [--conjugation-pair 0] [--parallel-tail-alignment 0] [--tail-alignment-pair 0] [--parallel-root-columns 0] [--root-columns-pair 0] [--auto-primary-bits 20] [--auto-primary-pair 0] [--native-auto-directory PATH] [--native-auto-pair 0] [--batch-product-trees 0] [--batch-trees-pair 0] [--thread-sweep 0] [--matrix-tiled 1] [--matrix-pair 0] [--matrix-dot8 1] [--matrix-dot8-pair 0] [--factor-short-ntt 1] [--factor-short-ntt-pair 0]");
 // Public primary40 fails this retained precision profile; keep it experimental.
 std::string native_auto_directory;
 uint32_t auto_primary_bits=20;
 uint64_t cache_gib=90,s2c_cache_gib=9;bool direct_rns=true,s2c_rns=true,tensor_rns=true,auto_rns=true,signed_basis=true,basis_pair=false,dd_wide_rns=true,dd_pair=false,dd_primary_pair=false,parallel_conjugation=false,conjugation_pair=false,parallel_tail_alignment=false,tail_alignment_pair=false,parallel_root_columns=false,root_columns_pair=false,auto_primary_pair=false,native_auto_pair=false,batch_product_trees=false,batch_trees_pair=false,thread_sweep=false,matrix_tiled=true,matrix_pair=false,matrix_dot8=true,matrix_dot8_pair=false,factor_short_ntt=true,factor_short_ntt_pair=false;
 for(int i=8;i<argc;i+=2){
  const std::string option=argv[i];
  if(option=="--native-auto-directory"){native_auto_directory=argv[i+1];continue;}
  const auto value=std::stoull(argv[i+1]);
  if(value>1024)throw std::invalid_argument("cache budget exceeds supported range");
  if(option=="--cache-gib")cache_gib=value;
  else if(option=="--s2c-cache-gib")s2c_cache_gib=value;
  else if(option=="--matrix-tiled"){if(value>1)throw std::invalid_argument("matrix tiled flag must be 0 or 1");matrix_tiled=value;}
  else if(option=="--factor-short-ntt"){if(value>1)throw std::invalid_argument("factor short NTT flag must be 0 or 1");factor_short_ntt=value;}
  else if(option=="--factor-short-ntt-pair"){if(value>1)throw std::invalid_argument("factor short NTT pair flag must be 0 or 1");factor_short_ntt_pair=value;}
  else if(option=="--matrix-dot8"){if(value>1)throw std::invalid_argument("matrix dot8 flag must be 0 or 1");matrix_dot8=value;}
  else if(option=="--matrix-dot8-pair"){if(value>1)throw std::invalid_argument("matrix dot8 pair flag must be 0 or 1");matrix_dot8_pair=value;}
  else if(option=="--matrix-pair"){if(value>1)throw std::invalid_argument("matrix pair flag must be 0 or 1");matrix_pair=value;}
  else if(option=="--thread-sweep"){if(value>1)throw std::invalid_argument("thread sweep flag must be 0 or 1");thread_sweep=value;}
  else if(option=="--dd-wide-rns"){if(value>1)throw std::invalid_argument("DD wide RNS flag must be 0 or 1");dd_wide_rns=value;}
  else if(option=="--batch-product-trees"){if(value>1)throw std::invalid_argument("batch product trees flag must be 0 or 1");batch_product_trees=value;}
  else if(option=="--batch-trees-pair"){if(value>1)throw std::invalid_argument("batch trees pair flag must be 0 or 1");batch_trees_pair=value;}
  else if(option=="--native-auto-pair"){if(value>1)throw std::invalid_argument("native auto pair flag must be 0 or 1");native_auto_pair=value;}
  else if(option=="--auto-primary-bits"){if(value!=20&&value!=40)throw std::invalid_argument("automorphism primary bits must be 20 or 40");auto_primary_bits=value;}
  else if(option=="--auto-primary-pair"){if(value>1)throw std::invalid_argument("automorphism primary pair flag must be 0 or 1");auto_primary_pair=value;}
  else if(option=="--parallel-root-columns"){if(value>1)throw std::invalid_argument("parallel root columns flag must be 0 or 1");parallel_root_columns=value;}
  else if(option=="--root-columns-pair"){if(value>1)throw std::invalid_argument("root columns pair flag must be 0 or 1");root_columns_pair=value;}
  else if(option=="--parallel-tail-alignment"){if(value>1)throw std::invalid_argument("parallel tail alignment flag must be 0 or 1");parallel_tail_alignment=value;}
  else if(option=="--tail-alignment-pair"){if(value>1)throw std::invalid_argument("tail alignment pair flag must be 0 or 1");tail_alignment_pair=value;}
  else if(option=="--parallel-conjugation"){if(value>1)throw std::invalid_argument("parallel conjugation flag must be 0 or 1");parallel_conjugation=value;}
  else if(option=="--conjugation-pair"){if(value>1)throw std::invalid_argument("conjugation pair flag must be 0 or 1");conjugation_pair=value;}
  else if(option=="--dd-primary-pair"){if(value>1)throw std::invalid_argument("DD primary pair flag must be 0 or 1");dd_primary_pair=value;}
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
 if((basis_pair||dd_pair||dd_primary_pair||conjugation_pair||tail_alignment_pair||root_columns_pair||auto_primary_pair||native_auto_pair||batch_trees_pair||thread_sweep||matrix_pair)&&requested!="tanh")throw std::invalid_argument("paired replay requires tanh");
 if(int(basis_pair)+int(dd_pair)+int(dd_primary_pair)+int(conjugation_pair)+int(tail_alignment_pair)+int(root_columns_pair)+int(auto_primary_pair)+int(native_auto_pair)+int(batch_trees_pair)+int(thread_sweep)+int(matrix_pair)+int(matrix_dot8_pair)+int(factor_short_ntt_pair)>1)throw std::invalid_argument("paired replay controls must be separate");
 if((parallel_conjugation||conjugation_pair)&&!direct_rns)throw std::invalid_argument("parallel factor conjugation requires direct RNS factors");
 if((auto_primary_pair||auto_primary_bits==40)&&!auto_rns)throw std::invalid_argument("wider public automorphisms require direct RNS gadgets");
 if(batch_product_trees&&(dd_pair||dd_primary_pair||root_columns_pair))throw std::invalid_argument("batched trees require a separate tree-backend comparison");
 if(native_auto_pair&&native_auto_directory.empty())throw std::invalid_argument("native auto pair requires a saved key directory");
 if(!native_auto_directory.empty()&&(!auto_rns||auto_primary_pair||auto_primary_bits!=20))throw std::invalid_argument("native auto mode requires RNS and the default public width selection");
 if(dd_primary_pair&&!dd_wide_rns)throw std::invalid_argument("DD primary pair requires wide RNS gadgets");
 if(factor_short_ntt_pair && requested=="all")throw std::invalid_argument("factor short NTT pair requires a single function");
 if(factor_short_ntt_pair&&!direct_rns)throw std::invalid_argument("factor short NTT pair requires direct RNS factors");
 if(matrix_dot8_pair && requested=="all")throw std::invalid_argument("matrix dot8 pair requires a single function");
 if(matrix_dot8_pair&&!direct_rns)throw std::invalid_argument("matrix dot8 comparison requires direct RNS factors");
 if(matrix_pair&&!direct_rns)throw std::invalid_argument("matrix pair requires direct RNS factors");
#ifdef _OPENMP
 if((thread_sweep||matrix_pair)&&omp_get_proc_bind()!=omp_proc_bind_close)
  throw std::invalid_argument("core replays require OMP_PROC_BIND=close with ordered singleton CPU places");
#endif
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
  uint64_t head,tail;ar(head,tail);if(head || tail)throw std::runtime_error("expected file-backed selectors");auto original=std::make_unique<CKKSGaloisKey<P,345>>();ar(*original);CKKSReduceGaloisKey<P,93>(key->s2c_input_galois,*original);CKKSReduceGaloisKey<P,56>(key->s2c_output_galois,*original);CKKSReduceAutoKey<P,178>(key->narrow.output_conjugation,(*original)[P::nbit]);}
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
 using File=CKKSNarrowSelectorFile<TableSource>;auto raw=std::make_shared<File>(argv[5],File::Existing{});
 uint64_t selector_retained_bytes=0;
 decltype(key->narrow.factor_source) coarse_matrix,tiled_matrix_source,serial_conjugation_source,parallel_conjugation_source,four_term_source,eight_term_source,scalar_short_source,simd_short_source;
 if(direct_rns){
  using Cache=CKKSNarrowFactorSelectorCacheRNS32<TableSource>;
  std::shared_ptr<const Cache> retained;
  if(cache_gib){retained=std::make_shared<Cache>(raw,cache_gib<<30,40,
   [&](bool tail,size_t group){std::cout<<"prepared direct RNS cache "<<(tail?"tail":"head")<<" group="<<group<<std::endl;});selector_retained_bytes=Cache::retained_bytes;}
  using Source=CKKSNarrowDeferredRescaleFactorSourceRNS32<Narrow,TableSource,TableSource>;
  if(factor_short_ntt_pair){
   scalar_short_source=std::make_shared<Source>(raw,retained,matrix_tiled,parallel_conjugation,matrix_dot8,false);
   simd_short_source=std::make_shared<Source>(raw,retained,matrix_tiled,parallel_conjugation,matrix_dot8,true);
   key->narrow.factor_source=scalar_short_source;
  }
  else if(matrix_dot8_pair){
   four_term_source=std::make_shared<Source>(raw,retained,matrix_tiled,parallel_conjugation,false,factor_short_ntt);
   eight_term_source=std::make_shared<Source>(raw,retained,matrix_tiled,parallel_conjugation,true,factor_short_ntt);
   key->narrow.factor_source=four_term_source;
  }
  else if(matrix_pair){coarse_matrix=std::make_shared<Source>(raw,retained,false,parallel_conjugation,matrix_dot8,factor_short_ntt);tiled_matrix_source=std::make_shared<Source>(raw,retained,true,parallel_conjugation,matrix_dot8,factor_short_ntt);key->narrow.factor_source=coarse_matrix;}
  else if(conjugation_pair){
   serial_conjugation_source=std::make_shared<Source>(raw,retained,matrix_tiled,false,matrix_dot8,factor_short_ntt);
   parallel_conjugation_source=std::make_shared<Source>(raw,retained,matrix_tiled,true,matrix_dot8,factor_short_ntt);
   key->narrow.factor_source=serial_conjugation_source;
  }
  else key->narrow.factor_source=std::make_shared<Source>(raw,retained,matrix_tiled,parallel_conjugation,matrix_dot8,factor_short_ntt);
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
 decltype(key->narrow.product_tree_evaluator) digit_tree,wide_tree,primary20_tree,primary40_tree;
 const auto prepare_tree=[&](bool wide){
  if(wide){
   if(tensor_rns)CKKSNarrowUseAuxiliaryNTTDD<48,20,true,true,true>(key->narrow,*master);
   else CKKSNarrowUseAuxiliaryNTTDD<48,20,true,false,true>(key->narrow,*master);
  }else{
   if(tensor_rns)CKKSNarrowUseAuxiliaryNTTDD<48,20,true,true>(key->narrow,*master);
   else CKKSNarrowUseAuxiliaryNTTDD<48,20,true>(key->narrow,*master);
  }
 };
 if(dd_pair){prepare_tree(false);digit_tree=key->narrow.product_tree_evaluator;prepare_tree(true);wide_tree=key->narrow.product_tree_evaluator;}
 else prepare_tree(dd_wide_rns);
 if(dd_primary_pair){
  primary20_tree=key->narrow.product_tree_evaluator;
  auto grouped=makeCKKSCompactDDRelinKeyWithPrimaryBase<40>(*master);
  if(tensor_rns)CKKSNarrowUseAuxiliaryNTTDD<48,20,true,true,true>(key->narrow,*grouped);
  else CKKSNarrowUseAuxiliaryNTTDD<48,20,true,false,true>(key->narrow,*grouped);
  primary40_tree=key->narrow.product_tree_evaluator;
  key->narrow.product_tree_evaluator=primary20_tree;
 }
 auto autoStart=std::chrono::steady_clock::now();
 decltype(key->narrow.automorphism_evaluator) primary20_auto,primary40_auto,native_auto;
 if(!native_auto_directory.empty()){
  using Cache=PreparedPublicAutoNTT<P,435,44,20,true,112>;
  auto cached=std::make_shared<std::map<uint32_t,std::shared_ptr<Cache>>>();
  std::set<uint32_t> exponents{2*P::n-1};
  for(std::size_t ell=0;ell<Narrow::tail_blocks_per_pack;++ell)exponents.insert(ckks_detail::rotationAutomorphismForSteps<P>(2*(ell+1)*Narrow::tail_width));
  for(auto exponent:exponents){
   std::ifstream file(std::filesystem::path(native_auto_directory)/(std::to_string(exponent)+".bin"),std::ios::binary);
   if(!file)throw std::runtime_error("saved native automorphism key missing");
   auto rows=CKKSReadNativeAutoRows<P,435,44>(file,exponent,origin,3.2);
   cached->emplace(exponent,std::make_shared<Cache>(*rows));
  }
  native_auto=narrow_detail::makePublicAutoEvaluator<Narrow,398,121>(cached);
 }
 if(native_auto_directory.empty()||native_auto_pair){
 auto packed=CKKSNarrowLoadPackedKey<PackedSource>(argv[7]);
 if(auto_primary_pair){
  using Cache40=PreparedPublicAutoNTT<P,435,40,20,true>;
  auto regrouped=std::make_shared<std::map<uint32_t,std::shared_ptr<Cache40>>>();
  CKKSNarrowUsePublicNTTAutomorphisms<398,20,20,true,121>(key->narrow,*packed,[&](const auto&rows){
   auto wider=makePublicAutoGadgetRowsWithPrimaryBase<40>(rows);
   regrouped->emplace(rows.exponent,std::make_shared<Cache40>(*wider));
  });
  primary20_auto=key->narrow.automorphism_evaluator;
  primary40_auto=narrow_detail::makePublicAutoEvaluator<Narrow,398,121>(regrouped);
 }
 else if(auto_primary_bits==40)CKKSNarrowUsePublicNTTAutomorphisms<398,40,20,true,121>(key->narrow,*packed);
 else if(auto_rns)CKKSNarrowUsePublicNTTAutomorphisms<398,20,20,true,121>(key->narrow,*packed);
 else CKKSNarrowUsePublicNTTAutomorphisms<398,20,20,false,121>(key->narrow,*packed);
 if(native_auto_pair)primary20_auto=key->narrow.automorphism_evaluator;
 }else key->narrow.automorphism_evaluator=native_auto;
 auto basis_raw=makeCKKSCompactDDRelinKeyAtLevel<213>(*master);auto basis_ntt=makeCKKSDDRelinKeyNTT<20>(*basis_raw);
 CKKSAuxiliaryDDRelinKeyNTT<P,165,48,20,20> basis_product(basis_ntt);
 const auto auto_setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-autoStart).count();

 auto input=std::make_unique<Refresh::Input>();auto slots=std::make_unique<CKKSSlotVector<P>>();const auto inputfile=root/"candidate-input.bin";
 const std::array<uint64_t,5> inputmeta{0x4c4c3134494e3031ULL,P::n,448,85,44};bool reused=std::filesystem::exists(inputfile);
 if(reused){std::ifstream f(inputfile,std::ios::binary);cereal::PortableBinaryInputArchive ar(f);std::array<uint64_t,5> meta;std::string stored;ar(meta,stored,*slots,*input);if(meta!=inputmeta || stored!=origin || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("retained input mismatch");}
 else throw std::runtime_error("retained candidate-input.bin required; no input encryption in this replay");
 const auto setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 auto functions=requested=="all"?std::vector<std::string>{"gelu","tanh","sigmoid"}:std::vector<std::string>{requested};
 if(basis_pair||dd_pair||dd_primary_pair||conjugation_pair||tail_alignment_pair||root_columns_pair||auto_primary_pair||native_auto_pair||batch_trees_pair||matrix_pair)functions.assign(4,"tanh");
 if(auto_primary_pair||native_auto_pair||batch_trees_pair)functions.assign(5,"tanh"); // Warm-up, then 20/40/40/20.
 if(factor_short_ntt_pair)functions.assign(5,requested); // Warm-up, then scalar/SIMD/SIMD/scalar.
 if(matrix_dot8_pair)functions.assign(5,requested); // Warm-up, then 4/8/8/4 terms.
 if(thread_sweep)functions.assign(cpu_profiles.size()+1,"tanh"); // Warm-up, then the symmetric core sweep.
 const auto sweep_tree_evaluator=key->narrow.product_tree_evaluator;
 std::size_t trial_index=0;
 for(const auto& function:functions) {
 const bool factor_short_ntt_warmup=factor_short_ntt_pair&&trial_index==0;
 const bool use_factor_short_ntt=factor_short_ntt_pair?(trial_index==2||trial_index==3):factor_short_ntt;
 if(factor_short_ntt_pair)key->narrow.factor_source=use_factor_short_ntt?simd_short_source:scalar_short_source;
 const bool matrix_dot8_warmup=matrix_dot8_pair&&trial_index==0;
 const bool use_matrix_dot8=matrix_dot8_pair?(trial_index==2||trial_index==3):matrix_dot8;
 if(matrix_dot8_pair)key->narrow.factor_source=use_matrix_dot8?eight_term_source:four_term_source;
 const bool batch_trees_warmup=batch_trees_pair&&trial_index==0;
 const bool use_batch_product_trees=batch_trees_pair?(trial_index==2||trial_index==3):batch_product_trees;
 key->narrow.batch_product_trees=use_batch_product_trees;
 const bool thread_sweep_warmup=thread_sweep&&trial_index==0;
 const auto sweep_profile_index=thread_sweep&&trial_index?trial_index-1:std::size_t{0};
 narrow_benchmark::CPURecorder tree_entries;
 if(thread_sweep)key->narrow.product_tree_evaluator=[base=sweep_tree_evaluator,&tree_entries](auto&result,auto factors,const auto&raw){
  tree_entries.record();base(result,std::move(factors),raw);
 };
 const bool native_auto_warmup=native_auto_pair&&trial_index==0;
 const bool use_native_auto=native_auto_pair?(trial_index==2||trial_index==3):!native_auto_directory.empty();
 if(native_auto_pair)key->narrow.automorphism_evaluator=use_native_auto?native_auto:primary20_auto;
 const bool auto_primary_warmup=auto_primary_pair&&trial_index==0;
 const bool use_auto_primary40=auto_primary_pair?(trial_index==2||trial_index==3):auto_primary_bits==40;
 if(auto_primary_pair)key->narrow.automorphism_evaluator=use_auto_primary40?primary40_auto:primary20_auto;
 const bool use_parallel_root_columns=root_columns_pair?(trial_index==1||trial_index==2):parallel_root_columns;
 key->narrow.parallel_root_columns=use_parallel_root_columns;
 const bool use_parallel_tail_alignment=tail_alignment_pair?(trial_index==1||trial_index==2):parallel_tail_alignment;
 key->narrow.parallel_public_tail_alignment=use_parallel_tail_alignment;
 const bool use_parallel_conjugation=conjugation_pair?(trial_index==1||trial_index==2):parallel_conjugation;
 if(conjugation_pair)key->narrow.factor_source=use_parallel_conjugation?parallel_conjugation_source:serial_conjugation_source;
 const bool use_matrix_tiled=matrix_pair?(trial_index==1||trial_index==2):matrix_tiled;
 const auto *active_profile=thread_sweep?&cpu_profiles[sweep_profile_index]:matrix_pair?&(use_matrix_tiled?cpu_profiles[2]:cpu_profiles.front()):nullptr;
 if(active_profile)narrow_benchmark::applyCPUProfile(*active_profile);
 if(matrix_pair)key->narrow.factor_source=use_matrix_tiled?tiled_matrix_source:coarse_matrix;
 const bool use_signed_basis=basis_pair?(trial_index==1||trial_index==2):signed_basis;
 const bool use_primary40=dd_primary_pair&&(trial_index==1||trial_index==2);
 if(dd_primary_pair)key->narrow.product_tree_evaluator=use_primary40?primary40_tree:primary20_tree;
 const bool use_dd_wide=dd_pair?(trial_index==1||trial_index==2):dd_wide_rns;
 if(dd_pair)key->narrow.product_tree_evaluator=use_dd_wide?wide_tree:digit_tree;
 const auto artifact=function+((basis_pair||dd_pair||dd_primary_pair||conjugation_pair||tail_alignment_pair||root_columns_pair||auto_primary_pair||native_auto_pair||batch_trees_pair||thread_sweep||matrix_pair||matrix_dot8_pair||factor_short_ntt_pair)?"-trial"+std::to_string(trial_index):"");
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
 std::vector<uint32_t> observed_roots;observed_roots.reserve(plan.frequencies().size());
 const CKKSNarrowRootObserver<Narrow,44> root_observer=[&](const auto&request,const auto&){observed_roots.push_back(request.frequency);};
 const bool check_root_observer=root_columns_pair||parallel_root_columns||batch_product_trees||batch_trees_pair;
 std::vector<Narrow::Output> roots;CKKSNarrowBlindRotateBatch<Narrow,44>(roots,requests,key->narrow,plan.frequencies().size(),check_root_observer?&root_observer:nullptr);
 if(check_root_observer&&observed_roots!=plan.frequencies())throw std::runtime_error("root observer order changed");
 auto weighted=std::make_unique<Basis::Output>();CKKSNarrowEvalSeriesFromBasisRootsPrepared<Basis>(*weighted,roots,plan,approx.coefficients,basis_product,
 [&](auto& result,const auto& input){key->narrow.automorphism_evaluator(result,input,121,2*P::n-1);},
 [&](auto& result,const auto& input){key->narrow.automorphism_evaluator(result,input,165,2*P::n-1);});const auto blind_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
 auto out=std::make_unique<Output>();CKKSScaleDown<P,121,72,28>(*out,*weighted);auto square=std::make_unique<Square>();auto postkey=makeCKKSCompactDDRelinKeyAtLevel<141>(*master);auto postntt=makeCKKSDDRelinKeyNTT<20>(*postkey);CKKSAuxiliaryDDRelinKeyNTT<P,93,48,20,20> auxiliary(postntt);CKKSMultRelinThenRescale<P>(*square,*out,*out,auxiliary);
 auto decoded=std::make_unique<CKKSSlotVector<P>>();ckksSlotDecrypt<P>(*decoded,*out,*dense);double error=0,amplitude=0;
 for(size_t i=0;i<P::n/2;++i){const double v=target((*slots)[i].real());error=narrow_benchmark::finiteMaximum(error,std::abs((*decoded)[i]-std::complex<double>{v,0}));amplitude=narrow_benchmark::finiteMaximum(amplitude,std::abs(v));}
 ckksSlotDecrypt<P>(*decoded,*square,*dense);double product_error=0;for(size_t i=0;i<P::n/2;++i){double v=target((*slots)[i].real());product_error=narrow_benchmark::finiteMaximum(product_error,std::abs((*decoded)[i]-std::complex<double>{v*v,0}));}
 const auto tolerance=std::ldexp(1.0,-23),product_tolerance=(2*amplitude+1)*tolerance+tolerance*tolerance;bool passed=error<tolerance && product_error<product_tolerance;
 {std::ofstream f(root/(artifact+"-outputs.bin"),std::ios::binary);cereal::PortableBinaryOutputArchive ar(f);const std::array<uint64_t,9> meta{0x4c4c313446425331ULL,P::n,448,85,44,48,93,44,49};ar(meta,function,*slots,*input,*prepared,*out,*square);}
 {std::ofstream f(root/(artifact+"-report.json"));f<<std::setprecision(17)<<"{\"function\":\""<<function<<"\",\"n\":"<<P::n<<",\"complex_input\":false,\"ntt_dd\":true,\"parallel_product_pairs\":"<<(parallel_pairs?"true":"false")<<",\"small_prime_factors\":"<<(small_prime_factors?"true":"false")<<",\"dd_small_prime_backend\":"<<(CKKSPreparedDDRelinKeyNTT<P,433,20,20>::uses_small_primes?"true":"false")<<",\"dd_convolution_primes\":"<<(use_dd_wide?(use_primary40?CKKSPreparedGadgetRNS32<P,433,40>::prime_count:CKKSPreparedGadgetRNS32<P,433,20>::prime_count):CKKSPreparedDDRelinKeyNTT<P,433,20,20>::prime_count)<<",\"dd_cache_chunk_bits\":"<<(use_dd_wide?80:20)<<",\"dd_primary_bits\":"<<(use_primary40?40:20)<<",\"basis_dd_primary_bits\":20,\"dd_secondary_bits\":20,\"dd_auxiliary_bits\":48,\"public_derived_ntt_automorphisms\":"<<(use_native_auto?"false":"true")<<",\"native_auto_keys_loaded\":"<<(use_native_auto?"true":"false")<<",\"native_auto_pair_warmup\":"<<(native_auto_warmup?"true":"false")<<",\"native_auto_pair_trial\":"<<((native_auto_pair&&!native_auto_warmup)?std::to_string(trial_index-1):"null")<<",\"conjugation_primary_bits\":"<<(use_native_auto?44:use_auto_primary40?40:20)<<",\"public_auto_gadget_rows\":"<<(use_native_auto?10:use_auto_primary40?11:22)<<",\"public_auto_rns_primes\":"<<(auto_rns?(use_native_auto?CKKSPreparedGadgetRNS32<P,435,44,112>::prime_count:use_auto_primary40?CKKSPreparedGadgetRNS32<P,435,40>::prime_count:CKKSPreparedGadgetRNS32<P,435,20>::prime_count):0)<<",\"auto_primary_pair_warmup\":"<<(auto_primary_warmup?"true":"false")<<",\"auto_primary_pair_trial\":"<<((auto_primary_pair&&!auto_primary_warmup)?std::to_string(trial_index-1):"null")<<",\"public_auto_chunk_bits\":"<<(use_native_auto?112:80)<<",\"basis_shares_public_conjugation_cache\":true,\"automorphism_setup_seconds\":"<<auto_setup_seconds<<",\"auto_rns\":"<<(auto_rns?"true":"false")<<",\"tensor_rns\":"<<(tensor_rns?"true":"false")<<",\"direct_rns_factors\":"<<(direct_rns?"true":"false")<<",\"selector_residue_rows\":"<<(direct_rns?17:36)<<",\"omp_threads\":"<<narrow_benchmark::benchmarkThreads()<<",\"thread_sweep_warmup\":"<<(thread_sweep_warmup?"true":"false")<<",\"tree_entry_cpu_ids\":"<<((thread_sweep&&!use_batch_product_trees)?narrow_benchmark::CPUProfileJSON(tree_entries.cpus()):"null")<<",\"thread_sweep_trial\":"<<((thread_sweep&&!thread_sweep_warmup)?std::to_string(trial_index-1):"null")<<",\"pinned_cpu_ids\":"<<(active_profile?narrow_benchmark::CPUProfileJSON(*active_profile):"null")<<",\"matrix_tiled\":"<<((direct_rns&&use_matrix_tiled)?"true":"false")<<",\"factor_short_ntt\":"<<(use_factor_short_ntt?"true":"false")<<",\"factor_short_ntt_pair_warmup\":"<<(factor_short_ntt_warmup?"true":"false")<<",\"factor_short_ntt_pair_trial\":"<<((factor_short_ntt_pair&&!factor_short_ntt_warmup)?std::to_string(trial_index-1):"null")<<",\"matrix_dot8\":"<<(use_matrix_dot8?"true":"false")<<",\"matrix_dot8_pair_warmup\":"<<(matrix_dot8_warmup?"true":"false")<<",\"matrix_dot8_pair_trial\":"<<((matrix_dot8_pair&&!matrix_dot8_warmup)?std::to_string(trial_index-1):"null")<<",\"matrix_pair_trial\":"<<(matrix_pair?std::to_string(trial_index):"null")<<",\"dd_wide_rns\":"<<(use_dd_wide?"true":"false")<<",\"batch_product_trees\":"<<(use_batch_product_trees?"true":"false")<<",\"batch_trees_pair_warmup\":"<<(batch_trees_warmup?"true":"false")<<",\"batch_trees_pair_trial\":"<<((batch_trees_pair&&!batch_trees_warmup)?std::to_string(trial_index-1):"null")<<",\"parallel_root_columns\":"<<(use_parallel_root_columns?"true":"false")<<",\"root_observer_order_validated\":"<<(check_root_observer?"true":"false")<<",\"root_columns_pair_trial\":"<<(root_columns_pair?std::to_string(trial_index):"null")<<",\"parallel_public_tail_alignment\":"<<(use_parallel_tail_alignment?"true":"false")<<",\"tail_alignment_pair_trial\":"<<(tail_alignment_pair?std::to_string(trial_index):"null")<<",\"parallel_factor_conjugation\":"<<(use_parallel_conjugation?"true":"false")<<",\"conjugation_pair_trial\":"<<(conjugation_pair?std::to_string(trial_index):"null")<<",\"dd_primary_pair_trial\":"<<(dd_primary_pair?std::to_string(trial_index):"null")<<",\"dd_pair_trial\":"<<(dd_pair?std::to_string(trial_index):"null")<<",\"basis_pair_trial\":"<<(basis_pair?std::to_string(trial_index):"null")<<",\"signed_frequency_basis\":"<<(plan.usesDifferences()?"true":"false")<<",\"basis_root_count\":"<<plan.frequencies().size()<<",\"frequency_basis\":true,\"selector_cache_retained_bytes\":"<<selector_retained_bytes<<",\"degree\":"<<degree<<",\"selector_matrix_digit_bits\":"<<(direct_rns?0:52)<<",\"selector_chunk_size\":40,\"fresh_private_keys_generated\":false,\"new_evaluation_keys_generated\":0,\"input_reused\":"<<(reused?"true":"false")<<",\"key_load_and_input_setup_seconds\":"<<setup_seconds<<",\"s2c_rns\":"<<(prepared_s2c_rns?"true":"false")<<",\"s2c_rns_primes\":"<<(prepared_s2c_rns?5:0)<<",\"s2c_public_setup_seconds\":"<<s2c_setup_seconds<<",\"s2c_prepared_direct_groups\":"<<(prepared_s2c_rns?prepared_s2c_rns->direct.groups.size():prepared_s2c?prepared_s2c->direct.groups.size():0)<<",\"s2c_prepared_conjugate_groups\":"<<(prepared_s2c_rns?prepared_s2c_rns->conjugate.groups.size():prepared_s2c?prepared_s2c->conjugate.groups.size():0)<<",\"s2c_retained_spectrum_bytes\":"<<s2c_retained_bytes<<",\"s2c_public_keys_prepared\":"<<(have_prepared_s2c?s2c_input_cache->prepared_keys+s2c_output_cache->prepared_keys:0)<<",\"s2c_seconds\":"<<s2c_seconds<<",\"blind_rotation_and_series_seconds\":"<<blind_seconds<<",\"functional_evaluation_seconds\":"<<s2c_seconds+blind_seconds<<",\"output_log_q\":93,\"output_log_delta\":44,\"post_product_log_q\":49,\"max_error\":"<<error<<",\"post_product_max_error\":"<<product_error<<",\"product_tolerance\":"<<product_tolerance<<",\"accuracy_passed\":"<<(passed?"true":"false")<<",\"security_validated\":false}";}
 std::cout<<artifact<<" error="<<error<<" square error="<<product_error<<" seconds="<<s2c_seconds+blind_seconds<<std::endl;if(!passed)throw std::runtime_error("functional accuracy gate failed");
 ++trial_index;
}
} catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
