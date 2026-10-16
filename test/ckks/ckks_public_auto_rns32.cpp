#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_auto_ntt.hpp>
#include <ckks/native_auto_file.hpp>
#include <sstream>
#include <iostream>
#include <random>
using namespace TFHEpp;
template<class P>struct DispatchSchedule{using Param=P;static constexpr std::uint32_t factor_log_q=385,output_log_q=165;};
template<std::uint32_t Nbit,std::uint32_t Bits=20,std::uint32_t ChunkBits=80>void check(){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;using T=typename P::T;
 constexpr std::uint32_t Q=435;
 using Rows=PublicAutoGadgetRows<P,Q,Bits>;
 auto rows=std::make_unique<Rows>();rows->exponent=2*P::n-1;
 std::mt19937_64 random(20261016);
 for(auto &r:rows->ciphertexts)for(auto &c:r)for(auto &v:c){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,Q>(v);}
 if constexpr(Bits==20){
  auto grouped=makePublicAutoGadgetRowsWithPrimaryBase<40>(*rows);
  if(grouped->exponent!=rows->exponent)throw std::runtime_error("public gadget regrouping changed exponent");
  for(std::size_t r=0;r<grouped->rows;++r){
   const auto power=(grouped->rows-r-1)*40;
   if(grouped->ciphertexts[r]!=rows->ciphertexts[rows->rows-1-power/20])throw std::runtime_error("public gadget regrouping changed encrypted row bytes");
  }
 }
 PreparedPublicAutoNTT<P,Q,Bits,20> reference(*rows);
 PreparedPublicAutoNTT<P,Q,Bits,20,true,ChunkBits> direct(*rows);
 auto cache_map=std::make_shared<std::map<uint32_t,std::shared_ptr<decltype(direct)>>>();
 cache_map->emplace(rows->exponent,std::make_shared<decltype(direct)>(direct));
 auto dispatch=narrow_detail::makePublicAutoEvaluator<DispatchSchedule<P>,398,121>(cache_map);
 auto input=std::make_unique<TRLWE<P>>();
 for(auto &c:*input)for(auto &v:c){for(auto&w:v.limb)w=random();}
 auto expected=std::make_unique<TRLWE<P>>(),actual=std::make_unique<TRLWE<P>>();
 auto saved=std::make_unique<TRLWE<P>>(*input);
 const auto compare=[&]<std::uint32_t DataQ>(){
  for(int trial=0;trial<5;++trial){
   *input=*saved;
   if(trial==3)for(auto&c:*input)c.fill({});
   if(trial==4)for(auto&c:*input)c.fill((T{1}<<(DataQ-1))+T{1});
   double old=0,now=0;
   const auto run_old=[&]{auto t=std::chrono::steady_clock::now();reference.template apply<DataQ>(*expected,*input);old=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
   const auto run_new=[&]{auto t=std::chrono::steady_clock::now();direct.template apply<DataQ>(*actual,*input);now=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
   if(trial%2){run_new();run_old();}else{run_old();run_new();}
   if(*actual!=*expected)throw std::runtime_error("wide-chunk RNS automorphism differs from digit NTT");
   std::cout<<"N="<<P::n<<" bits="<<Bits<<" data_q="<<DataQ<<" trial="<<trial<<" old="<<old<<" rns="<<now<<" old_bytes="<<reference.byteSize()<<" rns_bytes="<<direct.byteSize()<<" identical\n";
  }
  *input=*saved;
 };
 compare.template operator()<385>();compare.template operator()<398>();compare.template operator()<121>();
 if constexpr(Bits==44){
  const auto compare_batch=[&]<uint32_t DataQ>(){
   std::vector<TRLWE<P>> in(3,*saved),want(3),got(3);
   for(auto&c:in[1])c.fill(T{1}<<DataQ);
   for(auto&c:in[2])c.fill((T{1}<<(DataQ-1))+T{1});
   for(std::size_t j=0;j<in.size();++j)reference.template apply<DataQ>(want[j],in[j]);
   direct.template applyBatch<DataQ>(got,in);
   if(got!=want)throw std::runtime_error("batched native auto differs from digit NTT");
   auto alias=in;direct.template applyBatch<DataQ>(alias,alias);
   if(alias!=want)throw std::runtime_error("in-place native auto batch changed coefficients");
   std::vector<TRLWE<P>> empty;direct.template applyBatch<DataQ>(empty,empty);
   bool rejected=false;try{direct.template applyBatch<DataQ>(empty,in);}catch(const std::invalid_argument&){rejected=true;}
   if(!rejected)throw std::runtime_error("invalid native auto batch shape accepted");
  };
  compare_batch.template operator()<398>();compare_batch.template operator()<385>();compare_batch.template operator()<165>();compare_batch.template operator()<121>();
  auto batch_dispatch=narrow_detail::makePublicAutoBatchEvaluator<DispatchSchedule<P>,398,121>(cache_map);
  std::vector<TRLWE<P>> in(3,*saved),got(3);
  batch_dispatch(got,in,165,rows->exponent);reference.template apply<165>(*expected,*saved);
  for(const auto&v:got)if(v!=*expected)throw std::runtime_error("native auto batch dispatch changed level");
  bool rejected=false;try{batch_dispatch(got,in,122,rows->exponent);}catch(const std::invalid_argument&){rejected=true;}
  if(!rejected)throw std::runtime_error("unregistered native auto batch level accepted");
  rejected=false;try{batch_dispatch(got,in,165,1);}catch(const std::invalid_argument&){rejected=true;}
  if(!rejected)throw std::runtime_error("unregistered native auto batch exponent accepted");
  std::cout<<"native auto batches match digit NTT; alias, empty, shape and dispatch guards pass\n";
 }
 reference.template apply<165>(*expected,*input);
 dispatch(*actual,*input,165,rows->exponent);
 if(*actual!=*expected)throw std::runtime_error("shared root-level cache differs");
 reference.template apply<121>(*expected,*input);
 dispatch(*actual,*input,121,rows->exponent);
 if(*actual!=*expected)throw std::runtime_error("shared series-level cache differs");
 reference.template apply<156>(*expected,*input);
 auto intermediate_dispatch=narrow_detail::makePublicAutoEvaluator<DispatchSchedule<P>,398,121,156>(cache_map);
 intermediate_dispatch(*actual,*input,156,rows->exponent);
 if(*actual!=*expected)throw std::runtime_error("shared intermediate-level cache differs");
 bool missing_exponent=false;
 try{dispatch(*actual,*input,121,1);}catch(const std::invalid_argument&){missing_exponent=true;}
 if(!missing_exponent)throw std::runtime_error("unregistered public exponent accepted");
 bool rejected=false;
 try{narrow_detail::applyPublicAutoAtLevel<DispatchSchedule<P>,398>(*actual,*input,121,direct);}
 catch(const std::invalid_argument&){rejected=true;}
 if(!rejected)throw std::runtime_error("unregistered extra level accepted");
 rejected=false;
 try{narrow_detail::applyPublicAutoAtLevel<DispatchSchedule<P>,398,121>(*actual,*input,122,direct);}
 catch(const std::invalid_argument&){rejected=true;}
 if(!rejected)throw std::runtime_error("unexpected extra level accepted");
 if constexpr(Nbit<=4){
  // Compile and construct the existing default wide-primary backend too.
  PublicAutoGadgetRows<P,Q,64> legacy_rows;
  PreparedPublicAutoNTT<P,Q,64> legacy(legacy_rows);
  CKKSPreparedGadgetRNS32<P,Q,Bits,ChunkBits> gadget(*rows);
  auto digits=std::make_unique<std::array<Polynomial<P>,Rows::rows>>();
  ckks_detail::activeBaseDecomposePolynomialRows<P,Q,Bits,Rows::rows>(*digits,(*input)[0]);
  gadget.apply(*actual,(*input)[0]);for(auto&c:*expected)c.fill({});
  for(std::size_t r=0;r<Rows::rows;++r)for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)for(std::size_t j=0;j<P::n;++j){
   const auto v=(*digits)[r][i]*rows->ciphertexts[r][c][j];const auto k=(i+j)%P::n;
   (*expected)[c][k]=ckks_detail::reduceToLevel<P,Q>((*expected)[c][k]+(i+j<P::n?v:-v));
  }
  if(*actual!=*expected)throw std::runtime_error("full RNS gadget differs from schoolbook oracle");
 }
}
void checkNativeRows()
{
 using P=narrow_benchmark::BenchmarkParam<4,7>;constexpr uint32_t Q=435,Bits=44,exponent=2*P::n-1;
 Key<P> fixed{};for(std::size_t i=0;i<P::n;++i)fixed[i]=P::T::from_signed_i64(int(i%3)-1);
 auto rows=std::make_unique<PublicAutoGadgetRows<P,Q,Bits>>();CKKSGenerateNativeAutoGadgetRows(*rows,exponent,fixed);
 for(std::size_t row=0;row<rows->rows;++row){
  auto phase=rows->ciphertexts[row][1];Polynomial<P> expected{};
  for(std::size_t i=0;i<P::n;++i)for(std::size_t j=0;j<P::n;++j){
   const auto value=rows->ciphertexts[row][0][i]*fixed[j];phase[(i+j)%P::n]-=i+j<P::n?value:-value;
  }
  for(std::size_t i=0;i<P::n;++i){
   const auto position=i*exponent%(2*P::n);const auto value=fixed[i]<<((rows->rows-row-1)*Bits);
   expected[position%P::n]=position<P::n?value:-value;
  }
  for(std::size_t i=0;i<P::n;++i){
   const auto error=ckks_detail::reduceToLevel<P,Q>(phase[i]-expected[i]);
   const auto small=std::bit_cast<int64_t>(error.limb[0]);
   if(error!=ckks_detail::reduceToLevel<P,Q>(P::T::from_signed_i64(small))||small>128||small< -128)throw std::runtime_error("native gadget phase differs from independent schoolbook oracle");
  }
 }
 const std::string origin(64,'a');std::stringstream saved;
 CKKSWriteNativeAutoRows(saved,*rows,origin,3.2);const auto bytes=saved.str();
 auto loaded=CKKSReadNativeAutoRows<P,Q,Bits>(saved,exponent,origin,3.2);
 if(loaded->exponent!=rows->exponent||loaded->ciphertexts!=rows->ciphertexts)throw std::runtime_error("native archive changed encrypted rows");
 const auto reject=[&](std::string data,uint32_t target,const std::string&source,double sigma){
  bool rejected=false;try{std::stringstream input(data);CKKSReadNativeAutoRows<P,Q,Bits>(input,target,source,sigma);}catch(const std::exception&){rejected=true;}
  if(!rejected)throw std::runtime_error("invalid native archive accepted");
 };
 reject(bytes,1,origin,3.2);reject(bytes,exponent,std::string(64,'b'),3.2);reject(bytes,exponent,origin,6.4);
 reject(bytes.substr(0,bytes.size()-1),exponent,origin,3.2);reject(bytes+"x",exponent,origin,3.2);
 std::cout<<"native key phases match schoolbook; archive profile/origin/truncation guards pass\n";
}

int main(int argc,char**argv){if(argc>1&&std::string(argv[1])=="--native-batch-only"){check<4,44,112>();check<4,44,145>();check<9,44,145>();check<14,44,112>();check<14,44,145>();return 0;}if(argc>1&&std::string(argv[1])=="--native-chunks-only"){check<4,44,145>();check<9,44,145>();check<14,44,145>();return 0;}check<4>();check<9>();check<4,40>();check<9,40>();check<4,44,112>();check<9,44,112>();checkNativeRows();if(argc>1){check<14>();check<14,40>();}}
