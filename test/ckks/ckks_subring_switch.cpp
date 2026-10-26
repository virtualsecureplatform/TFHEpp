#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/subring_switch.hpp>
#include <sstream>
#include <random>

namespace {
using P=narrow_benchmark::BenchmarkParam<6,7>;
struct LoudP:P { static constexpr double ckks_min_noise_stddev=512.; };
struct L:TFHEpp::lvl1param { static constexpr std::uint32_t n=8,nbit=3,k=1; };
struct FullTarget:TFHEpp::lvl1param { static constexpr std::uint32_t n=1024,nbit=10,k=1; };
using K=TFHEpp::CKKSToTFHESubringSwitchKey32<P,L>;
using Input=TFHEpp::CKKSCiphertext<P,32,29>;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F> void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}require(rejected,"invalid input accepted");}
template<std::size_t N> auto multiply(const std::array<std::uint32_t,N>&a,const std::array<std::uint32_t,N>&b){
 std::array<std::uint32_t,N> out{};
 for(std::size_t i=0;i<N;++i)for(std::size_t j=0;j<N;++j){
  const auto product=std::uint32_t(std::uint64_t(a[i])*b[j]);
  if(i+j<N)out[i+j]+=product;else out[i+j-N]-=product;
 }return out;
}
std::uint32_t phase(const TFHEpp::TLWE<L>&ct,const TFHEpp::Key<L>&secret){
 auto out=ct[L::n];for(std::size_t i=0;i<L::n;++i)out-=ct[i]*secret[i];return out;
}
void arithmetic(){
 using namespace TFHEpp;
 auto key=std::make_unique<K>();std::mt19937_64 random(20261025);
 for(auto&r:key->ciphertexts)for(auto&c:r)for(auto&v:c)v=random();
 PreparedCKKSToTFHESubringSwitch32<P,L> prepared(*key);
 static_assert(decltype(prepared)::prime_count==2);
 require(prepared.byteSize()==2*16*2*P::n*4,"unexpected prepared cache size");
 for(int trial=0;trial<4;++trial){
  Input input;for(auto&c:input.ct)for(auto&v:c)v=P::T{random()};
  for(std::size_t i=0;i<P::n;++i){
   if(trial==1)input.ct[0][i]=P::T{};
   if(trial==2)input.ct[0][i]=P::T{std::uint32_t(1u<<(i%32))};
   if(trial==3)input.ct[0][i]=P::T{std::array<std::uint32_t,8>{0,1,2,3,0x7fffffff,0x80000000,0xfffffffe,0xffffffff}[i%8]};
  }
  std::array<std::array<std::uint32_t,P::n>,2> expected{};
  for(std::size_t r=0;r<16;++r){
   std::array<std::uint32_t,P::n> digits{};
   for(std::size_t i=0;i<P::n;++i){
    const auto a=std::uint32_t(input.ct[0][i].limb[0])+std::uint32_t{0xaaaaaaaa};
    digits[i]=std::uint32_t(int((a>>(30-2*r))&3)-2);
   }
   for(std::size_t c=0;c<2;++c){auto product=multiply(digits,key->ciphertexts[r][c]);for(std::size_t i=0;i<P::n;++i)expected[c][i]-=product[i];}
  }
  for(std::size_t i=0;i<P::n;++i)expected[1][i]+=std::uint32_t(input.ct[1][i].limb[0]);
  std::vector<TRLWE<L>> split;prepared.switchAndSplit(split,input);
  require(split.size()==K::stride,"split count");
  for(std::size_t r=0;r<K::stride;++r)for(std::size_t c=0;c<2;++c)for(std::size_t j=0;j<L::n;++j)
   require(split[r][c][j]==expected[c][r+K::stride*j],"NTT/split differs from integer oracle");
  for(auto count:{std::size_t{1},std::size_t{7},std::size_t{33},std::size_t{P::n}}){
   std::vector<TLWE<L>> out;prepared.apply(out,input,count);
   for(std::size_t i=0;i<count;++i){
    const auto r=i%K::stride,j=i/K::stride;
    require(out[i][L::n]==expected[1][i],"extracted body order");
    for(std::size_t k=0;k<L::n;++k){const auto value=k<=j?expected[0][r+K::stride*(j-k)]:std::uint32_t{0}-expected[0][r+K::stride*(L::n+j-k)];require(out[i][k]==value,"extracted mask order/wrap");}
   }
  }
  std::vector<TLWE<L>> out;rejects([&]{prepared.apply(out,input,0);});rejects([&]{prepared.apply(out,input,P::n+1);});
 }
 // Concurrent calls use immutable caches and independent scratch.
 Input zero{};std::array<std::vector<TLWE<L>>,2> outputs;
#pragma omp parallel for
 for(int i=0;i<2;++i)prepared.apply(outputs[i],zero,P::n);
 require(outputs[0]==outputs[1],"concurrent cache reuse differs");
}
void encrypted(){
 using namespace TFHEpp;Key<P> source{};Key<L> target{};
 for(std::size_t i=0;i<P::n;++i)source[i]=P::T::from_signed_i64(int(i%3)-1);
 for(std::size_t i=0;i<L::n;++i)target[i]=std::uint32_t(int(i%3)-1);
 auto key=std::make_unique<K>();CKKSToTFHESubringSwitchKeyGen32(*key,source,target);
 auto loud=std::make_unique<CKKSToTFHESubringSwitchKey32<LoudP,L>>();Key<LoudP> loud_source{};
 rejects([&]{CKKSToTFHESubringSwitchKeyGen32(*loud,loud_source,target);});
 PreparedCKKSToTFHESubringSwitch32<P,L> prepared(*key);
 std::array<std::uint32_t,P::n> a{},s{};std::mt19937 random(8124);
 for(std::size_t i=0;i<P::n;++i){a[i]=random();s[i]=source[i].limb[0];}
 auto as=multiply(a,s);Input input;
 for(std::size_t i=0;i<P::n;++i){input.ct[0][i]=P::T{a[i]};input.ct[1][i]=P::T{as[i]+std::uint32_t(i%2?1u<<29:0u-(1u<<29))};}
 std::vector<TLWE<L>> outputs;prepared.apply(outputs,input,P::n);
 for(std::size_t i=0;i<P::n;++i){const auto expected=std::uint32_t(i%2?1u<<29:0u-(1u<<29));auto error=std::bit_cast<std::int32_t>(phase(outputs[i],target)-expected);require(std::abs(double(error))<0.02575*4294967296.,"encrypted subring phase error");require(tlweSymDecrypt<L>(outputs[i],target)==bool(i%2),"Boolean extraction");}
 const std::string aorigin(64,'a'),borigin(64,'b');std::stringstream stream;
 CKKSWriteSubringBridgeKey<P,L>(stream,*key,1,aorigin,borigin);const auto bytes=stream.str();auto loaded=std::make_unique<K>();
 CKKSReadSubringBridgeKey<P,L>(stream,*loaded,1,aorigin,borigin);require(loaded->ciphertexts==key->ciphertexts,"archive roundtrip");
 auto read=[&](std::string data,std::uint64_t kind,std::string origin){std::stringstream file(data);CKKSReadSubringBridgeKey<P,L>(file,*loaded,kind,origin,borigin);};
 rejects([&]{read(bytes,2,aorigin);});rejects([&]{read(bytes,1,std::string(64,'c'));});rejects([&]{read(bytes.substr(0,bytes.size()-1),1,aorigin);});rejects([&]{read(bytes+"x",1,aorigin);});
 auto altered=[&](std::size_t field,std::uint64_t value){
  std::stringstream file;{cereal::PortableBinaryOutputArchive ar(file);auto profile=CKKSSubringBridgeProfile<P,L>(1);profile[field]=value;std::array<char,64>a{},b{};a.fill('a');b.fill('b');ar(profile,a,b,*key);}return file.str();
 };
 rejects([&]{read(altered(10,std::bit_cast<std::uint64_t>(128.)),1,aorigin);});
 rejects([&]{read(altered(3,L::n*2),1,aorigin);});
 rejects([&]{read(altered(0,0),1,aorigin);});
 rejects([&]{read(altered(0,0x53554252494e4731ULL),1,aorigin);});
 auto bad=source;bad[0]=P::T{2};rejects([&]{CKKSToTFHESubringSwitchKeyGen32(*key,bad,target);});
 auto direct=std::make_unique<CKKSToTFHEKeySwitchKey<P,L,2,16>>();
 rejects([&]{CKKSToTFHEKeyGen<P,L,2,16>(*direct,source,target,0.);});
 CKKSToTFHEKeyGen<P,L,2,16>(*direct,source,target,std::ldexp(256.,-32));
 TLWE<L> out;CKKSCoeffToTFHE<P,L,32,29,2,16>(out,input,0,*direct);require(!tlweSymDecrypt<L>(out,target),"explicit-noise direct switch");
}
void full_batch(){
 using namespace TFHEpp;using Source=narrow_benchmark::BenchmarkParam<14,7>;
 using Target=FullTarget;
 auto source=std::make_unique<Key<Source>>();Key<Target> target{};
 for(std::size_t i=0;i<Source::n;++i)(*source)[i]=Source::T::from_signed_i64(int(i%3)-1);
 for(std::size_t i=0;i<Target::n;++i)target[i]=std::uint32_t(int(i%3)-1);
 auto key=std::make_unique<CKKSToTFHESubringSwitchKey32<Source,Target>>();
 CKKSToTFHESubringSwitchKeyGen32(*key,*source,target);
 PreparedCKKSToTFHESubringSwitch32<Source,Target> prepared(*key);key.reset();
 auto message=std::make_unique<Polynomial<Source>>();
 for(std::size_t i=0;i<Source::n;++i)(*message)[i]=Source::T{std::uint32_t(i%2?1u<<29:0u-(1u<<29))};
 auto input_row=std::make_unique<CKKSSeededKeySwitchRow<Source,32>>();
 ckks_scheme_switch_detail::exactEncryptPolynomialAtLevel<Source,32>(*input_row,*message,*source,{0,0});
 auto input=std::make_unique<CKKSCiphertext<Source,32,29>>();
 ckks_detail::fillSeededUniformPolynomialAtLevel<Source,32>(input->ct[0],input_row->mask_seeds[0]);
 input->ct[1]=input_row->body;
 std::vector<TLWE<Target>> outputs;prepared.apply(outputs,*input,8192);
 for(std::size_t i=0;i<outputs.size();++i){
  const auto actual=tlweSymPhase<Target>(outputs[i],target),expected=std::uint32_t((*message)[i].limb[0]);
  require(std::abs(std::ldexp(double(std::bit_cast<std::int32_t>(actual-expected)),-32))<.025746934137548232,"full-batch phase error");
  require(tlweSymDecrypt<Target>(outputs[i],target)==bool(i%2),"full-batch bit mismatch");
 }
}
}
int main()try{arithmetic();encrypted();full_batch();std::cout<<"Subring switch integer, encryption, split, archive, cache and 8192-bit tests passed\n";}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
