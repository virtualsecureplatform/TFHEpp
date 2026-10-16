// Deterministic noncryptographic coefficient rows; no keys generated or encrypted.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_gadget_batch_rns32.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <iomanip>
#include <iostream>
using namespace TFHEpp;
template<class P,std::uint32_t Q,std::uint32_t Bits>struct Rows {
 static constexpr std::size_t rows=(Q+Bits-1)/Bits;
 std::array<TRLWE<P>,rows> ciphertexts;
};
template<class P,std::uint32_t Q,std::uint32_t Bits>auto fixtureRows(){
 auto rows=std::make_unique<Rows<P,Q,Bits>>();std::mt19937_64 random(20261027+Q+Bits);
 for(auto&r:rows->ciphertexts)for(auto&c:r)for(auto&v:c){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,Q>(v);}return rows;
}
template<class P,std::uint32_t Q,std::uint32_t Bits,std::uint32_t Lift>auto inputs(std::size_t size,bool edges){
 std::vector<Polynomial<P>> result(size);std::mt19937_64 random(20261028+Q+Bits);
 for(std::size_t j=0;j<size;++j)for(auto&v:result[j]){
  for(auto&w:v.limb)w=random();
  v=ckks_detail::reduceToLevel<P,Q>(v<<Lift);
  if(edges){
   if(j%7==1)v=typename P::T{1}<<Q;
   if(j%7==2)v=typename P::T{}-typename P::T{1};
   if(j%7==3)v=(typename P::T{1}<<(Q-1))+typename P::T{1};
   if(j%7==4)v=typename P::T{random()&((std::uint64_t{1}<<Bits)-1)};
   if(j%7==5)v=(typename P::T{1}<<Q)+typename P::T{1};
  }
 }
 return result;
}
template<class P,std::uint32_t Q,std::uint32_t Bits>void oracle(TRLWE<P>&out,const Polynomial<P>&input,const Rows<P,Q,Bits>&key){
 using I=boost::multiprecision::cpp_int;
 const auto integer=[](const auto&v){I r=0;for(std::size_t w=P::T::limbs;w-->0;){r<<=64;r+=v.limb[w];}return r;};
 const I modulus=I{1}<<Q,mask=modulus-1,half=modulus/2,base=I{1}<<Bits,dhalf=base/2;
 std::array<std::vector<I>,2> sums;for(auto&c:sums)c.resize(P::n);
 for(std::size_t i=0;i<P::n;++i){
  I value=integer(input[i])&mask;if(value>=half)value-=modulus;
  for(std::size_t r=0;r<key.rows;++r){
   I digit=(value+dhalf)%base;if(digit<0)digit+=base;digit-=dhalf;value=(value-digit)/base;
   const auto row=key.rows-1-r;
   for(std::size_t c=0;c<2;++c)for(std::size_t j=0;j<P::n;++j){const I product=digit*integer(key.ciphertexts[row][c][j]);sums[c][(i+j)%P::n]+=(i+j<P::n?product:-product);}
  }
 }
 for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i){
  I value=sums[c][i]&mask;out[c][i]={};
  for(auto&w:out[c][i].limb){w=static_cast<std::uint64_t>(value&((I{1}<<64)-1));value>>=64;}
 }
}
template<std::uint32_t Nbit,std::uint32_t Q,std::uint32_t Bits,std::uint32_t Chunk,std::uint32_t Lift>void check(){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;
 auto rows=fixtureRows<P,Q,Bits>();CKKSPreparedGadgetRNS32<P,Q,Bits,Chunk> candidate(*rows);
 CKKSPreparedGadgetRNS32<P,Q,Bits,Chunk> original(*rows);
 auto in=inputs<P,Q,Bits,Lift>(7,true);std::vector<TRLWE<P>> expected(in.size()),actual(in.size());
 for(std::size_t j=0;j<in.size();++j){
  original.apply(expected[j],in[j]);
  if constexpr(Nbit<=4){TRLWE<P> independent;oracle<P,Q,Bits>(independent,in[j],*rows);if(independent!=expected[j])throw std::runtime_error("single gadget differs from independent integer oracle");}
 }
 const auto compare=[&](std::size_t columns,std::size_t window){
  for(auto&r:actual)for(auto&c:r)c.fill(typename P::T{}-typename P::T{1});
  if(columns==1)CKKSApplyGadgetBatchRNS32<1>(actual,in,candidate,window);
  else if(columns==4)CKKSApplyGadgetBatchRNS32<4>(actual,in,candidate,window);
  else CKKSApplyGadgetBatchRNS32<8>(actual,in,candidate,window);
  if(actual!=expected)throw std::runtime_error("batched gadget changed exact coefficients");
 };
 for(auto columns:{1u,4u,8u})for(auto window:{3u,8u})compare(columns,window);
#pragma omp parallel num_threads(2)
 {
#pragma omp single
  CKKSApplyGadgetBatchRNS32<4>(actual,in,candidate,3);
 }
 if(actual!=expected)throw std::runtime_error("batch failed under caller's nesting policy");
 auto zero=in;for(auto&poly:zero)poly.fill(typename P::T{1}<<Q);
 CKKSApplyGadgetBatchRNS32<4>(actual,zero,candidate,8);
 for(const auto&r:actual)for(const auto&c:r)for(const auto&v:c)
  if(v!=typename P::T{})throw std::runtime_error("all-zero batch retained stale output");
 std::vector<TRLWE<P>> empty;std::vector<Polynomial<P>> no_input;CKKSApplyGadgetBatchRNS32<4>(empty,no_input,candidate);
 bool rejected=false;try{CKKSApplyGadgetBatchRNS32<4>(empty,in,candidate);}catch(const std::invalid_argument&){rejected=true;}if(!rejected)throw std::runtime_error("batch shape mismatch accepted");
 rejected=false;try{CKKSApplyGadgetBatchRNS32<4>(actual,in,candidate,0);}catch(const std::invalid_argument&){rejected=true;}if(!rejected)throw std::runtime_error("zero batch window accepted");
 std::cout<<"check N="<<P::n<<" q="<<Q<<" bits="<<Bits<<" chunk="<<Chunk<<" mixed active rows and zero clearing identical\n";
}
template<std::uint32_t Q,std::uint32_t Bits,std::uint32_t Chunk,std::uint32_t Lift>void benchmark(std::size_t size){
 using P=narrow_benchmark::BenchmarkParam<14,7>;
 auto rows=fixtureRows<P,Q,Bits>();CKKSPreparedGadgetRNS32<P,Q,Bits,Chunk> key(*rows);
 auto in=inputs<P,Q,Bits,Lift>(size,false);std::vector<TRLWE<P>> reference(size),out(size);
 [[maybe_unused]] const auto threads=std::min<int>(size,omp_get_max_threads());
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
 for(std::size_t j=0;j<size;++j)key.apply(reference[j],in[j]);
 const auto run=[&](int mode){
  const auto start=std::chrono::steady_clock::now();
  if(!mode){
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
   for(std::size_t j=0;j<size;++j)key.apply(out[j],in[j]);
  }else if(mode==1)CKKSApplyGadgetBatchRNS32<1>(out,in,key,32);
  else if(mode==4)CKKSApplyGadgetBatchRNS32<4>(out,in,key,32);
  else if(mode==8)CKKSApplyGadgetBatchRNS32<8>(out,in,key,32);
  else if(mode==48)CKKSApplyGadgetBatchRNS32<4>(out,in,key,8);
  else CKKSApplyGadgetBatchRNS32<4>(out,in,key,16);
  const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if(out!=reference)throw std::runtime_error("timed batch gadget mismatch");
  return seconds;
 };
 std::vector<int> modes{0,1,4,8};if(size==32){modes.push_back(48);modes.push_back(416);}
 for(auto mode:modes)run(mode);
 auto order=modes;std::reverse(modes.begin(),modes.end());order.insert(order.end(),modes.begin(),modes.end());
 for(auto mode:order)for(int trial=0;trial<3;++trial)std::cout<<"q="<<Q<<" bits="<<Bits<<" chunk="<<Chunk<<" lift="<<Lift<<" batch="<<size<<" mode="<<mode<<" trial="<<trial<<" seconds="<<run(mode)<<" identical\n";
}
int main(int argc,char**argv){
 std::cout<<std::setprecision(17);
 if(argc>1&&std::string(argv[1])=="--check"){
  check<4,93,4,16,8>();check<4,389,40,80,48>();check<4,345,40,115,48>();check<4,301,40,80,48>();check<4,435,44,112,37>();
  check<4,433,40,145,48>();check<4,435,44,145,37>();check<4,257,40,145,48>();
  check<9,433,40,145,48>();check<9,435,44,145,270>();
  check<14,389,40,80,48>();check<14,345,40,115,48>();check<14,301,40,80,48>();check<14,257,40,145,48>();check<14,435,44,112,37>();
  check<14,433,40,145,48>();check<14,435,44,145,37>();return 0;
 }
 for(std::size_t size:{4,8,32}){benchmark<433,40,145,48>(size);benchmark<435,44,145,37>(size);benchmark<435,44,145,270>(size);}
}
