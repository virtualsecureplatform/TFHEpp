#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_gadget_rns32.hpp>
#include <ckks/prepared_dd_ntt.hpp>
#include <ckks/key_projection.hpp>
#include <iostream>
using namespace TFHEpp;
template<std::uint32_t Q, std::uint32_t Bits> void checkSignedDigits() {
 using P=narrow_benchmark::BenchmarkParam<11,7>;
 constexpr auto rows=(Q+Bits-1)/Bits;
 if constexpr(rows*Bits<=P::T::digits) {
  auto input=std::make_unique<Polynomial<P>>();
  auto wide=std::make_unique<std::array<Polynomial<P>,rows>>();
  auto compact=std::make_unique<std::array<std::array<std::int64_t,P::n>,rows>>();
  std::mt19937_64 random(0xabc000+Q*100+Bits);
  for(int trial=0;trial<3;++trial) {
   for(auto&v:*input) { for(auto&w:v.limb)w=random();
    if(trial==1)v=ckks_detail::reduceToLevel<P,Q>(v);
   }
   // Zero, signed endpoints, modulus midpoint, and carry boundaries.
   (*input)[0]={};(*input)[1]=typename P::T{1};(*input)[2]=-typename P::T{1};
   (*input)[3]=typename P::T{1}<<(Q-1);
   for(std::size_t r=0;r<rows;++r)for(int d=-1;d<=1;++d)
    (*input)[4+r*3+d+1]=(typename P::T{1}<<(r*Bits+Bits-1))+P::T::from_signed_i64(d);
   ckks_detail::activeBaseDecomposePolynomialRows<P,Q,Bits,rows>(*wide,*input);
   ckks_detail::activeBaseDecomposeSignedRows<P,Q,Bits,rows>(*compact,*input);
   for(std::size_t i=0;i<P::n;++i) {
    typename P::T reconstructed{};
    for(std::size_t r=0;r<rows;++r) {
     if((*compact)[r][i]!=std::bit_cast<std::int64_t>((*wide)[r][i].limb[0]))
      throw std::runtime_error("compact signed digit differs from wide reference");
     reconstructed+=P::T::from_signed_i64((*compact)[r][i])<<((rows-r-1)*Bits);
    }
    if(ckks_detail::reduceToLevel<P,Q>(reconstructed)!=ckks_detail::reduceToLevel<P,Q>((*input)[i]))
     throw std::runtime_error("signed digits do not reconstruct input");
   }
  }
  std::cout<<"signed digits q="<<Q<<" bits="<<Bits<<" identical and reconstruct input\n";
 }
}
template<std::uint32_t Q> void checkSignedBases() {
 checkSignedDigits<Q,1>();checkSignedDigits<Q,20>();checkSignedDigits<Q,40>();
 checkSignedDigits<Q,44>();checkSignedDigits<Q,64>();
}
template<std::uint32_t Nbit,std::uint32_t DataQ>void check(){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;
 constexpr auto Q=DataQ+48;using Compact=CKKSCompactDDRelinKey<P,Q,20,20>;
 auto key=std::make_unique<Compact>();auto ordinary=std::make_unique<TRLWE<P>>();
 auto digits=std::make_unique<std::array<TRLWE<P>,Compact::bbar_rows>>();
 std::mt19937_64 random(20261017+Q);
 for(std::size_t r=0;r<Compact::primary_rows;++r){
  for(auto&c:*ordinary)for(auto&v:c){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,Q>(v);}
  ckks_detail::activeBaseDecomposeTRLWERows<P,Q,20,Compact::bbar_rows>(*digits,*ordinary);
  for(std::size_t d=0;d<Compact::bbar_rows;++d)for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)key->data[r][d][c][i]=Compact::pack((*digits)[d][c][i]);
 }
 CKKSPreparedDDRelinKeyNTT<P,Q,20,20> old(*key);CKKSPreparedGadgetRNS32<P,Q,20> direct(*key);
 auto input=std::make_unique<Polynomial<P>>();auto expected=std::make_unique<TRLWE<P>>(),actual=std::make_unique<TRLWE<P>>();
 for(auto&v:*input){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,Q>(v<<48);}
 for(int trial=0;trial<4;++trial){
  double old_time=0,new_time=0;
  const auto run_old=[&]{auto t=std::chrono::steady_clock::now();old.apply(*expected,*input);old_time=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
  const auto run_new=[&]{auto t=std::chrono::steady_clock::now();direct.apply(*actual,*input);new_time=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
  if(trial%2){run_new();run_old();}else{run_old();run_new();}
  if(*actual!=*expected)throw std::runtime_error("wide DD gadget differs from digit NTT");
  std::cout<<"N="<<P::n<<" data_q="<<DataQ<<" trial="<<trial<<" old="<<old_time<<" rns="<<new_time<<" old_bytes="<<old.byteSize()<<" rns_bytes="<<direct.byteSize()<<" identical\n";
 }
 auto grouped=makeCKKSCompactDDRelinKeyWithPrimaryBase<40>(*key);
 using Grouped=CKKSCompactDDRelinKey<P,Q,40,20>;
 for(std::size_t r=0;r<Grouped::primary_rows;++r){
  const auto power=(Grouped::primary_rows-r-1)*40;
  if(grouped->data[r]!=key->data[Compact::primary_rows-1-power/20])throw std::runtime_error("primary regrouping changed encrypted row bytes");
 }
 CKKSPreparedDDRelinKeyNTT<P,Q,40,20> grouped_reference(*grouped);
 CKKSPreparedGadgetRNS32<P,Q,40> grouped_direct(*grouped);
 grouped_reference.apply(*expected,*input);grouped_direct.apply(*actual,*input);
 if(*actual!=*expected)throw std::runtime_error("forty-bit grouped gadget differs from digit NTT");
 input->fill({});direct.apply(*actual,*input);for(auto&c:*actual)for(auto&v:c)if(v!=typename P::T{})throw std::runtime_error("zero input retained stale output");
}
int main(int argc,char**){checkSignedBases<85>();checkSignedBases<209>();checkSignedBases<257>();checkSignedBases<385>();checkSignedBases<433>();checkSignedBases<435>();checkSignedBases<448>();check<4,385>();check<9,385>();check<9,341>();check<9,297>();check<9,253>();check<9,209>();if(argc>1){check<14,385>();check<14,341>();check<14,297>();check<14,253>();check<14,209>();}}
