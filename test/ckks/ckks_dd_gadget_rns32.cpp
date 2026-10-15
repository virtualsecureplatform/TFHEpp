#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_gadget_rns32.hpp>
#include <ckks/prepared_dd_ntt.hpp>
#include <iostream>
using namespace TFHEpp;
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
 input->fill({});direct.apply(*actual,*input);for(auto&c:*actual)for(auto&v:c)if(v!=typename P::T{})throw std::runtime_error("zero input retained stale output");
}
int main(int argc,char**){check<4,385>();check<9,385>();check<9,341>();check<9,297>();check<9,253>();check<9,209>();if(argc>1){check<14,385>();check<14,341>();check<14,297>();check<14,253>();check<14,209>();}}
