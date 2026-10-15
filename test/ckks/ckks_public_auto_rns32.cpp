#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/public_auto_ntt.hpp>
#include <iostream>
#include <random>
using namespace TFHEpp;
template<std::uint32_t Nbit>void check(){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;using T=typename P::T;
 constexpr std::uint32_t Q=435,Bits=20;
 using Rows=PublicAutoGadgetRows<P,Q,Bits>;
 auto rows=std::make_unique<Rows>();rows->exponent=2*P::n-1;
 std::mt19937_64 random(20261016);
 for(auto &r:rows->ciphertexts)for(auto &c:r)for(auto &v:c){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,Q>(v);}
 PreparedPublicAutoNTT<P,Q,Bits,20> reference(*rows);
 PreparedPublicAutoNTT<P,Q,Bits,20,true> direct(*rows);
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
   std::cout<<"N="<<P::n<<" data_q="<<DataQ<<" trial="<<trial<<" old="<<old<<" rns="<<now<<" old_bytes="<<reference.byteSize()<<" rns_bytes="<<direct.byteSize()<<" identical\n";
  }
  *input=*saved;
 };
 compare.template operator()<385>();compare.template operator()<398>();compare.template operator()<121>();
 if constexpr(Nbit<=4){
  // Compile and construct the existing default wide-primary backend too.
  PublicAutoGadgetRows<P,Q,64> legacy_rows;
  PreparedPublicAutoNTT<P,Q,64> legacy(legacy_rows);
  CKKSPreparedGadgetRNS32<P,Q,Bits> gadget(*rows);
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
int main(int argc,char**){check<4>();check<9>();if(argc>1)check<14>();}
