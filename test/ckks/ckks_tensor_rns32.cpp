#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/tensor_rns32.hpp>
#include <iostream>
#include <random>
using namespace TFHEpp;
template<std::uint32_t Nbit,std::uint32_t Q>void check(){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;
 using T=typename P::T;
 auto a=std::make_unique<TRLWE<P>>(),b=std::make_unique<TRLWE<P>>();
 auto expected=std::make_unique<TRLWE3<P>>(),actual=std::make_unique<TRLWE3<P>>(),wide=std::make_unique<TRLWE3<P>>();
 std::mt19937_64 random(20261016+Q);
 for(int trial=0;trial<5;++trial){
  for(auto *input:{a.get(),b.get()})for(auto &component:*input)for(auto &v:component){
   for(auto &word:v.limb)word=random();
   if(trial==1)v=T{0}-T{1};
   if(trial==2)v=(T{1}<<(Q-1))+T{1};
   if(trial==3)v={};
   if(trial==4)v=(T{1}<<((Q+1)/2))-T{1};
  }
  auto t=std::chrono::steady_clock::now();
  CKKSTensorProductRescale<P,Q,Q,0>(*expected,*a,*b);
  const auto fft=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
  double rns=0,full=0;
  const auto compact=[&]{t=std::chrono::steady_clock::now();CKKSTensorProductRNS32<P,Q>(*actual,*a,*b);rns=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
  const auto control=[&]{t=std::chrono::steady_clock::now();CKKSTensorProductRNS32<P,Q,false>(*wide,*a,*b);full=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
  if(trial%2){compact();control();}else{control();compact();}
  if(*wide!=*actual)throw std::runtime_error("compact CRT tensor differs from full-width CRT");
  if(*expected!=*actual)throw std::runtime_error("tensor differs from FFT reference");
  if constexpr(Nbit<=4){
   auto oracle=std::make_unique<TRLWE3<P>>();
   const auto add=[&](std::size_t c,std::size_t x,std::size_t y){
    for(std::size_t i=0;i<P::n;++i)for(std::size_t j=0;j<P::n;++j){
     const auto v=(*a)[x][i]*(*b)[y][j];const auto k=(i+j)%P::n;
     (*oracle)[c][k]=ckks_detail::reduceToLevel<P,Q>((*oracle)[c][k]+(i+j<P::n?v:-v));
    }
   };
   add(0,0,1);add(0,1,0);add(1,1,1);add(2,0,0);
   if(*oracle!=*actual)throw std::runtime_error("tensor differs from schoolbook oracle");
  }
  std::cout<<"N="<<P::n<<" Q="<<Q<<" trial="<<trial<<" fft="<<fft<<" wide_crt="<<full<<" compact_crt="<<rns<<" identical\n";
 }
}
int main(int argc,char**){check<4,1>();check<4,64>();check<4,128>();check<4,256>();check<4,384>();check<4,93>();check<4,385>();check<4,440>();check<9,93>();check<9,165>();check<9,209>();check<9,253>();check<9,297>();check<9,341>();check<9,385>();if(argc>1){check<14,165>();check<14,209>();check<14,253>();check<14,297>();check<14,341>();check<14,385>();}}
