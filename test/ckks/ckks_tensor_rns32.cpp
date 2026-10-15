#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/tensor_rns32.hpp>
#include <iostream>
#include <random>
using namespace TFHEpp;
template<std::uint32_t Nbit,std::uint32_t Q>void check(){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;
 using T=typename P::T;
 auto a=std::make_unique<TRLWE<P>>(),b=std::make_unique<TRLWE<P>>();
 auto natural=std::make_unique<TRLWE3<P>>();
 auto expected=std::make_unique<TRLWE3<P>>(),actual=std::make_unique<TRLWE3<P>>(),wide=std::make_unique<TRLWE3<P>>(),garner=std::make_unique<TRLWE3<P>>(),scalar=std::make_unique<TRLWE3<P>>();
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
  double rns=0,full=0,legacy=0;
  const auto compact=[&]{t=std::chrono::steady_clock::now();CKKSTensorProductRNS32<P,Q>(*actual,*a,*b);rns=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
  const auto control=[&]{t=std::chrono::steady_clock::now();CKKSTensorProductRNS32<P,Q,false,false>(*wide,*a,*b);full=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
  const auto previous=[&]{t=std::chrono::steady_clock::now();CKKSTensorProductRNS32<P,Q,true,false>(*garner,*a,*b);legacy=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
  if(trial%2){compact();previous();control();}else{control();previous();compact();}
  CKKSTensorProductRNS32<P,Q,true,true,false>(*scalar,*a,*b);
  CKKSTensorProductRNS32<P,Q,true,true,true,false>(*natural,*a,*b);
  if(*natural!=*actual)throw std::runtime_error("bit-reversed tensor differs from natural-order transforms");
  if(*scalar!=*actual)throw std::runtime_error("batched RNS tensor encoding differs from scalar encoding");
  if(*garner!=*actual)throw std::runtime_error("direct CRT tensor differs from compact Garner CRT");
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
  std::cout<<"N="<<P::n<<" Q="<<Q<<" trial="<<trial<<" fft="<<fft<<" wide_crt="<<full<<" compact_garner="<<legacy<<" compact_direct="<<rns<<" identical\n";
 }
}
template<std::uint32_t Q>void benchmarkBatch(){
 using P=narrow_benchmark::BenchmarkParam<14,7>;
 auto a=std::make_unique<TRLWE<P>>(),b=std::make_unique<TRLWE<P>>();
 auto expected=std::make_unique<TRLWE3<P>>(),out=std::make_unique<TRLWE3<P>>();
 std::mt19937_64 random(20261019+Q);
 for(auto*input:{a.get(),b.get()})for(auto&c:*input)for(auto&v:c)for(auto&w:v.limb)w=random();
 CKKSTensorProductRNS32<P,Q,true,true,false>(*expected,*a,*b);
 for(bool batch:{false,true,true,false})for(int trial=0;trial<5;++trial){
  const auto start=std::chrono::steady_clock::now();
  if(batch)CKKSTensorProductRNS32<P,Q>(*out,*a,*b);
  else CKKSTensorProductRNS32<P,Q,true,true,false>(*out,*a,*b);
  const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if(*out!=*expected)throw std::runtime_error("timed batched encoding tensor mismatch");
  std::cout<<"tensor q="<<Q<<" batch="<<batch<<" trial="<<trial<<" seconds="<<seconds<<" identical\n";
 }
}
template<std::uint32_t Q>void benchmarkOrdering(){
 using P=narrow_benchmark::BenchmarkParam<14,7>;
 auto a=std::make_unique<TRLWE<P>>(),b=std::make_unique<TRLWE<P>>();
 auto reference=std::make_unique<TRLWE3<P>>(),out=std::make_unique<TRLWE3<P>>();
 std::mt19937_64 random(20261017+Q);
 for(auto*input:{a.get(),b.get()})for(auto&component:*input)for(auto&v:component)for(auto&w:v.limb)w=random();
 CKKSTensorProductRNS32<P,Q,true,true,true,false>(*reference,*a,*b);
 // Warm both schedules, then symmetric blocks of measured repetitions.
 for(bool order:{false,true}){
  if(order)CKKSTensorProductRNS32<P,Q>(*out,*a,*b);
  else CKKSTensorProductRNS32<P,Q,true,true,true,false>(*out,*a,*b);
  if(*out!=*reference)throw std::runtime_error("warm tensor order mismatch");
 }
 for(bool order:{false,true,true,false})for(int trial=0;trial<5;++trial){
  const auto start=std::chrono::steady_clock::now();
  if(order)CKKSTensorProductRNS32<P,Q>(*out,*a,*b);
  else CKKSTensorProductRNS32<P,Q,true,true,true,false>(*out,*a,*b);
  const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if(*out!=*reference)throw std::runtime_error("timed tensor order mismatch");
  std::cout<<"tensor q="<<Q<<" bit_reversed="<<order<<" trial="<<trial<<" seconds="<<seconds<<" identical\n";
 }
}
int main(int argc,char**){if(argc>3){benchmarkOrdering<165>();benchmarkOrdering<209>();benchmarkOrdering<253>();benchmarkOrdering<297>();benchmarkOrdering<341>();benchmarkOrdering<385>();return 0;}if(argc>2){benchmarkBatch<165>();benchmarkBatch<209>();benchmarkBatch<253>();benchmarkBatch<297>();benchmarkBatch<341>();benchmarkBatch<385>();return 0;}check<4,1>();check<4,64>();check<4,128>();check<4,256>();check<4,384>();check<4,93>();check<4,385>();check<4,440>();check<9,93>();check<9,165>();check<9,209>();check<9,253>();check<9,297>();check<9,341>();check<9,385>();if(argc>1){check<14,165>();check<14,209>();check<14,253>();check<14,297>();check<14,341>();check<14,385>();}}
