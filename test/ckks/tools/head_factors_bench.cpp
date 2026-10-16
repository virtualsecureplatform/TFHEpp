// Deterministic coefficient fixtures; no keys generated or encrypted.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <sys/resource.h>
#include <iostream>
#include <iomanip>
using namespace TFHEpp;
template<unsigned Nbit>void run(std::size_t columns){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;
 using Factor=CKKSCiphertext<P,385,44>;
 constexpr std::size_t Blocks=31;
 std::vector<std::vector<Factor>> source(Blocks,std::vector<Factor>(columns));
 std::mt19937_64 random(20261101+Nbit+columns);
 for(auto&group:source)for(auto&v:group)for(auto&c:v.ct)for(auto&w:c)for(auto&limb:w.limb)limb=random();
 const auto build=[&](bool append){
  std::vector<std::vector<Factor>> head(columns);
  for(auto&column:head){if(!append)column.resize(Blocks);column.reserve(Blocks+1);}
  for(std::size_t x=0;x<Blocks;++x)for(std::size_t c=0;c<columns;++c){
   if(append)head[c].push_back(source[x][c]);else head[c][x]=source[x][c];
  }
  for(std::size_t c=0;c<columns;++c){head[c].emplace_back();head[c].back().ct[1][0]=typename P::T{91+c};}
  return head;
 };
 for(int mode:{0,1,0,1,1,0})for(int trial=0;trial<3;++trial){
  rusage before{},after{};getrusage(RUSAGE_SELF,&before);
  auto start=std::chrono::steady_clock::now();auto actual=build(mode);
  auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  getrusage(RUSAGE_SELF,&after);
  if(actual.size()!=columns)throw std::runtime_error("head column count differs");
  for(std::size_t c=0;c<columns;++c){
   if(actual[c].size()!=Blocks+1||actual[c].capacity()<Blocks+1)throw std::runtime_error("head leaf shape differs");
   for(std::size_t x=0;x<Blocks;++x)if(actual[c][x].ct!=source[x][c].ct)throw std::runtime_error("head leaf order or bits differ");
   auto beta=std::make_unique<Factor>();beta->ct[1][0]=typename P::T{91+c};
   if(actual[c].back().ct!=beta->ct)throw std::runtime_error("head beta position or zero component differs");
  }
  std::cout<<"N="<<P::n<<" columns="<<columns<<" mode="<<mode<<" trial="<<trial<<" seconds="<<seconds<<" minor_faults="<<after.ru_minflt-before.ru_minflt<<" identical\n";
 }
 std::cout<<"removed_placeholder_bytes N="<<P::n<<" columns="<<columns<<" per_pass="<<Blocks*columns*sizeof(Factor)<<'\n';
}
int main(){std::cout<<std::setprecision(17);run<9>(0);run<9>(1);run<9>(4);run<14>(1);run<14>(4);run<14>(7);}
