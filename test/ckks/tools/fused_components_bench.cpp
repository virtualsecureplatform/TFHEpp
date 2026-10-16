// Modes: fused (default), eight-terms, hoist, point-major.
// fused=1 selects the candidate; fused=0 is its matched control.
// Public-residue CPU experiment; no keys or encryption. Build from repo root:
// c++ -O3 -march=native -funroll-loops -fopenmp -std=c++20 -Iinclude \
//   -Ithirdparties/cereal/include test/ckks/tools/fused_components_bench.cpp -o /tmp/fused-components-bench
#include <modular_rns32.hpp>
#include <ckks/narrow_matrix.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <omp.h>
using namespace TFHEpp;
constexpr std::size_t N=16384,Count=528,Primes=17;
struct Selector{std::array<std::vector<uint32_t>,Primes> spectra;};
using Tables=std::array<std::array<std::vector<uint32_t>,Primes>,2>;
using Sums=std::array<std::array<std::vector<uint32_t>,Primes>,4>;
template<bool Fused,bool Packed,std::size_t Terms,bool HoistFirst=false,bool HoistSecond=false>
void multiply(Sums&sums,const Tables&base,const std::vector<Selector>&cache,const std::array<size_t,Count>&exponents){
 constexpr std::size_t Columns=Packed?4:8;
#pragma omp parallel for collapse(2) schedule(dynamic,1)
 for(size_t prime=0;prime<Primes;++prime)for(size_t begin=0;begin<N;begin+=256){
  modular_ntt32::Montgomery arithmetic(modular_rns32::primes[prime].value);
  for(size_t local=0;local<Count;local+=Terms){
   std::array<const uint32_t*,Terms> selectors0,selectors1;
   if constexpr(HoistFirst || HoistSecond)for(size_t t=0;t<Terms;++t){
    selectors0[t]=cache[local+t].spectra[prime].data();
    if constexpr(HoistSecond)selectors1[t]=selectors0[t]+N;
   }
   for(size_t i=begin;i<begin+256;++i){
   std::array<const uint32_t*,Terms> left,right;
   std::array<uint32_t,Terms> values0,values1;
   for(size_t t=0;t<Terms;++t){
    const auto position=narrow_detail::narrowNTTPermutationIndex<N>(exponents[local+t],i)*8;
    left[t]=base[0][prime].data()+position;
    if constexpr(Packed)right[t]=left[t]+4;else right[t]=base[1][prime].data()+position;
    if constexpr(HoistFirst)values0[t]=selectors0[t][i];else values0[t]=cache[local+t].spectra[prime][i];
    if constexpr(HoistSecond)values1[t]=selectors1[t][i];else values1[t]=cache[local+t].spectra[prime][N+i];
   }
   auto*a0=sums[0][prime].data()+i*Columns;auto*b0=sums[2][prime].data()+i*Columns;
   auto*a1=sums[1][prime].data()+i*Columns;auto*b1=sums[3][prime].data()+i*Columns;
   if constexpr(Fused)arithmetic.addDotPairComponents<Terms,Packed>(a0,b0,a1,b1,left,right,values0,values1,Columns);
   else if constexpr(Packed){arithmetic.addPackedDotPair<Terms>(a0,b0,left,values0);arithmetic.addPackedDotPair<Terms>(a1,b1,left,values1);}
   else{arithmetic.addDotPair<Terms>(a0,b0,left,right,values0,Columns);arithmetic.addDotPair<Terms>(a1,b1,left,right,values1,Columns);}
  }
  }
  for(auto&component:sums)for(size_t i=begin*Columns;i<(begin+256)*Columns;++i)component[prime][i]=arithmetic.canonical(component[prime][i]);
 }
}
template<bool Packed,std::size_t Terms>
void pointMajor(Sums&sums,const Tables&base,const std::vector<Selector>&cache,const std::array<size_t,Count>&exponents){
 constexpr size_t Columns=Packed?4:8;
#pragma omp parallel for collapse(2) schedule(dynamic,1)
 for(size_t prime=0;prime<Primes;++prime)for(size_t begin=0;begin<N;begin+=256){
  modular_ntt32::Montgomery arithmetic(modular_rns32::primes[prime].value);
  std::array<const uint32_t*,Count> sel0,sel1;
  for(size_t t=0;t<Count;++t){sel0[t]=cache[t].spectra[prime].data();sel1[t]=sel0[t]+N;}
  for(size_t i=begin;i<begin+256;++i){
   alignas(32) std::array<std::array<uint32_t,Columns>,4> local_sums{};
   for(size_t local=0;local<Count;local+=Terms){
    std::array<const uint32_t*,Terms> left,right;
    std::array<uint32_t,Terms> v0,v1;
    for(size_t t=0;t<Terms;++t){
     auto position=narrow_detail::narrowNTTPermutationIndex<N>(exponents[local+t],i)*8;
     left[t]=base[0][prime].data()+position;
     if constexpr(Packed)right[t]=left[t]+4;else right[t]=base[1][prime].data()+position;
     v0[t]=sel0[local+t][i];v1[t]=sel1[local+t][i];
    }
    arithmetic.addDotPairComponents<Terms,Packed>(local_sums[0].data(),local_sums[2].data(),local_sums[1].data(),local_sums[3].data(),left,right,v0,v1,Columns);
   }
   for(size_t c=0;c<4;++c)for(size_t j=0;j<Columns;++j)sums[c][prime][i*Columns+j]=arithmetic.canonical(local_sums[c][j]);
  }
 }
}
int main(int argc,char**argv){
 const std::string mode=argc>1?argv[1]:"fused";
 if(mode!="fused" && mode!="eight-terms" && mode!="hoist" && mode!="point-major")throw std::invalid_argument("mode must be fused, eight-terms, hoist or point-major");
 std::mt19937_64 random(20261017);std::vector<Selector> cache(Count);std::array<size_t,Count> exponents;
 for(size_t y=0;y<Count;++y){exponents[y]=modular_ntt::power(5,(N/2-y)%(N/2),2*N);
  for(size_t p=0;p<Primes;++p){auto&v=cache[y].spectra[p];v.resize(2*N);for(auto&w:v)w=random()%modular_rns32::primes[p].value;}}
 for(size_t columns:(mode=="eight-terms"?std::vector<size_t>{8}:std::vector<size_t>{4,8})){
  Tables base;Sums sums;
  for(size_t t=0;t<2;++t)for(size_t p=0;p<Primes;++p){base[t][p].resize(t==1&&columns==4?0:N*8);for(auto&w:base[t][p])w=random()%modular_rns32::primes[p].value;}
  for(auto&component:sums)for(auto&v:component)v.resize(N*columns);
  const auto run=[&](bool fused){
   for(auto&component:sums)for(auto&v:component)std::fill(v.begin(),v.end(),0);
   const auto start=std::chrono::steady_clock::now();
   if(mode=="point-major"){
    if(columns==4){if(fused)pointMajor<true,8>(sums,base,cache,exponents);else multiply<true,true,8>(sums,base,cache,exponents);}
    else{if(fused)pointMajor<false,4>(sums,base,cache,exponents);else multiply<true,false,4>(sums,base,cache,exponents);}
   }else if(mode=="eight-terms"){
    if(fused)multiply<true,false,8>(sums,base,cache,exponents);else multiply<true,false,4>(sums,base,cache,exponents);
   }else if(mode=="hoist"){
    if(columns==4){if(fused)multiply<true,true,8,true,true>(sums,base,cache,exponents);else multiply<true,true,8,true,false>(sums,base,cache,exponents);}
    else{if(fused)multiply<true,false,4,true,true>(sums,base,cache,exponents);else multiply<true,false,4,true,false>(sums,base,cache,exponents);}
   }else{
   if(columns==4){if(fused)multiply<true,true,8>(sums,base,cache,exponents);else multiply<false,true,8>(sums,base,cache,exponents);}
   else{if(fused)multiply<true,false,4>(sums,base,cache,exponents);else multiply<false,false,4>(sums,base,cache,exponents);}
   }
   return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  };
  run(false);const auto reference=sums;run(true);if(sums!=reference)throw std::runtime_error("warm fused matrix mismatch");
  for(bool fused:{false,true,true,false})for(int trial=0;trial<5;++trial){
   const auto seconds=run(fused);if(sums!=reference)throw std::runtime_error("timed fused matrix mismatch");
   std::cout<<"mode "<<mode<<" columns "<<columns<<" fused "<<fused<<" trial "<<trial<<" seconds "<<seconds<<" identical\n";
  }
 }
}
