// Standalone CPU experiment; deterministic public residue fixtures, no keys.
// Build from the repository root:
// c++ -O3 -march=native -funroll-loops -fopenmp -std=c++20 -Iinclude \
//   -Ithirdparties/cereal/include test/ckks/tools/retained_matrix_bench.cpp -o /tmp/retained-matrix-bench
// Reports use mode0 for40-row chunks; modes32/64/128/256 consume all528 rows
// per point tile of that size. The arithmetic always uses four-term dots.
#include <modular_rns32.hpp>
#include <ckks/narrow_matrix.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <omp.h>
using namespace TFHEpp;
constexpr std::size_t N=16384,Count=528,Primes=17;
struct Selector {std::array<std::vector<uint32_t>,Primes> spectra;};
using Tables=std::array<std::array<std::vector<uint32_t>,Primes>,2>;
using Sums=std::array<std::array<std::vector<uint32_t>,Primes>,4>;
template<size_t Terms,size_t Tile=256,bool All=false> void multiply(Sums& sums,const Tables& base,const std::vector<Selector>& cache,
 const std::array<size_t,Count>& exponents,size_t stride){
for(size_t first_y=0;first_y<Count;first_y+=(All?Count:40)){
 const auto end_y=std::min(first_y+(All?Count:40),Count);
#pragma omp parallel for collapse(3) schedule(dynamic,1)
 for(size_t prime=0;prime<Primes;++prime)for(size_t component=0;component<2;++component)for(size_t begin=0;begin<N;begin+=Tile){
  modular_ntt32::Montgomery arithmetic(modular_rns32::primes[prime].value);
  auto*first=sums[component][prime].data();auto*second=sums[2+component][prime].data();
  for(size_t local=first_y;local<end_y;local+=Terms){
   std::array<const uint32_t*,Terms> selectors;
   for(size_t t=0;t<Terms;++t)selectors[t]=cache[local+t].spectra[prime].data()+component*N;
   for(size_t i=begin;i<begin+Tile;++i){
    std::array<const uint32_t*,Terms> left,right;std::array<uint32_t,Terms> values;
    for(size_t t=0;t<Terms;++t){auto position=narrow_detail::narrowNTTPermutationIndex<N>(exponents[local+t],i)*stride;
     left[t]=base[0][prime].data()+position;right[t]=base[1][prime].data()+position;values[t]=selectors[t][i];}
    arithmetic.addDotPair<Terms>(first+i*stride,second+i*stride,left,right,values,stride);
   }
  }
  for(size_t i=begin*stride;i<(begin+Tile)*stride;++i){first[i]=arithmetic.canonical(first[i]);second[i]=arithmetic.canonical(second[i]);}
 }
}
}
int main(){
 std::mt19937_64 random(20261017);std::vector<Selector> cache(Count);std::array<size_t,Count> exponents;
 for(size_t y=0;y<Count;++y){exponents[y]=modular_ntt::power(5,(N/2-y)%(N/2),2*N);
  for(size_t p=0;p<Primes;++p){auto&v=cache[y].spectra[p];v.resize(2*N);for(auto&w:v)w=random()%modular_rns32::primes[p].value;}}
 for(size_t stride:{8}){
  Tables base;Sums sums;
  for(auto&t:base)for(size_t p=0;p<Primes;++p){t[p].resize(N*stride);for(auto&v:t[p])v=random()%modular_rns32::primes[p].value;}
  for(auto&t:sums)for(auto&v:t)v.resize(N*stride);
  auto reference=sums;
  for(int mode:{0,256,128,64,32,32,64,128,256,0})for(int trial=0;trial<5;++trial){
   for(auto&t:sums)for(auto&v:t)std::fill(v.begin(),v.end(),0);
   auto start=std::chrono::steady_clock::now();
   if(mode==0)multiply<4>(sums,base,cache,exponents,stride);
   else if(mode==256)multiply<4,256,true>(sums,base,cache,exponents,stride);
   else if(mode==128)multiply<4,128,true>(sums,base,cache,exponents,stride);
   else if(mode==64)multiply<4,64,true>(sums,base,cache,exponents,stride);
   else multiply<4,32,true>(sums,base,cache,exponents,stride);
   double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
   if(mode==0&&trial==0)reference=sums;
   if(sums!=reference)throw std::runtime_error("matrix dot mismatch");
   std::cout<<"stride "<<stride<<" terms "<<mode<<" trial "<<trial<<" seconds "<<seconds<<" identical\n";
  }
 }
}
