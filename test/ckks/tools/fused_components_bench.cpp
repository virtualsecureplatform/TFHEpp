// Modes: fused (default), eight-terms, hoist, point-major, recurrence, tile32/64/128/512, traffic, slab.
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
#include <fstream>
#include <sstream>
#ifdef __linux__
#include <sys/mman.h>
#include <sys/prctl.h>
#include <cerrno>
#endif
using namespace TFHEpp;
constexpr std::size_t N=16384,Count=528,Primes=17;
struct Selector{std::array<std::vector<uint32_t>,Primes> spectra;};
using Tables=std::array<std::array<std::vector<uint32_t>,Primes>,2>;
using Sums=std::array<std::array<std::vector<uint32_t>,Primes>,4>;
template<bool Fused,bool Packed,std::size_t Terms,bool HoistFirst=false,bool HoistSecond=false,bool Recurrence=false,size_t Tile=256,class Row=Selector>
void multiply(Sums&sums,const Tables&base,const std::vector<Row>&cache,const std::array<size_t,Count>&exponents){
 constexpr std::size_t Columns=Packed?4:8;
#pragma omp parallel for collapse(2) schedule(dynamic,1)
 for(size_t prime=0;prime<Primes;++prime)for(size_t begin=0;begin<N;begin+=Tile){
  modular_ntt32::Montgomery arithmetic(modular_rns32::primes[prime].value);
  for(size_t local=0;local<Count;local+=Terms){
   std::array<const uint32_t*,Terms> selectors0,selectors1;
   if constexpr(HoistFirst || HoistSecond)for(size_t t=0;t<Terms;++t){
    selectors0[t]=cache[local+t].spectra[prime].data();
    if constexpr(HoistSecond)selectors1[t]=selectors0[t]+N;
   }
   std::array<size_t,Terms> positions;
   if constexpr(Recurrence)for(size_t t=0;t<Terms;++t)positions[t]=narrow_detail::narrowNTTPermutationIndex<N>(exponents[local+t],begin);
   for(size_t i=begin;i<begin+Tile;++i){
   std::array<const uint32_t*,Terms> left,right;
   std::array<uint32_t,Terms> values0,values1;
   for(size_t t=0;t<Terms;++t){
    size_t position;
    if constexpr(Recurrence){position=positions[t]*8;positions[t]=(positions[t]+exponents[local+t])&(N-1);}
    else position=narrow_detail::narrowNTTPermutationIndex<N>(exponents[local+t],i)*8;
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
  for(auto&component:sums)for(size_t i=begin*Columns;i<(begin+Tile)*Columns;++i)component[prime][i]=arithmetic.canonical(component[prime][i]);
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
// Read every selector word without any public-table arithmetic. These are
// traffic proxies, not lower bounds on actual DRAM traffic or kernel time.
template<bool Flat>uint32_t selectorReadProxy(const std::vector<Selector>&cache)
{
 std::vector<uint32_t> partial(Flat?Primes*Count:Primes*(N/256));
 struct XorWords {
#if defined(__AVX2__)
  __m256i vector=_mm256_setzero_si256();
#endif
  uint32_t scalar=0;
  void add(const uint32_t*data,size_t length){
   size_t i=0;
#if defined(__AVX2__)
   for(;i+8<=length;i+=8)vector=_mm256_xor_si256(vector,_mm256_loadu_si256(reinterpret_cast<const __m256i*>(data+i)));
#endif
   for(;i<length;++i)scalar^=data[i];
  }
  uint32_t finish(){
#if defined(__AVX2__)
   alignas(32) uint32_t words[8];_mm256_store_si256(reinterpret_cast<__m256i*>(words),vector);
   for(auto word:words)scalar^=word;
#endif
   return scalar;
  }
 };
 if constexpr(Flat){
#pragma omp parallel for collapse(2) schedule(static)
  for(size_t p=0;p<Primes;++p)for(size_t y=0;y<Count;++y){
   XorWords value;value.add(cache[y].spectra[p].data(),2*N);
   partial[p*Count+y]=value.finish();
  }
 }else{
#pragma omp parallel for collapse(2) schedule(dynamic,1)
  for(size_t p=0;p<Primes;++p)for(size_t begin=0;begin<N;begin+=256){
   XorWords value;
   for(size_t y=0;y<Count;++y)for(size_t c=0;c<2;++c)value.add(cache[y].spectra[p].data()+c*N+begin,256);
   partial[p*(N/256)+begin/256]=value.finish();
  }
 }
 uint32_t result=0;for(auto word:partial)result^=word;return result;
}
#ifdef __linux__
struct SelectorSlab {
 void*mapping=MAP_FAILED;size_t mapped_bytes=0;uint32_t*data=nullptr;int advice_error=0;
 SelectorSlab()=default;SelectorSlab(const SelectorSlab&)=delete;SelectorSlab&operator=(const SelectorSlab&)=delete;
 void allocate(bool huge,size_t stride=2*N){
  constexpr size_t alignment=2*1024*1024;
  const size_t bytes=(Count*stride*sizeof(uint32_t)+alignment-1)&~size_t(alignment-1);
  mapped_bytes=bytes+alignment;
  mapping=mmap(nullptr,mapped_bytes,PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
  if(mapping==MAP_FAILED)throw std::bad_alloc();
  data=reinterpret_cast<uint32_t*>((reinterpret_cast<uintptr_t>(mapping)+alignment-1)&~uintptr_t(alignment-1));
  if(madvise(data,bytes,huge?MADV_HUGEPAGE:MADV_NOHUGEPAGE))advice_error=errno;
 }
 ~SelectorSlab(){if(mapping!=MAP_FAILED)munmap(mapping,mapped_bytes);}
};
struct ScopedBenchmarkTHP {
 int original=prctl(PR_GET_THP_DISABLE,0,0,0,0),enable_error=0;
 ScopedBenchmarkTHP(){if(original==1 && prctl(PR_SET_THP_DISABLE,0,0,0,0))enable_error=errno;}
 ~ScopedBenchmarkTHP(){if(original==1 && !enable_error)prctl(PR_SET_THP_DISABLE,1,0,0,0);}
};
struct SelectorView {std::array<std::span<uint32_t>,Primes> spectra;};
size_t anonymousHugeKiB(){
 std::ifstream input("/proc/self/smaps");std::string line;size_t total=0;
 while(std::getline(input,line))if(line.starts_with("AnonHugePages:")){std::istringstream row(line.substr(14));size_t value=0;row>>value;total+=value;}
 return total;
}
#endif
int main(int argc,char**argv){
 const std::string mode=argc>1?argv[1]:"fused";
 if(mode!="fused" && mode!="eight-terms" && mode!="hoist" && mode!="point-major" && mode!="recurrence" && mode!="tile32" && mode!="tile64" && mode!="tile128" && mode!="tile512" && mode!="traffic" && mode!="slab")throw std::invalid_argument("unsupported benchmark mode");
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
   const auto tiled=[&]<size_t Tile>(){
    if(columns==4){if(fused)multiply<true,true,8,false,false,false,Tile>(sums,base,cache,exponents);else multiply<true,true,8>(sums,base,cache,exponents);}
    else{if(fused)multiply<true,false,4,false,false,false,Tile>(sums,base,cache,exponents);else multiply<true,false,4>(sums,base,cache,exponents);}
   };
   if(mode=="tile32")tiled.template operator()<32>();
   else if(mode=="tile64")tiled.template operator()<64>();
   else if(mode=="tile128")tiled.template operator()<128>();
   else if(mode=="tile512")tiled.template operator()<512>();
   else if(mode=="recurrence"){
    if(columns==4){if(fused)multiply<true,true,8,false,false,true>(sums,base,cache,exponents);else multiply<true,true,8>(sums,base,cache,exponents);}
    else{if(fused)multiply<true,false,4,false,false,true>(sums,base,cache,exponents);else multiply<true,false,4>(sums,base,cache,exponents);}
   }else if(mode=="point-major"){
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
  if(mode=="slab"){
#ifdef __linux__
   ScopedBenchmarkTHP thp;
   std::array<SelectorSlab,Primes> plain,huge,plain_padded,huge_padded;
   std::vector<SelectorView> heap_view(Count),plain_view(Count),huge_view(Count),plain_padded_view(Count),huge_padded_view(Count);
   constexpr size_t padded_stride=2*N+16;
   const auto before=anonymousHugeKiB();size_t advice_failures=0;
   for(size_t p=0;p<Primes;++p){
    plain[p].allocate(false);huge[p].allocate(true);
    plain_padded[p].allocate(false,padded_stride);huge_padded[p].allocate(true,padded_stride);
    advice_failures+=(plain[p].advice_error!=0)+(huge[p].advice_error!=0)+(plain_padded[p].advice_error!=0)+(huge_padded[p].advice_error!=0);
    for(size_t y=0;y<Count;++y){
     heap_view[y].spectra[p]=cache[y].spectra[p];
     plain_view[y].spectra[p]={plain[p].data+y*2*N,2*N};huge_view[y].spectra[p]={huge[p].data+y*2*N,2*N};
     plain_padded_view[y].spectra[p]={plain_padded[p].data+y*padded_stride,2*N};huge_padded_view[y].spectra[p]={huge_padded[p].data+y*padded_stride,2*N};
     std::copy_n(cache[y].spectra[p].data(),2*N,plain_padded_view[y].spectra[p].data());
     std::copy_n(cache[y].spectra[p].data(),2*N,huge_padded_view[y].spectra[p].data());
     std::copy_n(cache[y].spectra[p].data(),2*N,plain_view[y].spectra[p].data());
     std::copy_n(cache[y].spectra[p].data(),2*N,huge_view[y].spectra[p].data());
    }
   }
   const auto after=anonymousHugeKiB();
   std::cout<<"slab metadata columns="<<columns<<" anon_huge_before_kib="<<before<<" anon_huge_after_kib="<<after<<" inherited_thp_disabled="<<thp.original<<" thp_enable_error="<<thp.enable_error<<" advice_failures="<<advice_failures<<'\n';
   run(true);const auto reference=sums;
   const auto with_view=[&](const auto&view){
    for(auto&component:sums)for(auto&v:component)std::fill(v.begin(),v.end(),0);
    const auto start=std::chrono::steady_clock::now();
    if(columns==4)multiply<true,true,8>(sums,base,view,exponents);
    else multiply<true,false,4>(sums,base,view,exponents);
    const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if(sums!=reference)throw std::runtime_error("selector slab layout changed matrix residues");return seconds;
   };
   with_view(heap_view);with_view(plain_view);with_view(huge_view);with_view(plain_padded_view);with_view(huge_padded_view);
   for(int kind:{0,1,2,3,4,5,5,4,3,2,1,0})for(int trial=0;trial<5;++trial){
    const auto seconds=kind==0?run(true):kind==1?with_view(heap_view):kind==2?with_view(plain_view):kind==3?with_view(huge_view):kind==4?with_view(plain_padded_view):with_view(huge_padded_view);
    if(sums!=reference)throw std::runtime_error("original matrix changed in slab control");
    std::cout<<"slab columns="<<columns<<" kind="<<(kind==0?"original":kind==1?"heap-view":kind==2?"plain-slab":kind==3?"huge-slab":kind==4?"plain-padded":"huge-padded")<<" trial="<<trial<<" seconds="<<seconds<<" identical\n";
   }
   continue;
#else
   throw std::runtime_error("selector slab experiment requires Linux");
#endif
  }
  if(mode=="traffic"){
   uint32_t expected=0;for(const auto&row:cache)for(const auto&spectrum:row.spectra)for(auto word:spectrum)expected^=word;
   run(true);const auto reference=sums;
   if(selectorReadProxy<true>(cache)!=expected || selectorReadProxy<false>(cache)!=expected)throw std::runtime_error("selector traffic warmup checksum differs");
   constexpr uint64_t bytes=uint64_t(Count)*Primes*2*N*sizeof(uint32_t);
   for(int kind:{0,1,2,2,1,0})for(int trial=0;trial<5;++trial){
    double seconds;
    if(kind==0){seconds=run(true);if(sums!=reference)throw std::runtime_error("traffic comparison changed matrix residues");}
    else{
     const auto start=std::chrono::steady_clock::now();
     const auto checksum=kind==1?selectorReadProxy<false>(cache):selectorReadProxy<true>(cache);
     seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
     if(checksum!=expected)throw std::runtime_error("selector traffic checksum differs");
    }
    std::cout<<"traffic columns="<<columns<<" kind="<<(kind==0?"matrix":kind==1?"tiled-read":"flat-read")<<" trial="<<trial<<" seconds="<<seconds<<" selector_bytes="<<bytes<<" logical_GBps="<<bytes/seconds/1e9<<" checked\n";
   }
   continue;
  }
  run(false);const auto reference=sums;run(true);if(sums!=reference)throw std::runtime_error("warm fused matrix mismatch");
  for(bool fused:{false,true,true,false})for(int trial=0;trial<5;++trial){
   const auto seconds=run(fused);if(sums!=reference)throw std::runtime_error("timed fused matrix mismatch");
   std::cout<<"mode "<<mode<<" columns "<<columns<<" fused "<<fused<<" trial "<<trial<<" seconds "<<seconds<<" identical\n";
  }
 }
}
