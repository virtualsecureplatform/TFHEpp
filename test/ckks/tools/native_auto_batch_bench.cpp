// Complete native automorphisms on deterministic noncryptographic rows.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <ckks/public_auto_ntt.hpp>
#include <iomanip>
#include <iostream>
using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<14,7>;
using Cache=PreparedPublicAutoNTT<P,435,44,20,true,145>;
template<std::uint32_t DataQ>void benchmark(const Cache&key,std::size_t size){
 std::vector<TRLWE<P>> input(size),reference(size),out(size);std::mt19937_64 random(20261029+DataQ);
 for(auto&r:input)for(auto&c:r)for(auto&v:c)for(auto&w:v.limb)w=random();
 const auto threads=std::min<int>(size,omp_get_max_threads());
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
 for(std::size_t j=0;j<size;++j)key.template apply<DataQ>(reference[j],input[j]);
 const auto run=[&](int window){
  const auto start=std::chrono::steady_clock::now();
  if(!window){
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
   for(std::size_t j=0;j<size;++j)key.template apply<DataQ>(out[j],input[j]);
  }else key.template applyBatch<DataQ>(out,input,window);
  const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if(out!=reference)throw std::runtime_error("batched full automorphism changed exact ciphertext");
  return seconds;
 };
 for(int window:{0,4,8,32})run(window);
 for(int window:{0,4,8,32,32,8,4,0})for(int trial=0;trial<4;++trial)
  std::cout<<"data_q="<<DataQ<<" batch="<<size<<" window="<<window<<" trial="<<trial<<" seconds="<<run(window)<<" identical\n";
}
int main(){
 auto rows=std::make_unique<PublicAutoGadgetRows<P,435,44>>();rows->exponent=2*P::n-1;std::mt19937_64 random(20261029);
 for(auto&r:rows->ciphertexts)for(auto&c:r)for(auto&v:c){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,435>(v);}
 Cache key(*rows);std::cout<<std::setprecision(17);
 for(std::size_t size:{3,4,7,8}){benchmark<398>(key,size);benchmark<385>(key,size);benchmark<165>(key,size);}
}
