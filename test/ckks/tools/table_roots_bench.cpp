// Deterministic public phases; no keys generated or encrypted.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <cstring>
#include <iostream>
#include <iomanip>
using namespace TFHEpp;
template<unsigned Nbit>void check(){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;
 using S=CKKSNarrowBootstrapSchedule<P,48,427,44,42,31,true,13,true>;
 std::array<double,P::n> a;std::mt19937_64 random(20261030+Nbit);
 for(auto&v:a)v=std::ldexp(double(random()>>11),-53)-0.5;
 std::array<double,8> edges{0.0,-0.0,0.5,-0.5,std::nextafter(0.5,0.0),std::nextafter(-0.5,0.0),std::ldexp(1.0,-48),-std::ldexp(1.0,-48)};
 std::copy(edges.begin(),edges.end(),a.begin());
 for(std::uint32_t frequency:{1u,2u,7u,9u,10u,15u,40u,~0u}){
  auto roots=std::make_unique<narrow_detail::NarrowTableRoots<P>>(a,frequency);
  for(std::size_t x=0;x<S::blocks;++x){
   auto old=narrow_detail::headTable<S>(a,x,frequency),now=narrow_detail::headTable<S>(*roots,x);
   if(std::memcmp(old.data(),now.data(),sizeof(old)))throw std::runtime_error("cached head roots differ in bits");
  }
  for(std::size_t q=0;q<S::tail_packs;++q){
   auto old=narrow_detail::tailTable<S>(a,q,frequency),now=narrow_detail::tailTable<S>(*roots,q);
   if(std::memcmp(old.data(),now.data(),sizeof(old)))throw std::runtime_error("cached tail roots differ in bits");
  }
 }
 std::cout<<"check N="<<P::n<<" frequencies=8 all head/tail tables bit-identical\n";
}
void benchmark(){
 using P=narrow_benchmark::BenchmarkParam<14,7>;
 using S=CKKSNarrowBootstrapSchedule<P,48,427,44,42,31,true,13,true>;
 constexpr std::array<std::uint32_t,4> frequencies{7,9,10,15};
 std::array<double,P::n> a;std::mt19937_64 random(20261031);
 for(auto&v:a)v=std::ldexp(double(random()>>11),-53)-0.5;
 const auto build=[&](bool reuse){
  std::vector<narrow_detail::NarrowTableRoots<P>> roots;
  if(reuse){roots.reserve(frequencies.size());for(auto f:frequencies)roots.emplace_back(a,f);}
  std::vector<CKKSSlotVector<P>> tables;tables.reserve(frequencies.size()*(S::blocks+S::tail_packs));
  for(std::size_t column=0;column<frequencies.size();++column){
   for(std::size_t x=0;x<S::blocks;++x)tables.push_back(reuse?narrow_detail::headTable<S>(roots[column],x):narrow_detail::headTable<S>(a,x,frequencies[column]));
   for(std::size_t q=0;q<S::tail_packs;++q)tables.push_back(reuse?narrow_detail::tailTable<S>(roots[column],q):narrow_detail::tailTable<S>(a,q,frequencies[column]));
  }
  return tables;
 };
 const auto reference=build(false);
 for(int mode:{0,1,0,1,1,0})for(int trial=0;trial<3;++trial){
  auto start=std::chrono::steady_clock::now();auto actual=build(mode);
  auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if(actual.size()!=reference.size()||std::memcmp(actual.data(),reference.data(),actual.size()*sizeof(actual.front())))throw std::runtime_error("timed root tables differ in bits");
  std::cout<<"mode="<<mode<<" trial="<<trial<<" seconds="<<seconds<<" bit-identical\n";
 }
 std::cout<<"root_calls_per_column original="<<S::blocks*S::slots+S::blocks*2*S::tail_width<<" reused="<<2*P::n<<" extra_payload_bytes_per_column="<<sizeof(narrow_detail::NarrowTableRoots<P>)<<'\n';
}
int main(){std::cout<<std::setprecision(17);check<9>();check<14>();benchmark();}
