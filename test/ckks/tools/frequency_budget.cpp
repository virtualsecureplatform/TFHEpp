// Offline public frequency and approximation analysis; no keys or ciphertexts.
// Enumerate positive root sets in1..40 with one sum/difference per frequency.
// Build: c++ -O3 -std=c++20 -Iinclude -Ithirdparties/cereal/include
//        test/ckks/tools/frequency_budget.cpp -o /tmp/frequency-budget
// Reports coverage, degree40 truncation and fresh windowed Tanh Remez fits.
#include <ckks/narrow_approximation.hpp>
#include <bit>
#include <cstdint>
#include <iostream>
#include <iomanip>
using namespace TFHEpp;
template<std::size_t K> struct Search {
 std::array<unsigned,K> current{},best{};unsigned range=0;std::uint64_t tested=0;
 void visit(unsigned at,unsigned first,std::uint64_t mask){
  if(at==K){++tested;const auto r=unsigned(std::countr_one(mask));if(r>range || (r==range && current.back()<best.back())){range=r;best=current;}return;}
  for(unsigned v=first;v<=40-(K-at-1);++v){
   auto next=mask|(std::uint64_t{1}<<(v-1));
   if(2*v<=40)next|=std::uint64_t{1}<<(2*v-1);
   for(unsigned j=0;j<at;++j){if(v+current[j]<=40)next|=std::uint64_t{1}<<(v+current[j]-1);next|=std::uint64_t{1}<<(v-current[j]-1);}
   current[at]=v;visit(at+1,v+1,next);
  }
 }
 void run(){visit(0,1,0);std::cout<<"{\"roots\":"<<K<<",\"tested\":"<<tested<<",\"range\":"<<range<<",\"basis\":[";for(unsigned j=0;j<K;++j)std::cout<<(j?",":"")<<best[j];std::cout<<"]}\n";}
};
int main(){
 std::cout<<std::setprecision(17);Search<3>{}.run();Search<4>{}.run();Search<5>{}.run();Search<6>{}.run();
 const auto approx=CKKSBuildNarrowWindowedTrigonometricApproximation([](double x){return std::tanh(x);},[](double x){const auto y=std::tanh(x);return 1-y*y;},[](double x){const auto y=std::tanh(x);return -2*y*(1-y*y);},{1,1,2,4,40},4,16,40);
 std::array<double,41> tail{},sampled{};for(unsigned k=40;k>0;--k)tail[k-1]=tail[k]+2*std::abs(approx.coefficients[k]);
 for(unsigned i=0;i<=131072;++i){
  const long double x=-4+8.L*i/131072;const auto angle=2*std::acos(-1.L)*x/16;
  const auto z=std::complex<long double>(std::cos(angle),std::sin(angle));auto power=z;
  long double sum=approx.coefficients[0].real();const auto target=std::tanh(x);
  for(unsigned k=1;k<=40;++k){sum+=2*std::real(std::complex<long double>(approx.coefficients[k])*power);sampled[k]=std::max(sampled[k],double(std::abs(sum-target)));power*=z;}
 }
 std::cout<<"{\"full_bound\":"<<*approx.uniform_error_bound<<",\"l1\":"<<tail[0]<<",\"truncations\":[";
 for(unsigned k=20;k<=40;++k)std::cout<<(k>20?",":"")<<"{\"degree\":"<<k<<",\"tail\":"<<tail[k]<<",\"sampled\":"<<sampled[k]<<",\"triangle_bound\":"<<*approx.uniform_error_bound+tail[k]<<"}";
 std::cout<<"]}\n";
 for(unsigned degree:{28,30,32,34,36,38,40}){
  const auto a=CKKSBuildNarrowWindowedTrigonometricApproximation([](double x){return std::tanh(x);},[](double x){const auto y=std::tanh(x);return 1-y*y;},[](double x){const auto y=std::tanh(x);return -2*y*(1-y*y);},{1,1,2,4,40},4,16,degree);
  double norm=std::abs(a.coefficients[0]);for(unsigned k=1;k<=degree;++k)norm+=2*std::abs(a.coefficients[k]);
  std::cout<<"{\"remez_degree\":"<<degree<<",\"sampled\":"<<a.max_sampled_error<<",\"uniform_bound\":"<<*a.uniform_error_bound<<",\"l1\":"<<norm<<"}\n";
 }

}
