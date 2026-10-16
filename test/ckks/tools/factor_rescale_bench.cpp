// Deterministic coefficient fixtures; no key generation or encryption.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_factor_rns32.hpp>
#include <iostream>
#include <iomanip>
using namespace TFHEpp;
template<class P>struct SourceSchedule{using Param=P;static constexpr unsigned factor_log_q=398,log_delta=57;using Factor=CKKSCiphertext<P,398,57>;};
template<class P>struct TargetSchedule{using Param=P;static constexpr unsigned factor_log_q=385,log_delta=44;using Factor=CKKSCiphertext<P,385,44>;};
template<unsigned Nbit>void run(std::size_t columns){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;
 using Source=SourceSchedule<P>;using Target=TargetSchedule<P>;
 std::vector<typename Source::Factor> input(columns);
 std::mt19937_64 random(20261028);
 for(auto&v:input)for(auto&c:v.ct)for(auto&w:c){for(auto&limb:w.limb)limb=random();}
 std::vector<typename Target::Factor> expected,actual;
 narrow_detail::rescaleFactorBatchRNS32<Target,Source>(expected,input,false);
 // Exercise a reused output vector and a caller-owned outer team as well.
 actual.resize(columns+1);
#pragma omp parallel num_threads(2)
 {
#pragma omp single
  narrow_detail::rescaleFactorBatchRNS32<Target,Source>(actual,input,true);
 }
 if(actual.size()!=expected.size())throw std::runtime_error("factor projection shape differs");
 for(std::size_t j=0;j<columns;++j)if(actual[j].ct!=expected[j].ct)throw std::runtime_error("nested factor projection differs");
 for(bool allocate:{false,true})for(int mode:{0,1,1,0}){
  double seconds=0;
  for(int trial=0;trial<36;++trial){
   if(allocate)std::vector<typename Target::Factor>{}.swap(actual);
   auto start=std::chrono::steady_clock::now();
   narrow_detail::rescaleFactorBatchRNS32<Target,Source>(actual,input,mode);
   seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
   for(std::size_t j=0;j<columns;++j)if(actual[j].ct!=expected[j].ct)throw std::runtime_error("parallel factor projection differs");
  }
  std::cout<<"N="<<P::n<<" columns="<<columns<<" allocate="<<allocate<<" mode="<<mode<<" seconds="<<seconds<<" identical\n";
 }
}
int main(){std::cout<<std::setprecision(17);for(auto n:{0u,1u,3u,4u,7u,14u})run<9>(n);for(auto n:{1u,3u,4u,7u,14u})run<14>(n);}
