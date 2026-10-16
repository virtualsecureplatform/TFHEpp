// Deterministic public residues; no keys generated or encrypted.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_factor_rns32.hpp>
#include <ckks/factor_finish_rns32.hpp>
#include <iostream>
#include <iomanip>
using namespace TFHEpp;
template<class P>struct Schedule {using Param=P;static constexpr unsigned selector_log_q=440,plain_log_delta=42;};
template<class S,std::size_t Count,class CRT>
void reference(std::vector<TRLWE<typename S::Param>>&first,std::vector<TRLWE<typename S::Param>>&second,
 const std::vector<std::array<std::vector<uint32_t>,Count>>&sums,std::size_t columns,std::size_t stride,const CRT&crt){
 using P=typename S::Param;const auto&plans=narrow_detail::factorPlansRNS32<P,Count,false>();
 std::vector<std::array<TRLWE<P>,2>> accum(columns);
#pragma omp parallel for collapse(3) schedule(static)
 for(std::size_t t=0;t<2;++t)for(std::size_t c=0;c<=P::k;++c)for(std::size_t column=0;column<columns;++column){
  std::array<std::vector<uint32_t>,Count> spectra;
  for(std::size_t p=0;p<Count;++p){
   auto&v=spectra[p];v.resize(P::n);const auto&sum=sums[t*(P::k+1)+c][p];
   for(std::size_t i=0;i<P::n;++i)v[i]=sum[i*stride+column];
   plans[p].inverse(v);
  }
  std::array<std::span<uint32_t>,Count> residues;
  for(std::size_t p=0;p<Count;++p)residues[p]=spectra[p];
  std::vector<typename P::T> coefficients(P::n);crt.reconstruct(residues,coefficients);
  for(std::size_t i=0;i<P::n;++i)accum[column][t][c][i]=ckks_detail::reduceToLevel<P,S::selector_log_q>(coefficients[i]);
 }
 first.resize(columns);second.resize(columns);
#pragma omp parallel for collapse(3) schedule(dynamic,1)
 for(std::size_t column=0;column<columns;++column)for(std::size_t c=0;c<=P::k;++c)for(std::size_t begin=0;begin<P::n;begin+=256)
  for(std::size_t i=begin;i<std::min(begin+256,std::size_t(P::n));++i){
   first[column][c][i]=ckks_detail::roundedLevelRightShift<P,S::selector_log_q,S::plain_log_delta>(accum[column][0][c][i]);
   second[column][c][i]=ckks_detail::roundedLevelRightShift<P,S::selector_log_q,S::plain_log_delta>(accum[column][1][c][i]);
  }
}
template<unsigned Nbit>void run(std::size_t columns){
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;using S=Schedule<P>;constexpr std::size_t Count=17;
 const auto stride=columns==3?4:columns;
 std::mt19937_64 random(20261029+Nbit+columns);
 std::vector<std::array<std::vector<uint32_t>,Count>> sums(4);
 for(auto&row:sums)for(std::size_t p=0;p<Count;++p){row[p].resize(stride*P::n);for(auto&v:row[p])v=random()%modular_rns32::primes[p].value;}
 const modular_rns32::DirectCenteredCRT<P::T::limbs,Count> crt;
 const auto&plans=narrow_detail::factorPlansRNS32<P,Count,false>();
 std::vector<TRLWE<P>> ef,es;reference<S>(ef,es,sums,columns,stride,crt);
 if constexpr(Nbit<=9){
  const modular_rns32::CenteredCRT<P::T::limbs,Count> garner;std::vector<TRLWE<P>> gf,gs;
  reference<S>(gf,gs,sums,columns,stride,garner);
  if(gf!=ef||gs!=es)throw std::runtime_error("independent Garner factor finish differs");
 }
 std::vector<CKKSCiphertext<P,398,57>> wrapped;
 std::vector<TRLWE<P>> wrapped_second;
 narrow_detail::inverseRoundFactorBatchRNS32<S>(wrapped,wrapped_second,sums,columns,stride,plans,crt);
 if(wrapped_second!=es)throw std::runtime_error("wrapped factor second stream differs");
 for(std::size_t j=0;j<columns;++j)if(wrapped[j].ct!=ef[j])throw std::runtime_error("wrapped factor first stream differs");
 for(int mode:{0,1,0,1,1,0})for(int trial=0;trial<3;++trial){
  std::vector<TRLWE<P>> first,second;
  auto start=std::chrono::steady_clock::now();
  if(mode)narrow_detail::inverseRoundFactorBatchRNS32<S>(first,second,sums,columns,stride,plans,crt);
  else reference<S>(first,second,sums,columns,stride,crt);
  auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if(first!=ef||second!=es)throw std::runtime_error("factor finish changed coefficients");
  std::cout<<"N="<<P::n<<" columns="<<columns<<" mode="<<mode<<" trial="<<trial<<" seconds="<<seconds<<" identical\n";
 }
 std::vector<TRLWE<P>> first,second;
#pragma omp parallel num_threads(2)
 {
#pragma omp single
  narrow_detail::inverseRoundFactorBatchRNS32<S>(first,second,sums,columns,stride,plans,crt);
 }
 if(first!=ef||second!=es)throw std::runtime_error("nested factor finish differs");
}
int main(){std::cout<<std::setprecision(17);run<9>(0);run<9>(1);run<9>(3);run<9>(4);run<14>(3);run<14>(4);}
