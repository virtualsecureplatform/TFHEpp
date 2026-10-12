#pragma once
#include "public_auto_ntt.hpp"
#include <map>
namespace TFHEpp {
template<uint32_t HighFactorQ, class Schedule, class PackedSchedule>
inline void CKKSNarrowUsePublicNTTAutomorphisms(CKKSNarrowBootstrapKey<Schedule>& key,
    const CKKSNarrowPackedKey<PackedSchedule>& packed) {
 using P=typename Schedule::Param;constexpr auto Q=CKKSNarrowPackedKey<PackedSchedule>::log_q;
 static_assert(Q>HighFactorQ && Q>Schedule::factor_log_q && Q>Schedule::output_log_q);
 using Cache=PreparedPublicAutoNTT<P,Q>;
 auto prepared=std::make_shared<std::map<uint32_t,std::shared_ptr<Cache>>>();
 PublicPackedAutoDeriver<P,PackedSchedule> derive(packed);
 auto attach=[&](uint32_t exponent) {
  if(prepared->contains(exponent))return;
  auto [rows,steps]=derive.derive(exponent);
  prepared->emplace(exponent,std::make_shared<Cache>(*rows));
 };
 attach(2*P::n-1);
 for(size_t ell=0;ell<Schedule::tail_blocks_per_pack;++ell)
  attach(ckks_detail::rotationAutomorphismForSteps<P>(2*(ell+1)*Schedule::tail_width));
 key.automorphism_evaluator=[prepared](TRLWE<P>& result,const TRLWE<P>& input,uint32_t level,uint32_t exponent) {
  const auto found=prepared->find(exponent);if(found==prepared->end())throw std::invalid_argument("unprepared narrow automorphism exponent");
  const auto& cache=*found->second;
  if(level==Schedule::factor_log_q)cache.template apply<Schedule::factor_log_q>(result,input);
  else if(level==Schedule::output_log_q)cache.template apply<Schedule::output_log_q>(result,input);
  else if(level==HighFactorQ)cache.template apply<HighFactorQ>(result,input);
  else throw std::invalid_argument("unexpected narrow automorphism data level");
 };
 key.factor_conjugation_fft.reset();key.factor_conjugation_projected_fft.reset();
 key.output_conjugation_fft.reset();key.output_conjugation_projected_fft.reset();
 key.tail_rotations_fft.reset();
}
}
