#pragma once
#include "public_auto_ntt.hpp"
#include <map>
namespace TFHEpp {
namespace narrow_detail {
// Route each supported data level through one immutable public key cache.
// Optional series/intermediate levels reuse conjugation in reconstruction.
template<class Schedule,uint32_t HighFactorQ,uint32_t ExtraDataQ=0,uint32_t IntermediateDataQ=0,class Cache>
inline void applyPublicAutoAtLevel(TRLWE<typename Schedule::Param>&result,
    const TRLWE<typename Schedule::Param>&input,uint32_t level,const Cache&cache){
 if(level==Schedule::factor_log_q)cache.template apply<Schedule::factor_log_q>(result,input);
 else if(level==Schedule::output_log_q)cache.template apply<Schedule::output_log_q>(result,input);
 else if(level==HighFactorQ)cache.template apply<HighFactorQ>(result,input);
 else {
  if constexpr(IntermediateDataQ>0){
   if(level==IntermediateDataQ){cache.template apply<IntermediateDataQ>(result,input);return;}
  }
  if constexpr(ExtraDataQ>0){
   if(level==ExtraDataQ)cache.template apply<ExtraDataQ>(result,input);
   else throw std::invalid_argument("unexpected narrow automorphism data level");
  }else throw std::invalid_argument("unexpected narrow automorphism data level");
 }
}
// Reuse the same registered-level dispatcher for ordinary and regrouped maps.
template<class Schedule,uint32_t HighFactorQ,uint32_t ExtraDataQ,uint32_t IntermediateDataQ=0,class Cache>
inline auto makePublicAutoEvaluator(std::shared_ptr<std::map<uint32_t,std::shared_ptr<Cache>>> prepared)
{
 using P=typename Schedule::Param;
 return [prepared](TRLWE<P>&result,const TRLWE<P>&input,uint32_t level,uint32_t exponent){
  const auto found=prepared->find(exponent);
  if(found==prepared->end())throw std::invalid_argument("unprepared narrow automorphism exponent");
  applyPublicAutoAtLevel<Schedule,HighFactorQ,ExtraDataQ,IntermediateDataQ>(result,input,level,*found->second);
 };
}

}
// The optional setup observer can prepare another public row representation
// before temporary derived rows are released; no new encryption is needed.
template<uint32_t HighFactorQ, uint32_t Primary=64, uint32_t Secondary=32, bool DirectRNS=false, uint32_t ExtraDataQ=0,
         class Schedule, class PackedSchedule, class OnDerived=std::nullptr_t>
inline void CKKSNarrowUsePublicNTTAutomorphisms(CKKSNarrowBootstrapKey<Schedule>& key,
    const CKKSNarrowPackedKey<PackedSchedule>& packed,OnDerived on_derived=nullptr) {
 using P=typename Schedule::Param;constexpr auto Q=CKKSNarrowPackedKey<PackedSchedule>::log_q;
 static_assert(ExtraDataQ==0 || ExtraDataQ<Q);
 static_assert(Q>HighFactorQ && Q>Schedule::factor_log_q && Q>Schedule::output_log_q);
 using Cache=PreparedPublicAutoNTT<P,Q,Primary,Secondary,DirectRNS>;
 auto prepared=std::make_shared<std::map<uint32_t,std::shared_ptr<Cache>>>();
 PublicPackedAutoDeriver<P,PackedSchedule,Primary> derive(packed);
 auto attach=[&](uint32_t exponent) {
  if(prepared->contains(exponent))return;
  auto [rows,steps]=derive.derive(exponent);
  prepared->emplace(exponent,std::make_shared<Cache>(*rows));
  if constexpr(!std::is_same_v<OnDerived,std::nullptr_t>){const auto&public_rows=*rows;on_derived(public_rows);}
 };
 attach(2*P::n-1);
 for(size_t ell=0;ell<Schedule::tail_blocks_per_pack;++ell)
  attach(ckks_detail::rotationAutomorphismForSteps<P>(2*(ell+1)*Schedule::tail_width));
 key.automorphism_evaluator=narrow_detail::makePublicAutoEvaluator<Schedule,HighFactorQ,ExtraDataQ>(prepared);
 key.factor_conjugation_fft.reset();key.factor_conjugation_projected_fft.reset();
 key.output_conjugation_fft.reset();key.output_conjugation_projected_fft.reset();
 key.tail_rotations_fft.reset();
}
}
