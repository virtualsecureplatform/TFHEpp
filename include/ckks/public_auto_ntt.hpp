#pragma once
#include <ckks/narrow_packed_file.hpp>
#include <ckks/prepared_dd_ntt.hpp>
#include <ckks/prepared_gadget_rns32.hpp>
namespace TFHEpp {
// Public encryptions of gadget*s^exponent; never a relinearization key.
template<class P, uint32_t Q, uint32_t Bits=64>
struct PublicAutoGadgetRows {
 static_assert(P::k==1 && is_multilimb_uint_v<typename P::T>);
 static_assert(Bits>0 && Bits<=64 && ((Q+Bits-1)/Bits)*Bits<=P::T::digits);
 static constexpr size_t rows=(Q+Bits-1)/Bits;
 uint32_t exponent=1;
 std::array<TRLWE<P>,rows> ciphertexts{};
};
// Regroup retained encrypted gadget rows without deriving or encrypting again.
template<uint32_t NewBits,class P,uint32_t Q,uint32_t OldBits>
inline auto makePublicAutoGadgetRowsWithPrimaryBase(const PublicAutoGadgetRows<P,Q,OldBits>&source)
{
 static_assert(NewBits>=OldBits && NewBits%OldBits==0);
 using Target=PublicAutoGadgetRows<P,Q,NewBits>;
 auto result=std::make_unique<Target>();result->exponent=source.exponent;
 for(std::size_t row=0;row<Target::rows;++row){
  const auto power=(Target::rows-row-1)*NewBits;
  result->ciphertexts[row]=source.ciphertexts[source.rows-1-power/OldBits];
 }
 return result;
}

template<class TargetP, class PackedSchedule, uint32_t Bits=64>
class PublicPackedAutoDeriver {
 using SourceP=typename PackedSchedule::Param;
 static constexpr uint32_t Q=CKKSNarrowPackedKey<PackedSchedule>::log_q;
 static_assert(TargetP::n==SourceP::n && TargetP::k==SourceP::k);
 static_assert(std::is_same_v<typename TargetP::T,typename SourceP::T>);
 const CKKSNarrowPackedKey<PackedSchedule>& packed_;
 std::array<std::unique_ptr<CKKSPreparedAutoKeyFFT<SourceP,Q>>,SourceP::nbit> caches_;
public:
 explicit PublicPackedAutoDeriver(const CKKSNarrowPackedKey<PackedSchedule>& packed):packed_(packed){}
 auto derive(uint32_t exponent) {
  if(exponent==0 || exponent>=2*SourceP::n || !(exponent&1))throw std::invalid_argument("invalid public automorphism exponent");
  std::vector<size_t> depths;uint32_t product=1;
  for(size_t bit=1;bit<=SourceP::nbit;++bit)if(((product^exponent)>>bit)&1) {
   depths.push_back(SourceP::nbit-bit);
   product=(uint64_t(product)*(1+(uint32_t{1}<<bit)))%(2*SourceP::n);
  }
  if(product!=exponent)throw std::logic_error("automorphism factorization failed");
  for(auto depth:depths)if(!caches_[depth])caches_[depth]=makeCKKSPreparedAutoKeyFFT<SourceP,Q>(packed_.substitutions[depth]);
  auto result=std::make_unique<PublicAutoGadgetRows<TargetP,Q,Bits>>();result->exponent=exponent;
  auto current=std::make_unique<TRLWE<SourceP>>(),next=std::make_unique<TRLWE<SourceP>>();
  for(size_t row=0;row<result->rows;++row) {
   for(auto& c:*current)c.fill({});
   (*current)[0][0]=ckks_detail::reduceToLevel<SourceP,Q>(-typename SourceP::T{1}<<((result->rows-row-1)*Bits));
   for(auto depth:depths) {
    CKKSEvalAuto<SourceP,Q>(*next,*current,1+(SourceP::n>>depth),*caches_[depth]);
    current.swap(next);
   }
   for(size_t c=0;c<=TargetP::k;++c)for(size_t i=0;i<TargetP::n;++i)result->ciphertexts[row][c][i]=(*current)[c][i];
  }
  return std::pair{std::move(result),depths.size()};
 }
};
template<class P,uint32_t Q,uint32_t Bits=64,uint32_t Secondary=32,bool DirectRNS=false,uint32_t ChunkBits=80>
class PreparedPublicAutoNTT {
 using Cache=std::conditional_t<DirectRNS,CKKSPreparedGadgetRNS32<P,Q,Bits,ChunkBits>,CKKSPreparedDDRelinKeyNTT<P,Q,Bits,Secondary>>;
 uint32_t exponent_;
 std::shared_ptr<Cache> cache_;
public:
 explicit PreparedPublicAutoNTT(const PublicAutoGadgetRows<P,Q,Bits>& rows):exponent_(rows.exponent) {
  if(!exponent_ || exponent_>=2*P::n || !(exponent_&1))throw std::invalid_argument("invalid derived automorphism exponent");
  if constexpr(DirectRNS){cache_=std::make_shared<Cache>(rows);}
  else {
  // Reuse only the generic encrypted-gadget multiplication arithmetic.
  using Matrix=CKKSCompactDDRelinKey<P,Q,Bits,16>;
  static_assert(Matrix::primary_rows==rows.rows);
  auto matrix=std::make_unique<Matrix>();
  auto digits=std::make_unique<std::array<TRLWE<P>,Matrix::bbar_rows>>();
  for(size_t row=0;row<rows.rows;++row) {
   ckks_detail::activeBaseDecomposeTRLWERows<P,Q,16,Matrix::bbar_rows>(*digits,rows.ciphertexts[row]);
   for(size_t secondary=0;secondary<Matrix::bbar_rows;++secondary)
    for(size_t c=0;c<=P::k;++c)for(size_t i=0;i<P::n;++i)matrix->data[row][secondary][c][i]=Matrix::pack((*digits)[secondary][c][i]);
  }
  cache_=makeCKKSDDRelinKeyNTT<Secondary>(*matrix);
  }
 }
 size_t byteSize()const{return cache_->byteSize();}
 template<uint32_t DataQ> void apply(TRLWE<P>& result,const TRLWE<P>& input)const {
  static_assert(DataQ<Q);constexpr auto auxiliary=Q-DataQ;
  auto transformed=std::make_unique<Polynomial<P>>();auto body=std::make_unique<Polynomial<P>>();
  Automorphism<P>(*transformed,input[0],exponent_);Automorphism<P>(*body,input[1],exponent_);
  for(auto& coefficient:*transformed)coefficient=ckks_detail::reduceToLevel<P,DataQ>(coefficient)<<auxiliary;
  auto switched=std::make_unique<TRLWE<P>>();cache_->apply(*switched,*transformed);
  for(size_t i=0;i<P::n;++i) {
   result[0][i]=ckks_detail::reduceToLevel<P,DataQ>(-ckks_detail::roundedLevelRightShift<P,Q,auxiliary>((*switched)[0][i]));
   result[1][i]=ckks_detail::reduceToLevel<P,DataQ>((*body)[i]-ckks_detail::roundedLevelRightShift<P,Q,auxiliary>((*switched)[1][i]));
  }
 }
};
}
