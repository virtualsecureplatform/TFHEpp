#pragma once
#include "public_auto_ntt.hpp"

namespace TFHEpp {
// Direct evaluation-key setup from an existing secret. Replays load the saved
// public rows instead of calling this generator. The secret is never archived.
template<class P,uint32_t Q,uint32_t Bits>
inline void CKKSGenerateNativeAutoGadgetRows(PublicAutoGadgetRows<P,Q,Bits>&result,
    uint32_t exponent,const Key<P>&key,CKKSNoise noise={P::α,0})
{
 if(!exponent || exponent>=2*P::n || !(exponent&1))throw std::invalid_argument("invalid native automorphism exponent");
 auto secret=std::make_unique<Polynomial<P>>(),transformed=std::make_unique<Polynomial<P>>();
 for(std::size_t i=0;i<P::n;++i)(*secret)[i]=key[i];
 Automorphism<P>(*transformed,*secret,exponent);result.exponent=exponent;
 std::exception_ptr error;
#pragma omp parallel for schedule(static)
 for(std::size_t row=0;row<result.rows;++row){
  try{
   auto plaintext=std::make_unique<Polynomial<P>>();
   const auto power=(result.rows-row-1)*Bits;
   for(std::size_t i=0;i<P::n;++i)(*plaintext)[i]=ckks_detail::reduceToLevel<P,Q>((*transformed)[i]<<power);
   ckks_detail::encryptPolynomialAtLevel<P,Q>(result.ciphertexts[row],*plaintext,key,noise);
  }catch(...){
#pragma omp critical(ckks_native_auto_generation_error)
   {if(!error)error=std::current_exception();}
  }
 }
 if(error)std::rethrow_exception(error);
}

template<class P,uint32_t Q,uint32_t Bits>
inline auto CKKSNativeAutoFileMetadata(uint32_t exponent,double noise_stddev)
{
 return std::array<uint64_t,8>{0x4e4155544f303031ULL,P::n,P::T::digits,Q,Bits,
  PublicAutoGadgetRows<P,Q,Bits>::rows,exponent,std::bit_cast<uint64_t>(noise_stddev)};
}

// Fixed-size profile and origin fields precede public ciphertexts. No secret
// values or unbounded archive lengths occur in this format.
template<class P,uint32_t Q,uint32_t Bits>
inline void CKKSWriteNativeAutoRows(std::ostream&file,const PublicAutoGadgetRows<P,Q,Bits>&rows,
    const std::string&origin,double noise_stddev)
{
 if(origin.size()!=64)throw std::invalid_argument("native auto origin must be a SHA256 hex digest");
 std::array<char,64> source{};std::copy(origin.begin(),origin.end(),source.begin());
 cereal::PortableBinaryOutputArchive archive(file);
 archive(CKKSNativeAutoFileMetadata<P,Q,Bits>(rows.exponent,noise_stddev),source,rows.ciphertexts);
 if(!file)throw std::runtime_error("native auto archive write failed");
}

template<class P,uint32_t Q,uint32_t Bits>
inline auto CKKSReadNativeAutoRows(std::istream&file,uint32_t exponent,
    const std::string&origin,double noise_stddev)
{
 if(origin.size()!=64)throw std::invalid_argument("native auto origin must be a SHA256 hex digest");
 auto rows=std::make_unique<PublicAutoGadgetRows<P,Q,Bits>>();rows->exponent=exponent;
 cereal::PortableBinaryInputArchive archive(file);
 std::array<uint64_t,8> metadata{};std::array<char,64> source{};
 archive(metadata,source);
 if(metadata!=CKKSNativeAutoFileMetadata<P,Q,Bits>(exponent,noise_stddev) || std::string(source.begin(),source.end())!=origin)
  throw std::runtime_error("native auto archive profile/origin mismatch");
 archive(rows->ciphertexts);
 if(file.peek()!=std::char_traits<char>::eof())throw std::runtime_error("native auto archive has trailing bytes");
 return rows;
}
}
