#pragma once
#include "prepared_dd_ntt.hpp"
#include "prepared_gadget_batch_rns32.hpp"
#include "tensor_rns32.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
// Borrow typed operands from forests without gathering or copying ciphertexts.
// Outputs must be distinct and must not overlap any input in the batch.
template<class P,std::uint32_t Q,std::uint32_t Delta>
struct CKKSRNS32Product {
    CKKSCiphertext<P,Q-Delta,Delta>* output;
    const CKKSCiphertext<P,Q,Delta>* left;
    const CKKSCiphertext<P,Q,Delta>* right;
};
// Exact existing tensor, auxiliary rounding, linear addition and data rescale;
// only the scheduling and shared gadget application change.
template<class P, std::uint32_t Q, std::uint32_t Delta, std::uint32_t Auxiliary,
         std::uint32_t Bits, std::uint32_t Secondary, std::uint32_t Chunk>
inline void CKKSMultRelinThenRescaleBatchRNS32(
    std::type_identity_t<std::span<const CKKSRNS32Product<P,Q,Delta>>> products,
    const CKKSAuxiliaryDDRelinKeyNTT<P,Q,Auxiliary,Bits,Secondary,true,Chunk>& key,
    std::size_t window=8)
{
    static_assert(P::k==1 && Delta>0 && 2*Delta<Q);
    if(!window)throw std::invalid_argument("zero RNS32 product window");
    for(const auto&item:products)
        if(!item.output||!item.left||!item.right)
            throw std::invalid_argument("missing RNS32 product operand/output");
    constexpr auto KeyQ=Q+Auxiliary;
    // Every polynomial is overwritten by the next phase. Retain one bounded
    // workspace through all windows instead of zero-initializing it repeatedly.
    const auto capacity=std::min(window,products.size());
    std::vector<TRLWE3<P>> tensor(capacity);
    std::vector<Polynomial<P>> lifted(capacity);
    std::vector<TRLWE<P>> switched(capacity);
    std::vector<std::exception_ptr> errors(capacity);
    for(std::size_t offset=0;offset<products.size();offset+=window){
        const auto size=std::min(window,products.size()-offset);
#ifdef _OPENMP
        const auto threads=std::min<int>(size,omp_get_max_threads());
#else
        [[maybe_unused]] const auto threads=1;
#endif
        std::fill(errors.begin(),errors.end(),std::exception_ptr{});
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for(std::size_t j=0;j<size;++j){
            try{
                CKKSTensorProductRNS32<P,Q>(tensor[j],products[offset+j].left->ct,products[offset+j].right->ct);
                for(std::size_t i=0;i<P::n;++i)
                    lifted[j][i]=ckks_detail::reduceToLevel<P,KeyQ>(tensor[j][2][i]<<Auxiliary);
            }catch(...){errors[j]=std::current_exception();}
        }
        for(const auto&error:errors)if(error)std::rethrow_exception(error);
        CKKSApplyGadgetBatchRNS32<4>(std::span(switched).first(size),std::span<const Polynomial<P>>(lifted).first(size),*key.cache,window);
#pragma omp parallel for collapse(3) schedule(dynamic,1)
        for(std::size_t j=0;j<size;++j)for(std::size_t c=0;c<2;++c)
            for(std::size_t begin=0;begin<P::n;begin+=256)
                for(std::size_t i=begin;i<std::min(begin+256,std::size_t(P::n));++i){
                    auto value=ckks_detail::roundedLevelRightShift<P,KeyQ,Auxiliary>(switched[j][c][i]);
                    value=ckks_detail::reduceToLevel<P,Q>(value+tensor[j][c][i]);
                    products[offset+j].output->ct[c][i]=ckks_detail::roundedLevelRightShift<P,Q,Delta>(value);
                }
    }
}
} // namespace TFHEpp
#endif
