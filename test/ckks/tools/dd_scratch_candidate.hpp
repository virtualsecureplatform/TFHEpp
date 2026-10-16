#pragma once
#include <ckks/prepared_dd_batch_rns32.hpp>
#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp::dd_scratch_experiment {
// Frozen previous production helper for exact matched allocation comparisons.
// Exact existing tensor, auxiliary rounding, linear addition and data rescale;
// only the scheduling and shared gadget application change.
template<class P, std::uint32_t Q, std::uint32_t Delta, std::uint32_t Auxiliary,
         std::uint32_t Bits, std::uint32_t Secondary, std::uint32_t Chunk>
inline void multiplyReallocated(
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
    for(std::size_t offset=0;offset<products.size();offset+=window){
        const auto size=std::min(window,products.size()-offset);
        std::vector<TRLWE3<P>> tensor(size);
        std::vector<Polynomial<P>> lifted(size);
        std::vector<TRLWE<P>> switched(size);
#ifdef _OPENMP
        const auto threads=std::min<int>(size,omp_get_max_threads());
#else
        [[maybe_unused]] const auto threads=1;
#endif
        std::vector<std::exception_ptr> errors(size);
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for(std::size_t j=0;j<size;++j){
            try{
                CKKSTensorProductRNS32<P,Q>(tensor[j],products[offset+j].left->ct,products[offset+j].right->ct);
                for(std::size_t i=0;i<P::n;++i)
                    lifted[j][i]=ckks_detail::reduceToLevel<P,KeyQ>(tensor[j][2][i]<<Auxiliary);
            }catch(...){errors[j]=std::current_exception();}
        }
        for(const auto&error:errors)if(error)std::rethrow_exception(error);
        CKKSApplyGadgetBatchRNS32<4>(switched,lifted,*key.cache,window);
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
}
#endif
