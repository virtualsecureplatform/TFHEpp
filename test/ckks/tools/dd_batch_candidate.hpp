#pragma once
#include <ckks/prepared_dd_ntt.hpp>
#include <ckks/prepared_gadget_batch_rns32.hpp>
#include <ckks/tensor_rns32.hpp>

namespace TFHEpp::dd_batch_experiment {
// Exact existing tensor, auxiliary rounding, linear addition and data rescale;
// only the scheduling and shared gadget application change.
template<class P, std::uint32_t Q, std::uint32_t Delta, std::uint32_t Auxiliary,
         std::uint32_t Bits, std::uint32_t Chunk>
void multiply(std::span<CKKSCiphertext<P,Q-Delta,Delta>> output,
              std::span<const CKKSCiphertext<P,Q,Delta>> left,
              std::span<const CKKSCiphertext<P,Q,Delta>> right,
              const CKKSAuxiliaryDDRelinKeyNTT<P,Q,Auxiliary,Bits,20,true,Chunk>& key,
              std::size_t window=8)
{
    if(output.size()!=left.size()||left.size()!=right.size()||!window)
        throw std::invalid_argument("invalid DD multiply batch shape/window");
    if(left.empty())return;
    constexpr auto KeyQ=Q+Auxiliary;
    const auto size=left.size();
    std::vector<TRLWE3<P>> tensor(size);
    std::vector<Polynomial<P>> lifted(size);
    std::vector<TRLWE<P>> switched(size);
    const auto threads=std::min<int>(size,omp_get_max_threads());
    std::vector<std::exception_ptr> errors(size);
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
    for(std::size_t j=0;j<size;++j){
        try{
            CKKSTensorProductRNS32<P,Q>(tensor[j],left[j].ct,right[j].ct);
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
                output[j].ct[c][i]=ckks_detail::roundedLevelRightShift<P,Q,Delta>(value);
            }
}
}
