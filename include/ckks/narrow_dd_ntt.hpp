#pragma once
#include "narrow_bootstrapping.hpp"
#include "prepared_dd_ntt_chain.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
// Reuse a saved public wide-primary key without changing its encryption.
// The auxiliary lift protects switching precision; data rescale count stays
// equal to the existing before-rescale product tree's count.
template <std::uint32_t Auxiliary = 100, std::uint32_t Secondary = 32,
          class Schedule, std::uint32_t MasterQ, std::uint32_t Primary,
          std::uint32_t MasterSecondary>
inline void CKKSNarrowUseAuxiliaryNTTDD(CKKSNarrowBootstrapKey<Schedule> &key,
    const CKKSCompactDDRelinKey<typename Schedule::Param, MasterQ, Primary,
        MasterSecondary> &master)
{
    static_assert(Schedule::relinearize_before_rescale);
    using Chain=CKKSPreparedAuxiliaryDDRelinKeyChainNTT<typename Schedule::Param,
        Schedule::factor_log_q, Schedule::log_delta, Schedule::tree_depth,
        Auxiliary, Primary, Secondary>;
    auto prepared=std::make_shared<Chain>(master);
    key.product_tree_evaluator=[prepared](typename Schedule::TreeResult &result,
        std::vector<typename Schedule::Factor> factors,
        const typename Schedule::RelinKeys &raw) {
        narrow_detail::productTree<Schedule>(result,std::move(factors),raw,prepared.get());
    };
    key.relin_fft.reset();
    key.relin_projected_fft.reset();
}
} // namespace TFHEpp
#endif
