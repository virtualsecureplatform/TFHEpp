#pragma once
#include <ckks/narrow_factor_rns32.hpp>
#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp::factor_rescale_experiment {
// Frozen factor source from 0cb8d9b: final rescale remains serial.
template <class TargetSchedule,class SourceSchedule,class SelectorSchedule=SourceSchedule,bool DirectCRT=true>
class SerialProjectionFactorSourceRNS32 final : public CKKSNarrowFactorSource<TargetSchedule> {
    using P = typename TargetSchedule::Param;
    using Source = CKKSNarrowSelectorSource<SelectorSchedule>;
    static_assert(std::is_same_v<typename SourceSchedule::Selector,typename SelectorSchedule::Selector>);
    static_assert(std::is_same_v<P,typename SourceSchedule::Param>);
    static_assert(SourceSchedule::factor_log_q > TargetSchedule::factor_log_q);
    static constexpr auto drop = SourceSchedule::factor_log_q-TargetSchedule::factor_log_q;
    static_assert(SourceSchedule::log_delta-drop == TargetSchedule::log_delta);
    static_assert(SourceSchedule::plain_log_delta == TargetSchedule::plain_log_delta);
    static_assert(SourceSchedule::factor_key_log_q == TargetSchedule::factor_key_log_q);
    static_assert(SourceSchedule::blocks == TargetSchedule::blocks && SourceSchedule::blocks == SelectorSchedule::blocks);
    static_assert(SourceSchedule::block_width == SelectorSchedule::block_width);
    static_assert(SourceSchedule::head_slots == TargetSchedule::head_slots && SourceSchedule::head_slots == SelectorSchedule::head_slots);
    static_assert(SourceSchedule::tail_width == TargetSchedule::tail_width && SourceSchedule::tail_width == SelectorSchedule::tail_width);
    static_assert(SourceSchedule::tail_packs == TargetSchedule::tail_packs && SourceSchedule::tail_packs == SelectorSchedule::tail_packs);
    static_assert(SourceSchedule::tail_blocks_per_pack == TargetSchedule::tail_blocks_per_pack && SourceSchedule::tail_blocks_per_pack == SelectorSchedule::tail_blocks_per_pack);
    std::shared_ptr<const Source> source_;
    std::shared_ptr<const CKKSNarrowFactorSelectorCacheRNS32<SourceSchedule,SelectorSchedule>> retained_;
    bool tiled_matrix_,parallel_conjugation_,dot8_matrix_,vector_short_stages_,packed_matrix_,all_retained_rows_,joint_components_,parallel_table_ntt_,fused_components_,batch_table_encode_;
    struct FactorKeys {
        const CKKSGaloisKey<P,SourceSchedule::factor_key_log_q> &factor_galois;
        bool staged_fft_caches;
        const decltype(CKKSNarrowBootstrapKey<TargetSchedule>::automorphism_evaluator) &automorphism_evaluator;
        const decltype(CKKSNarrowBootstrapKey<SourceSchedule>::factor_conjugation_fft) &factor_conjugation_fft;
        const decltype(CKKSNarrowBootstrapKey<SourceSchedule>::factor_conjugation_projected_fft) &factor_conjugation_projected_fft;
        const decltype(CKKSNarrowBootstrapKey<TargetSchedule>::automorphism_batch_evaluator) &automorphism_batch_evaluator;
    };
public:
    explicit SerialProjectionFactorSourceRNS32(std::shared_ptr<const Source> source,
        std::shared_ptr<const CKKSNarrowFactorSelectorCacheRNS32<SourceSchedule,SelectorSchedule>> retained={},bool tiled_matrix=true,
        bool parallel_conjugation=false,bool dot8_matrix=true,bool vector_short_stages=false,bool packed_matrix=true,bool all_retained_rows=true,bool joint_components=false,bool parallel_table_ntt=false,bool fused_components=true,bool batch_table_encode=true)
        : source_(std::move(source)),retained_(std::move(retained)),tiled_matrix_(tiled_matrix),parallel_conjugation_(parallel_conjugation),dot8_matrix_(dot8_matrix),vector_short_stages_(vector_short_stages),packed_matrix_(packed_matrix),all_retained_rows_(all_retained_rows),joint_components_(joint_components),parallel_table_ntt_(parallel_table_ntt),fused_components_(fused_components),batch_table_encode_(batch_table_encode)
    {
        if (!source_ || !source_->hasRawChunks()) throw std::invalid_argument("direct RNS factor source requires raw chunks");
        if (retained_ && retained_->source()!=source_.get())
            throw std::invalid_argument("direct RNS cache belongs to another source");
    }
    std::size_t size(bool tail) const override { return source_->size(tail); }
    std::vector<typename TargetSchedule::Factor> evaluate(
        std::span<const CKKSSlotVector<P>> tables,bool tail,std::size_t group,
        const CKKSNarrowBootstrapKey<TargetSchedule> &key,
        const CKKSPreparedAutoKeyFFT<P,TargetSchedule::factor_key_log_q> *cache,
        CKKSNarrowFactorTimings *timings) const override
    {
        FactorKeys view{key.factor_galois,key.staged_fft_caches,key.automorphism_evaluator,
            key.factor_conjugation_fft,key.factor_conjugation_projected_fft,key.automorphism_batch_evaluator};
        auto factors = narrow_detail::packedFactorsBatchRNS32<SourceSchedule,DirectCRT>(tables,view,key.selector_ntt_chunk_size,
            [&](auto first,auto count,auto &raw) { source_->loadChunk(tail,group,first,count,raw); },cache,timings,
            source_->supportsRawChunkPrefetch(),retained_?retained_->group(tail,group):
                std::span<const narrow_detail::FactorSelectorRNS32<SourceSchedule>>{},tiled_matrix_,parallel_conjugation_,dot8_matrix_,vector_short_stages_,packed_matrix_,all_retained_rows_,joint_components_,parallel_table_ntt_,fused_components_,batch_table_encode_);
        std::vector<typename TargetSchedule::Factor> result(factors.size());
        for (std::size_t i = 0; i < factors.size(); ++i)
            CKKSScaleDown<P,SourceSchedule::factor_log_q,SourceSchedule::log_delta,drop>(result[i],factors[i]);
        return result;
    }
};
}
#endif
