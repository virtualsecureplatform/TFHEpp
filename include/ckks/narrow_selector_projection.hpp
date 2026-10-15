#pragma once

#include "narrow_bootstrapping.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Public coefficient reduction preserves the selector scale. No per-selector
// rounding is introduced before the weighted sum, and no keys are generated.
template<class TargetSchedule,class SourceSchedule>
class CKKSNarrowReducedSelectorSource final:public CKKSNarrowSelectorSource<TargetSchedule>{
    using P=typename TargetSchedule::Param;
    using Source=CKKSNarrowSelectorSource<SourceSchedule>;
    using Selector=typename TargetSchedule::Selector;
    static_assert(std::is_same_v<P,typename SourceSchedule::Param>);
    static_assert(TargetSchedule::selector_log_q<SourceSchedule::selector_log_q);
    static_assert(TargetSchedule::log_delta==SourceSchedule::log_delta);
    static_assert(TargetSchedule::input_log_q==SourceSchedule::input_log_q);
    static_assert(TargetSchedule::blocks==SourceSchedule::blocks && TargetSchedule::block_width==SourceSchedule::block_width);
    static_assert(TargetSchedule::tail_packs==SourceSchedule::tail_packs);
    std::shared_ptr<const Source> source_;
    static void convert(std::vector<Selector>&out,const std::vector<typename SourceSchedule::Selector>&raw){
        out.resize(raw.size());
        for(std::size_t i=0;i<raw.size();++i)CKKSLevelReduce<P>(out[i],raw[i]);
    }
public:
    explicit CKKSNarrowReducedSelectorSource(std::shared_ptr<const Source> source):source_(std::move(source)){
        if(!source_)throw std::invalid_argument("missing source for reduced selectors");
    }
    std::size_t size(bool tail)const override{return source_->size(tail);}
    bool hasRawChunks()const override{return source_->hasRawChunks();}
    bool supportsRawChunkPrefetch()const override{return source_->supportsRawChunkPrefetch();}
    void load(bool tail,std::size_t group,std::vector<Selector>&out)const override{
        std::vector<typename SourceSchedule::Selector> raw;source_->load(tail,group,raw);convert(out,raw);
    }
    void loadChunk(bool tail,std::size_t group,std::size_t first,std::size_t count,std::vector<Selector>&out)const override{
        std::vector<typename SourceSchedule::Selector> raw;source_->loadChunk(tail,group,first,count,raw);convert(out,raw);
    }
};

// Reuse higher-scale selectors without reconstructing a packed blind-rotation
// key. Conversion is public rounded rescaling, performed one requested chunk
// at a time. The selector geometry and encoded plaintext are unchanged.
// Prepared spectra belong to the original coefficients and cannot be reused.
template <class TargetSchedule, class SourceSchedule>
class CKKSNarrowRescaledSelectorSource final
    : public CKKSNarrowSelectorSource<TargetSchedule> {
    using P = typename TargetSchedule::Param;
    using Source = CKKSNarrowSelectorSource<SourceSchedule>;
    using Selector = typename TargetSchedule::Selector;
    static_assert(std::is_same_v<P, typename SourceSchedule::Param>);
    static_assert(SourceSchedule::selector_log_q > TargetSchedule::selector_log_q);
    static constexpr auto drop = SourceSchedule::selector_log_q - TargetSchedule::selector_log_q;
    static_assert(SourceSchedule::log_delta >= drop);
    static_assert(SourceSchedule::log_delta - drop == TargetSchedule::log_delta);
    static_assert(SourceSchedule::plain_log_delta == TargetSchedule::plain_log_delta);
    static_assert(SourceSchedule::input_log_q == TargetSchedule::input_log_q);
    static_assert(SourceSchedule::blocks == TargetSchedule::blocks);
    static_assert(SourceSchedule::block_width == TargetSchedule::block_width);
    static_assert(SourceSchedule::tail_packs == TargetSchedule::tail_packs);
    std::shared_ptr<const Source> source_;

    static void convert(std::vector<Selector> &result,
                        const std::vector<typename SourceSchedule::Selector> &raw)
    {
        result.resize(raw.size());
        for (std::size_t i = 0; i < raw.size(); ++i)
            CKKSScaleDown<P, SourceSchedule::selector_log_q,
                SourceSchedule::log_delta, drop>(result[i], raw[i]);
    }

public:
    explicit CKKSNarrowRescaledSelectorSource(std::shared_ptr<const Source> source)
        : source_(std::move(source))
    {
        if (!source_) throw std::invalid_argument("missing source for rescaled narrow selectors");
    }
    std::size_t size(bool tail) const override { return source_->size(tail); }
    bool hasRawChunks() const override { return source_->hasRawChunks(); }
    bool supportsRawChunkPrefetch() const override
    { return source_->supportsRawChunkPrefetch(); }
    void load(bool tail, std::size_t group, std::vector<Selector> &result) const override
    {
        std::vector<typename SourceSchedule::Selector> raw;
        source_->load(tail, group, raw);
        convert(result, raw);
    }
    void loadChunk(bool tail, std::size_t group, std::size_t first, std::size_t count,
                   std::vector<Selector> &result) const override
    {
        std::vector<typename SourceSchedule::Selector> raw;
        source_->loadChunk(tail, group, first, count, raw);
        convert(result, raw);
    }
};

// Compute with the retained high-scale selectors, then round once per factor.
// Rounding each selector before the weighted sum amplifies its rounding noise.
// The source and target must use the same conjugation key modulus; product-tree
// and tail-alignment evaluation continue at the target scale.
template <class TargetSchedule, class SourceSchedule, class SelectorSchedule = SourceSchedule>
class CKKSNarrowDeferredRescaleFactorSource final
    : public CKKSNarrowFactorSource<TargetSchedule> {
    using P = typename TargetSchedule::Param;
    // Selector archives describe their original reconstruction parameters.
    // Table encoding precision may change without changing selector bytes.
    using Source = CKKSNarrowSelectorSource<SelectorSchedule>;
    static_assert(std::is_same_v<typename SelectorSchedule::Selector,
                                 typename SourceSchedule::Selector>);
    static_assert(SelectorSchedule::input_log_q == SourceSchedule::input_log_q);
    static_assert(SelectorSchedule::blocks == SourceSchedule::blocks);
    static_assert(SelectorSchedule::block_width == SourceSchedule::block_width);
    static_assert(SelectorSchedule::head_slots == SourceSchedule::head_slots);
    static_assert(SelectorSchedule::tail_width == SourceSchedule::tail_width);
    static_assert(SelectorSchedule::tail_packs == SourceSchedule::tail_packs);
    static_assert(SelectorSchedule::tail_blocks_per_pack == SourceSchedule::tail_blocks_per_pack);
    static_assert(std::is_same_v<P, typename SourceSchedule::Param>);
    static_assert(SourceSchedule::factor_log_q > TargetSchedule::factor_log_q);
    static constexpr auto drop = SourceSchedule::factor_log_q - TargetSchedule::factor_log_q;
    static_assert(SourceSchedule::log_delta >= drop);
    static_assert(SourceSchedule::log_delta - drop == TargetSchedule::log_delta);
    static_assert(SourceSchedule::plain_log_delta == TargetSchedule::plain_log_delta);
    static_assert(SourceSchedule::input_log_q == TargetSchedule::input_log_q);
    static_assert(SourceSchedule::blocks == TargetSchedule::blocks);
    static_assert(SourceSchedule::block_width == TargetSchedule::block_width);
    static_assert(SourceSchedule::head_slots == TargetSchedule::head_slots);
    static_assert(SourceSchedule::tail_width == TargetSchedule::tail_width);
    static_assert(SourceSchedule::tail_packs == TargetSchedule::tail_packs);
    static_assert(SourceSchedule::tail_blocks_per_pack == TargetSchedule::tail_blocks_per_pack);
    static_assert(SourceSchedule::factor_key_log_q == TargetSchedule::factor_key_log_q);
    std::shared_ptr<const Source> source_;
    struct FactorKeys {
        const CKKSGaloisKey<P, SourceSchedule::factor_key_log_q> &factor_galois;
        bool staged_fft_caches;
        const decltype(CKKSNarrowBootstrapKey<TargetSchedule>::automorphism_evaluator) &automorphism_evaluator;
        const decltype(CKKSNarrowBootstrapKey<SourceSchedule>::factor_conjugation_fft) &factor_conjugation_fft;
        const decltype(CKKSNarrowBootstrapKey<SourceSchedule>::factor_conjugation_projected_fft) &factor_conjugation_projected_fft;
    };
public:
    explicit CKKSNarrowDeferredRescaleFactorSource(std::shared_ptr<const Source> source)
        : source_(std::move(source))
    {
        if (!source_ || !source_->hasRawChunks())
            throw std::invalid_argument("deferred factor rescaling requires raw selector chunks");
    }
    std::size_t size(bool tail) const override { return source_->size(tail); }
    std::vector<typename TargetSchedule::Factor> evaluate(
        std::span<const CKKSSlotVector<P>> tables, bool tail, std::size_t group,
        const CKKSNarrowBootstrapKey<TargetSchedule> &key,
        const CKKSPreparedAutoKeyFFT<P, TargetSchedule::factor_key_log_q> *cache,
        CKKSNarrowFactorTimings *timings) const override
    {
        if (key.selector_ntt_chunk_size == 0 || key.cache_mode == CKKSNarrowCacheMode::Reference)
            throw std::invalid_argument("deferred factor rescaling requires chunked matrix evaluation");
        FactorKeys view{key.factor_galois, key.staged_fft_caches, key.automorphism_evaluator,
            key.factor_conjugation_fft, key.factor_conjugation_projected_fft};
        const std::function<void(std::size_t, std::size_t,
            std::vector<typename SourceSchedule::Selector> &)> loader =
            [&](auto first, auto count, auto &chunk) {
                source_->loadChunk(tail, group, first, count, chunk);
            };
        auto factors = narrow_detail::packedFactorsBatch<SourceSchedule>(tables,
            {}, {}, 0, view, key.selector_ntt_chunk_size, cache, timings, {}, loader,
            source_->supportsRawChunkPrefetch());
        std::vector<typename TargetSchedule::Factor> result(factors.size());
        for (std::size_t i = 0; i < factors.size(); ++i)
            CKKSScaleDown<P, SourceSchedule::factor_log_q,
                SourceSchedule::log_delta, drop>(result[i], factors[i]);
        return result;
    }
};

// Opt-in correction for selectors produced by ExtractSelectorTerms. Their
// public slot masks have a small encoding gain error which repeats across the
// product tree. Correct that known gain in the existing transparent beta leaf;
// this requires neither secret information nor another multiplication level.
// Do not apply this wrapper to directly encrypted selectors or different masks.
template <class Schedule, class SelectorSchedule>
class CKKSNarrowPublicMaskNormalizedFactorSource final
    : public CKKSNarrowFactorSource<Schedule> {
    using P = typename Schedule::Param;
    static_assert(std::is_same_v<P, typename SelectorSchedule::Param>);
    static_assert(Schedule::blocks == SelectorSchedule::blocks);
    static_assert(Schedule::head_slots == SelectorSchedule::head_slots);
    static_assert(Schedule::tail_width == SelectorSchedule::tail_width);
    static_assert(Schedule::tail_blocks_per_pack == SelectorSchedule::tail_blocks_per_pack);
    static constexpr auto mask_log_q = SelectorSchedule::selector_log_q + SelectorSchedule::plain_log_delta;
    static_assert(mask_log_q <= std::numeric_limits<typename P::T>::digits);
    std::shared_ptr<const CKKSNarrowFactorSource<Schedule>> source_;
    CKKSSlotVector<P> correction_;
public:
    explicit CKKSNarrowPublicMaskNormalizedFactorSource(
        std::shared_ptr<const CKKSNarrowFactorSource<Schedule>> source)
        : source_(std::move(source))
    {
        if (!source_) throw std::invalid_argument("missing factor source for public mask normalization");
        std::vector<CKKSSlotVector<P>> masks(1 + Schedule::tail_blocks_per_pack);
        auto poly = std::make_unique<Polynomial<P>>();
        auto slots = std::make_unique<CKKSSlotVector<P>>();
        for (std::size_t j = 0; j < masks.size(); ++j) {
            slots->fill({});
            const auto begin = j ? (2 * (j - 1) + 1) * Schedule::tail_width : 0;
            const auto count = j ? Schedule::tail_width : Schedule::head_slots;
            std::fill_n(slots->begin() + begin, count, std::complex<double>{1, 0});
            ckksSlotEncode<P, mask_log_q, SelectorSchedule::plain_log_delta>(*poly, *slots);
            ckksSlotDecode<P, mask_log_q, SelectorSchedule::plain_log_delta>(masks[j], *poly);
        }
        for (std::size_t i = 0; i < Schedule::slots; ++i) {
            std::complex<double> gain{1, 0};
            for (std::size_t block = 0; block < Schedule::blocks; ++block) {
                if (i < Schedule::head_slots) gain *= masks[0][i];
                else {
                    const auto ell = block % Schedule::tail_blocks_per_pack;
                    const auto rotated = (i + 2 * (ell + 1) * Schedule::tail_width) % Schedule::slots;
                    gain *= masks[1 + ell][rotated];
                }
            }
            if (!std::isfinite(gain.real()) || !std::isfinite(gain.imag()) ||
                std::abs(gain) < 0.5 || std::abs(gain) > 2.0)
                throw std::invalid_argument("invalid public reconstruction mask gain");
            correction_[i] = 1.0 / gain;
        }
    }
    std::size_t size(bool tail) const override { return source_->size(tail); }
    void adjustBeta(CKKSSlotVector<P> &head, CKKSSlotVector<P> &tail) const override
    {
        source_->adjustBeta(head, tail);
        for (std::size_t i = 0; i < Schedule::slots; ++i) {
            head[i] *= correction_[i];
            tail[i] *= correction_[i];
        }
    }
    std::vector<typename Schedule::Factor> evaluate(
        std::span<const CKKSSlotVector<P>> tables, bool tail, std::size_t group,
        const CKKSNarrowBootstrapKey<Schedule> &key,
        const CKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q> *cache,
        CKKSNarrowFactorTimings *timings) const override
    { return source_->evaluate(tables, tail, group, key, cache, timings); }
};

}  // namespace TFHEpp
#endif
