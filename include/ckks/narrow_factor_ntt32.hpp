#pragma once

#include "narrow_selector_projection.hpp"
#include <modular_ntt32.hpp>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
namespace narrow_detail {

template <class P>
inline const auto &factorPlans32()
{
    static const std::array<modular_ntt32::NegacyclicNTTPlan,4> plans{
        modular_ntt32::NegacyclicNTTPlan(P::n,modular_ntt32::primes[0]),
        modular_ntt32::NegacyclicNTTPlan(P::n,modular_ntt32::primes[1]),
        modular_ntt32::NegacyclicNTTPlan(P::n,modular_ntt32::primes[2]),
        modular_ntt32::NegacyclicNTTPlan(P::n,modular_ntt32::primes[3])};
    return plans;
}

template <class Schedule>
struct FactorSelector32 {
    using P = typename Schedule::Param;
    static constexpr std::size_t digit_bits = 52;
    static constexpr std::size_t digits = (Schedule::selector_log_q+digit_bits-1)/digit_bits;
    static constexpr std::uint64_t mask = (std::uint64_t{1} << digit_bits)-1;
    static constexpr unsigned __int128 bound =
        ((static_cast<unsigned __int128>(1) << Schedule::plain_log_delta)+1)*mask*P::n*Schedule::block_width;
    static_assert(bound < modular_ntt32::crt_product/2,"selector convolution exceeds four-prime CRT range");
    std::array<std::vector<std::uint32_t>,4> spectra;
};

template <class Schedule>
inline void prepareFactorSelectors32(std::vector<FactorSelector32<Schedule>> &prepared,
                                    std::span<const typename Schedule::Selector> selectors)
{
    using P = typename Schedule::Param;
    using Cache = FactorSelector32<Schedule>;
    const auto &plans = factorPlans32<P>();
    prepared.resize(selectors.size());
#pragma omp parallel for schedule(static)
    for (std::size_t j = 0; j < selectors.size(); ++j)
        for (std::size_t prime = 0; prime < 4; ++prime) {
            auto &spectrum = prepared[j].spectra[prime];
            spectrum.resize((P::k+1)*Cache::digits*P::n);
            for (std::size_t component = 0; component <= P::k; ++component)
                for (std::size_t digit = 0; digit < Cache::digits; ++digit) {
                    auto row = std::span(spectrum).subspan((component*Cache::digits+digit)*P::n,P::n);
                    for (std::size_t i = 0; i < P::n; ++i)
                        row[i] = plans[prime].arithmetic().encode(
                            selectorDigit<P,Schedule::selector_log_q,Cache::digit_bits>(selectors[j].ct[component][i],digit));
                    plans[prime].forward(row);
                }
        }
}

// Same public tables, selector coefficients, torus recomposition and final
// rounding as packedFactorsBatch; only the exact convolution basis changes.
template <class Schedule, class FactorKeys, class Loader>
inline std::vector<typename Schedule::Factor> packedFactorsBatch32(
    std::span<const CKKSSlotVector<typename Schedule::Param>> tables,
    const FactorKeys &key, std::size_t chunk_size, Loader load,
    const CKKSPreparedAutoKeyFFT<typename Schedule::Param,Schedule::factor_key_log_q> *factor_fft,
    CKKSNarrowFactorTimings *timings, bool prefetch_raw = false)
{
    using P = typename Schedule::Param;
    using Cache = FactorSelector32<Schedule>;
    using Clock = std::chrono::steady_clock;
    if (!chunk_size) throw std::invalid_argument("small-prime selector chunk must be nonzero");
    const auto columns = tables.size();
    if (!columns) return {};
    const auto &plans = factorPlans32<P>();
    const modular_ntt32::FourPrimeCRT crt;
    auto start = Clock::now();
    std::array<std::array<std::vector<std::uint32_t>,4>,2> base;
    for (auto &signed_tables : base)
        for (auto &spectrum : signed_tables) spectrum.resize(columns*P::n);
    for (std::size_t t = 0; t < 2; ++t)
        for (std::size_t column = 0; column < columns; ++column) {
            CKKSSlotVector<P> transformed;
            for (std::size_t i = 0; i < Schedule::slots; ++i) {
                const auto z = tables[column][i];
                transformed[i] = t == 0 ? (z-std::complex<double>(0,1)*std::conj(z))/2.0
                                       : (std::conj(z)-std::complex<double>(0,1)*z)/2.0;
            }
            auto polynomial = std::make_unique<Polynomial<P>>();
            ckksSlotEncode<P,Schedule::selector_log_q,Schedule::plain_log_delta>(*polynomial,transformed);
            for (std::size_t prime = 0; prime < 4; ++prime) {
                const auto modulus = modular_ntt32::primes[prime].value;
                std::vector<std::uint32_t> spectrum(P::n);
                for (std::size_t i = 0; i < P::n; ++i) {
                    const auto [negative,magnitude] = ckks_detail::smallSignedMagnitude<P,Schedule::selector_log_q>((*polynomial)[i]);
                    if (static_cast<unsigned __int128>(magnitude)*Cache::mask*P::n*Schedule::block_width >= modular_ntt32::crt_product/2)
                        throw std::invalid_argument("small-prime table exceeds exact CRT range");
                    const auto value = magnitude % modulus;
                    spectrum[i] = negative && value ? modulus-value : value;
                }
                plans[prime].forward(spectrum);
                for (std::size_t i = 0; i < P::n; ++i) base[t][prime][i*columns+column] = spectrum[i];
            }
        }
    if (timings) timings->table_ntt_seconds += std::chrono::duration<double>(Clock::now()-start).count();
    std::vector<std::size_t> exponents(Schedule::block_width);
    for (std::size_t y = 0; y < exponents.size(); ++y)
        exponents[y] = modular_ntt::power(5,(Schedule::slots-y)%Schedule::slots,2*P::n);
    using Sums = std::array<std::vector<std::uint32_t>,4>;
    std::vector<Sums> sums(2*(P::k+1)*Cache::digits);
    auto index = [](std::size_t t,std::size_t component,std::size_t digit) {
        return (t*(P::k+1)+component)*Cache::digits+digit;
    };
    for (auto &row : sums) for (auto &spectrum : row) spectrum.resize(columns*P::n);
    std::future<std::vector<typename Schedule::Selector>> pending;
    const auto begin_read = [&](std::size_t first) {
        const auto count = std::min(chunk_size,Schedule::block_width-first);
        pending = std::async(std::launch::async,[&,first,count] {
            std::vector<typename Schedule::Selector> raw;
            load(first,count,raw); return raw;
        });
    };
    if (prefetch_raw) begin_read(0);
    for (std::size_t first_y = 0; first_y < Schedule::block_width;) {
        const auto count = std::min(chunk_size,Schedule::block_width-first_y);
        std::vector<typename Schedule::Selector> raw;
        start = Clock::now();
        if (prefetch_raw) raw = pending.get();
        else load(first_y,count,raw);
        if (raw.size() != count) throw std::invalid_argument("small-prime selector chunk count mismatch");
        if (timings) timings->selector_load_seconds += std::chrono::duration<double>(Clock::now()-start).count();
        start = Clock::now();
        std::vector<Cache> cache; prepareFactorSelectors32<Schedule>(cache,raw);
        std::vector<typename Schedule::Selector>{}.swap(raw);
        if (prefetch_raw && first_y+count < Schedule::block_width) begin_read(first_y+count);
        if (timings) timings->selector_ntt_seconds += std::chrono::duration<double>(Clock::now()-start).count();
        start = Clock::now();
#pragma omp parallel for collapse(3) schedule(static)
        for (std::size_t prime = 0; prime < 4; ++prime)
            for (std::size_t component = 0; component <= P::k; ++component)
                for (std::size_t digit = 0; digit < Cache::digits; ++digit) {
                    const auto &arithmetic = plans[prime].arithmetic();
                    auto *first = sums[index(0,component,digit)][prime].data();
                    auto *second = sums[index(1,component,digit)][prime].data();
                    for (std::size_t begin = 0; begin < P::n; begin += 256) {
                        const auto end = std::min(begin+256,std::size_t(P::n));
                        for (std::size_t local = 0; local < count; ++local) {
                            const auto *selector = cache[local].spectra[prime].data()+(component*Cache::digits+digit)*P::n;
                            for (std::size_t i = begin; i < end; ++i) {
                                const auto position = narrowNTTPermutationIndex<P::n>(exponents[first_y+local],i)*columns;
                                arithmetic.addPair(first+i*columns,second+i*columns,
                                    base[0][prime].data()+position,base[1][prime].data()+position,columns,selector[i]);
                            }
                        }
                        for (std::size_t i = begin*columns; i < end*columns; ++i) {
                            first[i] = arithmetic.canonical(first[i]); second[i] = arithmetic.canonical(second[i]);
                        }
                    }
                }
        if (timings) timings->matrix_seconds += std::chrono::duration<double>(Clock::now()-start).count();
        first_y += count;
    }
    for (auto &signed_tables : base) for (auto &spectrum : signed_tables) std::vector<std::uint32_t>{}.swap(spectrum);
    start = Clock::now();
    std::vector<std::array<TRLWE<P>,2>> accum(columns);
#pragma omp parallel for collapse(3) schedule(static)
    for (std::size_t t = 0; t < 2; ++t)
        for (std::size_t component = 0; component <= P::k; ++component)
            for (std::size_t column = 0; column < columns; ++column) {
                std::array<std::vector<std::uint32_t>,4> spectra;
                for (auto &spectrum : spectra) spectrum.resize(P::n);
                for (std::size_t digit = 0; digit < Cache::digits; ++digit) {
                    for (std::size_t prime = 0; prime < 4; ++prime) {
                        const auto &sum = sums[index(t,component,digit)][prime];
                        for (std::size_t i = 0; i < P::n; ++i) spectra[prime][i] = sum[i*columns+column];
                        plans[prime].inverse(spectra[prime]);
                    }
                    for (std::size_t i = 0; i < P::n; ++i)
                        accum[column][t][component][i] += ckks_detail::signedToLevel<P,Schedule::selector_log_q>(
                            crt.reconstructSigned({spectra[0][i],spectra[1][i],spectra[2][i],spectra[3][i]})) << (digit*Cache::digit_bits);
                }
            }
    sums.clear();
    if (timings) timings->inverse_crt_seconds += std::chrono::duration<double>(Clock::now()-start).count();
    start = Clock::now();
    std::vector<typename Schedule::Factor> result(columns);
    for (std::size_t column = 0; column < columns; ++column) {
        auto second = std::make_unique<typename Schedule::Factor>();
        auto conjugated = std::make_unique<typename Schedule::Factor>();
        for (std::size_t component = 0; component <= P::k; ++component)
            for (std::size_t i = 0; i < P::n; ++i) {
                result[column].ct[component][i] = ckks_detail::roundedLevelRightShift<P,Schedule::selector_log_q,Schedule::plain_log_delta>(accum[column][0][component][i]);
                second->ct[component][i] = ckks_detail::roundedLevelRightShift<P,Schedule::selector_log_q,Schedule::plain_log_delta>(accum[column][1][component][i]);
            }
        conjugate<Schedule,true>(conjugated->ct,second->ct,key,factor_fft);
        CKKSAddInPlace<P>(result[column],*conjugated);
    }
    if (timings) timings->conjugation_seconds += std::chrono::duration<double>(Clock::now()-start).count();
    return result;
}
} // namespace narrow_detail

template <class TargetSchedule,class SourceSchedule,class SelectorSchedule=SourceSchedule>
class CKKSNarrowDeferredRescaleFactorSource32 final : public CKKSNarrowFactorSource<TargetSchedule> {
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
    struct FactorKeys {
        const CKKSGaloisKey<P,SourceSchedule::factor_key_log_q> &factor_galois;
        bool staged_fft_caches;
        const decltype(CKKSNarrowBootstrapKey<TargetSchedule>::automorphism_evaluator) &automorphism_evaluator;
        const decltype(CKKSNarrowBootstrapKey<SourceSchedule>::factor_conjugation_fft) &factor_conjugation_fft;
        const decltype(CKKSNarrowBootstrapKey<SourceSchedule>::factor_conjugation_projected_fft) &factor_conjugation_projected_fft;
    };
public:
    explicit CKKSNarrowDeferredRescaleFactorSource32(std::shared_ptr<const Source> source) : source_(std::move(source))
    { if (!source_ || !source_->hasRawChunks()) throw std::invalid_argument("small-prime factor source requires raw chunks"); }
    std::size_t size(bool tail) const override { return source_->size(tail); }
    std::vector<typename TargetSchedule::Factor> evaluate(
        std::span<const CKKSSlotVector<P>> tables,bool tail,std::size_t group,
        const CKKSNarrowBootstrapKey<TargetSchedule> &key,
        const CKKSPreparedAutoKeyFFT<P,TargetSchedule::factor_key_log_q> *cache,
        CKKSNarrowFactorTimings *timings) const override
    {
        FactorKeys view{key.factor_galois,key.staged_fft_caches,key.automorphism_evaluator,
            key.factor_conjugation_fft,key.factor_conjugation_projected_fft};
        auto factors = narrow_detail::packedFactorsBatch32<SourceSchedule>(tables,view,key.selector_ntt_chunk_size,
            [&](auto first,auto count,auto &raw) { source_->loadChunk(tail,group,first,count,raw); },cache,timings,
            source_->supportsRawChunkPrefetch());
        std::vector<typename TargetSchedule::Factor> result(factors.size());
        for (std::size_t i = 0; i < factors.size(); ++i)
            CKKSScaleDown<P,SourceSchedule::factor_log_q,SourceSchedule::log_delta,drop>(result[i],factors[i]);
        return result;
    }
};
} // namespace TFHEpp
#endif
