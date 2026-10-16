#pragma once

#include "narrow_selector_projection.hpp"
#include <modular_rns32.hpp>
#include <exception>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
namespace narrow_detail {

template <class P,std::size_t Count=modular_rns32::prime_count,bool VectorShortStages=true>
inline const auto &factorPlansRNS32()
{
    static const auto plans=[] {
        std::vector<modular_ntt32::NegacyclicNTTPlan> result;
        result.reserve(Count);
        for(std::size_t p=0;p<Count;++p) result.emplace_back(P::n,modular_rns32::primes[p],VectorShortStages);
        return result;
    }();
    return plans;
}

template <class Schedule>
struct FactorSelectorRNS32 {
    using P = typename Schedule::Param;
    static constexpr std::size_t digits=1;
    static_assert(Schedule::plain_log_delta<64);
    static constexpr std::size_t prime_count=[] {
        if constexpr(requires{Schedule::selector_rns_prime_count;})return std::size_t(Schedule::selector_rns_prime_count);
        else return modular_rns32::prime_count;
    }();
    static_assert(prime_count>0 && prime_count<=modular_rns32::prime_count);
    static_assert(Schedule::selector_log_q+Schedule::plain_log_delta+P::nbit+
        std::bit_width(std::uint32_t(Schedule::block_width))+1<512);
    static constexpr auto bound=[] {
        const auto selector=MultiLimbUInt<8>{1}<<Schedule::selector_log_q;
        auto value=(selector<<Schedule::plain_log_delta)+selector;
        value=modular_rns32::multiplySmall(value,P::n);
        value=modular_rns32::multiplySmall(value,Schedule::block_width);
        return modular_rns32::multiplySmall(value,2);
    }();
    static_assert(bound<modular_rns32::basisProduct<prime_count>(),"direct selector convolution exceeds CRT range");
    std::array<std::vector<std::uint32_t>,prime_count> spectra;
};

template <class Schedule>
inline void prepareFactorSelectorsRNS32(std::vector<FactorSelectorRNS32<Schedule>> &prepared,
                                    std::span<const typename Schedule::Selector> selectors)
{
    using P = typename Schedule::Param;
    using Cache = FactorSelectorRNS32<Schedule>;
    const auto &plans = factorPlansRNS32<P,Cache::prime_count>();
    static const modular_rns32::CanonicalEncoder<P::T::limbs,Schedule::selector_log_q,Cache::prime_count> encoder;
    prepared.resize(selectors.size());
#pragma omp parallel for schedule(static)
    for (std::size_t j = 0; j < selectors.size(); ++j)
        for (std::size_t prime = 0; prime < Cache::prime_count; ++prime) {
            auto &spectrum = prepared[j].spectra[prime];
            spectrum.resize((P::k+1)*Cache::digits*P::n);
            for (std::size_t component = 0; component <= P::k; ++component)
                for (std::size_t digit = 0; digit < Cache::digits; ++digit) {
                    auto row = std::span(spectrum).subspan((component*Cache::digits+digit)*P::n,P::n);
                    for (std::size_t i = 0; i < P::n; ++i)
                        row[i] = plans[prime].arithmetic().encode(
                            encoder.residue(selectors[j].ct[component][i],prime));
                    plans[prime].forward(row);
                }
        }
}

// Same public tables, selector coefficients, torus modulus and final
// rounding as packedFactorsBatch; only the exact convolution basis changes.
template <class Schedule, bool DirectCRT=true, class FactorKeys, class Loader>
inline std::vector<typename Schedule::Factor> packedFactorsBatchRNS32(
    std::span<const CKKSSlotVector<typename Schedule::Param>> tables,
    const FactorKeys &key, std::size_t chunk_size, Loader load,
    const CKKSPreparedAutoKeyFFT<typename Schedule::Param,Schedule::factor_key_log_q> *factor_fft,
    CKKSNarrowFactorTimings *timings, bool prefetch_raw = false,
    std::span<const FactorSelectorRNS32<Schedule>> retained = {}, bool tiled_matrix = true,
    bool parallel_conjugation = false, bool dot8_matrix = true, bool vector_short_stages = false,
    bool packed_matrix = true,bool all_retained_rows = true,bool joint_components = false,bool parallel_table_ntt = false,bool fused_components = true)
{
    using P = typename Schedule::Param;
    using Cache = FactorSelectorRNS32<Schedule>;
    using Clock = std::chrono::steady_clock;
    if (!chunk_size) throw std::invalid_argument("direct RNS selector chunk must be nonzero");
    if (!retained.empty() && retained.size()!=Schedule::block_width)
        throw std::invalid_argument("retained direct RNS selector group size mismatch");
    if (!retained.empty()) prefetch_raw=false;
    const auto columns = tables.size();
    if (!columns) return {};
    // Three or five to seven real columns share one full SIMD row. Padding is public
    // zero workspace; only actual columns are inverse transformed and returned.
    const auto stride=columns==3?std::size_t{4}:columns>4 && columns<8?std::size_t{8}:columns;
    const bool packed_four=packed_matrix && stride==4;
    const auto &plans = vector_short_stages?factorPlansRNS32<P,Cache::prime_count,true>()
                                           :factorPlansRNS32<P,Cache::prime_count,false>();
    using CRT=std::conditional_t<DirectCRT && Cache::prime_count>=8,
        modular_rns32::DirectCenteredCRT<P::T::limbs,Cache::prime_count>,
        modular_rns32::CenteredCRT<P::T::limbs,Cache::prime_count>>;
    static const CRT crt;
    auto start = Clock::now();
    std::array<std::array<std::vector<std::uint32_t>,Cache::prime_count>,2> base;
    for(std::size_t t=0;t<2;++t)
        for(auto &spectrum:base[t])
            spectrum.resize(packed_four?(t==0?2*stride*P::n:0):stride*P::n);
    std::exception_ptr table_error;
    if(parallel_table_ntt){
        // Convert the wide torus coefficient once, rather than once per RNS
        // prime. Prime transforms then expose enough independent work even
        // when the number of requested columns is smaller than the team.
        using SmallCoefficient=std::pair<bool,std::uint64_t>;
        std::vector<std::vector<SmallCoefficient>> coefficients(2*columns);
#pragma omp parallel for collapse(2) schedule(static)
        for(std::size_t t=0;t<2;++t)for(std::size_t column=0;column<columns;++column){
            try{
                CKKSSlotVector<P> transformed;
                for(std::size_t i=0;i<Schedule::slots;++i){
                    const auto z=tables[column][i];
                    transformed[i]=t==0?(z-std::complex<double>(0,1)*std::conj(z))/2.0
                                       :(std::conj(z)-std::complex<double>(0,1)*z)/2.0;
                }
                auto polynomial=std::make_unique<Polynomial<P>>();
                ckksSlotEncode<P,Schedule::selector_log_q,Schedule::plain_log_delta>(*polynomial,transformed);
                auto&row=coefficients[t*columns+column];row.resize(P::n);
                for(std::size_t i=0;i<P::n;++i){
                    row[i]=ckks_detail::smallSignedMagnitude<P,Schedule::selector_log_q>((*polynomial)[i]);
                    if(row[i].second>(std::uint64_t{1}<<Schedule::plain_log_delta)+1)
                        throw std::invalid_argument("direct RNS table exceeds proved coefficient bound");
                }
            }catch(...){
#pragma omp critical(direct_rns_table_error)
                {if(!table_error)table_error=std::current_exception();}
            }
        }
        if(table_error)std::rethrow_exception(table_error);
#pragma omp parallel for collapse(3) schedule(static)
        for(std::size_t t=0;t<2;++t)for(std::size_t column=0;column<columns;++column)
            for(std::size_t prime=0;prime<Cache::prime_count;++prime){
                try{
                    const auto modulus=modular_rns32::primes[prime].value;
                    const auto&row=coefficients[t*columns+column];
                    std::vector<std::uint32_t> spectrum(P::n);
                    for(std::size_t i=0;i<P::n;++i){
                        const auto [negative,magnitude]=row[i];const auto value=magnitude%modulus;
                        spectrum[i]=negative&&value?modulus-value:value;
                    }
                    plans[prime].forward(spectrum);
                    for(std::size_t i=0;i<P::n;++i)
                        if(packed_four)base[0][prime][i*8+t*4+column]=spectrum[i];
                        else base[t][prime][i*stride+column]=spectrum[i];
                }catch(...){
#pragma omp critical(direct_rns_table_error)
                    {if(!table_error)table_error=std::current_exception();}
                }
            }
        if(table_error)std::rethrow_exception(table_error);
    }else{
#pragma omp parallel for collapse(2) schedule(static)
    for (std::size_t t = 0; t < 2; ++t)
        for (std::size_t column = 0; column < columns; ++column) {
            try {
                CKKSSlotVector<P> transformed;
                for (std::size_t i = 0; i < Schedule::slots; ++i) {
                    const auto z = tables[column][i];
                    transformed[i] = t == 0 ? (z-std::complex<double>(0,1)*std::conj(z))/2.0
                                           : (std::conj(z)-std::complex<double>(0,1)*z)/2.0;
                }
                auto polynomial = std::make_unique<Polynomial<P>>();
                ckksSlotEncode<P,Schedule::selector_log_q,Schedule::plain_log_delta>(*polynomial,transformed);
                for (std::size_t prime = 0; prime < Cache::prime_count; ++prime) {
                    const auto modulus = modular_rns32::primes[prime].value;
                    std::vector<std::uint32_t> spectrum(P::n);
                    for (std::size_t i = 0; i < P::n; ++i) {
                        const auto [negative,magnitude] = ckks_detail::smallSignedMagnitude<P,Schedule::selector_log_q>((*polynomial)[i]);
                        if (magnitude > (std::uint64_t{1}<<Schedule::plain_log_delta)+1)
                            throw std::invalid_argument("direct RNS table exceeds proved coefficient bound");
                        const auto value = magnitude % modulus;
                        spectrum[i] = negative && value ? modulus-value : value;
                    }
                    plans[prime].forward(spectrum);
                    for (std::size_t i = 0; i < P::n; ++i)
                        if(packed_four)base[0][prime][i*8+t*4+column]=spectrum[i];
                        else base[t][prime][i*stride+column]=spectrum[i];
                }
            }
            catch(...) {
#pragma omp critical(direct_rns_table_error)
                { if(!table_error)table_error=std::current_exception(); }
            }
        }
    if(table_error)std::rethrow_exception(table_error);
    }
    if (timings) timings->table_ntt_seconds += std::chrono::duration<double>(Clock::now()-start).count();
    std::vector<std::size_t> exponents(Schedule::block_width);
    for (std::size_t y = 0; y < exponents.size(); ++y)
        exponents[y] = modular_ntt::power(5,(Schedule::slots-y)%Schedule::slots,2*P::n);
    using Sums = std::array<std::vector<std::uint32_t>,Cache::prime_count>;
    std::vector<Sums> sums(2*(P::k+1)*Cache::digits);
    auto index = [](std::size_t t,std::size_t component,std::size_t digit) {
        return (t*(P::k+1)+component)*Cache::digits+digit;
    };
    for (auto &row : sums) for (auto &spectrum : row) spectrum.resize(stride*P::n);
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
        // Retained spectra are all available: finish each point tile while
        // its public tables are hot, instead of revisiting every prime tile
        // after each I/O-sized chunk. Streaming keeps its bounded workspace.
        const auto online_chunk=all_retained_rows && !retained.empty()?Schedule::block_width:chunk_size;
        const auto count = std::min(online_chunk,Schedule::block_width-first_y);
        std::vector<Cache> transient;
        std::span<const Cache> cache;
        if (!retained.empty()) cache=retained.subspan(first_y,count);
        else {
            std::vector<typename Schedule::Selector> raw;
            start = Clock::now();
            if (prefetch_raw) raw = pending.get();
            else load(first_y,count,raw);
            if (raw.size() != count) throw std::invalid_argument("direct RNS selector chunk count mismatch");
            if (timings) timings->selector_load_seconds += std::chrono::duration<double>(Clock::now()-start).count();
            start = Clock::now();
            prepareFactorSelectorsRNS32<Schedule>(transient,raw);
            std::vector<typename Schedule::Selector>{}.swap(raw);
            if (prefetch_raw && first_y+count < Schedule::block_width) begin_read(first_y+count);
            if (timings) timings->selector_ntt_seconds += std::chrono::duration<double>(Clock::now()-start).count();
            cache=transient;
        }
        start = Clock::now();
        const auto tile=[&]<bool Joint>(std::size_t prime,std::size_t component,std::size_t digit,std::size_t begin){
            const auto &arithmetic=plans[prime].arithmetic();
            auto *first=sums[index(0,component,digit)][prime].data();
            auto *second=sums[index(1,component,digit)][prime].data();
                        const auto end = std::min(begin+256,std::size_t(P::n));
                        std::size_t local=0;
                        const auto dot=[&]<std::size_t Terms,bool Packed=false>(){
                            std::array<const std::uint32_t*,Terms> selectors;
                            for(std::size_t t=0;t<Terms;++t)selectors[t]=cache[local+t].spectra[prime].data()+(component*Cache::digits+digit)*P::n;
                            for(std::size_t i=begin;i<end;++i){
                                std::array<const std::uint32_t*,Terms> left,right;
                                std::array<std::uint32_t,Terms> values;
                                for(std::size_t t=0;t<Terms;++t){
                                    const auto position=narrowNTTPermutationIndex<P::n>(exponents[first_y+local+t],i)*(Packed?8:stride);
                                    left[t]=base[0][prime].data()+position;
                                    if constexpr(!Packed)right[t]=base[1][prime].data()+position;
                                    values[t]=selectors[t][i];
                                }
                                if constexpr(Joint && P::k==1)if(fused_components){
                                    std::array<std::uint32_t,Terms> values1;
                                    for(std::size_t t=0;t<Terms;++t)values1[t]=cache[local+t].spectra[prime][Cache::digits*P::n+digit*P::n+i];
                                    arithmetic.template addDotPairComponents<Terms,Packed>(first+i*stride,second+i*stride,
                                        sums[index(0,1,digit)][prime].data()+i*stride,sums[index(1,1,digit)][prime].data()+i*stride,
                                        left,right,values,values1,stride);
                                    continue;
                                }
                                if constexpr(Packed)arithmetic.template addPackedDotPair<Terms>(first+i*4,second+i*4,left,values);
                                else arithmetic.template addDotPair<Terms>(first+i*stride,second+i*stride,left,right,values,stride);
                                // The public row addresses are shared by all ciphertext
                                // components; reuse them while their data is hot.
                                if constexpr(Joint)for(std::size_t c=component+1;c<=P::k;++c){
                                    for(std::size_t t=0;t<Terms;++t)
                                        values[t]=cache[local+t].spectra[prime][(c*Cache::digits+digit)*P::n+i];
                                    auto*a=sums[index(0,c,digit)][prime].data()+i*stride;
                                    auto*b=sums[index(1,c,digit)][prime].data()+i*stride;
                                    if constexpr(Packed)arithmetic.template addPackedDotPair<Terms>(a,b,left,values);
                                    else arithmetic.template addDotPair<Terms>(a,b,left,right,values,stride);
                                }
                            }
                        };
                        // Four columns pack both signed tables in one AVX2 vector.
                        // Wider rows retain four-term dots: eight terms were slower
                        // in the full-size matrix benchmark.
                        if(packed_four){
                            if(dot8_matrix)for(;local+8<=count;local+=8)dot.template operator()<8,true>();
                            for(;local+4<=count;local+=4)dot.template operator()<4,true>();
                        }else{
                            if(dot8_matrix && stride==4)
                                for(;local+8<=count;local+=8)dot.template operator()<8>();
                            for(;local+4<=count;local+=4)dot.template operator()<4>();
                        }
                        for (; local < count; ++local) {
                            const auto *selector = cache[local].spectra[prime].data()+(component*Cache::digits+digit)*P::n;
                            for (std::size_t i = begin; i < end; ++i) {
                                const auto position = narrowNTTPermutationIndex<P::n>(exponents[first_y+local],i)*(packed_four?8:stride);
                                for(std::size_t c=component;c<(Joint?P::k+1:component+1);++c)
                                    arithmetic.addPair(sums[index(0,c,digit)][prime].data()+i*stride,
                                        sums[index(1,c,digit)][prime].data()+i*stride,
                                        base[0][prime].data()+position,packed_four?base[0][prime].data()+position+4:
                                        base[1][prime].data()+position,stride,
                                        c==component?selector[i]:cache[local].spectra[prime][(c*Cache::digits+digit)*P::n+i]);
                            }
                        }
                        for (std::size_t i = begin*stride; i < end*stride; ++i) {
                            for(std::size_t c=component;c<(Joint?P::k+1:component+1);++c){
                                auto&a=sums[index(0,c,digit)][prime][i];auto&b=sums[index(1,c,digit)][prime][i];
                                a=arithmetic.canonical(a);b=arithmetic.canonical(b);
                            }
                        }
        };
        if(tiled_matrix){
#pragma omp parallel for collapse(4) schedule(dynamic,1)
            for(std::size_t prime=0;prime<Cache::prime_count;++prime)
                for(std::size_t component=0;component<(joint_components?1:P::k+1);++component)
                    for(std::size_t digit=0;digit<Cache::digits;++digit)
                        for(std::size_t begin=0;begin<P::n;begin+=256){
                            if(joint_components)tile.template operator()<true>(prime,component,digit,begin);
                            else tile.template operator()<false>(prime,component,digit,begin);
                        }
        }else{
#pragma omp parallel for collapse(3) schedule(static)
            for(std::size_t prime=0;prime<Cache::prime_count;++prime)
                for(std::size_t component=0;component<(joint_components?1:P::k+1);++component)
                    for(std::size_t digit=0;digit<Cache::digits;++digit)
                        for(std::size_t begin=0;begin<P::n;begin+=256){
                            if(joint_components)tile.template operator()<true>(prime,component,digit,begin);
                            else tile.template operator()<false>(prime,component,digit,begin);
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
                std::array<std::vector<std::uint32_t>,Cache::prime_count> spectra;
                for (std::size_t prime = 0; prime < Cache::prime_count; ++prime) {
                    auto &spectrum=spectra[prime];spectrum.resize(P::n);
                    const auto &sum=sums[index(t,component,0)][prime];
                    for(std::size_t i=0;i<P::n;++i)spectrum[i]=sum[i*stride+column];
                    plans[prime].inverse(spectrum);
                }
                std::array<std::span<std::uint32_t>,Cache::prime_count> residues;
                for(std::size_t prime=0;prime<Cache::prime_count;++prime)residues[prime]=spectra[prime];
                std::vector<typename P::T> coefficients(P::n);
                crt.reconstruct(residues,coefficients);
                for(std::size_t i=0;i<P::n;++i)
                    accum[column][t][component][i]=ckks_detail::reduceToLevel<P,Schedule::selector_log_q>(coefficients[i]);
            }
    sums.clear();
    if (timings) timings->inverse_crt_seconds += std::chrono::duration<double>(Clock::now()-start).count();
    start = Clock::now();
    std::vector<typename Schedule::Factor> result(columns);
    const auto finish_column = [&](std::size_t column) {
        auto second = std::make_unique<typename Schedule::Factor>();
        auto conjugated = std::make_unique<typename Schedule::Factor>();
        for (std::size_t component = 0; component <= P::k; ++component)
            for (std::size_t i = 0; i < P::n; ++i) {
                result[column].ct[component][i] = ckks_detail::roundedLevelRightShift<P,Schedule::selector_log_q,Schedule::plain_log_delta>(accum[column][0][component][i]);
                second->ct[component][i] = ckks_detail::roundedLevelRightShift<P,Schedule::selector_log_q,Schedule::plain_log_delta>(accum[column][1][component][i]);
            }
        if constexpr (Schedule::factor_log_q <= Schedule::factor_key_log_q)
            conjugate<Schedule,true>(conjugated->ct,second->ct,key,factor_fft);
        else {
            if (!key.automorphism_evaluator)
                throw std::invalid_argument("high-scale factors require a public automorphism backend");
            key.automorphism_evaluator(conjugated->ct,second->ct,Schedule::factor_log_q,2*P::n-1);
        }
        CKKSAddInPlace<P>(result[column],*conjugated);
    };
    // Public prepared caches may be shared by concurrent column evaluators.
    // Opt in only with a thread-safe callback; inner OpenMP work follows the
    // caller's nesting policy (one active level avoids oversubscription).
    if (parallel_conjugation && columns > 1) {
        std::exception_ptr error;
        [[maybe_unused]] const auto threads=columnThreads(columns);
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for (std::size_t column = 0; column < columns; ++column) {
            try { finish_column(column); }
            catch (...) {
#pragma omp critical(ckks_factor_conjugation_error)
                { if (!error) error = std::current_exception(); }
            }
        }
        if (error) std::rethrow_exception(error);
    }
    else for (std::size_t column = 0; column < columns; ++column) finish_column(column);
    if (timings) timings->conjugation_seconds += std::chrono::duration<double>(Clock::now()-start).count();
    return result;
}
} // namespace narrow_detail

// Key-dependent spectra only. The caller chooses a byte budget before any
// group is loaded; input tables and evaluated factors are never cached here.
template <class Schedule,class SelectorSchedule=Schedule>
class CKKSNarrowFactorSelectorCacheRNS32 {
    using P=typename Schedule::Param;
    using Source=CKKSNarrowSelectorSource<SelectorSchedule>;
    using Cache=narrow_detail::FactorSelectorRNS32<Schedule>;
    static_assert(std::is_same_v<typename Schedule::Selector,typename SelectorSchedule::Selector>);
    std::shared_ptr<const Source> source_;
    std::vector<std::vector<Cache>> groups_;
public:
    static constexpr std::uint64_t payload_bytes=std::uint64_t{Cache::prime_count}*sizeof(std::uint32_t)*
        (P::k+1)*Cache::digits*P::n*Schedule::block_width*(Schedule::blocks+Schedule::tail_packs);
    // Includes vector object overhead; allocation bookkeeping and temporary
    // preparation buffers remain additional bounded process workspace.
    static constexpr std::uint64_t retained_bytes=payload_bytes+
        sizeof(Cache)*std::uint64_t(Schedule::block_width)*(Schedule::blocks+Schedule::tail_packs)+
        sizeof(std::vector<Cache>)*(Schedule::blocks+Schedule::tail_packs);
    CKKSNarrowFactorSelectorCacheRNS32(std::shared_ptr<const Source> source,
        std::uint64_t budget_bytes,std::size_t chunk_size,
        std::function<void(bool,std::size_t)> observer={}) : source_(std::move(source))
    {
        if (!source_ || !source_->hasRawChunks() || !chunk_size || budget_bytes<retained_bytes)
            throw std::invalid_argument("invalid direct RNS cache source, chunk or memory budget");
        if (source_->size(false)!=Schedule::blocks*Schedule::block_width ||
            source_->size(true)!=Schedule::tail_packs*Schedule::block_width)
            throw std::invalid_argument("direct RNS cache source geometry mismatch");
        groups_.resize(Schedule::blocks+Schedule::tail_packs);
        for (bool tail:{false,true})
            for (std::size_t group=0;group<(tail?Schedule::tail_packs:Schedule::blocks);++group) {
                auto &destination=groups_[(tail?Schedule::blocks:0)+group];
                destination.resize(Schedule::block_width);
                for (std::size_t first=0;first<Schedule::block_width;) {
                    const auto count=std::min(chunk_size,Schedule::block_width-first);
                    std::vector<typename Schedule::Selector> raw;
                    source_->loadChunk(tail,group,first,count,raw);
                    if (raw.size()!=count) throw std::invalid_argument("direct RNS cache raw count mismatch");
                    std::vector<Cache> prepared;
                    narrow_detail::prepareFactorSelectorsRNS32<Schedule>(prepared,raw);
                    std::move(prepared.begin(),prepared.end(),destination.begin()+first);
                    first+=count;
                }
                if (observer) observer(tail,group);
            }
    }
    const Source *source() const { return source_.get(); }
    std::span<const Cache> group(bool tail,std::size_t index) const
    {
        if (index>=(tail?Schedule::tail_packs:Schedule::blocks))
            throw std::out_of_range("retained direct RNS group");
        return groups_[(tail?Schedule::blocks:0)+index];
    }
};

template <class TargetSchedule,class SourceSchedule,class SelectorSchedule=SourceSchedule,bool DirectCRT=true>
class CKKSNarrowDeferredRescaleFactorSourceRNS32 final : public CKKSNarrowFactorSource<TargetSchedule> {
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
    bool tiled_matrix_,parallel_conjugation_,dot8_matrix_,vector_short_stages_,packed_matrix_,all_retained_rows_,joint_components_,parallel_table_ntt_,fused_components_;
    struct FactorKeys {
        const CKKSGaloisKey<P,SourceSchedule::factor_key_log_q> &factor_galois;
        bool staged_fft_caches;
        const decltype(CKKSNarrowBootstrapKey<TargetSchedule>::automorphism_evaluator) &automorphism_evaluator;
        const decltype(CKKSNarrowBootstrapKey<SourceSchedule>::factor_conjugation_fft) &factor_conjugation_fft;
        const decltype(CKKSNarrowBootstrapKey<SourceSchedule>::factor_conjugation_projected_fft) &factor_conjugation_projected_fft;
    };
public:
    explicit CKKSNarrowDeferredRescaleFactorSourceRNS32(std::shared_ptr<const Source> source,
        std::shared_ptr<const CKKSNarrowFactorSelectorCacheRNS32<SourceSchedule,SelectorSchedule>> retained={},bool tiled_matrix=true,
        bool parallel_conjugation=false,bool dot8_matrix=true,bool vector_short_stages=false,bool packed_matrix=true,bool all_retained_rows=true,bool joint_components=false,bool parallel_table_ntt=false,bool fused_components=true)
        : source_(std::move(source)),retained_(std::move(retained)),tiled_matrix_(tiled_matrix),parallel_conjugation_(parallel_conjugation),dot8_matrix_(dot8_matrix),vector_short_stages_(vector_short_stages),packed_matrix_(packed_matrix),all_retained_rows_(all_retained_rows),joint_components_(joint_components),parallel_table_ntt_(parallel_table_ntt),fused_components_(fused_components)
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
            key.factor_conjugation_fft,key.factor_conjugation_projected_fft};
        auto factors = narrow_detail::packedFactorsBatchRNS32<SourceSchedule,DirectCRT>(tables,view,key.selector_ntt_chunk_size,
            [&](auto first,auto count,auto &raw) { source_->loadChunk(tail,group,first,count,raw); },cache,timings,
            source_->supportsRawChunkPrefetch(),retained_?retained_->group(tail,group):
                std::span<const narrow_detail::FactorSelectorRNS32<SourceSchedule>>{},tiled_matrix_,parallel_conjugation_,dot8_matrix_,vector_short_stages_,packed_matrix_,all_retained_rows_,joint_components_,parallel_table_ntt_,fused_components_);
        std::vector<typename TargetSchedule::Factor> result(factors.size());
        for (std::size_t i = 0; i < factors.size(); ++i)
            CKKSScaleDown<P,SourceSchedule::factor_log_q,SourceSchedule::log_delta,drop>(result[i],factors[i]);
        return result;
    }
};
} // namespace TFHEpp
#endif
