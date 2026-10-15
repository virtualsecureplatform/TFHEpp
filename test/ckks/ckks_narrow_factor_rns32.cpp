#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_factor_ntt32.hpp>
#include <ckks/narrow_factor_rns32.hpp>
#include <ckks/narrow_dd_ntt.hpp>
#include <iostream>
#include <atomic>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<9,7>;
using Source=CKKSNarrowBootstrapSchedule<P,48,440,57,42,31,false,0,true>;
using Target=CKKSNarrowBootstrapSchedule<P,48,427,44,42,31,true,13,true>;

struct ReducedSource:CKKSNarrowBootstrapSchedule<P,48,397,57,42,31,false,0,true>{static constexpr std::size_t selector_rns_prime_count=16;};
using LeanTarget=CKKSNarrowBootstrapSchedule<P,48,377,37,42,31,true,20,true>;

class TransparentSource final:public CKKSNarrowSelectorSource<Source> {
public:
    mutable std::size_t reads=0;
    std::size_t size(bool tail) const override
    { return (tail?Source::tail_packs:Source::blocks)*Source::block_width; }
    bool hasRawChunks() const override { return true; }
    void load(bool tail,std::size_t group,std::vector<Source::Selector>& result) const override
    { loadChunk(tail,group,0,Source::block_width,result); }
    void loadChunk(bool tail,std::size_t group,std::size_t first,std::size_t count,
                   std::vector<Source::Selector>& result) const override
    {
        ++reads;result.resize(count);
        for(std::size_t j=0;j<count;++j) {
            result[j]={};CKKSSlotVector<P> slots{};
            std::mt19937_64 random(20261015+group*1009+(tail?100000:0)+first+j);
            for(auto &coefficient:result[j].ct[0]) {
                for(auto &word:coefficient.limb)word=random();
                coefficient=ckks_detail::reduceToLevel<P,440>(coefficient);
            }
            if(!tail && first+j==group%Source::block_width)
                std::fill_n(slots.begin(),Source::head_slots,std::complex<double>{1,0});
            if(tail)for(std::size_t ell=0;ell<Source::tail_blocks_per_pack;++ell) {
                const auto block=group*Source::tail_blocks_per_pack+ell;
                if(block>=Source::blocks)break;
                if(first+j==block%Source::block_width)
                    std::fill_n(slots.begin()+(2*ell+1)*Source::tail_width,Source::tail_width,std::complex<double>{0,1});
            }
            ckksSlotEncode<P,440,57>(result[j].ct[1],slots);
        }
    }
};
template<class Chain>
struct FailingBatchKeys {
    template<std::size_t Level>
    auto get()const -> decltype(std::declval<const Chain&>().template get<Level>())
    {throw std::runtime_error("batch tree backend test failure");}
};

template<uint32_t Nbit,std::size_t H>
void checkBatchTrees()
{
    using BP=narrow_benchmark::BenchmarkParam<Nbit,7>;
    using S=CKKSNarrowBootstrapSchedule<BP,48,427,44,42,H,true,13,true>;
    using Chain=CKKSPreparedAuxiliaryDDRelinKeyChainNTT<BP,S::factor_log_q,44,S::tree_depth,48,20,20,true,true>;
    auto master=std::make_unique<CKKSCompactDDRelinKey<BP,433,20,20>>();
    std::mt19937_64 random(20261017+Nbit);
    for(auto&row:master->data)for(auto&digit:row)for(auto&component:digit)for(auto&value:component)value=int32_t(random()&((1U<<20)-1))-(1<<19);
    Chain prepared(*master);auto raw=std::make_unique<typename S::RelinKeys>();
    const auto maximum=std::size_t{1}<<S::tree_depth;
    for(std::size_t count:{0,1,7,14}){
        std::vector<std::vector<typename S::Factor>> factors(count);
        for(std::size_t tree=0;tree<count;++tree){
            const std::size_t lengths[]={maximum,maximum-1,maximum/2+1,1};
            factors[tree].resize(lengths[tree%4]);
            for(auto&factor:factors[tree])for(auto&component:factor.ct)for(auto&value:component){
                for(auto&word:value.limb)word=random();
                value=ckks_detail::reduceToLevel<BP,S::factor_log_q>(value);
            }
        }
        std::vector<typename S::TreeResult> expected(count),actual;
        for(std::size_t tree=0;tree<count;++tree)narrow_detail::productTree<S,0,Chain>(expected[tree],factors[tree],*raw,&prepared);
        narrow_detail::productTreeBatch<S,0,Chain>(actual,std::move(factors),*raw,&prepared);
        if(actual.size()!=expected.size())throw std::runtime_error("batched tree changed result count");
        for(std::size_t tree=0;tree<count;++tree)if(actual[tree].ct!=expected[tree].ct)throw std::runtime_error("batched tree changed paired arithmetic or odd-leaf carry");
    }
    if constexpr(use_multilimb_digit_fft_v<BP>){
        FailingBatchKeys<Chain> failing;bool rejected=false;
        std::vector<typename S::TreeResult> result;
        try{narrow_detail::productTreeBatch<S,0,FailingBatchKeys<Chain>>(result,
            std::vector<std::vector<typename S::Factor>>(2,std::vector<typename S::Factor>(2)),*raw,&failing);}
        catch(const std::runtime_error&error){rejected=std::string(error.what())=="batch tree backend test failure";}
        if(!rejected)throw std::runtime_error("batched tree did not propagate backend failure");
    }
    for(std::size_t length:{std::size_t{0},maximum+1}){
        bool rejected=false;std::vector<typename S::TreeResult> result;
        try{narrow_detail::productTreeBatch<S,0,Chain>(result,{std::vector<typename S::Factor>(length)},*raw,&prepared);}
        catch(const std::invalid_argument&){rejected=true;}
        if(!rejected)throw std::runtime_error("invalid batched tree shape accepted");
    }
}

void checkRootColumnDispatch()
{
    for(std::size_t columns:{0,1,4,7}){
        std::vector<std::size_t> expected(columns),actual(columns);
        narrow_detail::evaluateRootColumns(columns,false,[&](std::size_t column){expected[column]=column*column+17;});
        std::atomic<std::size_t> parallel_calls{0};
        std::atomic<int> team_size{0};
#ifdef _OPENMP
        const auto caller_threads=omp_get_max_threads();
#endif
        narrow_detail::evaluateRootColumns(columns,true,[&](std::size_t column){
#ifdef _OPENMP
            if(omp_in_parallel()){
                parallel_calls.fetch_add(1,std::memory_order_relaxed);
                team_size.store(omp_get_num_threads(),std::memory_order_relaxed);
            }
#endif
            actual[column]=column*column+17;
        });
        if(expected!=actual)throw std::runtime_error("root column dispatch changed indexed outputs");
#ifdef _OPENMP
        if(columns>1 && omp_get_max_threads()>1 && !parallel_calls.load())throw std::runtime_error("root column callbacks did not enter an active team");
        if(std::size_t(team_size.load())>columns)throw std::runtime_error("root column team exceeds available work");
        if(omp_get_max_threads()!=caller_threads)throw std::runtime_error("root column dispatch changed caller thread budget");
#endif
    }
    bool rejected=false;
    try{narrow_detail::evaluateRootColumns(7,true,[](std::size_t column){if(column==2)throw std::runtime_error("root column test failure");});}
    catch(const std::runtime_error& error){rejected=std::string(error.what())=="root column test failure";}
    if(!rejected)throw std::runtime_error("root column dispatch did not propagate an exception");
}

void checkTailAlignment()
{
    auto key=std::make_unique<CKKSNarrowBootstrapKey<Target>>();
    key->batch_product_trees=true;
    bool missing_rejected=false;std::vector<Target::Output> unused;
    try{CKKSNarrowBlindRotateBatch<Target,44>(unused,{},*key);}
    catch(const std::invalid_argument&error){missing_rejected=std::string(error.what())=="batched product trees require a prepared batch backend";}
    if(!missing_rejected)throw std::runtime_error("missing prepared batch backend accepted");
    key->batch_product_trees=false;
    std::atomic<std::size_t> parallel_calls{0};
    key->automorphism_evaluator=[&](auto& out,const auto& in,uint32_t level,uint32_t exponent){
        if(level!=Target::factor_log_q)throw std::runtime_error("incorrect public tail level");
#ifdef _OPENMP
        if(omp_in_parallel())parallel_calls.fetch_add(1,std::memory_order_relaxed);
#endif
        auto temporary=std::make_unique<TRLWE<P>>();
        for(std::size_t c=0;c<2;++c)Automorphism<P>((*temporary)[c],in[c],exponent);
        ckks_detail::reduceTRLWEToLevel<P,Target::factor_log_q>(*temporary);out=*temporary;
    };
    for(std::size_t columns:{0,1,4,7}){
        std::vector<std::vector<Target::Factor>> source(columns,std::vector<Target::Factor>(Target::blocks));
        std::vector<int> steps(columns*Target::blocks);
        for(std::size_t column=0;column<columns;++column)for(std::size_t block=0;block<Target::blocks;++block){
            steps[column*Target::blocks+block]=-int(Target::slots-2*(block%Target::tail_blocks_per_pack+1)*Target::tail_width);
            for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)
                source[column][block].ct[c][i]=P::T{1}<<((column+block+c+i)%Target::factor_log_q);
        }
        auto serial=source,parallel=source;
        key->parallel_public_tail_alignment=false;narrow_detail::alignPublicTail<Target>(serial,steps,*key);
        key->parallel_public_tail_alignment=true;narrow_detail::alignPublicTail<Target>(parallel,steps,*key);
        for(std::size_t column=0;column<columns;++column)for(std::size_t block=0;block<Target::blocks;++block)
            if(serial[column][block].ct!=parallel[column][block].ct)throw std::runtime_error("parallel tail alignment changed ciphertexts");
    }
#ifdef _OPENMP
    if(omp_get_max_threads()>1&&!parallel_calls.load())throw std::runtime_error("tail alignment callbacks did not enter an active team");
#endif
    // Archive only zero-valued public rows; no selector-key fixture generation.
    std::stringstream serial_payload,payload;
    key->parallel_public_tail_alignment=false;key->parallel_root_columns=false;key->batch_product_trees=false;
    {cereal::PortableBinaryOutputArchive archive(serial_payload);archive(*key);}
    key->parallel_public_tail_alignment=true;key->parallel_root_columns=true;key->batch_product_trees=true;
    key->product_tree_batch_evaluator=[](auto&,auto,const auto&){};
    {cereal::PortableBinaryOutputArchive archive(payload);archive(*key);}
    if(serial_payload.str()!=payload.str())throw std::runtime_error("runtime tail policy changed public archive bytes");
    auto loaded=std::make_unique<CKKSNarrowBootstrapKey<Target>>();
    loaded->parallel_public_tail_alignment=true;loaded->parallel_root_columns=true;loaded->batch_product_trees=true;loaded->product_tree_batch_evaluator=key->product_tree_batch_evaluator;
    {cereal::PortableBinaryInputArchive archive(payload);archive(*loaded);}
    if(loaded->parallel_public_tail_alignment || loaded->parallel_root_columns || loaded->batch_product_trees || loaded->product_tree_batch_evaluator)throw std::runtime_error("key archive retained a runtime parallel alignment policy");
    loaded->batch_product_trees=true;loaded->product_tree_batch_evaluator=key->product_tree_batch_evaluator;
    CKKSNarrowUseStagedFFT(*loaded,4);
    if(loaded->batch_product_trees||loaded->product_tree_batch_evaluator)throw std::runtime_error("staged policy retained an incompatible batch backend");
    std::vector<std::vector<Target::Factor>> tail(1,std::vector<Target::Factor>(Target::blocks));
    std::vector<int> steps(Target::blocks,1);
    key->automorphism_evaluator=[](auto&,const auto&,uint32_t,uint32_t){throw std::runtime_error("tail alignment test failure");};
    bool rejected=false;
    try{narrow_detail::alignPublicTail<Target>(tail,steps,*key);}
    catch(const std::runtime_error& error){rejected=std::string(error.what())=="tail alignment test failure";}
    if(!rejected)throw std::runtime_error("parallel tail alignment did not propagate an exception");
    steps.pop_back();rejected=false;
    try{narrow_detail::alignPublicTail<Target>(tail,steps,*key);}
    catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("malformed public tail batch accepted");
}

int main()
{
    checkTailAlignment();
    checkRootColumnDispatch();
    checkBatchTrees<4,4>();
    checkBatchTrees<9,31>();
    auto raw=std::make_shared<TransparentSource>();
    using Cache=CKKSNarrowFactorSelectorCache32<Source>;
    bool rejected=false;
    try { Cache insufficient(raw,Cache::retained_bytes-1,4); }
    catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected || raw->reads)throw std::runtime_error("cache budget did not reject before loading");
    auto cache=std::make_shared<Cache>(raw,Cache::retained_bytes,4);
    CKKSNarrowDeferredRescaleFactorSource32<Target,Source> streamed(raw),cached(raw,cache);
    using DirectCache=CKKSNarrowFactorSelectorCacheRNS32<Source>;
    rejected=false;const auto prior_reads=raw->reads;
    try { DirectCache bad(raw,DirectCache::retained_bytes-1,4); }
    catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected || raw->reads!=prior_reads)throw std::runtime_error("direct RNS budget guard failed");
    auto direct_cache=std::make_shared<DirectCache>(raw,DirectCache::retained_bytes,4);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> direct(raw,direct_cache),direct_streamed(raw);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> simd_short(raw,direct_cache,true,false,true,true);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> four_term(raw,direct_cache,true,false,false);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> unpacked(raw,direct_cache,true,false,true,false,false);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> joint(raw,direct_cache,true,false,true,false,true,true,true);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> joint_streamed(raw,{},true,false,true,false,true,true,true);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> chunked(raw,direct_cache,true,false,true,false,true,false);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source,Source,false> garner(raw,direct_cache);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> untiled(raw,direct_cache,false);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<Target,Source> parallel(raw,direct_cache,true,true);
    auto key=std::make_unique<CKKSNarrowBootstrapKey<Target>>();
    key->selector_ntt_chunk_size=4;
    std::atomic<std::size_t> parallel_callbacks{0};
    std::atomic<int> conjugation_team_size{0};
    key->automorphism_evaluator=[&](auto &result,const auto &input,uint32_t,uint32_t exponent) {
#ifdef _OPENMP
        if(omp_in_parallel()){
            parallel_callbacks.fetch_add(1,std::memory_order_relaxed);
            conjugation_team_size.store(omp_get_num_threads(),std::memory_order_relaxed);
        }
#endif
        for(std::size_t c=0;c<2;++c)Automorphism<P>(result[c],input[c],exponent);
    };
    std::vector<CKKSSlotVector<P>> tables(4);
    for(std::size_t j=0;j<tables.size();++j)
        for(std::size_t i=0;i<P::n/2;++i)tables[j][i]=std::polar(1.0,0.03*(j+1)*i);
    for(bool tail:{false,true}) {
        CKKSNarrowFactorTimings timing;
        const auto reference=streamed.evaluate(tables,tail,0,*key,nullptr,nullptr);
        const auto reads=raw->reads;
        const auto result=cached.evaluate(tables,tail,0,*key,nullptr,&timing);
        const auto exact=direct.evaluate(tables,tail,0,*key,nullptr,&timing);
        for(std::size_t j=0;j<tables.size();++j)
            if(exact[j].ct!=reference[j].ct)throw std::runtime_error("direct RNS factor differs from digit CRT");
        if(raw->reads!=reads || timing.selector_load_seconds || timing.selector_ntt_seconds)
            throw std::runtime_error("cached evaluation loaded or transformed raw selectors");
        const auto direct_uncached=direct_streamed.evaluate(tables,tail,0,*key,nullptr,nullptr);
        for(std::size_t j=0;j<tables.size();++j) {
            if(reference[j].ct!=result[j].ct)throw std::runtime_error("cached factor differs from streamed factor");
            if(direct_uncached[j].ct!=reference[j].ct)throw std::runtime_error("streamed direct RNS differs from digit CRT");
        }
    }
    // Exercise scalar tails, the four-column paired SIMD layout and the
    // nine-column vector/scalar boundary with non-multiple-of-four chunks.
    for(std::size_t columns:{1,2,3,4,5,6,7,9}){
        std::vector<CKKSSlotVector<P>> variants(columns);
        for(std::size_t c=0;c<columns;++c)for(std::size_t i=0;i<P::n/2;++i)variants[c][i]=std::polar(1.0,0.017*(c+1)*i);
        for(bool tail:{false,true}){
            key->selector_ntt_chunk_size=1;
            const auto reference=streamed.evaluate(variants,tail,0,*key,nullptr,nullptr);
            for(std::size_t chunk:{3,4,5,8,9,12,17,40}){
                key->selector_ntt_chunk_size=chunk;
                const auto actual=direct.evaluate(variants,tail,0,*key,nullptr,nullptr);
                const auto shared=joint.evaluate(variants,tail,0,*key,nullptr,nullptr);
                const auto shared_streamed=joint_streamed.evaluate(variants,tail,0,*key,nullptr,nullptr);
                for(std::size_t c=0;c<columns;++c)if(shared[c].ct!=actual[c].ct || shared_streamed[c].ct!=actual[c].ct)throw std::runtime_error("joint matrix components differ from separate components");
                const auto control=four_term.evaluate(variants,tail,0,*key,nullptr,nullptr);
                const auto separate=unpacked.evaluate(variants,tail,0,*key,nullptr,nullptr);
                const auto split=chunked.evaluate(variants,tail,0,*key,nullptr,nullptr);
                for(std::size_t c=0;c<columns;++c)if(split[c].ct!=actual[c].ct)throw std::runtime_error("retained whole-row tiles differ from chunked tiles");
                for(std::size_t c=0;c<columns;++c)if(separate[c].ct!=actual[c].ct)throw std::runtime_error("packed tables differ from separate table layout");
                const auto alternate=simd_short.evaluate(variants,tail,0,*key,nullptr,nullptr);
                for(std::size_t c=0;c<columns;++c)if(alternate[c].ct!=actual[c].ct)throw std::runtime_error("SIMD factor short stages differ from scalar short stages");
                const auto previous=garner.evaluate(variants,tail,0,*key,nullptr,nullptr);
                for(std::size_t c=0;c<columns;++c)if(previous[c].ct!=actual[c].ct)throw std::runtime_error("direct CRT factors differ from Garner CRT");
                for(std::size_t c=0;c<columns;++c)if(control[c].ct!=actual[c].ct)throw std::runtime_error("eight-term dots differ from four-term dots");
                conjugation_team_size.store(0);
                const auto concurrent=parallel.evaluate(variants,tail,0,*key,nullptr,nullptr);
                if(std::size_t(conjugation_team_size.load())>columns)throw std::runtime_error("factor conjugation team exceeds available columns");
                const auto coarse=untiled.evaluate(variants,tail,0,*key,nullptr,nullptr);
                for(std::size_t c=0;c<columns;++c)if(concurrent[c].ct!=reference[c].ct)throw std::runtime_error("parallel factor conjugations differ from digit CRT");
                for(std::size_t c=0;c<columns;++c)if(coarse[c].ct!=reference[c].ct)throw std::runtime_error("coarse RNS matrix differs from digit CRT");
                for(std::size_t c=0;c<columns;++c)if(actual[c].ct!=reference[c].ct)throw std::runtime_error("fused dot or remainder differs from single-term factors");
            }
        }
    }
#ifdef _OPENMP
    if(omp_get_max_threads()>1&&!parallel_callbacks.load())throw std::runtime_error("parallel factor callbacks did not enter an active OpenMP team");
#endif
    const auto original_auto=key->automorphism_evaluator;
    key->automorphism_evaluator=[](auto&,const auto&,uint32_t,uint32_t){throw std::runtime_error("conjugation test failure");};
    rejected=false;
    try{parallel.evaluate(tables,false,0,*key,nullptr,nullptr);}
    catch(const std::runtime_error& error){rejected=std::string(error.what())=="conjugation test failure";}
    if(!rejected)throw std::runtime_error("parallel conjugation did not propagate backend exception");
    key->automorphism_evaluator=original_auto;
    auto reduced=std::make_shared<CKKSNarrowReducedSelectorSource<ReducedSource,Source>>(raw);
    std::vector<Source::Selector> original_chunk;std::vector<ReducedSource::Selector> reduced_chunk;
    raw->loadChunk(false,0,0,3,original_chunk);reduced->loadChunk(false,0,0,3,reduced_chunk);
    for(std::size_t j=0;j<3;++j)for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)
        if(reduced_chunk[j].ct[c][i]!=ckks_detail::reduceToLevel<P,397>(original_chunk[j].ct[c][i]))throw std::runtime_error("selector reduction changed scale or coefficients");
    using ReducedCache=CKKSNarrowFactorSelectorCacheRNS32<ReducedSource>;
    auto reduced_cache=std::make_shared<ReducedCache>(reduced,ReducedCache::retained_bytes,4);
    CKKSNarrowDeferredRescaleFactorSourceRNS32<LeanTarget,ReducedSource> lean(reduced,reduced_cache);
    CKKSNarrowDeferredRescaleFactorSource32<LeanTarget,ReducedSource> lean_reference(reduced);
    auto lean_key=std::make_unique<CKKSNarrowBootstrapKey<LeanTarget>>();lean_key->selector_ntt_chunk_size=4;
    lean_key->automorphism_evaluator=[](auto&out,const auto&in,uint32_t,uint32_t exponent){for(std::size_t c=0;c<2;++c)Automorphism<P>(out[c],in[c],exponent);};
    for(bool tail:{false,true}){
        const auto expected=lean_reference.evaluate(tables,tail,0,*lean_key,nullptr,nullptr);
        const auto actual=lean.evaluate(tables,tail,0,*lean_key,nullptr,nullptr);
        for(std::size_t c=0;c<tables.size();++c)if(actual[c].ct!=expected[c].ct)throw std::runtime_error("sixteen-prime projected factors differ from independent digit CRT");
    }
    auto invalid=tables;
    std::fill(invalid[0].begin(),invalid[0].end(),std::complex<double>{8,0});
    rejected=false;
    try { direct.evaluate(invalid,false,0,*key,nullptr,nullptr); }
    catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected)throw std::runtime_error("unproved table coefficient bound accepted");
    std::cout<<"direct RNS and digit-CRT cached head/tail factors byte-identical; no online reads or selector transforms; retained bytes="<<DirectCache::retained_bytes<<'\n';
}
