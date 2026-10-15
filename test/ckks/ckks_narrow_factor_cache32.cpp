#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/narrow_factor_ntt32.hpp>
#include <iostream>

using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<9,7>;
using Source=CKKSNarrowBootstrapSchedule<P,48,440,57,42,31,false,0,true>;
using Target=CKKSNarrowBootstrapSchedule<P,48,427,44,42,31,true,13,true>;

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
int main()
{
    auto raw=std::make_shared<TransparentSource>();
    using Cache=CKKSNarrowFactorSelectorCache32<Source>;
    bool rejected=false;
    try { Cache insufficient(raw,Cache::retained_bytes-1,4); }
    catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected || raw->reads)throw std::runtime_error("cache budget did not reject before loading");
    auto cache=std::make_shared<Cache>(raw,Cache::retained_bytes,4);
    CKKSNarrowDeferredRescaleFactorSource32<Target,Source> streamed(raw),cached(raw,cache);
    auto key=std::make_unique<CKKSNarrowBootstrapKey<Target>>();
    key->selector_ntt_chunk_size=4;
    key->automorphism_evaluator=[](auto &result,const auto &input,uint32_t,uint32_t exponent) {
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
        if(raw->reads!=reads || timing.selector_load_seconds || timing.selector_ntt_seconds)
            throw std::runtime_error("cached evaluation loaded or transformed raw selectors");
        for(std::size_t j=0;j<tables.size();++j)
            if(reference[j].ct!=result[j].ct)throw std::runtime_error("cached factor differs from streamed factor");
    }
    std::cout<<"cached head/tail factors byte-identical; no online reads or selector transforms; retained bytes="<<Cache::retained_bytes<<'\n';
}
