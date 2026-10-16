#pragma once
#include "cmuxmem.hpp"
#include "circuitbootstrapping.hpp"
#include "gate.hpp"
#include <chrono>

namespace TFHEpp {
// Read-before-write RAM. A refresh backend receives ALL updated cell TLWEs
// and returns the same number under the same key and +/-mu convention.
// Coefficient zero alone is meaningful in each TRLWE cell. Inverse sample
// extraction therefore suffices to install a CKKS-refreshed TLWE; no extra
// TFHE bootstrap is hidden in the adapter.
template <class KSP, class CBP, class PKSP, class BRP,
          std::size_t AddressBits, std::size_t WordBits>
class RefreshableRAM {
public:
    using P=typename PKSP::targetP;
    static_assert(AddressBits>=2 && WordBits>0);
    static constexpr std::size_t cells=std::size_t{1}<<AddressBits;
    using Word=std::array<TLWE<P>,WordBits>;
    using Address=std::array<TLWE<P>,AddressBits>;
    using Memory=std::array<std::array<TRLWE<P>,cells>,WordBits>;
    struct Timing {double address=0,read=0,write=0,refresh=0;};
    Memory memory{};
    template <class Refresh>
    Word access(const Address &address,const TLWE<P> &write_enable,
                const Word &write_data,const EvalKey &ek,Refresh refresh,
                Timing *timing=nullptr)
    {
        using Clock=std::chrono::steady_clock;
        auto then=Clock::now();Timing times;
        auto selectors=std::make_unique<std::array<std::array<TRGSWFFT<P>,AddressBits>,2>>();
        for(std::size_t i=0;i<AddressBits;++i)
            CircuitBootstrappingWithInv<KSP,CBP,PKSP>((*selectors)[1][i],(*selectors)[0][i],address[i],ek);
        auto now=Clock::now();times.address=std::chrono::duration<double>(now-then).count();then=now;
        Word result{};
        for(std::size_t w=0;w<WordBits;++w){
            TRLWE<P> selected;
            RAMUX<P,AddressBits>(selected,(*selectors)[0],memory[w]);
            SampleExtractIndex<P>(result[w],selected,0);
        }
        now=Clock::now();times.read=std::chrono::duration<double>(now-then).count();then=now;
        std::vector<TLWE<P>> updated(cells*WordBits);
        for(std::size_t w=0;w<WordBits;++w){
            TRLWE<P> value;
            HomMUXwoSE<KSP,BRP>(value,write_enable,write_data[w],result[w],ek);
            for(std::size_t i=0;i<cells;++i){
                TRLWE<P> cell;
                RAMwriteBar<P,AddressBits>(cell,value,memory[w][i],i,*selectors);
                SampleExtractIndex<P>(updated[w*cells+i],cell,0);
            }
        }
        now=Clock::now();times.write=std::chrono::duration<double>(now-then).count();then=now;
        const auto refreshed=refresh(updated);
        if(refreshed.size()!=updated.size())throw std::runtime_error("RAM refresh returned wrong cell count");
        for(std::size_t w=0;w<WordBits;++w)for(std::size_t i=0;i<cells;++i)
            InvSampleExtractIndex<P>(memory[w][i],refreshed[w*cells+i],0);
        times.refresh=std::chrono::duration<double>(Clock::now()-then).count();
        if(timing)*timing=times;
        return result;
    }
};
} // namespace TFHEpp
