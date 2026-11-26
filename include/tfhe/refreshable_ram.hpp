#pragma once
#include "cmuxmem.hpp"
#include "circuitbootstrapping.hpp"
#include "gate.hpp"
#include <chrono>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

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
    // Circuit outside the refresh. The default is the fastest one.
    //   parallel: OpenMP over address bits, read-tree nodes and cells.
    //   demux:    route (new value - read value) to the addressed cell through
    //             a binary tree of external products, cells-1 products in
    //             place of cells*AddressBits CMUXes, and add it to every cell.
    // The original serial per-cell CMUX chain stays available: set both to
    // false, or run with TFHEPP_RAM_CIRCUIT=serial (or parallel, for the
    // chain under OpenMP) to change the default of every RAM in a program.
    struct Options {
        bool parallel=true,demux=true;
        static Options fromEnvironment()
        {
            const char *value=std::getenv("TFHEPP_RAM_CIRCUIT");
            const std::string name=value?value:"demux";
            if(name=="serial")return {false,false};
            if(name=="parallel")return {true,false};
            if(name=="demux")return {true,true};
            throw std::invalid_argument("TFHEPP_RAM_CIRCUIT: serial, parallel or demux");
        }
    };
    Options options=Options::fromEnvironment();
    Memory memory{};
    template <class Refresh>
    Word access(const Address &address,const TLWE<P> &write_enable,
                const Word &write_data,const EvalKey &ek,Refresh refresh,
                Timing *timing=nullptr)
    {
        using Clock=std::chrono::steady_clock;
        auto then=Clock::now();Timing times;
        const bool parallel=options.parallel;
        auto selectors=std::make_unique<std::array<std::array<TRGSWFFT<P>,AddressBits>,2>>();
#pragma omp parallel for schedule(dynamic) if(parallel)
        for(std::size_t i=0;i<AddressBits;++i)
            CircuitBootstrappingWithInv<KSP,CBP,PKSP>((*selectors)[1][i],(*selectors)[0][i],address[i],ek);
        auto now=Clock::now();times.address=std::chrono::duration<double>(now-then).count();then=now;
        Word result{};
        std::vector<TRLWE<P>> selected(WordBits);
        for(std::size_t w=0;w<WordBits;++w){
            if(parallel){
                std::vector<TRLWE<P>> tree(cells/2);
                for(std::size_t bit=0;bit<AddressBits;++bit){
                    const std::size_t nodes=cells>>(bit+1),stride=std::size_t{1}<<bit;
#pragma omp parallel for schedule(static) if(nodes>=64)
                    for(std::size_t i=0;i<nodes;++i){
                        if(bit==0)CMUXFFT<P>(tree[i],(*selectors)[0][0],memory[w][2*i],memory[w][2*i+1]);
                        else CMUXFFT<P>(tree[2*i*stride/2],(*selectors)[0][bit],
                            tree[2*i*stride/2],tree[(2*i+1)*stride/2]);
                    }
                }
                selected[w]=tree[0];
            }else RAMUX<P,AddressBits>(selected[w],(*selectors)[0],memory[w]);
            SampleExtractIndex<P>(result[w],selected[w],0);
        }
        now=Clock::now();times.read=std::chrono::duration<double>(now-then).count();then=now;
        std::vector<TLWE<P>> updated(cells*WordBits);
        for(std::size_t w=0;w<WordBits;++w){
            TRLWE<P> value;
            HomMUXwoSE<KSP,BRP>(value,write_enable,write_data[w],result[w],ek);
            if(options.demux){
                // tree[i] ends as (i == address) * (value - read).
                std::vector<TRLWE<P>> tree(cells);
                for(std::size_t k=0;k<=P::k;++k)for(std::size_t j=0;j<P::n;++j)
                    tree[0][k][j]=value[k][j]-selected[w][k][j];
                for(std::size_t bit=AddressBits;bit-->0;){
                    const std::size_t stride=std::size_t{1}<<bit,nodes=cells>>(bit+1);
#pragma omp parallel for schedule(static) if(parallel&&nodes>=64)
                    for(std::size_t i=0;i<nodes;++i){
                        auto &zero=tree[2*i*stride];auto &one=tree[(2*i+1)*stride];
                        ExternalProduct<P>(one,zero,(*selectors)[1][bit]);
                        for(std::size_t k=0;k<=P::k;++k)for(std::size_t j=0;j<P::n;++j)
                            zero[k][j]-=one[k][j];
                    }
                }
#pragma omp parallel for schedule(static) if(parallel)
                for(std::size_t i=0;i<cells;++i){
                    for(std::size_t k=0;k<=P::k;++k)for(std::size_t j=0;j<P::n;++j)
                        tree[i][k][j]+=memory[w][i][k][j];
                    SampleExtractIndex<P>(updated[w*cells+i],tree[i],0);
                }
                continue;
            }
#pragma omp parallel for schedule(static) if(parallel)
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
