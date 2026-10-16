#pragma once
#include "prepared_dd_ntt.hpp"
#include "key_projection.hpp"
#include "tensor_rns32.hpp"
#include "prepared_dd_batch_rns32.hpp"
#include <tuple>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Prepare every product-tree level from one saved public DD key. Projection
// changes the public representation only; this constructor never encrypts.
template <class P, std::uint32_t StartDataQ, std::uint32_t LogDelta,
          std::size_t Levels, std::uint32_t AuxiliaryLogScale = 100,
          std::uint32_t Bg = 64, std::uint32_t Bbar = 32, bool RNSTensor = false, bool WideRNS = false,
          class ChunkPolicy = CKKSUniformGadgetChunks<80>>
class CKKSPreparedAuxiliaryDDRelinKeyChainNTT {
    static_assert(Levels > 0 && LogDelta > 0);
    static_assert(StartDataQ > (Levels - 1) * LogDelta);
    template<std::uint32_t DataQ>
    using Key=CKKSAuxiliaryDDRelinKeyNTT<P,DataQ,AuxiliaryLogScale,Bg,Bbar,WideRNS,
        ChunkPolicy::template bits<P,DataQ+AuxiliaryLogScale,Bg>>;
    template <std::size_t... I>
    static auto type(std::index_sequence<I...>) -> std::tuple<
        Key<StartDataQ-I*LogDelta>...>;
    using Keys = decltype(type(std::make_index_sequence<Levels>{}));
    Keys keys;

    template<std::uint32_t DataQ,std::uint32_t MasterQ,std::uint32_t MasterBbar>
    static auto prepareLevel(const CKKSCompactDDRelinKey<P,MasterQ,Bg,MasterBbar>&master){
        using LevelKey=Key<DataQ>;
        auto projected=makeCKKSCompactDDRelinKeyAtLevel<DataQ+AuxiliaryLogScale>(master);
        if constexpr(WideRNS)return LevelKey(std::make_shared<typename LevelKey::Cache>(*projected));
        else return LevelKey(makeCKKSDDRelinKeyNTT<Bbar>(*projected));
    }
    template <std::uint32_t MasterQ, std::uint32_t MasterBbar, std::size_t... I>
    static Keys prepare(const CKKSCompactDDRelinKey<P, MasterQ, Bg, MasterBbar> &master,
                        std::index_sequence<I...>)
    {
        static_assert(StartDataQ + AuxiliaryLogScale <= MasterQ);
        return Keys{prepareLevel<StartDataQ-I*LogDelta>(master)...};
    }
public:
    template <std::uint32_t MasterQ, std::uint32_t MasterBbar>
    explicit CKKSPreparedAuxiliaryDDRelinKeyChainNTT(
        const CKKSCompactDDRelinKey<P, MasterQ, Bg, MasterBbar> &master)
        : keys(prepare(master, std::make_index_sequence<Levels>{})) {}

    template <std::size_t Level> decltype(auto) get() const
    {
        static_assert(Level < Levels);
        if constexpr(RNSTensor)
            return CKKSRelinKeyTensorRNS32<P,std::tuple_element_t<Level,Keys>>{std::get<Level>(keys)};
        else return (std::get<Level>(keys));
    }
    template<std::size_t Level,class Outputs,class Inputs>
    void multiplyBatch(Outputs&next,const Inputs&factors,
        const std::vector<std::pair<std::size_t,std::size_t>>&pairs)const
        requires(RNSTensor && WideRNS)
    {
        static_assert(Level<Levels);
        constexpr auto Q=StartDataQ-Level*LogDelta;
        std::vector<CKKSRNS32Product<P,Q,LogDelta>> work;
        work.reserve(pairs.size());
        for(const auto&[tree,pair]:pairs)
            work.push_back({&next[tree][pair],&factors[tree][2*pair],&factors[tree][2*pair+1]});
        CKKSMultRelinThenRescaleBatchRNS32<P,Q,LogDelta>(work,std::get<Level>(keys));
    }
    std::size_t byteSize() const
    {
        return std::apply([](const auto &...key) {
            return (key.cache->byteSize() + ...);
        }, keys);
    }
};
} // namespace TFHEpp
#endif
