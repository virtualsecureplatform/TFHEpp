#pragma once
#include "prepared_dd_ntt.hpp"
#include "key_projection.hpp"
#include "tensor_rns32.hpp"
#include <tuple>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Prepare every product-tree level from one saved public DD key. Projection
// changes the public representation only; this constructor never encrypts.
template <class P, std::uint32_t StartDataQ, std::uint32_t LogDelta,
          std::size_t Levels, std::uint32_t AuxiliaryLogScale = 100,
          std::uint32_t Bg = 64, std::uint32_t Bbar = 32, bool RNSTensor = false, bool WideRNS = false>
class CKKSPreparedAuxiliaryDDRelinKeyChainNTT {
    static_assert(Levels > 0 && LogDelta > 0);
    static_assert(StartDataQ > (Levels - 1) * LogDelta);
    template <std::size_t... I>
    static auto type(std::index_sequence<I...>) -> std::tuple<
        CKKSAuxiliaryDDRelinKeyNTT<P, StartDataQ - I * LogDelta,
            AuxiliaryLogScale, Bg, Bbar, WideRNS>...>;
    using Keys = decltype(type(std::make_index_sequence<Levels>{}));
    Keys keys;

    template<std::uint32_t DataQ,std::uint32_t MasterQ,std::uint32_t MasterBbar>
    static auto prepareLevel(const CKKSCompactDDRelinKey<P,MasterQ,Bg,MasterBbar>&master){
        using Key=CKKSAuxiliaryDDRelinKeyNTT<P,DataQ,AuxiliaryLogScale,Bg,Bbar,WideRNS>;
        auto projected=makeCKKSCompactDDRelinKeyAtLevel<DataQ+AuxiliaryLogScale>(master);
        if constexpr(WideRNS)return Key(std::make_shared<typename Key::Cache>(*projected));
        else return Key(makeCKKSDDRelinKeyNTT<Bbar>(*projected));
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
    std::size_t byteSize() const
    {
        return std::apply([](const auto &...key) {
            return (key.cache->byteSize() + ...);
        }, keys);
    }
};
} // namespace TFHEpp
#endif
