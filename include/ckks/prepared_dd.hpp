#pragma once

#include "ckks.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
namespace ckks_detail {
template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t Bgbit, std::uint32_t Bbarbit, class Seq>
struct PreparedDDChainTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t Bgbit, std::uint32_t Bbarbit, std::size_t... Is>
struct PreparedDDChainTuple<P, StartLogQ, LogDelta, Bgbit, Bbarbit,
                            std::index_sequence<Is...>> {
    using type = std::tuple<std::unique_ptr<CKKSDDRelinKeyFFT<P,
        StartLogQ - (Is + 1) * LogDelta, Bgbit, Bbarbit>>...>;
};
}  // namespace ckks_detail

// Derived public-key cache. Raw coefficient keys remain the serialized form.
template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t Depth, std::uint32_t Bgbit = P::Bgbit,
          std::uint32_t Bbarbit = P::B̅gbit>
struct CKKSPreparedDDRelinKeyChainFFT {
    using Tuple = typename ckks_detail::PreparedDDChainTuple<P, StartLogQ,
        LogDelta, Bgbit, Bbarbit, std::make_index_sequence<Depth>>::type;
    Tuple keys;
    template <std::size_t I>
    const auto &get() const { return *std::get<I>(keys); }
    std::size_t byteSize() const
    {
        return std::apply([](const auto &...key) {
            return (std::size_t{} + ... + sizeof(*key));
        }, keys);
    }
};

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t Depth, std::uint32_t Bgbit, std::uint32_t Bbarbit>
inline auto makeCKKSPreparedDDRelinKeyChainFFT(
    const CKKSDDRelinKeyChain<P, StartLogQ, LogDelta, Depth, Bgbit, Bbarbit> &raw)
{
    static_assert(use_multilimb_digit_fft_v<P>);
    // The FFT switch accumulates every primary-row convolution before its
    // inverse transform, so include that accumulation in the headroom check.
    static_assert(Bgbit + Bbarbit + P::nbit +
        std::bit_width(CKKSDDRelinPrimaryRowCountForLevel<P,
            StartLogQ - LogDelta, Bgbit>()) + 3 < std::numeric_limits<double>::digits);
    using Cache = CKKSPreparedDDRelinKeyChainFFT<P, StartLogQ, LogDelta,
        Depth, Bgbit, Bbarbit>;
    auto result = std::make_unique<Cache>();
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        ((std::get<Is>(result->keys) = makeCKKSDDRelinKeyFFT(
            raw.template get<Is>())), ...);
    }(std::make_index_sequence<Depth>{});
    return result;
}
}  // namespace TFHEpp
#endif
