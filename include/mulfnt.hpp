#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "params.hpp"

namespace TFHEpp {

template <class P>
inline constexpr std::uint64_t FNTMaxAccumulation()
{
    constexpr std::uint64_t nonce_digit = std::uint64_t{1} << (P::Bgₐbit - 1);
    constexpr std::uint64_t main_digit = std::uint64_t{1} << (P::Bgbit - 1);
    constexpr std::uint64_t max_digit =
        nonce_digit > main_digit ? nonce_digit : main_digit;
    constexpr std::uint64_t max_chunk =
        (std::uint64_t{1} << FNTChunkBits<P>) - 1;
    return (P::k * P::lₐ + P::l) * P::n * max_digit * max_chunk;
}

template <class P>
inline constexpr void FNTSupportCheck()
{
    static_assert(std::is_same_v<typename P::T, std::uint32_t>,
                  "FNT backend currently supports uint32_t torus parameters");
    static_assert(
        P::nbit <= 10,
        "FNT backend currently supports polynomial rings with nbit <= 10");
    static_assert(P::l̅ == 1 && P::l̅ₐ == 1,
                  "FNT backend currently supports standard decomposition only");
    static_assert(
        FNTMaxAccumulation<P>() < static_cast<std::uint64_t>(FNTpp::P / 2),
        "FNT chunk size is too large for exact accumulated products");
}

inline int64_t FNTAddMod(const int64_t a, const int64_t b)
{
    const int64_t res = a + b;
    return res >= FNTpp::P ? res - FNTpp::P : res;
}

inline int64_t FNTMulMod(const int64_t a, const int64_t b)
{
    return static_cast<int64_t>((static_cast<__int128_t>(a) * b) % FNTpp::P);
}

inline int64_t FNTSignedToMod(const int64_t value)
{
    return value < 0 ? value + FNTpp::P : value;
}

inline int64_t FNTModToSigned(const int64_t value)
{
    return value > FNTpp::P / 2 ? value - FNTpp::P : value;
}

template <class P>
inline std::uint32_t FNTChunkMask(const std::uint32_t chunk)
{
    constexpr std::uint32_t width = std::numeric_limits<typename P::T>::digits;
    const std::uint32_t shift = chunk * FNTChunkBits<P>;
    const std::uint32_t bits = std::min(FNTChunkBits<P>, width - shift);
    return bits == width ? std::numeric_limits<std::uint32_t>::max()
                         : ((std::uint32_t{1} << bits) - 1);
}

template <class P>
inline void TwistFNTDigit(PolynomialFNT<P> &res, const Polynomial<P> &poly)
{
    FNTSupportCheck<P>();

    std::array<int64_t, P::n> converted;
    for (std::uint32_t i = 0; i < P::n; i++)
        converted[i] = FNTSignedToMod(static_cast<int32_t>(poly[i]));
    FNTpp::TwistFNT<P::nbit>(res, converted);
}

template <class P>
inline void TwistFNTTorus(PolynomialFNTChunks<P> &res,
                          const Polynomial<P> &poly)
{
    FNTSupportCheck<P>();

    for (std::uint32_t chunk = 0; chunk < FNTChunkCount<P>; chunk++) {
        std::array<int64_t, P::n> converted;
        const std::uint32_t shift = chunk * FNTChunkBits<P>;
        const std::uint32_t mask = FNTChunkMask<P>(chunk);
        for (std::uint32_t i = 0; i < P::n; i++)
            converted[i] = (poly[i] >> shift) & mask;
        FNTpp::TwistFNT<P::nbit>(res[chunk], converted);
    }
}

template <class P>
inline void FNTMulAdd(PolynomialFNT<P> &acc, const PolynomialFNT<P> &a,
                      const PolynomialFNT<P> &b)
{
    FNTSupportCheck<P>();

    for (std::uint32_t i = 0; i < FNTTransformedSize<P>; i++)
        acc[i] = FNTAddMod(acc[i], FNTMulMod(a[i], b[i]));
}

template <class P>
inline void TwistFNTToTorus(Polynomial<P> &res,
                            const PolynomialFNTChunks<P> &chunks)
{
    FNTSupportCheck<P>();

    res = {};
    for (std::uint32_t chunk = 0; chunk < FNTChunkCount<P>; chunk++) {
        std::array<int64_t, P::n> converted;
        FNTpp::TwistIFNT<P::nbit>(converted, chunks[chunk]);

        const std::uint32_t shift = chunk * FNTChunkBits<P>;
        for (std::uint32_t i = 0; i < P::n; i++)
            res[i] += static_cast<typename P::T>(FNTModToSigned(converted[i]))
                      << shift;
    }
}

}  // namespace TFHEpp
