#pragma once

#include <array>
#include <cstdint>
#include <type_traits>

#include "params.hpp"

namespace TFHEpp {

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
inline void TwistFNTTorus(PolynomialFNT<P> &res, const Polynomial<P> &poly)
{
    FNTSupportCheck<P>();

    std::array<int64_t, P::n> converted;
    for (std::uint32_t i = 0; i < P::n; i++) converted[i] = poly[i];
    FNTpp::TwistFNT<P::nbit>(res, converted);
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
inline void TwistFNTToTorus(Polynomial<P> &res, const PolynomialFNT<P> &fnt)
{
    FNTSupportCheck<P>();

    std::array<int64_t, P::n> converted;
    FNTpp::TwistIFNT<P::nbit>(converted, fnt);
    // Reduce mod 2^32 by narrowing. In particular, the FNT residue 2^32 maps
    // back to torus zero.
    for (std::uint32_t i = 0; i < P::n; i++)
        res[i] = static_cast<typename P::T>(converted[i]);
}

}  // namespace TFHEpp
