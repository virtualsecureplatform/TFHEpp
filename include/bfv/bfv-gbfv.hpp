#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "tfhe/evalkeygens.hpp"
#include "tfhe/keyswitch.hpp"
#include "params.hpp"
#include "tfhe/trlwe.hpp"
#include "utils.hpp"

namespace TFHEpp::bfvboot::gbfv {

template <class P>
using TraceKey = std::vector<EvalAutoKey<P>>;

template <class P>
constexpr uint32_t TraceExponent(const int index)
{
    assert(index > 0 && index < static_cast<int>(P::nbit));
    return (P::n >> index) + 1;
}

template <class P>
void TraceKeyGen(TraceKey<P> &trace_key, const Key<P> &key,
                 const int start_index, const int end_index)
{
    assert(start_index > 0);
    assert(start_index <= end_index);
    assert(end_index < static_cast<int>(P::nbit));

    trace_key.resize(end_index - start_index + 1);
    for (int index = start_index; index <= end_index; index++)
        evalautokeygen<P>(trace_key[index - start_index],
                          TraceExponent<P>(index), key);
}

template <class P>
void EvaluateTrace(TRLWE<P> &res, const TRLWE<P> &ct,
                   const TraceKey<P> &trace_key, const int start_index,
                   const int end_index)
{
    assert(start_index > 0);
    assert(start_index <= end_index);
    assert(end_index < static_cast<int>(P::nbit));
    assert(trace_key.size() >=
           static_cast<std::size_t>(end_index - start_index + 1));

    res = ct;
    auto auto_ct = std::make_unique<TRLWE<P>>();
    for (int index = start_index; index <= end_index; index++) {
        EvalAuto<P>(*auto_ct, res, TraceExponent<P>(index),
                    trace_key[index - start_index]);
        TRLWEAdd<P>(res, res, *auto_ct);
    }
}

template <class P>
void EvaluateTrace(TRLWE<P> &res, const TRLWE<P> &ct,
                   const TraceKey<P> &trace_key, const int start_index)
{
    assert(!trace_key.empty());
    EvaluateTrace<P>(
        res, ct, trace_key, start_index,
        start_index + static_cast<int>(trace_key.size()) - 1);
}

template <class P>
void EvaluateTraceTorus(Polynomial<P> &res, const Polynomial<P> &poly,
                        const int start_index, const int end_index)
{
    assert(start_index > 0);
    assert(start_index <= end_index);
    assert(end_index < static_cast<int>(P::nbit));

    res = poly;
    auto aut = std::make_unique<Polynomial<P>>();
    for (int index = start_index; index <= end_index; index++) {
        Automorphism<P>(*aut, res, TraceExponent<P>(index));
        for (uint32_t i = 0; i < P::n; i++) res[i] += (*aut)[i];
    }
}

template <class P>
void AutomorphismPlainMod(std::array<uint64_t, P::n> &res,
                          const std::array<uint64_t, P::n> &poly,
                          const uint32_t d)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr uint32_t nmask = (uint32_t{1} << P::nbit) - 1;
    constexpr uint32_t signmask = uint32_t{1} << P::nbit;

    res = {};
    for (uint32_t i = 0; i < P::n; i++) {
        const uint32_t index = i * d;
        const uint32_t dst = index & nmask;
        const uint64_t coeff = poly[i] % t;
        if ((index & signmask) == 0) {
            res[dst] += coeff;
            res[dst] %= t;
        }
        else if (coeff != 0) {
            res[dst] += t - coeff;
            res[dst] %= t;
        }
    }
}

template <class P>
void EvaluateTracePlainMod(std::array<uint64_t, P::n> &res,
                           const std::array<uint64_t, P::n> &poly,
                           const int start_index, const int end_index)
{
    assert(start_index > 0);
    assert(start_index <= end_index);
    assert(end_index < static_cast<int>(P::nbit));

    res = poly;
    std::array<uint64_t, P::n> aut{};
    for (int index = start_index; index <= end_index; index++) {
        AutomorphismPlainMod<P>(aut, res, TraceExponent<P>(index));
        for (uint32_t i = 0; i < P::n; i++)
            res[i] = (res[i] + aut[i]) %
                     static_cast<uint64_t>(P::plain_modulus);
    }
}

namespace detail {

template <class T>
inline bool centered_negative(const T x)
{
    if constexpr (is_multilimb_uint_v<T>)
        return x.bit(std::numeric_limits<T>::digits - 1);
    else
        return ((x >> (std::numeric_limits<T>::digits - 1)) & T{1}) != T{0};
}

inline bool round_remainder_up(const uint64_t rem, const uint64_t divisor)
{
    return rem > divisor / 2 || (divisor % 2 == 0 && rem == divisor / 2);
}

template <class T>
inline T round_div_unsigned_magnitude(const T x, const uint64_t divisor)
{
    assert(divisor != 0);
    if (divisor == 1) return x;

    T quotient = x / divisor;
    uint64_t rem;
    if constexpr (is_multilimb_uint_v<T>)
        rem = x % divisor;
    else
        rem = static_cast<uint64_t>(x % divisor);

    if (round_remainder_up(rem, divisor)) quotient += T{1};
    return quotient;
}

template <class T>
inline T round_div_centered(const T x, const uint64_t divisor)
{
    if (!centered_negative<T>(x))
        return round_div_unsigned_magnitude<T>(x, divisor);

    const T magnitude = T{0} - x;
    const T quotient = round_div_unsigned_magnitude<T>(magnitude, divisor);
    return T{0} - quotient;
}

}  // namespace detail

template <class FromP, class ToP>
void ScaleAndRoundPlainMod(TRLWE<ToP> &res, const TRLWE<FromP> &ct)
{
    static_assert(FromP::n == ToP::n && FromP::k == ToP::k,
                  "plain modulus switching requires matching dimensions");
    static_assert(std::is_same_v<typename FromP::T, typename ToP::T>,
                  "plain modulus switching requires matching torus types");

    constexpr uint64_t old_t = static_cast<uint64_t>(FromP::plain_modulus);
    constexpr uint64_t new_t = static_cast<uint64_t>(ToP::plain_modulus);
    using T = typename ToP::T;

    // This mirrors the reference ScaleAndRoundCiphertext coefficient operation.
    // When increasing the plaintext modulus, an encrypted ciphertext preserves
    // the old plaintext residue class but cannot define a canonical high
    // p-adic digit: component wrap carries are hidden by the mask. Transparent
    // ciphertexts get the exact representative; encrypted ciphertexts are only
    // guaranteed to round-trip after scaling back down.
    if constexpr (old_t == new_t) {
        for (int c = 0; c <= static_cast<int>(ToP::k); c++)
            for (uint32_t i = 0; i < ToP::n; i++) res[c][i] = ct[c][i];
    }
    else if constexpr (old_t > new_t && old_t % new_t == 0) {
        constexpr uint64_t ratio = old_t / new_t;
        for (int c = 0; c <= static_cast<int>(ToP::k); c++)
            for (uint32_t i = 0; i < ToP::n; i++) res[c][i] = ct[c][i] * ratio;
    }
    else if constexpr (new_t > old_t && new_t % old_t == 0) {
        constexpr uint64_t divisor = new_t / old_t;
        for (int c = 0; c <= static_cast<int>(ToP::k); c++)
            for (uint32_t i = 0; i < ToP::n; i++)
                res[c][i] =
                    detail::round_div_centered<T>(ct[c][i], divisor);
    }
    else {
        static_assert(false_v<FromP>,
                      "unsupported non-integral plaintext modulus switch");
    }
}

}  // namespace TFHEpp::bfvboot::gbfv
