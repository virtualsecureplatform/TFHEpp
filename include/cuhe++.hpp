#pragma once

#include <bits/stdint-uintn.h>

#include <array>
#include <cassert>
#include <cstdint>

#include "INTorus.hpp"

namespace cuHEpp {

template <uint8_t bit>
uint32_t BitReverse(uint32_t in)
{
    if constexpr (bit > 1) {
        const uint32_t center = in & ((bit & 1) << (bit / 2));
        return (BitReverse<bit / 2>(in & ((1U << (bit / 2)) - 1))
                << (bit + 1) / 2) |
               center | BitReverse<bit / 2>(in >> ((bit + 1) / 2));
    }
    else {
        return in;
    }
}

// NTT implementation
// https://nufhe.readthedocs.io/en/latest/implementation_details.html
constexpr uint64_t W = 12037493425763644479ULL;

template <uint32_t Nbit>
inline std::array<std::array<INTorus, 1U << Nbit>, 2> TwistGen()
{
    constexpr uint32_t N = 1U << Nbit;

    std::array<std::array<INTorus, 1U << Nbit>, 2> twist;
    const INTorus w = INTorus(W).Pow(1U << (32 - Nbit - 1));
    twist[0][0] = twist[1][0] = INTorus(1, false);
    for (uint32_t i = 1; i < N; i++) twist[1][i] = twist[1][i - 1] * w;
    assert((twist[1][N - 1] * w).Pow(2).value == 1);
    twist[0][N - 1] = twist[1][N - 1] * w * w;
    for (uint32_t i = 2; i < N; i++) twist[0][N - i] = twist[0][N - i + 1] * w;
    assert((twist[0][1] * w).value == 1);
    return twist;
}

template <uint32_t Nbit>
inline std::array<std::array<INTorus, 1U << Nbit>, 2> TableGen()
{
    constexpr uint32_t N = 1U << Nbit;

    std::array<std::array<INTorus, N>, 2> table;
    const INTorus w = INTorus(W).Pow(1U << (32 - Nbit));
    table[0][0] = table[1][0] = INTorus(1, false);
    for (uint32_t i = 1; i < N; i++) table[1][i] = table[1][i - 1] * w;
    for (uint32_t i = 1; i < N; i++) table[0][i] = table[1][N - i];
    return table;
}

inline void ButterflyAdd(INTorus *const res, const uint32_t size)
{
    for (uint32_t index = 0; index < size / 2; index++) {
        const INTorus temp = res[index];
        res[index] += res[index + size / 2];
        res[index + size / 2] = temp - res[index + size / 2];
    }
}

template <uint32_t Nbit>
inline void TwiddleMul(INTorus *const res, const uint32_t size,
                       const uint32_t stride,
                       const std::array<INTorus, 1 << Nbit> &table)
{
    for (uint32_t index = 1; index < size; index++)
        res[index] *= table[stride * index];
}

template <uint8_t radixbit>
inline void INTTradixButterfly(INTorus *const res, const uint32_t size)
{
    static_assert(radixbit <= 6, "radix 64 is the maximum!");
    if constexpr (radixbit != 0) {
        ButterflyAdd(res, size);
        const uint32_t block = size >> radixbit;
        for (int i = 1; i < (1 << (radixbit - 1)); i++)
            for (int j = 0; j < block; j++)
                res[i * block + j + size / 2] = res[i * block + j + size / 2]
                                                << (3 * (i << (6 - radixbit)));
        INTTradixButterfly<radixbit - 1>(&res[0], size / 2);
        INTTradixButterfly<radixbit - 1>(&res[size / 2], size / 2);
    }
}

template <uint32_t Nbit, uint8_t radixbit>
inline void INTTradix(INTorus *const res, const uint32_t size,
                      const uint32_t num_block,
                      const std::array<INTorus, 1 << Nbit> &table)
{
    INTTradixButterfly<radixbit>(res, size);
    for (uint32_t i = 1; i < (1 << radixbit); i++)
        TwiddleMul<Nbit>(&res[i * (size >> radixbit)], size >> radixbit,
                         BitReverse<radixbit>(i) * num_block, table);
}

template <uint32_t Nbit, uint8_t radixbit>
inline void INTT(std::array<INTorus, 1 << Nbit> &res,
                 const std::array<INTorus, 1 << Nbit> &table)
{
    for (uint8_t sizebit = Nbit; sizebit > radixbit; sizebit -= radixbit) {
        const uint32_t size = 1U << sizebit;
        const uint32_t num_block = 1U << (Nbit - sizebit);
        for (uint32_t block = 0; block < num_block; block++)
            INTTradix<Nbit, radixbit>(&res[size * block], size, num_block,
                                      table);
    }
    constexpr uint8_t remainder = ((Nbit - 1) % radixbit) + 1;
    constexpr uint32_t size = 1U << remainder;
    constexpr uint32_t num_block = 1U << (Nbit - remainder);
    for (uint32_t block = 0; block < num_block; block++)
        INTTradixButterfly<remainder>(&res[size * block], size);
}

template <typename T = uint32_t, uint32_t Nbit>
inline void TwistMulInvert(std::array<INTorus, 1 << Nbit> &res,
                           const std::array<T, 1 << Nbit> &a,
                           const std::array<INTorus, 1 << Nbit> &twist)
{
    constexpr uint32_t N = 1 << Nbit;
    if constexpr (std::is_same_v<T, uint64_t>) {
        for (int i = 0; i < N; i++)
            res[i] = INTorus(static_cast<uint64_t>(
                                 (static_cast<__uint128_t>(a[i]) * P) >> 64),
                             true) *
                     twist[i];
    }
    else {
        for (int i = 0; i < N; i++) res[i] = INTorus(a[i], false) * twist[i];
    }
}

template <typename T, uint32_t Nbit>
void TwistINTT(std::array<INTorus, 1 << Nbit> &res,
               const std::array<T, 1 << Nbit> &a,
               const std::array<INTorus, 1 << Nbit> &table,
               const std::array<INTorus, 1 << Nbit> &twist)
{
    TwistMulInvert<T, Nbit>(res, a, twist);
    INTT<Nbit, 6>(res, table);
}

template <uint8_t radixbit>
inline void NTTradixButterfly(INTorus *const res, const uint32_t size)
{
    static_assert(radixbit <= 6, "radix 64 is the maximum!");
    if constexpr (radixbit != 0) {
        NTTradixButterfly<radixbit - 1>(&res[size / 2], size / 2);
        NTTradixButterfly<radixbit - 1>(&res[0], size / 2);
        const uint32_t block = size >> radixbit;
        if constexpr (radixbit != 1)
            for (int i = 1; i < (1 << (radixbit - 1)); i++)
                for (int j = 0; j < block; j++)
                    res[i * block + j + size / 2] =
                        res[i * block + j + size / 2]
                        << (3 * (64 - (i << (6 - radixbit))));
        ButterflyAdd(res, size);
    }
}

template <uint32_t Nbit, uint8_t radixbit>
inline void NTTradix(INTorus *const res, const uint32_t size,
                     const uint32_t num_block,
                     const std::array<INTorus, 1 << Nbit> &table)
{
    for (uint32_t i = 1; i < (1 << radixbit); i++)
        TwiddleMul<Nbit>(&res[i * (size >> radixbit)], size >> radixbit,
                         BitReverse<radixbit>(i) * num_block, table);
    NTTradixButterfly<radixbit>(res, size);
}

template <uint32_t Nbit, uint8_t radixbit>
void NTT(std::array<INTorus, 1 << Nbit> &res,
         const std::array<INTorus, 1 << Nbit> &table)
{
    constexpr uint8_t remainder = ((Nbit - 1) % radixbit) + 1;
    constexpr uint32_t size = 1U << remainder;
    constexpr uint32_t num_block = 1U << (Nbit - remainder);
    for (uint32_t block = 0; block < num_block; block++)
        NTTradixButterfly<remainder>(&res[size * block], size);
    for (uint8_t sizebit = remainder + radixbit; sizebit <= Nbit;
         sizebit += radixbit) {
        const uint32_t size = 1U << sizebit;
        const uint32_t num_block = 1U << (Nbit - sizebit);
        for (uint32_t block = 0; block < num_block; block++)
            NTTradix<Nbit, radixbit>(&res[size * block], size, num_block,
                                     table);
    }
}

template <typename T = uint32_t, uint32_t Nbit>
inline void TwistMulDirect(std::array<T, 1 << Nbit> &res,
                           const std::array<INTorus, 1 << Nbit> &a,
                           const std::array<INTorus, 1 << Nbit> &twist)
{
    const INTorus invN = InvPow2(Nbit);
    constexpr uint32_t N = 1 << Nbit;
    if constexpr (std::is_same_v<T, uint64_t>) {
        for (int i = 0; i < N; i++)
            res[i] = static_cast<T>(
                (static_cast<__uint128_t>((a[i] * twist[i] * invN).value)
                 << 64) /
                P);
    }
    else {
        for (int i = 0; i < N; i++)
            res[i] = static_cast<T>((a[i] * twist[i] * invN).value);
    }
}

template <typename T, uint32_t Nbit>
void TwistNTT(std::array<T, 1 << Nbit> &res, std::array<INTorus, 1 << Nbit> &a,
              const std::array<INTorus, 1 << Nbit> &table,
              const std::array<INTorus, 1 << Nbit> &twist)
{
    NTT<Nbit, 6>(a, table);
    TwistMulDirect<T, Nbit>(res, a, twist);
}

template <typename T, uint32_t Nbit>
void PolyMullvl1(std::array<T, 1 << Nbit> &res, std::array<T, 1 << Nbit> &a,
                 std::array<T, 1 << Nbit> &b,
                 const std::array<std::array<INTorus, 1 << Nbit>, 2> &table,
                 const std::array<std::array<INTorus, 1 << Nbit>, 2> &twist)
{
    std::array<INTorus, 1 << Nbit> ntta, nttb;
    TwistINTT<T, Nbit>(ntta, a, table[1], twist[1]);
    TwistINTT<T, Nbit>(nttb, b, table[1], twist[1]);
    for (int i = 0; i < (1U << Nbit); i++) ntta[i] *= nttb[i];
    TwistNTT<T, Nbit>(res, ntta, table[0], twist[0]);
}
}  // namespace cuHEpp