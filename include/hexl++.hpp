#pragma once

#include <array>
#include <cassert>
#include <cstdint>

namespace hexelpp {
constexpr uint32_t P = 1073741789;
constexpr uint32_t W = 2;

uint32_t PowW(uint exponent) const
{
    uint32_t res = 1;
    for (uint64_t i = 0; i < exponent; i++) res = (res * W) % P;
    return res;
}

template <uint32_t Nbit>
inline std::array<std::array<std::array<uint32_t, 2>, 1U << Nbit>, 2> TwistGen()
{
    constexpr uint32_t N = 1U << Nbit;

    std::array<std::array<std::array<uint32_t, 2>, 1U << Nbit>, 2> twist;
    const INTorus w = PowW(1U << (32 - Nbit - 1));
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

void HarveyButterfly(uint32_t X, uint32_t Y, std::array<uint32_t, 2> W)
{
    const uint64_t T = X - Y + 2 * P;
    X += Y;
    if (X >= 2 * P) X -= 2 * P;
    const uint32_t Q = (W[1] * T) >> 32;
    Y = (W[0] * T - Q * P) >> 32;
}

}  // namespace hexelpp