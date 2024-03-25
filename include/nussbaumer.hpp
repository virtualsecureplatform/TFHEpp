#pragma once
#include <span>

namespace Nussbaumer {

template <typename T, uint rbit>
inline void PolynomialMulByXai(const std::span<T, 1ull << rbit> res,
                               const size_t a)
{
    if (a == 0)
        return;
    else {
        constexpr size_t r = 1ull << rbit;
        std::array<T, r> temp;
        std::copy(res.begin(), res.end(), temp.begin());
        if (a < r) {
            for (int i = 0; i < a; i++) res[i] = -temp[i - a + r];
            for (int i = a; i < r; i++) res[i] = temp[i - a];
        }
        else {
            const size_t aa = a - r;
            for (int i = 0; i < aa; i++) res[i] = temp[i - aa + r];
            for (int i = aa; i < r; i++) res[i] = -temp[i - aa];
        }
    }
}

template <typename T, uint mbit, uint rbit>
void NussbaumerButterfly(const std::span<T, (1u << (rbit + mbit))> res)
{
    constexpr size_t m = 1ull << mbit;
    constexpr size_t r = 1ull << rbit;
    for (int i = 0; i < m / 2; i++)
        for (int j = 0; j < r; j++) {
            const T temp = res[i * r + j];
            res[i * r + j] += res[(i + m / 2) * r + j];
            res[(i + m / 2) * r + j] = temp - res[(i + m / 2) * r + j];
        }
    if constexpr (mbit != 1) {
        constexpr size_t stride = 1ull << (rbit - mbit);
        for (int i = 1; i < m / 2; i++)
            PolynomialMulByXai<T, rbit>(
                static_cast<std::span<T, r>>(res.subspan((i + m / 2) * r, r)),
                i * stride);
        NussbaumerButterfly<T, mbit - 1, rbit>(
            res.template subspan<0, m * r / 2>());
        NussbaumerButterfly<T, mbit - 1, rbit>(
            res.template subspan<m * r / 2, m * r / 2>());
    }
}

template <typename T, uint Nbit>
void NussbaumerTransform(std::span<T, (1ull << Nbit)> res)
{
    if constexpr (Nbit == 1) {
        const T temp = res[0];
        res[0] += res[1];
        res[1] = temp - res[1];
        return;
    }
    else {
        // initialize
        constexpr uint mbit = Nbit / 2;
        constexpr size_t m = 1ull << mbit;
        constexpr uint rbit = Nbit - mbit;
        constexpr size_t r = 1ull << rbit;
        std::array<T, (1ull << Nbit)> temp;
        std::copy(res.begin(), res.end(), temp.begin());
        // reorder
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < r; j++) res[i * r + j] = temp[m * j + i];
        }
        NussbaumerButterfly<T, mbit, rbit>(res);
        for (int i = 0; i < m; i++)
            NussbaumerTransform<T, rbit>(
                static_cast<std::span<T, r>>(res.subspan(i * r, r)));
    }
}

template <typename T, uint mbit, uint rbit>
void InverseNussbaumerButterfly(const std::span<T, (1u << (rbit + mbit))> res)
{
    constexpr size_t m = 1ull << mbit;
    constexpr size_t r = 1ull << rbit;
    if constexpr (mbit != 1) {
        constexpr size_t stride = 1ull << (rbit - mbit);
        InverseNussbaumerButterfly<T, mbit - 1, rbit>(
            res.template subspan<0, m * r / 2>());
        InverseNussbaumerButterfly<T, mbit - 1, rbit>(
            res.template subspan<m * r / 2, m * r / 2>());
        for (int i = 1; i < m / 2; i++)
            PolynomialMulByXai<T, rbit>(
                static_cast<std::span<T, r>>(res.subspan((i + m / 2) * r, r)),
                2 * r - i * stride);
    }
    for (int i = 0; i < m / 2; i++)
        for (int j = 0; j < r; j++) {
            const T temp = res[i * r + j];
            res[i * r + j] += res[(i + m / 2) * r + j];
            res[(i + m / 2) * r + j] = temp - res[(i + m / 2) * r + j];
        }
}

template <typename T, uint Nbit>
void InverseNussbaumerTransform(std::span<T, (1ull << Nbit)> res)
{
    if constexpr (Nbit == 1) {
        const T temp = res[0];
        res[0] += res[1];
        res[1] = temp - res[1];
        return;
    }
    else {
        // initialize
        constexpr uint mbit = Nbit / 2;
        constexpr size_t m = 1ull << mbit;
        constexpr uint rbit = Nbit - mbit;
        constexpr size_t r = 1ull << rbit;
        for (int i = 0; i < m; i++)
            InverseNussbaumerTransform<T, rbit>(
                static_cast<std::span<T, r>>(res.subspan(i * r, r)));
        InverseNussbaumerButterfly<T, mbit, rbit>(res);
        std::array<T, (1ull << Nbit)> temp;
        std::copy(res.begin(), res.end(), temp.begin());
        // reorder
        for (int i = 0; i < m; i++)
            for (int j = 0; j < r; j++) res[m * j + i] = temp[i * r + j];
    }
}
}  // namespace Nussbaumer