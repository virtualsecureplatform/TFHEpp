#pragma once

#include <bitset>

#include "detwfa.hpp"

namespace TFHEpp {
template <class P, uint32_t address_bit>
void RAMUX(TRLWE<P> &res, const array<TRGSWFFT<P>, address_bit> &invaddress,
           const array<TRLWE<P>, 1 << address_bit> &data)
{
    constexpr uint32_t num_trlwe = 1 << address_bit;
    array<TRLWE<P>, num_trlwe / 2> temp;

    for (uint32_t index = 0; index < num_trlwe / 2; index++) {
        CMUXFFT<P>(temp[index], invaddress[0], data[2 * index],
                   data[2 * index + 1]);
    }

    for (uint32_t bit = 0; bit < (address_bit - 2); bit++) {
        const uint32_t stride = 1 << bit;
        for (uint32_t index = 0; index < (num_trlwe >> (bit + 2)); index++) {
            CMUXFFT<P>(temp[(2 * index) * stride], invaddress[bit + 1],
                       temp[(2 * index) * stride],
                       temp[(2 * index + 1) * stride]);
        }
    }
    constexpr uint32_t stride = 1 << (address_bit - 2);
    CMUXFFT<P>(res, invaddress[address_bit - 1], temp[0], temp[stride]);
}

template <class P, uint32_t address_bit>
void RAMwriteBar(
    TRLWE<P> &res, const TRLWE<P> &writed, const TRLWE<P> &encmem,
    const int index,
    const std::array<std::array<TRGSWFFT<P>, address_bit>, 2> &bootedTRGSW)
{
    const std::bitset<address_bit> addressbitset(index);
    res = writed;
    for (int i = 0; i < address_bit; i++)
        CMUXFFT<P>(res, bootedTRGSW[addressbitset[i]][i], res, encmem);
}

// MUX tree for vertical packing
template <class P, uint32_t address_bit, uint32_t width_bit>
void UROMUX(TRLWE<P> &res, const array<TRGSWFFT<P>, address_bit> &invaddress,
            const array<TRLWE<P>, 1 << (address_bit - width_bit)> &data)
{
    constexpr uint32_t Ubit = address_bit - width_bit;
    constexpr uint32_t num_trlwe = 1 << (Ubit);
    array<TRLWE<P>, num_trlwe / 2> temp;

    for (int index = 0; index < num_trlwe / 2; index++) {
        CMUXFFT<P>(temp[index], invaddress[width_bit], data[2 * index],
                   data[2 * index + 1]);
    }

    for (int bit = 0; bit < (Ubit - 2); bit++) {
        const uint32_t stride = 1 << bit;
        for (uint32_t index = 0; index < (num_trlwe >> (bit + 2)); index++) {
            CMUXFFT<P>(
                temp[(2 * index) * stride], invaddress[width_bit + bit + 1],
                temp[(2 * index) * stride], temp[(2 * index + 1) * stride]);
        }
    }

    constexpr uint32_t stride = 1 << (Ubit - 2);
    CMUXFFT<P>(res, invaddress[address_bit - 1], temp[0], temp[stride]);
}

// MUX tree for horizontal packing
template <class P, uint32_t address_bit, uint32_t width_bit>
void LROMUX(vector<TLWE<P>> &res,
            const array<TRGSWFFT<P>, address_bit> &address,
            const TRLWE<P> &data)
{
    TRLWE<P> temp, acc;
    PolynomialMulByXaiMinusOne<P>(temp[0], data[0], 2 * P::n - (P::n >> 1));
    PolynomialMulByXaiMinusOne<P>(temp[1], data[1], 2 * P::n - (P::n >> 1));
    trgswfftExternalProduct<P>(temp, temp, address[width_bit - 1]);
    for (int i = 0; i < P::n; i++) {
        // initialize acc
        acc[0][i] = temp[0][i] + data[0][i];
        acc[1][i] = temp[1][i] + data[1][i];
    }

    for (uint32_t bit = 2; bit <= width_bit; bit++) {
        PolynomialMulByXaiMinusOne<P>(temp[0], acc[0],
                                      2 * P::n - (P::n >> bit));
        PolynomialMulByXaiMinusOne<P>(temp[1], acc[1],
                                      2 * P::n - (P::n >> bit));
        trgswfftExternalProduct<P>(temp, temp, address[width_bit - bit]);
        for (int i = 0; i < P::n; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }

    constexpr uint32_t word = 1 << (P::nbit - width_bit);
    for (int i = 0; i < word; i++) SampleExtractIndex<P>(res[i], acc, i);
}
}  // namespace TFHEpp