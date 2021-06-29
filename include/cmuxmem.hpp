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

template <class ksP, uint32_t address_bit, uint32_t width_bit>
void LROMUX(vector<TLWE<typename ksP::targetP>> &res,
            const array<TRGSWFFT<typename ksP::domainP>, address_bit> &address,
            const TRLWE<typename ksP::domainP> &data,
            const KeySwitchingKey<ksP> &ksk)
{
    TRLWE<typename ksP::domainP> temp, acc;
    PolynomialMulByXaiMinusOne<typename ksP::domainP>(
        temp[0], data[0], 2 * ksP::domainP::n - (ksP::domainP::n >> 1));
    PolynomialMulByXaiMinusOne<typename ksP::domainP>(
        temp[1], data[1], 2 * ksP::domainP::n - (ksP::domainP::n >> 1));
    trgswfftExternalProduct<typename ksP::domainP>(temp, temp,
                                                   address[width_bit - 1]);
    for (int i = 0; i < ksP::domainP::n; i++) {
        // initialize acc
        acc[0][i] = temp[0][i] + data[0][i];
        acc[1][i] = temp[1][i] + data[1][i];
    }

    for (uint32_t bit = 2; bit <= width_bit; bit++) {
        PolynomialMulByXaiMinusOne<typename ksP::domainP>(
            temp[0], acc[0], 2 * ksP::domainP::n - (ksP::domainP::n >> bit));
        PolynomialMulByXaiMinusOne<typename ksP::domainP>(
            temp[1], acc[1], 2 * ksP::domainP::n - (ksP::domainP::n >> bit));
        trgswfftExternalProduct<typename ksP::domainP>(
            temp, temp, address[width_bit - bit]);
        for (int i = 0; i < ksP::domainP::n; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }

    constexpr uint32_t word = 1 << (ksP::domainP::nbit - width_bit);
    array<TLWE<typename ksP::domainP>, word> reslvl1;
    for (int i = 0; i < word; i++)
        SampleExtractIndex<typename ksP::domainP>(reslvl1[i], acc, i);
    for (int i = 0; i < word; i++)
        IdentityKeySwitch<ksP>(res[i], reslvl1[i], ksk);
}
}  // namespace TFHEpp