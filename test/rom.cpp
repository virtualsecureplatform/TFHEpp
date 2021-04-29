#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

template <class P,uint32_t address_bit, uint32_t width_bit>
void UROMUX(TRLWE<P> &res,
            const array<TRGSWFFT<P>, address_bit> &invaddress,
            const array<TRLWE<P>, 1 << (address_bit - width_bit)> &data)
{
    constexpr uint32_t Ubit = address_bit - width_bit;
    constexpr uint32_t num_trlwe = 1 << (Ubit);
    array<TRLWE<P>, num_trlwe / 2> temp;

    for (uint32_t index = 0; index < num_trlwe / 2; index++) {
        CMUXFFT<P>(temp[index], invaddress[width_bit], data[2 * index],
                           data[2 * index + 1]);
    }

    for (uint32_t bit = 0; bit < (Ubit - 2); bit++) {
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

template <class ksP,uint32_t address_bit, uint32_t width_bit>
void LROMUX(vector<TLWE<typename ksP::targetP>> &res,
            const array<TRGSWFFT<typename ksP::domainP>, address_bit> &address,
            const TRLWE<typename ksP::domainP> &data,
            const KeySwitchingKey<ksP> &ksk)
{
    TRLWE<typename ksP::domainP> temp, acc;
    PolynomialMulByXaiMinusOnelvl1(temp[0], data[0],
                                   2 * ksP::domainP::n - (ksP::domainP::n >> 1));
    PolynomialMulByXaiMinusOnelvl1(temp[1], data[1],
                                   2 * ksP::domainP::n - (ksP::domainP::n >> 1));
    trgswfftExternalProduct<typename ksP::domainP>(temp, temp, address[width_bit - 1]);
    for (int i = 0; i < ksP::domainP::n; i++) {
        acc[0][i] = temp[0][i] + data[0][i];
        acc[1][i] = temp[1][i] + data[1][i];
    }

    for (uint32_t bit = 2; bit <= width_bit; bit++) {
        PolynomialMulByXaiMinusOnelvl1(
            temp[0], acc[0], 2 * ksP::domainP::n - (ksP::domainP::n >> bit));
        PolynomialMulByXaiMinusOnelvl1(
            temp[1], acc[1], 2 * ksP::domainP::n - (ksP::domainP::n >> bit));
        trgswfftExternalProduct<typename ksP::domainP>(temp, temp,
                                           address[width_bit - bit]);
        for (int i = 0; i < ksP::domainP::n; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }

    const uint32_t width = 1 << width_bit;
    array<TLWE<typename ksP::domainP>, width> reslvl1;
    for (int i = 0; i < width; i++)
        SampleExtractIndex<typename ksP::domainP>(reslvl1[i], acc, i);
    for (int i = 0; i < width; i++)
        IdentityKeySwitch<ksP>(res[i], reslvl1[i], ksk);
}

int main()
{
    constexpr uint32_t address_bit = 8;  // Address by words.
    constexpr uint32_t words_bit = 5;
    constexpr uint32_t width_bit =
        lvl1param::nbit -
        words_bit;  // log_2 of how many words are in one TRLWElvl1 message.
    static_assert(address_bit >= width_bit);
    constexpr uint32_t width = 1 << width_bit;
    constexpr uint32_t num_trlwe = 1 << (address_bit - width_bit);
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CloudKey *ck = new CloudKey(*sk);
    vector<array<uint8_t, lvl1param::n>> pmemory(num_trlwe);
    vector<array<uint32_t, lvl1param::n>> pmu(num_trlwe);
    vector<uint8_t> address(address_bit);
    vector<uint8_t> pres(width);

    for (array<uint8_t, lvl1param::n> &i : pmemory)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_trlwe; i++)
        for (int j = 0; j < lvl1param::n; j++)
            pmu[i][j] =
                pmemory[i][j]
                    ? 2 * lvl1param::μ
                    : -2 *
                          lvl1param::μ;  // This will increase noise torellance.
    for (uint8_t &p : address) p = binary(engine);

    array<TRGSWFFT<TFHEpp::lvl1param>, address_bit> bootedTGSW;
    vector<TLWE<lvl0param>> encaddress(address_bit);
    array<TRLWE<lvl1param>, num_trlwe> encmemory;
    vector<TLWE<lvl0param>> encres(width);

    encaddress = bootsSymEncrypt(address, *sk);
    for (int i = 0; i < num_trlwe; i++)
        encmemory[i] =
            trlweSymEncrypt<lvl1param>(pmu[i], lvl1param::α, (*sk).key.lvl1);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < words_bit; i++)
        CircuitBootstrappingFFT<lvl02param, lvl21param>(
            bootedTGSW[i], encaddress[i], (*ck).ck);
    for (int i = words_bit; i < address_bit; i++)
        CircuitBootstrappingFFTInv<lvl02param, lvl21param>(
            bootedTGSW[i], encaddress[i], (*ck).ck);
    TRLWE<lvl1param> encumemory;

    UROMUX<lvl1param,address_bit, width_bit>(encumemory, bootedTGSW, encmemory);
    LROMUX<lvl10param,address_bit, width_bit>(encres, bootedTGSW, encumemory,
                                   (*ck).gk.ksk);
    end = chrono::system_clock::now();

    pres = bootsSymDecrypt(encres, *sk);
    uint32_t uaddress = 0;
    uint32_t laddress = 0;
    for (int i = 0; i < (address_bit - width_bit); i++)
        uaddress += address[i + width_bit] << i;
    array<bool, lvl1param::n> umemory;
    umemory = trlweSymDecrypt<lvl1param>(encumemory, (*sk).key.lvl1);

    for (int i = 0; i < width_bit; i++)
        laddress += static_cast<uint32_t>(address[i]) << (i + words_bit);
    for (uint32_t i = 0; i < (1 << words_bit); i++)
        assert(static_cast<int>(pres[i]) ==
               static_cast<int>(pmemory[uaddress][laddress + i]));
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
}