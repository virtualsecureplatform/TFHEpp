#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

template <uint32_t address_bit, uint32_t width_bit>
void UROMUX(TRLWElvl1 &res, const array<TRGSWFFTlvl1, address_bit> &invaddress,
            const array<TRLWElvl1, 1 << (address_bit - width_bit)> &data)
{
    constexpr uint32_t Ubit = address_bit - width_bit;
    constexpr uint32_t num_trlwe = 1 << (Ubit);
    array<TRLWElvl1, num_trlwe / 2> temp;

    for (uint32_t index = 0; index < num_trlwe / 2; index++) {
        CMUXFFTlvl1(temp[index], invaddress[width_bit], data[2 * index],
                    data[2 * index + 1]);
    }

    for (uint32_t bit = 0; bit < (Ubit - 2); bit++) {
        const uint32_t stride = 1 << bit;
        for (uint32_t index = 0; index < (num_trlwe >> (bit + 2)); index++) {
            CMUXFFTlvl1(
                temp[(2 * index) * stride], invaddress[width_bit + bit + 1],
                temp[(2 * index) * stride], temp[(2 * index + 1) * stride]);
        }
    }

    constexpr uint32_t stride = 1 << (Ubit - 2);
    CMUXFFTlvl1(res, invaddress[address_bit - 1], temp[0], temp[stride]);
}

template <uint32_t address_bit, uint32_t width_bit>
void LROMUX(vector<TLWElvl0> &res,
            const array<TRGSWFFTlvl1, address_bit> &address,
            const TRLWElvl1 &data, const KeySwitchingKey &ksk)
{
    TRLWElvl1 temp, acc;
    PolynomialMulByXaiMinusOnelvl1(temp[0], data[0], 2 * DEF_N - (DEF_N >> 1));
    PolynomialMulByXaiMinusOnelvl1(temp[1], data[1], 2 * DEF_N - (DEF_N >> 1));
    trgswfftExternalProductlvl1(temp, temp, address[width_bit - 1]);
    for (int i = 0; i < DEF_N; i++) {
        acc[0][i] = temp[0][i] + data[0][i];
        acc[1][i] = temp[1][i] + data[1][i];
    }

    for (uint32_t bit = 2; bit <= width_bit; bit++) {
        PolynomialMulByXaiMinusOnelvl1(temp[0], acc[0],
                                       2 * DEF_N - (DEF_N >> bit));
        PolynomialMulByXaiMinusOnelvl1(temp[1], acc[1],
                                       2 * DEF_N - (DEF_N >> bit));
        trgswfftExternalProductlvl1(temp, temp, address[width_bit - bit]);
        for (int i = 0; i < DEF_N; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }

    const uint32_t width = 1 << width_bit;
    array<TLWElvl1, width> reslvl1;
    for (int i = 0; i < width; i++) SampleExtractIndexlvl1(reslvl1[i], acc, i);
    for (int i = 0; i < width; i++)
        IdentityKeySwitchlvl10(res[i], reslvl1[i], ksk);
}

int main()
{
    constexpr uint32_t address_bit = 8;  // Address by words.
    constexpr uint32_t words_bit = 5;
    constexpr uint32_t width_bit =
        DEF_Nbit -
        words_bit;  // log_2 of how many words are in one TRLWElvl1 message.
    static_assert(address_bit >= width_bit);
    const uint32_t width = 1 << width_bit;
    const uint32_t num_trlwe = 1 << (address_bit - width_bit);
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CloudKey *ck = new CloudKey(*sk);
    vector<array<uint8_t, DEF_N>> pmemory(num_trlwe);
    vector<array<uint32_t, DEF_N>> pmu(num_trlwe);
    vector<uint8_t> address(address_bit);
    vector<uint8_t> pres(width);

    for (array<uint8_t, DEF_N> &i : pmemory)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_trlwe; i++)
        for (int j = 0; j < DEF_N; j++)
            pmu[i][j] = pmemory[i][j] ? 2*DEF_μ : -2*DEF_μ; // This will increase noise torellance.
    for (uint8_t &p : address) p = binary(engine);

    array<TRGSWFFTlvl1, address_bit> bootedTGSW;
    vector<TLWElvl0> encaddress(address_bit);
    array<TRLWElvl1, num_trlwe> encmemory;
    vector<TLWElvl0> encres(width);

    encaddress = bootsSymEncrypt(address, *sk);
    for (int i = 0; i < num_trlwe; i++)
        encmemory[i] = trlweSymEncryptlvl1(pmu[i], DEF_αbk, (*sk).key.lvl1);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < words_bit; i++)
        CircuitBootstrappingFFT(bootedTGSW[i], encaddress[i], (*ck).ck);
    for (int i = words_bit; i < address_bit; i++)
        CircuitBootstrappingFFTInv(bootedTGSW[i], encaddress[i], (*ck).ck);
    TRLWElvl1 encumemory;

    UROMUX<address_bit, width_bit>(encumemory, bootedTGSW, encmemory);
    LROMUX<address_bit, width_bit>(encres, bootedTGSW, encumemory,
                                   (*ck).gk.ksk);
    end = chrono::system_clock::now();

    pres = bootsSymDecrypt(encres, *sk);
    uint32_t uaddress = 0;
    uint32_t laddress = 0;
    for (int i = 0; i < (address_bit - width_bit); i++)
        uaddress += address[i + width_bit] << i;
    array<bool, DEF_N> umemory;
    umemory = trlweSymDecryptlvl1(encumemory, (*sk).key.lvl1);

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