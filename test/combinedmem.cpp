#include <bitset>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

template <uint32_t address_bit, uint32_t words_bit>
void combUROMUX(TRLWE<lvl1param> &res,
                const array<TRGSWFFT<lvl1param>, address_bit> &invaddress,
                const array<TRLWE<lvl1param>, 1 << (address_bit + words_bit -
                                                    lvl1param::nbit)> &data)
{
    constexpr uint32_t width_bit =
        lvl1param::nbit -
        words_bit;  // log_2 of how many words are in one TRLWElvl1 message.
    constexpr uint32_t Ubit = address_bit - width_bit;
    constexpr uint32_t num_trlwe = 1 << (Ubit);
    array<TRLWE<lvl1param>, num_trlwe / 2> temp;

    for (uint32_t index = 0; index < num_trlwe / 2; index++) {
        CMUXFFT<lvl1param>(temp[index], invaddress[width_bit], data[2 * index],
                    data[2 * index + 1]);
    }

    for (uint32_t bit = 0; bit < (Ubit - 2); bit++) {
        const uint32_t stride = 1 << bit;
        for (uint32_t index = 0; index < (num_trlwe >> (bit + 2)); index++) {
            CMUXFFT<lvl1param>(
                temp[(2 * index) * stride], invaddress[width_bit + bit + 1],
                temp[(2 * index) * stride], temp[(2 * index + 1) * stride]);
        }
    }

    constexpr uint32_t stride = 1 << (Ubit - 2);
    CMUXFFT<lvl1param>(res, invaddress[address_bit - 1], temp[0], temp[stride]);
}

template <uint32_t address_bit, uint32_t words_bit>
void combLROMUX(array<TLWE<lvl0param>, 1U << words_bit> &res,
                const array<TRGSWFFT<lvl1param>, address_bit> &address,
                const TRLWE<lvl1param> &data,
                const KeySwitchingKey<lvl10param> &ksk)
{
    constexpr uint32_t width_bit =
        lvl1param::nbit -
        words_bit;  // log_2 of how many words are in one TRLWElvl1 message.
    TRLWE<lvl1param> temp, acc;
    PolynomialMulByXaiMinusOnelvl1(temp[0], data[0],
                                   2 * lvl1param::n - (lvl1param::n >> 1));
    PolynomialMulByXaiMinusOnelvl1(temp[1], data[1],
                                   2 * lvl1param::n - (lvl1param::n >> 1));
    trgswfftExternalProduct<lvl1param>(temp, temp, address[width_bit - 1]);
    for (int i = 0; i < lvl1param::n; i++) {
        acc[0][i] = temp[0][i] + data[0][i];
        acc[1][i] = temp[1][i] + data[1][i];
    }

    for (uint32_t bit = 2; bit <= width_bit; bit++) {
        PolynomialMulByXaiMinusOnelvl1(
            temp[0], acc[0], 2 * lvl1param::n - (lvl1param::n >> bit));
        PolynomialMulByXaiMinusOnelvl1(
            temp[1], acc[1], 2 * lvl1param::n - (lvl1param::n >> bit));
        trgswfftExternalProduct<lvl1param>(temp, temp, address[width_bit - bit]);
        for (int i = 0; i < lvl1param::n; i++) {
            acc[0][i] += temp[0][i];
            acc[1][i] += temp[1][i];
        }
    }

    const uint32_t width = 1 << width_bit;
    array<TLWE<lvl1param>, width> reslvl1;
    for (int i = 0; i < width; i++)
        SampleExtractIndex<lvl1param>(reslvl1[i], acc, i);
    for (int i = 0; i < width; i++)
        IdentityKeySwitch<lvl10param>(res[i], reslvl1[i], ksk);
}

template <uint32_t address_bit, uint32_t words_bit>
void combRAMUX(array<TRLWE<lvl1param>, 1U << words_bit> &res,
               const array<TRGSWFFT<lvl1param>, address_bit> &invaddress,
               const array<array<TRLWE<lvl1param>, 1 << address_bit>,
                           1U << words_bit> &data)
{
    constexpr uint32_t words = 1U << words_bit;
    constexpr uint32_t num_trlwe = 1 << address_bit;
    for (int i = 0; i < words; i++) {
        array<TRLWE<lvl1param>, num_trlwe / 2> temp;

        for (uint32_t index = 0; index < num_trlwe / 2; index++) {
            CMUXFFT<lvl1param>(temp[index], invaddress[0], data[i][2 * index],
                        data[i][2 * index + 1]);
        }

        for (uint32_t bit = 0; bit < (address_bit - 2); bit++) {
            const uint32_t stride = 1 << bit;
            for (uint32_t index = 0; index < (num_trlwe >> (bit + 2));
                 index++) {
                CMUXFFT<lvl1param>(temp[(2 * index) * stride], invaddress[bit + 1],
                            temp[(2 * index) * stride],
                            temp[(2 * index + 1) * stride]);
            }
        }
        constexpr uint32_t stride = 1 << (address_bit - 2);
        CMUXFFT<lvl1param>(res[i], invaddress[address_bit - 1], temp[0], temp[stride]);
    }
}

template <uint32_t address_bit, uint32_t words_bit>
void combWRAM(
    array<array<TRLWE<lvl1param>, 1U << address_bit>, 1U << words_bit> &encram,
    const array<array<TRGSWFFT<lvl1param>, address_bit>, 2> &address,
    const array<TRLWE<lvl1param>, 1U << words_bit> &encwritep,
    const GateKey &gk)
{
    constexpr uint32_t memsize = 1U << address_bit;
    constexpr uint32_t words = 1U << words_bit;

#pragma omp parallel for
    for (int i = 0; i < memsize; i++) {
        const bitset<address_bit> addressbitset(i);
        for (int j = 0; j < words; j++) {
            TRLWE<lvl1param> temp = encwritep[j];
            for (int k = 0; k < address_bit; k++)
                CMUXFFT<lvl1param>(temp, address[addressbitset[k]][k], temp,
                            encram[j][i]);
            TLWE<lvl1param> temp2;
            SampleExtractIndex<lvl1param>(temp2, temp, 0);
            TLWE<lvl0param> temp3;
            IdentityKeySwitch<lvl10param>(temp3, temp2, gk.ksk);
            GateBootstrappingTLWE2TRLWEFFT<lvl01param>(encram[j][i], temp3,
                                                gk.bkfftlvl01);
        }
    }
}

int main()
{
    constexpr uint32_t address_bit = 10;  // Address by bytes.
    constexpr uint32_t words_bit = 3;
    constexpr uint32_t words = 1U << words_bit;
    constexpr uint32_t memsize = 1 << address_bit;
    constexpr uint32_t numromtrlwe =
        1U << (address_bit - 1 + words_bit - lvl1param::nbit);
    constexpr uint32_t numramtrlwe = 1U << (address_bit - 1);
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CloudKey *ck = new CloudKey(*sk);
    vector<uint8_t> ramp(memsize / 2 * words);  // unit of memsize is byte(8bit)
    vector<uint8_t> romp(memsize / 2 * words);
    vector<array<array<uint32_t, lvl1param::n>, numramtrlwe>> ramu(words);
    vector<array<uint32_t, lvl1param::n>> romu(numromtrlwe);
    vector<uint8_t> address(address_bit);
    array<uint8_t, words> pres;
    array<uint8_t, words> writep;

    for (uint8_t wrflag = 0; wrflag <= 1; wrflag++) {
        for (uint8_t msbaddress = 0; msbaddress <= 1; msbaddress++) {
            for (uint8_t &p : address) p = binary(engine);
            address[address_bit - 1] = msbaddress;
            for (uint8_t &p : ramp) p = binary(engine);
            for (uint8_t &p : romp) p = binary(engine);

            for (int i = 0; i < words; i++) {
                for (int j = 0; j < numramtrlwe; j++) {
                    ramu[i][j] = {};
                    ramu[i][j][0] =
                        ramp[j * words + i] ? lvl1param::μ : -lvl1param::μ;
                }
            }
            for (int i = 0; i < numromtrlwe; i++) {
                for (int j = 0; j < lvl1param::n; j++) {
                    romu[i][j] = romp[i * lvl1param::n + j] ? lvl1param::μ
                                                            : -lvl1param::μ;
                }
            }

            uint32_t addressint = 0;
            for (int i = 0; i < address_bit - 1; i++)
                addressint += static_cast<uint32_t>(address[i]) << i;

            for (int i = 0; i < words; i++)
                writep[i] = ramp[addressint * words + i] > 0 ? 0 : 1;

            array<array<TRGSWFFT<lvl1param>, address_bit - 1>, 2> *bootedTGSW =
                new array<array<TRGSWFFT<lvl1param>, address_bit - 1>,
                          2>;  // MSB of address is evaluated by HomMUX, not
                               // CMUX.
            vector<TLWE<lvl0param>> encaddress(address_bit);
            array<TRLWE<lvl1param>, numromtrlwe> encrom;
            array<array<TRLWE<lvl1param>, numramtrlwe>, words> encram;

            array<TRLWE<lvl1param>, words> encramread;
            array<TLWE<lvl1param>, words> encramreadlvl1;
            array<TLWE<lvl0param>, words> encramreadres;

            TRLWE<lvl1param> encumemory;
            array<TLWE<lvl0param>, words> encromreadres;

            array<TLWE<lvl0param>, words> encreadres;

            TLWE<lvl0param> encwrflag;
            array<TLWE<lvl0param>, words> encwritep;
            array<TRLWE<lvl1param>, words> writed;

            encaddress = bootsSymEncrypt(address, *sk);
            for (int i = 0; i < numromtrlwe; i++)
                encrom[i] =
                    trlweSymEncryptlvl1(romu[i], lvl1param::α, (*sk).key.lvl1);
            for (int i = 0; i < words; i++)
                for (int j = 0; j < numramtrlwe; j++)
                    encram[i][j] = trlweSymEncryptlvl1(ramu[i][j], lvl1param::α,
                                                       (*sk).key.lvl1);

            encwrflag = tlweSymEncrypt<lvl0param>(
                (wrflag > 0) ? lvl0param::μ : -lvl0param::μ, lvl0param::α,
                (*sk).key.lvl0);
            for (int i = 0; i < words; i++)
                encwritep[i] = tlweSymEncrypt<lvl0param>(
                    writep[i] ? lvl0param::μ : -lvl0param::μ, lvl0param::α,
                    (*sk).key.lvl0);

            chrono::system_clock::time_point start, end;
            start = chrono::system_clock::now();
            // Addres CB
            for (int i = 0; i < address_bit - 1; i++) {
                CircuitBootstrappingFFTwithInvlvl01((*bootedTGSW)[1][i],
                                                    (*bootedTGSW)[0][i],
                                                    encaddress[i], (*ck).ck);
            }

            // Read
            combRAMUX<address_bit - 1, words_bit>(encramread, (*bootedTGSW)[0],
                                                  encram);
            for (int i = 0; i < words; i++) {
                SampleExtractIndex<lvl1param>(encramreadlvl1[i], encramread[i],
                                              0);
                IdentityKeySwitch<lvl10param>(encramreadres[i],
                                              encramreadlvl1[i], (*ck).gk.ksk);
            }

            combUROMUX<address_bit - 1, words_bit>(encumemory, (*bootedTGSW)[0],
                                                   encrom);
            combLROMUX<address_bit - 1, words_bit>(
                encromreadres, (*bootedTGSW)[1], encumemory, (*ck).gk.ksk);

            for (int i = 0; i < words; i++)
                HomMUX(encreadres[i], encaddress[address_bit - 1],
                       encramreadres[i], encromreadres[i], (*ck).gk);

            // Controll
            TLWE<lvl0param> cs;
            HomAND(cs, encwrflag, encaddress[address_bit - 1], (*ck).gk);
            for (int i = 0; i < words; i++)
                HomMUXwoSE(writed[i], cs, encwritep[i], encramreadres[i],
                           (*ck).gk);

            // Write
            combWRAM<address_bit - 1, words_bit>(encram, *bootedTGSW, writed,
                                                 (*ck).gk);

            end = chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      start)
                    .count();
            cout << elapsed << "ms" << endl;

            // test
            for (int i = 0; i < words; i++)
                pres[i] =
                    tlweSymDecrypt<lvl0param>(encreadres[i], (*sk).key.lvl0);

            for (int i = 0; i < words; i++)
                assert(static_cast<int>(ramp[addressint * words + i]) ==
                       static_cast<int>(tlweSymDecrypt<lvl0param>(
                           encramreadres[i], (*sk).key.lvl0)));
            for (int i = 0; i < words; i++)
                assert(static_cast<int>(romp[addressint * words + i]) ==
                       static_cast<int>(tlweSymDecrypt<lvl0param>(
                           encromreadres[i], (*sk).key.lvl0)));

            for (int i = 0; i < words; i++)
                assert(static_cast<int>(pres[i]) ==
                       static_cast<int>(address[address_bit - 1] > 0
                                            ? ramp[addressint * words + i]
                                            : romp[addressint * words + i]));

            array<array<bool, lvl1param::n>, words> pwriteres;
            for (int i = 0; i < words; i++)
                pwriteres[i] =
                    trlweSymDecryptlvl1(encram[i][addressint], (*sk).key.lvl1);

            cout << static_cast<int>(wrflag > 0) << ":"
                 << static_cast<int>(
                        tlweSymDecrypt<lvl0param>(encwrflag, (*sk).key.lvl0))
                 << endl;
            cout << static_cast<int>(address[address_bit - 1] > 0) << ":"
                 << static_cast<int>(tlweSymDecrypt<lvl0param>(
                        encaddress[address_bit - 1], (*sk).key.lvl0))
                 << endl;
            bool csp = ((wrflag > 0) & (address[address_bit - 1] > 0));
            cout << static_cast<int>(csp) << ":"
                 << static_cast<int>(
                        tlweSymDecrypt<lvl0param>(cs, (*sk).key.lvl0))
                 << endl;
            assert(static_cast<int>(tlweSymDecrypt<lvl0param>(
                       cs, (*sk).key.lvl0)) == static_cast<int>(csp));
            array<array<bool, lvl1param::n>, words> writedp;
            for (int i = 0; i < words; i++)
                writedp[i] = trlweSymDecryptlvl1(writed[i], (*sk).key.lvl1);

            for (int i = 0; i < words; i++)
                assert(static_cast<int>(writep[i]) ==
                       static_cast<int>(tlweSymDecrypt<lvl0param>(
                           encwritep[i], (*sk).key.lvl0)));
            for (int i = 0; i < words; i++)
                assert(static_cast<int>(writedp[i][0]) ==
                       static_cast<int>(csp ? writep[i]
                                            : ramp[addressint * words + i]));
            for (int i = 0; i < words; i++)
                assert(static_cast<int>(pwriteres[i][0]) ==
                       static_cast<int>(csp ? writep[i]
                                            : ramp[addressint * words + i]));
        }
    }
    cout << "Passed" << endl;
}