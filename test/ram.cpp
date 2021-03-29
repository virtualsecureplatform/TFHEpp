#include <bitset>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

template <uint32_t address_bit>
void RAMUX(TRLWE<lvl1param> &res,
           const array<TRGSWFFT<lvl1param>, address_bit> &invaddress,
           const array<TRLWE<lvl1param>, 1 << address_bit> &data)
{
    constexpr uint32_t num_trlwe = 1 << address_bit;
    array<TRLWE<lvl1param>, num_trlwe / 2> temp;

    for (uint32_t index = 0; index < num_trlwe / 2; index++) {
        CMUXFFT<lvl1param>(temp[index], invaddress[0], data[2 * index],
                           data[2 * index + 1]);
    }

    for (uint32_t bit = 0; bit < (address_bit - 2); bit++) {
        const uint32_t stride = 1 << bit;
        for (uint32_t index = 0; index < (num_trlwe >> (bit + 2)); index++) {
            CMUXFFT<lvl1param>(temp[(2 * index) * stride], invaddress[bit + 1],
                               temp[(2 * index) * stride],
                               temp[(2 * index + 1) * stride]);
        }
    }
    constexpr uint32_t stride = 1 << (address_bit - 2);
    CMUXFFT<lvl1param>(res, invaddress[address_bit - 1], temp[0], temp[stride]);
}

int main()
{
    constexpr uint32_t address_bit = 9;  // Address by bytes.
    constexpr uint32_t memsize = 1 << address_bit;
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CloudKey *ck = new CloudKey(*sk);
    vector<uint8_t> pmemory(memsize);
    vector<array<uint32_t, lvl1param::n>> pmu(memsize);
    vector<uint8_t> address(address_bit);
    uint8_t pres;
    uint8_t wrflag;
    uint8_t writep;

    for (uint8_t &p : pmemory) p = binary(engine);

    for (int i = 0; i < memsize; i++) {
        pmu[i] = {};
        pmu[i][0] = pmemory[i] ? lvl1param::μ : -lvl1param::μ;
    }
    for (uint8_t &p : address) p = binary(engine);
    uint32_t addressint = 0;
    for (int i = 0; i < address_bit; i++)
        addressint += static_cast<uint32_t>(address[i]) << i;

    wrflag = binary(engine);
    writep = pmemory[addressint] > 0 ? 0 : 1;

    array<array<TRGSWFFT<lvl1param>, address_bit>, 2> *bootedTGSW =
        new array<array<TRGSWFFT<lvl1param>, address_bit>, 2>;
    vector<TLWE<lvl0param>> encaddress(address_bit);
    array<TRLWE<lvl1param>, memsize> encmemory;
    TLWE<lvl1param> encreadreslvl1;
    TLWE<lvl0param> encreadres;
    TRLWE<lvl1param> encumemory;
    TLWE<lvl0param> cs;
    TLWE<lvl0param> c1;
    TRLWE<lvl1param> writed;

    encaddress = bootsSymEncrypt(address, *sk);
    for (int i = 0; i < memsize; i++)
        encmemory[i] =
            trlweSymEncrypt<lvl1param>(pmu[i], lvl1param::α, (*sk).key.lvl1);
    cs = tlweSymEncrypt<lvl0param>(wrflag ? lvl0param::μ : -lvl0param::μ,
                                   lvl0param::α, (*sk).key.lvl0);
    c1 = tlweSymEncrypt<lvl0param>(writep ? lvl0param::μ : -lvl0param::μ,
                                   lvl0param::α, (*sk).key.lvl0);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    // Addres CB
    for (int i = 0; i < address_bit; i++) {
        CircuitBootstrappingFFTwithInv<lvl02param, lvl21param>(
            (*bootedTGSW)[1][i], (*bootedTGSW)[0][i], encaddress[i], (*ck).ck);
    }

    // Read
    RAMUX<address_bit>(encumemory, (*bootedTGSW)[0], encmemory);
    SampleExtractIndex<lvl1param>(encreadreslvl1, encumemory, 0);
    IdentityKeySwitch<lvl10param>(encreadres, encreadreslvl1, (*ck).gk.ksk);

    // Write
    HomMUXwoSE(writed, cs, c1, encreadres, (*ck).gk);
    for (int i = 0; i < memsize; i++) {
        const bitset<address_bit> addressbitset(i);
        TRLWE<lvl1param> temp = writed;
        for (int j = 0; j < address_bit; j++)
            CMUXFFT<lvl1param>(temp, (*bootedTGSW)[addressbitset[j]][j], temp,
                               encmemory[i]);
        TLWE<lvl1param> temp2;
        SampleExtractIndex<lvl1param>(temp2, temp, 0);
        TLWE<lvl0param> temp3;
        IdentityKeySwitch<lvl10param>(temp3, temp2, (*ck).gk.ksk);
        GateBootstrappingTLWE2TRLWEFFT<lvl01param>(encmemory[i], temp3,
                                                   (*ck).gk.bkfftlvl01);
    }

    end = chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
    pres = tlweSymDecrypt<lvl0param>(encreadres, (*sk).key.lvl0);

    assert(static_cast<int>(pres) == static_cast<int>(pmemory[addressint]));

    array<bool, lvl1param::n> pwriteres =
        trlweSymDecrypt<lvl1param>(encmemory[addressint], (*sk).key.lvl1);
    assert(static_cast<int>(pwriteres[0]) ==
           static_cast<int>((wrflag > 0) ? writep : pmemory[addressint]));

    cout << "Passed" << endl;
}