#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

template <uint32_t address_bit, uint32_t block, uint32_t bitpos>
void writeMUX(
    array<TRLWE<lvl1param>, (1 << address_bit)> &res,
    const array<array<TRGSWFFT<lvl1param>, address_bit - 1>, 2> &address)
{
    trgswfftExternalProduct<lvl1param>(res[block + (1 << bitpos)], res[block],
                                       address[0][bitpos]);
    trgswfftExternalProduct<lvl1param>(res[block], res[block],
                                       address[1][bitpos]);
    if constexpr (bitpos > 0) {
        writeMUX<address_bit, block + (1 << bitpos), bitpos - 1>(res, address);
        writeMUX<address_bit, block, bitpos - 1>(res, address);
    }
}

int main()
{
    constexpr uint32_t address_bit = 9;  // Address by words.
    constexpr uint32_t memsize = 1 << address_bit;
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CloudKey *ck = new CloudKey(*sk);
    array<uint8_t, memsize> pmemory;
    array<array<typename lvl1param::T, lvl1param::n>, memsize> pmu = {};
    vector<uint8_t> address(address_bit);
    uint8_t pres;

    for (uint8_t &i : pmemory) i = binary(engine);
    for (int i = 0; i < memsize; i++)
        pmu[i][0] = pmemory[i] ? lvl0param::μ : -lvl0param::μ;
    for (uint8_t &p : address) p = binary(engine);

    array<array<TRGSWFFT<lvl1param>, address_bit - 1>, 2> *bootedTGSW =
        new array<array<TRGSWFFT<lvl1param>, address_bit - 1>, 2>;
    vector<TLWE<lvl0param>> encaddress(address_bit);
    array<TRLWE<lvl1param>, memsize> encmemory;
    TLWE<lvl0param> encres;
    TRLWE<lvl1param> datum;
    Polynomial<lvl1param> respoly = {};
    respoly[0] = pres ? lvl1param::μ : -lvl1param::μ;
    array<TRLWE<lvl1param>, memsize> trlweaddress;

    encaddress = bootsSymEncrypt(address, *sk);
    for (int i = 0; i < memsize; i++)
        encmemory[i] =
            trlweSymEncrypt<lvl1param>(pmu[i], lvl1param::α, (*sk).key.lvl1);
    datum = trlweSymEncrypt<lvl1param>(respoly, lvl1param::α, (*sk).key.lvl1);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < address_bit - 1; i++) {
        CircuitBootstrappingFFT<lvl02param, lvl21param>(
            (*bootedTGSW)[0][i], encaddress[i], (*ck).ck);
        for (int j = 0; j < 2 * lvl1param::l; j++)
            for (int k = 0; k < lvl1param::n; k++) {
                (*bootedTGSW)[1][i][j][0][k] = -(*bootedTGSW)[0][i][j][0][k];
                (*bootedTGSW)[1][i][j][1][k] = -(*bootedTGSW)[0][i][j][1][k];
            }
        for (int j = 0; j < lvl1param::l; j++)
            for (int k = 0; k < lvl1param::n / 2; k++) {
                const double doubleh = static_cast<double>(
                    1U << (32 - (j + 1) * lvl1param::Bgbit));
                (*bootedTGSW)[1][i][j][0][k] += doubleh;
                (*bootedTGSW)[1][i][j + lvl1param::l][1][k] += doubleh;
            }
    }
    TRLWE<lvl1param> msbaddress;
    GateBootstrappingTLWE2TRLWEFFT<lvl01param>(
        msbaddress, encaddress[address_bit - 1], (*ck).gk.bkfftlvl01);

    trlweaddress[memsize >> 1] = msbaddress;
    trlweaddress[memsize >> 1][1][0] += lvl1param::μ;
    for (int i = 0; i < lvl1param::n; i++) {
        trlweaddress[0][0][i] = -msbaddress[0][i];
        trlweaddress[0][1][i] = -msbaddress[1][i];
    }
    trlweaddress[0][1][0] += lvl1param::μ;

    writeMUX<address_bit, (memsize >> 1), address_bit - 2>(trlweaddress,
                                                           *bootedTGSW);
    writeMUX<address_bit, 0, address_bit - 2>(trlweaddress, *bootedTGSW);

    for (int i = 0; i < memsize; i++) {
        trlweaddress[i][1][0] -= lvl1param::μ;
        ExtractSwitchAndHomMUX(encmemory[i], trlweaddress[i], datum,
                               encmemory[i], (*ck).gk);
    }

    end = chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
    uint32_t intaddress = 0;
    for (int i = 0; i < address_bit; i++) intaddress += address[i] << i;
    assert(pres == trlweSymDecrypt<lvl1param>(encmemory[intaddress],
                                              (*sk).key.lvl1)[0]);
    cout << "Passed" << endl;
}