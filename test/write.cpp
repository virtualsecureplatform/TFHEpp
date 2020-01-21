#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

template <uint32_t address_bit, uint32_t block, uint32_t bitpos>
void writeMUX(array<TRLWElvl1, (1 << address_bit)> &res,
              const array<array<TRGSWFFTlvl1, address_bit - 1>, 2> &address)
{
    trgswfftExternalProductlvl1(res[block + (1 << bitpos)], res[block],
                                address[0][bitpos]);
    trgswfftExternalProductlvl1(res[block], res[block], address[1][bitpos]);
    if constexpr (bitpos > 0) {
        writeMUX<address_bit, block + (1 << bitpos), bitpos - 1>(res, address);
        writeMUX<address_bit, block, bitpos - 1>(res, address);
    }
}

int main()
{
    const uint32_t address_bit = 9;  // Address by words.
    const uint32_t memsize = 1 << address_bit;
    random_device seeder;
    default_random_engine engine(seeder());
    uniform_int_distribution<uint8_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CloudKey *ck = new CloudKey(*sk);
    array<uint8_t, memsize> pmemory;
    array<array<uint32_t, DEF_N>, memsize> pmu = {};
    vector<uint8_t> address(address_bit);
    uint8_t pres;

    for (uint8_t &i : pmemory) i = binary(engine);
    for (int i = 0; i < memsize; i++) pmu[i][0] = pmemory[i] ? DEF_μ : -DEF_μ;
    for (uint8_t &p : address) p = binary(engine);

    array<array<TRGSWFFTlvl1, address_bit - 1>, 2> *bootedTGSW =
        new array<array<TRGSWFFTlvl1, address_bit - 1>, 2>;
    vector<TLWElvl0> encaddress(address_bit);
    array<TRLWElvl1, memsize> encmemory;
    TLWElvl0 encres;
    TRLWElvl1 datum;
    Polynomiallvl1 respoly = {};
    respoly[0] = pres ? DEF_μ : -DEF_μ;
    array<TRLWElvl1, memsize> trlweaddress;

    encaddress = bootsSymEncrypt(address, *sk);
    for (int i = 0; i < memsize; i++)
        encmemory[i] = trlweSymEncryptlvl1(pmu[i], DEF_αbk, (*sk).key.lvl1);
    datum = trlweSymEncryptlvl1(respoly, DEF_αbk, (*sk).key.lvl1);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < address_bit - 1; i++) {
        CircuitBootstrappingFFT((*bootedTGSW)[0][i], encaddress[i], (*ck).ck);
        for (int j = 0; j < 2 * DEF_l; j++)
            for (int k = 0; k < DEF_N; k++) {
                (*bootedTGSW)[1][i][j][0][k] = -(*bootedTGSW)[0][i][j][0][k];
                (*bootedTGSW)[1][i][j][1][k] = -(*bootedTGSW)[0][i][j][1][k];
            }
        for (int j = 0; j < DEF_l; j++)
            for (int k = 0; k < DEF_N / 2; k++) {
                const double doubleh =
                    static_cast<double>(1U << (32 - (j + 1) * DEF_Bgbit));
                (*bootedTGSW)[1][i][j][0][k] += doubleh;
                (*bootedTGSW)[1][i][j + DEF_l][1][k] += doubleh;
            }
    }
    TRLWElvl1 msbaddress;
    GateBootstrappingTLWE2TRLWEFFTlvl01(msbaddress, encaddress[address_bit - 1],
                                        (*ck).gk);

    trlweaddress[memsize >> 1] = msbaddress;
    trlweaddress[memsize >> 1][1][0] += DEF_μ;
    for (int i = 0; i < DEF_N; i++) {
        trlweaddress[0][0][i] = -msbaddress[0][i];
        trlweaddress[0][1][i] = -msbaddress[1][i];
    }
    trlweaddress[0][1][0] += DEF_μ;

    writeMUX<address_bit, (memsize >> 1), address_bit - 2>(trlweaddress,
                                                           *bootedTGSW);
    writeMUX<address_bit, 0, address_bit - 2>(trlweaddress, *bootedTGSW);

    for (int i = 0; i < memsize; i++) {
        trlweaddress[i][1][0] -= DEF_μ;
        ExtractSwitchAndHomMUX(encmemory[i], trlweaddress[i], datum,
                               encmemory[i], (*ck).gk);
    }

    end = chrono::system_clock::now();
    uint32_t intaddress = 0;
    for (int i = 0; i < address_bit; i++) intaddress += address[i] << i;
    assert(pres ==
           trlweSymDecryptlvl1(encmemory[intaddress], (*sk).key.lvl1)[0]);
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed << "ms" << endl;
}