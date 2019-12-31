#include <cassert>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    const uint32_t num_test = 1000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    cout << "test p=1" << endl;
    for (int test = 0; test < num_test; test++) {
        lweKey key;

        array<bool, DEF_N> p;
        for (bool &i : p) i = (binary(engine) > 0);
        array<uint32_t, DEF_N> pmu;
        for (int i = 0; i < DEF_N; i++) pmu[i] = p[i] ? DEF_MU : -DEF_MU;
        TRLWElvl1 c = trlweSymEncryptlvl1(pmu, DEF_αbk, key.lvl1);

        TRGSWFFTlvl1 trgswfft = trgswfftSymEncryptlvl1(1, DEF_αbk, key.lvl1);
        trgswfftExternalProductlvl1(c, trgswfft);
        array<bool, DEF_N> p2 = trlweSymDecryptlvl1(c, key.lvl1);
        for (int i = 0; i < DEF_N; i++) assert(p[i] == p2[i]);
    }
    cout << "Passed" << endl;

    cout << "test p=-1" << endl;
    for (int test = 0; test < num_test; test++) {
        lweKey key;

        array<bool, DEF_N> p;
        for (bool &i : p) i = binary(engine) > 0;
        array<uint32_t, DEF_N> pmu;
        for (int i = 0; i < DEF_N; i++) pmu[i] = p[i] ? DEF_MU : -DEF_MU;
        TRLWElvl1 c = trlweSymEncryptlvl1(pmu, DEF_αbk, key.lvl1);

        TRGSWFFTlvl1 trgswfft = trgswfftSymEncryptlvl1(-1, DEF_αbk, key.lvl1);
        trgswfftExternalProductlvl1(c, trgswfft);
        array<bool, DEF_N> p2 = trlweSymDecryptlvl1(c, key.lvl1);
        for (int i = 0; i < DEF_N; i++) assert(p[i] == !p2[i]);
    }
    cout << "Passed" << endl;
}