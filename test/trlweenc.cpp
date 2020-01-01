#include <array>
#include <cassert>
#include <iostream>
#include <tfhe++.hpp>

using namespace TFHEpp;

int main()
{
    cout << "lvl1" << endl;
    const uint32_t num_test = 1000;
    for (int test = 0; test < num_test; test++) {
        random_device seed_gen;
        default_random_engine engine(seed_gen());
        uniform_int_distribution<uint32_t> binary(0, 1);

        lweKey key;
        array<bool, DEF_N> p;
        for (bool &i : p) i = binary(engine) > 0;
        array<uint32_t, DEF_N> pmu;
        for (int i = 0; i < DEF_N; i++) pmu[i] = p[i] ? DEF_μ : -DEF_μ;
        TRLWElvl1 c = trlweSymEncryptlvl1(pmu, DEF_αbk, key.lvl1);
        array<bool, DEF_N> p2 = trlweSymDecryptlvl1(c, key.lvl1);
        for (int i = 0; i < DEF_N; i++) assert(p[i] == p2[i]);
    }
    cout << "Passed" << endl;

    cout << "lvl2" << endl;
    for (int test = 0; test < num_test; test++) {
        random_device seed_gen;
        default_random_engine engine(seed_gen());
        uniform_int_distribution<uint32_t> binary(0, 1);

        lweKey key;
        array<bool, DEF_nbar> p;
        for (bool &i : p) i = binary(engine) > 0;
        array<uint64_t, DEF_nbar> pmu;
        for (int i = 0; i < DEF_nbar; i++) pmu[i] = p[i] ? DEF_μbar : -DEF_μbar;
        TRLWElvl2 c = trlweSymEncryptlvl2(pmu, DEF_αbk, key.lvl2);
        array<bool, DEF_nbar> p2 = trlweSymDecryptlvl2(c, key.lvl2);
        for (int i = 0; i < DEF_nbar; i++) assert(p[i] == p2[i]);
    }
    cout << "Passed" << endl;
}