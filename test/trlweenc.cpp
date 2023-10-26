#include <array>
#include <cassert>
#include <iostream>
#include <tfhe++.hpp>

using namespace TFHEpp;

int main()
{
    cout << "lvl1" << endl;
    constexpr uint32_t num_test = 1000;
    for (int test = 0; test < num_test; test++) {
        random_device seed_gen;
        default_random_engine engine(seed_gen());
        uniform_int_distribution<typename lvl1param::T> binary(0, 1);

        lweKey key;
        array<bool, lvl1param::n> p;
        for (bool &i : p) i = binary(engine) > 0;
        array<typename lvl1param::T, lvl1param::n> pmu;
        for (int i = 0; i < lvl1param::n; i++)
            pmu[i] = p[i] ? lvl1param::μ : -lvl1param::μ;
        TRLWE<TFHEpp::lvl1param> c = trlweSymEncrypt<lvl1param>(pmu, key.lvl1);
        // if constexpr(hasq<lvl1param>) for (int i = 0; i < lvl1param::n; i++)
        // if(c[lvl1param::k][i] >= lvl1param::q)
        // std::cout<<i<<":"<<c[lvl1param::k][i]<<std::endl;
        array<bool, lvl1param::n> p2 = trlweSymDecrypt<lvl1param>(c, key.lvl1);
        for (int i = 0; i < lvl1param::n; i++) assert(p[i] == p2[i]);
    }
    cout << "Passed" << endl;

    cout << "lvl2" << endl;
    for (int test = 0; test < num_test; test++) {
        random_device seed_gen;
        default_random_engine engine(seed_gen());
        uniform_int_distribution<uint32_t> binary(0, 1);

        lweKey key;
        array<bool, lvl2param::n> p;
        for (bool &i : p) i = binary(engine) > 0;
        array<typename lvl2param::T, lvl2param::n> pmu;
        for (int i = 0; i < lvl2param::n; i++)
            pmu[i] = p[i] ? lvl2param::μ : -lvl2param::μ;
        TRLWE<lvl2param> c = trlweSymEncrypt<lvl2param>(pmu, key.lvl2);
        array<bool, lvl2param::n> p2 = trlweSymDecrypt<lvl2param>(c, key.lvl2);
        for (int i = 0; i < lvl2param::n; i++) assert(p[i] == p2[i]);
    }
    cout << "Passed" << endl;
}