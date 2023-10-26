#include <cassert>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    constexpr uint32_t num_test = 1000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    cout << "test p=1" << endl;

    cout << "lvl1" << endl;
    for (int test = 0; test < num_test; test++) {
        lweKey key;

        array<bool, lvl1param::n> p;
        for (bool &i : p) i = (binary(engine) > 0);
        Polynomial<lvl1param> pmu;
        for (int i = 0; i < lvl1param::n; i++)
            pmu[i] = p[i] ? lvl1param::μ : -lvl1param::μ;
        TRLWE<lvl1param> c = trlweSymEncrypt<lvl1param>(pmu, key.lvl1);

        const Polynomial<TFHEpp::lvl1param> plainpoly = {
            static_cast<typename lvl1param::T>(1)};

        TRGSWFFT<lvl1param> trgswfft =
            trgswfftSymEncrypt<lvl1param>(plainpoly, key.lvl1);
        trgswfftExternalProduct<lvl1param>(c, c, trgswfft);
        array<bool, lvl1param::n> p2 = trlweSymDecrypt<lvl1param>(c, key.lvl1);
        for (int i = 0; i < lvl1param::n; i++) assert(p[i] == p2[i]);
    }
    cout << "Passed" << endl;

    cout << "lvl2" << endl;
    for (int test = 0; test < num_test; test++) {
        lweKey key;

        array<bool, lvl2param::n> p;
        for (bool &i : p) i = (binary(engine) > 0);
        Polynomial<lvl2param> pmu;
        for (int i = 0; i < lvl2param::n; i++)
            pmu[i] = p[i] ? lvl2param::μ : -lvl2param::μ;
        TRLWE<lvl2param> c = trlweSymEncrypt<lvl2param>(pmu, key.lvl2);

        const Polynomial<TFHEpp::lvl2param> plainpoly = {
            static_cast<typename lvl2param::T>(1)};

        TRGSWFFT<lvl2param> trgswfft =
            trgswfftSymEncrypt<lvl2param>(plainpoly, key.lvl2);
        trgswfftExternalProduct<lvl2param>(c, c, trgswfft);
        array<bool, lvl2param::n> p2 = trlweSymDecrypt<lvl2param>(c, key.lvl2);
        for (int i = 0; i < lvl2param::n; i++) assert(p[i] == p2[i]);
    }
    cout << "Passed" << endl;

    cout << "test p=-1" << endl;

    cout << "lvl1" << endl;
    for (int test = 0; test < num_test; test++) {
        lweKey key;

        array<bool, lvl1param::n> p;
        for (bool &i : p) i = binary(engine) > 0;
        array<typename TFHEpp::lvl1param::T, lvl1param::n> pmu;
        for (int i = 0; i < lvl1param::n; i++)
            pmu[i] = p[i] ? lvl1param::μ : -lvl1param::μ;
        TRLWE<lvl1param> c = trlweSymEncrypt<lvl1param>(pmu, key.lvl1);

        const Polynomial<TFHEpp::lvl1param> plainpoly = {
            static_cast<typename lvl1param::T>(-1)};

        TRGSWFFT<lvl1param> trgswfft =
            trgswfftSymEncrypt<lvl1param>(plainpoly, key.lvl1);
        trgswfftExternalProduct<lvl1param>(c, c, trgswfft);
        array<bool, lvl1param::n> p2 = trlweSymDecrypt<lvl1param>(c, key.lvl1);
        for (int i = 0; i < lvl1param::n; i++) assert(p[i] == !p2[i]);
    }
    cout << "Passed" << endl;

    cout << "lvl2" << endl;
    for (int test = 0; test < num_test; test++) {
        lweKey key;

        array<bool, lvl2param::n> p;
        for (bool &i : p) i = binary(engine) > 0;
        Polynomial<lvl2param> pmu;
        for (int i = 0; i < lvl2param::n; i++)
            pmu[i] = p[i] ? lvl2param::μ : -lvl2param::μ;
        TRLWE<lvl2param> c = trlweSymEncrypt<lvl2param>(pmu, key.lvl2);

        const Polynomial<TFHEpp::lvl2param> plainpoly = {
            static_cast<typename lvl2param::T>(-1)};

        TRGSWFFT<lvl2param> trgswfft =
            trgswfftSymEncrypt<lvl2param>(plainpoly, key.lvl2);
        trgswfftExternalProduct<lvl2param>(c, c, trgswfft);
        array<bool, lvl2param::n> p2 = trlweSymDecrypt<lvl2param>(c, key.lvl2);
        for (int i = 0; i < lvl2param::n; i++) assert(p[i] == !p2[i]);
    }
    cout << "Passed" << endl;
}