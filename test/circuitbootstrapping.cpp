#include <gperftools/profiler.h>

#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    constexpr uint32_t num_test = 10;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CircuitKey<lvl02param, lvl21param> *ck =
        new CircuitKey<lvl02param, lvl21param>(*sk);
    vector<array<uint8_t, lvl1param::n>> pa(num_test);
    vector<array<lvl1param::T, lvl1param::n>> pmu(num_test);
    vector<uint8_t> pones(num_test);
    array<bool, lvl1param::n> pres;
    for (array<uint8_t, lvl1param::n> &i : pa)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < lvl1param::n; j++)
            pmu[i][j] = pa[i][j] ? lvl1param::μ : -lvl1param::μ;
    for (int i = 0; i < num_test; i++) pones[i] = true;
    vector<TRLWE<lvl1param>> ca(num_test);
    vector<TLWE<lvl0param>> cones(num_test);
    vector<TRGSWFFT<lvl1param>> bootedTGSW(num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = trlweSymEncrypt<lvl1param>(pmu[i], lvl1param::α, sk->key.lvl1);
    cones = bootsSymEncrypt(pones, *sk);

    chrono::system_clock::time_point start, end;
    ProfilerStart("cb.prof");
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CircuitBootstrappingFFT<lvl02param, lvl21param>(bootedTGSW[test],
                                                        cones[test], *ck);
    }
    end = chrono::system_clock::now();
    ProfilerStop();
    for (int test = 0; test < num_test; test++) {
        trgswfftExternalProduct<lvl1param>(ca[test], ca[test],
                                           bootedTGSW[test]);
        pres = trlweSymDecrypt<lvl1param>(ca[test], sk->key.lvl1);
        for (int i = 0; i < lvl1param::n; i++) assert(pres[i] == pa[test][i]);
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}
