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
    CircuitKey<lvl02param, lvl22param> *ck =
        new CircuitKey<lvl02param, lvl22param>(*sk);
    vector<array<uint8_t, TFHEpp::lvl2param::n>> pa(num_test);
    vector<array<typename lvl2param::T, TFHEpp::lvl2param::n>> pmu(num_test);
    vector<uint8_t> pones(num_test);
    array<bool, TFHEpp::lvl2param::n> pres;
    for (array<uint8_t, TFHEpp::lvl2param::n> &i : pa)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < TFHEpp::lvl2param::n; j++)
            pmu[i][j] = pa[i][j] ? lvl2param::μ : -lvl2param::μ;
    for (int i = 0; i < num_test; i++) pones[i] = true;
    vector<TRLWE<lvl2param>> ca(num_test);
    vector<TLWE<lvl0param>> cones(num_test);
    vector<TRGSWFFT<lvl2param>> bootedTGSW(num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = trlweSymEncryptlvl2(pmu[i], TFHEpp::lvl2param::α, sk->key.lvl2);
    cones = bootsSymEncrypt(pones, *sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CircuitBootstrappingFFTlvl02(bootedTGSW[test], cones[test], *ck);
    }
    end = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        trgswfftExternalProductlvl2(ca[test], ca[test], bootedTGSW[test]);
        pres = trlweSymDecryptlvl2(ca[test], sk->key.lvl2);
        for (int i = 0; i < TFHEpp::lvl2param::n; i++)
            assert(pres[i] == pa[test][i]);
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}
