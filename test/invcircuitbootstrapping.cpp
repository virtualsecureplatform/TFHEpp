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

    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl21param;

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);
    vector<array<uint8_t, lvl1param::n>> pa(num_test);
    vector<array<typename TFHEpp::lvl1param::T, lvl1param::n>> pmu(num_test);
    vector<uint8_t> pzeros(num_test);
    array<bool, lvl1param::n> pres;
    for (array<uint8_t, lvl1param::n> &i : pa)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < lvl1param::n; j++)
            pmu[i][j] = pa[i][j] ? lvl1param::μ : -lvl1param::μ;
    for (int i = 0; i < num_test; i++) pzeros[i] = false;
    vector<TRLWE<lvl1param>,
           TFHEpp::AlignedAllocator<TFHEpp::TRLWE<TFHEpp::lvl1param>, 64>>
        ca(num_test);
    vector<TLWE<lvl1param>> czeros(num_test);
    vector<TRGSWFFT<lvl1param>,
           TFHEpp::AlignedAllocator<TFHEpp::TRGSWFFT<TFHEpp::lvl1param>, 64>>
        bootedTGSW(num_test);
    vector<TRGSWFFT<lvl1param>,
           TFHEpp::AlignedAllocator<TFHEpp::TRGSWFFT<TFHEpp::lvl1param>, 64>>
        invbootedTGSW(num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = trlweSymEncrypt<lvl1param>(pmu[i], sk->key.lvl1);
    czeros = bootsSymEncrypt(pzeros, *sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CircuitBootstrappingFFTwithInv<lvl10param, lvl02param, lvl21param>(
            bootedTGSW[test], invbootedTGSW[test], czeros[test], ek);
    }
    end = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        trgswfftExternalProduct<lvl1param>(ca[test], ca[test],
                                           invbootedTGSW[test]);
        pres = trlweSymDecrypt<lvl1param>(ca[test], sk->key.lvl1);
        for (int i = 0; i < lvl1param::n; i++) assert(pres[i] == pa[test][i]);
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}