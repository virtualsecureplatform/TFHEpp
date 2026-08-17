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
    using bkP = TFHEpp::cblvl02param;
    using privksP = TFHEpp::cblvl21param;

    TFHEpp::SecretKey* sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);
    using targetP = typename privksP::targetP;
    vector<array<uint8_t, targetP::n>> pa(num_test);
    vector<array<typename targetP::T, targetP::n>> pmu(num_test);
    vector<uint8_t> pzeros(num_test);
    array<bool, targetP::n> pres;
    for (array<uint8_t, targetP::n>& i : pa)
        for (uint8_t& p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < targetP::n; j++)
            pmu[i][j] = pa[i][j] ? targetP::μ : -targetP::μ;
    for (int i = 0; i < num_test; i++) pzeros[i] = false;
    vector<TRLWE<targetP>, TFHEpp::AlignedAllocator<TFHEpp::TRLWE<targetP>, 64>>
        ca(num_test);
    vector<TLWE<targetP>> czeros(num_test);
    vector<TRGSWFFT<targetP>,
           TFHEpp::AlignedAllocator<TFHEpp::TRGSWFFT<targetP>, 64>>
        bootedTGSW(num_test);
    vector<TRGSWFFT<targetP>,
           TFHEpp::AlignedAllocator<TFHEpp::TRGSWFFT<targetP>, 64>>
        invbootedTGSW(num_test);

    for (int i = 0; i < num_test; i++)
        trlweSymEncrypt<targetP>(ca[i], pmu[i],
                                 sk->key.getIndependent<targetP>());
    bootsSymEncrypt<targetP>(czeros, pzeros,
                             sk->key.getIndependent<targetP>());

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CircuitBootstrappingWithInv<lvl10param, cblvl02param, privksP>(
            bootedTGSW[test], invbootedTGSW[test], czeros[test], ek);
    }
    end = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        ExternalProduct<targetP>(ca[test], ca[test], invbootedTGSW[test]);
        pres = trlweSymDecrypt<targetP>(ca[test],
                                        sk->key.getIndependent<targetP>());
        for (int i = 0; i < targetP::n; i++) assert(pres[i] == pa[test][i]);
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}
