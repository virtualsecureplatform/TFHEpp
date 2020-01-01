#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    const uint32_t num_test = 10;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CloudKey *ck = new CloudKey(*sk);
    vector<array<bool,DEF_N>> pa(num_test);
    vector<array<uint32_t,DEF_N>> pmu(num_test);
    vector<bool> pones(num_test);
    array<bool,DEF_N> pres;
    for (array<bool,DEF_N> &i : pa) for(bool &p:i) p = binary(engine) > 0;
    for (int i = 0;i<num_test;i++) for(int j = 0;j<DEF_N;j++) pmu[i][j] = pa[i][j]?DEF_μ:-DEF_μ;
    for (int i = 0; i < num_test; i++) pones[i] = true;
    vector<TRLWElvl1> ca(num_test);
    vector<TLWElvl0> cones(num_test);
    vector<TRGSWFFTlvl1> bootedTGSW(num_test);

    for(int i = 0;i<num_test;i++) ca[i] = trlweSymEncryptlvl1(pmu[i], DEF_αbk,sk->key.lvl1);
    cones = bootsSymEncrypt(pones, *sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CircuitBootstrappingFFT(bootedTGSW[test],cones[test],*ck);
    }
    end = chrono::system_clock::now();
    
    for (int test = 0; test < num_test; test++) {
        trgswfftExternalProductlvl1(ca[test], bootedTGSW[test]);
        pres = trlweSymDecryptlvl1(ca[test],sk->key.lvl1);
        for (int i = 0; i < num_test; i++) assert(pres[i] == pa[test][i]);
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}