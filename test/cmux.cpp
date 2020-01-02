#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    const uint32_t num_test = 100;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    vector<int32_t> ps(num_test);
    vector<array<uint8_t, DEF_N>> p1(num_test);
    vector<array<uint8_t, DEF_N>> p0(num_test);
    
    vector<array<uint32_t, DEF_N>> pmu1(num_test);
    vector<array<uint32_t, DEF_N>> pmu0(num_test);
    array<bool, DEF_N> pres;

    for (int32_t &p : ps) p = binary(engine);
    for (array<uint8_t, DEF_N> &i : p1)
        for (uint8_t &p : i) p = binary(engine);
    for (array<uint8_t, DEF_N> &i : p0)
        for (uint8_t &p : i) p = binary(engine);

    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < DEF_N; j++) pmu1[i][j] = (p1[i][j]>0) ? DEF_μ : -DEF_μ;
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < DEF_N; j++) pmu0[i][j] = (p0[i][j]>0) ? DEF_μ : -DEF_μ;
    vector<TRGSWFFTlvl1> cs(num_test);
    vector<TRLWElvl1> c1(num_test);
    vector<TRLWElvl1> c0(num_test);
    vector<TRLWElvl1> cres(num_test);

    for (int i = 0; i < num_test; i++)
        cs[i] = trgswfftSymEncryptlvl1(ps[i], DEF_αbk, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        c1[i] = trlweSymEncryptlvl1(pmu1[i], DEF_αbk, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        c0[i] = trlweSymEncryptlvl1(pmu0[i], DEF_αbk, sk->key.lvl1);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CMUXFFTlvl1(cres[test],cs[test], c1[test], c0[test]);
    }
    end = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        pres = trlweSymDecryptlvl1(cres[test], sk->key.lvl1);
        for (int i = 0; i < DEF_N; i++) assert(pres[i] == ((ps[test]>0)?p1[test][i]:p0[test][i]) > 0);
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count();
    cout << elapsed / num_test << "μs" << endl;
}