#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    constexpr uint32_t num_test = 100;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    vector<int32_t> ps(num_test);
    vector<array<uint32_t, lvl1param::n>> p1(num_test);

    vector<array<uint32_t, lvl1param::n>> pmu1(num_test);
    array<bool, lvl1param::n> pres;

    for (int32_t &p : ps) p = binary(engine);
    for (array<uint32_t, lvl1param::n> &i : p1)
        for (uint32_t &p : i) p = binary(engine);

    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < lvl1param::n; j++)
            pmu1[i][j] = (p1[i][j] > 0) ? lvl1param::μ : -lvl1param::μ;
    vector<BootstrappingKeyElementFFT<lvl01param>> cs(num_test);
    vector<TRLWE<lvl1param>> c1(num_test);
    vector<TRLWE<lvl1param>> cres(num_test);

    for (int i = 0; i < num_test; i++) {
        Polynomial<TFHEpp::lvl1param> plainpoly = {};
        plainpoly[0] = ps[i];
        cs[i][TFHEpp::lvl0param::key_value_diff - 1] =
            trgswfftSymEncrypt<lvl1param>(plainpoly, lvl1param::α,
                                          sk->key.lvl1);
    }
    for (int i = 0; i < num_test; i++)
        c1[i] = trlweSymEncrypt<lvl1param>(pmu1[i], lvl1param::α, sk->key.lvl1);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CMUXFFTwithPolynomialMulByXaiMinusOne<lvl01param>(c1[test], cs[test],
                                                          -2);
    }
    end = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        pres = trlweSymDecrypt<lvl1param>(c1[test], sk->key.lvl1);
        TFHEpp::Polynomial<TFHEpp::lvl1param> polyres = pmu1[test];
        if (ps[test] == 1)
            TFHEpp::PolynomialMulByXai<lvl1param>(polyres, pmu1[test],
                                                  2 * lvl1param::n - 2);
        for (int i = 0; i < lvl1param::n; i++) {
            // std::cout<<i<<":"<<ps[test]<<":"<<pres[i]<<":"<<(static_cast<int>(polyres[i])>0?1:0)<<std::endl;
            assert(pres[i] == (static_cast<int>(polyres[i]) > 0 ? 1 : 0));
        }
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count();
    cout << elapsed / num_test << "μs" << endl;
}