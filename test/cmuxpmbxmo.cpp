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
    uniform_int_distribution<uint> exponentgen(0, 2 * TFHEpp::lvl1param::n - 1);

    SecretKey *sk = new SecretKey;
    std::vector<int32_t> ps(num_test);
    std::vector<array<typename TFHEpp::lvl1param::T, lvl1param::n>> p1(
        num_test);

    std::vector<std::array<typename TFHEpp::lvl1param::T, lvl1param::n>> pmu1(
        num_test);
    std::vector<uint> exponents(num_test);
    std::array<bool, TFHEpp::lvl1param::n> pres;

    for (int32_t &p : ps) p = binary(engine);
    for (uint &p : exponents) p = exponentgen(engine);
    for (array<typename TFHEpp::lvl1param::T, lvl1param::n> &i : p1)
        for (typename TFHEpp::lvl1param::T &p : i) p = binary(engine);

    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < lvl1param::n; j++)
            pmu1[i][j] = (p1[i][j] > 0) ? lvl1param::μ : -lvl1param::μ;
    std::vector<TFHEpp::BootstrappingKeyElementFFT<TFHEpp::lvl01param>,
                TFHEpp::AlignedAllocator<
                    TFHEpp::BootstrappingKeyElementFFT<TFHEpp::lvl01param>, 64>>
        cs(num_test);
    std::vector<TRLWE<lvl1param>> c1(num_test);
    std::vector<TRLWE<lvl1param>> cres(num_test);

    for (int i = 0; i < num_test; i++) {
        TFHEpp::Polynomial<TFHEpp::lvl1param> plainpoly = {};
        plainpoly[0] = ps[i];
        cs[i][TFHEpp::lvl0param::key_value_diff - 1] =
            trgswfftSymEncrypt<lvl1param>(plainpoly, sk->key.lvl1);
    }
    for (int i = 0; i < num_test; i++)
        c1[i] = trlweSymEncrypt<lvl1param>(pmu1[i], sk->key.lvl1);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CMUXFFTwithPolynomialMulByXaiMinusOne<lvl01param>(c1[test], cs[test],
                                                          exponents[test]);
    }
    end = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        pres = trlweSymDecrypt<lvl1param>(c1[test], sk->key.lvl1);
        TFHEpp::Polynomial<TFHEpp::lvl1param> polyres = pmu1[test];
        if (ps[test] == 1)
            TFHEpp::PolynomialMulByXai<lvl1param>(polyres, pmu1[test],
                                                  exponents[test]);
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