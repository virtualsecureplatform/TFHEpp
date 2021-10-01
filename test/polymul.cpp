#include <algorithm>
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
    uniform_int_distribution<uint32_t> Bgdist(0, lvl1param::Bg);
    uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);

    cout << "Start LVL1 test." << endl;
    for (int test = 0; test < num_test; test++) {
        Polynomial<lvl1param> a;
        for (uint32_t &i : a) i = Torus32dist(engine);
        PolynomialInFD<lvl1param> resfft;
        TFHEpp::TwistIFFT<lvl1param>(resfft, a);
        Polynomial<lvl1param> res;
        TFHEpp::TwistFFT<lvl1param>(res, resfft);
        for (int i = 0; i < lvl1param::n; i++)
            assert(abs(static_cast<int32_t>(a[i] - res[i])) <= 1);
    }
    cout << "FFT Passed" << endl;

    for (int test = 0; test < num_test; test++) {
        array<uint32_t, lvl1param::n> a;
        for (int i = 0; i < lvl1param::n; i++)
            a[i] = Bgdist(engine) - lvl1param::Bg / 2;
        for (uint32_t &i : a) i = Bgdist(engine) - lvl1param::Bg / 2;
        array<uint32_t, lvl1param::n> b;
        for (uint32_t &i : b) i = Torus32dist(engine);

        Polynomial<lvl1param> polymul;
        TFHEpp::PolyMul<lvl1param>(polymul, a, b);
        Polynomial<lvl1param> naieve = {};
        for (int i = 0; i < lvl1param::n; i++) {
            for (int j = 0; j <= i; j++)
                naieve[i] += static_cast<int32_t>(a[j]) * b[i - j];
            for (int j = i + 1; j < lvl1param::n; j++)
                naieve[i] -=
                    static_cast<int32_t>(a[j]) * b[lvl1param::n + i - j];
        }
        for (int i = 0; i < lvl1param::n; i++)
            assert(abs(static_cast<int32_t>(naieve[i] - polymul[i])) <= 1);
    }
    cout << "PolyMul Passed" << endl;

    uniform_int_distribution<uint64_t> Bgbardist(0, lvl2param::Bg);
    uniform_int_distribution<uint64_t> Torus64dist(0, UINT64_MAX);

    cout << "Start LVL2 test." << endl;
    for (int test = 0; test < num_test; test++) {
        Polynomial<lvl2param> a;
        for (uint64_t &i : a) i = Torus64dist(engine);
        PolynomialInFD<lvl2param> resfft;
        TFHEpp::TwistIFFT<lvl2param>(resfft, a);
        Polynomial<lvl2param> res;
        TFHEpp::TwistFFT<lvl2param>(res, resfft);
        for (int i = 0; i < lvl2param::n; i++)
            assert(abs(static_cast<int64_t>(a[i] - res[i])) <= (1 << 14));
    }
    cout << "FFT Passed" << endl;

    return 0;
}