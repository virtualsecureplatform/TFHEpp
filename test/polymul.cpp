#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    const uint32_t num_test = 1000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> Bgdist(0, DEF_Bg);
    uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);

    cout << "Start LVL1 test." << endl;
    for (volatile int test; test < num_test; test++) {
        Polynomiallvl1 a;
        for (uint32_t &i : a) i = Torus32dist(engine);
        PolynomialInFDlvl1 resfft;
        TFHEpp::TwistIFFTlvl1(resfft, a);
        Polynomiallvl1 res;
        TFHEpp::TwistFFTlvl1(res, resfft);
        for (int i = 0; i < DEF_N; i++)
            assert(abs(static_cast<int32_t>(a[i] - res[i])) <= 1);
    }
    cout << "FFT Passed" << endl;

    for (volatile int test; test < num_test; test++) {
        array<uint32_t, DEF_N> a;
        for (int i = 0; i < DEF_N; i++) a[i] = Bgdist(engine) - DEF_Bg / 2;
        for (uint32_t &i : a) i = Bgdist(engine) - DEF_Bg / 2;
        array<uint32_t, DEF_N> b;
        for (uint32_t &i : b) i = Torus32dist(engine);

        Polynomiallvl1 polymul;
        TFHEpp::PolyMullvl1(polymul, a, b);
        Polynomiallvl1 naieve = {};
        for (int i = 0; i < DEF_N; i++) {
            for (int j = 0; j <= i; j++)
                naieve[i] += static_cast<int32_t>(a[j]) * b[i - j];
            for (int j = i + 1; j < DEF_N; j++)
                naieve[i] -= static_cast<int32_t>(a[j]) * b[DEF_N + i - j];
        }
        for (int i = 0; i < DEF_N; i++)
            assert(abs(static_cast<int32_t>(naieve[i] - polymul[i])) <= 1);
    }
    cout << "PolyMul Passed" << endl;

    uniform_int_distribution<uint64_t> Bgbardist(0, DEF_Bgbar);
    uniform_int_distribution<uint64_t> Torus64dist(0, UINT64_MAX);

    cout << "Start LVL2 test." << endl;
    for (int test = 0; test < num_test; test++) {
        Polynomiallvl2 a;
        for (uint64_t &i : a) i = Torus64dist(engine);
        PolynomialInFDlvl2 resfft;
        TFHEpp::TwistIFFTlvl2(resfft, a);
        Polynomiallvl2 res;
        TFHEpp::TwistFFTlvl2(res, resfft);
        for (int i = 0; i < DEF_N; i++)
            assert(abs(static_cast<int64_t>(a[i] - res[i])) <= (1 << 14));
    }
    cout << "FFT Passed" << endl;

    return 0;
}