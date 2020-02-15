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

    std::cout << "Start LVL1 test." << endl;

    for (volatile int test; test < num_test; test++) {
        PolynomialInFDlvl1 a = {};
        Polynomiallvl1 b,res,naieve;
        for (uint32_t &i : b) i = Torus32dist(engine);
        uint32_t shift = Torus32dist(engine) & (2*DEF_N -1);

        Polynomiallvl1 polymul;
        PolynomialInFDlvl1 polymulfft,bfft;
        TwistIFFTlvl1(bfft,b);

        PolynomialMulByXaiSubFFTlvl1(polymulfft, bfft, a ,shift);
        TwistFFTlvl1(res,bfft);
        PolynomialMulByXailvl1(naieve,b,shift);
        for(int i = 0;i<DEF_N;i++)assert(naieve[i] == res[i]);
    }
    std::cout << "PolyMul Passed" << endl;

    return 0;
}