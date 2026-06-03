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

#ifdef USE_INTERLEAVED_FORMAT
    std::cout << "USE_INTERLEAVED_FORMAT" << std::endl;
#endif

    cout << "Start LVL1 test." << endl;
    for (int test = 0; test < num_test; test++) {
        Polynomial<lvl1param> a;
        for (typename TFHEpp::lvl1param::T &i : a) i = Torus32dist(engine);
        PolynomialInFD<lvl1param> resfft;
        TFHEpp::TwistIFFT<lvl1param>(resfft, a);
        Polynomial<lvl1param> res;
        TFHEpp::TwistFFT<lvl1param>(res, resfft);
        for (int i = 0; i < lvl1param::n; i++)
            assert(abs(static_cast<int32_t>(a[i] - res[i])) <= 1);
    }
    cout << "FFT Passed" << endl;

    for (int test = 0; test < num_test; test++) {
        alignas(64) array<typename TFHEpp::lvl1param::T, lvl1param::n> a;
        for (int i = 0; i < lvl1param::n; i++)
            a[i] = Bgdist(engine) - lvl1param::Bg / 2;
        for (typename TFHEpp::lvl1param::T &i : a)
            i = Bgdist(engine) - lvl1param::Bg / 2;
        alignas(64) array<typename TFHEpp::lvl1param::T, lvl1param::n> b;
        for (typename TFHEpp::lvl1param::T &i : b) i = Torus32dist(engine);

        alignas(64) Polynomial<lvl1param> polymul;
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

    std::cout << "PolyMulRescale Test" << std::endl;
    for (int test = 0; test < num_test; test++) {
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<typename TFHEpp::lvl1param::T> message(
            0, (1ULL << 32) - 1);

        TFHEpp::Polynomial<TFHEpp::lvl1param> p0, p1, pres, ptrue;
        for (typename TFHEpp::lvl1param::T &i : p0) i = message(engine);
        for (typename TFHEpp::lvl1param::T &i : p1) i = message(engine);

        TFHEpp::PolyMulRescaleUnsigned<TFHEpp::lvl1param>(pres, p0, p1);
        TFHEpp::PolyMulNaieveRescaleUnsigned<TFHEpp::lvl1param>(ptrue, p0, p1);

        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            // std::cout<<pres[i]<<":"<<ptrue[i]<<std::endl;
            assert(abs(static_cast<int>(pres[i] - ptrue[i])) <= 2);
        }
    }
    std::cout << "PolyMulRescaleUnsigned Passed" << std::endl;
    return 0;
}