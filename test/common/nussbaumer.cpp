#include <algorithm>
#include <cassert>
#include <iostream>
#include <nussbaumer.hpp>
#include <random>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> Bgdist(0, TFHEpp::lvl1param::Bg);
    std::uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);

    // std::cout << "Start LVL1 test." << std::endl;
    for (int test = 0; test < num_test; test++) {
        using T = uint64_t;
        std::array<T, TFHEpp::lvl1param::n> a, res;
        for (T &i : a) i = Torus32dist(engine);
        res = a;
        Nussbaumer::NussbaumerTransform<T, TFHEpp::lvl1param::nbit>(
            std::span{res});
        Nussbaumer::InverseNussbaumerTransform<T, TFHEpp::lvl1param::nbit>(
            std::span{res});
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            assert(abs(static_cast<int32_t>(
                           a[i] - res[i] / TFHEpp::lvl1param::n) <= 1));
    }
    std::cout << "Id Passed" << std::endl;

    // for (int test = 0; test < num_test; test++) {
    //     array<typename TFHEpp::lvl1param::T, lvl1param::n> a;
    //     for (int i = 0; i < lvl1param::n; i++)
    //         a[i] = Bgdist(engine) - lvl1param::Bg / 2;
    //     for (typename TFHEpp::lvl1param::T &i : a)
    //         i = Bgdist(engine) - lvl1param::Bg / 2;
    //     array<typename TFHEpp::lvl1param::T, lvl1param::n> b;
    //     for (typename TFHEpp::lvl1param::T &i : b) i = Torus32dist(engine);

    //     Polynomial<lvl1param> polymul;
    //     TFHEpp::PolyMul<lvl1param>(polymul, a, b);
    //     Polynomial<lvl1param> naieve = {};
    //     for (int i = 0; i < lvl1param::n; i++) {
    //         for (int j = 0; j <= i; j++)
    //             naieve[i] += static_cast<int32_t>(a[j]) * b[i - j];
    //         for (int j = i + 1; j < lvl1param::n; j++)
    //             naieve[i] -=
    //                 static_cast<int32_t>(a[j]) * b[lvl1param::n + i - j];
    //     }
    //     for (int i = 0; i < lvl1param::n; i++)
    //         assert(abs(static_cast<int32_t>(naieve[i] - polymul[i])) <= 1);
    // }
    // cout << "PolyMul Passed" << endl;

    // uniform_int_distribution<uint64_t> Bgbardist(0, lvl2param::Bg);
    // uniform_int_distribution<uint64_t> Torus64dist(0, UINT64_MAX);

    // cout << "Start LVL2 test." << endl;
    // for (int test = 0; test < num_test; test++) {
    //     Polynomial<lvl2param> a;
    //     for (uint64_t &i : a) i = Torus64dist(engine);
    //     PolynomialInFD<lvl2param> resfft;
    //     TFHEpp::TwistIFFT<lvl2param>(resfft, a);
    //     Polynomial<lvl2param> res;
    //     TFHEpp::TwistFFT<lvl2param>(res, resfft);
    //     for (int i = 0; i < lvl2param::n; i++)
    //         assert(abs(static_cast<int64_t>(a[i] - res[i])) <= (1 << 14));
    // }
    // cout << "FFT Passed" << endl;

    return 0;
}