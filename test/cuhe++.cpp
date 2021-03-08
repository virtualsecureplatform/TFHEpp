#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <tfhe++.hpp>
#include <cuhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> Bgdist(0, TFHEpp::lvl1param::Bg);
    std::uniform_int_distribution<uint32_t> Torus32dist(0, std::numeric_limits<typename TFHEpp::lvl1param::T>::max());

    std::array<std::array<cuHEpp::INTorus,TFHEpp::lvl1param::n>,2> twistlvl1;
    std::array<std::array<cuHEpp::INTorus,TFHEpp::lvl1param::n/2>,2> tablelvl1;
    twistlvl1 = cuHEpp::TwistGen<TFHEpp::lvl1param::nbit>();
    tablelvl1 = cuHEpp::TableGen<TFHEpp::lvl1param::nbit>();

    for(int i = 0;i<TFHEpp::lvl1param::n;i++) assert((twistlvl1[0][i]*twistlvl1[1][i]).value==1);
    for(int i = 1;i<TFHEpp::lvl1param::n;i++) assert(twistlvl1[0][i].value!=1);
    for(int i = 0;i<TFHEpp::lvl1param::n/2;i++) assert((tablelvl1[0][i]*tablelvl1[1][i]).value==1);
    for(int i = 1;i<TFHEpp::lvl1param::n/2;i++) assert(tablelvl1[0][i].value!=1);
    assert((tablelvl1[1][TFHEpp::lvl1param::n/2-1]*tablelvl1[1][1]).value == cuHEpp::P-1);

    std::cout << "Start NTT only test." << std::endl;
    for (int test; test < num_test; test++) {
        // std::array<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::n> a,res;
        std::array<cuHEpp::INTorus,TFHEpp::lvl1param::n> res;
        TFHEpp::Polynomial<TFHEpp::lvl1param> a;
        for (uint32_t &i : a) i = Torus32dist(engine);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) res[i] = cuHEpp::INTorus(a[i]);
        cuHEpp::INTT<TFHEpp::lvl1param::nbit,0>(res.begin(), tablelvl1[1]);
        cuHEpp::NTT<TFHEpp::lvl1param::nbit,0>(res.begin(), tablelvl1[0]);

        const cuHEpp::INTorus invN = cuHEpp::InvPow2(TFHEpp::lvl1param::nbit);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            res[i] *= invN;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            std::cout<<res[i].value<<":"<<a[i]<<std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            assert(a[i] == res[i].value);
    }
    std::cout << "Start NTT only test PASS" << std::endl;

    // std::cout << "Start LVL1 test." << std::endl;
    // for (int test; test < num_test; test++) {
    //     // std::array<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::n> a,res;
    //     TFHEpp::Polynomial<TFHEpp::lvl1param> a,res;
    //     for (uint32_t &i : a) i = Torus32dist(engine);
    //     std::array<cuHEpp::INTorus,TFHEpp::lvl1param::n> resntt;
    //     cuHEpp::TwistINTTlvl1<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::nbit>(resntt, a, tablelvl1[1], twistlvl1[1]);
    //     cuHEpp::TwistNTTlvl1<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::nbit>(res, resntt, tablelvl1[0], twistlvl1[0]);
    //     for (int i = 0; i < TFHEpp::lvl1param::n; i++)
    //         std::cout<<res[i]<<":"<<a[i]<<std::endl;
    //     for (int i = 0; i < TFHEpp::lvl1param::n; i++)
    //         assert(a[i] == res[i]);
    // }
    // cout << "NTT Passed" << endl;

    // for (int test; test < num_test; test++) {
    //     array<uint32_t, DEF_N> a;
    //     for (int i = 0; i < DEF_N; i++) a[i] = Bgdist(engine) - DEF_Bg / 2;
    //     for (uint32_t &i : a) i = Bgdist(engine) - DEF_Bg / 2;
    //     array<uint32_t, DEF_N> b;
    //     for (uint32_t &i : b) i = Torus32dist(engine);

    //     Polynomiallvl1 polymul;
    //     SPQLIOSpp::PolyMullvl1(polymul, a, b);
    //     Polynomiallvl1 naieve = {};
    //     for (int i = 0; i < DEF_N; i++) {
    //         for (int j = 0; j <= i; j++)
    //             naieve[i] += static_cast<int32_t>(a[j]) * b[i - j];
    //         for (int j = i + 1; j < DEF_N; j++)
    //             naieve[i] -= static_cast<int32_t>(a[j]) * b[DEF_N + i - j];
    //     }
    //     for (int i = 0; i < DEF_N; i++)
    //         assert(abs(static_cast<int32_t>(naieve[i] - polymul[i])) <= 1);
    // }
    // cout << "PolyMul Passed" << endl;

    // uniform_int_distribution<uint64_t> Bgbardist(0, DEF_Bgbar);
    // uniform_int_distribution<uint64_t> Torus64dist(0, UINT64_MAX);

    // cout << "Start LVL2 test." << endl;
    // for (int test = 0; test < num_test; test++) {
    //     Polynomiallvl2 a;
    //     for (uint64_t &i : a) i = Torus64dist(engine);
    //     PolynomialInFDlvl2 resfft;
    //     TFHEpp::TwistIFFTlvl2(resfft, a);
    //     Polynomiallvl2 res;
    //     TFHEpp::TwistFFTlvl2(res, resfft);
    //     for (int i = 0; i < DEF_N; i++)
    //         assert(abs(static_cast<int64_t>(a[i] - res[i])) <= (1 << 14));
    // }
    // cout << "FFT Passed" << endl;

    return 0;
}