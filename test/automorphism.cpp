#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> Bgdist(0, TFHEpp::lvl1param::Bg);
    std::uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);
    std::uniform_int_distribution<uint> ldist(1, TFHEpp::lvl1param::nbit + 1);

    std::cout << "Automorphism test" << std::endl;
    for (int test = 0; test < num_test; test++) {
        TFHEpp::Polynomial<TFHEpp::lvl1param> a, b;
        for (typename TFHEpp::lvl1param::T &i : a)
            i = Bgdist(engine) - TFHEpp::lvl1param::Bg / 2;
        for (typename TFHEpp::lvl1param::T &i : b) i = Torus32dist(engine);

        const uint d = (1U << ldist(engine)) + 1;

        TFHEpp::Polynomial<TFHEpp::lvl1param> beforeautomul, autoaftermul;
        TFHEpp::PolyMul<TFHEpp::lvl1param>(beforeautomul, a, b);
        TFHEpp::Automorphism<TFHEpp::lvl1param>(autoaftermul, beforeautomul, d);

        TFHEpp::Polynomial<TFHEpp::lvl1param> autoa, autob;
        TFHEpp::Automorphism<TFHEpp::lvl1param>(autoa, a, d);
        TFHEpp::Automorphism<TFHEpp::lvl1param>(autob, b, d);

        TFHEpp::Polynomial<TFHEpp::lvl1param> mulafterauto;
        TFHEpp::PolyMul<TFHEpp::lvl1param>(mulafterauto, autoa, autob);

        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            assert(abs(static_cast<int32_t>(mulafterauto[i] -
                                            autoaftermul[i])) <= 1);
    }
    std::cout << "PASS" << std::endl;
}