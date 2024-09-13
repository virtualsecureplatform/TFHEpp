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
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    std::uniform_int_distribution<uint> ldist(1, TFHEpp::lvl1param::nbit + 1);

    std::cout << "EvalAuto test" << std::endl;
    for (int test = 0; test < num_test; test++) {
        TFHEpp::SecretKey sk;
        TFHEpp::Polynomial<TFHEpp::lvl1param> pa, pmu, autokey;
        for (typename TFHEpp::lvl1param::T &i : pa) i = binary(engine) ? 1 : -1;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            pmu[i] =
                (pa[i] == 1) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;

        const uint d = (1U << ldist(engine)) + 1;

        TFHEpp::TRLWE<TFHEpp::lvl1param> ca =
            TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(pmu, sk.key.lvl1);

        TFHEpp::EvalAutoKey<TFHEpp::lvl1param> eak;
        TFHEpp::evalautokeygen<TFHEpp::lvl1param>(eak, d, sk.key.lvl1);

        TFHEpp::TRLWE<TFHEpp::lvl1param> cres;
        TFHEpp::EvalAuto<TFHEpp::lvl1param>(cres, ca, d, eak);

        TFHEpp::Polynomial<TFHEpp::lvl1param> autopoly;
        TFHEpp::Automorphism<TFHEpp::lvl1param>(autopoly, pa, d);

        std::array<bool, TFHEpp::lvl1param::n> pres =
            TFHEpp::trlweSymDecrypt<TFHEpp::lvl1param>(cres, sk.key.lvl1);

        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            assert(pres[i] == (autopoly[i] == 1));
    }
    std::cout << "PASS" << std::endl;

    return 0;
}