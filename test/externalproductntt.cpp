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

    std::cout << "test p=1" << std::endl;

    std::cout << "lvl1" << std::endl;
    for (int test = 0; test < num_test; test++) {
        TFHEpp::lweKey key;

        std::array<bool, TFHEpp::lvl1param::n> p;
        for (bool &i : p) i = (binary(engine) > 0);
        TFHEpp::Polynomial<TFHEpp::lvl1param> pmu;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            pmu[i] = p[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;
        TFHEpp::TRLWE<TFHEpp::lvl1param> c =
            TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(pmu, key.lvl1);

        const TFHEpp::Polynomial<TFHEpp::lvl1param> plainpoly = {
            static_cast<typename TFHEpp::lvl1param::T>(1)};

        TFHEpp::TRGSWNTT<TFHEpp::lvl1param> trgswntt =
            TFHEpp::trgswnttSymEncrypt<TFHEpp::lvl1param>(plainpoly, key.lvl1);
        TFHEpp::trgswnttExternalProduct<TFHEpp::lvl1param>(c, c, trgswntt);
        std::array<bool, TFHEpp::lvl1param::n> p2 =
            TFHEpp::trlweSymDecrypt<TFHEpp::lvl1param>(c, key.lvl1);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) assert(p[i] == p2[i]);
    }
    std::cout << "Passed" << std::endl;

    std::cout << "test p=-1" << std::endl;

    std::cout << "lvl1" << std::endl;
    for (int test = 0; test < num_test; test++) {
        TFHEpp::lweKey key;

        std::array<bool, TFHEpp::lvl1param::n> p;
        for (bool &i : p) i = binary(engine) > 0;
        std::array<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::n> pmu;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            pmu[i] = p[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;
        TFHEpp::TRLWE<TFHEpp::lvl1param> c =
            TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(pmu, key.lvl1);

        const TFHEpp::Polynomial<TFHEpp::lvl1param> plainpoly = {
            static_cast<typename TFHEpp::lvl1param::T>(-1)};

        TFHEpp::TRGSWNTT<TFHEpp::lvl1param> trgswntt =
            TFHEpp::trgswnttSymEncrypt<TFHEpp::lvl1param>(plainpoly, key.lvl1);
        TFHEpp::trgswnttExternalProduct<TFHEpp::lvl1param>(c, c, trgswntt);
        std::array<bool, TFHEpp::lvl1param::n> p2 =
            TFHEpp::trlweSymDecrypt<TFHEpp::lvl1param>(c, key.lvl1);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) assert(p[i] == !p2[i]);
    }
    std::cout << "Passed" << std::endl;
}