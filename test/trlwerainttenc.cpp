#include <array>
#include <cassert>
#include <iostream>
#include <tfhe++.hpp>

int main()
{
    std::cout << "lvl1" << std::endl;
    // static_assert(TFHEpp::hasq<TFHEpp::lvl1param>);
    // std::cout<< TFHEpp::lvl1param::q<<std::endl;
    constexpr uint32_t num_test = 1000;
    if constexpr (TFHEpp::hasq<TFHEpp::lvl1param>)
        for (int test = 0; test < num_test; test++) {
            std::random_device seed_gen;
            std::default_random_engine engine(seed_gen());
            std::uniform_int_distribution<typename TFHEpp::lvl1param::T> binary(
                0, 1);

            TFHEpp::lweKey key;
            std::array<bool, TFHEpp::lvl1param::n> p;
            for (bool &i : p) i = binary(engine) > 0;
            std::array<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::n>
                pmu = {};
            for (int i = 0; i < TFHEpp::lvl1param::n; i++)
                pmu[i] = p[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;
            TFHEpp::TRLWERAINTT<TFHEpp::lvl1param> craintt =
                TFHEpp::trlwerainttSymEncrypt<TFHEpp::lvl1param>(pmu, 3,
                                                                 key.lvl1);
            TFHEpp::TRLWE<TFHEpp::lvl1param> c;
            for (int k = 0; k <= TFHEpp::lvl1param::k; k++) {
                raintt::TwistNTT<typename TFHEpp::lvl1param::T,
                                 TFHEpp::lvl1param::nbit, false>(
                    c[k], craintt[k], (*TFHEpp::raintttable)[0],
                    (*TFHEpp::raintttwist)[0]);
                for (int i = 0; i < TFHEpp::lvl1param::n; i++)
                    c[k][i] = static_cast<int32_t>(c[k][i]) < 0
                                  ? c[k][i] + raintt::P
                                  : c[k][i];
            }
            std::array<bool, TFHEpp::lvl1param::n> p2 =
                TFHEpp::trlweSymDecrypt<TFHEpp::lvl1param>(c, key.lvl1);
            // for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            // std::cout<<test<<":"<<i<<":"<<c[TFHEpp::lvl1param::k][i]<<":"<<(p[i]?1:0)<<":"<<(p2[i]?1:0)<<std::endl;
            for (int i = 0; i < TFHEpp::lvl1param::n; i++)
                assert(p[i] == p2[i]);
        }
    else
        std::cout << "Nothing to do" << std::endl;
    std::cout << "Passed" << std::endl;
}