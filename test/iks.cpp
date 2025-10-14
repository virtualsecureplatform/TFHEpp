#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> tlwe(num_test);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl0param>> res(num_test);
    std::array<bool, num_test> p;
    for (int i = 0; i < num_test; i++) p[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++)
        tlwe[i] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            p[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            sk.key.get<TFHEpp::lvl1param>());

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl10param>(
            res[test], tlwe[test], ek.getiksk<TFHEpp::lvl10param>());
    }
    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
    for (int i = 0; i < num_test; i++) {
        bool p2 = TFHEpp::tlweSymDecrypt<TFHEpp::lvl0param>(
            res[i], sk.key.get<TFHEpp::lvl0param>());
        assert(p[i] == p2);
    }
    std::cout << "Passed" << std::endl;
}