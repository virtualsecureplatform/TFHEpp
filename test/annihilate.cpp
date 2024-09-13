#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    std::vector<uint8_t> pa(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> ca(num_test);
    std::vector<std::array<uint8_t, TFHEpp::lvl1param::n>> pin(num_test);

    std::vector<std::array<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::n>>
        pmu(num_test);

    for (std::array<uint8_t, TFHEpp::lvl1param::n> &i : pin)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < TFHEpp::lvl1param::n; j++)
            pmu[i][j] =
                (pin[i][j] > 0)
                    ? (TFHEpp::lvl1param::μ)
                    : -(TFHEpp::lvl1param::μ);

    std::vector<TFHEpp::TRLWE<TFHEpp::lvl1param>> cin(num_test);
    for (int i = 0; i < num_test; i++)
        cin[i] =
            TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(pmu[i], sk->key.lvl1);

    std::vector<TFHEpp::TRLWE<TFHEpp::lvl1param>> cres(num_test);

    std::unique_ptr<TFHEpp::AnnihilateKey<TFHEpp::lvl1param>> ahk(
        new TFHEpp::AnnihilateKey<TFHEpp::lvl1param>());
    TFHEpp::annihilatekeygen<TFHEpp::lvl1param>(*ahk, *sk);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::AnnihilateKeySwitching<TFHEpp::lvl1param>(cres[test], cin[test],
                                                          *ahk);
    }

    end = std::chrono::system_clock::now();
    std::vector<std::array<bool, TFHEpp::lvl1param::n>> pres(num_test);
    for (int i = 0; i < num_test; i++)
        pres[i] =
            TFHEpp::trlweSymDecrypt<TFHEpp::lvl1param>(cres[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++) assert(pres[i][0] == (pin[i][0] > 0));
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}