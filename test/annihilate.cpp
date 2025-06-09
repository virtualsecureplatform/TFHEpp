#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
import tfhepp;

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    using P = TFHEpp::lvl2param;

    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    std::vector<uint8_t> pa(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    std::vector<TFHEpp::TLWE<P>> ca(num_test);
    std::vector<std::array<uint8_t, P::n>> pin(num_test);

    std::vector<std::array<typename P::T, P::n>> pmu(num_test);

    for (std::array<uint8_t, P::n> &i : pin)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < P::n; j++)
            pmu[i][j] = (pin[i][j] > 0) ? (P::μ) : -(P::μ);

    std::vector<TFHEpp::TRLWE<P>> cin(num_test);
    for (int i = 0; i < num_test; i++)
        cin[i] = TFHEpp::trlweSymEncrypt<P>(pmu[i], sk->key.get<P>());

    std::vector<TFHEpp::TRLWE<P>> cres(num_test);

    std::unique_ptr<TFHEpp::AnnihilateKey<P>> ahk(
        new TFHEpp::AnnihilateKey<P>());
    TFHEpp::annihilatekeygen<P>(*ahk, *sk);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::AnnihilateKeySwitching<P>(cres[test], cin[test], *ahk);
    }

    end = std::chrono::system_clock::now();
    std::vector<std::array<bool, P::n>> pres(num_test);
    for (int i = 0; i < num_test; i++)
        pres[i] = TFHEpp::trlweSymDecrypt<P>(cres[i], sk->key.get<P>());
    for (int i = 0; i < num_test; i++) assert(pres[i][0] == (pin[i][0] > 0));
    // TFHEpp::Polynomial<P> phase =
    // TFHEpp::trlwePhase<P>(cres[0], sk->key.get<P>());
    // for (int i = 0; i < P::n; i++)
    // std::cout << static_cast<int64_t>(phase[i]) << ":";
    // std::cout << std::endl;
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}