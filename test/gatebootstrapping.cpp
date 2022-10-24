#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 10;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    using bkP = TFHEpp::lvl02param;
    using iksP = TFHEpp::lvl20param;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<bkP>(sk);
    ek.emplaceiksk<iksP>(sk);
    std::array<TFHEpp::TLWE<typename iksP::domainP>, num_test> tlwe, bootedtlwe;
    std::array<bool, num_test> p;
    for (int i = 0; i < num_test; i++) p[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++)
        tlwe[i] = TFHEpp::tlweSymEncrypt<typename iksP::domainP>(
            p[i] ? iksP::domainP::μ : -iksP::domainP::μ, iksP::domainP::α,
            sk.key.get<typename iksP::domainP>());

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::GateBootstrapping<iksP, bkP>(bootedtlwe[test], tlwe[test], ek);
    }

    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
    for (int i = 0; i < num_test; i++) {
        bool p2 = TFHEpp::tlweSymDecrypt<typename iksP::domainP>(
            bootedtlwe[i], sk.key.get<typename iksP::domainP>());
        assert(p[i] == p2);
    }
    std::cout << "Passed" << std::endl;
}