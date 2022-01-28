#include <cassert>
#include <iostream>
#include <random>
#include <chrono>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 10;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkntt<TFHEpp::lvl01param>(sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    std::array<TFHEpp::TLWE<TFHEpp::lvl1param>,num_test> tlwe,bootedtlwe;
    std::array<bool,num_test> p;
    for(int i = 0; i < num_test; i++) p[i] = binary(engine) > 0;
    for(int i = 0; i < num_test; i++) tlwe[i] =
            TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
                p[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
                TFHEpp::lvl1param::α, sk.key.lvl1);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::GateBootstrappingNTT(bootedtlwe[test], tlwe[test], ek);
    }

    end = std::chrono::system_clock::now();
    for(int i = 0; i < num_test; i++){
        bool p2 =
                TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(bootedtlwe[i], sk.key.lvl1);
        assert(p[i] == p2);
    }
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}