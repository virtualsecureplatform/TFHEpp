#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl1param;
    constexpr uint32_t num_test = 1000;
    constexpr uint l = 4;
    constexpr uint numtlwe = 1 << l;
    constexpr uint segment = 1 << (P::nbit - l);
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    constexpr uint32_t plain_modulus = 1 << (l + 1);
    std::uniform_int_distribution<uint32_t> intgen(0, plain_modulus - 1);

    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());

    std::vector<std::array<TFHEpp::TLWE<P>, numtlwe>> ca(num_test);

    std::vector<std::vector<uint8_t>> pin(num_test);
    for (std::vector<uint8_t> &i : pin) {
        i.resize(numtlwe);
        for (uint8_t &p : i) p = intgen(engine);
    }
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < numtlwe; j++)
            ca[i][j] =
                TFHEpp::tlweSymIntEncrypt<P, plain_modulus>(pin[i][j], *sk);

    std::vector<TFHEpp::TRLWE<P>> cres(num_test);

    std::unique_ptr<TFHEpp::AnnihilateKey<P>> ahk(
        new TFHEpp::AnnihilateKey<P>());
    TFHEpp::annihilatekeygen<P>(*ahk, *sk);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::TLWE2TablePacking<P, numtlwe>(cres[test], ca[test], *ahk);
    }

    end = std::chrono::system_clock::now();

    for (int i = 0; i < num_test; i++) {
        TFHEpp::Polynomial<P> pres;
        pres = TFHEpp::trlweSymIntDecrypt<P, plain_modulus>(cres[i], *sk);
        // for (int j = 0; j < numtlwe; j++)
        //     for (int k = 0; k < segment; k++)
        //         std::cout << static_cast<int64_t>(pres[j * segment + k]) <<
        //         ":" << static_cast<int64_t>(pin[i][j]) << std::endl;
        for (int j = 0; j < numtlwe; j++)
            for (int k = 0; k < segment; k++)
                assert(pres[j * segment + k] == pin[i][j]);
    }

    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}