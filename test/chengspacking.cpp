#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    constexpr uint l = 4;
    constexpr uint numtlwe = 1<<l;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());

    std::vector<std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>>> ca(num_test);

    std::vector<std::vector<uint8_t>> pin(num_test);
    for (std::vector<uint8_t> &i : pin){
        i.resize(numtlwe);
        for (uint8_t &p : i) p = binary(engine);
    }
    for (int i = 0; i < num_test; i++) ca[i] = TFHEpp::bootsSymEncrypt<TFHEpp::lvl1param>(pin[i], *sk);

    std::vector<TFHEpp::TRLWE<TFHEpp::lvl1param>> cres(num_test);

    std::unique_ptr<TFHEpp::AnnihilateKey<TFHEpp::lvl1param>> ahk(new TFHEpp::AnnihilateKey<TFHEpp::lvl1param>());
    TFHEpp::annihilatekeygen<TFHEpp::lvl1param>(*ahk, *sk);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::TLWE2TRLWEChengsPacking<TFHEpp::lvl1param>(cres[test], ca[test], *ahk);
    }

    end = std::chrono::system_clock::now();

    for (int i = 0; i < num_test; i++){
        std::array<bool, TFHEpp::lvl1param::n> pres = TFHEpp::trlweSymDecrypt<TFHEpp::lvl1param>(cres[i], sk->key.lvl1);
        for(int j = 0; j < numtlwe; j++) assert(pres[j*(TFHEpp::lvl1param::n>>l)] == (pin[i][j] > 0));
    }

    std::cout << "Passed" << std::endl;
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}