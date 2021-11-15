#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
    const uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey* sk = new TFHEpp::SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    std::vector<uint8_t> pa(num_test);
    std::vector<uint8_t> pb(num_test);
    std::vector<uint8_t> pand(num_test);
    std::vector<uint8_t> pnor(num_test);
    std::vector<uint8_t> pxor(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pb[i] = binary(engine) > 0;
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> ca(num_test);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> cb(num_test);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> carry(num_test);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> sum(num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            pa[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        cb[i] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            pb[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::HomHalfAdder<TFHEpp::lvl1param>(carry[test], sum[test], ca[test], cb[test], ek);
    }

    end = std::chrono::system_clock::now();
    for (int i = 0; i < num_test; i++)
        pand[i] =
            TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(carry[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        pxor[i] =
            TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(sum[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++) assert(pand[i] == (pa[i] & pb[i]));
    for (int i = 0; i < num_test; i++) assert(pxor[i] == (pa[i] ^ pb[i]));
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}
