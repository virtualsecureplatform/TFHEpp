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

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    for (int test = 0; test < num_test; test++) {
        bool p = binary(engine) > 0;
        TFHEpp::TLWE<TFHEpp::lvl1param> tlwe =
            TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
                p ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ, sk.key.lvl1);
        TFHEpp::TLWE<TFHEpp::lvl0param> res;
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl10param>(res, tlwe, *ek.iksklvl10);
        bool p2 = TFHEpp::tlweSymDecrypt<TFHEpp::lvl0param>(res, sk.key.lvl0);
        assert(p == p2);
    }
    std::cout << "Passed" << std::endl;
}