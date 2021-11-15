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
    std::vector<uint8_t> pc(num_test);
    std::vector<uint8_t> pcarry(num_test);
    std::vector<uint8_t> psum(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pb[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pc[i] = binary(engine) > 0;
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> ca(num_test);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> cb(num_test);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> cc(num_test);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> csum(num_test);
    std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>> ccarry(num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            pa[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        cb[i] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            pb[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        cc[i] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            pc[i] ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::Hom2BRFullAdder(ccarry[test], csum[test], ca[test], cb[test],
                                cc[test], ek);
    }

    end = std::chrono::system_clock::now();
    for (int i = 0; i < num_test; i++)
        pcarry[i] =
            TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(ccarry[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        psum[i] =
            TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(csum[i], sk->key.lvl1);
    // for (int i = 0; i < num_test; i++)
    // std::cout<<static_cast<int>(pcarry[i])<<":"<<static_cast<int>(psum[i])<<";"<<static_cast<int>(((pa[i]
    // + pb[i] + pc[i]) & 2) >> 1)<<":"<<static_cast<int>((pa[i] + pb[i] +
    // pc[i]) & 1)<<std::endl;
    for (int i = 0; i < num_test; i++)
        assert(pcarry[i] == (((pa[i] + pb[i] + pc[i]) & 2) >> 1));
    for (int i = 0; i < num_test; i++)
        assert(psum[i] == ((pa[i] + pb[i] + pc[i]) & 1));
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}
