#include <bits/stdint-uintn.h>
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
    std::uniform_int_distribution<uint32_t> indist(0, 3);

    TFHEpp::SecretKey* sk = new TFHEpp::SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    std::vector<uint8_t> pa(num_test);
    std::vector<uint8_t> pb(num_test);
    std::vector<uint8_t> pres(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = indist(engine);
    for (int i = 0; i < num_test; i++) pb[i] = indist(engine);
    for (int i = 0; i < num_test; i++) pres[i] = 0;
    std::vector<std::array<TFHEpp::TLWE<TFHEpp::lvl1param>,2>> ca(num_test);
    std::vector<std::array<TFHEpp::TLWE<TFHEpp::lvl1param>,2>> cb(num_test);
    std::vector<std::array<TFHEpp::TLWE<TFHEpp::lvl1param>,3>> cres(num_test);

    for (int i = 0; i < num_test; i++)
        for(int j = 0; j<2; j++){
            ca[i][j] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
                (pa[i]>>j)&1 ? TFHEpp::lvlMparam::μ : -TFHEpp::lvlMparam::μ,
                TFHEpp::lvl1param::α, sk->key.lvl1);
            TFHEpp::GateBootstrapping<TFHEpp::lvlM0param,TFHEpp::lvl0Mparam,TFHEpp::lvlMparam::μ>(ca[i][j],ca[i][j],ek);
        }
    for (int i = 0; i < num_test; i++)
        for(int j = 0; j<2; j++){
            cb[i][j] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
                (pb[i]>>j)&1 ? TFHEpp::lvlMparam::μ : -TFHEpp::lvlMparam::μ,
                TFHEpp::lvl1param::α, sk->key.lvl1);
            TFHEpp::GateBootstrapping<TFHEpp::lvlM0param,TFHEpp::lvl0Mparam,TFHEpp::lvlMparam::μ>(cb[i][j],cb[i][j],ek);
        }

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::TLWE<TFHEpp::lvl1param> ccarry;
        TFHEpp::HomHalfAdder<TFHEpp::lvlMparam>(ccarry, cres[test][0], ca[test][0], cb[test][0], ek);
        TFHEpp::HomFullAdder(cres[test][2], cres[test][1], ca[test][1], cb[test][1],
                                ccarry, ek);
    }

    end = std::chrono::system_clock::now();
    for (int i = 0; i < num_test; i++)
        for(int j = 0; j<3;j++)
        pres[i] +=
            static_cast<uint8_t>(TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(cres[i][j], sk->key.lvl1))<<j;
    // for (int i = 0; i < num_test; i++)
    // std::cout<<static_cast<int>(pcarry[i])<<":"<<static_cast<int>(psum[i])<<";"<<static_cast<int>(((pa[i]
    // + pb[i] + pc[i]) & 2) >> 1)<<":"<<static_cast<int>((pa[i] + pb[i] +
    // pc[i]) & 1)<<std::endl;
    for (int i = 0; i < num_test; i++)
        std::cout<<static_cast<int>(pres[i])<<':'<<pa[i] + pb[i]<<std::endl;
    for (int i = 0; i < num_test; i++)
        assert(pres[i] == (pa[i] + pb[i]) );

    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}
