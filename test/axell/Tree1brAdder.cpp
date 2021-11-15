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
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey* sk = new TFHEpp::SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    std::vector<std::array<uint8_t,6>> pin(num_test);
    std::vector<std::array<uint8_t,2>> pres(num_test);
    for (int i = 0; i < num_test; i++)for(int j = 0; j<6; j++) pin[i][j] = binary(engine);
    std::vector<std::array<TFHEpp::TLWE<TFHEpp::lvl1param>,6>> cin(num_test);
    std::vector<std::array<TFHEpp::TLWE<TFHEpp::lvl1param>,2>> cres(num_test);

    for (int i = 0; i < num_test; i++)
        for(int j = 0; j<6; j++)
            cin[i][j] = TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
                pin[i][j] ? TFHEpp::lvlMparam::μ : -TFHEpp::lvlMparam::μ,
                TFHEpp::lvl1param::α, sk->key.lvl1);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TFHEpp::TLWE<TFHEpp::lvl1param> ccarry,csum;
        TFHEpp::HomFullAdder(cres[test][1], cres[test][0], cin[test][0], cin[test][1],
                                cin[test][2], ek);
        TFHEpp::HomFullAdder(ccarry, csum, cin[test][3], cin[test][4],
                                cin[test][5], ek);
        TFHEpp::HomFullAdder(cres[test][1], cres[test][0],  cres[test][0], ccarry,
                                csum, ek);
    }

    end = std::chrono::system_clock::now();
    for (int i = 0; i < num_test; i++)
        for(int j = 0; j<2;j++)
        pres[i][j] +=
            static_cast<uint8_t>(TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(cres[i][j], sk->key.lvl1));
    // for (int i = 0; i < num_test; i++)
    // std::cout<<static_cast<int>(pcarry[i])<<":"<<static_cast<int>(psum[i])<<";"<<static_cast<int>(((pa[i]
    // + pb[i] + pc[i]) & 2) >> 1)<<":"<<static_cast<int>((pa[i] + pb[i] +
    // pc[i]) & 1)<<std::endl;
    // for (int i = 0; i < num_test; i++)
    //     std::cout<<static_cast<int>(pres[i])<<':'<<pa[i] + pb[i]<<std::endl;
    for (int i = 0; i < num_test; i++){
        const int a = (pin[i][0] + pin[i][1] + pin[i][2]) & 1;
        const int b = (pin[i][3] + pin[i][4] + pin[i][5]) & 1;
        const int c = ((pin[i][3] + pin[i][4] + pin[i][5])>>1) & 1;
        assert(pres[i][0] == ((a + b + c)&1));
        assert(pres[i][1] == (((a + b + c)>>1)&1));
    }

    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}
