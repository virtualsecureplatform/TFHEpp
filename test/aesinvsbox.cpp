#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using brP = TFHEpp::lvlh2param;
    using iksP = TFHEpp::lvl2hparam;
    using ahP = TFHEpp::AHlvl2param;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    constexpr uint32_t plain_modulus = 1 << (4 + 1);
    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    constexpr uint num_test = 1U << 8;
    // constexpr uint num_test = 1U<<4;
    std::vector<std::array<TFHEpp::TLWE<typename iksP::domainP>, 2>> cin(
        num_test);

    for (int i = 0; i < num_test; i++) {
        cin[i][0] =
            TFHEpp::tlweSymIntEncrypt<typename iksP::domainP, plain_modulus>(
                i & 0xF, *sk);
        cin[i][1] =
            TFHEpp::tlweSymIntEncrypt<typename iksP::domainP, plain_modulus>(
                (i >> 4), *sk);
    }
    std::vector<std::array<TFHEpp::TLWE<typename brP::targetP>, 2>> cres(
        num_test);
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(*sk);
    ek.emplaceiksk<iksP>(*sk);
    ek.emplaceahk<ahP>(*sk);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        std::cout << "test: " << test << std::endl;
        TFHEpp::AESInvSbox<iksP, brP, ahP>(cres[test], cin[test], ek);
    }

    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
    for (int i = 0; i < num_test; i++) {
        const uint8_t pres =
            (TFHEpp::tlweSymIntDecrypt<typename brP::targetP, plain_modulus>(
                 cres[i][1], *sk)
             << 4) +
            TFHEpp::tlweSymIntDecrypt<typename brP::targetP, plain_modulus>(
                cres[i][0], *sk);
        // std::cout << "test: " << i << " pres: " << (int)pres << " expected: "
        // << (int)inv_sbox[i>>4][i&0xF] << std::endl;
        assert(pres == inv_sbox[i >> 4][i & 0xF]);
    }
    std::cout << "Passed" << std::endl;
}