#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using brP = TFHEpp::lvl02param;
    using iksP = TFHEpp::lvl20param;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    constexpr uint num_test = 1U << 8;
    // constexpr uint num_test = 1U<<4;
    std::vector<std::array<TFHEpp::TLWE<typename iksP::domainP>, 8>> cin(
        num_test);

    for (int i = 0; i < num_test; i++) {
        for (int j = 0; j < 8; j++)
            cin[i][j] = TFHEpp::tlweSymEncrypt<typename iksP::domainP>(
                ((i >> j) & 0x1)
                    ? 1ULL << (std::numeric_limits<
                                   typename iksP::domainP::T>::digits -
                               2)
                    : -(1ULL << (std::numeric_limits<
                                     typename iksP::domainP::T>::digits -
                                 2)),
                *sk);
    }
    std::vector<std::array<TFHEpp::TLWE<typename brP::targetP>, 8>> cres(
        num_test);
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(*sk);
    ek.emplaceiksk<iksP>(*sk);
    ek.emplaceahk<typename brP::targetP>(*sk);
    ek.emplacecbsk<typename brP::targetP>(*sk);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        std::cout << "test: " << test << std::endl;
        TFHEpp::AESSboxROM<iksP, brP>(cres[test], cin[test], ek);
    }

    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
    for (int i = 0; i < num_test; i++) {
        uint8_t pres = 0;
        for (int j = 0; j < 8; j++) {
            if (TFHEpp::tlweSymDecrypt<typename brP::targetP>(cres[i][j], *sk))
                pres += (1 << j);
        }
        // std::cout << "test: " << i << " pres: " << (int)pres << " expected: "
        // << (int)inv_sbox[i>>4][i&0xF] << std::endl;
        assert(pres == sbox[i >> 4][i & 0xF]);
    }
    std::cout << "Passed" << std::endl;
}