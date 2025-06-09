#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl2param;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    constexpr uint num_test = 1000;
    std::vector<std::array<TFHEpp::TLWE<P>, 32>> cin(num_test);
    std::vector<std::array<uint8_t, 32>> pin(num_test);

    for (int i = 0; i < num_test; i++) {
        for (int j = 0; j < 32; j++) {
            pin[i][j] = binary(engine);
            cin[i][j] = TFHEpp::tlweSymEncrypt<P>(
                pin[i][j]
                    ? 1ULL << (std::numeric_limits<typename P::T>::digits - 2)
                    : -(1ULL
                        << (std::numeric_limits<typename P::T>::digits - 2)),
                *sk);
            cin[i][j][P::k * P::n] +=
                1ULL << (std::numeric_limits<typename P::T>::digits - 2);
        }
    }

    std::vector<std::array<TFHEpp::TLWE<P>, 32>> cres(num_test);
    std::vector<std::array<TFHEpp::TLWE<P>, 32>> cref(num_test);

    for (int test = 0; test < num_test; test++) {
        // std::cout << "test: " << test << std::endl;
        TFHEpp::MixColumn<P>(cref[test], cin[test]);
        TFHEpp::MixColumnDepth4<P>(cres[test], cin[test]);
        for (int j = 0; j < 32; j++) {
            cref[test][j][P::k * P::n] -=
                1ULL << (std::numeric_limits<typename P::T>::digits - 2);
            cres[test][j][P::k * P::n] -=
                1ULL << (std::numeric_limits<typename P::T>::digits - 2);
        }
    }

    for (int i = 0; i < num_test; i++) {
        uint32_t pincat = 0;
        uint32_t pres = 0;
        uint32_t pref = 0;
        for (int j = 0; j < 32; j++) {
            pres |= (static_cast<uint32_t>(
                        TFHEpp::tlweSymDecrypt<P>(cres[i][j], *sk)))
                    << j;
            pref |= (static_cast<uint32_t>(
                        TFHEpp::tlweSymDecrypt<P>(cref[i][j], *sk)))
                    << j;
            pincat |= static_cast<uint32_t>(pin[i][j]) << j;
        }
        // std::bitset<32> bpin(pincat);
        // std::bitset<32> bpres(pres);
        // std::bitset<32> bpref(pref);
        // std::cout << "i: " << i << std::endl;
        // std::cout << "pin: " << bpin << std::endl;
        // std::cout << "pres: " << bpres << std::endl;
        // std::cout << "pref: " << bpref << std::endl;
        assert(pres == pref);
    }
    std::cout << "PASS" << std::endl;
}