#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());

    using high2midP = TFHEpp::lvl31param;
    using mid2lowP = TFHEpp::lvl1hparam;
    using low2midP = TFHEpp::lvlh1param;
    constexpr auto numtest = 10;
    constexpr uint basebit = 4;
    constexpr uint64_t numdigits = 64 / 4;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    std::uniform_int_distribution<int32_t> lvl3gen(
        TFHEpp::lvl3param::key_value_min, TFHEpp::lvl3param::key_value_max);
    for (typename TFHEpp::lvl3param::T &i : sk.key.lvl3)
        i = lvl3gen(TFHEpp::generator);
    ek.emplacebkfft<low2midP>(sk);
    ek.emplaceiksk<mid2lowP>(sk);
    ek.emplaceiksk<high2midP>(sk);

    // Test input
    std::uniform_int_distribution<typename TFHEpp::lvl3param::T> messagegen(
        0, 2 * TFHEpp::lvl3param::plain_modulus - 1);
    std::uniform_int_distribution<typename TFHEpp::lvl3param::T> maskgen(
        0, std::numeric_limits<typename TFHEpp::lvl3param::T>::max());
    std::array<typename TFHEpp::lvl3param::T, numtest> plains{};
    for (typename TFHEpp::lvl3param::T &i : plains) {
        i = messagegen(engine);
    }
    std::array<TFHEpp::TLWE<TFHEpp::lvl3param>, numtest> ciphers{};
    for (uint i = 0; i < numtest; i++) {
        ciphers[i] = TFHEpp::tlweSymIntEncrypt<TFHEpp::lvl3param>(
            plains[i], TFHEpp::lvl3param::α, sk.key.lvl3);
        ciphers[i][TFHEpp::lvl3param::n] += maskgen(engine);
    }

    // Test output
    std::array<std::array<TFHEpp::TLWE<typename high2midP::targetP>, numdigits>,
               numtest>
        result_multiple{};

    // Convert TLWE
    for (uint i = 0; i < numtest; i++) {
        // converter.toLv1TLWE(ciphers.at(i), result_multiple.at(i));
        TFHEpp::HomDecomp<high2midP, mid2lowP, low2midP, basebit, numdigits>(
            result_multiple.at(i), ciphers.at(i), ek.getiksk<high2midP>(),
            ek.getiksk<mid2lowP>(), ek.getbkfft<low2midP>());
    }

    // Check the correctness of the results
    for (uint test = 0; test < numtest; test++) {
        uint64_t phase = TFHEpp::tlweSymPhase<typename high2midP::domainP>(
            ciphers.at(test), sk.key.lvl3);
        for (uint digit = 0; digit < numdigits; digit++) {
            int plainResult =
                TFHEpp::tlweSymIntDecrypt<typename high2midP::targetP,
                                          1U << basebit>(
                    result_multiple[test][digit],
                    sk.key.get<typename high2midP::targetP>());
            const uint64_t plainExpected =
                ((phase >> (basebit * digit)) & ((1ULL << basebit) - 1));
            plainResult = (plainResult + (1U << basebit)) % (1U << basebit);
            assert(plainExpected == plainResult);
        }
    }
    std::cout << "PASS" << std::endl;
}