#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <tfhe/gate.hpp>
#include <tfhe/key.hpp>
#include <tfhe/tlwe.hpp>
#include <type_traits>

int main()
{
#ifndef USE_BLOCK_BINARY
    std::cout << "Skipped" << std::endl;
    return 0;
#else
    constexpr uint32_t num_test = 8;
    using bkP = TFHEpp::lvl01param;
    using iksP = TFHEpp::lvl10param;

    static_assert(bkP::domainP::n == 630);
    static_assert(bkP::domainP::ell == 2);
    static_assert(
        std::is_same_v<typename iksP::domainP, typename bkP::targetP>);
    static_assert(
        std::is_same_v<typename iksP::targetP, typename bkP::domainP>);

    std::mt19937 engine(0);
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<bkP>(sk);
    ek.emplacesubiksk<iksP>(sk);

    const auto lvl0key = sk.key.getSubset<TFHEpp::lvl0param>();
    const auto lvl1key = sk.key.getSubset<TFHEpp::lvl1param>();
    for (int i = 0; i < TFHEpp::lvl0param::k * TFHEpp::lvl0param::n; i++) {
        if (static_cast<typename TFHEpp::lvl1param::T>(lvl0key[i]) !=
            lvl1key[i]) {
            std::cerr << "lvl0 key is not a prefix of lvl1 key" << std::endl;
            return 1;
        }
    }

    std::array<TFHEpp::TLWE<TFHEpp::lvl0param>, num_test> ca;
    std::array<TFHEpp::TLWE<TFHEpp::lvl0param>, num_test> cb;
    std::array<TFHEpp::TLWE<TFHEpp::lvl0param>, num_test> cres;
    std::array<bool, num_test> pa;
    std::array<bool, num_test> pb;

    for (int i = 0; i < num_test; i++) {
        pa[i] = binary(engine) > 0;
        pb[i] = binary(engine) > 0;
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
            ca[i], pa[i] ? TFHEpp::lvl0param::μ : -TFHEpp::lvl0param::μ, sk);
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
            cb[i], pb[i] ? TFHEpp::lvl0param::μ : -TFHEpp::lvl0param::μ, sk);
    }

    for (int i = 0; i < num_test; i++)
        TFHEpp::HomNAND<bkP, TFHEpp::lvl1param::μ, iksP>(cres[i], ca[i], cb[i],
                                                         ek);

    for (int i = 0; i < num_test; i++) {
        const bool decrypted =
            TFHEpp::tlweSymDecrypt<TFHEpp::lvl0param>(cres[i], sk);
        const bool expected = !(pa[i] & pb[i]);
        if (decrypted != expected) {
            std::cerr << "NAND decrypted " << decrypted << ", expected "
                      << expected << std::endl;
            return 1;
        }
    }

    std::cout << "Passed" << std::endl;
    return 0;
#endif
}
