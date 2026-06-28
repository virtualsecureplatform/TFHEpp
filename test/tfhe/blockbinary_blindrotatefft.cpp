#include <array>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>

#include <tfhe/gatebootstrapping.hpp>
#include <tfhe/key.hpp>
#include <tfhe/tlwe.hpp>
#include <tfhe/trlwe.hpp>

int main()
{
#ifndef USE_BLOCK_BINARY
    std::cout << "Skipped" << std::endl;
    return 0;
#else
    constexpr uint32_t num_test = 4;
    using bkP = TFHEpp::lvl01param;

    static_assert(bkP::domainP::n == 630);
    static_assert(bkP::domainP::ell == 2);
    static_assert(bkP::domainP::key_value_min == 0);
    static_assert(bkP::domainP::key_value_max == 1);

    std::mt19937 engine(0);
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey sk;
    auto bkfft =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFFT<bkP>>();
    TFHEpp::bkfftgen<bkP>(*bkfft, sk);

    const TFHEpp::Polynomial<typename bkP::targetP> testvector =
        TFHEpp::μpolygen<typename bkP::targetP, bkP::targetP::μ>();

    std::array<TFHEpp::TLWE<typename bkP::domainP>, num_test> tlwe;
    std::array<TFHEpp::TLWE<typename bkP::targetP>, num_test> bootedtlwe;
    std::array<bool, num_test> plain;

    for (int i = 0; i < num_test; i++) plain[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) {
        TFHEpp::tlweSymEncrypt<typename bkP::domainP>(
            tlwe[i], plain[i] ? bkP::domainP::μ : -bkP::domainP::μ,
            bkP::domainP::α, sk.key.get<typename bkP::domainP>());
    }

    for (int test = 0; test < num_test; test++) {
        alignas(64) TFHEpp::TRLWE<typename bkP::targetP> acc;
        TFHEpp::BlindRotate<bkP>(acc, tlwe[test], *bkfft, testvector);
        TFHEpp::SampleExtractIndex<typename bkP::targetP>(bootedtlwe[test],
                                                          acc, 0);
    }

    for (int i = 0; i < num_test; i++) {
        const bool decrypted = TFHEpp::tlweSymDecrypt<typename bkP::targetP>(
            bootedtlwe[i], sk.key.get<typename bkP::targetP>());
        if (plain[i] != decrypted) {
            std::cerr << "FFT BlindRotate decrypted " << decrypted
                      << ", expected " << plain[i] << std::endl;
            return 1;
        }
    }

    std::cout << "Passed" << std::endl;
    return 0;
#endif
}
