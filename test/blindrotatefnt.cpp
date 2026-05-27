#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>

#include <gatebootstrapping.hpp>
#include <key.hpp>
#include <tlwe.hpp>
#include <trlwe.hpp>

int main()
{
    constexpr uint32_t num_test = 2;
    using bkP = TFHEpp::lvl01param;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey sk;
    auto bkfnt =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFNT<bkP>>();

    auto start = std::chrono::system_clock::now();
    TFHEpp::bkfntgen<bkP>(*bkfnt, sk);
    auto keygen_end = std::chrono::system_clock::now();

    const TFHEpp::Polynomial<typename bkP::targetP> testvector =
        TFHEpp::μpolygen<typename bkP::targetP, bkP::targetP::μ>();

    std::array<TFHEpp::TLWE<typename bkP::domainP>, num_test> tlwe;
    std::array<TFHEpp::TLWE<typename bkP::targetP>, num_test> bootedtlwe;
    std::array<bool, num_test> p;

    for (int i = 0; i < num_test; i++) p[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++)
        TFHEpp::tlweSymEncrypt<typename bkP::domainP>(
            tlwe[i], p[i] ? bkP::domainP::μ : -bkP::domainP::μ,
            bkP::domainP::α, sk.key.get<typename bkP::domainP>());

    auto blindrotate_start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        alignas(64) TFHEpp::TRLWE<typename bkP::targetP> acc;
        TFHEpp::BlindRotate<bkP>(acc, tlwe[test], *bkfnt, testvector);
        TFHEpp::SampleExtractIndex<typename bkP::targetP>(bootedtlwe[test],
                                                          acc, 0);
    }
    auto end = std::chrono::system_clock::now();

    for (int i = 0; i < num_test; i++) {
        bool p2 = TFHEpp::tlweSymDecrypt<typename bkP::targetP>(
            bootedtlwe[i], sk.key.get<typename bkP::targetP>());
        if (p[i] != p2) {
            std::cerr << "FNT BlindRotate decrypted " << p2 << ", expected "
                      << p[i] << std::endl;
            return 1;
        }
    }

    const double keygen_elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(keygen_end -
                                                              start)
            .count();
    const double blindrotate_elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                              blindrotate_start)
            .count();

    std::cout << "Passed" << std::endl;
    std::cout << "keygen: " << keygen_elapsed << "ms" << std::endl;
    std::cout << "blindrotate: " << blindrotate_elapsed / num_test << "ms"
              << std::endl;
    return 0;
}
