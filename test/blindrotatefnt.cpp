#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>

#include <gatebootstrapping.hpp>
#include <key.hpp>
#include <trgsw.hpp>
#include <tlwe.hpp>
#include <trlwe.hpp>

int main()
{
    constexpr uint32_t num_test = 2;
    constexpr uint32_t num_active = 3;
    using bkP = TFHEpp::lvl01param;

    TFHEpp::SecretKey sk;
    const auto domainkey = sk.key.get<typename bkP::domainP>();
    const auto targetkey = sk.key.get<typename bkP::targetP>();

    std::array<uint32_t, num_active> active_index = {};
    uint32_t active_count = 0;
    for (uint32_t i = 0;
         i < bkP::domainP::k * bkP::domainP::n && active_count < num_active;
         i++)
        if (domainkey[i] != 0) active_index[active_count++] = i;
    if (active_count != num_active) {
        std::cerr << "Failed to find active secret-key coefficients"
                  << std::endl;
        return 1;
    }

    auto bkfnt =
        std::make_unique_for_overwrite<TFHEpp::BootstrappingKeyFNT<bkP>>();

    auto start = std::chrono::system_clock::now();
    for (uint32_t i = 0; i < num_active; i++) {
        TFHEpp::Polynomial<typename bkP::targetP> plainpoly = {};
        plainpoly[0] = domainkey[active_index[i]];
        TFHEpp::trgswSymEncrypt<typename bkP::targetP>(
            (*bkfnt)[active_index[i]], plainpoly, targetkey);
    }
    auto keygen_end = std::chrono::system_clock::now();

    const TFHEpp::Polynomial<typename bkP::targetP> testvector =
        TFHEpp::μpolygen<typename bkP::targetP, bkP::targetP::μ>();

    std::array<TFHEpp::TLWE<typename bkP::domainP>, num_test> tlwe;
    std::array<TFHEpp::TLWE<typename bkP::targetP>, num_test> bootedtlwe;
    const std::array<bool, num_test> p = {false, true};
    const std::array<uint32_t, num_active> rotation = {1, 7, 23};
    constexpr uint32_t modswitch_shift =
        std::numeric_limits<typename bkP::domainP::T>::digits - 1 -
        bkP::targetP::nbit;

    for (int test = 0; test < num_test; test++) {
        tlwe[test] = {};
        typename bkP::domainP::T b =
            p[test] ? bkP::domainP::μ : -bkP::domainP::μ;
        for (uint32_t i = 0; i < num_active; i++) {
            const typename bkP::domainP::T a = rotation[i]
                                               << modswitch_shift;
            tlwe[test][active_index[i]] = a;
            b += a * domainkey[active_index[i]];
        }
        tlwe[test][bkP::domainP::k * bkP::domainP::n] = b;
    }

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
        std::chrono::duration_cast<std::chrono::milliseconds>(end - keygen_end)
            .count();

    std::cout << "Passed" << std::endl;
    std::cout << "selected keygen: " << keygen_elapsed << "ms" << std::endl;
    std::cout << "blindrotate: " << blindrotate_elapsed / num_test << "ms"
              << std::endl;
    return 0;
}
