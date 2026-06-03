#include <array>
#include <cstdint>
#include <iostream>
#include <random>

#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 10;
    using P = TFHEpp::lvl1param;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::lweKey key;

    for (const int32_t message : {1, -1}) {
        TFHEpp::Polynomial<P> plain = {};
        plain[0] = static_cast<typename P::T>(message);

        TFHEpp::TRGSWFNT<P> trgswfnt;
        TFHEpp::trgswSymEncrypt<P>(trgswfnt, plain, key.get<P>());

        for (uint32_t test = 0; test < num_test; test++) {
            std::array<bool, P::n> p;
            TFHEpp::Polynomial<P> pmu;
            for (uint32_t i = 0; i < P::n; i++) {
                p[i] = binary(engine) > 0;
                pmu[i] = p[i] ? P::μ : -P::μ;
            }

            TFHEpp::TRLWE<P> c;
            TFHEpp::trlweSymEncrypt<P>(c, pmu, key.get<P>());
            TFHEpp::ExternalProduct<P>(c, c, trgswfnt);

            const std::array<bool, P::n> decrypted =
                TFHEpp::trlweSymDecrypt<P>(c, key.get<P>());
            for (uint32_t i = 0; i < P::n; i++) {
                const bool expected = message > 0 ? p[i] : !p[i];
                if (decrypted[i] != expected) {
                    std::cerr << "FNT ExternalProduct decrypted "
                              << decrypted[i] << ", expected " << expected
                              << " at test " << test << ", coefficient " << i
                              << std::endl;
                    return 1;
                }
            }
        }
    }

    std::cout << "Passed" << std::endl;
    return 0;
}
