#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <tfhe++.hpp>

int main()
{
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl21param;
    using targetP = typename privksP::targetP;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<bkP>(sk);
    ek.emplacesubiksk<privksP>(sk);
    ek.emplacesubprivksk4cb<privksP>(sk);

    std::array<targetP::T, targetP::n> message;
    for (std::size_t i = 0; i < message.size(); ++i)
        message[i] = (i & 1) ? targetP::μ : -targetP::μ;

    TFHEpp::TRLWE<targetP> encrypted;
    TFHEpp::trlweSymEncrypt<targetP>(encrypted, message,
                                     sk.key.getSubset<targetP>());

    for (const bool selector : {false, true}) {
        TFHEpp::TLWE<typename bkP::domainP> encrypted_selector;
        TFHEpp::tlweSymEncrypt<typename bkP::domainP>(
            encrypted_selector,
            selector ? bkP::domainP::μ : -bkP::domainP::μ,
            sk.key.getSubset<typename bkP::domainP>());

        TFHEpp::TRGSWFFT<targetP> normal, inverted;
        TFHEpp::CircuitBootstrappingSubsetWithInv<bkP, privksP>(
            normal, inverted, encrypted_selector, ek);

        for (const auto* control : {&normal, &inverted}) {
            auto product = encrypted;
            TFHEpp::ExternalProduct<targetP>(product, product, *control);
            const bool should_preserve = control == &normal ? selector
                                                             : !selector;
            if (should_preserve) {
                const auto decrypted = TFHEpp::trlweSymDecrypt<targetP>(
                    product, sk.key.getSubset<targetP>());
                for (std::size_t i = 0; i < decrypted.size(); ++i)
                    assert(decrypted[i] == static_cast<bool>(i & 1));
            }
            else {
                const auto phase = TFHEpp::trlwePhase<targetP>(
                    product, sk.key.getSubset<targetP>());
                for (const auto coefficient : phase) {
                    const auto centered = static_cast<std::int32_t>(coefficient);
                    assert(centered > -static_cast<std::int32_t>(targetP::μ));
                    assert(centered < static_cast<std::int32_t>(targetP::μ));
                }
            }
        }
    }

    std::cout << "Passed" << std::endl;
}
