#include <array>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>

#include <tfhe++.hpp>

int main()
{
    using bkP = TFHEpp::lvl01param;
    using DomainP = typename bkP::domainP;
    using TargetP = typename bkP::targetP;

    // Nine equal buckets give a 9-sparse secret instead of the usual 630
    // sequential CMUXs.  The production parameters and the regular TFHE FFT
    // bootstrapping key are intentionally unchanged.
    constexpr std::array<std::uint32_t, 10> blocks = {
        0, 70, 140, 210, 280, 350, 420, 490, 560, 630};

    std::mt19937 engine(17);
    TFHEpp::Key<DomainP> sparse_key;
    TFHEpp::structuredSparseBinaryKeyGen<DomainP>(sparse_key, blocks, engine);
    for (std::size_t block = 0; block + 1 < blocks.size(); block++) {
        std::uint32_t weight = 0;
        for (std::uint32_t i = blocks[block]; i < blocks[block + 1]; i++)
            weight += sparse_key[i];
        if (weight != 1) {
            std::cerr << "sparse key is not one-hot in block " << block
                      << '\n';
            return 1;
        }
    }
    TFHEpp::Key<TargetP> target_key;
    TFHEpp::keyGen<TargetP>(target_key);

    auto bk = std::make_unique<TFHEpp::BootstrappingKeyFFT<bkP>>();
    TFHEpp::bkfftgen<bkP>(*bk, sparse_key, target_key);

    std::bernoulli_distribution bit;
    for (int test = 0; test < 4; test++) {
        const bool plain = bit(engine);
        TFHEpp::TLWE<DomainP> input;
        TFHEpp::tlweSymEncrypt<DomainP>(
            input, plain ? DomainP::μ : -DomainP::μ, DomainP::α, sparse_key);

        TFHEpp::TLWE<TargetP> output;
        TFHEpp::StructuredSparseGateBootstrappingTLWE2TLWE<bkP>(
            output, input, *bk,
            TFHEpp::μpolygen<TargetP, TargetP::μ>(), blocks);
        if (TFHEpp::tlweSymDecrypt<TargetP>(output, target_key) != plain) {
            std::cerr << "structured sparse bootstrapping decrypted the wrong "
                         "bit\n";
            return 1;
        }
    }

    std::cout << "Passed" << std::endl;
}
