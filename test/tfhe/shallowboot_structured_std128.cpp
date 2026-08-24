#include <array>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>

#include <tfhe++.hpp>

int main()
{
    using bkP = TFHEpp::shallowboot::structured_binary_std128_gateparam;
    using DomainP = typename bkP::domainP;
    using TargetP = typename bkP::targetP;

    static_assert(DomainP::n == 1024);
    static_assert(DomainP::sparse_hamming_weight == 64);
    static_assert(bkP::blind_rotation_modulus == 512);

    constexpr auto blocks = [] {
        std::array<std::uint32_t, DomainP::sparse_hamming_weight + 1> result{};
        for (std::uint32_t i = 0; i <= DomainP::sparse_hamming_weight; i++)
            result[i] = i * (DomainP::n / DomainP::sparse_hamming_weight);
        return result;
    }();

    std::mt19937 engine(20261730);
    TFHEpp::Key<DomainP> source_key;
    TFHEpp::structuredSparseBinaryKeyGen<DomainP>(source_key, blocks, engine);
    TFHEpp::Key<TargetP> target_key;
    TFHEpp::keyGen<TargetP>(target_key);

    auto bootstrap_key =
        std::make_unique<TFHEpp::BootstrappingKeyFFT<bkP>>();
    TFHEpp::bkfftgen<bkP>(*bootstrap_key, source_key, target_key);

    std::bernoulli_distribution bit;
    for (int trial = 0; trial < 32; trial++) {
        const bool plain = bit(engine);
        TFHEpp::TLWE<DomainP> input;
        TFHEpp::tlweSymEncryptModQ<DomainP>(input, plain, source_key);
        TFHEpp::TLWE<TargetP> output;
        TFHEpp::StructuredSparseGateBootstrappingTLWE2TLWE<bkP>(
            output, input, *bootstrap_key,
            TFHEpp::ShallowBootBinaryIdentityTestVector<bkP, TargetP::μ>(),
            blocks);
        if (TFHEpp::tlweSymDecrypt<TargetP>(output, target_key) != plain) {
            std::cerr << "shallow structured bootstrap decrypted the wrong bit\n";
            return 1;
        }
    }

    std::cout << "Passed" << std::endl;
}
