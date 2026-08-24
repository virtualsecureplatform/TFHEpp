#include <chrono>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>

#include <tfhe++.hpp>

namespace {

template <class Function>
double AverageMilliseconds(const std::uint32_t repetitions, Function &&function)
{
    const auto begin = std::chrono::steady_clock::now();
    for (std::uint32_t iteration = 0; iteration < repetitions; iteration++)
        function();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - begin).count() /
           repetitions;
}

}  // namespace

int main()
{
#ifndef USE_BLOCK_BINARY
    std::cout << "Skipped: requires USE_BLOCK_BINARY" << std::endl;
    return 0;
#else
    using ShallowP = TFHEpp::shallowboot::structured_binary_std128_gateparam;
    using ShallowDomainP = typename ShallowP::domainP;
    using RingP = typename ShallowP::targetP;
    using BlockP = TFHEpp::lvl01param;
    using BlockDomainP = typename BlockP::domainP;

    constexpr std::uint32_t repetitions = 100;
    constexpr auto blocks = [] {
        std::array<std::uint32_t, ShallowDomainP::sparse_hamming_weight + 1>
            result{};
        for (std::uint32_t i = 0;
             i <= ShallowDomainP::sparse_hamming_weight; i++)
            result[i] =
                i * (ShallowDomainP::n / ShallowDomainP::sparse_hamming_weight);
        return result;
    }();

    std::mt19937 engine(20261730);
    TFHEpp::Key<ShallowDomainP> shallow_source_key;
    TFHEpp::structuredSparseBinaryKeyGen<ShallowDomainP>(shallow_source_key,
                                                         blocks, engine);
    TFHEpp::Key<BlockDomainP> block_source_key;
    TFHEpp::keyGen<BlockDomainP>(block_source_key);
    TFHEpp::Key<RingP> accumulator_key;
    TFHEpp::keyGen<RingP>(accumulator_key);

    auto shallow_bk = std::make_unique<TFHEpp::BootstrappingKeyFFT<ShallowP>>();
    TFHEpp::bkfftgen<ShallowP>(*shallow_bk, shallow_source_key,
                                accumulator_key);
    auto block_bk = std::make_unique<TFHEpp::BootstrappingKeyFFT<BlockP>>();
    TFHEpp::bkfftgen<BlockP>(*block_bk, block_source_key, accumulator_key);

    TFHEpp::TLWE<ShallowDomainP> shallow_input;
    TFHEpp::tlweSymEncryptModQ<ShallowDomainP>(shallow_input, true,
                                                shallow_source_key);
    TFHEpp::TLWE<BlockDomainP> block_input;
    TFHEpp::tlweSymEncrypt<BlockDomainP>(block_input, BlockDomainP::μ,
                                         BlockDomainP::α, block_source_key);
    const auto shallow_lut =
        TFHEpp::ShallowBootBinaryIdentityTestVector<ShallowP, RingP::μ>();
    const auto block_lut = TFHEpp::μpolygen<RingP, RingP::μ>();

    TFHEpp::TLWE<RingP> shallow_output;
    TFHEpp::TLWE<RingP> block_output;

    // Warm both paths and verify that timing excludes key generation.
    TFHEpp::StructuredSparseGateBootstrappingTLWE2TLWE<ShallowP>(
        shallow_output, shallow_input, *shallow_bk, shallow_lut, blocks);
    TFHEpp::GateBootstrappingTLWE2TLWE<BlockP>(block_output, block_input,
                                               *block_bk, block_lut);
    if (!TFHEpp::tlweSymDecrypt<RingP>(shallow_output, accumulator_key) ||
        !TFHEpp::tlweSymDecrypt<RingP>(block_output, accumulator_key)) {
        std::cerr << "benchmark setup produced an incorrect bootstrap"
                  << std::endl;
        return 1;
    }

    const double shallow_core_ms = AverageMilliseconds(repetitions, [&] {
        TFHEpp::StructuredSparseGateBootstrappingTLWE2TLWE<ShallowP>(
            shallow_output, shallow_input, *shallow_bk, shallow_lut, blocks);
    });
    const double block_core_ms = AverageMilliseconds(repetitions, [&] {
        TFHEpp::GateBootstrappingTLWE2TLWE<BlockP>(
            block_output, block_input, *block_bk, block_lut);
    });

    std::cout << "repetitions=" << repetitions << '\n'
              << "shallow_external_products="
              << ShallowDomainP::sparse_hamming_weight << '\n'
              << "blockbinary_external_products="
              << BlockDomainP::n / BlockDomainP::ell << '\n'
              << "shallow_core_ms=" << shallow_core_ms << '\n'
              << "blockbinary_core_ms=" << block_core_ms << '\n'
              << "core_speedup=" << block_core_ms / shallow_core_ms << '\n';
#endif
}
