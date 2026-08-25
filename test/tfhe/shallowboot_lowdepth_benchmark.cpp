#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include <tfhe/shallowboot_lowdepth.hpp>

template <class Function>
double AverageMilliseconds(const std::uint32_t repetitions, Function &&function)
{
    const auto begin = std::chrono::steady_clock::now();
    for (std::uint32_t repetition = 0; repetition < repetitions; repetition++)
        function();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - begin).count() /
           repetitions;
}

int main()
{
    using namespace TFHEpp::shallowboot::lowdepth;
    constexpr std::size_t h = 29;
    constexpr std::uint32_t q = 512;
    constexpr std::uint32_t repetitions = 20;
    const std::array<TFHEpp::modular_ntt::PrimeModulus, 2> primes = {
        TFHEpp::modular_ntt::wide_primes[0],
        TFHEpp::modular_ntt::wide_primes[2]};
    RNSRing high_ring(4096, primes);
    std::mt19937_64 engine(20261730);
    const RNSSecret binary_secret = BinaryNTTSecretGen(high_ring, engine);

    std::array<std::uint8_t, h> lwe_secret{};
    std::array<std::uint32_t, h> mask{};
    std::uint32_t dot = 0;
    for (std::size_t i = 0; i < h; i++) {
        lwe_secret[i] = 1;
        mask[i] = static_cast<std::uint32_t>((17 * i + 3) % q);
        dot += mask[i];
    }
    const RNSBootstrappingKey bootstrapping_key =
        BinaryNTTBootstrappingKeyGen(high_ring, binary_secret, lwe_secret,
                                     0.0, engine);
    std::vector<std::int64_t> small_coefficients(high_ring.degree());
    for (std::size_t i = 0; i < small_coefficients.size(); i++)
        small_coefficients[i] = static_cast<std::int64_t>(i % 3) - 1;
    const RNSSecret high_small_secret =
        QuadraticHintSmallSecretFromCoefficients(high_ring,
                                                 small_coefficients);
    const RNSKeySwitchKey boundary_key = BinaryNTTKeySwitchKeyGen(
        high_ring, binary_secret, high_small_secret, 4, 0.0, engine);
    constexpr TFHEpp::modular_ntt::PrimeModulus low_prime = {
        1125899906949121ULL, 14};
    Ring low_ring(4096, low_prime);
    const Secret low_small_secret = QuadraticHintSmallSecretFromCoefficients(
        low_ring, small_coefficients);
    const LWEKeySwitchKey output_key = LWEKeySwitchKeyGen(
        small_coefficients, lwe_secret, low_ring.modulus(), 4, 0.0, engine);
    const auto lut = Algorithm3BinarySignLUT(high_ring, q);
    const LWEPhase input = {
        std::vector<std::uint32_t>(mask.begin(), mask.end()),
        static_cast<std::uint32_t>((dot + 64) % q), q};

    TwoStageStats stats;
    const Ciphertext warm_accumulator = Algorithm3TwoStageBlindRotate(
        high_ring, binary_secret, bootstrapping_key, input, lut, 8, 4,
        boundary_key, high_small_secret, low_ring, low_small_secret, &stats);
    const LWEPhase warm_output = ModulusSwitch(
        LWEKeySwitch(SampleExtract(low_ring, warm_accumulator), output_key), q);
    if (LWEDecryptPhase(warm_output, lwe_secret) >= q / 2) {
        std::cerr << "Algorithm 3 benchmark setup failed" << std::endl;
        return 1;
    }

    Ciphertext accumulator;
    LWEPhase output;
    const double accumulator_ms = AverageMilliseconds(repetitions, [&] {
        accumulator = Algorithm3TwoStageBlindRotate(
            high_ring, binary_secret, bootstrapping_key, input, lut, 8, 4,
            boundary_key, high_small_secret, low_ring, low_small_secret);
    });
    const double output_ms = AverageMilliseconds(repetitions, [&] {
        output = ModulusSwitch(
            LWEKeySwitch(SampleExtract(low_ring, accumulator), output_key), q);
    });
    const double full_ms = AverageMilliseconds(repetitions, [&] {
        accumulator = Algorithm3TwoStageBlindRotate(
            high_ring, binary_secret, bootstrapping_key, input, lut, 8, 4,
            boundary_key, high_small_secret, low_ring, low_small_secret);
        output = ModulusSwitch(
            LWEKeySwitch(SampleExtract(low_ring, accumulator), output_key), q);
    });

    std::cout << "repetitions=" << repetitions << '\n'
              << "ring_degree=4096\n"
              << "lwe_weight=29\n"
              << "windows=8,4\n"
              << "high_product_layers="
              << stats.high_product.pointwise_multiplication_layers << '\n'
              << "modulus_boundaries=" << stats.modulus_boundaries << '\n'
              << "low_product_layers="
              << stats.low_product.pointwise_multiplication_layers << '\n'
              << "accumulator_ms=" << accumulator_ms << '\n'
              << "extract_keyswitch_modswitch_ms=" << output_ms << '\n'
              << "full_refresh_ms=" << full_ms << '\n';
}
