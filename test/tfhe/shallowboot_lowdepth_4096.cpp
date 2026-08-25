#include <array>
#include <cassert>
#include <cstdint>
#include <random>
#include <vector>

#include <tfhe/shallowboot_lowdepth.hpp>

int main()
{
    using namespace TFHEpp::shallowboot::lowdepth;
    constexpr std::size_t lwe_weight = 29;
    constexpr std::uint32_t lwe_modulus = 512;
    const std::array<TFHEpp::modular_ntt::PrimeModulus, 2> primes = {
        TFHEpp::modular_ntt::wide_primes[0],
        TFHEpp::modular_ntt::wide_primes[2]};
    RNSRing ring(4096, primes);
    std::mt19937_64 engine(20261730);
    const RNSSecret secret = BinaryNTTSecretGen(ring, engine);

    std::array<std::uint8_t, lwe_weight> lwe_secret{};
    std::array<std::uint32_t, lwe_weight> mask{};
    std::uint32_t dot_product = 0;
    for (std::size_t i = 0; i < lwe_weight; i++) {
        lwe_secret[i] = 1;
        mask[i] = static_cast<std::uint32_t>((17 * i + 3) % lwe_modulus);
        dot_product += mask[i];
    }
    constexpr std::uint32_t message_phase = 37;
    LWEPhase lwe = {std::vector<std::uint32_t>(mask.begin(), mask.end()),
                    static_cast<std::uint32_t>((dot_product + message_phase) %
                                               lwe_modulus),
                    lwe_modulus};

    const RNSBootstrappingKey bootstrapping_key =
        BinaryNTTBootstrappingKeyGen(ring, secret, lwe_secret, 0.0, engine);
    const unsigned __int128 high_modulus =
        static_cast<unsigned __int128>(ring[0].modulus()) * ring[1].modulus();
    const unsigned __int128 high_scale = high_modulus / 8;
    std::vector<std::vector<std::uint64_t>> test_vector(
        ring.levels(), std::vector<std::uint64_t>(ring.degree()));
    for (std::size_t level = 0; level < ring.levels(); level++) {
        test_vector[level][0] = static_cast<std::uint64_t>(
            high_scale % ring[level].modulus());
        ring[level].forward(test_vector[level]);
    }

    ProductTreeStats stats;
    const RNSCiphertext result = Algorithm3BlindRotate(
        ring, secret, bootstrapping_key, lwe, test_vector, &stats);
    const auto decrypted = DecryptNTT(ring, secret, result);
    constexpr std::size_t scale = 2 * 4096 / lwe_modulus;
    for (std::size_t level = 0; level < ring.levels(); level++) {
        auto expected = ring[level].monomialNTT(
            (lwe_modulus - message_phase) * scale);
        const std::uint64_t plaintext_scale = static_cast<std::uint64_t>(
            high_scale % ring[level].modulus());
        for (std::uint64_t &value : expected)
            value = TFHEpp::modular_ntt::multiply(
                value, plaintext_scale, ring[level].modulus());
        assert(decrypted[level] == expected);
    }
    assert(stats.pointwise_multiplication_layers == 5);
    assert(stats.pointwise_ciphertext_multiplications == lwe_weight);
    assert(stats.ntt_layers_inside_tree == 0);

    std::vector<std::int64_t> small_coefficients(ring.degree());
    for (std::size_t i = 0; i < small_coefficients.size(); i++)
        small_coefficients[i] = static_cast<std::int64_t>(i % 3) - 1;
    const RNSSecret small_secret =
        QuadraticHintSmallSecretFromCoefficients(ring, small_coefficients);
    const RNSKeySwitchKey boundary_key = BinaryNTTKeySwitchKeyGen(
        ring, secret, small_secret, 4, 0.0, engine);
    constexpr TFHEpp::modular_ntt::PrimeModulus low_prime = {
        1125899906949121ULL, 14};
    Ring low_ring(4096, low_prime);
    const Secret low_secret = QuadraticHintSmallSecretFromCoefficients(
        low_ring, small_coefficients);
    const Ciphertext low_result = Algorithm3ModulusBoundary(
        ring, secret, boundary_key, small_secret, low_ring, low_secret, result);
    const auto low_message = DecryptNTT(low_ring, low_secret, low_result);
    auto low_expected = low_ring.monomialNTT(
        (lwe_modulus - message_phase) * scale);
    for (std::uint64_t &value : low_expected)
        value = TFHEpp::modular_ntt::multiply(
            value, low_ring.modulus() / 8, low_ring.modulus());
    assert(low_message.size() == low_expected.size());

    TwoStageStats two_stage_stats;
    const Ciphertext two_stage_result = Algorithm3TwoStageBlindRotate(
        ring, secret, bootstrapping_key, lwe, test_vector, 8, 4,
        boundary_key, small_secret, low_ring, low_secret, &two_stage_stats);
    const auto two_stage_message =
        DecryptNTT(low_ring, low_secret, two_stage_result);
    assert(two_stage_message.size() == low_expected.size());
    assert(two_stage_stats.high_product.pointwise_multiplication_layers == 3);
    assert(two_stage_stats.high_product.pointwise_ciphertext_multiplications ==
           26);
    assert(two_stage_stats.modulus_boundaries == 4);
    assert(two_stage_stats.low_product.pointwise_multiplication_layers == 2);
    assert(two_stage_stats.low_product.pointwise_ciphertext_multiplications ==
           3);

    const LWEKeySwitchKey output_key = LWEKeySwitchKeyGen(
        small_coefficients, lwe_secret, low_ring.modulus(), 4, 0.0, engine);
    const auto sign_lut = Algorithm3BinarySignLUT(ring, lwe_modulus);
    for (const bool plain : {false, true}) {
        const std::uint32_t phase_message = plain ? 64 : lwe_modulus - 64;
        LWEPhase input = {
            std::vector<std::uint32_t>(mask.begin(), mask.end()),
            static_cast<std::uint32_t>((dot_product + phase_message) %
                                       lwe_modulus),
            lwe_modulus};
        const LWEPhase output = Algorithm3Bootstrap(
            ring, secret, bootstrapping_key, input, sign_lut, 8, 4,
            boundary_key, small_secret, low_ring, low_secret, output_key,
            lwe_modulus);
        const std::uint32_t output_phase = LWEDecryptPhase(output, lwe_secret);
        assert((output_phase < lwe_modulus / 2) == plain);
    }

    // h=37 passes the local sparse-LWE screen.  Its PBC has k=h+3=40
    // leaves; using 16 then 4 moves the additional layer above the modulus
    // boundary and leaves only two low-modulus product layers.
    constexpr std::size_t pbc_dimension = 1450;
    constexpr std::size_t pbc_weight = 37;
    std::vector<std::uint8_t> pbc_lwe_secret(pbc_dimension, 0);
    std::vector<std::uint32_t> pbc_support;
    pbc_support.reserve(pbc_weight);
    for (std::size_t i = 0; i < pbc_weight; i++) {
        const std::uint32_t index =
            static_cast<std::uint32_t>((47 * i + 11) % pbc_dimension);
        pbc_support.push_back(index);
        pbc_lwe_secret[index] = 1;
    }
    const PBCSchedule pbc_schedule = BuildPBCScheduleWithRetry(
        pbc_dimension, 40, 3, pbc_support, 20261730);
    const RNSPBCBootstrappingKey pbc_bootstrapping_key =
        PBCBootstrappingKeyGen(ring, secret, pbc_schedule, 0.0, engine);
    std::vector<std::uint32_t> pbc_mask(pbc_dimension);
    std::uint32_t pbc_dot_product = 0;
    for (std::size_t i = 0; i < pbc_dimension; i++) {
        pbc_mask[i] = static_cast<std::uint32_t>((101 * i + 17) % lwe_modulus);
        pbc_dot_product += pbc_mask[i] * pbc_lwe_secret[i];
    }
    const std::array<std::uint8_t, 1> pbc_output_secret = {0};
    const LWEKeySwitchKey pbc_output_key = LWEKeySwitchKeyGen(
        small_coefficients, pbc_output_secret, low_ring.modulus(), 4, 0.0,
        engine);
    for (const bool plain : {false, true}) {
        const std::uint32_t phase_message = plain ? 64 : lwe_modulus - 64;
        const LWEPhase pbc_input = {
            pbc_mask,
            static_cast<std::uint32_t>((pbc_dot_product + phase_message) %
                                       lwe_modulus),
            lwe_modulus};
        TwoStageStats pbc_stats;
        const LWEPhase pbc_output = Algorithm3PBCBootstrap(
            ring, secret, pbc_schedule, pbc_bootstrapping_key, pbc_input,
            sign_lut, 16, 4, boundary_key, small_secret, low_ring, low_secret,
            pbc_output_key, lwe_modulus, &pbc_stats);
        const std::uint32_t pbc_phase =
            LWEDecryptPhase(pbc_output, pbc_output_secret);
        assert((pbc_phase < lwe_modulus / 2) == plain);
        assert(pbc_stats.high_product.pointwise_multiplication_layers == 4);
        assert(pbc_stats.high_product.pointwise_ciphertext_multiplications ==
               37);
        assert(pbc_stats.modulus_boundaries == 3);
        assert(pbc_stats.low_product.pointwise_multiplication_layers == 2);
        assert(pbc_stats.low_product.pointwise_ciphertext_multiplications ==
               2);
    }

}
