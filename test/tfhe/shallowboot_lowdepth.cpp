#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <vector>

#include <tfhe/shallowboot_lowdepth.hpp>

int main()
{
    using namespace TFHEpp::shallowboot::lowdepth;
    Ring ring(32, TFHEpp::modular_ntt::wide_primes[0]);
    std::mt19937_64 engine(20261730);
    const Secret secret = BinaryNTTSecretGen(ring, engine);

    std::vector<std::uint64_t> initial_coeff(ring.degree());
    initial_coeff[0] = 7;
    std::vector<std::uint64_t> expected = initial_coeff;
    ring.forward(expected);

    std::vector<Ciphertext> factors;
    factors.push_back(TrivialEncryptNTT(ring, expected));
    for (std::size_t index = 0; index < 8; index++) {
        std::vector<std::uint64_t> bit_message(ring.degree());
        bit_message[index] = 1;
        ring.forward(bit_message);
        Ciphertext bit_ct = EncryptNTT(ring, secret, bit_message, 0.0, engine);
        const std::size_t exponent = 1 + 3 * index;
        factors.push_back(MonomialFactor(ring, bit_ct, exponent));

        std::vector<std::uint64_t> factor_plain =
            DecryptNTT(ring, secret, factors.back());
        for (std::size_t i = 0; i < ring.degree(); i++)
            expected[i] = TFHEpp::modular_ntt::multiply(
                expected[i], factor_plain[i], ring.modulus());
    }

    ProductTreeStats stats;
    const Ciphertext result = ProductTree(ring, secret, std::move(factors),
                                          &stats);
    const std::vector<std::uint64_t> decrypted =
        DecryptNTT(ring, secret, result);
    assert(decrypted == expected);
    assert(stats.pointwise_multiplication_layers == 4);
    assert(stats.pointwise_ciphertext_multiplications == 8);
    assert(stats.ntt_layers_inside_tree == 0);

    std::vector<std::uint64_t> accumulator_message(ring.degree());
    accumulator_message[0] = 11;
    ring.forward(accumulator_message);
    const Ciphertext accumulator = TrivialEncryptNTT(ring, accumulator_message);
    std::array<std::size_t, 7> exponents = {2, 5, 7, 11, 13, 17, 19};
    std::vector<Ciphertext> encrypted_bits;
    std::vector<std::uint64_t> blind_rotate_expected = accumulator_message;
    for (std::size_t index = 0; index < exponents.size(); index++) {
        std::vector<std::uint64_t> bit_message(ring.degree(),
                                               index % 2 == 0 ? 1 : 0);
        encrypted_bits.push_back(
            EncryptNTT(ring, secret, bit_message, 0.0, engine));
        const Ciphertext factor =
            MonomialFactor(ring, encrypted_bits.back(), exponents[index]);
        const auto factor_message = DecryptNTT(ring, secret, factor);
        for (std::size_t coefficient = 0; coefficient < ring.degree();
             coefficient++)
            blind_rotate_expected[coefficient] =
                TFHEpp::modular_ntt::multiply(
                    blind_rotate_expected[coefficient],
                    factor_message[coefficient], ring.modulus());
    }
    ProductTreeStats blind_rotate_stats;
    const Ciphertext blind_rotate = LowDepthBlindRotate(
        ring, secret, accumulator, encrypted_bits, exponents,
        &blind_rotate_stats);
    assert(DecryptNTT(ring, secret, blind_rotate) == blind_rotate_expected);
    assert(blind_rotate_stats.pointwise_multiplication_layers == 3);
    assert(blind_rotate_stats.pointwise_ciphertext_multiplications == 7);

    std::vector<Ciphertext> scheduled_inputs;
    for (std::size_t i = 0; i < 30; i++) {
        std::vector<std::uint64_t> message(ring.degree());
        message[0] = 1;
        ring.forward(message);
        scheduled_inputs.push_back(EncryptNTT(ring, secret, message, 0.0,
                                              engine));
    }
    constexpr std::array<std::uint32_t, 2> schedule = {8, 4};
    std::uint32_t boundary_calls = 0;
    ModulusScheduleStats schedule_stats;
    const Ciphertext scheduled_product = ScheduledProductTree(
        ring, secret, std::move(scheduled_inputs), schedule,
        [&boundary_calls](Ciphertext &) { boundary_calls++; }, &schedule_stats);
    const auto scheduled_message = DecryptNTT(ring, secret, scheduled_product);
    std::vector<std::uint64_t> scheduled_expected(ring.degree());
    scheduled_expected[0] = 1;
    ring.forward(scheduled_expected);
    const std::vector<std::uint64_t> scheduled_base = scheduled_expected;
    for (std::size_t power = 1; power < 30; power++)
        for (std::size_t coefficient = 0; coefficient < ring.degree();
             coefficient++)
            scheduled_expected[coefficient] = TFHEpp::modular_ntt::multiply(
                scheduled_expected[coefficient], scheduled_base[coefficient],
                ring.modulus());
    assert(scheduled_message == scheduled_expected);
    assert(schedule_stats.pointwise_multiplication_layers == 5);
    assert(schedule_stats.pointwise_ciphertext_multiplications == 29);
    assert(schedule_stats.modulus_boundaries == 4);
    assert(boundary_calls == 4);

    const std::array<std::uint8_t, 7> lwe_bits = {1, 0, 1, 0, 1, 0, 1};
    const BootstrappingKey bootstrapping_key = BinaryNTTBootstrappingKeyGen(
        ring, secret, lwe_bits, 0.0, engine);
    const Ciphertext keyed_blind_rotate = LowDepthBlindRotate(
        ring, secret, accumulator, bootstrapping_key, exponents);
    assert(DecryptNTT(ring, secret, keyed_blind_rotate) ==
           blind_rotate_expected);

    const Secret switched_secret = BinaryNTTSecretGen(ring, engine);
    const KeySwitchKey switch_key = BinaryNTTKeySwitchKeyGen(
        ring, secret, switched_secret, 4, 0.0, engine);
    const Ciphertext switched_ciphertext =
        BinaryNTTKeySwitch(ring, switch_key, keyed_blind_rotate);
    assert(DecryptNTT(ring, switched_secret, switched_ciphertext) ==
           blind_rotate_expected);

    const Secret small_secret =
        QuadraticHintSmallSecretGen(ring, 1, engine);
    const KeySwitchKey small_switch_key = BinaryNTTKeySwitchKeyGen(
        ring, secret, small_secret, 4, 0.0, engine);
    const Ciphertext small_switched =
        BinaryNTTKeySwitch(ring, small_switch_key, keyed_blind_rotate);
    assert(DecryptNTT(ring, small_secret, small_switched) ==
           blind_rotate_expected);

    constexpr TFHEpp::modular_ntt::PrimeModulus low_prime = {
        1125899906949121ULL, 14};
    Ring low_ring(32, low_prime);
    std::vector<std::int64_t> small_coefficients(ring.degree());
    for (std::size_t i = 0; i < small_coefficients.size(); i++)
        small_coefficients[i] = static_cast<std::int64_t>(i % 3) - 1;
    const Secret small_secret_high =
        QuadraticHintSmallSecretFromCoefficients(ring, small_coefficients);
    const Secret small_secret_low = QuadraticHintSmallSecretFromCoefficients(
        low_ring, small_coefficients);
    const KeySwitchKey boundary_key = BinaryNTTKeySwitchKeyGen(
        ring, secret, small_secret_high, 4, 0.0, engine);
    std::vector<std::uint64_t> high_message(ring.degree());
    high_message[0] = ring.modulus() / 8;
    ring.forward(high_message);
    const Ciphertext high_ciphertext =
        EncryptNTT(ring, secret, high_message, 0.0, engine);
    const Ciphertext low_ciphertext = Algorithm3ModulusBoundary(
        ring, secret, boundary_key, small_secret_high, low_ring,
        small_secret_low, high_ciphertext);
    auto low_message = DecryptNTT(low_ring, small_secret_low, low_ciphertext);
    low_ring.inverse(low_message);
    const auto expected_low = static_cast<std::int64_t>(low_ring.modulus() / 8);
    const auto actual_low = static_cast<std::int64_t>(
        TFHEpp::modular_ntt::centeredResidue(low_message[0], low_ring.modulus()));
    assert(std::llabs(actual_low - expected_low) < 128);

    const LWECiphertext extracted = SampleExtract(low_ring, low_ciphertext);
    const std::array<std::uint8_t, 5> output_lwe_secret = {1, 0, 1, 1, 0};
    const LWEKeySwitchKey lwe_switch_key = LWEKeySwitchKeyGen(
        small_coefficients, output_lwe_secret, low_ring.modulus(), 4, 0.0,
        engine);
    const LWECiphertext switched_lwe = LWEKeySwitch(extracted, lwe_switch_key);
    const LWEPhase modswitched_lwe = ModulusSwitch(switched_lwe, 64);
    const std::uint32_t phase =
        LWEDecryptPhase(modswitched_lwe, output_lwe_secret);
    assert(std::llabs(static_cast<std::int64_t>(phase) - 8) < 2);

    const std::array<std::uint8_t, 8> lwe_secret = {1, 0, 1, 1,
                                                     0, 1, 0, 1};
    const BootstrappingKey algorithm3_key = BinaryNTTBootstrappingKeyGen(
        ring, secret, lwe_secret, 0.0, engine);
    LWEPhase lwe;
    lwe.modulus = 64;
    lwe.a = {3, 8, 13, 21, 5, 34, 55, 2};
    std::uint32_t dot_product = 0;
    for (std::size_t i = 0; i < lwe.a.size(); i++)
        dot_product += lwe.a[i] * lwe_secret[i];
    constexpr std::uint32_t message_phase = 9;
    lwe.b = (dot_product + message_phase) % lwe.modulus;
    std::vector<std::uint64_t> test_vector(ring.degree());
    test_vector[0] = 1;
    ring.forward(test_vector);
    ProductTreeStats algorithm3_stats;
    const Ciphertext algorithm3_result = Algorithm3BlindRotate(
        ring, secret, algorithm3_key, lwe, test_vector, &algorithm3_stats);
    const auto algorithm3_message = DecryptNTT(ring, secret, algorithm3_result);
    const auto expected_monomial =
        ring.monomialNTT((lwe.modulus - message_phase) % lwe.modulus);
    assert(algorithm3_message == expected_monomial);
    assert(algorithm3_stats.pointwise_multiplication_layers == 4);
    assert(algorithm3_stats.pointwise_ciphertext_multiplications == 8);

    const std::array<std::uint8_t, 16> pbc_lwe_secret = {
        0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
    const std::array<std::uint32_t, 4> pbc_support = {1, 4, 9, 14};
    const PBCSchedule pbc_schedule = BuildPBCSchedule(
        16, 5, 3, pbc_support, 20261730);
    const PBCBootstrappingKey pbc_key = PBCBootstrappingKeyGen(
        ring, secret, pbc_schedule, 0.0, engine);
    LWEPhase pbc_lwe;
    pbc_lwe.modulus = 64;
    pbc_lwe.a = {3, 5, 7, 11, 13, 17, 19, 23,
                 29, 31, 37, 41, 43, 47, 53, 59};
    std::uint32_t pbc_dot = 0;
    for (std::size_t i = 0; i < pbc_lwe.a.size(); i++)
        pbc_dot += pbc_lwe.a[i] * pbc_lwe_secret[i];
    constexpr std::uint32_t pbc_message_phase = 6;
    pbc_lwe.b = (pbc_dot + pbc_message_phase) % pbc_lwe.modulus;
    ProductTreeStats pbc_stats;
    const Ciphertext pbc_result = Algorithm3PBCBlindRotate(
        ring, secret, pbc_schedule, pbc_key, pbc_lwe, test_vector,
        &pbc_stats);
    const auto pbc_message = DecryptNTT(ring, secret, pbc_result);
    assert(pbc_message ==
           ring.monomialNTT((pbc_lwe.modulus - pbc_message_phase) %
                            pbc_lwe.modulus));
    assert(pbc_stats.pointwise_multiplication_layers == 3);
    assert(pbc_stats.pointwise_ciphertext_multiplications == 4);

    const std::array<TFHEpp::modular_ntt::PrimeModulus, 2> rns_primes = {
        TFHEpp::modular_ntt::wide_primes[0],
        TFHEpp::modular_ntt::wide_primes[1]};
    RNSRing rns_ring(32, rns_primes);
    const RNSSecret rns_secret = BinaryNTTSecretGen(rns_ring, engine);
    std::vector<std::vector<std::uint64_t>> rns_message(
        rns_ring.levels(), std::vector<std::uint64_t>(rns_ring.degree()));
    for (std::size_t level = 0; level < rns_ring.levels(); level++) {
        rns_message[level][0] = 3;
        rns_ring[level].forward(rns_message[level]);
    }
    std::vector<RNSCiphertext> rns_inputs;
    for (std::size_t input = 0; input < 5; input++)
        rns_inputs.push_back(
            EncryptNTT(rns_ring, rns_secret, rns_message, 0.0, engine));
    ProductTreeStats rns_stats;
    const RNSCiphertext rns_product =
        ProductTree(rns_ring, rns_secret, std::move(rns_inputs), &rns_stats);
    const auto rns_decrypted = DecryptNTT(rns_ring, rns_secret, rns_product);
    for (std::size_t level = 0; level < rns_ring.levels(); level++) {
        std::vector<std::uint64_t> expected_rns = rns_message[level];
        for (std::size_t power = 1; power < 5; power++)
            for (std::size_t i = 0; i < rns_ring.degree(); i++)
                expected_rns[i] = TFHEpp::modular_ntt::multiply(
                    expected_rns[i], rns_message[level][i],
                    rns_ring[level].modulus());
        assert(rns_decrypted[level] == expected_rns);
    }
    assert(rns_stats.pointwise_multiplication_layers == 3);
    assert(rns_stats.pointwise_ciphertext_multiplications == 4);

    const std::array<std::uint8_t, 5> rns_lwe_secret = {1, 0, 1, 0, 1};
    const RNSBootstrappingKey rns_bootstrapping_key =
        BinaryNTTBootstrappingKeyGen(rns_ring, rns_secret, rns_lwe_secret,
                                     0.0, engine);
    LWEPhase rns_lwe;
    rns_lwe.modulus = 64;
    rns_lwe.a = {4, 9, 17, 33, 7};
    std::uint32_t rns_dot_product = 0;
    for (std::size_t i = 0; i < rns_lwe.a.size(); i++)
        rns_dot_product += rns_lwe.a[i] * rns_lwe_secret[i];
    constexpr std::uint32_t rns_message_phase = 12;
    rns_lwe.b = (rns_dot_product + rns_message_phase) % rns_lwe.modulus;
    std::vector<std::vector<std::uint64_t>> rns_test_vector(
        rns_ring.levels(), std::vector<std::uint64_t>(rns_ring.degree()));
    for (std::size_t level = 0; level < rns_ring.levels(); level++) {
        rns_test_vector[level][0] = 1;
        rns_ring[level].forward(rns_test_vector[level]);
    }
    ProductTreeStats rns_blind_rotate_stats;
    const RNSCiphertext rns_blind_rotate = Algorithm3BlindRotate(
        rns_ring, rns_secret, rns_bootstrapping_key, rns_lwe,
        rns_test_vector, &rns_blind_rotate_stats);
    const auto rns_blind_rotate_message =
        DecryptNTT(rns_ring, rns_secret, rns_blind_rotate);
    for (std::size_t level = 0; level < rns_ring.levels(); level++) {
        const auto expected = rns_ring[level].monomialNTT(
            (rns_lwe.modulus - rns_message_phase) % rns_lwe.modulus);
        assert(rns_blind_rotate_message[level] == expected);
    }
    assert(rns_blind_rotate_stats.pointwise_multiplication_layers == 3);
    assert(rns_blind_rotate_stats.pointwise_ciphertext_multiplications == 5);


    std::vector<std::int64_t> rns_small_coefficients(rns_ring.degree());
    for (std::size_t i = 0; i < rns_small_coefficients.size(); i++)
        rns_small_coefficients[i] = static_cast<std::int64_t>(i % 3) - 1;
    const RNSSecret rns_small_secret =
        QuadraticHintSmallSecretFromCoefficients(rns_ring,
                                                 rns_small_coefficients);
    const RNSKeySwitchKey rns_boundary_key = BinaryNTTKeySwitchKeyGen(
        rns_ring, rns_secret, rns_small_secret, 4, 0.0, engine);
    Ring rns_low_ring(32, low_prime);
    const Secret rns_small_low_secret =
        QuadraticHintSmallSecretFromCoefficients(rns_low_ring,
                                                 rns_small_coefficients);
    std::vector<std::uint64_t> qh_left(rns_low_ring.degree());
    std::vector<std::uint64_t> qh_right(rns_low_ring.degree());
    qh_left[0] = 7;
    qh_right[0] = 9;
    rns_low_ring.forward(qh_left);
    rns_low_ring.forward(qh_right);
    const Ciphertext qh_product = RelinearizationFreeMultiply(
        rns_low_ring, rns_small_low_secret,
        EncryptNTT(rns_low_ring, rns_small_low_secret, qh_left, 0.0, engine),
        EncryptNTT(rns_low_ring, rns_small_low_secret, qh_right, 0.0, engine));
    std::vector<std::uint64_t> qh_expected(rns_low_ring.degree());
    for (std::size_t i = 0; i < rns_low_ring.degree(); i++)
        qh_expected[i] = TFHEpp::modular_ntt::multiply(
            qh_left[i], qh_right[i], rns_low_ring.modulus());
    assert(DecryptNTT(rns_low_ring, rns_small_low_secret, qh_product) ==
           qh_expected);
    const unsigned __int128 rns_modulus =
        static_cast<unsigned __int128>(rns_ring[0].modulus()) *
        rns_ring[1].modulus();
    std::vector<std::vector<std::uint64_t>> rns_high_message(
        rns_ring.levels(), std::vector<std::uint64_t>(rns_ring.degree()));
    const unsigned __int128 encoded = rns_modulus / 8;
    for (std::size_t level = 0; level < rns_ring.levels(); level++) {
        rns_high_message[level][0] = static_cast<std::uint64_t>(
            encoded % rns_ring[level].modulus());
        rns_ring[level].forward(rns_high_message[level]);
    }
    const RNSCiphertext rns_high_ciphertext = EncryptNTT(
        rns_ring, rns_secret, rns_high_message, 0.0, engine);
    const Ciphertext rns_low_ciphertext = Algorithm3ModulusBoundary(
        rns_ring, rns_secret, rns_boundary_key, rns_small_secret,
        rns_low_ring, rns_small_low_secret, rns_high_ciphertext);
    auto rns_low_message =
        DecryptNTT(rns_low_ring, rns_small_low_secret, rns_low_ciphertext);
    rns_low_ring.inverse(rns_low_message);
    const auto rns_actual = static_cast<std::int64_t>(
        TFHEpp::modular_ntt::centeredResidue(rns_low_message[0],
                                             rns_low_ring.modulus()));
    assert(std::llabs(rns_actual -
                      static_cast<std::int64_t>(rns_low_ring.modulus() / 8)) <
           128);

    const RNSPBCBootstrappingKey rns_pbc_key = PBCBootstrappingKeyGen(
        rns_ring, rns_secret, pbc_schedule, 0.0, engine);
    TwoStageStats pbc_two_stage_stats;
    const Ciphertext pbc_two_stage = Algorithm3PBCTwoStageBlindRotate(
        rns_ring, rns_secret, pbc_schedule, rns_pbc_key, pbc_lwe,
        rns_high_message, 8, 4, rns_boundary_key, rns_small_secret,
        rns_low_ring, rns_small_low_secret, &pbc_two_stage_stats);
    const auto pbc_two_stage_message =
        DecryptNTT(rns_low_ring, rns_small_low_secret, pbc_two_stage);
    assert(pbc_two_stage_message.size() == rns_low_ring.degree());
    assert(pbc_two_stage_stats.high_product.pointwise_multiplication_layers ==
           3);
    assert(pbc_two_stage_stats.high_product.pointwise_ciphertext_multiplications ==
           4);
    assert(pbc_two_stage_stats.modulus_boundaries == 1);

    const LWEKeySwitchKey pbc_output_key = LWEKeySwitchKeyGen(
        rns_small_coefficients, pbc_lwe_secret, rns_low_ring.modulus(), 4,
        0.0, engine);
    const auto pbc_sign_lut = Algorithm3BinarySignLUT(rns_ring, 64);
    for (const bool plain : {false, true}) {
        LWEPhase pbc_input = pbc_lwe;
        const std::uint32_t message = plain ? 16 : 64 - 16;
        pbc_input.b = (pbc_dot + message) % pbc_input.modulus;
        const Ciphertext pbc_accumulator = Algorithm3PBCTwoStageBlindRotate(
            rns_ring, rns_secret, pbc_schedule, rns_pbc_key, pbc_input,
            pbc_sign_lut, 8, 4, rns_boundary_key, rns_small_secret,
            rns_low_ring, rns_small_low_secret);
        const LWECiphertext pbc_extracted =
            SampleExtract(rns_low_ring, pbc_accumulator);
        const LWEPhase pbc_output = ModulusSwitch(
            LWEKeySwitch(pbc_extracted, pbc_output_key), 64);
        const std::uint32_t pbc_phase =
            LWEDecryptPhase(pbc_output, pbc_lwe_secret);
        assert((pbc_phase < 32) == plain);
    }

    const CoefficientKeySwitchKey coefficient_switch_key =
        CoefficientKeySwitchKeyGen(ring, secret, small_secret, 8, 0.0,
                                   engine);
    const Ciphertext coefficient_switched = CoefficientKeySwitch(
        ring, coefficient_switch_key, keyed_blind_rotate);
    assert(DecryptNTT(ring, small_secret, coefficient_switched) ==
           blind_rotate_expected);
}
