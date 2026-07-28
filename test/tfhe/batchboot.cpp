#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>
#include <tfhe++.hpp>

namespace {

using RingP = TFHEpp::lvl1param;
constexpr std::uint32_t slots = 8;
using ModuleP = TFHEpp::BatchRingSwitchP<RingP, slots>;

TFHEpp::Polynomial<ModuleP> ModulePhase(
    const TFHEpp::TRLWE<ModuleP> &ciphertext,
    const TFHEpp::Key<ModuleP> &key)
{
    TFHEpp::Polynomial<ModuleP> phase = ciphertext[ModuleP::k];
    for (std::uint32_t component = 0; component < ModuleP::k; component++) {
        for (std::uint32_t output = 0; output < ModuleP::n; output++) {
            typename ModuleP::T product = 0;
            for (std::uint32_t j = 0; j <= output; j++)
                product += ciphertext[component][j] *
                           key[component * ModuleP::n + output - j];
            for (std::uint32_t j = output + 1; j < ModuleP::n; j++)
                product -= ciphertext[component][j] *
                           key[component * ModuleP::n + ModuleP::n + output -
                               j];
            phase[output] -= product;
        }
    }
    return phase;
}

void TestKeySerialization()
{
    auto original =
        std::make_unique<TFHEpp::BatchBootKey<ModuleP, RingP>>();
    original->components[0].negative_gaps.emplace_back(
        std::make_unique<TFHEpp::BatchEMPKey<RingP>>());
    original->components[0].final_positive =
        std::make_unique<TFHEpp::BatchEMPKey<RingP>>();

    std::stringstream stream;
    {
        cereal::PortableBinaryOutputArchive archive(stream);
        archive(*original);
    }

    auto restored =
        std::make_unique<TFHEpp::BatchBootKey<ModuleP, RingP>>();
    {
        cereal::PortableBinaryInputArchive archive(stream);
        archive(*restored);
    }
    assert(restored->components[0].negative_gaps.size() == 1);
    assert(restored->components[0].negative_gaps[0]);
    assert(restored->components[0].final_positive);
    assert(!restored->components[1].final_positive);
}

void TestFDConjugation()
{
    TFHEpp::Polynomial<RingP> input{};
    input[14] = RingP::μ;
    TFHEpp::PolynomialInFD<RingP> transformed;
    TFHEpp::TwistIFFT<RingP>(transformed, input);
    TFHEpp::PolynomialInFD<RingP> conjugated;
    TFHEpp::BatchConjugateInFD<RingP>(conjugated, transformed);

    TFHEpp::Polynomial<RingP> actual;
    TFHEpp::TwistFFT<RingP>(actual, conjugated);
    TFHEpp::Polynomial<RingP> expected;
    TFHEpp::Automorphism<RingP>(expected, input, 2 * RingP::n - 1);
    for (std::uint32_t i = 0; i < RingP::n; i++) {
        const auto error = static_cast<std::int32_t>(actual[i] - expected[i]);
        if (error < -2 || error > 2)
            std::cerr << "FD conjugation mismatch at coefficient " << i
                      << ": actual " << actual[i] << ", expected "
                      << expected[i] << '\n';
        assert(error >= -2);
        assert(error <= 2);
    }
}

void TestEMP()
{
    TFHEpp::Key<RingP> key{};
    key[7] = key[123] = 1;

    constexpr std::uint32_t emp_slots = 8;
    constexpr std::uint32_t exponent = 5;
    std::vector<TFHEpp::TRLWE<RingP>> original(emp_slots);
    std::vector<TFHEpp::Polynomial<RingP>> expected(emp_slots);
    for (std::uint32_t i = 0; i < emp_slots; i++) {
        TFHEpp::Polynomial<RingP> plaintext{};
        plaintext[11 + i] = RingP::μ;
        TFHEpp::trlweSymEncrypt<RingP>(original[i], plaintext, 0.0, key);
        expected[i] = plaintext;
    }

    auto emp_key = std::make_unique<TFHEpp::BatchEMPKey<RingP>>();
    TFHEpp::BatchEMPKeyGen<RingP>(*emp_key, exponent, emp_slots, key);
    for (const bool positive : {true, false}) {
        auto ciphertexts = original;
        auto oracle = expected;
        TFHEpp::BatchEMP<RingP>(ciphertexts, *emp_key, positive);

        std::vector<TFHEpp::Polynomial<RingP>> shifted(emp_slots);
        for (std::uint32_t output = 0; output < emp_slots; output++) {
            const std::uint32_t source =
                positive ? (output + emp_slots - exponent) % emp_slots
                         : (output + exponent) % emp_slots;
            const bool inverse = positive ? output < exponent
                                          : output + exponent >= emp_slots;
            if (inverse)
                TFHEpp::Automorphism<RingP>(shifted[output], oracle[source],
                                             2 * RingP::n - 1);
            else
                shifted[output] = oracle[source];
        }

        for (std::uint32_t i = 0; i < emp_slots; i++) {
            const auto phase = TFHEpp::trlwePhase<RingP>(ciphertexts[i], key);
            for (std::uint32_t j = 0; j < RingP::n; j++) {
                const auto error = static_cast<std::int32_t>(phase[j] -
                                                              shifted[i][j]);
                if (error <= -RingP::μ / 4 || error >= RingP::μ / 4)
                    std::cerr << "EMP mismatch for "
                              << (positive ? "positive" : "negative")
                              << " rotation, output " << i
                              << ", coefficient " << j << ", error "
                              << error << '\n';
                assert(error > -RingP::μ / 4);
                assert(error < RingP::μ / 4);
            }
        }
    }
}

void TestSparseFunctionalBatchBoot()
{
    constexpr std::uint32_t input_bits = 3;
    constexpr std::uint32_t stride = RingP::n / slots;

    TFHEpp::Key<RingP> sparse_key{};
    sparse_key[5] = 1;
    sparse_key[3 * stride + 5] = 1;
    sparse_key[700] = 1;

    TFHEpp::Key<RingP> accumulator_key{};
    for (std::uint32_t i = 0; i < RingP::n; i += 97)
        accumulator_key[i] = 1;

    const std::array<std::uint32_t, slots> messages{0, 1, 2, 3, 3, 2, 1, 0};
    TFHEpp::Polynomial<RingP> plaintext{};
    for (std::uint32_t i = 0; i < slots; i++)
        plaintext[i * stride] =
            static_cast<RingP::T>(messages[i])
            << (std::numeric_limits<RingP::T>::digits - input_bits);

    TFHEpp::TRLWE<RingP> source_ciphertext;
    TFHEpp::trlweSymEncrypt<RingP>(source_ciphertext, plaintext, 0.0,
                                   sparse_key);

    auto module_key =
        TFHEpp::BatchRingSwitchSecret<RingP, slots>(sparse_key);
    TFHEpp::TRLWE<ModuleP> module_ciphertext;
    TFHEpp::BatchRingSwitch<RingP, slots>(module_ciphertext,
                                          source_ciphertext);
    const auto extracted_phase = ModulePhase(module_ciphertext, module_key);
    for (std::uint32_t i = 0; i < slots; i++)
        assert(extracted_phase[i] == plaintext[i * stride]);

    auto boot_key =
        std::make_unique<TFHEpp::BatchBootKey<ModuleP, RingP>>();
    TFHEpp::BatchBootKeyGen<ModuleP, RingP>(
        *boot_key, module_key, accumulator_key);
    auto automorphism_keys =
        std::make_unique<TFHEpp::AnnihilateKey<RingP>>();
    TFHEpp::annihilatekeygen<RingP>(*automorphism_keys, accumulator_key);

    const auto identity =
        TFHEpp::MakeBatchBootTestVector<RingP, input_bits>(
            [](const std::uint32_t value) { return value; });
    std::vector<TFHEpp::TRLWE<RingP>> accumulators;
    TFHEpp::BatchBootAccumulate<ModuleP, RingP>(
        accumulators, module_ciphertext, identity.polynomial, *boot_key,
        identity.exponent_bias);
    for (std::uint32_t i = 0; i < slots; i++) {
        const auto accumulator_phase =
            TFHEpp::trlwePhase<RingP>(accumulators[i], accumulator_key);
        const auto actual =
            TFHEpp::BatchTorusDecode<input_bits>(accumulator_phase[0]);
        if (actual != messages[i]) {
            constexpr std::uint32_t rotation_modulus = 2 * RingP::n;
            auto modswitched = TFHEpp::BatchModSwitch<ModuleP>(
                module_ciphertext, rotation_modulus);
            std::int64_t phase_exponent =
                modswitched[ModuleP::k][i] + identity.exponent_bias;
            for (std::uint32_t component = 0; component < ModuleP::k;
                 component++) {
                for (std::uint32_t j = 0; j <= i; j++)
                    phase_exponent -=
                        modswitched[component][j] *
                        module_key[component * slots + i - j];
                for (std::uint32_t j = i + 1; j < slots; j++)
                    phase_exponent +=
                        modswitched[component][j] *
                        module_key[component * slots + slots + i - j];
            }
            phase_exponent %= rotation_modulus;
            if (phase_exponent < 0) phase_exponent += rotation_modulus;
            TFHEpp::Polynomial<RingP> oracle;
            TFHEpp::PolynomialMulByXai<RingP>(
                oracle, identity.polynomial,
                static_cast<std::uint32_t>(phase_exponent));
            std::cerr << "BatchBoot accumulator mismatch in slot " << i
                      << ": expected " << messages[i] << ", got " << actual
                      << ", phase exponent " << phase_exponent
                      << ", oracle value "
                      << TFHEpp::BatchTorusDecode<input_bits>(oracle[0])
                      << ", torus error "
                      << static_cast<std::int32_t>(accumulator_phase[0] -
                                                   oracle[0])
                      << '\n';
        }
        assert(actual == messages[i]);
    }

    TFHEpp::TRLWE<RingP> refreshed;
    TFHEpp::BatchRepack<RingP>(refreshed, accumulators,
                               *automorphism_keys);

    const auto phase = TFHEpp::trlwePhase<RingP>(refreshed, accumulator_key);
    for (std::uint32_t i = 0; i < slots; i++) {
        const std::uint32_t actual =
            TFHEpp::BatchTorusDecode<input_bits>(phase[i * stride]);
        if (actual != messages[i])
            std::cerr << "BatchBoot mismatch in slot " << i << ": expected "
                      << messages[i] << ", got " << actual
                      << ", torus error "
                      << static_cast<std::int32_t>(
                             phase[i * stride] -
                             (static_cast<RingP::T>(messages[i])
                              << (std::numeric_limits<RingP::T>::digits -
                                  input_bits)))
                      << '\n';
        assert(actual == messages[i]);
    }
}

void TestRLWEKeySwitch()
{
    TFHEpp::Key<RingP> source_key{};
    TFHEpp::Key<RingP> target_key{};
    source_key[3] = source_key[91] = 1;
    target_key[7] = target_key[123] = 1;

    TFHEpp::Polynomial<RingP> plaintext{};
    plaintext[0] = RingP::μ;
    plaintext[19] = -RingP::μ;
    TFHEpp::TRLWE<RingP> input;
    TFHEpp::trlweSymEncrypt<RingP>(input, plaintext, 0.0, source_key);

    auto key =
        std::make_unique<TFHEpp::BatchRLWEKeySwitchKey<RingP>>();
    TFHEpp::BatchRLWEKeySwitchKeyGen<RingP>(*key, source_key, target_key);
    TFHEpp::TRLWE<RingP> output;
    TFHEpp::BatchRLWEKeySwitch<RingP>(output, input, *key);

    const auto decrypted = TFHEpp::trlweSymDecrypt<RingP>(output, target_key);
    assert(decrypted[0]);
    assert(!decrypted[19]);
}

void TestExternalProductTree()
{
    TFHEpp::Key<RingP> key{};
    key[7] = key[123] = 1;

    TFHEpp::Polynomial<RingP> one{};
    one[0] = 1;
    std::vector<TFHEpp::TRGSWFFT<RingP>> controls(2);
    for (auto &control : controls)
        TFHEpp::trgswSymEncrypt<RingP>(control, one, key);

    std::array<TFHEpp::Polynomial<RingP>, 2> test_polynomials{};
    test_polynomials[1][0] = RingP::μ;
    auto automorphism_keys =
        std::make_unique<TFHEpp::AnnihilateKey<RingP>>();
    TFHEpp::annihilatekeygen<RingP>(*automorphism_keys, key);

    TFHEpp::TRLWE<RingP> output;
    TFHEpp::BatchExternalProductTree<RingP>(
        output, controls, test_polynomials, 2, *automorphism_keys);
    const auto phase = TFHEpp::trlwePhase<RingP>(output, key);
    const auto error = static_cast<std::int32_t>(phase[0] - RingP::μ);
    if (error <= -RingP::μ / 4 || error >= RingP::μ / 4)
        std::cerr << "External-product tree mismatch, signed error " << error
                  << '\n';
    assert(error > -RingP::μ / 4);
    assert(error < RingP::μ / 4);
}

void TestHalfCircuitBootstrap()
{
    using CircuitP = TFHEpp::lvl2param;
    constexpr std::uint32_t circuit_slots = 2;
    constexpr std::uint32_t input_bits = 4;
    constexpr std::uint32_t stride = CircuitP::n / circuit_slots;
    using CircuitModuleP =
        TFHEpp::BatchRingSwitchP<CircuitP, circuit_slots>;

    TFHEpp::Key<CircuitP> sparse_key{};
    sparse_key[stride + 13] = 1;
    TFHEpp::Key<CircuitP> accumulator_key{};
    accumulator_key[7] = accumulator_key[509] = 1;

    constexpr std::array<std::uint32_t, circuit_slots> messages{1, 3};
    TFHEpp::Polynomial<CircuitP> plaintext{};
    for (std::uint32_t i = 0; i < circuit_slots; i++)
        plaintext[i * stride] =
            static_cast<CircuitP::T>(messages[i])
            << (std::numeric_limits<CircuitP::T>::digits - input_bits);

    TFHEpp::TRLWE<CircuitP> source;
    TFHEpp::trlweSymEncrypt<CircuitP>(source, plaintext, 0.0, sparse_key);
    auto module_key =
        TFHEpp::BatchRingSwitchSecret<CircuitP, circuit_slots>(sparse_key);
    TFHEpp::TRLWE<CircuitModuleP> module_ciphertext;
    TFHEpp::BatchRingSwitch<CircuitP, circuit_slots>(module_ciphertext,
                                                     source);

    auto boot_key = std::make_unique<
        TFHEpp::BatchBootKey<CircuitModuleP, CircuitP>>();
    TFHEpp::BatchBootKeyGen<CircuitModuleP, CircuitP>(
        *boot_key, module_key, accumulator_key);
    auto automorphism_keys =
        std::make_unique<TFHEpp::AnnihilateKey<CircuitP>>();
    TFHEpp::annihilatekeygen<CircuitP>(*automorphism_keys, accumulator_key);
    auto scheme_switch_key =
        std::make_unique<TFHEpp::BatchSchemeSwitchKey<CircuitP>>();
    TFHEpp::BatchSchemeSwitchKeyGen<CircuitP>(*scheme_switch_key,
                                               accumulator_key);

    std::vector<TFHEpp::TRGSW<CircuitP>> circuit_ciphertexts;
    TFHEpp::BatchHalfCircuitBootstrap<CircuitModuleP, CircuitP>(
        circuit_ciphertexts, module_ciphertext, *boot_key,
        *automorphism_keys, *scheme_switch_key);
    assert(circuit_ciphertexts.size() == circuit_slots);

    TFHEpp::Polynomial<CircuitP> unit_plaintext{};
    unit_plaintext[0] = CircuitP::μ;
    TFHEpp::TRLWE<CircuitP> unit;
    TFHEpp::trlweSymEncrypt<CircuitP>(unit, unit_plaintext, 0.0,
                                      accumulator_key);
    for (std::uint32_t i = 0; i < circuit_slots; i++) {
        auto control = std::make_unique<TFHEpp::TRGSWFFT<CircuitP>>();
        TFHEpp::ApplyFFT2trgsw<CircuitP>(*control, circuit_ciphertexts[i]);
        TFHEpp::TRLWE<CircuitP> product;
        TFHEpp::ExternalProduct<CircuitP>(product, unit, *control);
        const auto phase =
            TFHEpp::trlwePhase<CircuitP>(product, accumulator_key);
        const std::uint32_t expected_exponent =
            (2 * CircuitP::n * messages[i]) / (1U << input_bits);
        const auto error = static_cast<std::int64_t>(
            phase[expected_exponent] - CircuitP::μ);
        if (error <= -static_cast<std::int64_t>(CircuitP::μ / 4) ||
            error >= static_cast<std::int64_t>(CircuitP::μ / 4))
            std::cerr << "Half.BatchCBS mismatch in slot " << i
                      << ", coefficient " << expected_exponent
                      << ", signed error " << error << '\n';
        assert(error > -static_cast<std::int64_t>(CircuitP::μ / 4));
        assert(error < static_cast<std::int64_t>(CircuitP::μ / 4));
    }

    std::array<TFHEpp::Polynomial<CircuitP>, 1> unary_test_polynomial{};
    unary_test_polynomial[0][0] = CircuitP::μ;
    TFHEpp::TRLWE<CircuitP> composed;
    TFHEpp::BatchCircuitBootstrap<CircuitModuleP, CircuitP>(
        composed, module_ciphertext, *boot_key, *automorphism_keys,
        *scheme_switch_key, 1, unary_test_polynomial, 2);
    const auto composed_phase =
        TFHEpp::trlwePhase<CircuitP>(composed, accumulator_key);
    const std::uint32_t expected_exponent =
        (2 * CircuitP::n * messages[0]) / (1U << input_bits);
    const auto composed_error = static_cast<std::int64_t>(
        composed_phase[expected_exponent] - CircuitP::μ);
    assert(composed_error > -static_cast<std::int64_t>(CircuitP::μ / 4));
    assert(composed_error < static_cast<std::int64_t>(CircuitP::μ / 4));
}

}  // namespace

int main()
{
    TestKeySerialization();
    TestFDConjugation();
    TestEMP();
    TestSparseFunctionalBatchBoot();
    TestRLWEKeySwitch();
    TestExternalProductTree();
    TestHalfCircuitBootstrap();
    std::cout << "BatchBoot passed" << std::endl;
}
