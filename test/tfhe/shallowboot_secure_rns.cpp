#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <random>
#include <iostream>
#include <vector>

#include <tfhe/shallowboot_lowdepth.hpp>

int main()
{
    using namespace TFHEpp::shallowboot::lowdepth;
    constexpr std::size_t degree = 8192;
    const std::array<RNSRing, 4> rings = {
        RNSRing(degree, secure_q151_primes),
        RNSRing(degree, secure_q69_primes),
        RNSRing(degree, secure_q52_primes),
        RNSRing(degree, secure_q36_primes)};
    std::vector<std::int64_t> small_secret(degree);
    for (std::size_t i = 0; i < degree; i++)
        small_secret[i] = static_cast<std::int64_t>(i % 3) - 1;
    std::array<RNSSecret, 4> secrets = {
        QuadraticHintSmallSecretFromCoefficients(rings[0], small_secret),
        QuadraticHintSmallSecretFromCoefficients(rings[1], small_secret),
        QuadraticHintSmallSecretFromCoefficients(rings[2], small_secret),
        QuadraticHintSmallSecretFromCoefficients(rings[3], small_secret)};

    std::mt19937_64 engine(20261730);
    {
        std::vector<std::vector<std::uint64_t>> message(
            rings[0].levels(), std::vector<std::uint64_t>(degree));
        for (std::size_t level = 0; level < rings[0].levels(); level++) {
            message[level][0] = 1;
            rings[0][level].forward(message[level]);
        }
        const RNSCiphertext fresh = EncryptLSBNTT(
            rings[0], secrets[0], message, 4, 0.75, engine);
        const RNSCiphertext boundary =
            ModulusSwitch(rings[0], rings[1], fresh, 4);
        auto boundary_phase = DecryptNTT(rings[1], secrets[1], boundary);
        rings[1][0].inverse(boundary_phase[0]);
        const std::uint64_t modulus = rings[1][0].modulus();
        const std::int64_t centered = boundary_phase[0][0] <= modulus / 2
            ? static_cast<std::int64_t>(boundary_phase[0][0])
            : -static_cast<std::int64_t>(modulus - boundary_phase[0][0]);
        if ((centered % 4 + 4) % 4 != 1) return 10;
        for (std::size_t i = 1; i < degree; i++) {
            const std::int64_t value = boundary_phase[0][i] <= modulus / 2
                ? static_cast<std::int64_t>(boundary_phase[0][i])
                : -static_cast<std::int64_t>(modulus - boundary_phase[0][i]);
            if ((value % 4 + 4) % 4 != 0) return 20;
        }
    }
    std::vector<RNSCiphertext> factors;
    factors.reserve(40);
    // A PBC bucket has at most 110 encrypted entries. Sampling one Gaussian
    // with their aggregate standard deviation is distribution-equivalent for
    // this phase/noise regression and avoids constructing a multi-GB key.
    const double aggregate_sigma = 0.75 * std::sqrt(110.0);
    for (std::size_t factor = 0; factor < 40; factor++) {
        std::vector<std::vector<std::uint64_t>> message(
            rings[0].levels(), std::vector<std::uint64_t>(degree));
        for (std::size_t level = 0; level < rings[0].levels(); level++) {
            message[level][0] = 1;
            rings[0][level].forward(message[level]);
        }
        factors.push_back(EncryptLSBNTT(
            rings[0], secrets[0], message, 4, aggregate_sigma, engine));
    }
    constexpr std::array<std::uint32_t, 4> windows = {8, 2, 2, 2};
    {
        std::vector<RNSCiphertext> diagnostic = factors;
        for (std::size_t stage = 0; stage < rings.size(); stage++) {
            std::vector<RNSCiphertext> outputs;
            for (std::size_t first = 0; first < diagnostic.size();
                 first += windows[stage]) {
                const std::size_t last = std::min<std::size_t>(
                    first + windows[stage], diagnostic.size());
                std::vector<RNSCiphertext> group;
                for (std::size_t i = first; i < last; i++)
                    group.push_back(diagnostic[i]);
                RNSCiphertext product = ProductTree(
                    rings[stage], secrets[stage], std::move(group));
                if (stage + 1 < rings.size() &&
                    !IsSameRNSBasis(rings[stage], rings[stage + 1]))
                    product = ModulusSwitch(
                        rings[stage], rings[stage + 1], product, 4);
                outputs.push_back(std::move(product));
            }
            const RNSRing &check_ring =
                stage + 1 < rings.size() ? rings[stage + 1] : rings[stage];
            const RNSSecret &check_secret =
                stage + 1 < rings.size() ? secrets[stage + 1] : secrets[stage];
            for (const RNSCiphertext &output : outputs) {
                auto phase_check = DecryptNTT(check_ring, check_secret, output);
                check_ring[0].inverse(phase_check[0]);
                const std::uint64_t modulus = check_ring[0].modulus();
                for (std::size_t i = 0; i < degree; i++) {
                    const std::int64_t value = phase_check[0][i] <= modulus / 2
                        ? static_cast<std::int64_t>(phase_check[0][i])
                        : -static_cast<std::int64_t>(
                              modulus - phase_check[0][i]);
                    const std::uint32_t wanted = i == 0 ? 1 : 0;
                    if (static_cast<std::uint32_t>((value % 4 + 4) % 4) !=
                        wanted)
                        return static_cast<int>(21 + stage);
                }
            }
            diagnostic = std::move(outputs);
        }
    }
    MultiStageStats stats;
    const RNSCiphertext result = Algorithm3QHMultiStageProductRNS(
        rings, secrets, std::move(factors), windows, 4, &stats);
    const auto phase = DecryptNTT(rings.back(), secrets.back(), result);
    std::vector<std::uint64_t> coefficients = phase.front();
    rings.back()[0].inverse(coefficients);
    // BGV/LSB boundaries preserve the selector polynomial modulo t=2.
    const std::uint64_t final_modulus = rings.back()[0].modulus();
    const auto parity = [final_modulus](const std::uint64_t residue) {
        const std::int64_t centered = residue <= final_modulus / 2
            ? static_cast<std::int64_t>(residue)
            : -static_cast<std::int64_t>(final_modulus - residue);
        return static_cast<std::uint32_t>((centered % 4 + 4) % 4);
    };
    if (parity(coefficients[0]) != 1) return 11;
    for (std::size_t i = 1; i < degree; i++)
        if (parity(coefficients[i]) != 0) return 12;
    if (stats.products.size() != 4 ||
        stats.products[0].pointwise_multiplication_layers != 3 ||
        stats.products[1].pointwise_multiplication_layers != 1 ||
        stats.products[2].pointwise_multiplication_layers != 1 ||
        stats.products[3].pointwise_multiplication_layers != 1 ||
        stats.modulus_boundaries != std::vector<std::uint32_t>({5, 3, 2}))
        return 13;

    // Finish the intended Algorithm-3 encoding path: modular inverse converts
    // LSB selectors to MSB without scaling their error, then a plaintext LUT
    // rotation produces an ordinary extractable TFHE accumulator.
    const RNSCiphertext msb_selector = LSBToMSB(rings.back(), result, 4);
    const auto lut = Algorithm3BinarySignPlaintextLUT(rings.back(), 512);
    const RNSCiphertext accumulator =
        MultiplyByPlaintextNTT(rings.back(), msb_selector, lut);
    const LWECiphertext extracted =
        SampleExtract(rings.back()[0], accumulator.residues[0]);
    const LWEPhase output = ModulusSwitch(extracted, 512);
    // Coefficient zero of the generated sign LUT is +Q/4 after multiplying
    // its -1 plaintext coefficient by the converted selector's -Q/4 phase.
    // The aggregate PBC noise fixture above models up to 110 encrypted
    // entries. Allow its deterministic 112-bin phase budget around Q/4.
    if (output.b < 16 || output.b > 240) {
        std::cerr << "unexpected sign-LUT phase: " << output.b << std::endl;
        return 14;
    }
    auto opposite_lut = lut;
    for (std::size_t level = 0; level < rings.back().levels(); level++)
        for (std::uint64_t &value : opposite_lut[level])
            value = TFHEpp::modular_ntt::negate(
                value, rings.back()[level].modulus());
    const RNSCiphertext opposite_accumulator =
        MultiplyByPlaintextNTT(rings.back(), msb_selector, opposite_lut);
    const LWEPhase opposite_output = ModulusSwitch(
        SampleExtract(rings.back()[0], opposite_accumulator.residues[0]), 512);
    if (opposite_output.b < 272 || opposite_output.b > 496) return 16;

    // Exercise the public PBC frontend and complete bootstrap wrapper. A
    // compact identity schedule keeps this integration test small; the first
    // half above separately exercises secure-parameter aggregate noise.
    PBCSchedule schedule;
    schedule.source_dimension = 40;
    schedule.bucket_count = 40;
    schedule.copies = 1;
    schedule.bucket_indices.resize(40);
    schedule.selected_source.resize(40);
    for (std::uint32_t i = 0; i < 40; i++) {
        schedule.bucket_indices[i] = {i};
        schedule.selected_source[i] = static_cast<std::int32_t>(i);
    }
    const RNSPBCBootstrappingKey pbc_key = PBCBootstrappingKeyGenLSB(
        rings[0], secrets[0], schedule, 4, 0.75, engine);
    const LWEPhase input = {std::vector<std::uint32_t>(40, 0), 0, 512};
    const std::array<std::uint8_t, 1> target_secret = {0};
    const LWEKeySwitchKey output_key = LWEKeySwitchKeyGen(
        small_secret, target_secret, 1U << 12, 3, 0.0,
        engine);
    MultiStageStats wrapper_stats;
    const LWEPhase wrapped = Algorithm3PBCQHMultiStageBootstrap(
        rings, secrets, schedule, pbc_key, input, lut, windows, output_key,
        512, 1U << 12, 4, &wrapper_stats);
    if (wrapped.a.size() != 1 || wrapped.b < 64 || wrapped.b > 192) {
        std::cerr << "wrapped phase=" << wrapped.b << '\n';
        return 15;
    }

    // Signed base-4 decomposition is algebraically exact even though its
    // noisy secure-parameter trial currently lacks sufficient margin.
    const LWEKeySwitchKey balanced_key = LWEKeySwitchKeyGen(
        small_secret, target_secret, 1U << 15, 2, 0.0, engine, 8, true);
    LWECiphertext probe;
    probe.modulus = 1U << 15;
    probe.a.resize(degree);
    probe.b = probe.modulus / 4;
    for (std::size_t i = 0; i < degree; i++) {
        probe.a[i] = (17 * i + 5) % probe.modulus;
        const std::int64_t key_value = small_secret[i];
        const std::uint64_t product = key_value >= 0
            ? probe.a[i] * static_cast<std::uint64_t>(key_value) % probe.modulus
            : TFHEpp::modular_ntt::negate(
                  probe.a[i] * static_cast<std::uint64_t>(-key_value) %
                      probe.modulus,
                  probe.modulus);
        probe.b = TFHEpp::modular_ntt::add(
            probe.b, product, probe.modulus);
    }
    const LWECiphertext balanced_output = LWEKeySwitch(probe, balanced_key);
    if (balanced_output.b != probe.modulus / 4) return 17;
}
