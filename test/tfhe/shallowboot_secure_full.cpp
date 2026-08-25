#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <string_view>
#include <variant>
#include <vector>

#include <tfhe/shallowboot_lowdepth.hpp>

int main(const int argc, const char *const argv[])
{
    using namespace TFHEpp::shallowboot::lowdepth;
    using clock = std::chrono::steady_clock;
    constexpr std::size_t degree = 8192;
    const bool seeded = argc > 1 && std::string_view(argv[1]) == "seeded";
    const bool structured = seeded ||
        (argc > 1 && std::string_view(argv[1]) == "structured");
    const std::uint32_t source_n = structured ? 1024 : 1450;
    const std::uint32_t source_h = structured ? 64 : 37;
    constexpr std::uint32_t output_n = 1450;
    constexpr std::uint32_t lwe_modulus = 512;
    const std::array<RNSRing, 4> rings = {
        RNSRing(degree, secure_q151_primes),
        RNSRing(degree, secure_q69_primes),
        RNSRing(degree, secure_q52_primes),
        RNSRing(degree, secure_q36_primes)};
    std::vector<std::int64_t> ring_key(degree);
    for (std::size_t i = 0; i < degree; i++)
        ring_key[i] = static_cast<std::int64_t>(i % 3) - 1;
    const std::array<RNSSecret, 4> secrets = {
        QuadraticHintSmallSecretFromCoefficients(rings[0], ring_key),
        QuadraticHintSmallSecretFromCoefficients(rings[1], ring_key),
        QuadraticHintSmallSecretFromCoefficients(rings[2], ring_key),
        QuadraticHintSmallSecretFromCoefficients(rings[3], ring_key)};

    std::vector<std::uint8_t> source_key(source_n, 0);
    PBCSchedule schedule;
    if (structured) {
        schedule.source_dimension = source_n;
        schedule.bucket_count = source_h;
        schedule.copies = 1;
        schedule.bucket_indices.resize(source_h);
        schedule.selected_source.resize(source_h);
        const std::uint32_t width = source_n / source_h;
        for (std::uint32_t block = 0; block < source_h; block++) {
            for (std::uint32_t offset = 0; offset < width; offset++)
                schedule.bucket_indices[block].push_back(block * width + offset);
            const std::uint32_t selected = block * width + 3;
            schedule.selected_source[block] = selected;
            source_key[selected] = 1;
        }
    }
    else {
        std::vector<std::uint32_t> support;
        for (std::uint32_t i = 0; i < source_h; i++) {
            const std::uint32_t index = (47 * i + 11) % source_n;
            source_key[index] = 1;
            support.push_back(index);
        }
        schedule = BuildPBCScheduleWithRetry(
            source_n, source_h + 3, 3, support, 20261730);
    }
    std::vector<std::uint8_t> output_lwe_key(output_n, 0);
    for (std::uint32_t i = 0; i < 60; i++)
        output_lwe_key[(53 * i + 7) % output_n] = 1;
    std::mt19937_64 engine(20261730);
    const auto keygen_start = clock::now();
#ifdef USE_BLAKE3
    using PBCKeyVariant = std::variant<RNSPBCBootstrappingKey,
                                       SeededRNSPBCBootstrappingKey>;
    PBCKeyVariant bootstrapping_key = seeded
        ? PBCKeyVariant(PBCBootstrappingKeyGenSeededLSB(
              rings[0], secrets[0], schedule, 4, 0.75, engine))
        : PBCKeyVariant(PBCBootstrappingKeyGenLSB(
              rings[0], secrets[0], schedule, 4, 0.75, engine));
#else
    const RNSPBCBootstrappingKey bootstrapping_key = PBCBootstrappingKeyGenLSB(
        rings[0], secrets[0], schedule, 4, 0.75, engine);
#endif
    const LWEKeySwitchKey output_key = LWEKeySwitchKeyGen(
        ring_key, output_lwe_key, 1U << 15, 1, 3.2, engine, 15);
    const auto keygen_end = clock::now();

    std::vector<std::uint32_t> mask(source_n);
    std::uniform_int_distribution<std::uint32_t> uniform(0, lwe_modulus - 1);
    std::uint32_t dot = 0;
    for (std::size_t i = 0; i < mask.size(); i++) {
        mask[i] = uniform(engine);
        dot = (dot + mask[i] * source_key[i]) % lwe_modulus;
    }
    // A representative rounded source error at the specified sigma=3.2.
    constexpr std::uint32_t source_error = 3;
    const LWEPhase input = {
        mask, static_cast<std::uint32_t>(
                  (dot + lwe_modulus / 4 + source_error) % lwe_modulus),
        lwe_modulus};
    const LWEPhase opposite_input = {
        mask, static_cast<std::uint32_t>(
                  (dot + 3 * lwe_modulus / 4 + source_error) % lwe_modulus),
        lwe_modulus};
    const auto lut = Algorithm3BinarySignPlaintextLUT(rings.back(),
                                                       lwe_modulus);
    constexpr std::array<std::uint32_t, 4> windows = {8, 2, 2, 2};
    MultiStageStats stats;
    const auto bootstrap = [&](const LWEPhase &value,
                               MultiStageStats *statistics = nullptr) {
#ifdef USE_BLAKE3
        return std::visit(
            [&](const auto &key) {
                return Algorithm3PBCQHMultiStageBootstrap(
                    rings, secrets, schedule, key, value, lut, windows,
                    output_key, lwe_modulus, 1U << 15, 4, statistics);
            },
            bootstrapping_key);
#else
        return Algorithm3PBCQHMultiStageBootstrap(
            rings, secrets, schedule, bootstrapping_key, value, lut, windows,
            output_key, lwe_modulus, 1U << 15, 4, statistics);
#endif
    };
    const auto bootstrap_start = clock::now();
    const LWEPhase output = bootstrap(input, &stats);
    const LWEPhase opposite_output = bootstrap(opposite_input);
    const auto bootstrap_end = clock::now();
    const std::uint32_t phase = LWEDecryptPhase(output, output_lwe_key);
    if (phase < 320 || phase > 448) {
        std::cerr << "unexpected output phase=" << phase << '\n';
        return 1;
    }
    const std::uint32_t opposite_phase =
        LWEDecryptPhase(opposite_output, output_lwe_key);
    if (opposite_phase < 64 || opposite_phase > 192) {
        std::cerr << "unexpected opposite output phase=" << opposite_phase
                  << '\n';
        return 1;
    }
    const double keygen_seconds =
        std::chrono::duration<double>(keygen_end - keygen_start).count();
    const double bootstrap_ms =
        std::chrono::duration<double, std::milli>(bootstrap_end -
                                                  bootstrap_start).count() /
        2;
    std::size_t entry_count = schedule.bucket_count;
    for (const auto &bucket : schedule.bucket_indices)
        entry_count += bucket.size();
    const double words_per_entry = seeded
        ? static_cast<double>(rings[0].levels() * degree) + 4
        : static_cast<double>(rings[0].levels() * 2 * degree);
    const double pbc_key_gib = static_cast<double>(entry_count) *
        words_per_entry * sizeof(std::uint64_t) /
        (1024.0 * 1024.0 * 1024.0);
    std::cout << "profile="
              << (seeded ? "structured-seeded"
                         : structured ? "structured" : "general")
              << '\n'
              << "pbc_entries=" << entry_count << '\n'
              << "estimated_pbc_key_gib=" << pbc_key_gib << '\n'
              << "keygen_seconds=" << keygen_seconds << '\n'
              << "bootstrap_ms=" << bootstrap_ms << '\n'
              << "output_phase=" << phase << '\n'
              << "opposite_output_phase=" << opposite_phase << '\n';
}
