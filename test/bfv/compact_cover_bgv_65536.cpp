#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <random>
#include <vector>

#include "bfv/compact-cover-bgv.hpp"

int main()
{
    namespace cc = TFHEpp::compact_cover_bgv;
    namespace ntt = TFHEpp::modular_ntt;

    static_assert(cc::degree == 65536);
    static_assert(cc::max_frontier_width == 368);
    static_assert(cc::automorphism_key_count == 362);

    const auto switch_exponents = cc::thinScheduleSwitchExponents();
    cc::FrontierManifest switch_manifest{switch_exponents.size(), 1,
                                         switch_exponents};
    switch_manifest.validate();
    if (cc::generatedSubgroupSize(switch_exponents) != cc::degree ||
        cc::ThinSchedule65536::baby_product != 182 ||
        cc::ThinSchedule65536::giant_product != 181 ||
        cc::ThinSchedule65536::peak_live_ciphertexts != 368) {
        std::cerr << "FAIL concrete thin-BGV schedule" << std::endl;
        return 1;
    }

    const auto estimate = cc::estimateResources(368, 15, 4);
    const std::size_t expected_frontier_bytes =
        2ULL * 65536 * 368 * 15 * sizeof(std::uint64_t);
    if (estimate.ciphertext_bytes != expected_frontier_bytes) {
        std::cerr << "FAIL compact-cover resource estimate" << std::endl;
        return 1;
    }
    std::cout << "N=65536 width=368 fifteen-limb ciphertext "
              << static_cast<double>(estimate.ciphertext_bytes) / (1ULL << 30)
              << " GiB" << std::endl;

    for (const auto prime : ntt::degree65536_primes) {
        if ((prime.value - 1) % (2 * cc::degree) != 0 ||
            (prime.value - 1) % (cc::ThinSchedule65536::plaintext_prime *
                                 static_cast<std::uint64_t>(
                                     cc::ThinSchedule65536::plaintext_prime)) !=
                0) {
            std::cerr << "FAIL large-degree NTT prime" << std::endl;
            return 1;
        }
    }

    // Real N=65536 negacyclic transform round trip.
    cc::FrontierElement transform_value(1, 1);
    auto transform_slice = transform_value.slice(0, 0);
    std::mt19937_64 engine(0x434f4d50414354ULL);
    const auto modulus = ntt::degree65536_primes[0].value;
    for (auto &value : transform_slice) value = engine() % modulus;
    const std::vector<std::uint64_t> original(transform_slice.begin(),
                                              transform_slice.end());
    cc::forwardNTT(transform_value);
    cc::inverseNTT(transform_value);
    if (!std::equal(original.begin(), original.end(),
                    transform_slice.begin())) {
        std::cerr << "FAIL N=65536 HEXL NTT round trip" << std::endl;
        return 1;
    }

    // Match the Lean backend definition: natural slot k is evaluation at
    // psi^(2k+1). A sparse polynomial makes this check linear-time while still
    // detecting backend-order and twist-direction mismatches.
    std::fill(transform_slice.begin(), transform_slice.end(), 0);
    transform_slice[0] = 7;
    transform_slice[1] = 11;
    transform_slice[17] = 3;
    transform_slice[cc::degree - 1] = 5;
    cc::forwardNTT(transform_value);
    const auto generator = ntt::degree65536_primes[0].primitive_root;
    const auto psi = ntt::power(
        generator, (modulus - 1) / (2 * cc::degree), modulus);
    for (const std::size_t slot :
         {std::size_t{0}, std::size_t{1}, std::size_t{17},
          std::size_t{32768}, std::size_t{65535}}) {
        const auto point = ntt::power(psi, 2 * slot + 1, modulus);
        auto expected = UINT64_C(7);
        expected = ntt::add(expected, ntt::multiply(11, point, modulus),
                            modulus);
        expected = ntt::add(
            expected,
            ntt::multiply(3, ntt::power(point, 17, modulus), modulus),
            modulus);
        expected = ntt::add(
            expected,
            ntt::multiply(5, ntt::power(point, cc::degree - 1, modulus),
                          modulus),
            modulus);
        if (transform_slice[slot] != expected) {
            std::cerr << "FAIL mathematical NTT ordering at slot " << slot
                      << std::endl;
            return 1;
        }
    }

    // Materialize and relabel the actual width-368, one-limb frontier.  The
    // fifteen-limb profile is opt-in because it consumes 5.39 GiB per
    // ciphertext.
    cc::FrontierElement source(368, 1);
    cc::FrontierElement target(368, 1);
    for (std::size_t label = 0; label < source.width(); ++label) {
        auto values = source.slice(0, label);
        for (std::size_t slot = 0; slot < cc::degree; ++slot)
            values[slot] = (label * 17 + slot * 3 + 1) % modulus;
    }
    std::vector<cc::RelabelEntry> mapping(368);
    for (std::size_t label = 0; label < mapping.size(); ++label)
        mapping[label] = {367 - label, label % 2 == 0 ? 1U : 3U};
    cc::relabel(target, source, mapping);
    for (const std::size_t label :
         {std::size_t{0}, std::size_t{1}, std::size_t{183}, std::size_t{367}}) {
        const auto input = source.slice(0, mapping[label].source_index);
        const auto output = target.slice(0, label);
        for (const std::size_t slot :
             {std::size_t{0}, std::size_t{1}, std::size_t{32768},
              std::size_t{65535}}) {
            if (output[slot] != input[cc::automorphismSlot(
                                    mapping[label].automorphism, slot)]) {
                std::cerr << "FAIL width-368 relabel" << std::endl;
                return 1;
            }
        }
    }

    // Small streamed-key fixture using the production N=65536 slice size.
    const auto root = std::filesystem::temp_directory_path() /
                      "tfhepp_compact_cover_bgv_65536";
    std::filesystem::create_directories(root);
    const auto key_path = root / "transition.key";
    constexpr std::size_t rows = 2;
    constexpr std::size_t width = 2;
    {
        cc::SeededTransitionKeyWriter writer(key_path, width, 1, rows,
                                             16ULL << 20);
        for (std::size_t row = 0; row < rows; ++row)
            for (std::size_t label = 0; label < width; ++label)
                writer.writeSeed({0x1000 + row * width + label, 2, 3, 4});
        std::vector<std::uint64_t> body(cc::degree);
        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t label = 0; label < width; ++label) {
                for (std::size_t slot = 0; slot < body.size(); ++slot)
                    body[slot] = (row + 3 * label + slot) % modulus;
                writer.writeBodySlice(body);
            }
        }
    }
    cc::SeededTransitionKeyProvider provider(key_path);
    std::vector<std::uint64_t> mask(cc::degree), body(cc::degree);
    provider.expandMaskSlice(1, 0, 1, mask);
    provider.loadBodySlice(1, 0, 1, body);
    if (mask[0] >= modulus || body[12345] != (1 + 3 + 12345) % modulus) {
        std::cerr << "FAIL streamed transition-key slice" << std::endl;
        return 1;
    }

    // Exact phase test for a streamed limb-local gadget transition.  The
    // fixture is noiseless so equality can be checked in every NTT slot.
    const auto algebra_key_path = root / "transition-algebra.key";
    constexpr unsigned gadget_bits = 16;
    const std::size_t gadget_rows = cc::gadgetRows(modulus, gadget_bits);
    std::vector<std::uint64_t> source_secret(cc::degree);
    std::vector<std::uint64_t> target_secret(cc::degree);
    for (std::size_t slot = 0; slot < cc::degree; ++slot) {
        source_secret[slot] = engine() % modulus;
        target_secret[slot] = engine() % modulus;
    }
    std::vector<cc::TransitionSeed> algebra_seeds(gadget_rows);
    for (std::size_t row = 0; row < gadget_rows; ++row)
        algebra_seeds[row] = {UINT64_C(0xabc000) + row, 11, 12, 13};
    {
        cc::SeededTransitionKeyWriter writer(algebra_key_path, 1, 1,
                                             gadget_rows, 8ULL << 20);
        for (const auto seed : algebra_seeds) writer.writeSeed(seed);
        std::vector<std::uint64_t> row_mask(cc::degree);
        std::vector<std::uint64_t> row_body(cc::degree);
        for (std::size_t row = 0; row < gadget_rows; ++row) {
            const auto gadget = (UINT64_C(1) << (row * gadget_bits)) % modulus;
            cc::expandSeededUniformSlice(algebra_seeds[row], modulus, row_mask);
            for (std::size_t slot = 0; slot < cc::degree; ++slot) {
                row_body[slot] = ntt::subtract(
                    ntt::multiply(row_mask[slot], target_secret[slot], modulus),
                    ntt::multiply(gadget, source_secret[slot], modulus),
                    modulus);
            }
            writer.writeBodySlice(row_body);
        }
    }
    cc::SeededTransitionKeyProvider algebra_key(algebra_key_path);
    cc::FrontierCiphertext transition_input(1, 1);
    cc::FrontierCiphertext transition_output(1, 1);
    auto input_mask = transition_input.mask.slice(0, 0);
    auto input_body = transition_input.body.slice(0, 0);
    for (auto &coefficient : input_mask) coefficient = engine() % modulus;
    cc::forwardNTT(transition_input.mask);
    std::vector<std::uint64_t> message(cc::degree);
    for (std::size_t slot = 0; slot < cc::degree; ++slot) {
        message[slot] = engine() % modulus;
        input_body[slot] = ntt::add(
            ntt::multiply(input_mask[slot], source_secret[slot], modulus),
            message[slot], modulus);
    }
    cc::applyTransition(transition_output, transition_input, algebra_key,
                        gadget_bits);
    const auto switched_mask = transition_output.mask.slice(0, 0);
    const auto switched_body = transition_output.body.slice(0, 0);
    for (std::size_t slot = 0; slot < cc::degree; ++slot) {
        const auto phase = ntt::subtract(
            switched_body[slot],
            ntt::multiply(switched_mask[slot], target_secret[slot], modulus),
            modulus);
        if (phase != message[slot]) {
            std::cerr << "FAIL streamed transition phase at slot " << slot
                      << std::endl;
            return 1;
        }
    }

    // Certified full-CRT decomposition with stored-uniform masks.  Two limbs
    // are enough to detect the per-limb decomposition bug: input coefficients
    // are sampled independently in both residues and reconstructed modulo the
    // complete product before taking each balanced digit.
    const auto uniform_path = root / "uniform-full-modulus.key";
    constexpr std::size_t uniform_limbs = 2;
    constexpr std::size_t uniform_rows = 5;
    cc::FullModulusGadget full_gadget(uniform_limbs, uniform_rows);
    std::vector<std::vector<std::uint64_t>> full_source_secret(
        uniform_limbs, std::vector<std::uint64_t>(cc::degree));
    std::vector<std::vector<std::uint64_t>> full_target_secret(
        uniform_limbs, std::vector<std::uint64_t>(cc::degree));
    {
        cc::UniformTransitionKeyWriter writer(uniform_path, 1, uniform_limbs,
                                              uniform_rows, 32ULL << 20);
        std::vector<std::uint64_t> row_mask(cc::degree), row_body(cc::degree);
        for (std::size_t row = 0; row < uniform_rows; ++row) {
            for (std::size_t limb = 0; limb < uniform_limbs; ++limb) {
                const auto prime = ntt::degree65536_primes[limb].value;
                if (row == 0) {
                    for (std::size_t slot = 0; slot < cc::degree; ++slot) {
                        full_source_secret[limb][slot] = engine() % prime;
                        full_target_secret[limb][slot] = engine() % prime;
                    }
                }
                boost::multiprecision::cpp_int gadget_power = 1;
                for (std::size_t i = 0; i < row; ++i)
                    gadget_power *= full_gadget.base();
                const auto gadget_residue =
                    cc::FullModulusGadget::residue(gadget_power, prime);
                for (std::size_t slot = 0; slot < cc::degree; ++slot) {
                    row_mask[slot] = engine() % prime;
                    row_body[slot] = ntt::subtract(
                        ntt::multiply(row_mask[slot],
                                      full_target_secret[limb][slot], prime),
                        ntt::multiply(gadget_residue,
                                      full_source_secret[limb][slot], prime),
                        prime);
                }
                writer.writeSlice(row_mask, row_body);
            }
        }
        writer.finish();
    }
    cc::UniformTransitionKeyProvider uniform_key(uniform_path);
    cc::FrontierCiphertext full_input(1, uniform_limbs);
    cc::FrontierCiphertext full_output(1, uniform_limbs);
    std::vector<std::vector<std::uint64_t>> full_message(
        uniform_limbs, std::vector<std::uint64_t>(cc::degree));
    for (std::size_t limb = 0; limb < uniform_limbs; ++limb) {
        const auto prime = ntt::degree65536_primes[limb].value;
        auto mask_slice = full_input.mask.slice(limb, 0);
        for (auto &coefficient : mask_slice) coefficient = engine() % prime;
    }
    cc::forwardNTT(full_input.mask);
    for (std::size_t limb = 0; limb < uniform_limbs; ++limb) {
        const auto prime = ntt::degree65536_primes[limb].value;
        const auto mask_slice = full_input.mask.slice(limb, 0);
        auto body_slice = full_input.body.slice(limb, 0);
        for (std::size_t slot = 0; slot < cc::degree; ++slot) {
            full_message[limb][slot] = engine() % prime;
            body_slice[slot] =
                ntt::add(ntt::multiply(mask_slice[slot],
                                       full_source_secret[limb][slot], prime),
                         full_message[limb][slot], prime);
        }
    }
    cc::applyFullModulusTransition(full_output, full_input, uniform_key);
    for (std::size_t limb = 0; limb < uniform_limbs; ++limb) {
        const auto prime = ntt::degree65536_primes[limb].value;
        const auto mask_slice = full_output.mask.slice(limb, 0);
        const auto body_slice = full_output.body.slice(limb, 0);
        for (std::size_t slot = 0; slot < cc::degree; ++slot) {
            const auto phase = ntt::subtract(
                body_slice[slot],
                ntt::multiply(mask_slice[slot], full_target_secret[limb][slot],
                              prime),
                prime);
            if (phase != full_message[limb][slot]) {
                std::cerr << "FAIL full-modulus transition phase at limb "
                          << limb << " slot " << slot << std::endl;
                return 1;
            }
        }
    }
    std::filesystem::remove_all(root);

    bool guard_triggered = false;
    try {
        cc::SeededTransitionKeyWriter impossible(
            std::filesystem::temp_directory_path() / "too_large.key", 368, 15,
            4, 1);
    }
    catch (const std::runtime_error &) {
        guard_triggered = true;
    }
    if (!guard_triggered) {
        std::cerr << "FAIL transition-key resource guard" << std::endl;
        return 1;
    }

    std::cout << "compact-cover N=65536 substrate PASS" << std::endl;
    return 0;
}
