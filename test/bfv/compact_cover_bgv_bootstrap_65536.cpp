#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>

#include "bfv/compact-cover-bgv-bootstrap.hpp"

int main()
{
    namespace bgv = TFHEpp::compact_cover_bgv;
    using Parameters = bgv::CompactBGV65536Parameters;
    static_assert(Parameters::ring_degree == 65536);
    static_assert(Parameters::rns_limbs == 20);
    static_assert(Parameters::gadget_digits == 23);

    std::mt19937_64 engine(UINT64_C(0x4342475642543635));
    bgv::CompactBGVKeyMaterial key;
    bgv::CompactBGVKeyGen(key, engine);
    if (!bgv::CompactBGVQuadraticHintValid(key)) {
        std::cerr << "FAIL compact BGV quadratic hint" << std::endl;
        return 1;
    }

    for (const std::uint64_t message :
         {UINT64_C(0), UINT64_C(1), Parameters::plaintext_prime - 1,
          UINT64_C(12345)}) {
        bgv::CompactBGVScalarCiphertext ciphertext;
        bgv::CompactBGVEncryptScalar(ciphertext, message, key.secret, engine);
        if (bgv::CompactBGVDecryptScalar(ciphertext, key.secret) != message) {
            std::cerr << "FAIL compact BGV scalar encryption" << std::endl;
            return 1;
        }
    }

    bgv::CompactBGVScalarCiphertext left;
    bgv::CompactBGVScalarCiphertext right;
    bgv::CompactBGVScalarCiphertext product;
    bgv::CompactBGVEncryptScalar(left, 17, key.secret, engine);
    bgv::CompactBGVEncryptScalar(right, 29, key.secret, engine);
    bgv::CompactBGVMultiply(product, left, right, key.quadratic_hint);
    if (bgv::CompactBGVDecryptScalar(product, key.secret) != 17 * 29) {
        std::cerr << "FAIL compact BGV quadratic-hint multiplication"
                  << std::endl;
        return 1;
    }

    const auto root = std::filesystem::temp_directory_path() /
                      "tfhepp_compact_bgv_bootstrap_65536";
    bgv::CompactBGVBootstrapKeyDirectoryOptions options{root, 12ULL << 30,
                                                        true};
    bgv::CompactBGVBootstrapKeyGenToDirectory(options, key, engine);
    if (!bgv::CompactBGVBootstrapKeyDirectoryComplete(root) ||
        !bgv::CompactBGVBootstrapKeyDirectoryManifestMatches(root)) {
        std::cerr << "FAIL compact BGV bootstrap directory" << std::endl;
        return 1;
    }
    const auto before =
        std::filesystem::last_write_time(bgv::CompactBGVBootstrapKeyPath(root));
    bgv::CompactBGVBootstrapKeyGenToDirectory(options, key, engine);
    const auto after =
        std::filesystem::last_write_time(bgv::CompactBGVBootstrapKeyPath(root));
    if (before != after) {
        std::cerr << "FAIL compact BGV resumable key generation" << std::endl;
        return 1;
    }

    bgv::CompactBGVBootstrapFilesystemKeyProvider provider(root);
    bgv::CompactBGVEraseMasterWitness(key);
    if (std::any_of(key.binary_ntt_witness.storage().begin(),
                    key.binary_ntt_witness.storage().end(),
                    [](const auto value) { return value != 0; })) {
        std::cerr << "FAIL compact BGV master-witness erasure" << std::endl;
        return 1;
    }
    bgv::CompactBGVScalarCiphertext input_high;
    bgv::CompactBGVScalarCiphertext input(Parameters::plaintext_prime, 1);
    bgv::CompactBGVScalarCiphertext refreshed;
    bgv::CompactBGVBootstrapTimings timings;
    bgv::CompactBGVEncryptScalar(input_high, 4242, key.secret, engine);
    bgv::CompactBGVModSwitch(input, input_high);
    const auto low_modulus = TFHEpp::modular_ntt::degree65536_primes[0].value;
    const auto injected_error =
        5 * (low_modulus / Parameters::plaintext_square);
    const auto injected_phase = TFHEpp::modular_ntt::multiply(
        Parameters::plaintext_prime, injected_error, low_modulus);
    for (auto &slot : input.value.body.slice(0, 0))
        slot = TFHEpp::modular_ntt::add(slot, injected_phase, low_modulus);
    if (bgv::CompactBGVDecryptScalar(input, key.secret) != 4242) {
        std::cerr << "FAIL compact BGV low-level modulus switch" << std::endl;
        return 1;
    }
    bgv::CompactBGVBootstrap(refreshed, input, provider, &timings);
    if (bgv::CompactBGVDecryptScalar(refreshed, key.secret) != 4242) {
        std::cerr << "FAIL compact BGV full scalar bootstrap" << std::endl;
        return 1;
    }

    // Audit the same public stages independently.  The secret is used only by
    // this test to measure centered phase errors; evaluation remains secret-free.
    bgv::CompactBGVScalarCiphertext audited_lifted(
        Parameters::plaintext_square, Parameters::rns_limbs);
    bgv::CompactBGVPhaseLiftTransition(audited_lifted, input, provider);
    const auto measured_error = [&](const auto &ciphertext) {
        return bgv::CompactBGVMaxErrorLog2(
            ciphertext, key.secret,
            bgv::CompactBGVDecryptPolynomial(ciphertext, key.secret));
    };
    if (measured_error(audited_lifted) >= 52.0) {
        std::cerr << "FAIL compact BGV phase-lift error bound" << std::endl;
        return 1;
    }
    bgv::CompactBGVScalarCiphertext audited_trace(
        Parameters::plaintext_square, Parameters::rns_limbs);
    bgv::CompactBGVTraceProjectConstant(audited_trace, audited_lifted,
                                        provider);
    if (audited_trace.limbs() != 18 ||
        measured_error(audited_trace) >= 38.0) {
        std::cerr << "FAIL compact BGV trace error bound" << std::endl;
        return 1;
    }
    const auto audited_polynomial =
        TFHEpp::digitext::GetLowestDigitRemovalPolynomialOverRange(
            Parameters::plaintext_prime, Parameters::digit_error_bound);
    std::uint64_t polynomial_hash = UINT64_C(1469598103934665603);
    for (const auto coefficient : audited_polynomial)
        for (unsigned byte = 0; byte < 8; ++byte) {
            polynomial_hash ^=
                static_cast<std::uint8_t>(coefficient >> (8 * byte));
            polynomial_hash *= UINT64_C(1099511628211);
        }
    if (polynomial_hash != Parameters::digit_polynomial_fnv1a) {
        std::cerr << "FAIL compact BGV digit-polynomial manifest"
                  << std::endl;
        return 1;
    }
    bgv::CompactBGVScalarCiphertext audited_removed(
        Parameters::plaintext_square, Parameters::rns_limbs);
    bgv::CompactBGVPolynomialEval(audited_removed, audited_polynomial,
                                  audited_trace, provider.quadraticHint());
    if (audited_removed.limbs() != 10 ||
        measured_error(audited_removed) >= 155.0) {
        std::cerr << "FAIL compact BGV digit-removal error bound" << std::endl;
        return 1;
    }
    bgv::CompactBGVScalarCiphertext audited_divided(
        Parameters::plaintext_prime, audited_removed.limbs());
    bgv::CompactBGVExactDivideByPlaintextPrime(audited_divided,
                                               audited_removed);
    if (bgv::CompactBGVDecryptScalar(audited_divided, key.secret) != 4242) {
        std::cerr << "FAIL compact BGV audited exact division" << std::endl;
        return 1;
    }
    if (refreshed.limbs() != 10) {
        std::cerr << "FAIL compact BGV certified output level" << std::endl;
        return 1;
    }
    if (timings.lift_milliseconds <= 0 ||
        timings.transition_milliseconds <= 0 ||
        timings.trace_milliseconds <= 0 ||
        timings.digit_removal_milliseconds <= 0 ||
        timings.divide_milliseconds <= 0) {
        std::cerr << "FAIL compact BGV bootstrap timings" << std::endl;
        return 1;
    }

    bgv::CompactBGVScalarCiphertext refreshed_low(Parameters::plaintext_prime,
                                                  1);
    bgv::CompactBGVModSwitch(refreshed_low, refreshed);
    bgv::CompactBGVScalarCiphertext refreshed_twice;
    bgv::CompactBGVBootstrap(refreshed_twice, refreshed_low, provider);
    if (bgv::CompactBGVDecryptScalar(refreshed_twice, key.secret) != 4242) {
        std::cerr << "FAIL compact BGV bootstrap of bootstrap" << std::endl;
        return 1;
    }

    bgv::CompactBGVScalarCiphertext refreshed_add_level(
        Parameters::plaintext_prime, 1);
    bgv::CompactBGVScalarCiphertext refreshed_twice_add_level(
        Parameters::plaintext_prime, 1);
    bgv::CompactBGVModSwitch(refreshed_add_level, refreshed);
    bgv::CompactBGVModSwitch(refreshed_twice_add_level, refreshed_twice);
    bgv::CompactBGVScalarCiphertext refreshed_sum(
        Parameters::plaintext_prime, 1);
    bgv::CompactBGVAdd(refreshed_sum, refreshed_add_level,
                       refreshed_twice_add_level);
    const auto refreshed_sum_message =
        (UINT64_C(4242) + UINT64_C(4242)) % Parameters::plaintext_prime;
    if (bgv::CompactBGVDecryptScalar(refreshed_sum, key.secret) !=
        refreshed_sum_message) {
        std::cerr << "FAIL compact BGV addition of refreshed outputs"
                  << std::endl;
        return 1;
    }
    bgv::CompactBGVScalarCiphertext bootstrapped_sum;
    bgv::CompactBGVBootstrap(bootstrapped_sum, refreshed_sum, provider);
    if (bgv::CompactBGVDecryptScalar(bootstrapped_sum, key.secret) !=
        refreshed_sum_message) {
        std::cerr << "FAIL compact BGV bootstrap after addition" << std::endl;
        return 1;
    }

    bgv::CompactBGVScalarCiphertext refreshed_mul_level(
        Parameters::plaintext_prime, 2);
    bgv::CompactBGVScalarCiphertext refreshed_twice_mul_level(
        Parameters::plaintext_prime, 2);
    bgv::CompactBGVModSwitch(refreshed_mul_level, refreshed);
    bgv::CompactBGVModSwitch(refreshed_twice_mul_level, refreshed_twice);
    const auto refreshed_product = bgv::CompactBGVMultiplyAndDrop(
        refreshed_mul_level, refreshed_twice_mul_level,
        provider.quadraticHint());
    if (refreshed_product.limbs() != 1) {
        std::cerr << "FAIL compact BGV multiplication closure level"
                  << std::endl;
        return 1;
    }
    const auto refreshed_product_message =
        (UINT64_C(4242) * UINT64_C(4242)) % Parameters::plaintext_prime;
    if (bgv::CompactBGVDecryptScalar(refreshed_product, key.secret) !=
        refreshed_product_message) {
        std::cerr << "FAIL compact BGV multiplication of refreshed outputs"
                  << std::endl;
        return 1;
    }
    bgv::CompactBGVScalarCiphertext bootstrapped_product;
    bgv::CompactBGVBootstrap(bootstrapped_product, refreshed_product, provider);
    if (bgv::CompactBGVDecryptScalar(bootstrapped_product, key.secret) !=
        refreshed_product_message) {
        std::cerr << "FAIL compact BGV bootstrap after multiplication"
                  << std::endl;
        return 1;
    }

    bool budget_guard = false;
    try {
        bgv::CompactBGVBootstrapKeyDirectoryOptions too_small{
            root / "too-small", 1, false};
        bgv::CompactBGVBootstrapKeyGenToDirectory(too_small, key, engine);
    }
    catch (const std::runtime_error &) {
        budget_guard = true;
    }
    if (!budget_guard) {
        std::cerr << "FAIL compact BGV key disk-budget guard" << std::endl;
        return 1;
    }

    {
        std::fstream manifest(bgv::CompactBGVBootstrapManifestPath(root),
                              std::ios::binary | std::ios::in | std::ios::out);
        std::uint64_t bad_magic = 0;
        manifest.write(reinterpret_cast<const char *>(&bad_magic),
                       sizeof(bad_magic));
    }
    if (bgv::CompactBGVBootstrapKeyDirectoryManifestMatches(root)) {
        std::cerr << "FAIL compact BGV corrupted manifest rejection"
                  << std::endl;
        return 1;
    }

    const auto resources = bgv::CompactBGVEstimateResources();
    std::cout << "compact BGV N=65536 scalar bootstrap PASS key="
              << static_cast<double>(resources.bootstrap_key_bytes) /
                     (1ULL << 20)
              << " MiB transition=" << timings.transition_milliseconds << " ms"
              << " trace=" << timings.trace_milliseconds << " ms"
              << " digit=" << timings.digit_removal_milliseconds << " ms"
              << std::endl;
    std::filesystem::remove_all(root);
    return 0;
}
