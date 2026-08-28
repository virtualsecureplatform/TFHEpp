#include <cstdint>
#include <iostream>
#include <random>

#include "bfv/regular-cover-bgv.hpp"

int main()
{
    // X^2+1 splits at ±psi modulo 65537.  Its Galois group is C2 and the
    // nontrivial automorphism swaps the two split coordinates.
    using Scheme = TFHEpp::regular_cover_bgv::Scheme<2, 2, 65537, 2>;
    std::mt19937_64 engine(0x5243424756ULL);

    const std::uint64_t psi = Scheme::pow(3, (65537 - 1) / 4);
    if (Scheme::mul(psi, psi) != 65536) {
        std::cerr << "FAIL degree-two cyclotomic root" << std::endl;
        return 1;
    }
    Scheme::Element operational_secret{};
    std::uniform_int_distribution<int> coefficient(-1, 1);
    for (auto &component : operational_secret) {
        const std::uint64_t c0 = Scheme::fromSigned(coefficient(engine));
        const std::uint64_t c1 = Scheme::fromSigned(coefficient(engine));
        component[0] = Scheme::add(c0, Scheme::mul(c1, psi));
        component[1] = Scheme::sub(c0, Scheme::mul(c1, psi));
    }
    const auto key =
        Scheme::keygenWithOperationalSecret(engine, operational_secret);
    if (!Scheme::checkQuadraticHint(key)) {
        std::cerr << "FAIL quadratic hint" << std::endl;
        return 1;
    }

    // Completion must produce a binary witness in every cover/NTT coordinate.
    for (const auto &component : key.completion.witness)
        for (const auto value : component)
            if (value != 0 && value != 1) {
                std::cerr << "FAIL non-binary completed witness" << std::endl;
                return 1;
            }

    // Base-2 decomposition uses at most 17 rows for q=65537.  With row errors
    // in [-1,1] and input error bound 2, the homomorphic phase stays in the
    // interpolation interval [-20,20].
    const auto digit_polynomial = Scheme::interpolationPolynomial(20);
    const auto boot_key =
        Scheme::makeBootKey(key.operational_secret, engine, 2, 1);

    bool saw_nonzero_boot_error = false;
    std::uint64_t gadget = 1;
    for (const auto &row : boot_key.rows) {
        const auto row_phase = Scheme::phase(row, key.operational_secret);
        const auto message = Scheme::scalarMul(key.operational_secret, gadget);
        const auto error = Scheme::elementSub(row_phase, message);
        for (const auto &component : error) {
            for (const auto value : component) {
                const std::int64_t centered = Scheme::centered(value);
                if (centered % 2 != 0 || centered < -2 || centered > 2) {
                    std::cerr << "FAIL boot-key error law" << std::endl;
                    return 1;
                }
                saw_nonzero_boot_error |= centered != 0;
            }
        }
        gadget *= 2;
    }
    if (!saw_nonzero_boot_error) {
        std::cerr << "FAIL boot key unexpectedly noiseless" << std::endl;
        return 1;
    }

    for (std::uint64_t message = 0; message < 2; ++message) {
        for (int trial = 0; trial < 8; ++trial) {
            auto ciphertext = Scheme::encryptScalar(
                message, key.operational_secret, engine, 2);
            if (Scheme::decryptScalar(ciphertext, key.operational_secret) !=
                message) {
                std::cerr << "FAIL initial decrypt" << std::endl;
                return 1;
            }

            // Reuse the same noisy boot key twice under the same operational
            // secret.
            ciphertext =
                Scheme::refresh(ciphertext, boot_key, key.operational_secret,
                                key.hint, digit_polynomial, engine, 1);
            if (Scheme::decryptScalar(ciphertext, key.operational_secret) !=
                message) {
                std::cerr << "FAIL first refresh" << std::endl;
                return 1;
            }
            ciphertext =
                Scheme::refresh(ciphertext, boot_key, key.operational_secret,
                                key.hint, digit_polynomial, engine, 1);
            if (Scheme::decryptScalar(ciphertext, key.operational_secret) !=
                message) {
                std::cerr << "FAIL second refresh" << std::endl;
                return 1;
            }
        }
    }

    // Quadratic-hint multiplication retains a two-component ciphertext.
    const auto one =
        Scheme::encryptScalar(1, key.operational_secret, engine, 0);
    const auto product = Scheme::quadraticHintMul(one, one, key.hint);
    if (Scheme::decryptScalar(product, key.operational_secret) != 1) {
        std::cerr << "FAIL quadratic-hint multiply" << std::endl;
        return 1;
    }

    // The lifted action preserves diagonal scalar plaintexts.
    const auto rotated = Scheme::liftedAutomorphism(one, 1);
    if (Scheme::decryptScalar(rotated, Scheme::liftedAutomorphism(
                                           key.operational_secret, 1)) != 1) {
        std::cerr << "FAIL lifted automorphism" << std::endl;
        return 1;
    }

    // The actual same-key lifted automorphism is followed by a gadget switch
    // back from sigma(z) to z.
    const auto switch_key =
        Scheme::makeSwitchKey(key.operational_secret, 1, 16, engine);
    const auto switched = Scheme::applySwitchKey(one, switch_key);
    if (Scheme::decryptScalar(switched, key.operational_secret) != 1) {
        std::cerr << "FAIL lifted automorphism key switch" << std::endl;
        return 1;
    }

    // Exercise the exact ordinary-base-RLWE -> invariant-cover -> simultaneous
    // automorphism compiler used by the security reduction.
    std::array<Scheme::BaseCiphertext, Scheme::group_size> source_rows{};
    std::array<Scheme::Base, Scheme::group_size> source_errors{};
    const Scheme::Element invariant = Scheme::elementMul(
        key.completion.multiplier,
        Scheme::elementSub(key.completion.witness, key.completion.offset));
    const Scheme::Base bottom_secret = invariant[0];
    std::uniform_int_distribution<std::uint64_t> uniform(0, 65536);
    for (std::size_t h = 0; h < Scheme::group_size; ++h) {
        for (std::size_t j = 0; j < Scheme::split_degree; ++j) {
            source_rows[h].mask[j] = uniform(engine);
            source_errors[h][j] =
                Scheme::fromSigned(static_cast<std::int64_t>((j + h) % 3) - 1);
            source_rows[h].body[j] = Scheme::add(
                Scheme::mul(source_rows[h].mask[j], bottom_secret[j]),
                source_errors[h][j]);
        }
    }
    const auto source_cover = Scheme::assembleCoverSamples(source_rows);
    std::array<Scheme::Element, Scheme::group_size> weights{};
    weights[1] = Scheme::constant(7);
    const auto compiled =
        Scheme::compileWitnessRow(source_cover, key.completion, weights);
    const auto compiled_phase = Scheme::phase(compiled, key.completion.witness);
    const auto assembled_error = Scheme::assembleCover(source_errors);
    const auto expected_phase = Scheme::elementAdd(
        assembled_error,
        Scheme::elementMul(
            weights[1], Scheme::liftedAutomorphism(key.completion.witness, 1)));
    if (compiled_phase != expected_phase) {
        std::cerr << "FAIL regular-cover security compiler" << std::endl;
        return 1;
    }

    std::cout << "regular-cover BGV prototype PASS" << std::endl;
    std::cout << "  R=F_q[X]/(X^2+1), |Gamma|=2, q=65537, t=2" << std::endl;
    std::cout << "  same boot key reused for two refreshes" << std::endl;
    std::cout << "  lifted automorphism switch and security compiler checked"
              << std::endl;
    return 0;
}
