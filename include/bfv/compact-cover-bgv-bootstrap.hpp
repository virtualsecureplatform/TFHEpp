#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "bfv/bfv-digitext.hpp"
#include "bfv/compact-cover-bgv.hpp"

namespace TFHEpp::compact_cover_bgv {

struct CompactBGV65536Parameters {
    static constexpr std::size_t ring_degree = degree;
    static constexpr std::size_t rns_limbs = 20;
    static constexpr std::size_t gadget_digits = 23;
    static constexpr std::size_t phase_lift_gadget_digits = 2;
    static constexpr std::size_t trace_key_count = 16;
    static constexpr std::uint64_t digit_error_bound = 23;
    static constexpr std::size_t secret_weight = 32;
    static constexpr unsigned cbd_eta = 20;
    static constexpr std::uint64_t plaintext_prime = 65537;
    static constexpr std::uint64_t plaintext_square =
        plaintext_prime * plaintext_prime;
    static constexpr std::uint32_t certificate_version = 8;
    static constexpr std::uint64_t bootstrap_manifest_magic =
        UINT64_C(0x4342475642543031);
    static constexpr const char *scalar_gate_manifest_sha256 =
        "9584b90e526fc67ca85c4ea1b6cea004ca64b30b70e4e6609d0961c7e6144843";
    static constexpr const char *certificate_sha256 =
        "155a50174b509e170a649674d38c55c450f475e43ae56d32a8865c3459eaab7e";
    static constexpr std::uint64_t digit_polynomial_fnv1a =
        UINT64_C(0xa4beb97e03c045ab);
    static constexpr const char *general_gate_manifest_sha256 =
        "6d376c090ad655badcb23b00da1e16c4429c83d3f20c3195df93eef5bf94049e";
};

struct CompactBGVSecretKey {
    std::vector<std::int8_t> coefficients;
    FrontierElement ntt;

    CompactBGVSecretKey()
        : coefficients(degree), ntt(1, CompactBGV65536Parameters::rns_limbs)
    {
    }
};

struct CompactBGVQuadraticHint {
    FrontierElement alpha;
    FrontierElement beta;
    FrontierElement u;
    FrontierElement v;

    CompactBGVQuadraticHint()
        : alpha(1, CompactBGV65536Parameters::rns_limbs),
          beta(1, CompactBGV65536Parameters::rns_limbs),
          u(1, CompactBGV65536Parameters::rns_limbs),
          v(1, CompactBGV65536Parameters::rns_limbs)
    {
    }
};

struct CompactBGVKeyMaterial {
    CompactBGVSecretKey secret;
    FrontierElement binary_ntt_witness;
    CompactBGVQuadraticHint quadratic_hint;

    CompactBGVKeyMaterial()
        : binary_ntt_witness(1, CompactBGV65536Parameters::rns_limbs)
    {
    }
};

struct CompactBGVScalarCiphertext {
    FrontierCiphertext value;
    std::uint64_t plaintext_modulus{};
    double log2_variance{};

    explicit CompactBGVScalarCiphertext(
        const std::uint64_t plaintext =
            CompactBGV65536Parameters::plaintext_prime,
        const std::size_t limbs = CompactBGV65536Parameters::rns_limbs)
        : value(1, limbs), plaintext_modulus(plaintext)
    {
    }

    std::size_t limbs() const { return value.mask.limbs(); }
};

struct CompactBGVBootstrapTimings {
    double lift_milliseconds{};
    double transition_milliseconds{};
    double trace_milliseconds{};
    double digit_removal_milliseconds{};
    double divide_milliseconds{};
};

struct CompactBGVResourceEstimate {
    std::size_t ciphertext_bytes{};
    std::size_t bootstrap_key_bytes{};
    std::size_t transition_scratch_bytes{};
};

inline std::size_t CompactBGVQuadraticHintBytes();
inline std::size_t CompactBGVTraceKeyLimbs(std::size_t index);

inline CompactBGVResourceEstimate CompactBGVEstimateResources()
{
    const auto phase_estimate =
        estimateResources(1, CompactBGV65536Parameters::rns_limbs,
                          CompactBGV65536Parameters::phase_lift_gadget_digits);
    std::size_t key_bytes = sizeof(TransitionKeyHeader) +
                            phase_estimate.uniform_key_bytes +
                            CompactBGVQuadraticHintBytes();
    for (std::size_t index = 0;
         index < CompactBGV65536Parameters::trace_key_count; ++index) {
        const auto estimate =
            estimateResources(1, CompactBGVTraceKeyLimbs(index),
                              CompactBGV65536Parameters::gadget_digits);
        key_bytes += sizeof(TransitionKeyHeader) + estimate.uniform_key_bytes;
    }
    return {phase_estimate.ciphertext_bytes, key_bytes,
            checkedProduct({8, CompactBGV65536Parameters::rns_limbs, degree,
                            sizeof(std::uint64_t)})};
}

template <class Engine>
inline std::uint64_t sampleUniformMod(Engine &engine,
                                      const std::uint64_t modulus)
{
    const std::uint64_t threshold = -modulus % modulus;
    std::uint64_t value;
    do value = static_cast<std::uint64_t>(engine());
    while (value < threshold);
    return value % modulus;
}

template <class Engine>
inline std::int64_t sampleCBD20(Engine &engine)
{
    unsigned positive = 0;
    unsigned negative = 0;
    for (unsigned bit = 0; bit < CompactBGV65536Parameters::cbd_eta; ++bit) {
        positive += static_cast<unsigned>(engine() & 1U);
        negative += static_cast<unsigned>(engine() & 1U);
    }
    return static_cast<std::int64_t>(positive) -
           static_cast<std::int64_t>(negative);
}

inline std::uint64_t signedResidue(const std::int64_t value,
                                   const std::uint64_t modulus)
{
    if (value >= 0) return static_cast<std::uint64_t>(value) % modulus;
    const auto magnitude = static_cast<std::uint64_t>(-(value + 1)) + 1;
    return modular_ntt::negate(magnitude % modulus, modulus);
}

template <class Engine>
inline void CompactBGVKeyGen(CompactBGVKeyMaterial &key, Engine &engine)
{
    std::vector<std::size_t> support(degree);
    std::iota(support.begin(), support.end(), 0);
    std::shuffle(support.begin(), support.end(), engine);
    std::fill(key.secret.coefficients.begin(), key.secret.coefficients.end(),
              0);
    for (std::size_t index = 0;
         index < CompactBGV65536Parameters::secret_weight; ++index)
        key.secret.coefficients[support[index]] = (engine() & 1U) == 0 ? -1 : 1;

    for (std::size_t limb = 0; limb < CompactBGV65536Parameters::rns_limbs;
         ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        auto secret = key.secret.ntt.slice(limb, 0);
        for (std::size_t coefficient = 0; coefficient < degree; ++coefficient)
            secret[coefficient] =
                signedResidue(key.secret.coefficients[coefficient], modulus);
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        plan.forward(secret);

        auto witness = key.binary_ntt_witness.slice(limb, 0);
        auto alpha = key.quadratic_hint.alpha.slice(limb, 0);
        auto beta = key.quadratic_hint.beta.slice(limb, 0);
        auto u = key.quadratic_hint.u.slice(limb, 0);
        auto v = key.quadratic_hint.v.slice(limb, 0);
        for (std::size_t slot = 0; slot < degree; ++slot) {
            if (limb == 0)
                witness[slot] = engine() & 1U;
            else
                witness[slot] = key.binary_ntt_witness.slice(0, 0)[slot];
            do alpha[slot] = sampleUniformMod(engine, modulus);
            while (alpha[slot] == 0);
            beta[slot] = modular_ntt::add(
                modular_ntt::multiply(alpha[slot], witness[slot], modulus),
                secret[slot], modulus);
            u[slot] = modular_ntt::subtract(
                modular_ntt::add(beta[slot], beta[slot], modulus), alpha[slot],
                modulus);
            v[slot] = modular_ntt::subtract(
                modular_ntt::multiply(alpha[slot], beta[slot], modulus),
                modular_ntt::multiply(beta[slot], beta[slot], modulus),
                modulus);
        }
    }
}

// Production entry point. The engine-taking overload remains available only
// for deterministic conformance tests and explicitly managed deployments.
inline void CompactBGVKeyGen(CompactBGVKeyMaterial &key)
{
    CompactBGVKeyGen(key, TFHEpp::generator);
}

inline bool CompactBGVQuadraticHintValid(const CompactBGVKeyMaterial &key)
{
    for (std::size_t limb = 0; limb < CompactBGV65536Parameters::rns_limbs;
         ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto secret = key.secret.ntt.slice(limb, 0);
        const auto u = key.quadratic_hint.u.slice(limb, 0);
        const auto v = key.quadratic_hint.v.slice(limb, 0);
        for (std::size_t slot = 0; slot < degree; ++slot)
            if (modular_ntt::multiply(secret[slot], secret[slot], modulus) !=
                modular_ntt::add(
                    modular_ntt::multiply(u[slot], secret[slot], modulus),
                    v[slot], modulus))
                return false;
    }
    return true;
}

inline void CompactBGVEraseMasterWitness(CompactBGVKeyMaterial &key)
{
    std::fill(key.binary_ntt_witness.storage().begin(),
              key.binary_ntt_witness.storage().end(), 0);
}

template <class Engine>
inline void CompactBGVEncryptScalar(CompactBGVScalarCiphertext &result,
                                    const std::uint64_t message,
                                    const CompactBGVSecretKey &secret,
                                    Engine &engine)
{
    if (result.plaintext_modulus !=
            CompactBGV65536Parameters::plaintext_prime &&
        result.plaintext_modulus != CompactBGV65536Parameters::plaintext_square)
        throw std::invalid_argument("unsupported compact BGV plaintext");
    std::vector<std::int64_t> error(degree);
    for (auto &coefficient : error) coefficient = sampleCBD20(engine);
    if (result.limbs() > secret.ntt.limbs())
        throw std::invalid_argument("compact BGV secret has too few limbs");
    for (std::size_t limb = 0; limb < result.limbs(); ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        auto mask = result.value.mask.slice(limb, 0);
        auto body = result.value.body.slice(limb, 0);
        std::vector<std::uint64_t> error_ntt(degree);
        for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
            mask[coefficient] = sampleUniformMod(engine, modulus);
            error_ntt[coefficient] = signedResidue(error[coefficient], modulus);
        }
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        plan.forward(mask);
        plan.forward(error_ntt);
        const auto secret_ntt = secret.ntt.slice(limb, 0);
        const auto reduced_message = message % result.plaintext_modulus;
        for (std::size_t slot = 0; slot < degree; ++slot) {
            body[slot] = modular_ntt::add(
                modular_ntt::add(modular_ntt::multiply(
                                     mask[slot], secret_ntt[slot], modulus),
                                 reduced_message, modulus),
                modular_ntt::multiply(result.plaintext_modulus % modulus,
                                      error_ntt[slot], modulus),
                modulus);
        }
    }
}

inline void CompactBGVEncryptScalar(CompactBGVScalarCiphertext &result,
                                    const std::uint64_t message,
                                    const CompactBGVSecretKey &secret)
{
    CompactBGVEncryptScalar(result, message, secret, TFHEpp::generator);
}

inline std::vector<std::uint64_t> CompactBGVDecryptPolynomial(
    const CompactBGVScalarCiphertext &ciphertext,
    const CompactBGVSecretKey &secret)
{
    const std::size_t limbs = ciphertext.limbs();
    if (limbs > secret.ntt.limbs())
        throw std::invalid_argument("compact BGV secret has too few limbs");
    std::vector<std::vector<std::uint64_t>> phase(
        limbs, std::vector<std::uint64_t>(degree));
    for (std::size_t limb = 0; limb < limbs; ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto mask = ciphertext.value.mask.slice(limb, 0);
        const auto body = ciphertext.value.body.slice(limb, 0);
        const auto secret_ntt = secret.ntt.slice(limb, 0);
        for (std::size_t slot = 0; slot < degree; ++slot)
            phase[limb][slot] = modular_ntt::subtract(
                body[slot],
                modular_ntt::multiply(mask[slot], secret_ntt[slot], modulus),
                modulus);
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        plan.inverse(phase[limb]);
    }
    FullModulusGadget crt(limbs, CompactBGV65536Parameters::gadget_digits);
    std::vector<std::uint64_t> result(degree), residues(limbs);
    for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
        for (std::size_t limb = 0; limb < limbs; ++limb)
            residues[limb] = phase[limb][coefficient];
        auto value = crt.reconstructCentered(residues);
        boost::multiprecision::cpp_int reduced =
            value % ciphertext.plaintext_modulus;
        if (reduced < 0) reduced += ciphertext.plaintext_modulus;
        result[coefficient] = reduced.convert_to<std::uint64_t>();
    }
    return result;
}

inline std::vector<boost::multiprecision::cpp_int>
CompactBGVCenteredPhaseCoefficients(
    const CompactBGVScalarCiphertext &ciphertext,
    const CompactBGVSecretKey &secret)
{
    const auto limbs = ciphertext.limbs();
    std::vector<std::vector<std::uint64_t>> phase(
        limbs, std::vector<std::uint64_t>(degree));
    for (std::size_t limb = 0; limb < limbs; ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto mask = ciphertext.value.mask.slice(limb, 0);
        const auto body = ciphertext.value.body.slice(limb, 0);
        const auto secret_ntt = secret.ntt.slice(limb, 0);
        for (std::size_t slot = 0; slot < degree; ++slot)
            phase[limb][slot] = modular_ntt::subtract(
                body[slot],
                modular_ntt::multiply(mask[slot], secret_ntt[slot], modulus),
                modulus);
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        plan.inverse(phase[limb]);
    }
    FullModulusGadget crt(limbs, CompactBGV65536Parameters::gadget_digits);
    std::vector<boost::multiprecision::cpp_int> result(degree);
    std::vector<std::uint64_t> residues(limbs);
    for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
        for (std::size_t limb = 0; limb < limbs; ++limb)
            residues[limb] = phase[limb][coefficient];
        result[coefficient] = crt.reconstructCentered(residues);
    }
    return result;
}

inline double CompactBGVMaxErrorLog2(
    const CompactBGVScalarCiphertext &ciphertext,
    const CompactBGVSecretKey &secret,
    const std::vector<std::uint64_t> &plaintext)
{
    if (plaintext.size() != degree)
        throw std::invalid_argument("compact BGV error plaintext size");
    const auto phase = CompactBGVCenteredPhaseCoefficients(ciphertext, secret);
    boost::multiprecision::cpp_int maximum = 0;
    for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
        boost::multiprecision::cpp_int difference =
            phase[coefficient] - plaintext[coefficient];
        if (difference % ciphertext.plaintext_modulus != 0)
            throw std::runtime_error(
                "compact BGV phase does not match expected plaintext");
        boost::multiprecision::cpp_int error =
            difference / ciphertext.plaintext_modulus;
        if (error < 0) error = -error;
        maximum = std::max(maximum, error);
    }
    if (maximum == 0) return -std::numeric_limits<double>::infinity();
    return static_cast<double>(boost::multiprecision::msb(maximum));
}

inline std::uint64_t CompactBGVDecryptScalar(
    const CompactBGVScalarCiphertext &ciphertext,
    const CompactBGVSecretKey &secret, const bool require_scalar = true)
{
    const auto polynomial = CompactBGVDecryptPolynomial(ciphertext, secret);
    if (require_scalar &&
        std::any_of(polynomial.begin() + 1, polynomial.end(),
                    [](const auto value) { return value != 0; }))
        throw std::runtime_error(
            "compact BGV ciphertext is not a scalar plaintext");
    return polynomial[0];
}

// OpenFHE/HElib-style BGV LSB modulus reduction on a nested prefix basis.
// The correction is a plaintext multiple, so every dropped prime preserves
// the message while scaling down the old invariant error.
inline void CompactBGVModSwitch(CompactBGVScalarCiphertext &result,
                                const CompactBGVScalarCiphertext &input)
{
    if (result.plaintext_modulus != input.plaintext_modulus ||
        result.limbs() == 0 || result.limbs() >= input.limbs())
        throw std::invalid_argument("invalid compact BGV modulus switch");
    std::vector<std::vector<std::uint64_t>> masks(
        input.limbs(), std::vector<std::uint64_t>(degree));
    std::vector<std::vector<std::uint64_t>> bodies(
        input.limbs(), std::vector<std::uint64_t>(degree));
    for (std::size_t limb = 0; limb < input.limbs(); ++limb) {
        std::copy(input.value.mask.slice(limb, 0).begin(),
                  input.value.mask.slice(limb, 0).end(), masks[limb].begin());
        std::copy(input.value.body.slice(limb, 0).begin(),
                  input.value.body.slice(limb, 0).end(), bodies[limb].begin());
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        plan.inverse(masks[limb]);
        plan.inverse(bodies[limb]);
    }

    for (std::size_t active = input.limbs(); active > result.limbs();
         --active) {
        const std::size_t dropped = active - 1;
        const auto ql = modular_ntt::degree65536_primes[dropped].value;
        const auto t_inverse =
            modular_ntt::invert(input.plaintext_modulus % ql, ql);
        const auto neg_t_inverse = modular_ntt::negate(t_inverse, ql);
        std::vector<std::uint64_t> ql_inverse(dropped);
        for (std::size_t limb = 0; limb < dropped; ++limb) {
            const auto modulus = modular_ntt::degree65536_primes[limb].value;
            ql_inverse[limb] = modular_ntt::invert(ql % modulus, modulus);
        }
        for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
            for (unsigned component = 0; component < 2; ++component) {
                const auto dropped_value = component == 0
                                               ? masks[dropped][coefficient]
                                               : bodies[dropped][coefficient];
                const auto delta_residue =
                    modular_ntt::multiply(dropped_value, neg_t_inverse, ql);
                const std::int64_t delta =
                    delta_residue <= ql / 2
                        ? static_cast<std::int64_t>(delta_residue)
                        : -static_cast<std::int64_t>(ql - delta_residue);
                for (std::size_t limb = 0; limb < dropped; ++limb) {
                    const auto modulus =
                        modular_ntt::degree65536_primes[limb].value;
                    const auto delta_mod = signedResidue(delta, modulus);
                    const auto correction = modular_ntt::multiply(
                        input.plaintext_modulus % modulus, delta_mod, modulus);
                    auto &value = component == 0 ? masks[limb][coefficient]
                                                 : bodies[limb][coefficient];
                    value = modular_ntt::multiply(
                        modular_ntt::add(value, correction, modulus),
                        ql_inverse[limb], modulus);
                }
            }
        }
    }

    for (std::size_t limb = 0; limb < result.limbs(); ++limb) {
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        plan.forward(masks[limb]);
        plan.forward(bodies[limb]);
        std::copy(masks[limb].begin(), masks[limb].end(),
                  result.value.mask.slice(limb, 0).begin());
        std::copy(bodies[limb].begin(), bodies[limb].end(),
                  result.value.body.slice(limb, 0).begin());
    }
}

inline void CompactBGVMultiply(CompactBGVScalarCiphertext &result,
                               const CompactBGVScalarCiphertext &left,
                               const CompactBGVScalarCiphertext &right,
                               const CompactBGVQuadraticHint &hint)
{
    if (left.plaintext_modulus != right.plaintext_modulus ||
        result.plaintext_modulus != left.plaintext_modulus)
        throw std::invalid_argument("compact BGV plaintext mismatch");
    if (left.limbs() != right.limbs() || result.limbs() != left.limbs() ||
        hint.u.limbs() < left.limbs())
        throw std::invalid_argument("compact BGV RNS level mismatch");
    for (std::size_t limb = 0; limb < left.limbs(); ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto a1 = left.value.mask.slice(limb, 0);
        const auto b1 = left.value.body.slice(limb, 0);
        const auto a2 = right.value.mask.slice(limb, 0);
        const auto b2 = right.value.body.slice(limb, 0);
        const auto u = hint.u.slice(limb, 0);
        const auto v = hint.v.slice(limb, 0);
        auto output_mask = result.value.mask.slice(limb, 0);
        auto output_body = result.value.body.slice(limb, 0);
        for (std::size_t slot = 0; slot < degree; ++slot) {
            const auto masks =
                modular_ntt::multiply(a1[slot], a2[slot], modulus);
            output_mask[slot] = modular_ntt::subtract(
                modular_ntt::add(
                    modular_ntt::multiply(a1[slot], b2[slot], modulus),
                    modular_ntt::multiply(a2[slot], b1[slot], modulus),
                    modulus),
                modular_ntt::multiply(masks, u[slot], modulus), modulus);
            output_body[slot] = modular_ntt::add(
                modular_ntt::multiply(b1[slot], b2[slot], modulus),
                modular_ntt::multiply(masks, v[slot], modulus), modulus);
        }
    }
}

inline void CompactBGVAdd(CompactBGVScalarCiphertext &result,
                          const CompactBGVScalarCiphertext &left,
                          const CompactBGVScalarCiphertext &right)
{
    if (left.plaintext_modulus != right.plaintext_modulus ||
        result.plaintext_modulus != left.plaintext_modulus ||
        left.limbs() != right.limbs() || result.limbs() != left.limbs())
        throw std::invalid_argument("compact BGV addition shape mismatch");
    add(result.value.mask, left.value.mask, right.value.mask);
    add(result.value.body, left.value.body, right.value.body);
}

inline void CompactBGVScalarMultiply(CompactBGVScalarCiphertext &result,
                                     const CompactBGVScalarCiphertext &input,
                                     const std::uint64_t scalar)
{
    if (result.plaintext_modulus != input.plaintext_modulus ||
        result.limbs() != input.limbs())
        throw std::invalid_argument(
            "compact BGV scalar multiplication shape mismatch");
    for (std::size_t limb = 0; limb < input.limbs(); ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        std::uint64_t plaintext_reduced = scalar % input.plaintext_modulus;
        const std::int64_t centered =
            plaintext_reduced <= input.plaintext_modulus / 2
                ? static_cast<std::int64_t>(plaintext_reduced)
                : -static_cast<std::int64_t>(input.plaintext_modulus -
                                             plaintext_reduced);
        const auto reduced = signedResidue(centered, modulus);
        for (std::size_t slot = 0; slot < degree; ++slot) {
            result.value.mask.slice(limb, 0)[slot] = modular_ntt::multiply(
                input.value.mask.slice(limb, 0)[slot], reduced, modulus);
            result.value.body.slice(limb, 0)[slot] = modular_ntt::multiply(
                input.value.body.slice(limb, 0)[slot], reduced, modulus);
        }
    }
}

inline void CompactBGVAddConstant(CompactBGVScalarCiphertext &value,
                                  const std::uint64_t constant)
{
    for (std::size_t limb = 0; limb < value.limbs(); ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto reduced = constant % modulus;
        auto body = value.value.body.slice(limb, 0);
        for (auto &slot : body) slot = modular_ntt::add(slot, reduced, modulus);
    }
}

inline CompactBGVScalarCiphertext CompactBGVZero(
    const std::uint64_t plaintext_modulus, const std::size_t limbs)
{
    return CompactBGVScalarCiphertext(plaintext_modulus, limbs);
}

inline CompactBGVScalarCiphertext CompactBGVAtLevel(
    const CompactBGVScalarCiphertext &input, const std::size_t limbs)
{
    if (limbs == 0 || limbs > input.limbs())
        throw std::invalid_argument("invalid compact BGV target level");
    if (limbs == input.limbs()) return input;
    CompactBGVScalarCiphertext result(input.plaintext_modulus, limbs);
    CompactBGVModSwitch(result, input);
    return result;
}

inline CompactBGVScalarCiphertext CompactBGVAddAligned(
    const CompactBGVScalarCiphertext &left,
    const CompactBGVScalarCiphertext &right)
{
    if (left.plaintext_modulus != right.plaintext_modulus)
        throw std::invalid_argument(
            "compact BGV aligned-add plaintext mismatch");
    const auto limbs = std::min(left.limbs(), right.limbs());
    const auto aligned_left = CompactBGVAtLevel(left, limbs);
    const auto aligned_right = CompactBGVAtLevel(right, limbs);
    CompactBGVScalarCiphertext result(left.plaintext_modulus, limbs);
    CompactBGVAdd(result, aligned_left, aligned_right);
    return result;
}

inline CompactBGVScalarCiphertext CompactBGVMultiplyAndDrop(
    const CompactBGVScalarCiphertext &left,
    const CompactBGVScalarCiphertext &right,
    const CompactBGVQuadraticHint &hint)
{
    if (left.plaintext_modulus != right.plaintext_modulus)
        throw std::invalid_argument(
            "compact BGV multiply-and-drop plaintext mismatch");
    const auto input_limbs = std::min(left.limbs(), right.limbs());
    if (input_limbs < 2)
        throw std::runtime_error("compact BGV multiplication exhausted levels");
    const auto aligned_left = CompactBGVAtLevel(left, input_limbs);
    const auto aligned_right = CompactBGVAtLevel(right, input_limbs);
    CompactBGVScalarCiphertext product(left.plaintext_modulus, input_limbs);
    CompactBGVMultiply(product, aligned_left, aligned_right, hint);
    CompactBGVScalarCiphertext result(left.plaintext_modulus, input_limbs - 1);
    CompactBGVModSwitch(result, product);
    return result;
}

struct CompactBGVBootstrapKeyDirectoryOptions {
    std::filesystem::path root;
    std::uintmax_t maximum_bytes = std::numeric_limits<std::uintmax_t>::max();
    bool resume = true;
};

struct CompactBGVBootstrapManifest {
    std::uint64_t magic{};
    std::uint32_t version{};
    std::uint32_t reserved{};
    std::uint64_t degree_value{};
    std::uint64_t limbs{};
    std::uint64_t rows{};
    std::uint64_t plaintext_prime{};
    std::uint64_t plaintext_square{};
    std::uint64_t digit_polynomial_fnv1a{};
    std::uint64_t key_file_checksum{};
    std::uint64_t quadratic_hint_checksum{};
    std::array<std::uint64_t, CompactBGV65536Parameters::trace_key_count>
        trace_key_checksums{};
    std::array<char, 65> gate_manifest_sha256{};
    std::array<char, 65> certificate_sha256{};
};

inline CompactBGVBootstrapManifest CompactBGVExpectedBootstrapManifest()
{
    CompactBGVBootstrapManifest manifest{
        CompactBGV65536Parameters::bootstrap_manifest_magic,
        CompactBGV65536Parameters::certificate_version,
        0,
        degree,
        CompactBGV65536Parameters::rns_limbs,
        CompactBGV65536Parameters::phase_lift_gadget_digits,
        CompactBGV65536Parameters::plaintext_prime,
        CompactBGV65536Parameters::plaintext_square,
        CompactBGV65536Parameters::digit_polynomial_fnv1a,
        0,
        0};
    std::copy_n(CompactBGV65536Parameters::scalar_gate_manifest_sha256, 64,
                manifest.gate_manifest_sha256.begin());
    std::copy_n(CompactBGV65536Parameters::certificate_sha256, 64,
                manifest.certificate_sha256.begin());
    return manifest;
}

inline std::filesystem::path CompactBGVBootstrapManifestPath(
    const std::filesystem::path &root)
{
    return root / "manifest.bin";
}

inline std::filesystem::path CompactBGVBootstrapKeyPath(
    const std::filesystem::path &root)
{
    return root / "scalar_bootstrap.key";
}

inline std::filesystem::path CompactBGVQuadraticHintPath(
    const std::filesystem::path &root)
{
    return root / "quadratic_hint.bin";
}

inline std::filesystem::path CompactBGVTraceKeyPath(
    const std::filesystem::path &root, const std::size_t index)
{
    if (index >= CompactBGV65536Parameters::trace_key_count)
        throw std::out_of_range("compact BGV trace-key index");
    return root / ("trace_" + std::to_string(index) + ".key");
}

inline std::size_t CompactBGVTraceKeyLimbs(const std::size_t index)
{
    if (index >= CompactBGV65536Parameters::trace_key_count)
        throw std::out_of_range("compact BGV trace-key level index");
    return index < 8 ? CompactBGV65536Parameters::rns_limbs
                     : CompactBGV65536Parameters::rns_limbs - 1;
}

inline std::array<std::uint32_t, CompactBGV65536Parameters::trace_key_count>
CompactBGVTraceExponents()
{
    std::array<std::uint32_t, CompactBGV65536Parameters::trace_key_count>
        result{};
    for (std::size_t index = 0; index + 1 < result.size(); ++index)
        result[index] = powModCyclotomic(5, std::size_t{1} << index);
    result.back() = 2 * degree - 1;
    return result;
}

inline std::uint64_t CompactBGVBootstrapKeyChecksum(
    const std::filesystem::path &path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
        throw std::runtime_error("cannot checksum compact BGV bootstrap key");
    std::uint64_t checksum = UINT64_C(1469598103934665603);
    std::array<char, 1 << 16> buffer{};
    while (stream) {
        stream.read(buffer.data(), buffer.size());
        const auto count = stream.gcount();
        for (std::streamsize i = 0; i < count; ++i) {
            checksum ^= static_cast<std::uint8_t>(buffer[i]);
            checksum *= UINT64_C(1099511628211);
        }
    }
    return checksum;
}

inline std::size_t CompactBGVQuadraticHintBytes()
{
    return checkedProduct({4, degree, CompactBGV65536Parameters::rns_limbs,
                           sizeof(std::uint64_t)});
}

inline void CompactBGVWriteQuadraticHint(const std::filesystem::path &path,
                                         const CompactBGVQuadraticHint &hint)
{
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream)
        throw std::runtime_error("cannot create compact BGV quadratic hint");
    for (const FrontierElement *element :
         {&hint.alpha, &hint.beta, &hint.u, &hint.v}) {
        if (element->width() != 1 ||
            element->limbs() != CompactBGV65536Parameters::rns_limbs)
            throw std::invalid_argument("compact BGV quadratic hint shape");
        stream.write(reinterpret_cast<const char *>(element->storage().data()),
                     static_cast<std::streamsize>(element->bytes()));
    }
    if (!stream)
        throw std::runtime_error("failed writing compact BGV quadratic hint");
}

inline void CompactBGVReadQuadraticHint(const std::filesystem::path &path,
                                        CompactBGVQuadraticHint &hint)
{
    if (std::filesystem::file_size(path) != CompactBGVQuadraticHintBytes())
        throw std::runtime_error("invalid compact BGV quadratic hint size");
    std::ifstream stream(path, std::ios::binary);
    for (FrontierElement *element : {&hint.alpha, &hint.beta, &hint.u, &hint.v})
        stream.read(reinterpret_cast<char *>(element->storage().data()),
                    static_cast<std::streamsize>(element->bytes()));
    if (!stream)
        throw std::runtime_error("failed reading compact BGV quadratic hint");
}

inline bool CompactBGVBootstrapKeyDirectoryManifestMatches(
    const std::filesystem::path &root)
{
    std::ifstream stream(CompactBGVBootstrapManifestPath(root),
                         std::ios::binary);
    CompactBGVBootstrapManifest actual{};
    stream.read(reinterpret_cast<char *>(&actual), sizeof(actual));
    const auto expected = CompactBGVExpectedBootstrapManifest();
    const auto key_path = CompactBGVBootstrapKeyPath(root);
    const auto hint_path = CompactBGVQuadraticHintPath(root);
    if (!(stream && std::filesystem::exists(key_path) &&
          std::filesystem::exists(hint_path) &&
          actual.magic == expected.magic &&
          actual.version == expected.version &&
          actual.degree_value == expected.degree_value &&
          actual.limbs == expected.limbs && actual.rows == expected.rows &&
          actual.plaintext_prime == expected.plaintext_prime &&
          actual.plaintext_square == expected.plaintext_square &&
          actual.digit_polynomial_fnv1a == expected.digit_polynomial_fnv1a &&
          actual.gate_manifest_sha256 == expected.gate_manifest_sha256 &&
          actual.certificate_sha256 == expected.certificate_sha256 &&
          actual.key_file_checksum == CompactBGVBootstrapKeyChecksum(key_path)))
        return false;
    if (actual.quadratic_hint_checksum !=
        CompactBGVBootstrapKeyChecksum(hint_path))
        return false;
    for (std::size_t index = 0;
         index < CompactBGV65536Parameters::trace_key_count; ++index) {
        const auto path = CompactBGVTraceKeyPath(root, index);
        if (!std::filesystem::exists(path) ||
            actual.trace_key_checksums[index] !=
                CompactBGVBootstrapKeyChecksum(path))
            return false;
    }
    return true;
}

inline bool CompactBGVBootstrapKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    if (!CompactBGVBootstrapKeyDirectoryManifestMatches(root) ||
        !std::filesystem::exists(CompactBGVBootstrapKeyPath(root)))
        return false;
    const auto estimate =
        estimateResources(1, CompactBGV65536Parameters::rns_limbs,
                          CompactBGV65536Parameters::phase_lift_gadget_digits);
    if (std::filesystem::file_size(CompactBGVBootstrapKeyPath(root)) !=
        sizeof(TransitionKeyHeader) + estimate.uniform_key_bytes)
        return false;
    if (std::filesystem::file_size(CompactBGVQuadraticHintPath(root)) !=
        CompactBGVQuadraticHintBytes())
        return false;
    for (std::size_t index = 0;
         index < CompactBGV65536Parameters::trace_key_count; ++index) {
        const auto trace_estimate =
            estimateResources(1, CompactBGVTraceKeyLimbs(index),
                              CompactBGV65536Parameters::gadget_digits);
        if (std::filesystem::file_size(CompactBGVTraceKeyPath(root, index)) !=
            sizeof(TransitionKeyHeader) + trace_estimate.uniform_key_bytes)
            return false;
    }
    return true;
}

inline bool CompactBGVUniformKeyFileMatches(const std::filesystem::path &path,
                                            const std::size_t limbs,
                                            const std::size_t rows)
{
    try {
        if (!std::filesystem::exists(path)) return false;
        UniformTransitionKeyProvider provider(path);
        return provider.header().width == 1 &&
               provider.header().limbs == limbs &&
               provider.header().rows == rows;
    }
    catch (const std::exception &) {
        return false;
    }
}

template <class Engine>
inline void CompactBGVWriteTransitionKey(const std::filesystem::path &path,
                                         const std::size_t rows,
                                         const FullModulusGadget &gadget,
                                         const FrontierElement &source_secret,
                                         const FrontierElement &target_secret,
                                         Engine &engine,
                                         const std::uintmax_t maximum_bytes)
{
    if (source_secret.width() != 1 || target_secret.width() != 1 ||
        source_secret.limbs() != target_secret.limbs() || gadget.rows() != rows)
        throw std::invalid_argument("invalid compact BGV transition secrets");
    const auto limbs = target_secret.limbs();
    UniformTransitionKeyWriter writer(path, 1, limbs, rows, maximum_bytes);
    std::vector<std::int64_t> error(degree);
    std::vector<std::uint64_t> mask(degree), body(degree), error_ntt(degree);
    boost::multiprecision::cpp_int gadget_power = 1;
    for (std::size_t row = 0; row < rows; ++row) {
        for (auto &coefficient : error) coefficient = sampleCBD20(engine);
        for (std::size_t limb = 0; limb < limbs; ++limb) {
            const auto modulus = modular_ntt::degree65536_primes[limb].value;
            for (std::size_t coefficient = 0; coefficient < degree;
                 ++coefficient) {
                mask[coefficient] = sampleUniformMod(engine, modulus);
                error_ntt[coefficient] =
                    signedResidue(error[coefficient], modulus);
            }
            modular_ntt::NegacyclicNTTPlan plan(
                degree, modular_ntt::degree65536_primes[limb]);
            plan.forward(mask);
            plan.forward(error_ntt);
            const auto gadget_residue =
                FullModulusGadget::residue(gadget_power, modulus);
            const auto source = source_secret.slice(limb, 0);
            const auto target = target_secret.slice(limb, 0);
            for (std::size_t slot = 0; slot < degree; ++slot) {
                body[slot] = modular_ntt::add(
                    modular_ntt::subtract(
                        modular_ntt::multiply(mask[slot], target[slot],
                                              modulus),
                        modular_ntt::multiply(gadget_residue, source[slot],
                                              modulus),
                        modulus),
                    modular_ntt::multiply(
                        CompactBGV65536Parameters::plaintext_square % modulus,
                        error_ntt[slot], modulus),
                    modulus);
            }
            writer.writeSlice(mask, body);
        }
        gadget_power *= gadget.base();
    }
    writer.finish();
}

template <class Engine>
inline void CompactBGVBootstrapKeyGenToDirectory(
    const CompactBGVBootstrapKeyDirectoryOptions &options,
    const CompactBGVKeyMaterial &key_material, Engine &engine)
{
    const auto &secret = key_material.secret;
    if (options.root.empty())
        throw std::invalid_argument("compact BGV key directory is empty");
    if (options.resume && CompactBGVBootstrapKeyDirectoryComplete(options.root))
        return;
    std::filesystem::create_directories(options.root);
    const auto required_bytes =
        CompactBGVEstimateResources().bootstrap_key_bytes;
    if (required_bytes > options.maximum_bytes)
        throw std::runtime_error(
            "compact BGV evaluation key exceeds configured disk budget");
    if (required_bytes > std::filesystem::space(options.root).available)
        throw std::runtime_error(
            "insufficient disk space for compact BGV evaluation key");
    const auto final_path = CompactBGVBootstrapKeyPath(options.root);
    const auto partial_path = final_path.string() + ".partial";
    const FullModulusGadget gadget(
        1, CompactBGV65536Parameters::phase_lift_gadget_digits);
    if (!CompactBGVUniformKeyFileMatches(
            final_path, CompactBGV65536Parameters::rns_limbs,
            CompactBGV65536Parameters::phase_lift_gadget_digits)) {
        CompactBGVWriteTransitionKey(
            partial_path, CompactBGV65536Parameters::phase_lift_gadget_digits,
            gadget, secret.ntt, secret.ntt, engine, options.maximum_bytes);
        std::filesystem::rename(partial_path, final_path);
    }

    const auto trace_exponents = CompactBGVTraceExponents();
    const std::array<RelabelEntry, 1> mapping{{{0, 1}}};
    for (std::size_t index = 0; index < trace_exponents.size(); ++index) {
        const auto limbs = CompactBGVTraceKeyLimbs(index);
        const FullModulusGadget trace_gadget(
            limbs, CompactBGV65536Parameters::gadget_digits);
        FrontierElement target_secret(1, limbs);
        for (std::size_t limb = 0; limb < limbs; ++limb)
            std::copy(secret.ntt.slice(limb, 0).begin(),
                      secret.ntt.slice(limb, 0).end(),
                      target_secret.slice(limb, 0).begin());
        FrontierElement automorphed_secret(1, limbs);
        auto actual_mapping = mapping;
        actual_mapping[0].automorphism = trace_exponents[index];
        relabel(automorphed_secret, target_secret, actual_mapping);
        const auto trace_final = CompactBGVTraceKeyPath(options.root, index);
        const auto trace_partial = trace_final.string() + ".partial";
        if (!CompactBGVUniformKeyFileMatches(
                trace_final, limbs, CompactBGV65536Parameters::gadget_digits)) {
            CompactBGVWriteTransitionKey(
                trace_partial, CompactBGV65536Parameters::gadget_digits,
                trace_gadget, automorphed_secret, target_secret, engine,
                options.maximum_bytes);
            std::filesystem::rename(trace_partial, trace_final);
        }
    }

    auto manifest = CompactBGVExpectedBootstrapManifest();
    manifest.key_file_checksum = CompactBGVBootstrapKeyChecksum(final_path);
    for (std::size_t index = 0; index < trace_exponents.size(); ++index)
        manifest.trace_key_checksums[index] = CompactBGVBootstrapKeyChecksum(
            CompactBGVTraceKeyPath(options.root, index));
    const auto hint_final = CompactBGVQuadraticHintPath(options.root);
    const auto hint_partial = hint_final.string() + ".partial";
    CompactBGVWriteQuadraticHint(hint_partial, key_material.quadratic_hint);
    std::filesystem::rename(hint_partial, hint_final);
    manifest.quadratic_hint_checksum =
        CompactBGVBootstrapKeyChecksum(hint_final);
    const auto manifest_partial =
        CompactBGVBootstrapManifestPath(options.root).string() + ".partial";
    {
        std::ofstream stream(manifest_partial,
                             std::ios::binary | std::ios::trunc);
        stream.write(reinterpret_cast<const char *>(&manifest),
                     sizeof(manifest));
        if (!stream)
            throw std::runtime_error(
                "failed writing compact BGV bootstrap manifest");
    }
    std::filesystem::rename(manifest_partial,
                            CompactBGVBootstrapManifestPath(options.root));
}

inline void CompactBGVBootstrapKeyGenToDirectory(
    const CompactBGVBootstrapKeyDirectoryOptions &options,
    const CompactBGVKeyMaterial &key_material)
{
    CompactBGVBootstrapKeyGenToDirectory(options, key_material,
                                         TFHEpp::generator);
}

class CompactBGVBootstrapFilesystemKeyProvider {
public:
    explicit CompactBGVBootstrapFilesystemKeyProvider(
        const std::filesystem::path &root)
        : provider_(CompactBGVBootstrapKeyPath(root))
    {
        if (!CompactBGVBootstrapKeyDirectoryComplete(root))
            throw std::runtime_error(
                "incomplete compact BGV bootstrap key directory");
        trace_.reserve(CompactBGV65536Parameters::trace_key_count);
        for (std::size_t index = 0;
             index < CompactBGV65536Parameters::trace_key_count; ++index)
            trace_.push_back(std::make_unique<UniformTransitionKeyProvider>(
                CompactBGVTraceKeyPath(root, index)));
        CompactBGVReadQuadraticHint(CompactBGVQuadraticHintPath(root), hint_);
    }

    const UniformTransitionKeyProvider &transition() const { return provider_; }
    const UniformTransitionKeyProvider &trace(const std::size_t index) const
    {
        if (index >= trace_.size())
            throw std::out_of_range("compact BGV trace provider index");
        return *trace_[index];
    }
    const CompactBGVQuadraticHint &quadraticHint() const { return hint_; }

private:
    UniformTransitionKeyProvider provider_;
    std::vector<std::unique_ptr<UniformTransitionKeyProvider>> trace_;
    CompactBGVQuadraticHint hint_;
};

inline void CompactBGVPhaseLiftTransition(
    CompactBGVScalarCiphertext &result, const CompactBGVScalarCiphertext &input,
    const CompactBGVBootstrapFilesystemKeyProvider &key,
    CompactBGVBootstrapTimings *timings = nullptr)
{
    using Clock = std::chrono::steady_clock;
    if (input.plaintext_modulus != CompactBGV65536Parameters::plaintext_prime ||
        input.limbs() != 1 ||
        result.plaintext_modulus !=
            CompactBGV65536Parameters::plaintext_square ||
        result.limbs() != CompactBGV65536Parameters::rns_limbs)
        throw std::invalid_argument(
            "compact BGV phase lift requires one input limb and full output");
    const auto lift_start = Clock::now();
    const auto low_modulus = modular_ntt::degree65536_primes[0].value;
    std::vector<std::uint64_t> low_mask(input.value.mask.slice(0, 0).begin(),
                                        input.value.mask.slice(0, 0).end());
    std::vector<std::uint64_t> low_body(input.value.body.slice(0, 0).begin(),
                                        input.value.body.slice(0, 0).end());
    modular_ntt::NegacyclicNTTPlan low_plan(degree,
                                            modular_ntt::degree65536_primes[0]);
    low_plan.inverse(low_mask);
    low_plan.inverse(low_body);
    const auto p = CompactBGV65536Parameters::plaintext_prime;
    std::vector<boost::multiprecision::cpp_int> remaining(degree);
    std::vector<std::int64_t> centered_body(degree);
    for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
        const auto scaled_mask =
            modular_ntt::multiply(low_mask[coefficient], p, low_modulus);
        const auto scaled_body =
            modular_ntt::multiply(low_body[coefficient], p, low_modulus);
        remaining[coefficient] =
            scaled_mask <= low_modulus / 2
                ? boost::multiprecision::cpp_int{scaled_mask}
                : -boost::multiprecision::cpp_int{low_modulus - scaled_mask};
        centered_body[coefficient] =
            scaled_body <= low_modulus / 2
                ? static_cast<std::int64_t>(scaled_body)
                : -static_cast<std::int64_t>(low_modulus - scaled_body);
    }
    const auto transition_start = Clock::now();
    const FullModulusGadget gadget(
        1, CompactBGV65536Parameters::phase_lift_gadget_digits);
    std::vector<std::vector<std::uint64_t>> digits(
        CompactBGV65536Parameters::rns_limbs,
        std::vector<std::uint64_t>(degree));
    std::vector<std::uint64_t> key_mask(degree), key_body(degree),
        product(degree);
    for (std::size_t limb = 0; limb < CompactBGV65536Parameters::rns_limbs;
         ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        auto output_mask = result.value.mask.slice(limb, 0);
        auto output_body = result.value.body.slice(limb, 0);
        std::fill(output_mask.begin(), output_mask.end(), 0);
        for (std::size_t coefficient = 0; coefficient < degree; ++coefficient)
            output_body[coefficient] =
                signedResidue(centered_body[coefficient], modulus);
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        plan.forward(output_body);
    }
    for (std::size_t row = 0;
         row < CompactBGV65536Parameters::phase_lift_gadget_digits; ++row) {
        for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
            const auto digit = gadget.takeBalancedDigit(
                remaining[coefficient],
                row + 1 == CompactBGV65536Parameters::phase_lift_gadget_digits);
            for (std::size_t limb = 0;
                 limb < CompactBGV65536Parameters::rns_limbs; ++limb)
                digits[limb][coefficient] = FullModulusGadget::residue(
                    digit, modular_ntt::degree65536_primes[limb].value);
        }
        for (std::size_t limb = 0; limb < CompactBGV65536Parameters::rns_limbs;
             ++limb) {
            const auto modulus = modular_ntt::degree65536_primes[limb].value;
            modular_ntt::NegacyclicNTTPlan plan(
                degree, modular_ntt::degree65536_primes[limb]);
            plan.forward(digits[limb]);
            key.transition().loadSlice(row, limb, 0, key_mask, key_body);
            addProductMod(result.value.mask.slice(limb, 0), digits[limb],
                          key_mask, modulus, product);
            addProductMod(result.value.body.slice(limb, 0), digits[limb],
                          key_body, modulus, product);
        }
    }
    const auto finish = Clock::now();
    if (timings != nullptr) {
        const auto milliseconds = [](const auto duration) {
            return std::chrono::duration<double, std::milli>(duration).count();
        };
        timings->lift_milliseconds =
            milliseconds(transition_start - lift_start);
        timings->transition_milliseconds =
            milliseconds(finish - transition_start);
        timings->divide_milliseconds = 0;
    }
}

inline void CompactBGVTraceProjectConstant(
    CompactBGVScalarCiphertext &result, const CompactBGVScalarCiphertext &input,
    const CompactBGVBootstrapFilesystemKeyProvider &keys)
{
    if (input.plaintext_modulus !=
            CompactBGV65536Parameters::plaintext_square ||
        input.limbs() != CompactBGV65536Parameters::rns_limbs)
        throw std::invalid_argument("invalid compact BGV trace input");
    CompactBGVScalarCiphertext current = input;
    const auto exponents = CompactBGVTraceExponents();
    const std::array<RelabelEntry, 1> identity_mapping{{{0, 1}}};
    for (std::size_t index = 0; index < exponents.size(); ++index) {
        auto mapping = identity_mapping;
        mapping[0].automorphism = exponents[index];
        CompactBGVScalarCiphertext automorphed(input.plaintext_modulus,
                                               current.limbs());
        relabel(automorphed.value.mask, current.value.mask, mapping);
        relabel(automorphed.value.body, current.value.body, mapping);
        CompactBGVScalarCiphertext switched(input.plaintext_modulus,
                                            current.limbs());
        applyFullModulusTransition(switched.value, automorphed.value,
                                   keys.trace(index));
        auto next = CompactBGVAddAligned(current, switched);
        current = std::move(next);
        if (index == 7 || index + 1 == exponents.size()) {
            CompactBGVScalarCiphertext reduced(input.plaintext_modulus,
                                               current.limbs() - 1);
            CompactBGVModSwitch(reduced, current);
            current = std::move(reduced);
        }
    }
    // 65536^(-1) mod 65537^2.
    constexpr std::uint64_t trace_normalizer = UINT64_C(4295032831);
    CompactBGVScalarCiphertext normalized(input.plaintext_modulus,
                                          current.limbs());
    CompactBGVScalarMultiply(normalized, current, trace_normalizer);
    result = std::move(normalized);
}

inline void CompactBGVPolynomialEval(
    CompactBGVScalarCiphertext &result,
    const std::vector<std::uint64_t> &coefficients,
    const CompactBGVScalarCiphertext &input,
    const CompactBGVQuadraticHint &hint)
{
    if (input.plaintext_modulus != CompactBGV65536Parameters::plaintext_square)
        throw std::invalid_argument("invalid compact BGV polynomial input");
    if (coefficients.empty()) {
        result = CompactBGVZero(input.plaintext_modulus, input.limbs());
        return;
    }
    const std::size_t polynomial_degree = coefficients.size() - 1;
    if (polynomial_degree == 0) {
        result = CompactBGVZero(input.plaintext_modulus, input.limbs());
        CompactBGVAddConstant(result, coefficients[0]);
        return;
    }
    bool odd_polynomial = coefficients[0] == 0;
    for (std::size_t power = 2; power < coefficients.size(); power += 2)
        odd_polynomial = odd_polynomial && coefficients[power] == 0;
    if (odd_polynomial) {
        auto square = CompactBGVMultiplyAndDrop(input, input, hint);
        std::vector<std::uint64_t> reduced;
        reduced.reserve((coefficients.size() + 1) / 2);
        for (std::size_t power = 1; power < coefficients.size(); power += 2)
            reduced.push_back(coefficients[power]);
        CompactBGVScalarCiphertext inner(input.plaintext_modulus,
                                         square.limbs());
        CompactBGVPolynomialEval(inner, reduced, square, hint);
        result = CompactBGVMultiplyAndDrop(inner, input, hint);
        return;
    }
    constexpr std::size_t baby_count = 3;
    constexpr std::size_t giant_count = 6;
    if (polynomial_degree > baby_count * (std::size_t{1} << giant_count))
        throw std::invalid_argument(
            "compact BGV digit polynomial exceeds fixed BSGS reach");

    const auto scale = [](const CompactBGVScalarCiphertext &value,
                          const std::uint64_t scalar) {
        CompactBGVScalarCiphertext scaled(value.plaintext_modulus,
                                          value.limbs());
        CompactBGVScalarMultiply(scaled, value, scalar);
        return scaled;
    };

    std::vector<CompactBGVScalarCiphertext> baby;
    baby.reserve(baby_count + 1);
    baby.emplace_back(input.plaintext_modulus, input.limbs());
    baby.push_back(input);
    baby.push_back(CompactBGVMultiplyAndDrop(baby[1], baby[1], hint));
    baby.push_back(CompactBGVMultiplyAndDrop(baby[2], baby[1], hint));

    std::vector<CompactBGVScalarCiphertext> giant;
    giant.reserve(giant_count);
    giant.push_back(baby[baby_count]);
    for (std::size_t index = 1; index < giant_count; ++index)
        giant.push_back(CompactBGVMultiplyAndDrop(giant[index - 1],
                                                  giant[index - 1], hint));

    std::function<CompactBGVScalarCiphertext(std::size_t, std::size_t,
                                             std::size_t)>
        evaluate = [&](const std::size_t offset, const std::size_t length,
                       const std::size_t level) {
            CompactBGVScalarCiphertext accumulator(input.plaintext_modulus,
                                                   input.limbs());
            if (length == 0) return accumulator;
            CompactBGVAddConstant(accumulator, coefficients[offset]);
            if (length == 1) return accumulator;
            if (level == 0) {
                const auto last = std::min(length - 1, baby_count);
                for (std::size_t power = 1; power <= last; ++power) {
                    if (coefficients[offset + power] == 0) continue;
                    const auto term =
                        scale(baby[power], coefficients[offset + power]);
                    accumulator = CompactBGVAddAligned(accumulator, term);
                }
                return accumulator;
            }
            const std::size_t split =
                std::min(length, baby_count * (std::size_t{1} << (level - 1)));
            auto lower = evaluate(offset, split, level - 1);
            if (split == length) return lower;
            auto upper = evaluate(offset + split, length - split, level - 1);
            const auto product =
                CompactBGVMultiplyAndDrop(upper, giant[level - 1], hint);
            return CompactBGVAddAligned(lower, product);
        };
    result = evaluate(0, coefficients.size(), giant_count);
}

inline void CompactBGVExactDivideByPlaintextPrime(
    CompactBGVScalarCiphertext &result, const CompactBGVScalarCiphertext &input)
{
    if (input.plaintext_modulus !=
            CompactBGV65536Parameters::plaintext_square ||
        result.plaintext_modulus !=
            CompactBGV65536Parameters::plaintext_prime ||
        result.limbs() != input.limbs())
        throw std::invalid_argument("invalid compact BGV exact division");
    for (std::size_t limb = 0; limb < input.limbs(); ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto inverse = modular_ntt::invert(
            CompactBGV65536Parameters::plaintext_prime % modulus, modulus);
        for (std::size_t slot = 0; slot < degree; ++slot) {
            result.value.mask.slice(limb, 0)[slot] = modular_ntt::multiply(
                input.value.mask.slice(limb, 0)[slot], inverse, modulus);
            result.value.body.slice(limb, 0)[slot] = modular_ntt::multiply(
                input.value.body.slice(limb, 0)[slot], inverse, modulus);
        }
    }
}

inline void CompactBGVBootstrap(
    CompactBGVScalarCiphertext &result, const CompactBGVScalarCiphertext &input,
    const CompactBGVBootstrapFilesystemKeyProvider &keys,
    CompactBGVBootstrapTimings *timings = nullptr)
{
    using Clock = std::chrono::steady_clock;
    if (input.limbs() != 1 ||
        input.plaintext_modulus != CompactBGV65536Parameters::plaintext_prime)
        throw std::invalid_argument(
            "compact BGV bootstrap requires low input and full output");
    CompactBGVScalarCiphertext lifted(
        CompactBGV65536Parameters::plaintext_square,
        CompactBGV65536Parameters::rns_limbs);
    CompactBGVBootstrapTimings phase_timings;
    CompactBGVPhaseLiftTransition(lifted, input, keys, &phase_timings);
    const auto trace_start = Clock::now();
    CompactBGVScalarCiphertext projected(
        CompactBGV65536Parameters::plaintext_square,
        CompactBGV65536Parameters::rns_limbs);
    CompactBGVTraceProjectConstant(projected, lifted, keys);
    const auto digit_start = Clock::now();
    const auto polynomial = digitext::GetLowestDigitRemovalPolynomialOverRange(
        CompactBGV65536Parameters::plaintext_prime,
        CompactBGV65536Parameters::digit_error_bound);
    CompactBGVScalarCiphertext removed(
        CompactBGV65536Parameters::plaintext_square,
        CompactBGV65536Parameters::rns_limbs);
    CompactBGVPolynomialEval(removed, polynomial, projected,
                             keys.quadraticHint());
    const auto divide_start = Clock::now();
    CompactBGVScalarCiphertext divided(
        CompactBGV65536Parameters::plaintext_prime, removed.limbs());
    CompactBGVExactDivideByPlaintextPrime(divided, removed);
    result = std::move(divided);
    const auto finish = Clock::now();
    if (timings != nullptr) {
        const auto milliseconds = [](const auto duration) {
            return std::chrono::duration<double, std::milli>(duration).count();
        };
        timings->lift_milliseconds = phase_timings.lift_milliseconds;
        timings->transition_milliseconds =
            phase_timings.transition_milliseconds;
        timings->trace_milliseconds = milliseconds(digit_start - trace_start);
        timings->digit_removal_milliseconds =
            milliseconds(divide_start - digit_start);
        timings->divide_milliseconds = milliseconds(finish - divide_start);
    }
}

}  // namespace TFHEpp::compact_cover_bgv
