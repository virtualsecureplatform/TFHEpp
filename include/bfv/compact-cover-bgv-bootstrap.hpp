#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "bfv/compact-cover-bgv.hpp"

namespace TFHEpp::compact_cover_bgv {

struct CompactBGV65536Parameters {
    static constexpr std::size_t ring_degree = degree;
    static constexpr std::size_t rns_limbs = 15;
    static constexpr std::size_t gadget_digits = 5;
    static constexpr std::size_t secret_weight = 32;
    static constexpr unsigned cbd_eta = 20;
    static constexpr std::uint64_t plaintext_prime = 65537;
    static constexpr std::uint64_t plaintext_square =
        plaintext_prime * plaintext_prime;
    static constexpr std::uint32_t certificate_version = 1;
    static constexpr std::uint64_t bootstrap_manifest_magic =
        UINT64_C(0x4342475642543031);
    static constexpr const char *scalar_gate_manifest_sha256 =
        "209b8826908383c92bce2ea41f27eda9febbf250e69e2c5541f5c18b76e454f0";
    static constexpr const char *certificate_sha256 =
        "143be2a21f31cd5bb2bb3d6359b81f1103a05329639fd42152c3503669dd39b1";
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
            CompactBGV65536Parameters::plaintext_prime)
        : value(1, CompactBGV65536Parameters::rns_limbs),
          plaintext_modulus(plaintext)
    {
    }
};

struct CompactBGVBootstrapTimings {
    double lift_milliseconds{};
    double transition_milliseconds{};
    double divide_milliseconds{};
};

struct CompactBGVResourceEstimate {
    std::size_t ciphertext_bytes{};
    std::size_t bootstrap_key_bytes{};
    std::size_t transition_scratch_bytes{};
};

inline CompactBGVResourceEstimate CompactBGVEstimateResources()
{
    const auto estimate =
        estimateResources(1, CompactBGV65536Parameters::rns_limbs,
                          CompactBGV65536Parameters::gadget_digits);
    return {estimate.ciphertext_bytes,
            sizeof(TransitionKeyHeader) + estimate.uniform_key_bytes,
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
    for (std::size_t limb = 0; limb < CompactBGV65536Parameters::rns_limbs;
         ++limb) {
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

inline std::vector<std::uint64_t> CompactBGVDecryptPolynomial(
    const CompactBGVScalarCiphertext &ciphertext,
    const CompactBGVSecretKey &secret)
{
    constexpr std::size_t limbs = CompactBGV65536Parameters::rns_limbs;
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

inline void CompactBGVMultiply(CompactBGVScalarCiphertext &result,
                               const CompactBGVScalarCiphertext &left,
                               const CompactBGVScalarCiphertext &right,
                               const CompactBGVQuadraticHint &hint)
{
    if (left.plaintext_modulus != right.plaintext_modulus ||
        result.plaintext_modulus != left.plaintext_modulus)
        throw std::invalid_argument("compact BGV plaintext mismatch");
    for (std::size_t limb = 0; limb < CompactBGV65536Parameters::rns_limbs;
         ++limb) {
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
    std::uint64_t key_file_checksum{};
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
        CompactBGV65536Parameters::gadget_digits,
        CompactBGV65536Parameters::plaintext_prime,
        CompactBGV65536Parameters::plaintext_square,
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

inline bool CompactBGVBootstrapKeyDirectoryManifestMatches(
    const std::filesystem::path &root)
{
    std::ifstream stream(CompactBGVBootstrapManifestPath(root),
                         std::ios::binary);
    CompactBGVBootstrapManifest actual{};
    stream.read(reinterpret_cast<char *>(&actual), sizeof(actual));
    const auto expected = CompactBGVExpectedBootstrapManifest();
    const auto key_path = CompactBGVBootstrapKeyPath(root);
    return stream && std::filesystem::exists(key_path) &&
           actual.magic == expected.magic &&
           actual.version == expected.version &&
           actual.degree_value == expected.degree_value &&
           actual.limbs == expected.limbs && actual.rows == expected.rows &&
           actual.plaintext_prime == expected.plaintext_prime &&
           actual.plaintext_square == expected.plaintext_square &&
           actual.gate_manifest_sha256 == expected.gate_manifest_sha256 &&
           actual.certificate_sha256 == expected.certificate_sha256 &&
           actual.key_file_checksum == CompactBGVBootstrapKeyChecksum(key_path);
}

inline bool CompactBGVBootstrapKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    if (!CompactBGVBootstrapKeyDirectoryManifestMatches(root) ||
        !std::filesystem::exists(CompactBGVBootstrapKeyPath(root)))
        return false;
    const auto estimate =
        estimateResources(1, CompactBGV65536Parameters::rns_limbs,
                          CompactBGV65536Parameters::gadget_digits);
    return std::filesystem::file_size(CompactBGVBootstrapKeyPath(root)) ==
           sizeof(TransitionKeyHeader) + estimate.uniform_key_bytes;
}

template <class Engine>
inline void CompactBGVBootstrapKeyGenToDirectory(
    const CompactBGVBootstrapKeyDirectoryOptions &options,
    const CompactBGVSecretKey &secret, Engine &engine)
{
    if (options.root.empty())
        throw std::invalid_argument("compact BGV key directory is empty");
    if (options.resume && CompactBGVBootstrapKeyDirectoryComplete(options.root))
        return;
    std::filesystem::create_directories(options.root);
    const auto final_path = CompactBGVBootstrapKeyPath(options.root);
    const auto partial_path = final_path.string() + ".partial";
    const FullModulusGadget gadget(CompactBGV65536Parameters::rns_limbs,
                                   CompactBGV65536Parameters::gadget_digits);
    {
        UniformTransitionKeyWriter writer(
            partial_path, 1, CompactBGV65536Parameters::rns_limbs,
            CompactBGV65536Parameters::gadget_digits, options.maximum_bytes);
        std::vector<std::int64_t> error(degree);
        std::vector<std::uint64_t> mask(degree), body(degree),
            error_ntt(degree);
        boost::multiprecision::cpp_int gadget_power = 1;
        for (std::size_t row = 0;
             row < CompactBGV65536Parameters::gadget_digits; ++row) {
            for (auto &coefficient : error) coefficient = sampleCBD20(engine);
            for (std::size_t limb = 0;
                 limb < CompactBGV65536Parameters::rns_limbs; ++limb) {
                const auto modulus =
                    modular_ntt::degree65536_primes[limb].value;
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
                const auto secret_ntt = secret.ntt.slice(limb, 0);
                for (std::size_t slot = 0; slot < degree; ++slot) {
                    body[slot] = modular_ntt::add(
                        modular_ntt::subtract(
                            modular_ntt::multiply(mask[slot], secret_ntt[slot],
                                                  modulus),
                            modular_ntt::multiply(gadget_residue,
                                                  secret_ntt[slot], modulus),
                            modulus),
                        modular_ntt::multiply(
                            CompactBGV65536Parameters::plaintext_square %
                                modulus,
                            error_ntt[slot], modulus),
                        modulus);
                }
                writer.writeSlice(mask, body);
            }
            gadget_power *= gadget.base();
        }
        writer.finish();
    }
    std::filesystem::rename(partial_path, final_path);
    auto manifest = CompactBGVExpectedBootstrapManifest();
    manifest.key_file_checksum = CompactBGVBootstrapKeyChecksum(final_path);
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

class CompactBGVBootstrapFilesystemKeyProvider {
public:
    explicit CompactBGVBootstrapFilesystemKeyProvider(
        const std::filesystem::path &root)
        : provider_(CompactBGVBootstrapKeyPath(root))
    {
        if (!CompactBGVBootstrapKeyDirectoryComplete(root))
            throw std::runtime_error(
                "incomplete compact BGV bootstrap key directory");
    }

    const UniformTransitionKeyProvider &transition() const { return provider_; }

private:
    UniformTransitionKeyProvider provider_;
};

inline void CompactBGVBootstrap(
    CompactBGVScalarCiphertext &result, const CompactBGVScalarCiphertext &input,
    const CompactBGVBootstrapFilesystemKeyProvider &key,
    CompactBGVBootstrapTimings *timings = nullptr)
{
    using Clock = std::chrono::steady_clock;
    if (input.plaintext_modulus != CompactBGV65536Parameters::plaintext_prime)
        throw std::invalid_argument(
            "compact BGV bootstrap input must use plaintext p");
    CompactBGVScalarCiphertext lifted(
        CompactBGV65536Parameters::plaintext_square);
    const auto lift_start = Clock::now();
    for (std::size_t limb = 0; limb < CompactBGV65536Parameters::rns_limbs;
         ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto scalar =
            CompactBGV65536Parameters::plaintext_prime % modulus;
        for (std::size_t slot = 0; slot < degree; ++slot) {
            lifted.value.mask.slice(limb, 0)[slot] = modular_ntt::multiply(
                input.value.mask.slice(limb, 0)[slot], scalar, modulus);
            lifted.value.body.slice(limb, 0)[slot] = modular_ntt::multiply(
                input.value.body.slice(limb, 0)[slot], scalar, modulus);
        }
    }
    const auto transition_start = Clock::now();
    CompactBGVScalarCiphertext switched(
        CompactBGV65536Parameters::plaintext_square);
    applyFullModulusTransition(switched.value, lifted.value, key.transition());

    const auto divide_start = Clock::now();
    result.plaintext_modulus = CompactBGV65536Parameters::plaintext_prime;
    for (std::size_t limb = 0; limb < CompactBGV65536Parameters::rns_limbs;
         ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto inverse = modular_ntt::invert(
            CompactBGV65536Parameters::plaintext_prime % modulus, modulus);
        for (std::size_t slot = 0; slot < degree; ++slot) {
            result.value.mask.slice(limb, 0)[slot] = modular_ntt::multiply(
                switched.value.mask.slice(limb, 0)[slot], inverse, modulus);
            result.value.body.slice(limb, 0)[slot] = modular_ntt::multiply(
                switched.value.body.slice(limb, 0)[slot], inverse, modulus);
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
            milliseconds(divide_start - transition_start);
        timings->divide_milliseconds = milliseconds(finish - divide_start);
    }
}

}  // namespace TFHEpp::compact_cover_bgv
