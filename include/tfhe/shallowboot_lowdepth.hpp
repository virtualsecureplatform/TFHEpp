#pragma once

// Experimental implementation of Algorithm 3 from Jain--Lin--Liu--Saha,
// ePrint 2026/1730.  This header deliberately uses an explicit prime-modulus
// NTT ring rather than TFHEpp's power-of-two torus types: Algorithm 3 needs a
// Binary-NTT secret and ciphertexts that remain in NTT form across a product
// tree. It implements PBC factor compilation, multi-stage nested-RNS BGV
// boundaries, modular-inverse LSB-to-MSB conversion, LUT evaluation,
// extraction, and raw LWE key switching. The secure executable parameter is
// exposed through this experimental API; TFHEpp TLWE/EvalKey serialization is
// intentionally separate.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#ifdef USE_BLAKE3
#include "../BLAKE3PRNG.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../modular_ntt.hpp"

namespace TFHEpp::shallowboot::lowdepth {

class Ring {
public:
    Ring(const std::size_t degree, const modular_ntt::PrimeModulus modulus)
        : degree_(degree), modulus_(modulus.value), plan_(degree, modulus)
    {
    }

    std::size_t degree() const { return degree_; }
    std::uint64_t modulus() const { return modulus_; }

    void forward(std::vector<std::uint64_t> &values) const
    {
        requireSize(values);
        plan_.forward(values);
    }

    void inverse(std::vector<std::uint64_t> &values) const
    {
        requireSize(values);
        plan_.inverse(values);
    }

    std::vector<std::uint64_t> monomialNTT(const std::size_t exponent) const
    {
        std::vector<std::uint64_t> result(degree_);
        const std::size_t reduced = exponent % (2 * degree_);
        if (reduced < degree_)
            result[reduced] = 1;
        else
            result[reduced - degree_] = modulus_ - 1;
        forward(result);
        return result;
    }

    const std::vector<std::uint64_t> &cachedMonomialNTT(
        const std::size_t exponent) const
    {
        const std::size_t reduced = exponent % (2 * degree_);
        const auto found = monomial_cache_.find(reduced);
        if (found != monomial_cache_.end()) return found->second;
        auto [inserted, unused] = monomial_cache_.emplace(
            reduced, monomialNTT(reduced));
        return inserted->second;
    }

    void requireSize(const std::vector<std::uint64_t> &values) const
    {
        if (values.size() != degree_)
            throw std::invalid_argument("Binary-NTT vector has wrong degree");
    }

private:
    std::size_t degree_;
    std::uint64_t modulus_;
    modular_ntt::NegacyclicNTTPlan plan_;
    mutable std::unordered_map<std::size_t, std::vector<std::uint64_t>>
        monomial_cache_;
};

// The secret is binary in NTT slots, equivalently s^2=s.  The general
// quadratic-hint form s^2=u*s+v is retained in the representation so the
// multiplication formula matches Algorithm 3 lines 14--17 exactly.
struct Secret {
    std::vector<std::uint64_t> slots;
    std::vector<std::uint64_t> quadratic_u;
    std::vector<std::uint64_t> quadratic_v;
};

struct Ciphertext {
    std::vector<std::uint64_t> a;
    std::vector<std::uint64_t> b;
};

struct ProductTreeStats {
    std::uint32_t pointwise_multiplication_layers = 0;
    std::uint32_t pointwise_ciphertext_multiplications = 0;
    std::uint32_t ntt_layers_inside_tree = 0;
};

struct ModulusScheduleStats : ProductTreeStats {
    std::uint32_t modulus_boundaries = 0;
};

struct BootstrappingKey {
    std::vector<Ciphertext> encrypted_lwe_bits;
};

struct KeySwitchKey {
    std::uint32_t basebit = 0;
    std::vector<Ciphertext> digits;
};

// A conventional RLWE key switch, in contrast to the slotwise NTT helper
// below.  Algorithm 3's modulus boundary explicitly performs INTT before its
// first key switch: decomposing the coefficient representation keeps the
// error sampled in `EncryptNTT` small after the result returns to NTT form.
// The key is larger (one RLWE ciphertext per coefficient and digit), but it
// is the error-correct representation required at this boundary.
struct CoefficientKeySwitchKey {
    std::uint32_t basebit = 0;
    std::vector<std::vector<Ciphertext>> entries;
};

struct LWEPhase {
    std::vector<std::uint32_t> a;
    std::uint32_t b;
    std::uint32_t modulus;
};

struct LWECiphertext {
    std::vector<std::uint64_t> a;
    std::uint64_t b = 0;
    std::uint64_t modulus = 0;
};

struct LWEKeySwitchKey {
    std::uint32_t basebit = 0;
    bool balanced = false;
    std::vector<std::vector<LWECiphertext>> entries;
};

inline void validate(const Ring &ring, const Secret &secret);

// A compact, deterministic PBC encoder for Algorithm 2/3. Every source LWE
// coordinate is copied into `copies` public buckets. Key generation matches
// the sparse support to distinct buckets; retry with a new salt if matching
// fails. This is the same operational contract as the paper's PBC.Ecd/Sched
// interface, without claiming its asymptotic failure bound.
struct PBCSchedule {
    std::uint32_t source_dimension = 0;
    std::uint32_t bucket_count = 0;
    std::uint32_t copies = 0;
    std::uint64_t salt = 0;
    std::vector<std::vector<std::uint32_t>> bucket_indices;
    std::vector<std::int32_t> selected_source;
};

inline std::uint64_t PBCHash(const std::uint32_t index,
                             const std::uint32_t copy,
                             const std::uint64_t salt)
{
    std::uint64_t value = static_cast<std::uint64_t>(index) +
                          0x9e3779b97f4a7c15ULL * (copy + 1) + salt;
    value ^= value >> 30;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27;
    value *= 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

inline PBCSchedule BuildPBCSchedule(
    const std::uint32_t source_dimension, const std::uint32_t bucket_count,
    const std::uint32_t copies, const std::span<const std::uint32_t> support,
    const std::uint64_t salt)
{
    if (source_dimension == 0 || bucket_count < support.size() || copies == 0)
        throw std::invalid_argument("invalid PBC schedule parameters");
    PBCSchedule result;
    result.source_dimension = source_dimension;
    result.bucket_count = bucket_count;
    result.copies = copies;
    result.salt = salt;
    result.bucket_indices.resize(bucket_count);
    std::vector<std::vector<std::uint32_t>> candidates(source_dimension);
    for (std::uint32_t index = 0; index < source_dimension; index++) {
        for (std::uint32_t copy = 0; copy < copies; copy++) {
            const std::uint32_t bucket = static_cast<std::uint32_t>(
                PBCHash(index, copy, salt) % bucket_count);
            if (std::find(candidates[index].begin(), candidates[index].end(),
                          bucket) == candidates[index].end()) {
                candidates[index].push_back(bucket);
                result.bucket_indices[bucket].push_back(index);
            }
        }
    }
    for (const std::uint32_t index : support)
        if (index >= source_dimension)
            throw std::invalid_argument("PBC support index is out of range");

    std::vector<std::int32_t> bucket_owner(bucket_count, -1);
    const auto augment = [&](const auto &self, const std::uint32_t support_pos,
                             std::vector<bool> &seen) -> bool {
        const std::uint32_t index = support[support_pos];
        for (const std::uint32_t bucket : candidates[index]) {
            if (seen[bucket]) continue;
            seen[bucket] = true;
            if (bucket_owner[bucket] < 0 ||
                self(self, static_cast<std::uint32_t>(bucket_owner[bucket]),
                     seen)) {
                bucket_owner[bucket] = static_cast<std::int32_t>(support_pos);
                return true;
            }
        }
        return false;
    };
    for (std::uint32_t support_pos = 0; support_pos < support.size();
         support_pos++) {
        std::vector<bool> seen(bucket_count);
        if (!augment(augment, support_pos, seen))
            throw std::runtime_error(
                "PBC matching failed; retry with a different salt");
    }
    result.selected_source.assign(bucket_count, -1);
    for (std::uint32_t bucket = 0; bucket < bucket_count; bucket++)
        if (bucket_owner[bucket] >= 0)
            result.selected_source[bucket] =
                static_cast<std::int32_t>(support[bucket_owner[bucket]]);
    return result;
}

inline PBCSchedule BuildPBCScheduleWithRetry(
    const std::uint32_t source_dimension, const std::uint32_t bucket_count,
    const std::uint32_t copies, const std::span<const std::uint32_t> support,
    const std::uint64_t initial_salt, const std::uint32_t attempts = 128)
{
    for (std::uint32_t attempt = 0; attempt < attempts; attempt++) {
        try {
            return BuildPBCSchedule(source_dimension, bucket_count, copies,
                                    support, initial_salt + attempt);
        }
        catch (const std::runtime_error &) {
            // Matching failures are expected PBC retries; malformed input is
            // rejected by BuildPBCSchedule before this point.
        }
    }
    throw std::runtime_error("PBC matching failed after all schedule retries");
}

struct PBCBootstrappingKey {
    std::vector<std::vector<Ciphertext>> entries;
    std::vector<Ciphertext> dummies;
};

template <class URBG>
PBCBootstrappingKey PBCBootstrappingKeyGen(
    const Ring &ring, const Secret &secret, const PBCSchedule &schedule,
    const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, secret);
    PBCBootstrappingKey result;
    result.entries.resize(schedule.bucket_count);
    result.dummies.reserve(schedule.bucket_count);
    for (std::size_t bucket = 0; bucket < schedule.bucket_indices.size();
         bucket++) {
        for (const std::uint32_t source : schedule.bucket_indices[bucket]) {
            std::vector<std::uint64_t> message(ring.degree());
            message[0] = schedule.selected_source[bucket] ==
                                 static_cast<std::int32_t>(source)
                             ? 1
                             : 0;
            ring.forward(message);
            result.entries[bucket].push_back(
                EncryptNTT(ring, secret, message, coefficient_noise_stddev,
                           engine));
        }
        std::vector<std::uint64_t> dummy(ring.degree());
        dummy[0] = schedule.selected_source[bucket] < 0 ? 1 : 0;
        ring.forward(dummy);
        result.dummies.push_back(
            EncryptNTT(ring, secret, dummy, coefficient_noise_stddev, engine));
    }
    return result;
}

inline Ciphertext AddCiphertexts(const Ring &ring, const Ciphertext &left,
                                 const Ciphertext &right)
{
    ring.requireSize(left.a);
    ring.requireSize(left.b);
    ring.requireSize(right.a);
    ring.requireSize(right.b);
    Ciphertext result;
    result.a.resize(ring.degree());
    result.b.resize(ring.degree());
#ifdef USE_HEXL
    intel::hexl::EltwiseAddMod(result.a.data(), left.a.data(), right.a.data(),
                               ring.degree(), ring.modulus());
    intel::hexl::EltwiseAddMod(result.b.data(), left.b.data(), right.b.data(),
                               ring.degree(), ring.modulus());
#else
    for (std::size_t i = 0; i < ring.degree(); i++) {
        result.a[i] = modular_ntt::add(left.a[i], right.a[i], ring.modulus());
        result.b[i] = modular_ntt::add(left.b[i], right.b[i], ring.modulus());
    }
#endif
    return result;
}

inline Ciphertext MultiplyByMonomial(const Ring &ring,
                                     const Ciphertext &ciphertext,
                                     const std::size_t exponent)
{
    ring.requireSize(ciphertext.a);
    ring.requireSize(ciphertext.b);
    const auto monomial = ring.monomialNTT(exponent);
    Ciphertext result;
    result.a.resize(ring.degree());
    result.b.resize(ring.degree());
#ifdef USE_HEXL
    intel::hexl::EltwiseMultMod(result.a.data(), monomial.data(),
                                ciphertext.a.data(), ring.degree(),
                                ring.modulus(), 1);
    intel::hexl::EltwiseMultMod(result.b.data(), monomial.data(),
                                ciphertext.b.data(), ring.degree(),
                                ring.modulus(), 1);
#else
    for (std::size_t i = 0; i < ring.degree(); i++) {
        result.a[i] = modular_ntt::multiply(monomial[i], ciphertext.a[i],
                                             ring.modulus());
        result.b[i] = modular_ntt::multiply(monomial[i], ciphertext.b[i],
                                             ring.modulus());
    }
#endif
    return result;
}

inline void validate(const Ring &ring, const Secret &secret)
{
    ring.requireSize(secret.slots);
    ring.requireSize(secret.quadratic_u);
    ring.requireSize(secret.quadratic_v);
    const std::uint64_t modulus = ring.modulus();
    for (std::size_t i = 0; i < ring.degree(); i++) {
        if (secret.slots[i] >= modulus || secret.quadratic_u[i] >= modulus ||
            secret.quadratic_v[i] >= modulus)
            throw std::invalid_argument("Binary-NTT secret is not reduced");
        const std::uint64_t square = modular_ntt::multiply(
            secret.slots[i], secret.slots[i], modulus);
        const std::uint64_t hinted = modular_ntt::add(
            modular_ntt::multiply(secret.quadratic_u[i], secret.slots[i],
                                  modulus),
            secret.quadratic_v[i], modulus);
        if (square != hinted)
            throw std::invalid_argument("Binary-NTT quadratic hint is invalid");
    }
}

template <class URBG>
Secret BinaryNTTSecretGen(const Ring &ring, URBG &engine)
{
    Secret result;
    result.slots.resize(ring.degree());
    result.quadratic_u.assign(ring.degree(), 1);
    result.quadratic_v.assign(ring.degree(), 0);
    std::bernoulli_distribution bit;
    for (std::uint64_t &slot : result.slots) slot = bit(engine) ? 1 : 0;
    return result;
}

// Algorithm 3 may switch from a Binary-NTT secret to a coefficient-small
// quadratic-hint secret before modulus switching. Publishing v=s^2 (u=0)
// gives the required s^2=u*s+v relation without relinearization.
template <class URBG>
Secret QuadraticHintSmallSecretGen(const Ring &ring,
                                   const std::uint32_t coefficient_bound,
                                   URBG &engine)
{
    if (coefficient_bound == 0)
        throw std::invalid_argument("small-secret bound must be positive");
    std::uniform_int_distribution<std::int64_t> coefficient(
        -static_cast<std::int64_t>(coefficient_bound),
        static_cast<std::int64_t>(coefficient_bound));
    std::vector<std::uint64_t> secret_coefficients(ring.degree());
    for (std::uint64_t &value : secret_coefficients) {
        const std::int64_t sampled = coefficient(engine);
        value = sampled >= 0
                    ? static_cast<std::uint64_t>(sampled)
                    : modular_ntt::negate(
                          static_cast<std::uint64_t>(-sampled), ring.modulus());
    }
    ring.forward(secret_coefficients);
    Secret result;
    result.slots = std::move(secret_coefficients);
    result.quadratic_u.assign(ring.degree(), 0);
    result.quadratic_v.resize(ring.degree());
    for (std::size_t i = 0; i < ring.degree(); i++)
        result.quadratic_v[i] = modular_ntt::multiply(
            result.slots[i], result.slots[i], ring.modulus());
    return result;
}

inline Secret QuadraticHintSmallSecretFromCoefficients(
    const Ring &ring, const std::span<const std::int64_t> coefficients)
{
    if (coefficients.size() != ring.degree())
        throw std::invalid_argument("small secret has wrong coefficient count");
    std::vector<std::uint64_t> values(ring.degree());
    for (std::size_t i = 0; i < ring.degree(); i++)
        values[i] = coefficients[i] >= 0
                        ? static_cast<std::uint64_t>(coefficients[i]) %
                              ring.modulus()
                        : modular_ntt::negate(
                              static_cast<std::uint64_t>(-coefficients[i]) %
                                  ring.modulus(),
                              ring.modulus());
    ring.forward(values);
    Secret result;
    result.slots = std::move(values);
    result.quadratic_u.assign(ring.degree(), 0);
    result.quadratic_v.resize(ring.degree());
    for (std::size_t i = 0; i < ring.degree(); i++)
        result.quadratic_v[i] = modular_ntt::multiply(
            result.slots[i], result.slots[i], ring.modulus());
    return result;
}

template <class URBG>
Secret QuadraticHintSmallSecretWithRandomHint(
    const Ring &ring, const std::span<const std::int64_t> coefficients,
    URBG &engine)
{
    Secret result = QuadraticHintSmallSecretFromCoefficients(ring, coefficients);
    std::uniform_int_distribution<std::uint64_t> uniform(0, ring.modulus() - 1);
    std::vector<std::uint64_t> u_coefficients(ring.degree());
    for (std::uint64_t &coefficient : u_coefficients)
        coefficient = uniform(engine);
    ring.forward(u_coefficients);
    result.quadratic_u = std::move(u_coefficients);
    for (std::size_t i = 0; i < ring.degree(); i++)
        result.quadratic_v[i] = modular_ntt::subtract(
            modular_ntt::multiply(result.slots[i], result.slots[i],
                                  ring.modulus()),
            modular_ntt::multiply(result.quadratic_u[i], result.slots[i],
                                  ring.modulus()),
            ring.modulus());
    return result;
}

inline Ciphertext TrivialEncryptNTT(const Ring &ring,
                                    const std::vector<std::uint64_t> &message)
{
    ring.requireSize(message);
    Ciphertext result;
    result.a.assign(ring.degree(), 0);
    result.b = message;
    return result;
}

template <class URBG>
BootstrappingKey BinaryNTTBootstrappingKeyGen(
    const Ring &ring, const Secret &secret,
    const std::span<const std::uint8_t> lwe_secret_bits,
    const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, secret);
    BootstrappingKey result;
    result.encrypted_lwe_bits.reserve(lwe_secret_bits.size());
    for (const std::uint8_t bit : lwe_secret_bits) {
        if (bit > 1)
            throw std::invalid_argument("Algorithm 3 requires binary LWE bits");
        std::vector<std::uint64_t> message(ring.degree());
        message[0] = bit;
        ring.forward(message);
        result.encrypted_lwe_bits.push_back(
            EncryptNTT(ring, secret, message, coefficient_noise_stddev,
                       engine));
    }
    return result;
}

template <class URBG>
KeySwitchKey BinaryNTTKeySwitchKeyGen(
    const Ring &ring, const Secret &old_secret, const Secret &new_secret,
    const std::uint32_t basebit, const double coefficient_noise_stddev,
    URBG &engine)
{
    validate(ring, old_secret);
    validate(ring, new_secret);
    if (basebit == 0 || basebit > 16)
        throw std::invalid_argument("Binary-NTT key-switch base is unsupported");
    const std::uint32_t digit_count = (64 + basebit - 1) / basebit;
    const std::uint64_t base = std::uint64_t{1} << basebit;
    KeySwitchKey result;
    result.basebit = basebit;
    result.digits.reserve(digit_count);
    std::uint64_t scale = 1;
    for (std::uint32_t digit = 0; digit < digit_count; digit++) {
        std::vector<std::uint64_t> message(ring.degree());
        for (std::size_t i = 0; i < ring.degree(); i++)
            message[i] = modular_ntt::multiply(scale, old_secret.slots[i],
                                                ring.modulus());
        result.digits.push_back(EncryptNTT(ring, new_secret, message,
                                           coefficient_noise_stddev, engine));
        scale = modular_ntt::multiply(scale, base, ring.modulus());
    }
    return result;
}

template <class URBG>
CoefficientKeySwitchKey CoefficientKeySwitchKeyGen(
    const Ring &ring, const Secret &old_secret, const Secret &new_secret,
    const std::uint32_t basebit, const double coefficient_noise_stddev,
    URBG &engine)
{
    validate(ring, old_secret);
    validate(ring, new_secret);
    if (basebit == 0 || basebit > 32)
        throw std::invalid_argument(
            "coefficient key-switch base is unsupported");
    const std::uint32_t digit_count = (64 + basebit - 1) / basebit;
    const std::uint64_t base = std::uint64_t{1} << basebit;
    CoefficientKeySwitchKey result;
    result.basebit = basebit;
    result.entries.resize(ring.degree());
    for (std::size_t coefficient = 0; coefficient < ring.degree(); coefficient++) {
        result.entries[coefficient].reserve(digit_count);
        std::uint64_t scale = 1;
        for (std::uint32_t digit = 0; digit < digit_count; digit++) {
            const std::vector<std::uint64_t> monomial =
                ring.monomialNTT(coefficient);
            std::vector<std::uint64_t> message(ring.degree());
            for (std::size_t slot = 0; slot < ring.degree(); slot++)
                message[slot] = modular_ntt::multiply(
                    modular_ntt::multiply(scale, monomial[slot],
                                          ring.modulus()),
                    old_secret.slots[slot], ring.modulus());
            result.entries[coefficient].push_back(EncryptNTT(
                ring, new_secret, message, coefficient_noise_stddev, engine));
            scale = modular_ntt::multiply(scale, base, ring.modulus());
        }
    }
    return result;
}

inline Ciphertext BinaryNTTKeySwitch(const Ring &ring, const KeySwitchKey &key,
                                     const Ciphertext &input)
{
    if (key.basebit == 0 || key.digits.empty())
        throw std::invalid_argument("Binary-NTT key-switch key is empty");
    ring.requireSize(input.a);
    ring.requireSize(input.b);
    for (const Ciphertext &digit : key.digits) {
        ring.requireSize(digit.a);
        ring.requireSize(digit.b);
    }
    const std::uint64_t modulus = ring.modulus();
    const std::uint64_t mask = (std::uint64_t{1} << key.basebit) - 1;
    Ciphertext result;
    result.a.assign(ring.degree(), 0);
    result.b = input.b;
    for (std::size_t i = 0; i < ring.degree(); i++) {
        std::uint64_t value = input.a[i];
        for (std::size_t digit = 0; digit < key.digits.size(); digit++) {
            const std::uint64_t coefficient = value & mask;
            value >>= key.basebit;
            if (coefficient == 0) continue;
            result.a[i] = modular_ntt::subtract(
                result.a[i],
                modular_ntt::multiply(coefficient, key.digits[digit].a[i],
                                      modulus),
                modulus);
            result.b[i] = modular_ntt::subtract(
                result.b[i],
                modular_ntt::multiply(coefficient, key.digits[digit].b[i],
                                      modulus),
                modulus);
        }
    }
    return result;
}

inline Ciphertext CoefficientKeySwitch(const Ring &ring,
                                       const CoefficientKeySwitchKey &key,
                                       const Ciphertext &input)
{
    if (key.basebit == 0 || key.entries.size() != ring.degree())
        throw std::invalid_argument("invalid coefficient key-switch key");
    ring.requireSize(input.a);
    ring.requireSize(input.b);
    const std::uint64_t mask = key.basebit == 64
        ? std::numeric_limits<std::uint64_t>::max()
        : (std::uint64_t{1} << key.basebit) - 1;
    std::vector<std::uint64_t> coefficients = input.a;
    ring.inverse(coefficients);
    Ciphertext result;
    result.a.assign(ring.degree(), 0);
    result.b = input.b;
    for (std::size_t coefficient = 0; coefficient < ring.degree(); coefficient++) {
        if (key.entries[coefficient].empty())
            throw std::invalid_argument("empty coefficient key-switch entry");
        std::uint64_t value = coefficients[coefficient];
        for (std::size_t digit = 0; digit < key.entries[coefficient].size();
             digit++) {
            const std::uint64_t decomposition = value & mask;
            value >>= key.basebit;
            if (decomposition == 0) continue;
            const Ciphertext &entry = key.entries[coefficient][digit];
            for (std::size_t slot = 0; slot < ring.degree(); slot++) {
                result.a[slot] = modular_ntt::subtract(
                    result.a[slot], modular_ntt::multiply(
                        decomposition, entry.a[slot], ring.modulus()),
                    ring.modulus());
                result.b[slot] = modular_ntt::subtract(
                    result.b[slot], modular_ntt::multiply(
                        decomposition, entry.b[slot], ring.modulus()),
                    ring.modulus());
            }
        }
    }
    return result;
}

inline std::uint64_t RoundModulusSwitchCoefficient(
    const __int128 value, const std::uint64_t source_modulus,
    const std::uint64_t target_modulus)
{
    const bool negative = value < 0;
    const unsigned __int128 magnitude = negative
        ? static_cast<unsigned __int128>(-(value + 1)) + 1
        : static_cast<unsigned __int128>(value);
    const std::uint64_t rounded = static_cast<std::uint64_t>(
        (magnitude * target_modulus + source_modulus / 2) / source_modulus);
    const std::uint64_t reduced = rounded % target_modulus;
    return negative ? modular_ntt::negate(reduced, target_modulus) : reduced;
}

// Coefficient-domain modulus switch.  It is valid when the caller has first
// switched to a secret with the same small coefficient representation in both
// rings; Binary-NTT slot secrets alone cannot be modulus-switched this way.
inline Ciphertext ModulusSwitch(const Ring &source, const Ring &target,
                                const Ciphertext &input)
{
    if (source.degree() != target.degree())
        throw std::invalid_argument("modulus switch changes ring degree");
    source.requireSize(input.a);
    source.requireSize(input.b);
    Ciphertext result;
    result.a = input.a;
    result.b = input.b;
    source.inverse(result.a);
    source.inverse(result.b);
    for (std::size_t i = 0; i < source.degree(); i++) {
        result.a[i] = RoundModulusSwitchCoefficient(
            modular_ntt::centeredResidue(result.a[i], source.modulus()),
            source.modulus(), target.modulus());
        result.b[i] = RoundModulusSwitchCoefficient(
            modular_ntt::centeredResidue(result.b[i], source.modulus()),
            source.modulus(), target.modulus());
    }
    target.forward(result.a);
    target.forward(result.b);
    return result;
}

inline Ciphertext Algorithm3ModulusBoundary(
    const Ring &source_ring, const Secret &binary_ntt_secret,
    const KeySwitchKey &switch_key, const Secret &small_secret_source,
    const Ring &target_ring, const Secret &small_secret_target,
    const Ciphertext &input)
{
    validate(source_ring, binary_ntt_secret);
    validate(source_ring, small_secret_source);
    validate(target_ring, small_secret_target);
    const Ciphertext switched =
        BinaryNTTKeySwitch(source_ring, switch_key, input);
    const Ciphertext reduced =
        ModulusSwitch(source_ring, target_ring, switched);
    return reduced;
}

inline Ciphertext Algorithm3ModulusBoundary(
    const Ring &source_ring, const Secret &binary_ntt_secret,
    const CoefficientKeySwitchKey &switch_key,
    const Secret &small_secret_source, const Ring &target_ring,
    const Secret &small_secret_target, const Ciphertext &input)
{
    validate(source_ring, binary_ntt_secret);
    validate(source_ring, small_secret_source);
    validate(target_ring, small_secret_target);
    return ModulusSwitch(source_ring, target_ring,
                         CoefficientKeySwitch(source_ring, switch_key, input));
}

// Standard sample extraction for b-a*s RLWE convention. The output secret is
// the coefficient representation of the small RLWE secret.
inline LWECiphertext SampleExtract(const Ring &ring, const Ciphertext &input)
{
    ring.requireSize(input.a);
    ring.requireSize(input.b);
    LWECiphertext result;
    result.a.resize(ring.degree());
    result.modulus = ring.modulus();
    std::vector<std::uint64_t> a = input.a;
    std::vector<std::uint64_t> b = input.b;
    ring.inverse(a);
    ring.inverse(b);
    result.a[0] = a[0];
    for (std::size_t i = 1; i < ring.degree(); i++)
        result.a[i] = modular_ntt::negate(a[ring.degree() - i], ring.modulus());
    result.b = b[0];
    return result;
}

template <class URBG>
LWECiphertext LWEEncrypt(const std::uint64_t message,
                         const std::span<const std::uint8_t> secret,
                         const std::uint64_t modulus,
                         const double noise_stddev, URBG &engine)
{
    if (secret.empty() || modulus < 2 || noise_stddev < 0)
        throw std::invalid_argument("invalid LWE encryption parameters");
    std::uniform_int_distribution<std::uint64_t> uniform(0, modulus - 1);
    std::normal_distribution<double> gaussian(0.0, noise_stddev);
    LWECiphertext result;
    result.a.resize(secret.size());
    result.modulus = modulus;
    result.b = message % modulus;
    for (std::size_t i = 0; i < secret.size(); i++) {
        if (secret[i] > 1)
            throw std::invalid_argument("LWE key switch expects binary target key");
        result.a[i] = uniform(engine);
        result.b = modular_ntt::add(
            result.b, secret[i] == 0 ? 0 : result.a[i], modulus);
    }
    const std::int64_t error =
        static_cast<std::int64_t>(std::llround(gaussian(engine)));
    const std::uint64_t encoded_error = error >= 0
        ? static_cast<std::uint64_t>(error) % modulus
        : modular_ntt::negate(static_cast<std::uint64_t>(-error) % modulus,
                              modulus);
    result.b = modular_ntt::add(result.b, encoded_error, modulus);
    return result;
}

template <class URBG>
LWEKeySwitchKey LWEKeySwitchKeyGen(
    const std::span<const std::int64_t> source_secret,
    const std::span<const std::uint8_t> target_secret,
    const std::uint64_t modulus, const std::uint32_t basebit,
    const double noise_stddev, URBG &engine,
    const std::uint32_t requested_digits = 0,
    const bool balanced = false)
{
    if (source_secret.empty() || basebit == 0 || basebit > 16)
        throw std::invalid_argument("invalid LWE key-switch parameters");
    const std::uint32_t digits = requested_digits == 0
        ? (64 + basebit - 1) / basebit
        : requested_digits;
    LWEKeySwitchKey result;
    result.basebit = basebit;
    result.balanced = balanced;
    result.entries.resize(source_secret.size());
    const std::uint64_t base = std::uint64_t{1} << basebit;
    for (std::size_t i = 0; i < source_secret.size(); i++) {
        result.entries[i].reserve(digits);
        std::uint64_t digit_scale = 1;
        for (std::uint32_t digit = 0; digit < digits; digit++) {
            const std::int64_t raw = source_secret[i];
            const std::uint64_t message = raw >= 0
                ? modular_ntt::multiply(static_cast<std::uint64_t>(raw),
                                         digit_scale, modulus)
                : modular_ntt::negate(
                      modular_ntt::multiply(static_cast<std::uint64_t>(-raw),
                                             digit_scale, modulus),
                      modulus);
            result.entries[i].push_back(
                LWEEncrypt(message, target_secret, modulus, noise_stddev,
                           engine));
            digit_scale = modular_ntt::multiply(digit_scale, base, modulus);
        }
    }
    return result;
}

inline LWECiphertext LWEKeySwitch(const LWECiphertext &input,
                                  const LWEKeySwitchKey &key)
{
    if (input.modulus < 2 || key.basebit == 0 ||
        key.entries.size() != input.a.size())
        throw std::invalid_argument("invalid LWE key-switch input");
    const std::uint64_t modulus = input.modulus;
    const std::uint64_t mask = (std::uint64_t{1} << key.basebit) - 1;
    const std::uint64_t base = std::uint64_t{1} << key.basebit;
    const auto next_digit = [&](std::uint64_t &value) {
        std::int64_t digit = static_cast<std::int64_t>(value & mask);
        if (key.balanced && digit >= static_cast<std::int64_t>(base / 2)) {
            digit -= static_cast<std::int64_t>(base);
            value = (value + base) >> key.basebit;
        }
        else
            value >>= key.basebit;
        return digit;
    };
    const bool power_of_two_modulus =
        (modulus & (modulus - 1)) == 0;
    const std::uint64_t modulus_mask = modulus - 1;
    const auto multiply_mod = [&](const std::uint64_t left,
                                  const std::uint64_t right) {
        return power_of_two_modulus
            ? (left * right) & modulus_mask
            : modular_ntt::multiply(left, right, modulus);
    };
    const std::size_t target_dimension = key.entries.front().front().a.size();
    LWECiphertext result;
    result.a.assign(target_dimension, 0);
    result.b = input.b;
    result.modulus = modulus;
#ifdef _OPENMP
    const int thread_count = omp_get_max_threads();
    if (thread_count > 1 && input.a.size() > 1) {
        std::vector<std::vector<std::uint64_t>> partial_a(
            thread_count, std::vector<std::uint64_t>(target_dimension));
        std::vector<std::uint64_t> partial_b(thread_count, 0);
#pragma omp parallel
        {
            const int thread = omp_get_thread_num();
#pragma omp for
            for (std::int64_t signed_i = 0;
                 signed_i < static_cast<std::int64_t>(input.a.size());
                 signed_i++) {
                const std::size_t i = static_cast<std::size_t>(signed_i);
                std::uint64_t value = input.a[i];
                for (std::size_t digit = 0;
                     digit < key.entries[i].size(); digit++) {
                    const std::int64_t digit_value = next_digit(value);
                    if (digit_value == 0) continue;
                    const std::uint64_t coefficient = digit_value >= 0
                        ? static_cast<std::uint64_t>(digit_value)
                        : static_cast<std::uint64_t>(-digit_value);
                    const LWECiphertext &entry = key.entries[i][digit];
                    for (std::size_t j = 0; j < target_dimension; j++) {
                        const std::uint64_t product =
                            multiply_mod(coefficient, entry.a[j]);
                        partial_a[thread][j] = digit_value >= 0
                            ? modular_ntt::subtract(partial_a[thread][j],
                                                    product, modulus)
                            : modular_ntt::add(partial_a[thread][j], product,
                                               modulus);
                    }
                    const std::uint64_t b_product =
                        multiply_mod(coefficient, entry.b);
                    partial_b[thread] = digit_value >= 0
                        ? modular_ntt::subtract(partial_b[thread], b_product,
                                                modulus)
                        : modular_ntt::add(partial_b[thread], b_product,
                                           modulus);
                }
            }
        }
        for (int thread = 0; thread < thread_count; thread++) {
            for (std::size_t j = 0; j < target_dimension; j++)
                result.a[j] = modular_ntt::add(
                    result.a[j], partial_a[thread][j], modulus);
            result.b = modular_ntt::add(result.b, partial_b[thread], modulus);
        }
        return result;
    }
#endif
    for (std::size_t i = 0; i < input.a.size(); i++) {
        std::uint64_t value = input.a[i];
        for (std::size_t digit = 0; digit < key.entries[i].size(); digit++) {
            const std::int64_t digit_value = next_digit(value);
            if (digit_value == 0) continue;
            const std::uint64_t coefficient = digit_value >= 0
                ? static_cast<std::uint64_t>(digit_value)
                : static_cast<std::uint64_t>(-digit_value);
            const LWECiphertext &entry = key.entries[i][digit];
            for (std::size_t j = 0; j < target_dimension; j++) {
                const std::uint64_t product =
                    multiply_mod(coefficient, entry.a[j]);
                result.a[j] = digit_value >= 0
                    ? modular_ntt::subtract(result.a[j], product, modulus)
                    : modular_ntt::add(result.a[j], product, modulus);
            }
            const std::uint64_t b_product =
                multiply_mod(coefficient, entry.b);
            result.b = digit_value >= 0
                ? modular_ntt::subtract(result.b, b_product, modulus)
                : modular_ntt::add(result.b, b_product, modulus);
        }
    }
    return result;
}

inline LWECiphertext ModulusSwitchCiphertext(
    const LWECiphertext &input, const std::uint64_t target_modulus)
{
    if (target_modulus < 2)
        throw std::invalid_argument("LWE modulus switch requires modulus >= 2");
    LWECiphertext result;
    result.modulus = target_modulus;
    result.a.resize(input.a.size());
    for (std::size_t i = 0; i < input.a.size(); i++)
        result.a[i] = RoundModulusSwitchCoefficient(
            modular_ntt::centeredResidue(input.a[i], input.modulus),
            input.modulus, target_modulus);
    result.b = RoundModulusSwitchCoefficient(
        modular_ntt::centeredResidue(input.b, input.modulus), input.modulus,
        target_modulus);
    return result;
}

inline LWEPhase ModulusSwitch(const LWECiphertext &input,
                              const std::uint32_t target_modulus)
{
    const LWECiphertext switched =
        ModulusSwitchCiphertext(input, target_modulus);
    LWEPhase result;
    result.modulus = target_modulus;
    result.a.reserve(switched.a.size());
    for (const std::uint64_t value : switched.a)
        result.a.push_back(static_cast<std::uint32_t>(value));
    result.b = static_cast<std::uint32_t>(switched.b);
    return result;
}

inline std::uint32_t LWEDecryptPhase(
    const LWEPhase &ciphertext, const std::span<const std::uint8_t> secret)
{
    if (ciphertext.a.size() != secret.size() || ciphertext.modulus < 2)
        throw std::invalid_argument("invalid LWE decryption input");
    std::uint32_t phase = ciphertext.b % ciphertext.modulus;
    for (std::size_t i = 0; i < secret.size(); i++) {
        if (secret[i] > 1)
            throw std::invalid_argument("LWE secret is not binary");
        if (secret[i] != 0)
            phase = static_cast<std::uint32_t>(
                modular_ntt::subtract(phase, ciphertext.a[i],
                                      ciphertext.modulus));
    }
    return phase;
}

template <class URBG>
Ciphertext EncryptNTT(const Ring &ring, const Secret &secret,
                      const std::vector<std::uint64_t> &message,
                      const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, secret);
    ring.requireSize(message);
    if (coefficient_noise_stddev < 0)
        throw std::invalid_argument("noise standard deviation must be nonnegative");

    const std::uint64_t modulus = ring.modulus();
    Ciphertext result;
    result.a.resize(ring.degree());
    result.b.resize(ring.degree());
    std::uniform_int_distribution<std::uint64_t> uniform(0, modulus - 1);
    std::normal_distribution<double> gaussian(0.0, coefficient_noise_stddev);
    std::vector<std::uint64_t> error(ring.degree());
    for (std::size_t i = 0; i < ring.degree(); i++) {
        result.a[i] = uniform(engine);
        const std::int64_t sampled =
            static_cast<std::int64_t>(std::llround(gaussian(engine)));
        error[i] = sampled >= 0
                       ? static_cast<std::uint64_t>(sampled) % modulus
                       : modular_ntt::negate(
                             static_cast<std::uint64_t>(-sampled) % modulus,
                             modulus);
    }
    ring.forward(error);
    for (std::size_t i = 0; i < ring.degree(); i++) {
        const std::uint64_t masked = modular_ntt::multiply(
            result.a[i], secret.slots[i], modulus);
        result.b[i] = modular_ntt::add(
            message[i], modular_ntt::add(masked, error[i], modulus), modulus);
    }
    return result;
}

inline std::vector<std::uint64_t> DecryptNTT(const Ring &ring,
                                             const Secret &secret,
                                             const Ciphertext &ciphertext)
{
    validate(ring, secret);
    ring.requireSize(ciphertext.a);
    ring.requireSize(ciphertext.b);
    std::vector<std::uint64_t> result(ring.degree());
    for (std::size_t i = 0; i < ring.degree(); i++)
        result[i] = modular_ntt::subtract(
            ciphertext.b[i],
            modular_ntt::multiply(ciphertext.a[i], secret.slots[i],
                                  ring.modulus()),
            ring.modulus());
    return result;
}

// Algorithm 3, Mult.  No decomposition or relinearization occurs here: all
// operations are pointwise in the retained NTT representation.
inline Ciphertext RelinearizationFreeMultiply(const Ring &ring,
                                              const Secret &secret,
                                              const Ciphertext &left,
                                              const Ciphertext &right)
{
    validate(ring, secret);
    ring.requireSize(left.a);
    ring.requireSize(left.b);
    ring.requireSize(right.a);
    ring.requireSize(right.b);

    const std::uint64_t modulus = ring.modulus();
    Ciphertext result;
    result.a.resize(ring.degree());
    result.b.resize(ring.degree());
#ifdef USE_HEXL
    // Keep all six tensor/hint products in HEXL's vector kernels.  These are
    // the hot operations in every level of the Algorithm-3 product tree.
    std::vector<std::uint64_t> x1(ring.degree());
    std::vector<std::uint64_t> scratch(ring.degree());
    intel::hexl::EltwiseMultMod(x1.data(), left.a.data(), right.a.data(),
                                ring.degree(), modulus, 1);
    intel::hexl::EltwiseMultMod(result.a.data(), left.a.data(), right.b.data(),
                                ring.degree(), modulus, 1);
    intel::hexl::EltwiseMultMod(scratch.data(), left.b.data(), right.a.data(),
                                ring.degree(), modulus, 1);
    intel::hexl::EltwiseAddMod(result.a.data(), result.a.data(), scratch.data(),
                               ring.degree(), modulus);
    intel::hexl::EltwiseMultMod(scratch.data(), x1.data(),
                                secret.quadratic_u.data(), ring.degree(),
                                modulus, 1);
    intel::hexl::EltwiseSubMod(result.a.data(), result.a.data(), scratch.data(),
                               ring.degree(), modulus);
    intel::hexl::EltwiseMultMod(result.b.data(), left.b.data(), right.b.data(),
                                ring.degree(), modulus, 1);
    intel::hexl::EltwiseMultMod(scratch.data(), x1.data(),
                                secret.quadratic_v.data(), ring.degree(),
                                modulus, 1);
    intel::hexl::EltwiseAddMod(result.b.data(), result.b.data(), scratch.data(),
                               ring.degree(), modulus);
#else
    for (std::size_t i = 0; i < ring.degree(); i++) {
        const std::uint64_t x1 = modular_ntt::multiply(left.a[i], right.a[i],
                                                        modulus);
        const std::uint64_t x2 = modular_ntt::multiply(left.a[i], right.b[i],
                                                        modulus);
        const std::uint64_t x3 = modular_ntt::multiply(left.b[i], right.a[i],
                                                        modulus);
        const std::uint64_t x4 = modular_ntt::multiply(left.b[i], right.b[i],
                                                        modulus);
        // TFHEpp's low-depth ciphertext convention is b-a*s=m.  Substitute
        // s^2=u*s+v into (b1-a1*s)(b2-a2*s): the resulting two components
        // are a'=x2+x3-x1*u and b'=x4+x1*v.  This is the sign-adjusted form
        // of Algorithm 3's Mult routine.
        result.a[i] = modular_ntt::subtract(
            modular_ntt::add(x2, x3, modulus),
            modular_ntt::multiply(x1, secret.quadratic_u[i], modulus),
            modulus);
        result.b[i] = modular_ntt::add(
            x4, modular_ntt::multiply(x1, secret.quadratic_v[i], modulus),
            modulus);
    }
#endif
    return result;
}

// Multiply a collection in a balanced tree.  Every level is independent and
// can therefore be dispatched in parallel.  No NTT or inverse NTT is called
// after the input ciphertexts enter this routine.
inline Ciphertext ProductTree(const Ring &ring, const Secret &secret,
                              std::vector<Ciphertext> inputs,
                              ProductTreeStats *stats = nullptr)
{
    if (inputs.empty())
        throw std::invalid_argument("Binary-NTT product tree requires an input");
    for (const Ciphertext &input : inputs) {
        ring.requireSize(input.a);
        ring.requireSize(input.b);
    }

    ProductTreeStats local_stats;
    while (inputs.size() > 1) {
        const std::size_t pairs = inputs.size() / 2;
        std::vector<Ciphertext> next((inputs.size() + 1) / 2);
#pragma omp parallel for if (pairs > 1)
        for (std::int64_t pair = 0; pair < static_cast<std::int64_t>(pairs);
             pair++)
            next[static_cast<std::size_t>(pair)] =
                RelinearizationFreeMultiply(
                    ring, secret, inputs[2 * static_cast<std::size_t>(pair)],
                    inputs[2 * static_cast<std::size_t>(pair) + 1]);
        if ((inputs.size() & 1U) != 0)
            next.back() = std::move(inputs.back());
        local_stats.pointwise_multiplication_layers++;
        local_stats.pointwise_ciphertext_multiplications +=
            static_cast<std::uint32_t>(pairs);
        inputs = std::move(next);
    }
    if (stats != nullptr) *stats = local_stats;
    return std::move(inputs.front());
}

// Homomorphically form 1 + (X^a - 1)s for a ciphertext encrypting s.  This
// is the low-depth analogue of the prepared bootstrapping-key entries in
// Algorithm 3, line 23.
inline Ciphertext MonomialFactor(const Ring &ring, const Ciphertext &bit_ct,
                                 const std::size_t exponent)
{
    ring.requireSize(bit_ct.a);
    ring.requireSize(bit_ct.b);
    const std::uint64_t modulus = ring.modulus();
    std::vector<std::uint64_t> monomial = ring.monomialNTT(exponent);
    Ciphertext result;
    result.a.resize(ring.degree());
    result.b.resize(ring.degree());
    for (std::size_t i = 0; i < ring.degree(); i++) {
        const std::uint64_t multiplier =
            modular_ntt::subtract(monomial[i], 1, modulus);
        result.a[i] = modular_ntt::multiply(multiplier, bit_ct.a[i], modulus);
        result.b[i] = modular_ntt::add(
            modular_ntt::multiply(multiplier, bit_ct.b[i], modulus), 1,
            modulus);
    }
    return result;
}

// Algorithm 3's low-depth blind-rotation core after PBC/structured-key
// preparation: `encrypted_bits[i]` encrypts the i-th selected LWE-key bit and
// exponents[i] is the corresponding (already modulus-switched) LWE mask.
// The initial accumulator and all prepared factors are multiplied in one
// balanced tree, rather than through the sequential RGSW accumulator chain.
inline Ciphertext LowDepthBlindRotate(
    const Ring &ring, const Secret &secret, const Ciphertext &initial,
    const std::span<const Ciphertext> encrypted_bits,
    const std::span<const std::size_t> exponents,
    ProductTreeStats *stats = nullptr)
{
    if (encrypted_bits.size() != exponents.size())
        throw std::invalid_argument(
            "low-depth blind rotation has mismatched factors and exponents");
    std::vector<Ciphertext> factors;
    factors.reserve(encrypted_bits.size() + 1);
    factors.push_back(initial);
    for (std::size_t i = 0; i < encrypted_bits.size(); i++)
        factors.push_back(MonomialFactor(ring, encrypted_bits[i], exponents[i]));
    return ProductTree(ring, secret, std::move(factors), stats);
}

inline Ciphertext PBCPreparedFactor(
    const Ring &ring, const PBCSchedule &schedule,
    const PBCBootstrappingKey &bootstrapping_key,
    const std::uint32_t bucket, const std::span<const std::uint32_t> lwe_a,
    const std::size_t exponent_scale)
{
    if (bucket >= schedule.bucket_count || lwe_a.size() != schedule.source_dimension ||
        bootstrapping_key.entries.size() != schedule.bucket_count ||
        bootstrapping_key.dummies.size() != schedule.bucket_count ||
        bootstrapping_key.entries[bucket].size() !=
            schedule.bucket_indices[bucket].size())
        throw std::invalid_argument("invalid PBC prepared-factor inputs");
    Ciphertext result = bootstrapping_key.dummies[bucket];
    for (std::size_t slot = 0;
         slot < schedule.bucket_indices[bucket].size(); slot++) {
        const std::uint32_t source = schedule.bucket_indices[bucket][slot];
        result = AddCiphertexts(
            ring, result,
            MultiplyByMonomial(ring, bootstrapping_key.entries[bucket][slot],
                               lwe_a[source] * exponent_scale));
    }
    return result;
}

inline Ciphertext Algorithm3PBCBlindRotate(
    const Ring &ring, const Secret &secret, const PBCSchedule &schedule,
    const PBCBootstrappingKey &bootstrapping_key, const LWEPhase &lwe,
    const std::vector<std::uint64_t> &test_vector,
    ProductTreeStats *stats = nullptr)
{
    if (lwe.modulus == 0 || (2 * ring.degree()) % lwe.modulus != 0 ||
        lwe.a.size() != schedule.source_dimension || lwe.b >= lwe.modulus)
        throw std::invalid_argument("invalid Algorithm 3 PBC LWE input");
    ring.requireSize(test_vector);
    const std::size_t exponent_scale = 2 * ring.degree() / lwe.modulus;
    const std::size_t initial_exponent =
        (lwe.modulus - lwe.b) % lwe.modulus;
    Ciphertext initial = TrivialEncryptNTT(ring, test_vector);
    initial = MultiplyByMonomial(ring, initial,
                                 initial_exponent * exponent_scale);
    // Keep exactly `bucket_count` leaves: the accumulator is public, so it
    // can be absorbed into one PBC factor.  Thus the paper's k=h+3=32
    // schedule fits precisely into windows 8 then 4.
    std::vector<Ciphertext> factors;
    factors.reserve(schedule.bucket_count);
    for (std::uint32_t bucket = 0; bucket < schedule.bucket_count; bucket++) {
        Ciphertext factor = PBCPreparedFactor(ring, schedule,
                                               bootstrapping_key, bucket,
                                               lwe.a, exponent_scale);
        if (bucket == 0)
            for (std::size_t i = 0; i < ring.degree(); i++) {
                factor.a[i] = modular_ntt::multiply(
                    factor.a[i], initial.b[i], ring.modulus());
                factor.b[i] = modular_ntt::multiply(
                    factor.b[i], initial.b[i], ring.modulus());
            }
        factors.push_back(std::move(factor));
    }
    return ProductTree(ring, secret, std::move(factors), stats);
}

// Evaluate Algorithm 3's window schedule.  A boundary callback is invoked on
// every group before the following window; a concrete deployment supplies the
// INTT -> key switch/modulus switch -> NTT transition there.  Keeping the
// callback explicit prevents an identity modulus change from being mistaken
// for the paper's noise-control step.
template <class Boundary>
Ciphertext ScheduledProductTree(const Ring &ring, const Secret &secret,
                                std::vector<Ciphertext> inputs,
                                const std::span<const std::uint32_t> windows,
                                Boundary &&boundary,
                                ModulusScheduleStats *stats = nullptr)
{
    if (inputs.empty() || windows.empty())
        throw std::invalid_argument(
            "scheduled product tree requires inputs and a schedule");
    ModulusScheduleStats local_stats;
    for (const std::uint32_t width : windows) {
        if (width < 2 || (width & (width - 1)) != 0)
            throw std::invalid_argument(
                "Algorithm 3 multiplication windows must be powers of two");
        std::vector<Ciphertext> next;
        next.reserve((inputs.size() + width - 1) / width);
        std::uint32_t stage_depth = 0;
        for (std::size_t first = 0; first < inputs.size(); first += width) {
            const std::size_t last = std::min<std::size_t>(first + width,
                                                            inputs.size());
            std::vector<Ciphertext> group;
            group.reserve(last - first);
            for (std::size_t i = first; i < last; i++)
                group.push_back(std::move(inputs[i]));
            ProductTreeStats group_stats;
            next.push_back(
                ProductTree(ring, secret, std::move(group), &group_stats));
            stage_depth = std::max(stage_depth,
                                   group_stats.pointwise_multiplication_layers);
            local_stats.pointwise_ciphertext_multiplications +=
                group_stats.pointwise_ciphertext_multiplications;
        }
        local_stats.pointwise_multiplication_layers += stage_depth;
        if (next.size() == 1) {
            if (stats != nullptr) *stats = local_stats;
            return std::move(next.front());
        }
        for (Ciphertext &ciphertext : next) boundary(ciphertext);
        local_stats.modulus_boundaries +=
            static_cast<std::uint32_t>(next.size());
        inputs = std::move(next);
    }
    if (inputs.size() != 1)
        throw std::invalid_argument(
            "Algorithm 3 schedule does not reduce to one ciphertext");
    if (stats != nullptr) *stats = local_stats;
    return std::move(inputs.front());
}

inline Ciphertext LowDepthBlindRotate(
    const Ring &ring, const Secret &secret, const Ciphertext &initial,
    const BootstrappingKey &bootstrapping_key,
    const std::span<const std::size_t> exponents,
    ProductTreeStats *stats = nullptr)
{
    return LowDepthBlindRotate(ring, secret, initial,
                               bootstrapping_key.encrypted_lwe_bits,
                               exponents, stats);
}

// Algorithm 3's blind-rotation product, specialized to a binary LWE secret.
// The LWE modulus must divide 2N so every LWE coefficient maps to an exact
// negacyclic-ring monomial.  The caller supplies the encrypted binary secret
// bits as `bootstrapping_key`; the LWE phase itself remains public input.
inline Ciphertext Algorithm3BlindRotate(
    const Ring &ring, const Secret &secret, const BootstrappingKey &bootstrapping_key,
    const LWEPhase &lwe, const std::vector<std::uint64_t> &test_vector,
    ProductTreeStats *stats = nullptr)
{
    if (lwe.modulus == 0 || (2 * ring.degree()) % lwe.modulus != 0)
        throw std::invalid_argument(
            "Algorithm 3 requires an LWE modulus dividing 2N");
    if (lwe.a.size() != bootstrapping_key.encrypted_lwe_bits.size())
        throw std::invalid_argument(
            "Algorithm 3 LWE and bootstrapping-key dimensions differ");
    ring.requireSize(test_vector);
    const std::size_t scale = 2 * ring.degree() / lwe.modulus;
    const std::size_t initial_exponent =
        (lwe.modulus - (lwe.b % lwe.modulus)) % lwe.modulus;
    std::vector<std::uint64_t> initial_message = test_vector;
    // Apply X^{-b}; all subsequent factors contribute X^{a_i*s_i}.
    const std::vector<std::uint64_t> initial_monomial =
        ring.monomialNTT(initial_exponent * scale);
    for (std::size_t i = 0; i < ring.degree(); i++)
        initial_message[i] = modular_ntt::multiply(
            initial_message[i], initial_monomial[i], ring.modulus());
    const Ciphertext initial = TrivialEncryptNTT(ring, initial_message);

    std::vector<std::size_t> exponents(lwe.a.size());
    for (std::size_t i = 0; i < lwe.a.size(); i++)
        exponents[i] = (lwe.a[i] % lwe.modulus) * scale;
    return LowDepthBlindRotate(ring, secret, initial, bootstrapping_key,
                               exponents, stats);
}

// The concrete Algorithm 3 parameters need a modulus above one machine word.
// RNSRing keeps the same Binary-NTT secret in every prime component and gives
// the product tree a two-prime (about 124-bit) execution path.  The primes are
// intentionally explicit; callers may choose a smaller product once a proved
// modulus/noise schedule is fixed.
class RNSRing {
public:
    RNSRing(const std::size_t degree,
            const std::span<const modular_ntt::PrimeModulus> moduli)
    {
        if (moduli.empty())
            throw std::invalid_argument("RNS ring requires a prime modulus");
        rings_.reserve(moduli.size());
        primes_.assign(moduli.begin(), moduli.end());
        for (const auto modulus : moduli) rings_.emplace_back(degree, modulus);
    }

    std::size_t degree() const { return rings_.front().degree(); }
    std::size_t levels() const { return rings_.size(); }
    const Ring &operator[](const std::size_t level) const
    {
        return rings_.at(level);
    }
    const modular_ntt::PrimeModulus &prime(const std::size_t level) const
    {
        return primes_.at(level);
    }

private:
    std::vector<Ring> rings_;
    std::vector<modular_ntt::PrimeModulus> primes_;
};

// Concrete N=8192 NTT primes for qh_ss_source_screened_rns. Their products
// have bit lengths 151, 69, 52, and 36 respectively. Each level is a prefix
// of the preceding RNS basis, as required by BGV modulus switching. Each prime is 1 modulo
// 2N, and the second field is a primitive root modulo that prime.
inline constexpr std::array<modular_ntt::PrimeModulus, 6>
    secure_q151_primes{{
        {34359754753ULL, 5},
        {65537ULL, 3},
        {147457ULL, 10},
        {268582913ULL, 3},
        {134250497ULL, 3},
        {134348801ULL, 3},
    }};
inline constexpr std::array<modular_ntt::PrimeModulus, 3>
    secure_q69_primes{{
        {34359754753ULL, 5},
        {65537ULL, 3},
        {147457ULL, 10},
    }};
inline constexpr std::array<modular_ntt::PrimeModulus, 2>
    secure_q52_primes{{{34359754753ULL, 5}, {65537ULL, 3}}};
inline constexpr std::array<modular_ntt::PrimeModulus, 1>
    secure_q36_primes{{{34359754753ULL, 5}}};


struct RNSSecret {
    std::vector<Secret> residues;
};

struct RNSCiphertext {
    std::vector<Ciphertext> residues;
};

struct RNSBootstrappingKey {
    std::vector<RNSCiphertext> encrypted_lwe_bits;
};

struct RNSPBCBootstrappingKey {
    std::vector<std::vector<RNSCiphertext>> entries;
    std::vector<RNSCiphertext> dummies;
};

#ifdef USE_BLAKE3
struct SeededRNSCiphertext {
    std::array<std::uint8_t, 32> a_seed{};
    std::vector<std::vector<std::uint64_t>> b_residues;
};

struct SeededRNSPBCBootstrappingKey {
    std::vector<std::vector<SeededRNSCiphertext>> entries;
    std::vector<SeededRNSCiphertext> dummies;
};
#endif

struct RNSKeySwitchKey {
    std::vector<KeySwitchKey> residues;
};

struct RNSCoefficientKeySwitchKey {
    std::vector<CoefficientKeySwitchKey> residues;
};

struct TwoStageStats {
    ProductTreeStats high_product;
    ProductTreeStats low_product;
    std::uint32_t modulus_boundaries = 0;
};

struct MultiStageStats {
    std::vector<ProductTreeStats> products;
    std::vector<std::uint32_t> modulus_boundaries;
};

template <class URBG>
RNSSecret BinaryNTTSecretGen(const RNSRing &ring, URBG &engine)
{
    std::vector<std::uint64_t> bits(ring.degree());
    std::bernoulli_distribution bit;
    for (std::uint64_t &slot : bits) slot = bit(engine) ? 1 : 0;

    RNSSecret result;
    result.residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++) {
        result.residues[level].slots = bits;
        result.residues[level].quadratic_u.assign(ring.degree(), 1);
        result.residues[level].quadratic_v.assign(ring.degree(), 0);
    }
    return result;
}

inline RNSSecret QuadraticHintSmallSecretFromCoefficients(
    const RNSRing &ring, const std::span<const std::int64_t> coefficients)
{
    RNSSecret result;
    result.residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++)
        result.residues[level] = QuadraticHintSmallSecretFromCoefficients(
            ring[level], coefficients);
    return result;
}

template <class URBG>
RNSSecret QuadraticHintSmallSecretGen(
    const RNSRing &ring, const std::span<const std::int64_t> coefficients,
    URBG &engine)
{
    RNSSecret result = QuadraticHintSmallSecretFromCoefficients(ring, coefficients);
    for (std::size_t level = 0; level < ring.levels(); level++) {
        std::uniform_int_distribution<std::uint64_t> uniform(
            0, ring[level].modulus() - 1);
        std::vector<std::uint64_t> u_coefficients(ring.degree());
        for (std::uint64_t &coefficient : u_coefficients)
            coefficient = uniform(engine);
        ring[level].forward(u_coefficients);
        result.residues[level].quadratic_u = std::move(u_coefficients);
        for (std::size_t i = 0; i < ring.degree(); i++)
            result.residues[level].quadratic_v[i] = modular_ntt::subtract(
                modular_ntt::multiply(result.residues[level].slots[i],
                                      result.residues[level].slots[i],
                                      ring[level].modulus()),
                modular_ntt::multiply(result.residues[level].quadratic_u[i],
                                      result.residues[level].slots[i],
                                      ring[level].modulus()),
                ring[level].modulus());
    }
    return result;
}

inline void validate(const RNSRing &ring, const RNSSecret &secret)
{
    if (secret.residues.size() != ring.levels())
        throw std::invalid_argument("RNS Binary-NTT secret has wrong levels");
    for (std::size_t level = 0; level < ring.levels(); level++)
        validate(ring[level], secret.residues[level]);
}

inline RNSCiphertext TrivialEncryptNTT(
    const RNSRing &ring,
    const std::vector<std::vector<std::uint64_t>> &messages)
{
    if (messages.size() != ring.levels())
        throw std::invalid_argument("RNS message has wrong levels");
    RNSCiphertext result;
    result.residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++)
        result.residues[level] = TrivialEncryptNTT(ring[level], messages[level]);
    return result;
}

template <class URBG>
RNSCiphertext EncryptNTT(
    const RNSRing &ring, const RNSSecret &secret,
    const std::vector<std::vector<std::uint64_t>> &messages,
    const double coefficient_noise_stddev, URBG &engine,
    const std::uint64_t error_multiplier = 1)
{
    validate(ring, secret);
    if (messages.size() != ring.levels())
        throw std::invalid_argument("RNS message has wrong levels");
    if (coefficient_noise_stddev < 0)
        throw std::invalid_argument("noise standard deviation must be nonnegative");
    // The RNS residues encode one RLWE ciphertext modulo the product of the
    // primes.  In particular its coefficient error must be shared across
    // residues; independently sampled errors CRT-reconstruct to a value of
    // size Q and invalidate modulus switching.
    std::normal_distribution<double> gaussian(0.0, coefficient_noise_stddev);
    std::vector<std::int64_t> shared_error(ring.degree());
    for (std::int64_t &coefficient : shared_error)
        coefficient = static_cast<std::int64_t>(std::llround(gaussian(engine)));
    RNSCiphertext result;
    result.residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++) {
        ring[level].requireSize(messages[level]);
        const std::uint64_t modulus = ring[level].modulus();
        std::uniform_int_distribution<std::uint64_t> uniform(0, modulus - 1);
        std::vector<std::uint64_t> error(ring.degree());
        for (std::size_t coefficient = 0; coefficient < ring.degree();
             coefficient++) {
            const std::int64_t value = shared_error[coefficient] *
                                       static_cast<std::int64_t>(error_multiplier);
            error[coefficient] = value >= 0
                ? static_cast<std::uint64_t>(value) % modulus
                : modular_ntt::negate(static_cast<std::uint64_t>(-value) %
                                      modulus, modulus);
        }
        ring[level].forward(error);
        Ciphertext &ciphertext = result.residues[level];
        ciphertext.a.resize(ring.degree());
        ciphertext.b.resize(ring.degree());
        for (std::size_t slot = 0; slot < ring.degree(); slot++) {
            ciphertext.a[slot] = uniform(engine);
            ciphertext.b[slot] = modular_ntt::add(
                messages[level][slot],
                modular_ntt::add(
                    modular_ntt::multiply(ciphertext.a[slot],
                                          secret.residues[level].slots[slot],
                                          modulus),
                    error[slot], modulus),
                modulus);
        }
    }
    return result;
}

template <class URBG>
RNSCiphertext EncryptLSBNTT(
    const RNSRing &ring, const RNSSecret &secret,
    const std::vector<std::vector<std::uint64_t>> &messages,
    const std::uint64_t plaintext_modulus,
    const double coefficient_noise_stddev, URBG &engine)
{
    if (plaintext_modulus < 2)
        throw std::invalid_argument("LSB plaintext modulus must be at least two");
    return EncryptNTT(ring, secret, messages, coefficient_noise_stddev,
                      engine, plaintext_modulus);
}

template <class URBG>
RNSBootstrappingKey BinaryNTTBootstrappingKeyGen(
    const RNSRing &ring, const RNSSecret &secret,
    const std::span<const std::uint8_t> lwe_secret_bits,
    const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, secret);
    RNSBootstrappingKey result;
    result.encrypted_lwe_bits.reserve(lwe_secret_bits.size());
    for (const std::uint8_t bit : lwe_secret_bits) {
        if (bit > 1)
            throw std::invalid_argument("Algorithm 3 requires binary LWE bits");
        std::vector<std::vector<std::uint64_t>> message(
            ring.levels(), std::vector<std::uint64_t>(ring.degree()));
        for (std::size_t level = 0; level < ring.levels(); level++) {
            message[level][0] = bit;
            ring[level].forward(message[level]);
        }
        result.encrypted_lwe_bits.push_back(
            EncryptNTT(ring, secret, message, coefficient_noise_stddev,
                       engine));
    }
    return result;
}

template <class URBG>
RNSPBCBootstrappingKey PBCBootstrappingKeyGen(
    const RNSRing &ring, const RNSSecret &secret, const PBCSchedule &schedule,
    const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, secret);
    RNSPBCBootstrappingKey result;
    result.entries.resize(schedule.bucket_count);
    result.dummies.reserve(schedule.bucket_count);
    for (std::size_t bucket = 0; bucket < schedule.bucket_indices.size();
         bucket++) {
        for (const std::uint32_t source : schedule.bucket_indices[bucket]) {
            std::vector<std::vector<std::uint64_t>> message(
                ring.levels(), std::vector<std::uint64_t>(ring.degree()));
            for (std::size_t level = 0; level < ring.levels(); level++) {
                message[level][0] = schedule.selected_source[bucket] ==
                                         static_cast<std::int32_t>(source)
                                     ? 1
                                     : 0;
                ring[level].forward(message[level]);
            }
            result.entries[bucket].push_back(
                EncryptNTT(ring, secret, message, coefficient_noise_stddev,
                           engine));
        }
        std::vector<std::vector<std::uint64_t>> dummy(
            ring.levels(), std::vector<std::uint64_t>(ring.degree()));
        for (std::size_t level = 0; level < ring.levels(); level++) {
            dummy[level][0] = schedule.selected_source[bucket] < 0 ? 1 : 0;
            ring[level].forward(dummy[level]);
        }
        result.dummies.push_back(
            EncryptNTT(ring, secret, dummy, coefficient_noise_stddev, engine));
    }
    return result;
}

template <class URBG>
RNSPBCBootstrappingKey PBCBootstrappingKeyGenLSB(
    const RNSRing &ring, const RNSSecret &secret, const PBCSchedule &schedule,
    const std::uint64_t plaintext_modulus,
    const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, secret);
    if (plaintext_modulus < 2)
        throw std::invalid_argument("LSB plaintext modulus must be at least two");
    RNSPBCBootstrappingKey result;
    result.entries.resize(schedule.bucket_count);
    result.dummies.reserve(schedule.bucket_count);
    for (std::size_t bucket = 0; bucket < schedule.bucket_indices.size(); bucket++) {
        for (const std::uint32_t source : schedule.bucket_indices[bucket]) {
            std::vector<std::vector<std::uint64_t>> message(
                ring.levels(), std::vector<std::uint64_t>(ring.degree()));
            for (std::size_t level = 0; level < ring.levels(); level++) {
                message[level][0] = schedule.selected_source[bucket] ==
                                            static_cast<std::int32_t>(source)
                                        ? 1
                                        : 0;
                ring[level].forward(message[level]);
            }
            result.entries[bucket].push_back(EncryptLSBNTT(
                ring, secret, message, plaintext_modulus,
                coefficient_noise_stddev, engine));
        }
        std::vector<std::vector<std::uint64_t>> dummy(
            ring.levels(), std::vector<std::uint64_t>(ring.degree()));
        for (std::size_t level = 0; level < ring.levels(); level++) {
            dummy[level][0] = schedule.selected_source[bucket] < 0 ? 1 : 0;
            ring[level].forward(dummy[level]);
        }
        result.dummies.push_back(EncryptLSBNTT(
            ring, secret, dummy, plaintext_modulus,
            coefficient_noise_stddev, engine));
    }
    return result;
}

#ifdef USE_BLAKE3
template <class URBG>
SeededRNSCiphertext EncryptSeededLSBNTT(
    const RNSRing &ring, const RNSSecret &secret,
    const std::vector<std::vector<std::uint64_t>> &messages,
    const std::uint64_t plaintext_modulus,
    const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, secret);
    if (messages.size() != ring.levels() || plaintext_modulus < 2)
        throw std::invalid_argument("invalid seeded LSB encryption input");
    SeededRNSCiphertext result;
    for (std::uint8_t &byte : result.a_seed)
        byte = static_cast<std::uint8_t>(engine());
    BLAKE3PRNG::BLAKE3PRNG<std::uint64_t> mask_engine(result.a_seed);
    std::normal_distribution<double> gaussian(0.0, coefficient_noise_stddev);
    std::vector<std::int64_t> shared_error(ring.degree());
    for (std::int64_t &coefficient : shared_error)
        coefficient = static_cast<std::int64_t>(std::llround(gaussian(engine))) *
                      static_cast<std::int64_t>(plaintext_modulus);
    result.b_residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++) {
        ring[level].requireSize(messages[level]);
        const std::uint64_t modulus = ring[level].modulus();
        std::uniform_int_distribution<std::uint64_t> uniform(0, modulus - 1);
        std::vector<std::uint64_t> error(ring.degree());
        for (std::size_t coefficient = 0; coefficient < ring.degree(); coefficient++) {
            const std::int64_t value = shared_error[coefficient];
            error[coefficient] = value >= 0
                ? static_cast<std::uint64_t>(value) % modulus
                : modular_ntt::negate(
                      static_cast<std::uint64_t>(-value) % modulus, modulus);
        }
        ring[level].forward(error);
        auto &b = result.b_residues[level];
        b.resize(ring.degree());
        for (std::size_t slot = 0; slot < ring.degree(); slot++) {
            const std::uint64_t a = uniform(mask_engine);
            b[slot] = modular_ntt::add(
                messages[level][slot],
                modular_ntt::add(
                    modular_ntt::multiply(a, secret.residues[level].slots[slot],
                                          modulus),
                    error[slot], modulus),
                modulus);
        }
    }
    return result;
}

template <class URBG>
SeededRNSPBCBootstrappingKey PBCBootstrappingKeyGenSeededLSB(
    const RNSRing &ring, const RNSSecret &secret, const PBCSchedule &schedule,
    const std::uint64_t plaintext_modulus,
    const double coefficient_noise_stddev, URBG &engine)
{
    SeededRNSPBCBootstrappingKey result;
    result.entries.resize(schedule.bucket_count);
    result.dummies.reserve(schedule.bucket_count);
    for (std::size_t bucket = 0; bucket < schedule.bucket_indices.size(); bucket++) {
        for (const std::uint32_t source : schedule.bucket_indices[bucket]) {
            std::vector<std::vector<std::uint64_t>> message(
                ring.levels(), std::vector<std::uint64_t>(ring.degree()));
            for (std::size_t level = 0; level < ring.levels(); level++) {
                message[level][0] = schedule.selected_source[bucket] ==
                                            static_cast<std::int32_t>(source)
                                        ? 1
                                        : 0;
                ring[level].forward(message[level]);
            }
            result.entries[bucket].push_back(EncryptSeededLSBNTT(
                ring, secret, message, plaintext_modulus,
                coefficient_noise_stddev, engine));
        }
        std::vector<std::vector<std::uint64_t>> dummy(
            ring.levels(), std::vector<std::uint64_t>(ring.degree()));
        for (std::size_t level = 0; level < ring.levels(); level++) {
            dummy[level][0] = schedule.selected_source[bucket] < 0 ? 1 : 0;
            ring[level].forward(dummy[level]);
        }
        result.dummies.push_back(EncryptSeededLSBNTT(
            ring, secret, dummy, plaintext_modulus,
            coefficient_noise_stddev, engine));
    }
    return result;
}
#endif

template <class URBG>
RNSKeySwitchKey BinaryNTTKeySwitchKeyGen(
    const RNSRing &ring, const RNSSecret &old_secret,
    const RNSSecret &new_secret, const std::uint32_t basebit,
    const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, old_secret);
    validate(ring, new_secret);
    RNSKeySwitchKey result;
    result.residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++)
        result.residues[level] = BinaryNTTKeySwitchKeyGen(
            ring[level], old_secret.residues[level],
            new_secret.residues[level], basebit, coefficient_noise_stddev,
            engine);
    return result;
}

template <class URBG>
RNSCoefficientKeySwitchKey CoefficientKeySwitchKeyGen(
    const RNSRing &ring, const RNSSecret &old_secret,
    const RNSSecret &new_secret, const std::uint32_t basebit,
    const double coefficient_noise_stddev, URBG &engine)
{
    validate(ring, old_secret);
    validate(ring, new_secret);
    if (basebit == 0 || basebit > 32)
        throw std::invalid_argument(
            "RNS coefficient key-switch base is unsupported");
    const std::uint32_t digit_count = (128 + basebit - 1) / basebit;
    const std::uint64_t base = std::uint64_t{1} << basebit;
    RNSCoefficientKeySwitchKey result;
    result.residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++) {
        result.residues[level].basebit = basebit;
        result.residues[level].entries.resize(ring.degree());
        for (auto &entry : result.residues[level].entries)
            entry.reserve(digit_count);
    }
    for (std::size_t coefficient = 0; coefficient < ring.degree(); coefficient++) {
        std::vector<std::vector<std::uint64_t>> messages(
            ring.levels(), std::vector<std::uint64_t>(ring.degree()));
        std::vector<std::uint64_t> scales(ring.levels(), 1);
        for (std::uint32_t digit = 0; digit < digit_count; digit++) {
            for (std::size_t level = 0; level < ring.levels(); level++) {
                const std::uint64_t modulus = ring[level].modulus();
                const auto monomial = ring[level].monomialNTT(coefficient);
                for (std::size_t slot = 0; slot < ring.degree(); slot++)
                    messages[level][slot] = modular_ntt::multiply(
                        modular_ntt::multiply(scales[level], monomial[slot],
                                              modulus),
                        old_secret.residues[level].slots[slot], modulus);
            }
            const RNSCiphertext encrypted = EncryptNTT(
                ring, new_secret, messages, coefficient_noise_stddev, engine);
            for (std::size_t level = 0; level < ring.levels(); level++)
                result.residues[level].entries[coefficient].push_back(
                    encrypted.residues[level]);
            for (std::size_t level = 0; level < ring.levels(); level++)
                scales[level] = modular_ntt::multiply(
                    scales[level], base, ring[level].modulus());
        }
    }
    return result;
}

inline RNSCiphertext RelinearizationFreeMultiply(
    const RNSRing &ring, const RNSSecret &secret, const RNSCiphertext &left,
    const RNSCiphertext &right)
{
    validate(ring, secret);
    if (left.residues.size() != ring.levels() ||
        right.residues.size() != ring.levels())
        throw std::invalid_argument("RNS ciphertext has wrong levels");
    RNSCiphertext result;
    result.residues.resize(ring.levels());
#pragma omp parallel for if (ring.levels() > 1)
    for (std::int64_t level = 0;
         level < static_cast<std::int64_t>(ring.levels()); level++)
        result.residues[static_cast<std::size_t>(level)] =
            RelinearizationFreeMultiply(
                ring[static_cast<std::size_t>(level)],
                secret.residues[static_cast<std::size_t>(level)],
                left.residues[static_cast<std::size_t>(level)],
                right.residues[static_cast<std::size_t>(level)]);
    return result;
}

inline RNSCiphertext BinaryNTTKeySwitch(const RNSRing &ring,
                                        const RNSKeySwitchKey &key,
                                        const RNSCiphertext &input)
{
    if (key.residues.size() != ring.levels() ||
        input.residues.size() != ring.levels())
        throw std::invalid_argument("RNS key-switch has wrong levels");
    RNSCiphertext result;
    result.residues.resize(ring.levels());
#pragma omp parallel for if (ring.levels() > 1)
    for (std::int64_t level = 0;
         level < static_cast<std::int64_t>(ring.levels()); level++)
        result.residues[static_cast<std::size_t>(level)] = BinaryNTTKeySwitch(
            ring[static_cast<std::size_t>(level)],
            key.residues[static_cast<std::size_t>(level)],
            input.residues[static_cast<std::size_t>(level)]);
    return result;
}

inline RNSCiphertext CoefficientKeySwitch(
    const RNSRing &ring, const RNSCoefficientKeySwitchKey &key,
    const RNSCiphertext &input)
{
    if (ring.levels() != 2 || key.residues.size() != ring.levels() ||
        input.residues.size() != ring.levels())
        throw std::invalid_argument("RNS coefficient key-switch has wrong levels");
    const std::uint32_t basebit = key.residues[0].basebit;
    if (basebit == 0 || basebit > 32 ||
        key.residues[1].basebit != basebit)
        throw std::invalid_argument("invalid RNS coefficient key-switch base");
    const std::uint32_t digit_count = (128 + basebit - 1) / basebit;
    const __int128 base = static_cast<__int128>(std::uint64_t{1} << basebit);
    const __int128 half_base = base / 2;
    std::array<std::vector<std::uint64_t>, 2> coefficients = {
        input.residues[0].a, input.residues[1].a};
    ring[0].inverse(coefficients[0]);
    ring[1].inverse(coefficients[1]);
    const modular_ntt::TwoPrimeCRT crt(ring.prime(0), ring.prime(1));
    RNSCiphertext result;
    result.residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++) {
        result.residues[level].a.assign(ring.degree(), 0);
        result.residues[level].b = input.residues[level].b;
    }
    for (std::size_t coefficient = 0; coefficient < ring.degree(); coefficient++) {
        __int128 value = crt.reconstructSigned(
            coefficients[0][coefficient], coefficients[1][coefficient]);
        for (std::uint32_t digit = 0; digit < digit_count; digit++) {
            __int128 decomposition = value % base;
            if (decomposition > half_base) decomposition -= base;
            if (decomposition < -half_base) decomposition += base;
            value = (value - decomposition) / base;
            if (decomposition == 0) continue;
            for (std::size_t level = 0; level < ring.levels(); level++) {
                if (key.residues[level].entries.size() != ring.degree() ||
                    key.residues[level].entries[coefficient].size() !=
                        digit_count)
                    throw std::invalid_argument(
                        "invalid RNS coefficient key-switch entry");
                const std::uint64_t modulus = ring[level].modulus();
                const std::uint64_t scalar = decomposition >= 0
                    ? static_cast<std::uint64_t>(decomposition) % modulus
                    : modular_ntt::negate(
                          static_cast<std::uint64_t>(-decomposition) % modulus,
                          modulus);
                const Ciphertext &entry =
                    key.residues[level].entries[coefficient][digit];
                for (std::size_t slot = 0; slot < ring.degree(); slot++) {
                    result.residues[level].a[slot] = modular_ntt::subtract(
                        result.residues[level].a[slot],
                        modular_ntt::multiply(scalar, entry.a[slot], modulus),
                        modulus);
                    result.residues[level].b[slot] = modular_ntt::subtract(
                        result.residues[level].b[slot],
                        modular_ntt::multiply(scalar, entry.b[slot], modulus),
                        modulus);
                }
            }
        }
        if (value != 0)
            throw std::invalid_argument(
                "RNS coefficient key-switch decomposition overflow");
    }
    return result;
}

inline RNSCiphertext MonomialFactor(const RNSRing &ring,
                                    const RNSCiphertext &bit_ct,
                                    const std::size_t exponent)
{
    if (bit_ct.residues.size() != ring.levels())
        throw std::invalid_argument("RNS ciphertext has wrong levels");
    RNSCiphertext result;
    result.residues.resize(ring.levels());
#pragma omp parallel for if (ring.levels() > 1)
    for (std::int64_t level = 0;
         level < static_cast<std::int64_t>(ring.levels()); level++)
        result.residues[static_cast<std::size_t>(level)] = MonomialFactor(
            ring[static_cast<std::size_t>(level)],
            bit_ct.residues[static_cast<std::size_t>(level)], exponent);
    return result;
}

inline RNSCiphertext AddCiphertexts(const RNSRing &ring,
                                    const RNSCiphertext &left,
                                    const RNSCiphertext &right)
{
    if (left.residues.size() != ring.levels() ||
        right.residues.size() != ring.levels())
        throw std::invalid_argument("RNS addition has wrong levels");
    RNSCiphertext result;
    result.residues.resize(ring.levels());
    // PBC factor construction performs many small additions.  Keep the two
    // RNS residues deterministic here; parallelism belongs to the product
    // tree, where the independent work is substantially larger.
    for (std::size_t level = 0; level < ring.levels(); level++)
        result.residues[level] = AddCiphertexts(
            ring[level], left.residues[level], right.residues[level]);
    return result;
}

inline RNSCiphertext MultiplyByMonomial(const RNSRing &ring,
                                        const RNSCiphertext &ciphertext,
                                        const std::size_t exponent)
{
    if (ciphertext.residues.size() != ring.levels())
        throw std::invalid_argument("RNS monomial multiply has wrong levels");
    RNSCiphertext result;
    result.residues.resize(ring.levels());
#pragma omp parallel for if (ring.levels() > 1)
    for (std::int64_t level = 0;
         level < static_cast<std::int64_t>(ring.levels()); level++)
        result.residues[static_cast<std::size_t>(level)] = MultiplyByMonomial(
            ring[static_cast<std::size_t>(level)],
            ciphertext.residues[static_cast<std::size_t>(level)], exponent);
    return result;
}

inline RNSCiphertext PBCPreparedFactor(
    const RNSRing &ring, const PBCSchedule &schedule,
    const RNSPBCBootstrappingKey &bootstrapping_key,
    const std::uint32_t bucket, const std::span<const std::uint32_t> lwe_a,
    const std::size_t exponent_scale)
{
    if (bucket >= schedule.bucket_count || lwe_a.size() != schedule.source_dimension ||
        bootstrapping_key.entries.size() != schedule.bucket_count ||
        bootstrapping_key.dummies.size() != schedule.bucket_count ||
        bootstrapping_key.entries[bucket].size() !=
            schedule.bucket_indices[bucket].size())
        throw std::invalid_argument("invalid RNS PBC prepared-factor inputs");
    RNSCiphertext result = bootstrapping_key.dummies[bucket];
#ifdef USE_HEXL
    std::vector<std::uint64_t> product(ring.degree());
#endif
    for (std::size_t slot = 0;
         slot < schedule.bucket_indices[bucket].size(); slot++) {
        const std::uint32_t source = schedule.bucket_indices[bucket][slot];
        const RNSCiphertext &entry = bootstrapping_key.entries[bucket][slot];
        const std::size_t exponent = lwe_a[source] * exponent_scale;
        for (std::size_t level = 0; level < ring.levels(); level++) {
            const std::uint64_t modulus = ring[level].modulus();
            const auto &monomial = ring[level].cachedMonomialNTT(exponent);
#ifdef USE_HEXL
            intel::hexl::EltwiseMultMod(
                product.data(), entry.residues[level].a.data(),
                monomial.data(), ring.degree(), modulus, 1);
            intel::hexl::EltwiseAddMod(
                result.residues[level].a.data(),
                result.residues[level].a.data(), product.data(),
                ring.degree(), modulus);
            intel::hexl::EltwiseMultMod(
                product.data(), entry.residues[level].b.data(),
                monomial.data(), ring.degree(), modulus, 1);
            intel::hexl::EltwiseAddMod(
                result.residues[level].b.data(),
                result.residues[level].b.data(), product.data(),
                ring.degree(), modulus);
#else
            for (std::size_t i = 0; i < ring.degree(); i++) {
                result.residues[level].a[i] = modular_ntt::add(
                    result.residues[level].a[i],
                    modular_ntt::multiply(entry.residues[level].a[i],
                                          monomial[i], modulus),
                    modulus);
                result.residues[level].b[i] = modular_ntt::add(
                    result.residues[level].b[i],
                    modular_ntt::multiply(entry.residues[level].b[i],
                                          monomial[i], modulus),
                    modulus);
            }
#endif
        }
    }
    return result;
}

#ifdef USE_BLAKE3
inline RNSCiphertext ExpandSeededRNSCiphertext(
    const RNSRing &ring, const SeededRNSCiphertext &seeded)
{
    if (seeded.b_residues.size() != ring.levels())
        throw std::invalid_argument("seeded RNS ciphertext has wrong levels");
    BLAKE3PRNG::BLAKE3PRNG<std::uint64_t> mask_engine(seeded.a_seed);
    RNSCiphertext result;
    result.residues.resize(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++) {
        ring[level].requireSize(seeded.b_residues[level]);
        std::uniform_int_distribution<std::uint64_t> uniform(
            0, ring[level].modulus() - 1);
        result.residues[level].a.resize(ring.degree());
        for (std::uint64_t &value : result.residues[level].a)
            value = uniform(mask_engine);
        result.residues[level].b = seeded.b_residues[level];
    }
    return result;
}

inline RNSCiphertext PBCPreparedFactor(
    const RNSRing &ring, const PBCSchedule &schedule,
    const SeededRNSPBCBootstrappingKey &bootstrapping_key,
    const std::uint32_t bucket, const std::span<const std::uint32_t> lwe_a,
    const std::size_t exponent_scale)
{
    if (bucket >= schedule.bucket_count ||
        lwe_a.size() != schedule.source_dimension ||
        bootstrapping_key.entries.size() != schedule.bucket_count ||
        bootstrapping_key.dummies.size() != schedule.bucket_count ||
        bootstrapping_key.entries[bucket].size() !=
            schedule.bucket_indices[bucket].size())
        throw std::invalid_argument("invalid seeded RNS PBC factor inputs");
    RNSCiphertext result =
        ExpandSeededRNSCiphertext(ring, bootstrapping_key.dummies[bucket]);
    std::vector<std::uint64_t> a(ring.degree());
#ifdef USE_HEXL
    std::vector<std::uint64_t> product(ring.degree());
#endif
    for (std::size_t slot = 0;
         slot < schedule.bucket_indices[bucket].size(); slot++) {
        const std::uint32_t source = schedule.bucket_indices[bucket][slot];
        const SeededRNSCiphertext &entry =
            bootstrapping_key.entries[bucket][slot];
        if (entry.b_residues.size() != ring.levels())
            throw std::invalid_argument("seeded RNS PBC entry has wrong levels");
        BLAKE3PRNG::BLAKE3PRNG<std::uint64_t> mask_engine(entry.a_seed);
        const std::size_t exponent = lwe_a[source] * exponent_scale;
        for (std::size_t level = 0; level < ring.levels(); level++) {
            const std::uint64_t modulus = ring[level].modulus();
            std::uniform_int_distribution<std::uint64_t> uniform(0, modulus - 1);
            for (std::uint64_t &value : a) value = uniform(mask_engine);
            const auto &monomial = ring[level].cachedMonomialNTT(exponent);
#ifdef USE_HEXL
            intel::hexl::EltwiseMultMod(
                product.data(), a.data(), monomial.data(), ring.degree(),
                modulus, 1);
            intel::hexl::EltwiseAddMod(
                result.residues[level].a.data(),
                result.residues[level].a.data(), product.data(),
                ring.degree(), modulus);
            intel::hexl::EltwiseMultMod(
                product.data(), entry.b_residues[level].data(),
                monomial.data(), ring.degree(), modulus, 1);
            intel::hexl::EltwiseAddMod(
                result.residues[level].b.data(),
                result.residues[level].b.data(), product.data(),
                ring.degree(), modulus);
#else
            for (std::size_t i = 0; i < ring.degree(); i++) {
                result.residues[level].a[i] = modular_ntt::add(
                    result.residues[level].a[i],
                    modular_ntt::multiply(a[i], monomial[i], modulus),
                    modulus);
                result.residues[level].b[i] = modular_ntt::add(
                    result.residues[level].b[i],
                    modular_ntt::multiply(entry.b_residues[level][i],
                                          monomial[i], modulus),
                    modulus);
            }
#endif
        }
    }
    return result;
}
#endif

inline RNSCiphertext ProductTree(const RNSRing &ring, const RNSSecret &secret,
                                 std::vector<RNSCiphertext> inputs,
                                 ProductTreeStats *stats = nullptr)
{
    if (inputs.empty())
        throw std::invalid_argument("RNS product tree requires an input");
    ProductTreeStats local_stats;
    while (inputs.size() > 1) {
        const std::size_t pairs = inputs.size() / 2;
        std::vector<RNSCiphertext> next((inputs.size() + 1) / 2);
#pragma omp parallel for if (pairs > 1)
        for (std::int64_t pair = 0; pair < static_cast<std::int64_t>(pairs);
             pair++)
            next[static_cast<std::size_t>(pair)] =
                RelinearizationFreeMultiply(
                    ring, secret, inputs[2 * static_cast<std::size_t>(pair)],
                    inputs[2 * static_cast<std::size_t>(pair) + 1]);
        if ((inputs.size() & 1U) != 0)
            next.back() = std::move(inputs.back());
        local_stats.pointwise_multiplication_layers++;
        local_stats.pointwise_ciphertext_multiplications +=
            static_cast<std::uint32_t>(pairs);
        inputs = std::move(next);
    }
    if (stats != nullptr) *stats = local_stats;
    return std::move(inputs.front());
}

inline RNSCiphertext LowDepthBlindRotate(
    const RNSRing &ring, const RNSSecret &secret, const RNSCiphertext &initial,
    const RNSBootstrappingKey &bootstrapping_key,
    const std::span<const std::size_t> exponents,
    ProductTreeStats *stats = nullptr)
{
    if (bootstrapping_key.encrypted_lwe_bits.size() != exponents.size())
        throw std::invalid_argument(
            "low-depth RNS blind rotation has mismatched factors and exponents");
    std::vector<RNSCiphertext> factors;
    factors.reserve(exponents.size() + 1);
    factors.push_back(initial);
    for (std::size_t i = 0; i < exponents.size(); i++)
        factors.push_back(MonomialFactor(
            ring, bootstrapping_key.encrypted_lwe_bits[i], exponents[i]));
    return ProductTree(ring, secret, std::move(factors), stats);
}

inline std::vector<std::vector<std::uint64_t>> DecryptNTT(
    const RNSRing &ring, const RNSSecret &secret,
    const RNSCiphertext &ciphertext)
{
    validate(ring, secret);
    if (ciphertext.residues.size() != ring.levels())
        throw std::invalid_argument("RNS ciphertext has wrong levels");
    std::vector<std::vector<std::uint64_t>> result(ring.levels());
    for (std::size_t level = 0; level < ring.levels(); level++)
        result[level] = DecryptNTT(ring[level], secret.residues[level],
                                   ciphertext.residues[level]);
    return result;
}

inline RNSCiphertext Algorithm3BlindRotate(
    const RNSRing &ring, const RNSSecret &secret,
    const RNSBootstrappingKey &bootstrapping_key, const LWEPhase &lwe,
    const std::vector<std::vector<std::uint64_t>> &test_vector,
    ProductTreeStats *stats = nullptr)
{
    if (lwe.modulus == 0 || (2 * ring.degree()) % lwe.modulus != 0)
        throw std::invalid_argument(
            "Algorithm 3 requires an LWE modulus dividing 2N");
    if (lwe.a.size() != bootstrapping_key.encrypted_lwe_bits.size() ||
        test_vector.size() != ring.levels())
        throw std::invalid_argument("Algorithm 3 RNS input has wrong dimensions");
    const std::size_t scale = 2 * ring.degree() / lwe.modulus;
    const std::size_t initial_exponent =
        (lwe.modulus - (lwe.b % lwe.modulus)) % lwe.modulus;
    std::vector<std::vector<std::uint64_t>> initial_message = test_vector;
    for (std::size_t level = 0; level < ring.levels(); level++) {
        ring[level].requireSize(initial_message[level]);
        const auto monomial =
            ring[level].monomialNTT(initial_exponent * scale);
        for (std::size_t i = 0; i < ring.degree(); i++)
            initial_message[level][i] = modular_ntt::multiply(
                initial_message[level][i], monomial[i], ring[level].modulus());
    }
    const RNSCiphertext initial = TrivialEncryptNTT(ring, initial_message);
    std::vector<std::size_t> exponents(lwe.a.size());
    for (std::size_t i = 0; i < lwe.a.size(); i++)
        exponents[i] = (lwe.a[i] % lwe.modulus) * scale;
    return LowDepthBlindRotate(ring, secret, initial, bootstrapping_key,
                               exponents, stats);
}

inline std::uint64_t RoundRNSModulusSwitchCoefficient(
    const __int128 value, const unsigned __int128 source_modulus,
    const std::uint64_t target_modulus)
{
    using boost::multiprecision::cpp_int;
    const bool negative = value < 0;
    const unsigned __int128 magnitude = negative
        ? static_cast<unsigned __int128>(-(value + 1)) + 1
        : static_cast<unsigned __int128>(value);
    const cpp_int numerator = cpp_int(magnitude) * target_modulus +
                              cpp_int(source_modulus) / 2;
    const std::uint64_t rounded = static_cast<std::uint64_t>(
        numerator / cpp_int(source_modulus));
    const std::uint64_t reduced = rounded % target_modulus;
    return negative ? modular_ntt::negate(reduced, target_modulus) : reduced;
}

// Exact CRT modulus switch used at an Algorithm-3 window boundary.
// The high RNS secret must already have been switched to a coefficient-small
// quadratic-hint secret corresponding to `target_secret`.
inline Ciphertext ModulusSwitch(const RNSRing &source, const Ring &target,
                                const RNSCiphertext &input)
{
    using boost::multiprecision::cpp_int;
    if (source.levels() < 2 || source.degree() != target.degree() ||
        input.residues.size() != source.levels())
        throw std::invalid_argument(
            "RNS modulus switch requires matching equal-degree levels");
    std::vector<Ciphertext> coefficients = input.residues;
    cpp_int product = 1;
    for (std::size_t level = 0; level < source.levels(); level++) {
        source[level].inverse(coefficients[level].a);
        source[level].inverse(coefficients[level].b);
        product *= source[level].modulus();
    }
    const auto reconstruct = [&](const bool b_component,
                                 const std::size_t coefficient) {
        cpp_int value = b_component ? coefficients[0].b[coefficient]
                                    : coefficients[0].a[coefficient];
        cpp_int accumulated_modulus = source[0].modulus();
        for (std::size_t level = 1; level < source.levels(); level++) {
            const std::uint64_t prime = source[level].modulus();
            const std::uint64_t current = static_cast<std::uint64_t>(
                (value % prime).convert_to<std::uint64_t>());
            const std::uint64_t wanted = b_component
                ? coefficients[level].b[coefficient]
                : coefficients[level].a[coefficient];
            const std::uint64_t difference = modular_ntt::subtract(
                wanted, current, prime);
            const std::uint64_t accumulated_mod_prime =
                static_cast<std::uint64_t>(
                    (accumulated_modulus % prime).convert_to<std::uint64_t>());
            const std::uint64_t quotient = modular_ntt::multiply(
                difference, modular_ntt::invert(accumulated_mod_prime, prime),
                prime);
            value += accumulated_modulus * quotient;
            accumulated_modulus *= prime;
        }
        if (value > product / 2) value -= product;
        return value;
    };
    const auto scale = [&](const cpp_int &value) {
        const bool negative = value < 0;
        const cpp_int magnitude = negative ? -value : value;
        const cpp_int rounded =
            (magnitude * target.modulus() + product / 2) / product;
        const std::uint64_t residue = static_cast<std::uint64_t>(
            (rounded % target.modulus()).convert_to<std::uint64_t>());
        return negative ? modular_ntt::negate(residue, target.modulus())
                        : residue;
    };
    Ciphertext result;
    result.a.resize(target.degree());
    result.b.resize(target.degree());
    for (std::size_t i = 0; i < target.degree(); i++) {
        result.a[i] = scale(reconstruct(false, i));
        result.b[i] = scale(reconstruct(true, i));
    }
    target.forward(result.a);
    target.forward(result.b);
    return result;
}

inline bool IsNestedRNSPrefix(const RNSRing &source, const RNSRing &target)
{
    if (target.levels() >= source.levels() ||
        target.degree() != source.degree())
        return false;
    for (std::size_t level = 0; level < target.levels(); level++)
        if (source.prime(level).value != target.prime(level).value)
            return false;
    return true;
}

inline bool IsSameRNSBasis(const RNSRing &left, const RNSRing &right)
{
    if (left.degree() != right.degree() || left.levels() != right.levels())
        return false;
    for (std::size_t level = 0; level < left.levels(); level++)
        if (left.prime(level).value != right.prime(level).value)
            return false;
    return true;
}

inline RNSCiphertext ModulusSwitchNestedLSB(
    const RNSRing &source, const RNSRing &target,
    const RNSCiphertext &input, const std::uint64_t plaintext_modulus)
{
    if (!IsNestedRNSPrefix(source, target) || plaintext_modulus < 2 ||
        input.residues.size() != source.levels())
        throw std::invalid_argument("invalid nested RNS LSB modulus switch");
    std::vector<Ciphertext> coefficients = input.residues;
    for (std::size_t level = 0; level < source.levels(); level++) {
        source[level].inverse(coefficients[level].a);
        source[level].inverse(coefficients[level].b);
    }
    // OpenFHE-style BGV ModReduce: discard one tower at a time. The correction
    // lives in the discarded prime, is centered when mapped to surviving
    // towers, and makes division by q_l exact. No CRT interpolation is needed.
    for (std::size_t active = source.levels(); active > target.levels();
         active--) {
        const std::size_t dropped = active - 1;
        const std::uint64_t ql = source[dropped].modulus();
        const std::uint64_t t_inverse = modular_ntt::invert(
            plaintext_modulus % ql, ql);
        const std::uint64_t neg_t_inverse =
            modular_ntt::negate(t_inverse, ql);
        std::vector<std::uint64_t> ql_inverse(dropped);
        for (std::size_t level = 0; level < dropped; level++)
            ql_inverse[level] = modular_ntt::invert(
                ql % source[level].modulus(), source[level].modulus());
#pragma omp parallel for if (source.degree() > 1)
        for (std::int64_t signed_coefficient = 0;
             signed_coefficient < static_cast<std::int64_t>(source.degree());
             signed_coefficient++) {
            const std::size_t coefficient =
                static_cast<std::size_t>(signed_coefficient);
            for (std::size_t component = 0; component < 2; component++) {
                const std::uint64_t dropped_value = component == 0
                    ? coefficients[dropped].a[coefficient]
                    : coefficients[dropped].b[coefficient];
                const std::uint64_t delta_residue = modular_ntt::multiply(
                    dropped_value, neg_t_inverse, ql);
                const std::int64_t delta = delta_residue <= ql / 2
                    ? static_cast<std::int64_t>(delta_residue)
                    : -static_cast<std::int64_t>(ql - delta_residue);
                for (std::size_t level = 0; level < dropped; level++) {
                    const std::uint64_t prime = source[level].modulus();
                    const std::uint64_t delta_mod = delta >= 0
                        ? static_cast<std::uint64_t>(delta) % prime
                        : modular_ntt::negate(
                              static_cast<std::uint64_t>(-delta) % prime,
                              prime);
                    const std::uint64_t correction = modular_ntt::multiply(
                        plaintext_modulus % prime, delta_mod, prime);
                    std::uint64_t &value = component == 0
                        ? coefficients[level].a[coefficient]
                        : coefficients[level].b[coefficient];
                    value = modular_ntt::multiply(
                        modular_ntt::add(value, correction, prime),
                        ql_inverse[level], prime);
                }
            }
        }
    }
    RNSCiphertext result;
    result.residues.assign(coefficients.begin(),
                           coefficients.begin() + target.levels());
    for (std::size_t level = 0; level < target.levels(); level++) {
        target[level].forward(result.residues[level].a);
        target[level].forward(result.residues[level].b);
    }
    return result;
}

inline RNSCiphertext ModulusSwitch(const RNSRing &source,
                                   const RNSRing &target,
                                   const RNSCiphertext &input,
                                   const std::uint64_t plaintext_modulus = 0)
{
    using boost::multiprecision::cpp_int;
    if (source.levels() == 0 || target.levels() == 0 ||
        source.degree() != target.degree() ||
        input.residues.size() != source.levels() || plaintext_modulus == 1)
        throw std::invalid_argument("invalid RNS-to-RNS modulus switch");
    if (plaintext_modulus != 0 && IsNestedRNSPrefix(source, target))
        return ModulusSwitchNestedLSB(
            source, target, input, plaintext_modulus);
    std::vector<Ciphertext> coefficients = input.residues;
    cpp_int source_modulus = 1;
    cpp_int target_modulus = 1;
    for (std::size_t level = 0; level < source.levels(); level++) {
        source[level].inverse(coefficients[level].a);
        source[level].inverse(coefficients[level].b);
        source_modulus *= source[level].modulus();
    }
    for (std::size_t level = 0; level < target.levels(); level++)
        target_modulus *= target[level].modulus();
    cpp_int dropped_modulus = 0;
    cpp_int plaintext_inverse = 0;
    if (plaintext_modulus != 0) {
        if (source_modulus % target_modulus != 0)
            throw std::invalid_argument(
                "LSB modulus switch requires a nested RNS modulus chain");
        dropped_modulus = source_modulus / target_modulus;
        // Extended Euclid over cpp_int for t^{-1} modulo the dropped basis.
        cpp_int old_r = plaintext_modulus;
        cpp_int r = dropped_modulus;
        cpp_int old_s = 1;
        cpp_int s = 0;
        while (r != 0) {
            const cpp_int quotient = old_r / r;
            const cpp_int next_r = old_r - quotient * r;
            old_r = r;
            r = next_r;
            const cpp_int next_s = old_s - quotient * s;
            old_s = s;
            s = next_s;
        }
        if (old_r != 1)
            throw std::invalid_argument(
                "plaintext modulus is not invertible in dropped RNS basis");
        plaintext_inverse = (old_s % dropped_modulus + dropped_modulus) %
                            dropped_modulus;
    }
    const auto reconstruct = [&](const bool b_component,
                                 const std::size_t coefficient) {
        cpp_int value = b_component ? coefficients[0].b[coefficient]
                                    : coefficients[0].a[coefficient];
        cpp_int accumulated = source[0].modulus();
        for (std::size_t level = 1; level < source.levels(); level++) {
            const std::uint64_t prime = source[level].modulus();
            const std::uint64_t current =
                (value % prime).convert_to<std::uint64_t>();
            const std::uint64_t wanted = b_component
                ? coefficients[level].b[coefficient]
                : coefficients[level].a[coefficient];
            const std::uint64_t difference =
                modular_ntt::subtract(wanted, current, prime);
            const std::uint64_t accumulated_mod =
                (accumulated % prime).convert_to<std::uint64_t>();
            const std::uint64_t quotient = modular_ntt::multiply(
                difference, modular_ntt::invert(accumulated_mod, prime), prime);
            value += accumulated * quotient;
            accumulated *= prime;
        }
        if (plaintext_modulus == 0 && value > source_modulus / 2)
            value -= source_modulus;
        return value;
    };
    RNSCiphertext result;
    result.residues.resize(target.levels());
    for (auto &residue : result.residues) {
        residue.a.resize(target.degree());
        residue.b.resize(target.degree());
    }
    for (std::size_t i = 0; i < target.degree(); i++)
        for (std::size_t component = 0; component < 2; component++) {
            const cpp_int value = reconstruct(component == 1, i);
            cpp_int signed_rounded;
            if (plaintext_modulus == 0) {
                const bool negative = value < 0;
                const cpp_int magnitude = negative ? -value : value;
                const cpp_int rounded =
                    (magnitude * target_modulus + source_modulus / 2) /
                    source_modulus;
                signed_rounded = negative ? -rounded : rounded;
            }
            else {
                // Appendix A.1's BGV/LSB switch. Select k modulo the dropped
                // basis so value-t*k is exactly divisible by Q/Q'.
                const cpp_int k =
                    (value % dropped_modulus) * plaintext_inverse %
                    dropped_modulus;
                signed_rounded =
                    (value - cpp_int(plaintext_modulus) * k) /
                    dropped_modulus;
            }
            for (std::size_t level = 0; level < target.levels(); level++) {
                const std::uint64_t prime = target[level].modulus();
                cpp_int reduced = signed_rounded % prime;
                if (reduced < 0) reduced += prime;
                const std::uint64_t residue =
                    reduced.convert_to<std::uint64_t>();
                if (component == 0)
                    result.residues[level].a[i] = residue;
                else
                    result.residues[level].b[i] = residue;
            }
        }
    for (std::size_t level = 0; level < target.levels(); level++) {
        target[level].forward(result.residues[level].a);
        target[level].forward(result.residues[level].b);
    }
    return result;
}

inline Ciphertext Algorithm3ModulusBoundary(
    const RNSRing &source_ring, const RNSSecret &binary_ntt_secret,
    const RNSKeySwitchKey &switch_key, const RNSSecret &small_secret_source,
    const Ring &target_ring, const Secret &small_secret_target,
    const RNSCiphertext &input)
{
    validate(source_ring, binary_ntt_secret);
    validate(source_ring, small_secret_source);
    validate(target_ring, small_secret_target);
    const RNSCiphertext switched =
        BinaryNTTKeySwitch(source_ring, switch_key, input);
    return ModulusSwitch(source_ring, target_ring, switched);
}

inline Ciphertext Algorithm3ModulusBoundary(
    const RNSRing &source_ring, const RNSSecret &binary_ntt_secret,
    const RNSCoefficientKeySwitchKey &switch_key,
    const RNSSecret &small_secret_source, const Ring &target_ring,
    const Secret &small_secret_target, const RNSCiphertext &input)
{
    validate(source_ring, binary_ntt_secret);
    validate(source_ring, small_secret_source);
    validate(target_ring, small_secret_target);
    const RNSCiphertext switched =
        CoefficientKeySwitch(source_ring, switch_key, input);
    return ModulusSwitch(source_ring, target_ring, switched);
}

template <class BoundaryKey>
inline Ciphertext Algorithm3TwoStageProduct(
    const RNSRing &high_ring, const RNSSecret &binary_ntt_secret,
    std::vector<RNSCiphertext> factors, const std::uint32_t high_window,
    const std::uint32_t low_window, const BoundaryKey &boundary_key,
    const RNSSecret &high_small_secret, const Ring &low_ring,
    const Secret &low_small_secret, TwoStageStats *stats = nullptr,
    const KeySwitchKey *low_restore_key = nullptr,
    const Secret *low_product_secret = nullptr)
{
    if (factors.empty() || high_window < 2 || low_window < 2 ||
        (high_window & (high_window - 1)) != 0 ||
        (low_window & (low_window - 1)) != 0)
        throw std::invalid_argument("Algorithm 3 windows must be powers of two");
    if ((low_restore_key == nullptr) != (low_product_secret == nullptr))
        throw std::invalid_argument(
            "Algorithm 3 low restore key and secret must be supplied together");
    if (low_product_secret != nullptr) validate(low_ring, *low_product_secret);
    TwoStageStats local_stats;
    std::vector<Ciphertext> low_inputs;
    low_inputs.reserve((factors.size() + high_window - 1) / high_window);
    for (std::size_t first = 0; first < factors.size(); first += high_window) {
        const std::size_t last = std::min<std::size_t>(
            first + high_window, factors.size());
        std::vector<RNSCiphertext> group;
        group.reserve(last - first);
        for (std::size_t i = first; i < last; i++)
            group.push_back(std::move(factors[i]));
        ProductTreeStats group_stats;
        const RNSCiphertext high_product = ProductTree(
            high_ring, binary_ntt_secret, std::move(group), &group_stats);
        local_stats.high_product.pointwise_multiplication_layers = std::max(
            local_stats.high_product.pointwise_multiplication_layers,
            group_stats.pointwise_multiplication_layers);
        local_stats.high_product.pointwise_ciphertext_multiplications +=
            group_stats.pointwise_ciphertext_multiplications;
        Ciphertext boundary = Algorithm3ModulusBoundary(
            high_ring, binary_ntt_secret, boundary_key, high_small_secret,
            low_ring, low_small_secret, high_product);
        if (low_restore_key != nullptr)
            boundary = BinaryNTTKeySwitch(low_ring, *low_restore_key, boundary);
        low_inputs.push_back(std::move(boundary));
        local_stats.modulus_boundaries++;
    }
    if (low_inputs.size() > low_window)
        throw std::invalid_argument(
            "Algorithm 3 low window cannot combine all boundary outputs");
    const Secret &product_secret = low_product_secret == nullptr
        ? low_small_secret
        : *low_product_secret;
    const Ciphertext result = ProductTree(low_ring, product_secret,
                                          std::move(low_inputs),
                                          &local_stats.low_product);
    if (stats != nullptr) *stats = local_stats;
    return result;
}

// Concrete two-stage realization of Algorithm 3's Section 7.2 shape.  The
// high stage multiplies groups of `high_window` prepared factors in RNS under
// the Binary-NTT secret. Each group then crosses the key-switch/modulus-switch
// boundary and the low stage multiplies the resulting ciphertexts under the
// coefficient-small quadratic-hint secret. For the paper schedule use 8,4.
template <class BoundaryKey>
inline Ciphertext Algorithm3TwoStageBlindRotate(
    const RNSRing &high_ring, const RNSSecret &binary_ntt_secret,
    const RNSBootstrappingKey &bootstrapping_key, const LWEPhase &lwe,
    const std::vector<std::vector<std::uint64_t>> &test_vector,
    const std::uint32_t high_window, const std::uint32_t low_window,
    const BoundaryKey &boundary_key, const RNSSecret &high_small_secret,
    const Ring &low_ring, const Secret &low_small_secret,
    TwoStageStats *stats = nullptr,
    const KeySwitchKey *low_restore_key = nullptr,
    const Secret *low_product_secret = nullptr)
{
    if (lwe.modulus == 0 || (2 * high_ring.degree()) % lwe.modulus != 0 ||
        lwe.a.size() != bootstrapping_key.encrypted_lwe_bits.size() ||
        test_vector.size() != high_ring.levels())
        throw std::invalid_argument("Algorithm 3 two-stage input is invalid");

    const std::size_t scale = 2 * high_ring.degree() / lwe.modulus;
    const std::size_t initial_exponent =
        (lwe.modulus - (lwe.b % lwe.modulus)) % lwe.modulus;
    std::vector<std::vector<std::uint64_t>> initial_message = test_vector;
    for (std::size_t level = 0; level < high_ring.levels(); level++) {
        high_ring[level].requireSize(initial_message[level]);
        const auto monomial =
            high_ring[level].monomialNTT(initial_exponent * scale);
        for (std::size_t i = 0; i < high_ring.degree(); i++)
            initial_message[level][i] = modular_ntt::multiply(
                initial_message[level][i], monomial[i],
                high_ring[level].modulus());
    }
    std::vector<RNSCiphertext> factors;
    factors.reserve(lwe.a.size() + 1);
    factors.push_back(TrivialEncryptNTT(high_ring, initial_message));
    for (std::size_t i = 0; i < lwe.a.size(); i++)
        factors.push_back(MonomialFactor(
            high_ring, bootstrapping_key.encrypted_lwe_bits[i],
            (lwe.a[i] % lwe.modulus) * scale));

    return Algorithm3TwoStageProduct(
        high_ring, binary_ntt_secret, std::move(factors), high_window,
        low_window, boundary_key, high_small_secret, low_ring,
        low_small_secret, stats, low_restore_key, low_product_secret);
}

template <class BoundaryKey>
inline Ciphertext Algorithm3PBCTwoStageBlindRotate(
    const RNSRing &high_ring, const RNSSecret &binary_ntt_secret,
    const PBCSchedule &schedule,
    const RNSPBCBootstrappingKey &bootstrapping_key, const LWEPhase &lwe,
    const std::vector<std::vector<std::uint64_t>> &test_vector,
    const std::uint32_t high_window, const std::uint32_t low_window,
    const BoundaryKey &boundary_key, const RNSSecret &high_small_secret,
    const Ring &low_ring, const Secret &low_small_secret,
    TwoStageStats *stats = nullptr,
    const KeySwitchKey *low_restore_key = nullptr,
    const Secret *low_product_secret = nullptr)
{
    if (lwe.modulus == 0 || (2 * high_ring.degree()) % lwe.modulus != 0 ||
        lwe.a.size() != schedule.source_dimension ||
        test_vector.size() != high_ring.levels())
        throw std::invalid_argument("invalid Algorithm 3 PBC two-stage input");
    const std::size_t scale = 2 * high_ring.degree() / lwe.modulus;
    const std::size_t initial_exponent =
        (lwe.modulus - (lwe.b % lwe.modulus)) % lwe.modulus;
    std::vector<std::vector<std::uint64_t>> initial_message = test_vector;
    for (std::size_t level = 0; level < high_ring.levels(); level++) {
        const auto monomial =
            high_ring[level].monomialNTT(initial_exponent * scale);
        for (std::size_t i = 0; i < high_ring.degree(); i++)
            initial_message[level][i] = modular_ntt::multiply(
                initial_message[level][i], monomial[i],
                high_ring[level].modulus());
    }
    const RNSCiphertext initial =
        TrivialEncryptNTT(high_ring, initial_message);
    std::vector<RNSCiphertext> factors;
    factors.reserve(schedule.bucket_count);
    for (std::uint32_t bucket = 0; bucket < schedule.bucket_count; bucket++) {
        RNSCiphertext factor = PBCPreparedFactor(
            high_ring, schedule, bootstrapping_key, bucket, lwe.a, scale);
        if (bucket == 0)
            for (std::size_t level = 0; level < high_ring.levels(); level++)
                for (std::size_t i = 0; i < high_ring.degree(); i++) {
                    const std::uint64_t modulus = high_ring[level].modulus();
                    factor.residues[level].a[i] = modular_ntt::multiply(
                        factor.residues[level].a[i],
                        initial.residues[level].b[i], modulus);
                    factor.residues[level].b[i] = modular_ntt::multiply(
                        factor.residues[level].b[i],
                        initial.residues[level].b[i], modulus);
                }
        factors.push_back(std::move(factor));
    }
    return Algorithm3TwoStageProduct(
        high_ring, binary_ntt_secret, std::move(factors), high_window,
        low_window, boundary_key, high_small_secret, low_ring,
        low_small_secret, stats, low_restore_key, low_product_secret);
}

// D2/QH-SS realization of the Algorithm-3 boundary.  Unlike the D1
// Binary-NTT route it does not need a key switch before modulus switching:
// the coefficient-small QH secret is shared across RNS residues and remains
// short after modulus reduction.
inline Ciphertext Algorithm3QHTwoStageProduct(
    const RNSRing &high_ring, const RNSSecret &high_qh_secret,
    std::vector<RNSCiphertext> factors, const std::uint32_t high_window,
    const std::uint32_t low_window, const Ring &low_ring,
    const Secret &low_qh_secret, TwoStageStats *stats = nullptr)
{
    if (factors.empty() || high_window < 2 || low_window < 2 ||
        (high_window & (high_window - 1)) != 0 ||
        (low_window & (low_window - 1)) != 0)
        throw std::invalid_argument("Algorithm 3 windows must be powers of two");
    validate(high_ring, high_qh_secret);
    validate(low_ring, low_qh_secret);
    TwoStageStats local_stats;
    std::vector<Ciphertext> low_inputs;
    low_inputs.reserve((factors.size() + high_window - 1) / high_window);
    for (std::size_t first = 0; first < factors.size(); first += high_window) {
        const std::size_t last = std::min<std::size_t>(
            first + high_window, factors.size());
        std::vector<RNSCiphertext> group;
        group.reserve(last - first);
        for (std::size_t i = first; i < last; i++)
            group.push_back(std::move(factors[i]));
        ProductTreeStats group_stats;
        const RNSCiphertext high_product = ProductTree(
            high_ring, high_qh_secret, std::move(group), &group_stats);
        local_stats.high_product.pointwise_multiplication_layers = std::max(
            local_stats.high_product.pointwise_multiplication_layers,
            group_stats.pointwise_multiplication_layers);
        local_stats.high_product.pointwise_ciphertext_multiplications +=
            group_stats.pointwise_ciphertext_multiplications;
        low_inputs.push_back(ModulusSwitch(high_ring, low_ring, high_product));
        local_stats.modulus_boundaries++;
    }
    if (low_inputs.size() > low_window)
        throw std::invalid_argument(
            "Algorithm 3 low window cannot combine all boundary outputs");
    const Ciphertext result = ProductTree(low_ring, low_qh_secret,
                                          std::move(low_inputs),
                                          &local_stats.low_product);
    if (stats != nullptr) *stats = local_stats;
    return result;
}

inline RNSCiphertext Algorithm3QHTwoStageProductRNS(
    const RNSRing &high_ring, const RNSSecret &high_qh_secret,
    std::vector<RNSCiphertext> factors, const std::uint32_t high_window,
    const std::uint32_t low_window, const RNSRing &low_ring,
    const RNSSecret &low_qh_secret, TwoStageStats *stats = nullptr)
{
    if (factors.empty() || high_window < 2 || low_window < 2 ||
        (high_window & (high_window - 1)) != 0 ||
        (low_window & (low_window - 1)) != 0)
        throw std::invalid_argument("Algorithm 3 windows must be powers of two");
    validate(high_ring, high_qh_secret);
    validate(low_ring, low_qh_secret);
    TwoStageStats local_stats;
    std::vector<RNSCiphertext> low_inputs;
    for (std::size_t first = 0; first < factors.size(); first += high_window) {
        const std::size_t last = std::min<std::size_t>(
            first + high_window, factors.size());
        std::vector<RNSCiphertext> group;
        for (std::size_t i = first; i < last; i++)
            group.push_back(std::move(factors[i]));
        ProductTreeStats group_stats;
        const RNSCiphertext high_product = ProductTree(
            high_ring, high_qh_secret, std::move(group), &group_stats);
        local_stats.high_product.pointwise_multiplication_layers = std::max(
            local_stats.high_product.pointwise_multiplication_layers,
            group_stats.pointwise_multiplication_layers);
        local_stats.high_product.pointwise_ciphertext_multiplications +=
            group_stats.pointwise_ciphertext_multiplications;
        low_inputs.push_back(ModulusSwitch(high_ring, low_ring, high_product));
        local_stats.modulus_boundaries++;
    }
    if (low_inputs.size() > low_window)
        throw std::invalid_argument("low RNS window cannot combine all inputs");
    const RNSCiphertext result = ProductTree(
        low_ring, low_qh_secret, std::move(low_inputs),
        &local_stats.low_product);
    if (stats != nullptr) *stats = local_stats;
    return result;
}

// General QH-SS form of Algorithm 3's delayed-modulus-switch schedule. Each
// window is multiplied independently at its current modulus; all resulting
// ciphertexts are then switched to the next ring. The last window must reduce
// the collection to one ciphertext. This realizes schedules such as the
// screened 8,2,2,2 tree with three RNS boundaries.
inline RNSCiphertext Algorithm3QHMultiStageProductRNS(
    const std::span<const RNSRing> rings,
    const std::span<const RNSSecret> secrets,
    std::vector<RNSCiphertext> factors,
    const std::span<const std::uint32_t> windows,
    const std::uint64_t plaintext_modulus = 2,
    MultiStageStats *stats = nullptr)
{
    if (rings.empty() || rings.size() != secrets.size() ||
        rings.size() != windows.size() || factors.empty())
        throw std::invalid_argument("invalid Algorithm 3 multi-stage schedule");
    const std::size_t degree = rings.front().degree();
    for (std::size_t stage = 0; stage < rings.size(); stage++) {
        if (rings[stage].degree() != degree || windows[stage] < 2 ||
            (windows[stage] & (windows[stage] - 1)) != 0)
            throw std::invalid_argument("invalid Algorithm 3 multi-stage ring/window");
        validate(rings[stage], secrets[stage]);
    }
    MultiStageStats local;
    local.products.resize(rings.size());
    local.modulus_boundaries.assign(rings.size() - 1, 0);
    for (std::size_t stage = 0; stage < rings.size(); stage++) {
        std::vector<RNSCiphertext> outputs;
        outputs.reserve((factors.size() + windows[stage] - 1) /
                        windows[stage]);
        for (std::size_t first = 0; first < factors.size();
             first += windows[stage]) {
            const std::size_t last = std::min<std::size_t>(
                first + windows[stage], factors.size());
            std::vector<RNSCiphertext> group;
            group.reserve(last - first);
            for (std::size_t i = first; i < last; i++)
                group.push_back(std::move(factors[i]));
            ProductTreeStats group_stats;
            RNSCiphertext product = ProductTree(
                rings[stage], secrets[stage], std::move(group), &group_stats);
            local.products[stage].pointwise_multiplication_layers = std::max(
                local.products[stage].pointwise_multiplication_layers,
                group_stats.pointwise_multiplication_layers);
            local.products[stage].pointwise_ciphertext_multiplications +=
                group_stats.pointwise_ciphertext_multiplications;
            if (stage + 1 < rings.size()) {
                if (!IsSameRNSBasis(rings[stage], rings[stage + 1])) {
                    product = ModulusSwitch(
                        rings[stage], rings[stage + 1], product,
                        plaintext_modulus);
                    local.modulus_boundaries[stage]++;
                }
            }
            outputs.push_back(std::move(product));
        }
        factors = std::move(outputs);
    }
    if (factors.size() != 1)
        throw std::invalid_argument(
            "Algorithm 3 multi-stage windows did not reduce to one output");
    if (stats != nullptr) *stats = std::move(local);
    return std::move(factors.front());
}

template <class PBCKey>
inline RNSCiphertext Algorithm3PBCQHMultiStageSelectorProductRNS(
    const std::span<const RNSRing> rings,
    const std::span<const RNSSecret> secrets,
    const PBCSchedule &schedule,
    const PBCKey &bootstrapping_key,
    const std::span<const std::uint32_t> lwe_a,
    const std::uint32_t lwe_modulus,
    const std::span<const std::uint32_t> windows,
    const std::uint64_t plaintext_modulus = 2,
    MultiStageStats *stats = nullptr)
{
    if (rings.empty() || lwe_modulus == 0 ||
        (2 * rings.front().degree()) % lwe_modulus != 0 ||
        lwe_a.size() != schedule.source_dimension)
        throw std::invalid_argument(
            "invalid Algorithm 3 multi-stage PBC selector input");
    const RNSRing &high_ring = rings.front();
    const std::size_t scale = 2 * high_ring.degree() / lwe_modulus;
    // Populate the shared monomial caches before buckets read them in
    // parallel. There are at most lwe_modulus distinct spectra per prime.
    for (const std::uint32_t value : lwe_a)
        for (std::size_t level = 0; level < high_ring.levels(); level++)
            static_cast<void>(
                high_ring[level].cachedMonomialNTT(value * scale));
    std::vector<RNSCiphertext> factors(schedule.bucket_count);
#pragma omp parallel for if (schedule.bucket_count > 1)
    for (std::int64_t bucket = 0;
         bucket < static_cast<std::int64_t>(schedule.bucket_count); bucket++)
        factors[static_cast<std::size_t>(bucket)] = PBCPreparedFactor(
            high_ring, schedule, bootstrapping_key,
            static_cast<std::uint32_t>(bucket), lwe_a, scale);
    return Algorithm3QHMultiStageProductRNS(
        rings, secrets, std::move(factors), windows, plaintext_modulus, stats);
}

// Convert BGV/LSB encoding m+t*e to an MSB encoding without multiplying the
// error by the MSB scale. For t=2 and odd Q, t^{-1}(m+te) has message at Q/2
// and error e. More generally the plaintext is multiplied by -Q^{-1} mod t.
inline RNSCiphertext LSBToMSB(const RNSRing &ring,
                              const RNSCiphertext &input,
                              const std::uint64_t plaintext_modulus = 2)
{
    if (plaintext_modulus < 2 || input.residues.size() != ring.levels())
        throw std::invalid_argument("invalid LSB-to-MSB conversion input");
    RNSCiphertext result = input;
    for (std::size_t level = 0; level < ring.levels(); level++) {
        const std::uint64_t modulus = ring[level].modulus();
        const std::uint64_t inverse = modular_ntt::invert(
            plaintext_modulus % modulus, modulus);
        for (std::size_t i = 0; i < ring.degree(); i++) {
            result.residues[level].a[i] = modular_ntt::multiply(
                result.residues[level].a[i], inverse, modulus);
            result.residues[level].b[i] = modular_ntt::multiply(
                result.residues[level].b[i], inverse, modulus);
        }
    }
    return result;
}

inline RNSCiphertext MultiplyByPlaintextNTT(
    const RNSRing &ring, const RNSCiphertext &input,
    const std::vector<std::vector<std::uint64_t>> &plaintext)
{
    if (input.residues.size() != ring.levels() ||
        plaintext.size() != ring.levels())
        throw std::invalid_argument("invalid RNS plaintext multiplication");
    RNSCiphertext result = input;
    for (std::size_t level = 0; level < ring.levels(); level++) {
        ring[level].requireSize(plaintext[level]);
#ifdef USE_HEXL
        intel::hexl::EltwiseMultMod(
            result.residues[level].a.data(), input.residues[level].a.data(),
            plaintext[level].data(), ring.degree(), ring[level].modulus(), 1);
        intel::hexl::EltwiseMultMod(
            result.residues[level].b.data(), input.residues[level].b.data(),
            plaintext[level].data(), ring.degree(), ring[level].modulus(), 1);
#else
        for (std::size_t i = 0; i < ring.degree(); i++) {
            result.residues[level].a[i] = modular_ntt::multiply(
                input.residues[level].a[i], plaintext[level][i],
                ring[level].modulus());
            result.residues[level].b[i] = modular_ntt::multiply(
                input.residues[level].b[i], plaintext[level][i],
                ring[level].modulus());
        }
#endif
    }
    return result;
}

inline std::vector<std::vector<std::uint64_t>>
Algorithm3BinarySignPlaintextLUT(const RNSRing &ring,
                                const std::uint32_t lwe_modulus)
{
    if (lwe_modulus == 0 || lwe_modulus > 2 * ring.degree() ||
        (2 * ring.degree()) % lwe_modulus != 0)
        throw std::invalid_argument("invalid Algorithm 3 plaintext LUT modulus");
    const std::size_t exponent_scale = 2 * ring.degree() / lwe_modulus;
    std::vector<std::vector<std::uint64_t>> result(
        ring.levels(), std::vector<std::uint64_t>(ring.degree()));
    for (std::uint32_t exponent = 0; exponent < lwe_modulus; exponent++) {
        const bool positive = exponent >= lwe_modulus / 2;
        const std::size_t ring_exponent = exponent * exponent_scale;
        for (std::size_t level = 0; level < ring.levels(); level++) {
            const std::uint64_t modulus = ring[level].modulus();
            const std::uint64_t desired = positive ? 1 : modulus - 1;
            if (ring_exponent == 0)
                result[level][0] = desired;
            else if (ring_exponent < ring.degree())
                result[level][ring.degree() - ring_exponent] =
                    modular_ntt::negate(desired, modulus);
            else if (ring_exponent == ring.degree())
                result[level][0] = modular_ntt::negate(desired, modulus);
            else
                result[level][2 * ring.degree() - ring_exponent] = desired;
        }
    }
    for (std::size_t level = 0; level < ring.levels(); level++)
        ring[level].forward(result[level]);
    return result;
}

template <class PBCKey>
inline RNSCiphertext Algorithm3PBCQHMultiStageBlindRotateRNS(
    const std::span<const RNSRing> rings,
    const std::span<const RNSSecret> secrets,
    const PBCSchedule &schedule,
    const PBCKey &bootstrapping_key, const LWEPhase &lwe,
    const std::vector<std::vector<std::uint64_t>> &plaintext_lut,
    const std::span<const std::uint32_t> windows,
    const std::uint64_t plaintext_modulus = 4,
    MultiStageStats *stats = nullptr)
{
    RNSCiphertext selectors = Algorithm3PBCQHMultiStageSelectorProductRNS(
        rings, secrets, schedule, bootstrapping_key, lwe.a, lwe.modulus,
        windows, plaintext_modulus, stats);
    const RNSRing &final_ring = rings.back();
    RNSCiphertext msb =
        LSBToMSB(final_ring, selectors, plaintext_modulus);
    std::vector<std::vector<std::uint64_t>> rotated = plaintext_lut;
    const std::size_t exponent_scale = 2 * final_ring.degree() / lwe.modulus;
    const std::size_t exponent =
        (lwe.modulus - lwe.b % lwe.modulus) * exponent_scale;
    for (std::size_t level = 0; level < final_ring.levels(); level++) {
        final_ring[level].requireSize(rotated[level]);
        const auto monomial = final_ring[level].monomialNTT(exponent);
        for (std::size_t i = 0; i < final_ring.degree(); i++)
            rotated[level][i] = modular_ntt::multiply(
                rotated[level][i], monomial[i], final_ring[level].modulus());
    }
    return MultiplyByPlaintextNTT(final_ring, msb, rotated);
}

template <class PBCKey>
inline LWEPhase Algorithm3PBCQHMultiStageBootstrap(
    const std::span<const RNSRing> rings,
    const std::span<const RNSSecret> secrets,
    const PBCSchedule &schedule,
    const PBCKey &bootstrapping_key, const LWEPhase &input,
    const std::vector<std::vector<std::uint64_t>> &plaintext_lut,
    const std::span<const std::uint32_t> windows,
    const LWEKeySwitchKey &output_key_switch,
    const std::uint32_t output_modulus,
    const std::uint64_t key_switch_modulus = 1U << 14,
    const std::uint64_t plaintext_modulus = 4,
    MultiStageStats *stats = nullptr)
{
    if (rings.empty() || rings.back().levels() != 1)
        throw std::invalid_argument(
            "Algorithm 3 output extraction requires a single-prime final ring");
    const RNSCiphertext accumulator =
        Algorithm3PBCQHMultiStageBlindRotateRNS(
            rings, secrets, schedule, bootstrapping_key, input,
            plaintext_lut, windows, plaintext_modulus, stats);
    const LWECiphertext extracted =
        SampleExtract(rings.back()[0], accumulator.residues[0]);
    const LWECiphertext boundary =
        ModulusSwitchCiphertext(extracted, key_switch_modulus);
    return ModulusSwitch(LWEKeySwitch(boundary, output_key_switch),
                         output_modulus);
}

inline Ciphertext Algorithm3PBCQHTwoStageBlindRotate(
    const RNSRing &high_ring, const RNSSecret &high_qh_secret,
    const PBCSchedule &schedule,
    const RNSPBCBootstrappingKey &bootstrapping_key, const LWEPhase &lwe,
    const std::vector<std::vector<std::uint64_t>> &test_vector,
    const std::uint32_t high_window, const std::uint32_t low_window,
    const Ring &low_ring, const Secret &low_qh_secret,
    TwoStageStats *stats = nullptr)
{
    if (lwe.modulus == 0 || (2 * high_ring.degree()) % lwe.modulus != 0 ||
        lwe.a.size() != schedule.source_dimension ||
        test_vector.size() != high_ring.levels())
        throw std::invalid_argument("invalid Algorithm 3 QH PBC input");
    const std::size_t scale = 2 * high_ring.degree() / lwe.modulus;
    std::vector<std::vector<std::uint64_t>> initial_message = test_vector;
    const std::size_t initial_exponent =
        (lwe.modulus - (lwe.b % lwe.modulus)) % lwe.modulus;
    for (std::size_t level = 0; level < high_ring.levels(); level++) {
        const auto monomial =
            high_ring[level].monomialNTT(initial_exponent * scale);
        for (std::size_t i = 0; i < high_ring.degree(); i++)
            initial_message[level][i] = modular_ntt::multiply(
                initial_message[level][i], monomial[i],
                high_ring[level].modulus());
    }
    const RNSCiphertext initial =
        TrivialEncryptNTT(high_ring, initial_message);
    std::vector<RNSCiphertext> factors;
    factors.reserve(schedule.bucket_count);
    for (std::uint32_t bucket = 0; bucket < schedule.bucket_count; bucket++) {
        RNSCiphertext factor = PBCPreparedFactor(
            high_ring, schedule, bootstrapping_key, bucket, lwe.a, scale);
        if (bucket == 0)
            for (std::size_t level = 0; level < high_ring.levels(); level++)
                for (std::size_t i = 0; i < high_ring.degree(); i++) {
                    const std::uint64_t modulus = high_ring[level].modulus();
                    factor.residues[level].a[i] = modular_ntt::multiply(
                        factor.residues[level].a[i],
                        initial.residues[level].b[i], modulus);
                    factor.residues[level].b[i] = modular_ntt::multiply(
                        factor.residues[level].b[i],
                        initial.residues[level].b[i], modulus);
                }
        factors.push_back(std::move(factor));
    }
    return Algorithm3QHTwoStageProduct(high_ring, high_qh_secret,
                                        std::move(factors), high_window,
                                        low_window, low_ring, low_qh_secret,
                                        stats);
}

inline RNSCiphertext Algorithm3PBCQHTwoStageBlindRotateRNS(
    const RNSRing &high_ring, const RNSSecret &high_qh_secret,
    const PBCSchedule &schedule,
    const RNSPBCBootstrappingKey &bootstrapping_key, const LWEPhase &lwe,
    const std::vector<std::vector<std::uint64_t>> &test_vector,
    const std::uint32_t high_window, const std::uint32_t low_window,
    const RNSRing &low_ring, const RNSSecret &low_qh_secret,
    TwoStageStats *stats = nullptr)
{
    if (lwe.modulus == 0 || (2 * high_ring.degree()) % lwe.modulus != 0 ||
        lwe.a.size() != schedule.source_dimension ||
        test_vector.size() != high_ring.levels())
        throw std::invalid_argument("invalid Algorithm 3 QH RNS PBC input");
    const std::size_t scale = 2 * high_ring.degree() / lwe.modulus;
    std::vector<std::vector<std::uint64_t>> initial_message = test_vector;
    const std::size_t initial_exponent =
        (lwe.modulus - (lwe.b % lwe.modulus)) % lwe.modulus;
    for (std::size_t level = 0; level < high_ring.levels(); level++) {
        const auto monomial =
            high_ring[level].monomialNTT(initial_exponent * scale);
        for (std::size_t i = 0; i < high_ring.degree(); i++)
            initial_message[level][i] = modular_ntt::multiply(
                initial_message[level][i], monomial[i],
                high_ring[level].modulus());
    }
    const RNSCiphertext initial =
        TrivialEncryptNTT(high_ring, initial_message);
    std::vector<RNSCiphertext> factors;
    factors.reserve(schedule.bucket_count);
    for (std::uint32_t bucket = 0; bucket < schedule.bucket_count; bucket++) {
        RNSCiphertext factor = PBCPreparedFactor(
            high_ring, schedule, bootstrapping_key, bucket, lwe.a, scale);
        if (bucket == 0)
            for (std::size_t level = 0; level < high_ring.levels(); level++)
                for (std::size_t i = 0; i < high_ring.degree(); i++) {
                    const std::uint64_t modulus = high_ring[level].modulus();
                    factor.residues[level].a[i] = modular_ntt::multiply(
                        factor.residues[level].a[i],
                        initial.residues[level].b[i], modulus);
                    factor.residues[level].b[i] = modular_ntt::multiply(
                        factor.residues[level].b[i],
                        initial.residues[level].b[i], modulus);
                }
        factors.push_back(std::move(factor));
    }
    return Algorithm3QHTwoStageProductRNS(
        high_ring, high_qh_secret, std::move(factors), high_window,
        low_window, low_ring, low_qh_secret, stats);
}

inline std::vector<std::vector<std::uint64_t>> Algorithm3BinarySignLUT(
    const RNSRing &ring, const std::uint32_t lwe_modulus)
{
    if (lwe_modulus == 0 || lwe_modulus > 2 * ring.degree() ||
        (2 * ring.degree()) % lwe_modulus != 0)
        throw std::invalid_argument("invalid Algorithm 3 binary LUT modulus");
    boost::multiprecision::cpp_int high_modulus = 1;
    for (std::size_t level = 0; level < ring.levels(); level++)
        high_modulus *= ring[level].modulus();
    const boost::multiprecision::cpp_int plaintext_scale = high_modulus / 8;
    const std::size_t exponent_scale = 2 * ring.degree() / lwe_modulus;
    std::vector<std::vector<std::uint64_t>> result(
        ring.levels(), std::vector<std::uint64_t>(ring.degree()));
    for (std::uint32_t exponent = 0; exponent < lwe_modulus; exponent++) {
        const bool positive = exponent >= lwe_modulus / 2;
        const std::size_t ring_exponent = exponent * exponent_scale;
        for (std::size_t level = 0; level < ring.levels(); level++) {
            const std::uint64_t modulus = ring[level].modulus();
            const std::uint64_t encoded = (plaintext_scale % modulus)
                                              .convert_to<std::uint64_t>();
            const std::uint64_t desired = positive
                ? encoded
                : modular_ntt::negate(encoded, modulus);
            if (ring_exponent == 0)
                result[level][0] = desired;
            else if (ring_exponent < ring.degree())
                result[level][ring.degree() - ring_exponent] =
                    modular_ntt::negate(desired, modulus);
            else if (ring_exponent == ring.degree())
                result[level][0] = modular_ntt::negate(desired, modulus);
            else
                result[level][2 * ring.degree() - ring_exponent] = desired;
        }
    }
    for (std::size_t level = 0; level < ring.levels(); level++)
        ring[level].forward(result[level]);
    return result;
}

template <class BoundaryKey>
inline LWEPhase Algorithm3Bootstrap(
    const RNSRing &high_ring, const RNSSecret &binary_ntt_secret,
    const RNSBootstrappingKey &bootstrapping_key, const LWEPhase &input,
    const std::vector<std::vector<std::uint64_t>> &test_vector,
    const std::uint32_t high_window, const std::uint32_t low_window,
    const BoundaryKey &boundary_key, const RNSSecret &high_small_secret,
    const Ring &low_ring, const Secret &low_small_secret,
    const LWEKeySwitchKey &output_key_switch,
    const std::uint32_t output_modulus, TwoStageStats *stats = nullptr)
{
    const Ciphertext accumulator = Algorithm3TwoStageBlindRotate(
        high_ring, binary_ntt_secret, bootstrapping_key, input, test_vector,
        high_window, low_window, boundary_key, high_small_secret, low_ring,
        low_small_secret, stats);
    const LWECiphertext extracted = SampleExtract(low_ring, accumulator);
    const LWECiphertext switched = LWEKeySwitch(extracted, output_key_switch);
    return ModulusSwitch(switched, output_modulus);
}

// PBC variant of Algorithm 3.  PBC compiles a general sparse source key into
// `schedule.bucket_count` one-hot factors before the same two-stage product,
// extraction, and output LWE key switch used by the structured path.
template <class BoundaryKey>
inline LWEPhase Algorithm3PBCBootstrap(
    const RNSRing &high_ring, const RNSSecret &binary_ntt_secret,
    const PBCSchedule &schedule,
    const RNSPBCBootstrappingKey &bootstrapping_key,
    const LWEPhase &input,
    const std::vector<std::vector<std::uint64_t>> &test_vector,
    const std::uint32_t high_window, const std::uint32_t low_window,
    const BoundaryKey &boundary_key, const RNSSecret &high_small_secret,
    const Ring &low_ring, const Secret &low_small_secret,
    const LWEKeySwitchKey &output_key_switch,
    const std::uint32_t output_modulus, TwoStageStats *stats = nullptr)
{
    const Ciphertext accumulator = Algorithm3PBCTwoStageBlindRotate(
        high_ring, binary_ntt_secret, schedule, bootstrapping_key, input,
        test_vector, high_window, low_window, boundary_key,
        high_small_secret, low_ring, low_small_secret, stats);
    const LWECiphertext extracted = SampleExtract(low_ring, accumulator);
    const LWECiphertext switched = LWEKeySwitch(extracted, output_key_switch);
    return ModulusSwitch(switched, output_modulus);
}

inline LWEPhase Algorithm3PBCQHBootstrap(
    const RNSRing &high_ring, const RNSSecret &high_qh_secret,
    const PBCSchedule &schedule,
    const RNSPBCBootstrappingKey &bootstrapping_key, const LWEPhase &input,
    const std::vector<std::vector<std::uint64_t>> &test_vector,
    const std::uint32_t high_window, const std::uint32_t low_window,
    const Ring &low_ring, const Secret &low_qh_secret,
    const LWEKeySwitchKey &output_key_switch,
    const std::uint32_t output_modulus, TwoStageStats *stats = nullptr)
{
    const Ciphertext accumulator = Algorithm3PBCQHTwoStageBlindRotate(
        high_ring, high_qh_secret, schedule, bootstrapping_key, input,
        test_vector, high_window, low_window, low_ring, low_qh_secret, stats);
    const LWECiphertext extracted = SampleExtract(low_ring, accumulator);
    return ModulusSwitch(LWEKeySwitch(extracted, output_key_switch),
                         output_modulus);
}

}  // namespace TFHEpp::shallowboot::lowdepth
