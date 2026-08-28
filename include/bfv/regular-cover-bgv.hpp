#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace TFHEpp::regular_cover_bgv {

// Proof-aligned executable prototype.  The base ring is stored in split/NTT
// coordinates, so multiplication is coordinatewise and the cyclic Galois
// action is an explicit coordinate permutation.
template <std::size_t N, std::size_t GroupSize, std::uint64_t Modulus,
          std::uint64_t PlaintextModulus>
struct Scheme {
    static_assert(N > 0 && GroupSize > 0 && N % GroupSize == 0);
    static_assert(Modulus > 2 && PlaintextModulus >= 2);
    static constexpr std::size_t split_degree = N;
    static constexpr std::size_t group_size = GroupSize;

    using Base = std::array<std::uint64_t, N>;
    using Element = std::array<Base, GroupSize>;

    struct Ciphertext {
        Element mask{};
        Element body{};
    };

    struct BaseCiphertext {
        Base mask{};
        Base body{};
    };

    struct SwitchKey {
        std::size_t automorphism{};
        std::uint64_t gadget_base{};
        std::vector<Ciphertext> rows{};
    };

    struct Completion {
        Element multiplier{};
        Element offset{};
        Element witness{};
    };

    struct PublicHint {
        Element alpha{};
        Element beta{};
        Element u{};
        Element v{};
    };

    struct KeyMaterial {
        Element operational_secret{};
        Completion completion{};
        PublicHint hint{};
    };

    static constexpr std::uint64_t add(const std::uint64_t left,
                                       const std::uint64_t right)
    {
        const std::uint64_t sum = left + right;
        return sum >= Modulus ? sum - Modulus : sum;
    }

    static constexpr std::uint64_t sub(const std::uint64_t left,
                                       const std::uint64_t right)
    {
        return left >= right ? left - right : left + Modulus - right;
    }

    static constexpr std::uint64_t neg(const std::uint64_t value)
    {
        return value == 0 ? 0 : Modulus - value;
    }

    static constexpr std::uint64_t mul(const std::uint64_t left,
                                       const std::uint64_t right)
    {
        return static_cast<std::uint64_t>(
            (static_cast<unsigned __int128>(left) * right) % Modulus);
    }

    static std::uint64_t pow(std::uint64_t value, std::uint64_t exponent)
    {
        std::uint64_t result = 1;
        while (exponent != 0) {
            if ((exponent & 1) != 0) result = mul(result, value);
            value = mul(value, value);
            exponent >>= 1;
        }
        return result;
    }

    static std::uint64_t inverse(const std::uint64_t value)
    {
        if (value == 0) throw std::invalid_argument("zero is not a unit");
        // The prototype profiles use a prime modulus.
        return pow(value, Modulus - 2);
    }

    static std::uint64_t fromSigned(const std::int64_t value)
    {
        if (value >= 0) return static_cast<std::uint64_t>(value) % Modulus;
        const std::uint64_t magnitude =
            static_cast<std::uint64_t>(-value) % Modulus;
        return magnitude == 0 ? 0 : Modulus - magnitude;
    }

    static std::int64_t centered(const std::uint64_t value)
    {
        return value <= Modulus / 2
                   ? static_cast<std::int64_t>(value)
                   : -static_cast<std::int64_t>(Modulus - value);
    }

    static Element zero() { return {}; }

    static Element constant(const std::uint64_t value)
    {
        Element result{};
        for (auto &component : result) component.fill(value % Modulus);
        return result;
    }

    static Element elementAdd(const Element &left, const Element &right)
    {
        Element result{};
        for (std::size_t h = 0; h < GroupSize; ++h)
            for (std::size_t j = 0; j < N; ++j)
                result[h][j] = add(left[h][j], right[h][j]);
        return result;
    }

    static Element elementSub(const Element &left, const Element &right)
    {
        Element result{};
        for (std::size_t h = 0; h < GroupSize; ++h)
            for (std::size_t j = 0; j < N; ++j)
                result[h][j] = sub(left[h][j], right[h][j]);
        return result;
    }

    static Element elementNeg(const Element &value)
    {
        Element result{};
        for (std::size_t h = 0; h < GroupSize; ++h)
            for (std::size_t j = 0; j < N; ++j) result[h][j] = neg(value[h][j]);
        return result;
    }

    static Element elementMul(const Element &left, const Element &right)
    {
        Element result{};
        for (std::size_t h = 0; h < GroupSize; ++h)
            for (std::size_t j = 0; j < N; ++j)
                result[h][j] = mul(left[h][j], right[h][j]);
        return result;
    }

    static Element scalarMul(const Element &value, const std::uint64_t scalar)
    {
        return elementMul(value, constant(scalar));
    }

    // Base automorphism action: cyclically permute split coordinates.
    static Base baseAutomorphism(const Base &value, const std::size_t sigma)
    {
        Base result{};
        for (std::size_t j = 0; j < N; ++j)
            result[j] = value[(j + N - (sigma % N)) % N];
        return result;
    }

    // Regular lifted action: sigma(x[sigma^-1 h]).
    static Element liftedAutomorphism(const Element &value,
                                      const std::size_t sigma)
    {
        Element result{};
        const std::size_t shift = sigma % GroupSize;
        for (std::size_t h = 0; h < GroupSize; ++h) {
            const std::size_t source = (h + GroupSize - shift) % GroupSize;
            result[h] = baseAutomorphism(value[source], sigma);
        }
        return result;
    }

    static Element fixedEmbedding(const Base &value)
    {
        Element result{};
        for (std::size_t h = 0; h < GroupSize; ++h)
            result[h] = baseAutomorphism(value, h);
        return result;
    }

    static Element assembleCover(const std::array<Base, GroupSize> &values)
    {
        Element result{};
        for (std::size_t h = 0; h < GroupSize; ++h)
            result[h] = baseAutomorphism(values[h], h);
        return result;
    }

    static Ciphertext assembleCoverSamples(
        const std::array<BaseCiphertext, GroupSize> &samples)
    {
        std::array<Base, GroupSize> masks{}, bodies{};
        for (std::size_t h = 0; h < GroupSize; ++h) {
            masks[h] = samples[h].mask;
            bodies[h] = samples[h].body;
        }
        return {assembleCover(masks), assembleCover(bodies)};
    }

    template <class Engine>
    static Completion makeCompletion(Engine &engine)
    {
        std::uniform_int_distribution<int> bit(0, 1);
        Base bottom{};
        for (auto &entry : bottom)
            entry = static_cast<std::uint64_t>(bit(engine));

        Completion completion{};
        const Element invariant = fixedEmbedding(bottom);
        for (std::size_t h = 0; h < GroupSize; ++h) {
            for (std::size_t j = 0; j < N; ++j) {
                // Orbit representative is d=j-h.  The identity component has
                // p=0, every other component gets one independent coin.
                const bool coin = h == 0 ? false : bit(engine) != 0;
                completion.multiplier[h][j] = coin ? Modulus - 1 : 1;
                completion.offset[h][j] = coin ? 1 : 0;
                completion.witness[h][j] =
                    add(mul(completion.multiplier[h][j], invariant[h][j]),
                        completion.offset[h][j]);
            }
        }
        return completion;
    }

    template <class Engine>
    static KeyMaterial keygenWithOperationalSecret(
        Engine &engine, const Element &operational_secret)
    {
        std::uniform_int_distribution<std::uint64_t> unit(1, Modulus - 1);
        KeyMaterial key{};
        key.completion = makeCompletion(engine);
        key.operational_secret = operational_secret;
        for (std::size_t h = 0; h < GroupSize; ++h) {
            for (std::size_t j = 0; j < N; ++j) {
                const std::uint64_t z = operational_secret[h][j];
                const std::uint64_t alpha = unit(engine);
                const std::uint64_t beta =
                    add(mul(alpha, key.completion.witness[h][j]), z);
                key.hint.alpha[h][j] = alpha;
                key.hint.beta[h][j] = beta;
                key.hint.u[h][j] = sub(add(beta, beta), alpha);
                key.hint.v[h][j] = sub(mul(alpha, beta), mul(beta, beta));
            }
        }
        return key;
    }

    template <class Engine>
    static KeyMaterial keygen(Engine &engine, const int secret_bound = 1)
    {
        if (secret_bound < 0)
            throw std::invalid_argument("negative secret bound");
        std::uniform_int_distribution<int> small(-secret_bound, secret_bound);
        Element operational_secret{};
        for (std::size_t h = 0; h < GroupSize; ++h) {
            for (std::size_t j = 0; j < N; ++j)
                operational_secret[h][j] = fromSigned(small(engine));
        }
        return keygenWithOperationalSecret(engine, operational_secret);
    }

    static bool checkQuadraticHint(const KeyMaterial &key)
    {
        return elementMul(key.operational_secret, key.operational_secret) ==
               elementAdd(elementMul(key.hint.u, key.operational_secret),
                          key.hint.v);
    }

    static Element phase(const Ciphertext &ciphertext, const Element &secret)
    {
        return elementSub(ciphertext.body, elementMul(ciphertext.mask, secret));
    }

    static Ciphertext compileWitnessRow(
        const Ciphertext &source, const Completion &completion,
        const std::array<Element, GroupSize> &automorphism_weights)
    {
        Element multiplier_shift{};
        Element offset_shift{};
        for (std::size_t sigma = 0; sigma < GroupSize; ++sigma) {
            multiplier_shift = elementAdd(
                multiplier_shift,
                elementMul(automorphism_weights[sigma],
                           liftedAutomorphism(completion.multiplier, sigma)));
            offset_shift = elementAdd(
                offset_shift,
                elementMul(automorphism_weights[sigma],
                           liftedAutomorphism(completion.offset, sigma)));
        }
        // H is involutive in the exact Binary-NTT completion.
        const Element target_mask = elementMul(
            elementSub(source.mask, multiplier_shift), completion.multiplier);
        const Element target_body = elementAdd(
            elementAdd(source.body, elementMul(target_mask, completion.offset)),
            offset_shift);
        return {target_mask, target_body};
    }

    template <class Engine>
    static Ciphertext encryptElement(const Element &message,
                                     const Element &secret, Engine &engine,
                                     const int error_bound = 1)
    {
        std::uniform_int_distribution<std::uint64_t> uniform(0, Modulus - 1);
        std::uniform_int_distribution<int> error(-error_bound, error_bound);
        Ciphertext result{};
        for (std::size_t h = 0; h < GroupSize; ++h) {
            for (std::size_t j = 0; j < N; ++j) {
                result.mask[h][j] = uniform(engine);
                const std::uint64_t scaled_error =
                    mul(PlaintextModulus % Modulus, fromSigned(error(engine)));
                result.body[h][j] = add(
                    add(mul(result.mask[h][j], secret[h][j]), message[h][j]),
                    scaled_error);
            }
        }
        return result;
    }

    template <class Engine>
    static Ciphertext encryptScalar(const std::uint64_t message,
                                    const Element &secret, Engine &engine,
                                    const int error_bound = 1)
    {
        return encryptElement(constant(message % PlaintextModulus), secret,
                              engine, error_bound);
    }

    static std::uint64_t decryptScalar(const Ciphertext &ciphertext,
                                       const Element &secret)
    {
        const Element decoded_phase = phase(ciphertext, secret);
        const std::int64_t value = centered(decoded_phase[0][0]);
        const std::int64_t modulus =
            static_cast<std::int64_t>(PlaintextModulus);
        const std::int64_t reduced = ((value % modulus) + modulus) % modulus;
        return static_cast<std::uint64_t>(reduced);
    }

    static Ciphertext ciphertextAdd(const Ciphertext &left,
                                    const Ciphertext &right)
    {
        return {elementAdd(left.mask, right.mask),
                elementAdd(left.body, right.body)};
    }

    static Ciphertext addConstant(const Ciphertext &ciphertext,
                                  const std::uint64_t value)
    {
        Ciphertext result = ciphertext;
        result.body = elementAdd(result.body, constant(value));
        return result;
    }

    static Ciphertext plaintextMul(const Ciphertext &ciphertext,
                                   const Element &plaintext)
    {
        return {elementMul(ciphertext.mask, plaintext),
                elementMul(ciphertext.body, plaintext)};
    }

    static Ciphertext quadraticHintMul(const Ciphertext &left,
                                       const Ciphertext &right,
                                       const PublicHint &hint)
    {
        const Element a1b2 = elementMul(left.mask, right.body);
        const Element a2b1 = elementMul(right.mask, left.body);
        const Element a1a2 = elementMul(left.mask, right.mask);
        Ciphertext result{};
        result.mask =
            elementSub(elementAdd(a1b2, a2b1), elementMul(a1a2, hint.u));
        result.body = elementAdd(elementMul(left.body, right.body),
                                 elementMul(a1a2, hint.v));
        return result;
    }

    static Ciphertext liftedAutomorphism(const Ciphertext &ciphertext,
                                         const std::size_t sigma)
    {
        return {liftedAutomorphism(ciphertext.mask, sigma),
                liftedAutomorphism(ciphertext.body, sigma)};
    }

    template <class Engine>
    static SwitchKey makeSwitchKey(const Element &secret,
                                   const std::size_t sigma,
                                   const std::uint64_t gadget_base,
                                   Engine &engine, const int error_bound = 1)
    {
        if (gadget_base < 2) throw std::invalid_argument("invalid gadget base");
        SwitchKey key{sigma % GroupSize, gadget_base, {}};
        const Element automorphed_secret =
            liftedAutomorphism(secret, key.automorphism);
        std::uint64_t gadget = 1;
        while (gadget < Modulus) {
            key.rows.push_back(
                encryptElement(scalarMul(automorphed_secret, gadget), secret,
                               engine, error_bound));
            if (gadget > (Modulus - 1) / gadget_base) break;
            gadget *= gadget_base;
        }
        return key;
    }

    static Ciphertext applySwitchKey(const Ciphertext &ciphertext,
                                     const SwitchKey &key)
    {
        const Ciphertext automorphed =
            liftedAutomorphism(ciphertext, key.automorphism);
        Element mask_sum{};
        Element body_sum{};
        Element remaining = automorphed.mask;
        for (const auto &row : key.rows) {
            Element digit{};
            for (std::size_t h = 0; h < GroupSize; ++h) {
                for (std::size_t j = 0; j < N; ++j) {
                    digit[h][j] = remaining[h][j] % key.gadget_base;
                    remaining[h][j] /= key.gadget_base;
                }
            }
            mask_sum = elementAdd(mask_sum, elementMul(digit, row.mask));
            body_sum = elementAdd(body_sum, elementMul(digit, row.body));
        }
        return {elementNeg(mask_sum), elementSub(automorphed.body, body_sum)};
    }

    static std::vector<std::uint64_t> interpolationPolynomial(
        const int error_bound)
    {
        if (error_bound < 0)
            throw std::invalid_argument("negative error bound");
        const int minimum = -static_cast<int>(PlaintextModulus) * error_bound;
        const int maximum = static_cast<int>(PlaintextModulus) * error_bound +
                            static_cast<int>(PlaintextModulus) - 1;
        const std::size_t count =
            static_cast<std::size_t>(maximum - minimum + 1);
        if (count >= Modulus)
            throw std::invalid_argument("interpolation domain wraps modulus");

        std::vector<std::uint64_t> polynomial(count, 0);
        for (int point = minimum; point <= maximum; ++point) {
            std::vector<std::uint64_t> basis{1};
            std::uint64_t denominator = 1;
            const std::uint64_t xpoint = fromSigned(point);
            for (int other = minimum; other <= maximum; ++other) {
                if (other == point) continue;
                const std::uint64_t xother = fromSigned(other);
                std::vector<std::uint64_t> next(basis.size() + 1, 0);
                for (std::size_t degree = 0; degree < basis.size(); ++degree) {
                    next[degree] =
                        add(next[degree], mul(basis[degree], neg(xother)));
                    next[degree + 1] = add(next[degree + 1], basis[degree]);
                }
                basis = std::move(next);
                denominator = mul(denominator, sub(xpoint, xother));
            }
            const int reduced = ((point % static_cast<int>(PlaintextModulus)) +
                                 static_cast<int>(PlaintextModulus)) %
                                static_cast<int>(PlaintextModulus);
            const std::uint64_t scale =
                mul(static_cast<std::uint64_t>(reduced), inverse(denominator));
            for (std::size_t degree = 0; degree < basis.size(); ++degree)
                polynomial[degree] =
                    add(polynomial[degree], mul(scale, basis[degree]));
        }
        return polynomial;
    }

    static Ciphertext polynomialEval(
        const std::vector<std::uint64_t> &polynomial, const Ciphertext &input,
        const PublicHint &hint)
    {
        if (polynomial.empty()) return {};
        Ciphertext accumulator{};
        accumulator.body = constant(polynomial.back());
        for (std::size_t index = polynomial.size() - 1; index > 0; --index) {
            accumulator = quadraticHintMul(accumulator, input, hint);
            accumulator = addConstant(accumulator, polynomial[index - 1]);
        }
        return accumulator;
    }

    template <class Engine>
    static SwitchKey makeBootKey(const Element &secret, Engine &engine,
                                 const std::uint64_t gadget_base = 2,
                                 const int error_bound = 1)
    {
        // Identity-automorphism gadget rows encrypt g_i*z directly with fresh
        // uniform masks and genuine small errors.
        return makeSwitchKey(secret, 0, gadget_base, engine, error_bound);
    }

    template <class Engine>
    static Ciphertext refresh(
        const Ciphertext &input, const SwitchKey &boot_key,
        const Element &secret, const PublicHint &hint,
        const std::vector<std::uint64_t> &digit_polynomial, Engine &engine,
        const int fresh_error_bound = 1)
    {
        // Identity key switching is homomorphic decryption: the output phase
        // is B-A*z plus the bounded gadget-weighted boot-key errors.
        Ciphertext encrypted_phase = applySwitchKey(input, boot_key);
        Ciphertext refreshed =
            polynomialEval(digit_polynomial, encrypted_phase, hint);
        const Ciphertext fresh_zero =
            encryptScalar(0, secret, engine, fresh_error_bound);
        return ciphertextAdd(refreshed, fresh_zero);
    }
};

}  // namespace TFHEpp::regular_cover_bgv
