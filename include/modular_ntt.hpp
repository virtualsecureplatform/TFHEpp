#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace TFHEpp::modular_ntt {

struct PrimeModulus {
    std::uint64_t value;
    std::uint64_t primitive_root;
};

// These primes support power-of-two transforms through length 2^12 or
// greater and a length-17 transform.  They are deliberately scheme-neutral:
// callers select one, two, or three primes according to their proved exact
// coefficient bound.
inline constexpr std::array<PrimeModulus, 3> wide_primes{{
    {4611686018426953729ULL, 11},
    {4611686018426884097ULL, 5},
    {4611686018426257409ULL, 3},
}};

inline std::uint64_t add(const std::uint64_t lhs, const std::uint64_t rhs,
                         const std::uint64_t modulus)
{
    // Every supported modulus is below 2^63, so this addition cannot wrap.
    const std::uint64_t sum = lhs + rhs;
    return sum >= modulus ? sum - modulus : sum;
}

inline std::uint64_t subtract(const std::uint64_t lhs, const std::uint64_t rhs,
                              const std::uint64_t modulus)
{
    return lhs >= rhs ? lhs - rhs : lhs + modulus - rhs;
}

inline std::uint64_t negate(const std::uint64_t value,
                            const std::uint64_t modulus)
{
    return value == 0 ? 0 : modulus - value;
}

inline std::uint64_t reduceWide(const unsigned __int128 value,
                                const std::uint64_t modulus)
{
    // The exact GL primes are pseudo-Mersenne: p = 2^62-c with c < 2^21.
    // For any 128-bit x, fold the high 62-bit part twice using
    // 2^62 = c (mod p).  The second fold is below 2p, so one subtraction
    // canonicalizes it.  This avoids the compiler's much slower 128-by-64
    // division helper.  Retain the generic reduction for callers that provide
    // a different modulus.
    constexpr std::uint64_t base = std::uint64_t{1} << 62;
    constexpr std::uint64_t mask = base - 1;
    if (modulus < base && base - modulus < (std::uint64_t{1} << 21)) {
        const std::uint64_t complement = base - modulus;
        const unsigned __int128 first =
            (value & mask) + (value >> 62) * complement;
        const std::uint64_t second = static_cast<std::uint64_t>(
            (first & mask) + (first >> 62) * complement);
        return second >= modulus ? second - modulus : second;
    }
    return static_cast<std::uint64_t>(value % modulus);
}

inline std::uint64_t multiply(const std::uint64_t lhs, const std::uint64_t rhs,
                              const std::uint64_t modulus)
{
    return reduceWide(static_cast<unsigned __int128>(lhs) * rhs, modulus);
}

struct ShoupMultiplier {
    std::uint64_t value = 0;
    std::uint64_t quotient = 0;

    ShoupMultiplier() = default;

    ShoupMultiplier(const std::uint64_t constant, const std::uint64_t modulus)
        : value(constant)
    {
        if (modulus == 0 || constant >= modulus ||
            modulus >= (std::uint64_t{1} << 63))
            throw std::invalid_argument("invalid Shoup multiplier");
        quotient = static_cast<std::uint64_t>(
            (static_cast<unsigned __int128>(constant) << 64) / modulus);
    }

    // input must be the canonical residue in [0, modulus).
    std::uint64_t apply(const std::uint64_t input,
                        const std::uint64_t modulus) const
    {
        std::uint64_t result = applyLazy(input, modulus);
        if (result >= modulus) result -= modulus;
        return result;
    }

    // For modulus < 2^62 this also accepts input in [0, 2*modulus) and
    // returns the congruent representative in that same lazy range.
    std::uint64_t applyLazy(const std::uint64_t input,
                            const std::uint64_t modulus) const
    {
        const std::uint64_t estimate = static_cast<std::uint64_t>(
            (static_cast<unsigned __int128>(input) * quotient) >> 64);
        return input * value - estimate * modulus;
    }
};

inline std::uint64_t addLazy(const std::uint64_t lhs, const std::uint64_t rhs,
                             const std::uint64_t twice_modulus)
{
    const std::uint64_t sum = lhs + rhs;
    return sum >= twice_modulus ? sum - twice_modulus : sum;
}

inline std::uint64_t subtractLazy(const std::uint64_t lhs,
                                  const std::uint64_t rhs,
                                  const std::uint64_t twice_modulus)
{
    return lhs >= rhs ? lhs - rhs : lhs + twice_modulus - rhs;
}

inline std::uint64_t power(std::uint64_t base, std::uint64_t exponent,
                           const std::uint64_t modulus)
{
    std::uint64_t result = 1;
    while (exponent != 0) {
        if ((exponent & 1U) != 0) result = multiply(result, base, modulus);
        base = multiply(base, base, modulus);
        exponent >>= 1;
    }
    return result;
}

inline std::uint64_t invert(const std::uint64_t value,
                            const std::uint64_t modulus)
{
    if (value == 0) throw std::invalid_argument("zero has no modular inverse");
    return power(value, modulus - 2, modulus);
}

constexpr bool isPowerOfTwo(const std::size_t value)
{
    return value != 0 && (value & (value - 1)) == 0;
}

class Radix2NTTPlan {
public:
    Radix2NTTPlan(const std::size_t size, const PrimeModulus prime)
        : size_(size), modulus_(prime.value)
    {
        if (!isPowerOfTwo(size_) || size_ < 2)
            throw std::invalid_argument("NTT size must be a power of two");
        if (modulus_ >= (std::uint64_t{1} << 63) || (modulus_ - 1) % size_ != 0)
            throw std::invalid_argument("modulus does not support NTT size");

        inverse_size_ =
            ShoupMultiplier(invert(size_ % modulus_, modulus_), modulus_);
        for (std::size_t i = 1, reversed = 0; i < size_; i++) {
            std::size_t bit = size_ >> 1;
            for (; (reversed & bit) != 0; bit >>= 1) reversed ^= bit;
            reversed ^= bit;
            if (i < reversed) bit_reversal_swaps_.emplace_back(i, reversed);
        }
        for (std::size_t length = 2; length <= size_; length <<= 1) {
            const std::uint64_t root =
                power(prime.primitive_root, (modulus_ - 1) / length, modulus_);
            if (power(root, length, modulus_) != 1 ||
                power(root, length / 2, modulus_) == 1)
                throw std::invalid_argument("invalid primitive NTT root");
            const std::uint64_t inverse_root = invert(root, modulus_);
            const std::size_t half = length >> 1;
            forward_twiddles_.emplace_back(half);
            inverse_twiddles_.emplace_back(half);
            forward_twiddles_.back()[0] = ShoupMultiplier(1, modulus_);
            inverse_twiddles_.back()[0] = ShoupMultiplier(1, modulus_);
            for (std::size_t j = 1; j < half; j++) {
                forward_twiddles_.back()[j] = ShoupMultiplier(
                    multiply(forward_twiddles_.back()[j - 1].value, root,
                             modulus_),
                    modulus_);
                inverse_twiddles_.back()[j] = ShoupMultiplier(
                    multiply(inverse_twiddles_.back()[j - 1].value,
                             inverse_root, modulus_),
                    modulus_);
            }
        }
    }

    std::size_t size() const { return size_; }
    std::uint64_t modulus() const { return modulus_; }

    void forward(const std::span<std::uint64_t> values) const
    {
        transform(values, false);
    }

    void inverse(const std::span<std::uint64_t> values) const
    {
        transform(values, true, true);
    }

    // Convolution plans can combine this normalization with a following
    // coefficient-wise multiplier and save one modular product per value.
    void inverseUnscaled(const std::span<std::uint64_t> values) const
    {
        transform(values, true, false);
    }

private:
    void transform(const std::span<std::uint64_t> values, const bool invert,
                   const bool normalize = false) const
    {
        if (values.size() != size_)
            throw std::invalid_argument("NTT input has the wrong size");

        for (const auto &[first, second] : bit_reversal_swaps_)
            std::swap(values[first], values[second]);

        if (modulus_ < (std::uint64_t{1} << 62)) {
            const std::uint64_t twice_modulus = 2 * modulus_;
            std::size_t stage = 0;
            for (std::size_t length = 2; length <= size_;
                 length <<= 1, stage++) {
                const auto &twiddles = invert ? inverse_twiddles_[stage]
                                              : forward_twiddles_[stage];
                const std::size_t half = length >> 1;
                for (std::size_t block = 0; block < size_; block += length)
                    for (std::size_t j = 0; j < half; j++) {
                        const std::uint64_t even = values[block + j];
                        const std::uint64_t odd =
                            j == 0 ? values[block + half]
                                   : twiddles[j].applyLazy(
                                         values[block + j + half], modulus_);
                        values[block + j] = addLazy(even, odd, twice_modulus);
                        values[block + j + half] =
                            subtractLazy(even, odd, twice_modulus);
                    }
            }
            for (std::uint64_t &value : values)
                if (value >= modulus_) value -= modulus_;
            if (normalize)
                for (std::uint64_t &value : values)
                    value = inverse_size_.apply(value, modulus_);
            return;
        }

        std::size_t stage = 0;
        for (std::size_t length = 2; length <= size_; length <<= 1, stage++) {
            const auto &twiddles =
                invert ? inverse_twiddles_[stage] : forward_twiddles_[stage];
            const std::size_t half = length >> 1;
            for (std::size_t block = 0; block < size_; block += length) {
                for (std::size_t j = 0; j < half; j++) {
                    const std::uint64_t even = values[block + j];
                    const std::uint64_t odd =
                        j == 0 ? values[block + half]
                               : twiddles[j].apply(values[block + j + half],
                                                   modulus_);
                    values[block + j] = add(even, odd, modulus_);
                    values[block + j + half] = subtract(even, odd, modulus_);
                }
            }
        }

        if (normalize)
            for (std::uint64_t &value : values)
                value = inverse_size_.apply(value, modulus_);
    }

    std::size_t size_;
    std::uint64_t modulus_;
    ShoupMultiplier inverse_size_{};
    std::vector<std::pair<std::size_t, std::size_t>> bit_reversal_swaps_;
    std::vector<std::vector<ShoupMultiplier>> forward_twiddles_;
    std::vector<std::vector<ShoupMultiplier>> inverse_twiddles_;
};

// Compile-time-sized radix-2 plan for small transforms embedded in another
// fixed algorithm (notably Rader's p=17 convolution). Flat twiddle tables and
// compile-time stages let the compiler remove the dynamic loop and vector
// machinery used by the general plan.
template <std::size_t Size>
class FixedRadix2NTTPlan {
public:
    explicit FixedRadix2NTTPlan(const PrimeModulus prime)
        : modulus_(prime.value), inverse_size_(0)
    {
        static_assert(isPowerOfTwo(Size) && Size >= 2);
        if (modulus_ >= (std::uint64_t{1} << 63) || (modulus_ - 1) % Size != 0)
            throw std::invalid_argument("modulus does not support NTT size");

        inverse_size_ = invert(Size % modulus_, modulus_);
        for (std::size_t length = 2; length <= Size; length <<= 1) {
            const std::uint64_t root =
                power(prime.primitive_root, (modulus_ - 1) / length, modulus_);
            if (power(root, length, modulus_) != 1 ||
                power(root, length / 2, modulus_) == 1)
                throw std::invalid_argument("invalid primitive NTT root");
            const std::uint64_t inverse_root = invert(root, modulus_);
            const std::size_t half = length >> 1;
            const std::size_t offset = half - 1;
            forward_twiddles_[offset] = ShoupMultiplier(1, modulus_);
            inverse_twiddles_[offset] = ShoupMultiplier(1, modulus_);
            for (std::size_t j = 1; j < half; j++) {
                forward_twiddles_[offset + j] = ShoupMultiplier(
                    multiply(forward_twiddles_[offset + j - 1].value, root,
                             modulus_),
                    modulus_);
                inverse_twiddles_[offset + j] = ShoupMultiplier(
                    multiply(inverse_twiddles_[offset + j - 1].value,
                             inverse_root, modulus_),
                    modulus_);
            }
        }
    }

    void forward(std::array<std::uint64_t, Size> &values) const
    {
        transform<false, false>(values);
    }

    void inverse(std::array<std::uint64_t, Size> &values) const
    {
        transform<true, true>(values);
    }

    void inverseUnscaled(std::array<std::uint64_t, Size> &values) const
    {
        transform<true, false>(values);
    }

private:
    static constexpr std::array<std::size_t, Size> bitReversal()
    {
        std::array<std::size_t, Size> result{};
        for (std::size_t i = 0; i < Size; i++) {
            std::size_t input = i;
            std::size_t reversed = 0;
            for (std::size_t bits = Size; bits > 1; bits >>= 1) {
                reversed = (reversed << 1) | (input & 1U);
                input >>= 1;
            }
            result[i] = reversed;
        }
        return result;
    }

    template <bool Invert, std::size_t Length>
    void applyStages(std::array<std::uint64_t, Size> &values) const
    {
        constexpr std::size_t half = Length >> 1;
        constexpr std::size_t offset = half - 1;
        const auto &twiddles = Invert ? inverse_twiddles_ : forward_twiddles_;
        for (std::size_t block = 0; block < Size; block += Length)
            for (std::size_t j = 0; j < half; j++) {
                const std::uint64_t even = values[block + j];
                const std::uint64_t odd =
                    j == 0 ? values[block + half]
                           : twiddles[offset + j].apply(
                                 values[block + j + half], modulus_);
                values[block + j] = add(even, odd, modulus_);
                values[block + j + half] = subtract(even, odd, modulus_);
            }
        if constexpr (Length < Size) applyStages<Invert, 2 * Length>(values);
    }

    template <bool Invert, bool Normalize>
    void transform(std::array<std::uint64_t, Size> &values) const
    {
        static constexpr auto bit_reversal = bitReversal();
        for (std::size_t i = 1; i < Size; i++)
            if (i < bit_reversal[i])
                std::swap(values[i], values[bit_reversal[i]]);
        applyStages<Invert, 2>(values);
        if constexpr (Normalize)
            for (std::uint64_t &value : values)
                value = multiply(value, inverse_size_, modulus_);
    }

    std::uint64_t modulus_;
    std::uint64_t inverse_size_;
    std::array<ShoupMultiplier, Size - 1> forward_twiddles_{};
    std::array<ShoupMultiplier, Size - 1> inverse_twiddles_{};
};

class NegacyclicNTTPlan {
public:
    NegacyclicNTTPlan(const std::size_t size, const PrimeModulus prime)
        : cyclic_(size, prime), modulus_(prime.value)
    {
        if ((modulus_ - 1) % (2 * size) != 0)
            throw std::invalid_argument(
                "modulus does not support negacyclic NTT size");
        const std::uint64_t psi =
            power(prime.primitive_root, (modulus_ - 1) / (2 * size), modulus_);
        if (power(psi, size, modulus_) != modulus_ - 1)
            throw std::invalid_argument("invalid negacyclic NTT root");

        forward_twist_.resize(size);
        inverse_twist_.resize(size);
        const std::uint64_t inverse_psi = invert(psi, modulus_);
        const std::uint64_t inverse_size = invert(size % modulus_, modulus_);
        forward_twist_[0] = ShoupMultiplier(1, modulus_);
        inverse_twist_[0] = ShoupMultiplier(1, modulus_);
        for (std::size_t i = 1; i < size; i++) {
            forward_twist_[i] = ShoupMultiplier(
                multiply(forward_twist_[i - 1].value, psi, modulus_), modulus_);
            inverse_twist_[i] = ShoupMultiplier(
                multiply(inverse_twist_[i - 1].value, inverse_psi, modulus_),
                modulus_);
        }
        for (ShoupMultiplier &twist : inverse_twist_)
            twist = ShoupMultiplier(
                multiply(twist.value, inverse_size, modulus_), modulus_);
    }

    std::size_t size() const { return cyclic_.size(); }
    std::uint64_t modulus() const { return modulus_; }

    void forward(const std::span<std::uint64_t> values) const
    {
        if (values.size() != size())
            throw std::invalid_argument(
                "negacyclic NTT input has the wrong size");
        for (std::size_t i = 0; i < size(); i++)
            values[i] = forward_twist_[i].apply(values[i], modulus_);
        cyclic_.forward(values);
    }

    void inverse(const std::span<std::uint64_t> values) const
    {
        if (values.size() != size())
            throw std::invalid_argument(
                "negacyclic NTT input has the wrong size");
        cyclic_.inverseUnscaled(values);
        for (std::size_t i = 0; i < size(); i++)
            values[i] = inverse_twist_[i].apply(values[i], modulus_);
    }

private:
    Radix2NTTPlan cyclic_;
    std::uint64_t modulus_;
    std::vector<ShoupMultiplier> forward_twist_;
    std::vector<ShoupMultiplier> inverse_twist_;
};

constexpr std::uint32_t smallPowerMod(std::uint32_t base,
                                      std::uint32_t exponent,
                                      const std::uint32_t modulus)
{
    std::uint64_t result = 1;
    std::uint64_t power_value = base % modulus;
    while (exponent != 0) {
        if ((exponent & 1U) != 0) result = (result * power_value) % modulus;
        power_value = (power_value * power_value) % modulus;
        exponent >>= 1;
    }
    return static_cast<std::uint32_t>(result);
}

constexpr std::uint32_t primitiveRootOfPrime(const std::uint32_t prime)
{
    const std::uint32_t order = prime - 1;
    std::array<std::uint32_t, 32> factors{};
    std::size_t factor_count = 0;
    std::uint32_t remaining = order;
    for (std::uint32_t factor = 2; factor * factor <= remaining; factor++) {
        if (remaining % factor != 0) continue;
        factors[factor_count++] = factor;
        while (remaining % factor == 0) remaining /= factor;
    }
    if (remaining > 1) factors[factor_count++] = remaining;

    for (std::uint32_t candidate = 2; candidate < prime; candidate++) {
        bool primitive = true;
        for (std::size_t i = 0; i < factor_count; i++)
            if (smallPowerMod(candidate, order / factors[i], prime) == 1) {
                primitive = false;
                break;
            }
        if (primitive) return candidate;
    }
    return 0;
}

// Rader's algorithm converts a prime-length DFT into one cyclic convolution
// of length PrimeLength-1.  PrimeLength-1 must be a power of two so the common
// radix-2 plan can execute that convolution.
template <std::size_t PrimeLength>
class RaderNTTPlan {
public:
    static constexpr std::size_t convolution_size = PrimeLength - 1;

    explicit RaderNTTPlan(const PrimeModulus prime)
        : prime_(prime),
          convolution_(prime),
          root_(power(prime.primitive_root, (prime.value - 1) / PrimeLength,
                      prime.value)),
          inverse_prime_(invert(PrimeLength, prime.value))
    {
        static_assert(PrimeLength >= 3);
        static_assert(isPowerOfTwo(convolution_size));
        if ((prime.value - 1) % PrimeLength != 0 ||
            power(root_, PrimeLength, prime.value) != 1 || root_ == 1)
            throw std::invalid_argument(
                "modulus does not support the prime-length NTT");

        constexpr std::uint32_t generator = primitiveRootOfPrime(PrimeLength);
        static_assert(generator != 0);
        powers_[0] = 1;
        for (std::size_t i = 1; i < convolution_size; i++)
            powers_[i] = static_cast<std::uint32_t>(
                (static_cast<std::uint64_t>(powers_[i - 1]) * generator) %
                PrimeLength);

        buildKernel(forward_kernel_, root_, 1);
        buildKernel(inverse_kernel_, invert(root_, prime.value),
                    inverse_prime_);
    }

    std::uint64_t modulus() const { return prime_.value; }
    std::uint64_t root() const { return root_; }

    void forward(std::array<std::uint64_t, PrimeLength> &values) const
    {
        apply(values, forward_kernel_, 1);
    }

    void inverse(std::array<std::uint64_t, PrimeLength> &values) const
    {
        apply(values, inverse_kernel_, inverse_prime_);
    }

private:
    void buildKernel(std::array<ShoupMultiplier, convolution_size> &kernel,
                     const std::uint64_t transform_root,
                     const std::uint64_t output_scale)
    {
        std::array<std::uint64_t, convolution_size> coefficients{};
        for (std::size_t t = 0; t < convolution_size; t++) {
            const std::size_t inverse_index =
                (convolution_size - t) % convolution_size;
            coefficients[t] =
                power(transform_root, powers_[inverse_index], prime_.value);
        }
        convolution_.forward(coefficients);
        const std::uint64_t inverse_convolution_size =
            invert(convolution_size % prime_.value, prime_.value);
        const std::uint64_t kernel_scale =
            multiply(inverse_convolution_size, output_scale, prime_.value);
        for (std::size_t t = 0; t < convolution_size; t++)
            kernel[t] = ShoupMultiplier(
                multiply(coefficients[t], kernel_scale, prime_.value),
                prime_.value);
    }

    void apply(std::array<std::uint64_t, PrimeLength> &values,
               const std::array<ShoupMultiplier, convolution_size> &kernel,
               const std::uint64_t constant_scale) const
    {
        const std::uint64_t zero = values[0];
        std::uint64_t sum = zero;
        std::array<std::uint64_t, convolution_size> permuted{};
        for (std::size_t t = 0; t < convolution_size; t++) {
            permuted[t] = values[powers_[t]];
            sum = add(sum, permuted[t], prime_.value);
        }
        convolution_.forward(permuted);
        for (std::size_t t = 0; t < convolution_size; t++)
            permuted[t] = kernel[t].apply(permuted[t], prime_.value);
        convolution_.inverseUnscaled(permuted);

        const std::uint64_t scaled_zero =
            constant_scale == 1 ? zero
                                : multiply(zero, constant_scale, prime_.value);
        values[0] = constant_scale == 1
                        ? sum
                        : multiply(sum, constant_scale, prime_.value);
        for (std::size_t m = 0; m < convolution_size; m++)
            values[powers_[m]] =
                add(scaled_zero,
                    permuted[(convolution_size - m) % convolution_size],
                    prime_.value);
    }

    PrimeModulus prime_;
    FixedRadix2NTTPlan<convolution_size> convolution_;
    std::uint64_t root_;
    std::uint64_t inverse_prime_;
    std::array<std::uint32_t, convolution_size> powers_{};
    std::array<ShoupMultiplier, convolution_size> forward_kernel_{};
    std::array<ShoupMultiplier, convolution_size> inverse_kernel_{};
};

// Evaluation/interpolation for Z_q[W]/Phi_p(W).  Coefficients are represented
// canonically with degree below p-1; transform values are ordered by the
// numeric nonzero p-th roots omega^1, ..., omega^(p-1).
template <std::size_t PrimeLength>
class PrimeCyclotomicNTTPlan {
public:
    static constexpr std::size_t dimension = PrimeLength - 1;

    explicit PrimeCyclotomicNTTPlan(const PrimeModulus prime)
        : rader_(prime), modulus_(prime.value), root_(rader_.root())
    {
        std::uint64_t root_power = root_;
        for (auto &weight : missing_value_weights_) {
            weight = ShoupMultiplier(root_power, modulus_);
            root_power = multiply(root_power, root_, modulus_);
        }
    }

    std::uint64_t modulus() const { return modulus_; }

    void forward(std::array<std::uint64_t, dimension> &values) const
    {
        std::array<std::uint64_t, PrimeLength> full{};
        for (std::size_t i = 0; i < dimension; i++) full[i] = values[i];
        rader_.forward(full);
        for (std::size_t k = 1; k < PrimeLength; k++) values[k - 1] = full[k];
    }

    void inverse(std::array<std::uint64_t, dimension> &values) const
    {
        std::array<std::uint64_t, PrimeLength> full{};
        std::uint64_t value_at_one = 0;
        for (std::size_t k = 1; k < PrimeLength; k++) {
            full[k] = values[k - 1];
            value_at_one =
                add(value_at_one,
                    missing_value_weights_[k - 1].apply(full[k], modulus_),
                    modulus_);
        }
        // The missing evaluation at W=1 is selected so that the interpolated
        // W^(p-1) coefficient is zero, i.e. the canonical Phi_p representative.
        full[0] = negate(value_at_one, modulus_);
        rader_.inverse(full);
        assert(full[dimension] == 0);
        for (std::size_t i = 0; i < dimension; i++) values[i] = full[i];
    }

private:
    RaderNTTPlan<PrimeLength> rader_;
    std::uint64_t modulus_;
    std::uint64_t root_;
    std::array<ShoupMultiplier, dimension> missing_value_weights_{};
};

class TwoPrimeCRT {
public:
    TwoPrimeCRT(const PrimeModulus first, const PrimeModulus second)
        : first_(first.value),
          second_(second.value),
          product_(static_cast<unsigned __int128>(first_) * second_),
          first_inverse_mod_second_(invert(first_ % second_, second_))
    {
    }

    unsigned __int128 modulusProduct() const { return product_; }

    __int128 reconstructSigned(const std::uint64_t first_residue,
                               const std::uint64_t second_residue) const
    {
        const std::uint64_t first_mod_second =
            first_ < second_ ? first_residue
            : first_ - second_ < second_
                ? (first_residue >= second_ ? first_residue - second_
                                            : first_residue)
                : first_residue % second_;
        const std::uint64_t delta =
            subtract(second_residue, first_mod_second, second_);
        const std::uint64_t quotient =
            multiply(delta, first_inverse_mod_second_, second_);
        const unsigned __int128 value =
            first_residue + static_cast<unsigned __int128>(first_) * quotient;
        if (value <= product_ / 2) return static_cast<__int128>(value);
        return -static_cast<__int128>(product_ - value);
    }

private:
    std::uint64_t first_;
    std::uint64_t second_;
    unsigned __int128 product_;
    std::uint64_t first_inverse_mod_second_;
};

// Reconstruct values that are only slightly wider than the first two-prime
// product without materializing the roughly 186-bit three-prime product.  A
// caller-proved |value| < 2^126 can differ from the centered two-prime lift by
// at most four multiples of q0*q1; rejecting a larger quotient keeps every
// intermediate within signed __int128.
class ThreePrimeCRT {
public:
    ThreePrimeCRT(const PrimeModulus first, const PrimeModulus second,
                  const PrimeModulus third)
        : first_two_(first, second),
          third_(third.value),
          first_two_product_(first_two_.modulusProduct()),
          product_inverse_mod_third_(invert(
              static_cast<std::uint64_t>(first_two_product_ % third_), third_))
    {
    }

    __int128 reconstructSignedBounded(const std::uint64_t first_residue,
                                      const std::uint64_t second_residue,
                                      const std::uint64_t third_residue) const
    {
        const __int128 lower =
            first_two_.reconstructSigned(first_residue, second_residue);
        const std::uint64_t lower_residue = signedResidue(lower, third_);
        const std::uint64_t delta =
            subtract(third_residue, lower_residue, third_);
        const std::uint64_t quotient_residue =
            multiply(delta, product_inverse_mod_third_, third_);
        const __int128 quotient =
            quotient_residue <= third_ / 2
                ? static_cast<__int128>(quotient_residue)
                : -static_cast<__int128>(third_ - quotient_residue);
        if (quotient < -4 || quotient > 4)
            throw std::overflow_error(
                "three-prime CRT value exceeds the signed 126-bit bound");
        return lower + static_cast<__int128>(first_two_product_) * quotient;
    }

private:
    static std::uint64_t signedResidue(const __int128 value,
                                       const std::uint64_t modulus)
    {
        if (value >= 0)
            return static_cast<std::uint64_t>(
                static_cast<unsigned __int128>(value) % modulus);
        const unsigned __int128 magnitude =
            static_cast<unsigned __int128>(-(value + 1)) + 1;
        const std::uint64_t residue =
            static_cast<std::uint64_t>(magnitude % modulus);
        return residue == 0 ? 0 : modulus - residue;
    }

    TwoPrimeCRT first_two_;
    std::uint64_t third_;
    unsigned __int128 first_two_product_;
    std::uint64_t product_inverse_mod_third_;
};

inline __int128 centeredResidue(const std::uint64_t residue,
                                const std::uint64_t modulus)
{
    if (residue <= modulus / 2) return residue;
    return -static_cast<__int128>(modulus - residue);
}

}  // namespace TFHEpp::modular_ntt
