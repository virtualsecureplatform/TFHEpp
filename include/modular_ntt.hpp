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

// These primes support power-of-two transforms through length 2^11 or
// greater and a length-17 transform.  They are deliberately scheme-neutral:
// callers may use either one prime or the pair according to their exact
// coefficient bound.
inline constexpr std::array<PrimeModulus, 2> wide_primes{{
    {4611686018426953729ULL, 11},
    {4611686018426884097ULL, 5},
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

inline std::uint64_t multiply(const std::uint64_t lhs, const std::uint64_t rhs,
                              const std::uint64_t modulus)
{
    return static_cast<std::uint64_t>(
        (static_cast<unsigned __int128>(lhs) * rhs) % modulus);
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
        : size_(size), modulus_(prime.value), inverse_size_(0)
    {
        if (!isPowerOfTwo(size_) || size_ < 2)
            throw std::invalid_argument("NTT size must be a power of two");
        if (modulus_ >= (std::uint64_t{1} << 63) || (modulus_ - 1) % size_ != 0)
            throw std::invalid_argument("modulus does not support NTT size");

        inverse_size_ = invert(size_ % modulus_, modulus_);
        for (std::size_t length = 2; length <= size_; length <<= 1) {
            const std::uint64_t root =
                power(prime.primitive_root, (modulus_ - 1) / length, modulus_);
            if (power(root, length, modulus_) != 1 ||
                power(root, length / 2, modulus_) == 1)
                throw std::invalid_argument("invalid primitive NTT root");
            forward_stage_roots_.push_back(root);
            inverse_stage_roots_.push_back(invert(root, modulus_));
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
        transform(values, true);
    }

private:
    void transform(const std::span<std::uint64_t> values,
                   const bool invert) const
    {
        if (values.size() != size_)
            throw std::invalid_argument("NTT input has the wrong size");

        for (std::size_t i = 1, j = 0; i < size_; i++) {
            std::size_t bit = size_ >> 1;
            for (; (j & bit) != 0; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap(values[i], values[j]);
        }

        std::size_t stage = 0;
        for (std::size_t length = 2; length <= size_; length <<= 1, stage++) {
            const std::uint64_t stage_root = invert
                                                 ? inverse_stage_roots_[stage]
                                                 : forward_stage_roots_[stage];
            const std::size_t half = length >> 1;
            for (std::size_t block = 0; block < size_; block += length) {
                std::uint64_t twiddle = 1;
                for (std::size_t j = 0; j < half; j++) {
                    const std::uint64_t even = values[block + j];
                    const std::uint64_t odd =
                        multiply(values[block + j + half], twiddle, modulus_);
                    values[block + j] = add(even, odd, modulus_);
                    values[block + j + half] = subtract(even, odd, modulus_);
                    twiddle = multiply(twiddle, stage_root, modulus_);
                }
            }
        }

        if (invert)
            for (std::uint64_t &value : values)
                value = multiply(value, inverse_size_, modulus_);
    }

    std::size_t size_;
    std::uint64_t modulus_;
    std::uint64_t inverse_size_;
    std::vector<std::uint64_t> forward_stage_roots_;
    std::vector<std::uint64_t> inverse_stage_roots_;
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
        forward_twist_[0] = inverse_twist_[0] = 1;
        for (std::size_t i = 1; i < size; i++) {
            forward_twist_[i] = multiply(forward_twist_[i - 1], psi, modulus_);
            inverse_twist_[i] =
                multiply(inverse_twist_[i - 1], inverse_psi, modulus_);
        }
    }

    std::size_t size() const { return cyclic_.size(); }
    std::uint64_t modulus() const { return modulus_; }

    void forward(const std::span<std::uint64_t> values) const
    {
        if (values.size() != size())
            throw std::invalid_argument(
                "negacyclic NTT input has the wrong size");
        for (std::size_t i = 0; i < size(); i++)
            values[i] = multiply(values[i], forward_twist_[i], modulus_);
        cyclic_.forward(values);
    }

    void inverse(const std::span<std::uint64_t> values) const
    {
        if (values.size() != size())
            throw std::invalid_argument(
                "negacyclic NTT input has the wrong size");
        cyclic_.inverse(values);
        for (std::size_t i = 0; i < size(); i++)
            values[i] = multiply(values[i], inverse_twist_[i], modulus_);
    }

private:
    Radix2NTTPlan cyclic_;
    std::uint64_t modulus_;
    std::vector<std::uint64_t> forward_twist_;
    std::vector<std::uint64_t> inverse_twist_;
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
          convolution_(convolution_size, prime),
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

        buildKernel(forward_kernel_, root_);
        buildKernel(inverse_kernel_, invert(root_, prime.value));
    }

    std::uint64_t modulus() const { return prime_.value; }
    std::uint64_t root() const { return root_; }

    void forward(std::array<std::uint64_t, PrimeLength> &values) const
    {
        apply(values, forward_kernel_);
    }

    void inverse(std::array<std::uint64_t, PrimeLength> &values) const
    {
        apply(values, inverse_kernel_);
        for (std::uint64_t &value : values)
            value = multiply(value, inverse_prime_, prime_.value);
    }

private:
    void buildKernel(std::array<std::uint64_t, convolution_size> &kernel,
                     const std::uint64_t transform_root)
    {
        for (std::size_t t = 0; t < convolution_size; t++) {
            const std::size_t inverse_index =
                (convolution_size - t) % convolution_size;
            kernel[t] =
                power(transform_root, powers_[inverse_index], prime_.value);
        }
        convolution_.forward(kernel);
    }

    void apply(std::array<std::uint64_t, PrimeLength> &values,
               const std::array<std::uint64_t, convolution_size> &kernel) const
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
            permuted[t] = multiply(permuted[t], kernel[t], prime_.value);
        convolution_.inverse(permuted);

        values[0] = sum;
        for (std::size_t m = 0; m < convolution_size; m++)
            values[powers_[m]] =
                add(zero, permuted[(convolution_size - m) % convolution_size],
                    prime_.value);
    }

    PrimeModulus prime_;
    Radix2NTTPlan convolution_;
    std::uint64_t root_;
    std::uint64_t inverse_prime_;
    std::array<std::uint32_t, convolution_size> powers_{};
    std::array<std::uint64_t, convolution_size> forward_kernel_{};
    std::array<std::uint64_t, convolution_size> inverse_kernel_{};
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
        std::uint64_t root_power = root_;
        std::uint64_t value_at_one = 0;
        for (std::size_t k = 1; k < PrimeLength; k++) {
            full[k] = values[k - 1];
            value_at_one =
                add(value_at_one, multiply(full[k], root_power, modulus_),
                    modulus_);
            root_power = multiply(root_power, root_, modulus_);
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
        const std::uint64_t delta =
            subtract(second_residue, first_residue % second_, second_);
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

inline __int128 centeredResidue(const std::uint64_t residue,
                                const std::uint64_t modulus)
{
    if (residue <= modulus / 2) return residue;
    return -static_cast<__int128>(modulus - residue);
}

}  // namespace TFHEpp::modular_ntt
