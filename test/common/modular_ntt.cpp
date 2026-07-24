#include <array>
#include <cstdint>
#include <iostream>
#include <modular_ntt.hpp>
#include <random>
#include <vector>

namespace {

using namespace TFHEpp::modular_ntt;

std::uint64_t signedResidue(const __int128 value, const std::uint64_t modulus)
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

bool checkMultiply(const PrimeModulus prime)
{
    const auto reference = [&](const std::uint64_t lhs,
                               const std::uint64_t rhs) {
        return static_cast<std::uint64_t>(static_cast<unsigned __int128>(lhs) *
                                          rhs % prime.value);
    };
    const std::array<std::uint64_t, 5> edges{0, 1, prime.value / 2,
                                             prime.value - 2, prime.value - 1};
    for (const auto lhs : edges)
        for (const auto rhs : edges)
            if (multiply(lhs, rhs, prime.value) != reference(lhs, rhs))
                return false;

    std::mt19937_64 rng(0x4d554c5449504c59ULL ^ prime.value);
    for (std::size_t sample = 0; sample < 100000; sample++) {
        const std::uint64_t lhs = rng() % prime.value;
        const std::uint64_t rhs = rng() % prime.value;
        if (multiply(lhs, rhs, prime.value) != reference(lhs, rhs))
            return false;
    }
    return true;
}

bool checkShoup(const PrimeModulus prime)
{
    const std::array<std::uint64_t, 5> edges{0, 1, prime.value / 2,
                                             prime.value - 2, prime.value - 1};
    for (const auto constant : edges) {
        const ShoupMultiplier multiplier(constant, prime.value);
        for (const auto input : edges)
            if (multiplier.apply(input, prime.value) !=
                multiply(input, constant, prime.value))
                return false;
    }

    std::mt19937_64 rng(0x53484f55504d554cULL ^ prime.value);
    for (std::size_t sample = 0; sample < 100000; sample++) {
        const std::uint64_t constant = rng() % prime.value;
        const std::uint64_t input = rng() % prime.value;
        const ShoupMultiplier multiplier(constant, prime.value);
        if (multiplier.apply(input, prime.value) !=
            multiply(input, constant, prime.value))
            return false;
        const std::uint64_t lazy_input = rng() % (2 * prime.value);
        const std::uint64_t lazy =
            multiplier.applyLazy(lazy_input, prime.value);
        if (lazy >= 2 * prime.value ||
            (lazy >= prime.value ? lazy - prime.value : lazy) !=
                multiply(lazy_input % prime.value, constant, prime.value))
            return false;
    }
    return true;
}

bool checkRadix2(const PrimeModulus prime)
{
    constexpr std::size_t size = 1024;
    Radix2NTTPlan plan(size, prime);
    std::mt19937_64 rng(0x4e5454);
    std::vector<std::uint64_t> values(size);
    for (auto &value : values) value = rng() % prime.value;
    const auto expected = values;
    plan.forward(values);
    plan.inverse(values);
    return values == expected;
}

bool checkNegacyclic(const PrimeModulus prime)
{
    constexpr std::size_t size = 16;
    NegacyclicNTTPlan plan(size, prime);
    std::array<std::uint64_t, size> lhs{};
    std::array<std::uint64_t, size> rhs{};
    std::array<std::uint64_t, size> expected{};
    for (std::size_t i = 0; i < size; i++) {
        lhs[i] = (7 * i + 3) % prime.value;
        rhs[i] = (11 * i + 5) % prime.value;
    }
    const auto lhs_coefficients = lhs;
    for (std::size_t i = 0; i < size; i++)
        for (std::size_t j = 0; j < size; j++) {
            const std::uint64_t product = multiply(lhs[i], rhs[j], prime.value);
            const std::size_t degree = i + j;
            const std::size_t output = degree % size;
            expected[output] =
                degree < size
                    ? add(expected[output], product, prime.value)
                    : subtract(expected[output], product, prime.value);
        }

    plan.forward(lhs);
    const std::uint64_t psi = power(
        prime.primitive_root, (prime.value - 1) / (2 * size), prime.value);
    for (std::size_t point = 0; point < size; point++) {
        const std::uint64_t root =
            power(psi, 2 * point + 1, prime.value);
        std::uint64_t expected_value = 0;
        std::uint64_t root_power = 1;
        for (std::size_t coefficient = 0; coefficient < size; coefficient++) {
            expected_value = add(
                expected_value,
                multiply(lhs_coefficients[coefficient], root_power,
                         prime.value),
                prime.value);
            root_power = multiply(root_power, root, prime.value);
        }
        if (lhs[point] != expected_value) return false;
    }
    plan.forward(rhs);
    for (std::size_t i = 0; i < size; i++)
        lhs[i] = multiply(lhs[i], rhs[i], prime.value);
    plan.inverse(lhs);
    return lhs == expected;
}

bool checkRader(const PrimeModulus prime)
{
    constexpr std::size_t p = 17;
    RaderNTTPlan<p> plan(prime);
    std::array<std::uint64_t, p> values{};
    for (std::size_t i = 0; i < p; i++) values[i] = 13 * i + 9;
    const auto input = values;
    std::array<std::uint64_t, p> expected{};
    for (std::size_t k = 0; k < p; k++)
        for (std::size_t j = 0; j < p; j++)
            expected[k] =
                add(expected[k],
                    multiply(input[j], power(plan.root(), j * k, prime.value),
                             prime.value),
                    prime.value);

    plan.forward(values);
    if (values != expected) return false;
    plan.inverse(values);
    return values == input;
}

bool checkCyclotomic(const PrimeModulus prime)
{
    constexpr std::size_t p = 17;
    constexpr std::size_t dimension = p - 1;
    PrimeCyclotomicNTTPlan<p> plan(prime);
    std::array<std::uint64_t, dimension> lhs{};
    std::array<std::uint64_t, dimension> rhs{};
    for (std::size_t i = 0; i < dimension; i++) {
        lhs[i] = 5 * i + 1;
        rhs[i] = 9 * i + 2;
    }
    const auto original = lhs;
    plan.forward(lhs);
    plan.inverse(lhs);
    if (lhs != original) return false;

    lhs = original;
    std::array<std::uint64_t, 2 * dimension - 1> unreduced{};
    for (std::size_t i = 0; i < dimension; i++)
        for (std::size_t j = 0; j < dimension; j++)
            unreduced[i + j] =
                add(unreduced[i + j], multiply(lhs[i], rhs[j], prime.value),
                    prime.value);
    std::array<std::uint64_t, dimension> expected{};
    for (std::size_t i = 0; i < dimension; i++) expected[i] = unreduced[i];
    for (std::size_t degree = p; degree < unreduced.size(); degree++)
        expected[degree - p] =
            add(expected[degree - p], unreduced[degree], prime.value);
    for (std::size_t i = 0; i < dimension; i++)
        expected[i] = subtract(expected[i], unreduced[dimension], prime.value);

    plan.forward(lhs);
    plan.forward(rhs);
    for (std::size_t i = 0; i < dimension; i++)
        lhs[i] = multiply(lhs[i], rhs[i], prime.value);
    plan.inverse(lhs);
    return lhs == expected;
}

bool checkCRT()
{
    const TwoPrimeCRT crt(wide_primes[0], wide_primes[1]);
    const std::array<__int128, 5> values{
        0,
        123456789,
        -987654321,
        (static_cast<__int128>(1) << 110) + 7654321,
        -((static_cast<__int128>(1) << 113) + 1234567),
    };
    for (const __int128 value : values)
        if (crt.reconstructSigned(signedResidue(value, wide_primes[0].value),
                                  signedResidue(value, wide_primes[1].value)) !=
            value)
            return false;

    std::mt19937_64 rng(0x43525452414e444fULL);
    for (std::size_t sample = 0; sample < 100000; sample++) {
        const std::uint64_t first = rng() % wide_primes[0].value;
        const std::uint64_t second = rng() % wide_primes[1].value;
        const __int128 reconstructed = crt.reconstructSigned(first, second);
        if (signedResidue(reconstructed, wide_primes[0].value) != first ||
            signedResidue(reconstructed, wide_primes[1].value) != second)
            return false;
    }

    const TwoPrimeCRT generic({17, 3}, {5, 2});
    const TwoPrimeCRT ascending({5, 2}, {17, 3});
    for (__int128 value = -42; value <= 42; value++) {
        const std::uint64_t mod17 = signedResidue(value, 17);
        const std::uint64_t mod5 = signedResidue(value, 5);
        if (generic.reconstructSigned(mod17, mod5) != value ||
            ascending.reconstructSigned(mod5, mod17) != value)
            return false;
    }
    return true;
}

bool checkThreePrimeCRT()
{
    const ThreePrimeCRT crt(wide_primes[0], wide_primes[1], wide_primes[2]);
    const std::array<__int128, 5> values{
        0,
        123456789,
        -987654321,
        (static_cast<__int128>(1) << 125) + 7654321,
        -((static_cast<__int128>(1) << 125) + 1234567),
    };
    for (const __int128 value : values)
        if (crt.reconstructSignedBounded(
                signedResidue(value, wide_primes[0].value),
                signedResidue(value, wide_primes[1].value),
                signedResidue(value, wide_primes[2].value)) != value)
            return false;
    return true;
}

}  // namespace

int main()
{
    for (const auto prime : TFHEpp::modular_ntt::wide_primes) {
        if (!checkMultiply(prime) || !checkShoup(prime) ||
            !checkRadix2(prime) || !checkNegacyclic(prime) ||
            !checkRader(prime) || !checkCyclotomic(prime)) {
            std::cerr << "modular NTT regression failed for " << prime.value
                      << std::endl;
            return 1;
        }
    }
    if (!checkCRT() || !checkThreePrimeCRT()) {
        std::cerr << "CRT regression failed" << std::endl;
        return 1;
    }
    std::cout << "modular NTT regression passed" << std::endl;
    return 0;
}
