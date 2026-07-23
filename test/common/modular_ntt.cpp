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
    return true;
}

}  // namespace

int main()
{
    for (const auto prime : TFHEpp::modular_ntt::wide_primes) {
        if (!checkRadix2(prime) || !checkNegacyclic(prime) ||
            !checkRader(prime) || !checkCyclotomic(prime)) {
            std::cerr << "modular NTT regression failed for " << prime.value
                      << std::endl;
            return 1;
        }
    }
    if (!checkCRT()) {
        std::cerr << "CRT regression failed" << std::endl;
        return 1;
    }
    std::cout << "modular NTT regression passed" << std::endl;
    return 0;
}
