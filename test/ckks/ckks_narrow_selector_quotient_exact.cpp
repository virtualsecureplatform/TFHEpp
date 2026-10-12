#include <ckks/narrow_matrix.hpp>
#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

int main()
{
    using namespace TFHEpp;
    std::mt19937_64 random(20262119);
    std::vector<std::uint64_t> moduli{1, 2, 3, 16, 257,
        (std::uint64_t{1} << 62) - 1, (std::uint64_t{1} << 63) - 1};
    for (const auto &prime : modular_ntt::degree65536_primes) moduli.push_back(prime.value);
    std::size_t checked = 0, products = 0;
    for (const auto modulus : moduli) {
        const narrow_detail::SelectorMultiplierFactory factory(modulus);
        for (std::size_t sample = 0; sample < 100006; ++sample) {
            const auto selector = sample < 6
                ? std::array<std::uint64_t, 6>{0, modulus - 1, modulus / 2,
                    (modulus - 1) / 2, modulus / 4, modulus - modulus / 4 - 1}[sample]
                : random() % modulus;
            const auto fixed = factory.make(selector);
            const auto expected = (static_cast<unsigned __int128>(selector) << 64) / modulus;
            if (fixed.value != selector || fixed.quotient != expected)
                throw std::runtime_error("selector reciprocal quotient differs from exact division");
            ++checked;
            for (const std::size_t columns : {1, 3, 4, 5, 8, 18}) {
                std::vector<std::uint64_t> first(columns), second(columns), left(columns), right(columns);
                for (std::size_t i = 0; i < columns; ++i) {
                    first[i] = random() % modulus; second[i] = random() % modulus;
                    left[i] = random() % modulus; right[i] = random() % modulus;
                }
                auto a = first, b = second;
                narrow_detail::addSelectorTablePairPrepared(first.data(), second.data(), left.data(), right.data(), columns, fixed, modulus);
                for (std::size_t i = 0; i < columns; ++i) {
                    if (first[i] != (static_cast<unsigned __int128>(selector) * left[i] + a[i]) % modulus ||
                        second[i] != (static_cast<unsigned __int128>(selector) * right[i] + b[i]) % modulus)
                        throw std::runtime_error("prepared selector accumulation differs from exact modular arithmetic");
                    products += 2;
                }
            }
        }
        try { factory.make(modulus); throw std::runtime_error("noncanonical selector accepted"); }
        catch (const std::invalid_argument &) {}
    }
    for (const auto modulus : {0ULL, 1ULL << 63}) {
        try { narrow_detail::SelectorMultiplierFactory invalid(modulus); throw std::runtime_error("invalid modulus accepted"); }
        catch (const std::invalid_argument &) {}
    }
    std::cout << "{\"exact_quotients_checked\":" << checked
              << ",\"independent_modular_products_checked\":" << products << ",\"passed\":true}\n";

    // Paired matrix trials compare the original division-based vector path
    // with the factory, using the same selectors, tables and permutations.
    constexpr std::size_t n = 16384, rows = 8, count = 16, columns = 8;
    // Keep the modulus runtime-selected as in the prime loop in the kernel;
    // a compile-time divisor lets the compiler precompute its own reciprocal.
    volatile std::uint64_t runtime_modulus = modular_ntt::degree65536_primes[0].value;
    const auto modulus = runtime_modulus;
    const narrow_detail::SelectorMultiplierFactory factory(modulus);
    std::vector<std::uint64_t> selectors(rows * count * n), tables(2 * n * columns);
    for (auto &value : selectors) value = random() % modulus;
    for (auto &value : tables) value = random() % modulus;
    std::vector<std::uint64_t> reference;
    for (const bool reciprocal : {false, true, true, false, false, true}) {
        std::vector<std::uint64_t> sums(2 * rows * n * columns);
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t repetition = 0; repetition < 2; ++repetition) {
#pragma omp parallel for schedule(static)
            for (std::size_t row = 0; row < rows; ++row)
                for (std::size_t y = 0; y < count; ++y) {
                    const auto exponent = modular_ntt::power(5, (n / 2 - y) % (n / 2), 2 * n);
                    for (std::size_t i = 0; i < n; ++i) {
                        const auto destination = (row * n + i) * columns;
                        const auto point = narrow_detail::narrowNTTPermutationIndex<n>(exponent, i) * columns;
                        const auto selector = selectors[(row * count + y) * n + i];
                        if (reciprocal)
                            narrow_detail::addSelectorTablePairPrepared(sums.data() + destination,
                                sums.data() + rows * n * columns + destination,
                                tables.data() + point, tables.data() + n * columns + point,
                                columns, factory.make(selector), modulus);
                        else
                            narrow_detail::addSelectorTablePairVector(sums.data() + destination,
                                sums.data() + rows * n * columns + destination,
                                tables.data() + point, tables.data() + n * columns + point,
                                columns, selector, modulus);
                    }
                }
        }
        const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (reference.empty()) reference = sums;
        else if (sums != reference) throw std::runtime_error("reciprocal matrix changed exact sums");
        std::cout << "{\"reciprocal\":" << (reciprocal ? "true" : "false")
                  << ",\"seconds\":" << seconds << ",\"exact_coefficients\":" << sums.size() << "}\n";
    }
}
