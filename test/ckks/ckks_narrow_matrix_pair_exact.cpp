#include <ckks/narrow_matrix.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using namespace TFHEpp;

int main()
{
    std::mt19937_64 random(20262119);
    std::size_t checked = 0;
    for (const auto &prime : {modular_ntt::degree65536_primes[0], modular_ntt::degree65536_primes[1],
                             modular_ntt::PrimeModulus{(std::uint64_t{1} << 62) - 1, 1},
                             modular_ntt::PrimeModulus{(std::uint64_t{1} << 63) - 1, 1}}) {
        const auto modulus = prime.value;
        for (const std::size_t columns : {1, 3, 4, 5, 8, 12, 18, 37}) {
            std::vector<std::uint64_t> first(columns), second(columns), left(columns), right(columns);
            auto expected_first = first, expected_second = second;
            for (std::size_t point = 0; point < 16384; point++) {
                for (std::size_t column = 0; column < columns; column++) {
                    first[column] = random() % modulus;
                    second[column] = random() % modulus;
                }
                expected_first = first;
                expected_second = second;
                for (std::size_t row = 0; row < 17; row++) {
                    const auto selector = row == 0 ? 0 : row == 1 ? modulus-1 : random() % modulus;
                    for (std::size_t column = 0; column < columns; column++) {
                        left[column] = row < 2 ? modulus-1 : random() % modulus;
                        right[column] = row < 2 ? 0 : random() % modulus;
                        expected_first[column] = (static_cast<unsigned __int128>(selector)*left[column]+expected_first[column]) % modulus;
                        expected_second[column] = (static_cast<unsigned __int128>(selector)*right[column]+expected_second[column]) % modulus;
                    }
                    narrow_detail::addSelectorTablePairVector(first.data(), second.data(), left.data(), right.data(), columns, selector, modulus);
                    if (first != expected_first || second != expected_second)
                        throw std::runtime_error("independent paired modular accumulation mismatch");
                    checked += 2*columns;
                }
            }
        }
    }
    std::cout << "{\"independent_modular_coefficients_checked\":" << checked << ",\"passed\":true}" << std::endl;

    // Match the full-ring 48-bit selector matrix: 8 digits, 2 ciphertext
    // components, 32 selector rows, 12 frequency columns and 2 signed tables.
    constexpr std::size_t n = 16384, rows = 16, count = 32, columns = 12;
    const auto modulus = modular_ntt::degree65536_primes[0].value;
    std::vector<std::uint64_t> selectors(rows*count*n), tables(2*n*columns);
    for (auto &value : selectors) value = random()%modulus;
    for (auto &value : tables) value = random()%modulus;
    std::vector<std::size_t> permutations(count*n);
    for (std::size_t y = 0; y < count; y++) {
        const auto exponent = modular_ntt::power(5, (n/2-y)%(n/2), 2*n);
        for (std::size_t i = 0; i < n; i++)
            permutations[y*n+i] = ((exponent*(2*i+1))%(2*n)-1)/2;
    }
    std::vector<std::uint64_t> expected;
    for (const bool vector : {false, true, true, false, false, true, true, false}) {
        std::vector<std::uint64_t> sums(2*rows*n*columns);
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t repetition = 0; repetition < 8; ++repetition) {
#pragma omp parallel for schedule(static)
            for (std::size_t row = 0; row < rows; row++)
                for (std::size_t y = 0; y < count; y++)
                    for (std::size_t i = 0; i < n; i++) {
                        const auto index = (row*n+i)*columns;
                        const auto table_index = permutations[y*n+i]*columns;
                        if (vector)
                            narrow_detail::addSelectorTablePairVector(sums.data()+index, sums.data()+rows*n*columns+index,
                                tables.data()+table_index, tables.data()+n*columns+table_index,
                                columns, selectors[(row*count+y)*n+i], modulus);
                        else
                            narrow_detail::addSelectorTablePair(sums.data()+index, sums.data()+rows*n*columns+index,
                                tables.data()+table_index, tables.data()+n*columns+table_index,
                                columns, selectors[(row*count+y)*n+i], modulus);
                    }
        }
        const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if (expected.empty()) expected = sums;
        else if (expected != sums) throw std::runtime_error("paired full matrix differs from original");
        std::cout << "{\"vector\":" << (vector ? "true" : "false") << ",\"seconds\":" << seconds
                  << ",\"exact_coefficients\":" << sums.size() << "}" << std::endl;
    }
}
