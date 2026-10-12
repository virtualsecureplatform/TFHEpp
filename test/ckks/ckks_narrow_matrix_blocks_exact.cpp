#include <ckks/narrow_matrix.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

template <std::size_t N, std::size_t Tile>
void check(std::size_t count, std::size_t columns, std::size_t rows, std::size_t repetitions)
{
    using namespace TFHEpp;
    std::mt19937_64 random(20262119 + N + columns + count);
    for (std::size_t prime = 0; prime < 4; ++prime) {
        const auto modulus = prime < 2 ? modular_ntt::degree65536_primes[prime].value
                                      : prime == 2 ? modular_ntt::wide_primes[0].value
                                                   : (std::uint64_t{1} << 61) - 1;
        const narrow_detail::SelectorMultiplierFactory factory(modulus);
        std::vector<std::uint64_t> selectors(rows * count * N), tables(2 * N * columns);
        std::vector<std::size_t> exponents(count);
        for (auto &value : selectors) value = random() % modulus;
        for (auto &value : tables) value = random() % modulus;
        // Exercise maximum canonical residues and repeated additions
        // near the lazy range boundary, in addition to random coefficients.
        for (std::size_t i = 0; i < selectors.size(); i += N)
            selectors[i] = modulus - 1;
        std::fill_n(tables.begin(), columns, modulus - 1);
        std::fill_n(tables.begin() + N * columns, columns, modulus - 1);
        for (std::size_t y = 0; y < count; ++y)
            exponents[y] = modular_ntt::power(5, (N / 2 - y % (N / 2)) % (N / 2), 2 * N);
        std::vector<std::uint64_t> reference;
        for (const bool blocked : {false, true, true, false}) {
            std::vector<std::uint64_t> sums(2 * rows * N * columns, modulus - 1);
            const auto start = std::chrono::steady_clock::now();
            for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
#pragma omp parallel for schedule(static)
                for (std::size_t row = 0; row < rows; ++row) {
                    const auto get_row = [&](std::size_t local) {
                        return narrow_detail::SelectorMatrixRow{
                            selectors.data() + (row * count + local) * N, exponents[local]};
                    };
                    auto *first = sums.data() + row * N * columns;
                    auto *second = sums.data() + (rows + row) * N * columns;
                    if (blocked)
                        narrow_detail::addSelectorMatrixRowsPair<N, Tile>(first, second,
                            tables.data(), tables.data() + N * columns,
                            columns, count, get_row, factory, modulus);
                    else
                        narrow_detail::addSelectorMatrixRowsPair<N, N, false>(first, second,
                            tables.data(), tables.data() + N * columns,
                            columns, count, get_row, factory, modulus);
                }
            }
            const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (reference.empty()) reference = sums;
            else if (reference != sums) throw std::runtime_error("matrix blocking changed modular sums");
            if constexpr (N == 512) {
                for (std::size_t row = 0; row < rows; ++row)
                    for (std::size_t i = 0; i < N; ++i)
                        for (std::size_t column = 0; column < columns; ++column)
                            for (std::size_t table = 0; table < 2; ++table) {
                                std::uint64_t expected = modulus - 1;
                                for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
                                    for (std::size_t y = 0; y < count; ++y) {
                                        const auto point = ((exponents[y] * (2 * i + 1)) % (2 * N) - 1) / 2;
                                        expected = (static_cast<unsigned __int128>(selectors[(row * count + y) * N + i]) *
                                            tables[(table * N + point) * columns + column] + expected) % modulus;
                                    }
                                if (sums[((table * rows + row) * N + i) * columns + column] != expected)
                                    throw std::runtime_error("matrix blocks differ from independent modular oracle");
                            }
            }
            std::cout << "{\"n\":" << N << ",\"prime\":" << prime
                      << ",\"columns\":" << columns << ",\"selector_rows\":" << count
                      << ",\"tile\":" << (blocked ? Tile : N) << ",\"seconds\":" << seconds
                      << ",\"exact_coefficients\":" << sums.size() << ",\"passed\":true}" << std::endl;
        }
    }
}

int main()
{
    check<512, 37>(3, 5, 3, 2); // Partial final block and scalar SIMD tail.
    check<32768, 256>(4, 8, 22, 8);
    check<32768, 256>(40, 18, 22, 2);
    check<32768, 256>(4, 22, 2, 2); // Combined-root SIMD tail.
}
