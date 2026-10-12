#include <ckks/narrow_matrix.hpp>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

template <std::size_t N>
std::size_t check()
{
    using namespace TFHEpp;
    std::mt19937_64 random(20262119 + N);
    std::size_t checked = 0;
    for (const auto &prime : modular_ntt::degree65536_primes) {
        const modular_ntt::NegacyclicNTTPlan plan(N, prime);
        std::vector<std::uint64_t> polynomial(N);
        for (auto &coefficient : polynomial) coefficient = random() % prime.value;
        auto base = polynomial;
        plan.forward(base);
        // Include wraparound rotations and both signs of the automorphism.
        for (const std::size_t rotation : {0UL, 1UL, N / 31, N / 2 - 1, N / 2, N - 1}) {
            const auto positive = modular_ntt::power(5, rotation, 2 * N);
            for (const auto exponent : {positive, 2 * N - positive}) {
                std::vector<std::uint64_t> rotated(N);
                // Independent coefficient-domain substitution X -> X^exponent
                // in Z_q[X]/(X^N+1), followed by a fresh transform.
                for (std::size_t j = 0; j < N; ++j) {
                    const auto degree = (j * exponent) % (2 * N);
                    rotated[degree % N] = degree < N || polynomial[j] == 0
                        ? polynomial[j] : prime.value - polynomial[j];
                }
                plan.forward(rotated);
                for (std::size_t i = 0; i < N; ++i) {
                    if (rotated[i] != base[narrow_detail::narrowNTTPermutationIndex<N>(exponent, i)])
                        throw std::runtime_error("implicit NTT permutation differs from coefficient automorphism");
                    ++checked;
                }
            }
        }
    }
    return checked;
}

int main()
{
    const auto checked = check<512>() + check<8192>() + check<16384>() + check<32768>();
    std::cout << "{\"independent_ntt_points_checked\":" << checked << ",\"passed\":true}\n";
}
