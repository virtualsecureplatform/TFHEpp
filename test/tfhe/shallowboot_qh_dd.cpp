#include <array>
#include <cassert>
#include <cstdint>
#include <random>
#include <vector>

#include <bfv/bfv++.hpp>
#include <bfv/bfv-slots.hpp>

struct ShallowDDParam : TFHEpp::lvl3simdparam {
    static const inline double α = std::pow(2.0, -105);
    static constexpr std::uint64_t simd_modulus = 114689;
    static constexpr std::uint64_t simd_psi = 80720;
    static constexpr std::uint64_t simd_psi_inv = 7887;
    static constexpr std::uint64_t simd_n_inv = 114661;
};


int main()
{
    using namespace TFHEpp;
    using P = ShallowDDParam;
    std::mt19937_64 engine(20261730);
    std::uniform_int_distribution<std::uint64_t> word;

    Key<P> key;
    Polynomial<P> secret;
    for (std::size_t i = 0; i < P::n; i++) {
        const bool positive = (engine() & 1U) != 0;
        key[i] = positive ? typename P::T{1} : static_cast<typename P::T>(-1);
        secret[i] = key[i];
    }

    Polynomial<P> quadratic_u;
    for (typename P::T &coefficient : quadratic_u)
        coefficient = (static_cast<typename P::T>(word(engine)) << 64) |
                      word(engine);
    Polynomial<P> secret_square;
    Polynomial<P> u_times_secret;
    Polynomial<P> quadratic_v;
    PolyMul<P>(secret_square, secret, secret);
    PolyMul<P>(u_times_secret, quadratic_u, secret);
    for (std::size_t i = 0; i < P::n; i++)
        quadratic_v[i] = secret_square[i] - u_times_secret[i];

    std::array<std::uint64_t, P::n> left_slots;
    std::array<std::uint64_t, P::n> right_slots;
    left_slots.fill(2);
    right_slots.fill(3);
    TRLWE<P> left;
    TRLWE<P> right;
    trlweSlotEncrypt<P>(left, left_slots, key);
    trlweSlotEncrypt<P>(right, right_slots, key);

    TRLWE<P> product;
    TRLWEQuadraticHintMultFullDD<P>(product, left, right, quadratic_u,
                                    quadratic_v);
    std::array<std::uint64_t, P::n> result_slots;
    trlweSlotDecrypt<P>(result_slots, product, key);
    for (const std::uint64_t value : result_slots) assert(value == 6);

    std::vector<TRLWE<P>> factors(4);
    for (std::size_t factor = 0; factor < factors.size(); factor++) {
        std::array<std::uint64_t, P::n> slots;
        slots.fill(static_cast<std::uint64_t>(factor + 2));
        trlweSlotEncrypt<P>(factors[factor], slots, key);
    }
    std::uint32_t layers = 0;
    std::uint32_t multiplications = 0;
    const TRLWE<P> tree_product = TRLWEQuadraticHintProductTreeFullDD<P>(
        std::move(factors), quadratic_u, quadratic_v, &layers,
        &multiplications);
    assert(layers == 2);
    assert(multiplications == 3);
    trlweSlotDecrypt<P>(result_slots, tree_product, key);
    for (const std::uint64_t value : result_slots) assert(value == 120);

    auto encrypt_coefficient_bit = [&](const bool bit) {
        Polynomial<P> encoded{};
        if (bit) encoded[0] = P::delta_int;
        TRLWE<P> ciphertext;
        trlweSymEncrypt<P>(ciphertext, encoded, key);
        return ciphertext;
    };
    const std::array<bool, 4> selector_bits = {true, false, true, true};
    const std::array<std::uint32_t, 4> exponents = {3, 5, 7, 11};
    std::vector<TRLWE<P>> monomial_factors(4);
    for (std::size_t i = 0; i < monomial_factors.size(); i++) {
        const TRLWE<P> encrypted_bit = encrypt_coefficient_bit(selector_bits[i]);
        TRLWEMonomialFactorBFV<P>(monomial_factors[i], encrypted_bit,
                                  exponents[i]);
        const Polynomial<P> factor_phase = trlwePhase<P>(monomial_factors[i], key);
        const std::uint32_t expected_index = selector_bits[i] ? exponents[i] : 0;
        for (std::size_t coefficient = 0; coefficient < P::n; coefficient++) {
            const std::uint64_t decoded = bfvDecodeCoeff<P>(factor_phase[coefficient]);
            assert(decoded == (coefficient == expected_index ? 1U : 0U));
        }
    }
    const TRLWE<P> monomial_product = TRLWEQuadraticHintProductTreeFullDD<P>(
        std::move(monomial_factors), quadratic_u, quadratic_v);
    const Polynomial<P> monomial_phase = trlwePhase<P>(monomial_product, key);
    constexpr std::uint32_t expected_exponent = 3 + 7 + 11;
    for (std::size_t coefficient = 0; coefficient < P::n; coefficient++) {
        const std::uint64_t decoded = bfvDecodeCoeff<P>(monomial_phase[coefficient]);
        assert(decoded == (coefficient == expected_exponent ? 1U : 0U));
    }

    constexpr std::size_t bucket_size = 8;
    constexpr std::size_t selected_slot = 3;
    constexpr std::uint32_t bucket_scale = 2;
    std::vector<TRLWE<P>> selected_entries(bucket_size);
    std::vector<std::uint32_t> bucket_masks(bucket_size);
    for (std::size_t i = 0; i < bucket_size; i++) {
        selected_entries[i] = encrypt_coefficient_bit(i == selected_slot);
        bucket_masks[i] = static_cast<std::uint32_t>(2 * i + 1);
    }
    const TRLWE<P> encrypted_zero = encrypt_coefficient_bit(false);
    const TRLWE<P> selected_bucket = TRLWEPBCPreparedFactorBFV<P>(
        selected_entries, encrypted_zero, bucket_masks, bucket_scale);
    const Polynomial<P> selected_bucket_phase = trlwePhase<P>(selected_bucket, key);
    const std::uint32_t selected_exponent =
        bucket_masks[selected_slot] * bucket_scale;
    for (std::size_t coefficient = 0; coefficient < P::n; coefficient++) {
        const std::uint64_t decoded =
            bfvDecodeCoeff<P>(selected_bucket_phase[coefficient]);
        assert(decoded == (coefficient == selected_exponent ? 1U : 0U));
    }

    std::vector<TRLWE<P>> empty_entries(bucket_size);
    for (TRLWE<P> &entry : empty_entries)
        entry = encrypt_coefficient_bit(false);
    const TRLWE<P> encrypted_one = encrypt_coefficient_bit(true);
    const TRLWE<P> dummy_bucket = TRLWEPBCPreparedFactorBFV<P>(
        empty_entries, encrypted_one, bucket_masks, bucket_scale);
    const Polynomial<P> dummy_bucket_phase = trlwePhase<P>(dummy_bucket, key);
    for (std::size_t coefficient = 0; coefficient < P::n; coefficient++) {
        const std::uint64_t decoded =
            bfvDecodeCoeff<P>(dummy_bucket_phase[coefficient]);
        assert(decoded == (coefficient == 0 ? 1U : 0U));
    }




}
