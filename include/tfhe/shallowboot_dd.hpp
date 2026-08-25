#pragma once

#include <cstdint>

#include <boost/multiprecision/cpp_int.hpp>

#include <bfv/bfv++.hpp>
#include <bfv/bfv-slots.hpp>
#include <tfhe/shallowboot_lowdepth.hpp>

namespace TFHEpp::shallowboot::dd {

template <class P, std::uint32_t ModulusBits>
inline Polynomial<P> PolynomialMulPowerOfTwo(const Polynomial<P> &left,
                                             const Polynomial<P> &right)
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>);
    static_assert(ModulusBits > 0 && ModulusBits < 128);
    constexpr __uint128_t mask =
        (static_cast<__uint128_t>(1) << ModulusBits) - 1;
    Polynomial<P> result{};
    for (std::size_t i = 0; i < P::n; i++) {
        __uint128_t accumulator = 0;
        for (std::size_t j = 0; j <= i; j++)
            accumulator += left[j] * right[i - j];
        for (std::size_t j = i + 1; j < P::n; j++)
            accumulator -= left[j] * right[P::n + i - j];
        result[i] = accumulator & mask;
    }
    return result;
}

template <class HighP, class LowP, std::uint32_t LowBits>
inline TRLWE<LowP> ModulusBoundaryToPowerOfTwoDD(const TRLWE<HighP> &high)
{
    static_assert(std::is_same_v<typename HighP::T, __uint128_t> &&
                  std::is_same_v<typename LowP::T, __uint128_t>);
    static_assert(HighP::n == LowP::n && HighP::k == 1 && LowP::k == 1);
    static_assert(LowBits > 0 && LowBits < 128);
    constexpr std::uint32_t shift = 128 - LowBits;
    constexpr __uint128_t rounding = static_cast<__uint128_t>(1)
                                      << (shift - 1);
    constexpr __uint128_t mask =
        (static_cast<__uint128_t>(1) << LowBits) - 1;
    TRLWE<LowP> result;
    for (std::size_t component = 0; component <= LowP::k; component++)
        for (std::size_t i = 0; i < LowP::n; i++)
            result[component][i] = ((high[component][i] + rounding) >> shift) &
                                   mask;
    return result;
}

template <class LowP, std::uint32_t LowBits>
inline void QuadraticHintMultiplyPowerOfTwoDD(
    TRLWE<LowP> &result, const TRLWE<LowP> &left,
    const TRLWE<LowP> &right, const Polynomial<LowP> &quadratic_u,
    const Polynomial<LowP> &quadratic_v)
{
    static_assert(std::is_same_v<typename LowP::T, __uint128_t> &&
                  LowP::k == 1);
    constexpr __uint128_t mask =
        (static_cast<__uint128_t>(1) << LowBits) - 1;
    TRLWE3<LowP> tensor;
    TRLWEMultWithoutRelinearizationFullDD<LowP>(tensor, left, right);
    Polynomial<LowP> square{};
    for (std::size_t i = 0; i < LowP::n; i++) {
        tensor[0][i] &= mask;
        tensor[1][i] &= mask;
        square[i] = tensor[2][i] & mask;
    }
    const Polynomial<LowP> square_u =
        PolynomialMulPowerOfTwo<LowP, LowBits>(square, quadratic_u);
    const Polynomial<LowP> square_v =
        PolynomialMulPowerOfTwo<LowP, LowBits>(square, quadratic_v);
    for (std::size_t i = 0; i < LowP::n; i++) {
        result[0][i] = (tensor[0][i] - square_u[i]) & mask;
        result[1][i] = (tensor[1][i] + square_v[i]) & mask;
    }
}

template <class LowP, std::uint32_t LowBits>
inline TRLWE<LowP> QuadraticHintProductTreePowerOfTwoDD(
    std::vector<TRLWE<LowP>> inputs, const Polynomial<LowP> &quadratic_u,
    const Polynomial<LowP> &quadratic_v)
{
    if (inputs.empty())
        throw std::invalid_argument("low DD product tree is empty");
    while (inputs.size() > 1) {
        const std::size_t pairs = inputs.size() / 2;
        std::vector<TRLWE<LowP>> next((inputs.size() + 1) / 2);
#pragma omp parallel for if (pairs > 1)
        for (std::int64_t pair = 0;
             pair < static_cast<std::int64_t>(pairs); pair++)
            QuadraticHintMultiplyPowerOfTwoDD<LowP, LowBits>(
                next[static_cast<std::size_t>(pair)],
                inputs[2 * static_cast<std::size_t>(pair)],
                inputs[2 * static_cast<std::size_t>(pair) + 1], quadratic_u,
                quadratic_v);
        if ((inputs.size() & 1U) != 0)
            next.back() = std::move(inputs.back());
        inputs = std::move(next);
    }
    return std::move(inputs.front());
}

// Algorithm-3 D2/QH-SS boundary: the high secret is already coefficient-small,
// so BFV modulus switching can be applied directly.  Convert the resulting
// coefficient arrays into the explicit low-prime NTT representation used by
// the remaining shallow product levels.
template <class HighP>
inline lowdepth::Ciphertext ModulusBoundaryToLowNTT(
    const TRLWE<HighP> &high, const lowdepth::Ring &low_ring)
{
    if (HighP::n != low_ring.degree())
        throw std::invalid_argument("DD boundary changes ring degree");
    TRLWE<HighP> reduced;
    TFHEpp::ModulusSwitch<HighP>(reduced, high, low_ring.modulus());
    lowdepth::Ciphertext result;
    result.a.resize(low_ring.degree());
    result.b.resize(low_ring.degree());
    for (std::size_t i = 0; i < low_ring.degree(); i++) {
        result.a[i] = static_cast<std::uint64_t>(reduced[0][i]);
        result.b[i] = static_cast<std::uint64_t>(reduced[HighP::k][i]);
    }
    low_ring.forward(result.a);
    low_ring.forward(result.b);
    return result;
}

template <class HighP>
inline lowdepth::RNSCiphertext ModulusBoundaryToRNS(
    const TRLWE<HighP> &high, const lowdepth::RNSRing &target)
{
    using boost::multiprecision::uint256_t;
    static_assert(std::is_same_v<typename HighP::T, __uint128_t> &&
                  HighP::k == 1);
    if (target.levels() != 2 || HighP::n != target.degree())
        throw std::invalid_argument(
            "DD-to-RNS boundary requires two equal-degree primes");
    const unsigned __int128 target_modulus =
        static_cast<unsigned __int128>(target[0].modulus()) *
        target[1].modulus();
    lowdepth::RNSCiphertext result;
    result.residues.resize(2);
    for (auto &residue : result.residues) {
        residue.a.resize(target.degree());
        residue.b.resize(target.degree());
    }
    for (std::size_t i = 0; i < target.degree(); i++) {
        for (std::size_t component = 0; component < 2; component++) {
            const __uint128_t value = high[component][i];
            const uint256_t numerator =
                uint256_t(value) * target_modulus +
                (uint256_t(1) << 127);
            const unsigned __int128 rounded =
                static_cast<unsigned __int128>(numerator >> 128) %
                target_modulus;
            for (std::size_t level = 0; level < 2; level++) {
                const std::uint64_t residue = static_cast<std::uint64_t>(
                    rounded % target[level].modulus());
                if (component == 0)
                    result.residues[level].a[i] = residue;
                else
                    result.residues[level].b[i] = residue;
            }
        }
    }
    for (std::size_t level = 0; level < 2; level++) {
        target[level].forward(result.residues[level].a);
        target[level].forward(result.residues[level].b);
    }
    return result;
}

inline lowdepth::RNSCiphertext QuadraticHintBFVMultiplyRNS(
    const lowdepth::RNSRing &ring, const lowdepth::RNSSecret &secret,
    const lowdepth::RNSCiphertext &left,
    const lowdepth::RNSCiphertext &right,
    const std::uint64_t plaintext_modulus)
{
    using boost::multiprecision::int256_t;
    if (ring.levels() != 2 || plaintext_modulus < 2 ||
        left.residues.size() != 2 || right.residues.size() != 2)
        throw std::invalid_argument("invalid CRT-wide BFV/QH multiplication");
    lowdepth::validate(ring, secret);
    const modular_ntt::TwoPrimeCRT crt(ring.prime(0), ring.prime(1));
    const unsigned __int128 modulus = crt.modulusProduct();
    const unsigned __int128 delta = modulus / plaintext_modulus;
    struct Components {
        std::array<std::vector<std::uint64_t>, 2> a;
        std::array<std::vector<std::uint64_t>, 2> b;
    } lhs, rhs;
    for (std::size_t level = 0; level < 2; level++) {
        lhs.a[level] = left.residues[level].a;
        lhs.b[level] = left.residues[level].b;
        rhs.a[level] = right.residues[level].a;
        rhs.b[level] = right.residues[level].b;
        ring[level].inverse(lhs.a[level]);
        ring[level].inverse(lhs.b[level]);
        ring[level].inverse(rhs.a[level]);
        ring[level].inverse(rhs.b[level]);
    }
    auto reconstruct = [&](const std::array<std::vector<std::uint64_t>, 2> &x) {
        std::vector<unsigned __int128> result(ring.degree());
        for (std::size_t i = 0; i < ring.degree(); i++) {
            const __int128 signed_value =
                crt.reconstructSigned(x[0][i], x[1][i]);
            result[i] = signed_value >= 0
                ? static_cast<unsigned __int128>(signed_value)
                : modulus -
                      (static_cast<unsigned __int128>(-(signed_value + 1)) + 1);
        }
        return result;
    };
    const auto la = reconstruct(lhs.a);
    const auto lb = reconstruct(lhs.b);
    const auto ra = reconstruct(rhs.a);
    const auto rb = reconstruct(rhs.b);
    auto convolution = [&](const std::vector<unsigned __int128> &x,
                           const std::vector<unsigned __int128> &y) {
        std::vector<int256_t> result(ring.degree());
        for (std::size_t i = 0; i < ring.degree(); i++) {
            int256_t accumulator = 0;
            for (std::size_t j = 0; j <= i; j++)
                accumulator += int256_t(x[j]) * y[i - j];
            for (std::size_t j = i + 1; j < ring.degree(); j++)
                accumulator -= int256_t(x[j]) *
                               y[ring.degree() + i - j];
            result[i] = accumulator;
        }
        return result;
    };
    auto aa = convolution(la, ra);
    auto cross = convolution(la, rb);
    const auto other_cross = convolution(lb, ra);
    auto bb = convolution(lb, rb);
    for (std::size_t i = 0; i < ring.degree(); i++)
        cross[i] += other_cross[i];
    auto scale_and_reduce = [&](const int256_t &value) {
        const bool negative = value < 0;
        const int256_t magnitude = negative ? -value : value;
        int256_t rounded = (magnitude + int256_t(delta) / 2) / delta;
        rounded %= modulus;
        unsigned __int128 residue =
            static_cast<unsigned __int128>(rounded);
        return negative && residue != 0 ? modulus - residue : residue;
    };
    std::array<std::vector<std::uint64_t>, 2> square;
    lowdepth::RNSCiphertext result;
    result.residues.resize(2);
    for (std::size_t level = 0; level < 2; level++) {
        square[level].resize(ring.degree());
        result.residues[level].a.resize(ring.degree());
        result.residues[level].b.resize(ring.degree());
    }
    for (std::size_t i = 0; i < ring.degree(); i++) {
        const unsigned __int128 square_value = scale_and_reduce(aa[i]);
        const unsigned __int128 cross_value = scale_and_reduce(cross[i]);
        const unsigned __int128 bb_value = scale_and_reduce(bb[i]);
        for (std::size_t level = 0; level < 2; level++) {
            const std::uint64_t prime = ring[level].modulus();
            square[level][i] = static_cast<std::uint64_t>(square_value % prime);
            result.residues[level].a[i] =
                static_cast<std::uint64_t>(cross_value % prime);
            result.residues[level].b[i] =
                static_cast<std::uint64_t>(bb_value % prime);
        }
    }
    for (std::size_t level = 0; level < 2; level++) {
        ring[level].forward(square[level]);
        ring[level].forward(result.residues[level].a);
        ring[level].forward(result.residues[level].b);
        const std::uint64_t prime = ring[level].modulus();
        for (std::size_t i = 0; i < ring.degree(); i++) {
            result.residues[level].a[i] = modular_ntt::subtract(
                result.residues[level].a[i],
                modular_ntt::multiply(square[level][i],
                    secret.residues[level].quadratic_u[i], prime), prime);
            result.residues[level].b[i] = modular_ntt::add(
                result.residues[level].b[i],
                modular_ntt::multiply(square[level][i],
                    secret.residues[level].quadratic_v[i], prime), prime);
        }
    }
    return result;
}

inline lowdepth::RNSCiphertext QuadraticHintBFVProductTreeRNS(
    const lowdepth::RNSRing &ring, const lowdepth::RNSSecret &secret,
    std::vector<lowdepth::RNSCiphertext> inputs,
    const std::uint64_t plaintext_modulus)
{
    if (inputs.empty())
        throw std::invalid_argument("RNS BFV/QH product tree is empty");
    while (inputs.size() > 1) {
        const std::size_t pairs = inputs.size() / 2;
        std::vector<lowdepth::RNSCiphertext> next((inputs.size() + 1) / 2);
#pragma omp parallel for if (pairs > 1)
        for (std::int64_t pair = 0;
             pair < static_cast<std::int64_t>(pairs); pair++)
            next[static_cast<std::size_t>(pair)] =
                QuadraticHintBFVMultiplyRNS(
                    ring, secret, inputs[2 * static_cast<std::size_t>(pair)],
                    inputs[2 * static_cast<std::size_t>(pair) + 1],
                    plaintext_modulus);
        if ((inputs.size() & 1U) != 0)
            next.back() = std::move(inputs.back());
        inputs = std::move(next);
    }
    return std::move(inputs.front());
}

inline lowdepth::RNSCiphertext QuadraticHintBGVMultiplyRNS(
    const lowdepth::RNSRing &ring, const lowdepth::RNSSecret &secret,
    const lowdepth::RNSCiphertext &left,
    const lowdepth::RNSCiphertext &right,
    const std::uint64_t plaintext_modulus)
{
    lowdepth::RNSCiphertext result =
        lowdepth::RelinearizationFreeMultiply(ring, secret, left, right);
    const unsigned __int128 modulus =
        static_cast<unsigned __int128>(ring[0].modulus()) *
        ring[1].modulus();
    const unsigned __int128 delta = modulus / plaintext_modulus;
    for (std::size_t level = 0; level < ring.levels(); level++) {
        const std::uint64_t prime = ring[level].modulus();
        const std::uint64_t inverse_delta = modular_ntt::invert(
            static_cast<std::uint64_t>(delta % prime), prime);
        for (std::size_t i = 0; i < ring.degree(); i++) {
            result.residues[level].a[i] = modular_ntt::multiply(
                result.residues[level].a[i], inverse_delta, prime);
            result.residues[level].b[i] = modular_ntt::multiply(
                result.residues[level].b[i], inverse_delta, prime);
        }
    }
    return result;
}

inline lowdepth::RNSCiphertext QuadraticHintBGVProductTreeRNS(
    const lowdepth::RNSRing &ring, const lowdepth::RNSSecret &secret,
    std::vector<lowdepth::RNSCiphertext> inputs,
    const std::uint64_t plaintext_modulus)
{
    if (inputs.empty())
        throw std::invalid_argument("RNS BGV/QH product tree is empty");
    while (inputs.size() > 1) {
        const std::size_t pairs = inputs.size() / 2;
        std::vector<lowdepth::RNSCiphertext> next((inputs.size() + 1) / 2);
#pragma omp parallel for if (pairs > 1)
        for (std::int64_t pair = 0;
             pair < static_cast<std::int64_t>(pairs); pair++)
            next[static_cast<std::size_t>(pair)] =
                QuadraticHintBGVMultiplyRNS(
                    ring, secret, inputs[2 * static_cast<std::size_t>(pair)],
                    inputs[2 * static_cast<std::size_t>(pair) + 1],
                    plaintext_modulus);
        if ((inputs.size() & 1U) != 0)
            next.back() = std::move(inputs.back());
        inputs = std::move(next);
    }
    return std::move(inputs.front());
}

inline lowdepth::Ciphertext QuadraticHintBFVMultiply(
    const lowdepth::Ring &ring, const lowdepth::Secret &secret,
    const lowdepth::Ciphertext &left, const lowdepth::Ciphertext &right,
    const std::uint64_t plaintext_modulus)
{
    if (plaintext_modulus < 2 || plaintext_modulus >= ring.modulus())
        throw std::invalid_argument("invalid low BFV plaintext modulus");
    const std::uint64_t modulus = ring.modulus();
    const std::uint64_t delta = modulus / plaintext_modulus;
    std::array<std::vector<std::uint64_t>, 4> coefficient_inputs = {
        left.a, left.b, right.a, right.b};
    for (auto &polynomial : coefficient_inputs) ring.inverse(polynomial);
    auto convolution = [&](const std::vector<std::uint64_t> &lhs,
                           const std::vector<std::uint64_t> &rhs) {
        using boost::multiprecision::int256_t;
        std::vector<int256_t> result(ring.degree());
        for (std::size_t i = 0; i < ring.degree(); i++) {
            int256_t accumulator = 0;
            for (std::size_t j = 0; j <= i; j++)
                accumulator += int256_t(
                    modular_ntt::centeredResidue(lhs[j], modulus)) *
                    modular_ntt::centeredResidue(rhs[i - j], modulus);
            for (std::size_t j = i + 1; j < ring.degree(); j++)
                accumulator -= int256_t(
                    modular_ntt::centeredResidue(lhs[j], modulus)) *
                    modular_ntt::centeredResidue(
                        rhs[ring.degree() + i - j], modulus);
            result[i] = accumulator;
        }
        return result;
    };
    auto aa = convolution(coefficient_inputs[0], coefficient_inputs[2]);
    auto ab = convolution(coefficient_inputs[0], coefficient_inputs[3]);
    auto ba = convolution(coefficient_inputs[1], coefficient_inputs[2]);
    auto bb = convolution(coefficient_inputs[1], coefficient_inputs[3]);
    auto rounded_residue = [&](const boost::multiprecision::int256_t &value) {
        using boost::multiprecision::int256_t;
        const bool negative = value < 0;
        const int256_t magnitude = negative ? -value : value;
        const std::uint64_t rounded = static_cast<std::uint64_t>(
            ((magnitude + delta / 2) / delta) % modulus);
        return negative ? modular_ntt::negate(rounded, modulus) : rounded;
    };
    lowdepth::Ciphertext square;
    lowdepth::Ciphertext result;
    square.a.resize(ring.degree());
    result.a.resize(ring.degree());
    result.b.resize(ring.degree());
    for (std::size_t i = 0; i < ring.degree(); i++) {
        square.a[i] = rounded_residue(aa[i]);
        result.a[i] = rounded_residue(ab[i] + ba[i]);
        result.b[i] = rounded_residue(bb[i]);
    }
    ring.forward(square.a);
    ring.forward(result.a);
    ring.forward(result.b);
    for (std::size_t i = 0; i < ring.degree(); i++) {
        result.a[i] = modular_ntt::subtract(
            result.a[i],
            modular_ntt::multiply(square.a[i], secret.quadratic_u[i], modulus),
            modulus);
        result.b[i] = modular_ntt::add(
            result.b[i],
            modular_ntt::multiply(square.a[i], secret.quadratic_v[i], modulus),
            modulus);
    }
    return result;
}

inline lowdepth::Ciphertext QuadraticHintBFVProductTree(
    const lowdepth::Ring &ring, const lowdepth::Secret &secret,
    std::vector<lowdepth::Ciphertext> inputs,
    const std::uint64_t plaintext_modulus,
    std::uint32_t *layers = nullptr,
    std::uint32_t *multiplications = nullptr)
{
    if (inputs.empty())
        throw std::invalid_argument("low BFV product tree is empty");
    std::uint32_t local_layers = 0;
    std::uint32_t local_multiplications = 0;
    while (inputs.size() > 1) {
        const std::size_t pairs = inputs.size() / 2;
        std::vector<lowdepth::Ciphertext> next((inputs.size() + 1) / 2);
#pragma omp parallel for if (pairs > 1)
        for (std::int64_t pair = 0;
             pair < static_cast<std::int64_t>(pairs); pair++)
            next[static_cast<std::size_t>(pair)] = QuadraticHintBFVMultiply(
                ring, secret, inputs[2 * static_cast<std::size_t>(pair)],
                inputs[2 * static_cast<std::size_t>(pair) + 1],
                plaintext_modulus);
        if ((inputs.size() & 1U) != 0)
            next.back() = std::move(inputs.back());
        inputs = std::move(next);
        local_layers++;
        local_multiplications += static_cast<std::uint32_t>(pairs);
    }
    if (layers != nullptr) *layers = local_layers;
    if (multiplications != nullptr) *multiplications = local_multiplications;
    return std::move(inputs.front());
}

// Low-prime BFV/QH multiplication using TFHEpp's existing FullDD scale-and-
// round implementation. LowP::delta_int must equal floor(q_low/t). Inputs are
// converted from the low-prime NTT representation to coefficient TRLWE, the
// tensor is scaled by FullDD, then reduced modulo q_low before QH collapse.
template <class LowP>
inline lowdepth::Ciphertext QuadraticHintBFVMultiplyFullDD(
    const lowdepth::Ring &ring, const lowdepth::Secret &secret,
    const lowdepth::Ciphertext &left, const lowdepth::Ciphertext &right)
{
    static_assert(LowP::k == 1 && std::is_same_v<typename LowP::T, __uint128_t>);
    if (LowP::n != ring.degree())
        throw std::invalid_argument("low FullDD ring degree mismatch");
    auto to_trlwe = [&](const lowdepth::Ciphertext &input) {
        TRLWE<LowP> result;
        std::vector<std::uint64_t> a = input.a;
        std::vector<std::uint64_t> b = input.b;
        ring.inverse(a);
        ring.inverse(b);
        for (std::size_t i = 0; i < ring.degree(); i++) {
            result[0][i] = a[i];
            result[1][i] = b[i];
        }
        return result;
    };
    const TRLWE<LowP> left_trlwe = to_trlwe(left);
    const TRLWE<LowP> right_trlwe = to_trlwe(right);
    TRLWE3<LowP> tensor;
    TRLWEMultWithoutRelinearizationFullDD<LowP>(tensor, left_trlwe,
                                                right_trlwe);
    const std::uint64_t modulus = ring.modulus();
    std::array<std::vector<std::uint64_t>, 3> reduced;
    for (auto &component : reduced) component.resize(ring.degree());
    for (std::size_t component = 0; component < 3; component++)
        for (std::size_t i = 0; i < ring.degree(); i++) {
            const __int128_t centered = static_cast<__int128_t>(tensor[component][i]);
            __int128_t residue = centered % static_cast<__int128_t>(modulus);
            if (residue < 0) residue += modulus;
            reduced[component][i] = static_cast<std::uint64_t>(residue);
        }
    for (auto &component : reduced) ring.forward(component);
    lowdepth::Ciphertext result;
    result.a.resize(ring.degree());
    result.b.resize(ring.degree());
    for (std::size_t i = 0; i < ring.degree(); i++) {
        result.a[i] = modular_ntt::subtract(
            reduced[0][i],
            modular_ntt::multiply(reduced[2][i], secret.quadratic_u[i],
                                  modulus),
            modulus);
        result.b[i] = modular_ntt::add(
            reduced[1][i],
            modular_ntt::multiply(reduced[2][i], secret.quadratic_v[i],
                                  modulus),
            modulus);
    }
    return result;
}

}  // namespace TFHEpp::shallowboot::dd
