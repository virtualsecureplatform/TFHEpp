#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "ckks/ckks.hpp"
#include "params.hpp"
#include "tfhe/key.hpp"

// The approximate GL implementation uses the active-level and multi-limb
// arithmetic introduced by the CKKS implementation.  Those facilities are
// available with TFHEpp's default 128-bit parameter family.
#ifdef TFHEPP_DEFAULT_128BIT_PARAMS

namespace TFHEpp {

namespace gl_detail {

consteval bool isPrime(const std::uint32_t value)
{
    if (value < 2) return false;
    for (std::uint32_t divisor = 2; divisor * divisor <= value; divisor++)
        if (value % divisor == 0) return false;
    return true;
}

consteval bool isPowerOfTwo(const std::uint32_t value)
{
    return value != 0 && (value & (value - 1)) == 0;
}

constexpr std::uint32_t powMod(std::uint32_t base, std::uint32_t exponent,
                               const std::uint32_t modulus)
{
    std::uint64_t result = 1;
    std::uint64_t power = base % modulus;
    while (exponent != 0) {
        if ((exponent & 1U) != 0) result = (result * power) % modulus;
        power = (power * power) % modulus;
        exponent >>= 1;
    }
    return static_cast<std::uint32_t>(result);
}

constexpr std::uint32_t primitiveRootPrime(const std::uint32_t prime)
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
        for (std::size_t i = 0; i < factor_count; i++) {
            if (powMod(candidate, order / factors[i], prime) == 1) {
                primitive = false;
                break;
            }
        }
        if (primitive) return candidate;
    }
    return 0;
}

template <class P>
inline std::int64_t torusToSignedSmall(const typename P::T value)
{
    if constexpr (is_multilimb_uint_v<typename P::T>) {
        return multilimb_to_signed_i64(value);
    }
    else {
        static_assert(std::is_same_v<typename P::T, __uint128_t>);
        const __int128_t signed_value = static_cast<__int128_t>(value);
        if (signed_value < std::numeric_limits<std::int64_t>::min() ||
            signed_value > std::numeric_limits<std::int64_t>::max())
            throw std::overflow_error("GL digit product exceeds int64_t");
        return static_cast<std::int64_t>(signed_value);
    }
}

template <class T>
inline void addSigned(T &destination, const T value, const bool negative)
{
    if (negative)
        destination -= value;
    else
        destination += value;
}

}  // namespace gl_detail

// GL works over
//   Z_q[I,X,Y,W]/(I^2+1, X^n-I, Y^n-I, Phi_p(W)).
// A Y coefficient is an RLWE polynomial of degree 2*n*phi(p).  BaseP::n is
// therefore the dimension of one (I,X,W) slice, not the matrix dimension.
template <class BaseP, std::uint32_t MatrixDimension,
          std::uint32_t CyclotomicOrder, std::uint32_t AuxiliaryLogQ = 32>
struct GLParameter {
    using baseP = BaseP;
    using T = typename BaseP::T;

    static constexpr std::uint32_t matrix_dimension = MatrixDimension;
    static constexpr std::uint32_t cyclotomic_order = CyclotomicOrder;
    static constexpr std::uint32_t phi = CyclotomicOrder - 1;
    static constexpr std::uint32_t slice_dimension = 2 * MatrixDimension * phi;
    static constexpr std::uint32_t auxiliary_log_q = AuxiliaryLogQ;
    static constexpr std::uint32_t primitive_root =
        gl_detail::primitiveRootPrime(CyclotomicOrder);

    static_assert(BaseP::k == 1,
                  "GL currently uses one RLWE mask polynomial per Y slice");
    static_assert(gl_detail::isPowerOfTwo(MatrixDimension),
                  "the GL X/Y dimension must be a power of two");
    static_assert(gl_detail::isPrime(CyclotomicOrder) &&
                      (CyclotomicOrder & 1U) != 0,
                  "this GL implementation currently supports odd prime p");
    static_assert(BaseP::n == slice_dimension,
                  "BaseP::n must equal 2*n*phi(p)");
    static_assert(primitive_root != 0);
    static_assert(AuxiliaryLogQ > 0,
                  "GL key switching requires a nonzero auxiliary modulus");
    static_assert(ckks_detail::supported_torus_v<BaseP>,
                  "GL requires the CKKS 128-bit or multi-limb torus backend");
};

// The parameter shape used in the GL paper: n=256, p=17 and slice dimension
// 8192.  lvl4param supplies the matching 128-bit storage/DD backend.  This is
// a compile-time shape baseline; the coefficient kernels below prioritize
// correctness and are not the paper's optimized full-size NTT kernels.
using GL256p17Parameter = GLParameter<lvl4param, 256, 17>;

// GL's reference parameter set uses coefficient-domain Gaussian noise with
// standard deviation 3.2.  CKKSNoise stores a modular standard deviation, so
// convert an absolute coefficient standard deviation at the active level.
template <std::uint32_t LogQ>
inline CKKSNoise GLNoiseAtLevel(const double coefficient_stdev = 3.2)
{
    return {std::ldexp(coefficient_stdev, -static_cast<int>(LogQ)), 0};
}

template <class GLP>
using GLBasePolynomial = Polynomial<typename GLP::baseP>;

template <class GLP>
class GLPolynomial {
public:
    using BasePolynomial = GLBasePolynomial<GLP>;

    GLPolynomial() : coefficients_(GLP::matrix_dimension) {}

    BasePolynomial &operator[](const std::size_t y) { return coefficients_[y]; }
    const BasePolynomial &operator[](const std::size_t y) const
    {
        return coefficients_[y];
    }

    constexpr std::size_t size() const { return GLP::matrix_dimension; }
    auto begin() { return coefficients_.begin(); }
    auto end() { return coefficients_.end(); }
    auto begin() const { return coefficients_.begin(); }
    auto end() const { return coefficients_.end(); }

private:
    std::vector<BasePolynomial> coefficients_;
};

template <class GLP>
struct GLCiphertextData {
    // Paper convention: body + mask*s decrypts the ciphertext.
    static constexpr std::size_t body = 0;
    static constexpr std::size_t mask = 1;

    std::array<GLPolynomial<GLP>, 2> data{};

    GLPolynomial<GLP> &operator[](const std::size_t component)
    {
        return data[component];
    }
    const GLPolynomial<GLP> &operator[](const std::size_t component) const
    {
        return data[component];
    }
};

template <class GLP>
using GLBaseCiphertextData = std::array<GLBasePolynomial<GLP>, 2>;

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
struct GLCiphertext {
    using Parameter = GLP;
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t log_budget = LogQ - LogDelta;

    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename GLP::baseP::T>::digits);
    static_assert(LogDelta < LogQ);

    GLCiphertextData<GLP> ct{};

    GLPolynomial<GLP> &operator[](const std::size_t component)
    {
        return ct[component];
    }
    const GLPolynomial<GLP> &operator[](const std::size_t component) const
    {
        return ct[component];
    }
};

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
struct GLPlaintext {
    using Parameter = GLP;
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;

    GLPolynomial<GLP> poly{};
};

template <class GLP>
class GLMatrixBatch {
public:
    using value_type = std::complex<double>;

    GLMatrixBatch()
        : values_(static_cast<std::size_t>(GLP::phi) * GLP::matrix_dimension *
                  GLP::matrix_dimension)
    {
    }

    value_type &operator()(const std::size_t batch, const std::size_t row,
                           const std::size_t column)
    {
        return values_[index(batch, row, column)];
    }
    const value_type &operator()(const std::size_t batch, const std::size_t row,
                                 const std::size_t column) const
    {
        return values_[index(batch, row, column)];
    }

    constexpr std::size_t batchCount() const { return GLP::phi; }
    constexpr std::size_t matrixDimension() const
    {
        return GLP::matrix_dimension;
    }

private:
    static std::size_t index(const std::size_t batch, const std::size_t row,
                             const std::size_t column)
    {
        if (batch >= GLP::phi || row >= GLP::matrix_dimension ||
            column >= GLP::matrix_dimension)
            throw std::out_of_range("GL matrix-batch index out of range");
        return (batch * GLP::matrix_dimension + row) * GLP::matrix_dimension +
               column;
    }

    std::vector<value_type> values_;
};

namespace gl_detail {

template <class GLP>
constexpr std::size_t baseIndex(const std::uint32_t gaussian,
                                const std::uint32_t x, const std::uint32_t w)
{
    return (static_cast<std::size_t>(w) * GLP::matrix_dimension + x) * 2 +
           gaussian;
}

template <class GLP>
inline void clear(GLBasePolynomial<GLP> &poly)
{
    poly.fill(typename GLP::T{0});
}

template <class GLP>
inline void clear(GLPolynomial<GLP> &poly)
{
    for (auto &slice : poly) clear<GLP>(slice);
}

template <class GLP, std::uint32_t LogQ>
inline void reduce(GLBasePolynomial<GLP> &poly)
{
    using P = typename GLP::baseP;
    for (auto &coefficient : poly)
        coefficient = ckks_detail::reduceToLevel<P, LogQ>(coefficient);
}

template <class GLP, std::uint32_t LogQ>
inline void reduce(GLPolynomial<GLP> &poly)
{
    for (auto &slice : poly) reduce<GLP, LogQ>(slice);
}

template <class GLP, std::uint32_t InputLogQ, std::uint32_t Shift>
inline typename GLP::T divideRoundLevel(typename GLP::T value)
{
    using P = typename GLP::baseP;
    using T = typename GLP::T;
    static_assert(Shift > 0);
    static_assert(InputLogQ > Shift);
    static_assert(InputLogQ <= std::numeric_limits<T>::digits);

    value = ckks_detail::reduceToLevel<P, InputLogQ>(value);
    const bool negative = (value & (T{1} << (InputLogQ - 1))) != T{0};
    T magnitude;
    if (!negative)
        magnitude = value;
    else if constexpr (InputLogQ == std::numeric_limits<T>::digits)
        magnitude = T{0} - value;
    else
        magnitude = (T{1} << InputLogQ) - value;

    magnitude += T{1} << (Shift - 1);
    const T quotient = magnitude >> Shift;
    return ckks_detail::reduceToLevel<P, InputLogQ - Shift>(
        negative ? T{0} - quotient : quotient);
}

template <class GLP, std::uint32_t InputLogQ, std::uint32_t Shift>
inline void divideRoundLevel(GLPolynomial<GLP> &poly)
{
    for (auto &slice : poly)
        for (auto &coefficient : slice)
            coefficient = divideRoundLevel<GLP, InputLogQ, Shift>(coefficient);
}

template <class GLP>
inline void addInPlace(GLBasePolynomial<GLP> &destination,
                       const GLBasePolynomial<GLP> &term)
{
    for (std::size_t i = 0; i < destination.size(); i++)
        destination[i] += term[i];
}

template <class GLP>
inline void addInPlace(GLPolynomial<GLP> &destination,
                       const GLPolynomial<GLP> &term)
{
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++)
        addInPlace<GLP>(destination[y], term[y]);
}

template <class GLP>
inline void multiplyByIInPlace(GLBasePolynomial<GLP> &poly)
{
    for (std::uint32_t w = 0; w < GLP::phi; w++) {
        for (std::uint32_t x = 0; x < GLP::matrix_dimension; x++) {
            const std::size_t real = baseIndex<GLP>(0, x, w);
            const std::size_t imag = baseIndex<GLP>(1, x, w);
            const typename GLP::T old_real = poly[real];
            poly[real] = typename GLP::T{0} - poly[imag];
            poly[imag] = old_real;
        }
    }
}

template <class GLP>
inline void addWProduct(GLBasePolynomial<GLP> &result,
                        const std::uint32_t gaussian, const std::uint32_t x,
                        const std::uint32_t w_degree,
                        const typename GLP::T value, const bool negative)
{
    if (w_degree < GLP::phi) {
        addSigned(result[baseIndex<GLP>(gaussian, x, w_degree)], value,
                  negative);
        return;
    }
    if (w_degree == GLP::phi) {
        // For prime p, Phi_p(W)=1+W+...+W^(p-1).
        for (std::uint32_t w = 0; w < GLP::phi; w++)
            addSigned(result[baseIndex<GLP>(gaussian, x, w)], value, !negative);
        return;
    }
    // The unreduced product has degree at most 2p-4.  W^p=1 therefore
    // handles every remaining degree in one step.
    addSigned(
        result[baseIndex<GLP>(gaussian, x, w_degree - GLP::cyclotomic_order)],
        value, negative);
}

template <class GLP>
inline void baseMultiply(GLBasePolynomial<GLP> &result,
                         const GLBasePolynomial<GLP> &lhs,
                         const GLBasePolynomial<GLP> &rhs)
{
    clear<GLP>(result);
    for (std::uint32_t wa = 0; wa < GLP::phi; wa++) {
        for (std::uint32_t xa = 0; xa < GLP::matrix_dimension; xa++) {
            for (std::uint32_t ia = 0; ia < 2; ia++) {
                const auto a = lhs[baseIndex<GLP>(ia, xa, wa)];
                if (a == typename GLP::T{0}) continue;
                for (std::uint32_t wb = 0; wb < GLP::phi; wb++) {
                    for (std::uint32_t xb = 0; xb < GLP::matrix_dimension;
                         xb++) {
                        for (std::uint32_t ib = 0; ib < 2; ib++) {
                            const auto b = rhs[baseIndex<GLP>(ib, xb, wb)];
                            if (b == typename GLP::T{0}) continue;

                            std::uint32_t x = xa + xb;
                            std::uint32_t gaussian = ia + ib;
                            bool negative = false;
                            if (x >= GLP::matrix_dimension) {
                                x -= GLP::matrix_dimension;
                                gaussian++;
                            }
                            if (gaussian >= 2) {
                                gaussian -= 2;
                                negative = !negative;
                            }
                            addWProduct<GLP>(result, gaussian, x, wa + wb,
                                             a * b, negative);
                        }
                    }
                }
            }
        }
    }
}

template <class GLP>
inline void polynomialMultiply(GLPolynomial<GLP> &result,
                               const GLPolynomial<GLP> &lhs,
                               const GLPolynomial<GLP> &rhs)
{
    clear<GLP>(result);
    auto product = std::make_unique<GLBasePolynomial<GLP>>();
    for (std::uint32_t ya = 0; ya < GLP::matrix_dimension; ya++) {
        for (std::uint32_t yb = 0; yb < GLP::matrix_dimension; yb++) {
            baseMultiply<GLP>(*product, lhs[ya], rhs[yb]);
            std::uint32_t y = ya + yb;
            if (y >= GLP::matrix_dimension) {
                y -= GLP::matrix_dimension;
                multiplyByIInPlace<GLP>(*product);
            }
            addInPlace<GLP>(result[y], *product);
        }
    }
}

template <class GLP>
inline void inverseWMonomial(
    std::array<std::pair<std::uint32_t, bool>, GLP::phi> &terms,
    std::size_t &term_count, const std::uint32_t exponent)
{
    term_count = 0;
    if (exponent == 0) {
        terms[term_count++] = {0, false};
    }
    else if (exponent == 1) {
        // W^-1 = W^(p-1) = -(1+...+W^(p-2)).
        for (std::uint32_t w = 0; w < GLP::phi; w++)
            terms[term_count++] = {w, true};
    }
    else {
        terms[term_count++] = {GLP::cyclotomic_order - exponent, false};
    }
}

// The paper's trace product a (*) b.  This applies
// (I,X,Y,W)->(-I,Y^-1,Z^-1,W^-1) to b and selects the zero Z coefficient.
template <class GLP>
inline void traceProduct(GLPolynomial<GLP> &result,
                         const GLPolynomial<GLP> &lhs,
                         const GLPolynomial<GLP> &rhs)
{
    clear<GLP>(result);
    std::array<std::pair<std::uint32_t, bool>, GLP::phi> inverse_w{};
    for (std::uint32_t z = 0; z < GLP::matrix_dimension; z++) {
        for (std::uint32_t wa = 0; wa < GLP::phi; wa++) {
            for (std::uint32_t xa = 0; xa < GLP::matrix_dimension; xa++) {
                for (std::uint32_t ia = 0; ia < 2; ia++) {
                    const auto a = lhs[z][baseIndex<GLP>(ia, xa, wa)];
                    if (a == typename GLP::T{0}) continue;
                    for (std::uint32_t wb = 0; wb < GLP::phi; wb++) {
                        std::size_t inverse_w_count = 0;
                        inverseWMonomial<GLP>(inverse_w, inverse_w_count, wb);
                        for (std::uint32_t xb = 0; xb < GLP::matrix_dimension;
                             xb++) {
                            for (std::uint32_t ib = 0; ib < 2; ib++) {
                                const auto b =
                                    rhs[z][baseIndex<GLP>(ib, xb, wb)];
                                if (b == typename GLP::T{0}) continue;

                                std::uint32_t y = 0;
                                std::uint32_t gaussian = ia + ib;
                                bool negative = ib != 0;  // I -> -I
                                if (xb != 0) {
                                    y = GLP::matrix_dimension - xb;
                                    gaussian++;  // Y^-xb = -I*Y^(n-xb)
                                    negative = !negative;
                                }
                                if (gaussian >= 2) {
                                    gaussian -= 2;
                                    negative = !negative;
                                }

                                for (std::size_t term = 0;
                                     term < inverse_w_count; term++) {
                                    const auto [w, inverse_negative] =
                                        inverse_w[term];
                                    addWProduct<GLP>(
                                        result[y], gaussian, xa, wa + w, a * b,
                                        negative != inverse_negative);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <class GLP>
inline GLBasePolynomial<GLP> keyPolynomial(const Key<typename GLP::baseP> &key)
{
    GLBasePolynomial<GLP> result{};
    for (std::uint32_t i = 0; i < GLP::baseP::n; i++) result[i] = key[i];
    return result;
}

template <class GLP>
inline void liftBase(GLPolynomial<GLP> &result,
                     const GLBasePolynomial<GLP> &base)
{
    clear<GLP>(result);
    result[0] = base;
}

// Map s(I,X,W) to s(-I,Y^-1,W^-1), represented in the extended ring.
template <class GLP>
inline void conjugateKeyToY(GLPolynomial<GLP> &result,
                            const GLBasePolynomial<GLP> &key)
{
    clear<GLP>(result);
    std::array<std::pair<std::uint32_t, bool>, GLP::phi> inverse_w{};
    for (std::uint32_t w = 0; w < GLP::phi; w++) {
        std::size_t inverse_w_count = 0;
        inverseWMonomial<GLP>(inverse_w, inverse_w_count, w);
        for (std::uint32_t x = 0; x < GLP::matrix_dimension; x++) {
            for (std::uint32_t gaussian = 0; gaussian < 2; gaussian++) {
                const auto value = key[baseIndex<GLP>(gaussian, x, w)];
                if (value == typename GLP::T{0}) continue;
                std::uint32_t out_gaussian = gaussian;
                std::uint32_t y = 0;
                bool negative = gaussian != 0;
                if (x != 0) {
                    y = GLP::matrix_dimension - x;
                    out_gaussian++;
                    negative = !negative;
                }
                if (out_gaussian >= 2) {
                    out_gaussian -= 2;
                    negative = !negative;
                }
                for (std::size_t term = 0; term < inverse_w_count; term++) {
                    const auto [out_w, inverse_negative] = inverse_w[term];
                    addSigned(result[y][baseIndex<GLP>(out_gaussian, 0, out_w)],
                              value, negative != inverse_negative);
                }
            }
        }
    }
}

template <class GLP>
inline void addReducedWMonomial(GLPolynomial<GLP> &result,
                                const std::uint32_t y,
                                const std::uint32_t gaussian,
                                const std::uint32_t x, const std::uint32_t w,
                                const typename GLP::T value,
                                const bool negative)
{
    if (w < GLP::phi) {
        addSigned(result[y][baseIndex<GLP>(gaussian, x, w)], value, negative);
        return;
    }

    // W^(p-1) = -(1 + W + ... + W^(p-2)) for prime p.
    for (std::uint32_t reduced_w = 0; reduced_w < GLP::phi; reduced_w++)
        addSigned(result[y][baseIndex<GLP>(gaussian, x, reduced_w)], value,
                  !negative);
}

template <class GLP>
inline void polynomialAutomorphism(GLPolynomial<GLP> &result,
                                   const GLPolynomial<GLP> &input,
                                   const std::uint32_t i_multiplier,
                                   const std::uint32_t x_multiplier,
                                   const std::uint32_t y_multiplier,
                                   const std::uint32_t w_multiplier,
                                   const bool swap_xy = false)
{
    constexpr std::uint32_t n = GLP::matrix_dimension;
    constexpr std::uint32_t four_n = 4 * n;
    constexpr std::uint32_t p = GLP::cyclotomic_order;
    clear<GLP>(result);

    for (std::uint32_t y = 0; y < n; y++) {
        const std::uint32_t mapped_y = static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(y_multiplier) * y) % four_n);
        const std::uint32_t y_power = mapped_y / n;
        const std::uint32_t y_remainder = mapped_y % n;
        for (std::uint32_t w = 0; w < GLP::phi; w++) {
            const std::uint32_t mapped_w = static_cast<std::uint32_t>(
                (static_cast<std::uint64_t>(w_multiplier) * w) % p);
            for (std::uint32_t x = 0; x < n; x++) {
                const std::uint32_t mapped_x = static_cast<std::uint32_t>(
                    (static_cast<std::uint64_t>(x_multiplier) * x) % four_n);
                const std::uint32_t x_power = mapped_x / n;
                const std::uint32_t x_remainder = mapped_x % n;
                const std::uint32_t output_x =
                    swap_xy ? y_remainder : x_remainder;
                const std::uint32_t output_y =
                    swap_xy ? x_remainder : y_remainder;

                for (std::uint32_t gaussian = 0; gaussian < 2; gaussian++) {
                    const auto value = input[y][baseIndex<GLP>(gaussian, x, w)];
                    if (value == typename GLP::T{0}) continue;
                    const std::uint32_t gaussian_power =
                        (i_multiplier * gaussian + x_power + y_power) & 3U;
                    addReducedWMonomial<GLP>(
                        result, output_y, gaussian_power & 1U, output_x,
                        mapped_w, value, gaussian_power >= 2);
                }
            }
        }
    }
}

template <class GLP>
inline void baseAutomorphism(GLBasePolynomial<GLP> &result,
                             const GLBasePolynomial<GLP> &input,
                             const std::uint32_t x_multiplier,
                             const std::uint32_t w_multiplier)
{
    GLPolynomial<GLP> lifted;
    GLPolynomial<GLP> transformed;
    liftBase<GLP>(lifted, input);
    polynomialAutomorphism<GLP>(transformed, lifted, x_multiplier & 3U,
                                x_multiplier, 1, w_multiplier);
    result = transformed[0];
}

template <class GLP>
inline std::uint32_t batchExponent(const std::uint32_t batch)
{
    return powMod(GLP::primitive_root, batch, GLP::cyclotomic_order);
}

template <class GLP>
inline std::complex<long double> xRoot(const std::uint32_t slot)
{
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    const std::uint32_t modulus = 4 * GLP::matrix_dimension;
    const std::uint32_t exponent = powMod(5, slot, modulus);
    const long double angle = 2.0L * pi * exponent / modulus;
    return {std::cos(angle), std::sin(angle)};
}

template <class GLP>
inline std::complex<long double> wRoot(const std::uint32_t exponent)
{
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    const long double angle = 2.0L * pi * exponent / GLP::cyclotomic_order;
    return {std::cos(angle), std::sin(angle)};
}

}  // namespace gl_detail

template <class GLP>
inline void GLTraceProduct(GLPolynomial<GLP> &result,
                           const GLPolynomial<GLP> &lhs,
                           const GLPolynomial<GLP> &rhs)
{
    gl_detail::traceProduct<GLP>(result, lhs, rhs);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLEncode(GLPlaintext<GLP, LogQ, LogDelta> &plaintext,
                     const GLMatrixBatch<GLP> &matrices)
{
    using P = typename GLP::baseP;
    static_assert(LogDelta < 63,
                  "the reference GL encoder rounds through int64 precision");

    constexpr std::uint32_t n = GLP::matrix_dimension;
    constexpr std::uint32_t p = GLP::cyclotomic_order;
    constexpr std::uint32_t phi = GLP::phi;
    const long double scale = std::ldexp(1.0L, LogDelta);

    std::array<std::complex<long double>, n> x_roots{};
    for (std::uint32_t j = 0; j < n; j++) x_roots[j] = gl_detail::xRoot<GLP>(j);

    std::vector<std::complex<long double>> w_values(
        static_cast<std::size_t>(n) * n * phi);
    auto w_value = [&](const std::uint32_t row, const std::uint32_t column,
                       const std::uint32_t w) -> std::complex<long double> & {
        return w_values[(static_cast<std::size_t>(row) * n + column) * phi + w];
    };

    // Invert the W canonical embedding.  Values are supplied at all nonzero
    // p-th roots.  Choosing the degree-(p-2) representative determines the
    // missing value at W=1.
    for (std::uint32_t row = 0; row < n; row++) {
        for (std::uint32_t column = 0; column < n; column++) {
            std::complex<long double> at_one = 0;
            for (std::uint32_t batch = 0; batch < phi; batch++) {
                const std::uint32_t exponent =
                    gl_detail::batchExponent<GLP>(batch);
                at_one -= static_cast<std::complex<long double>>(
                              matrices(batch, row, column)) *
                          gl_detail::wRoot<GLP>(exponent);
            }
            for (std::uint32_t w = 0; w < phi; w++) {
                std::complex<long double> coefficient = at_one;
                for (std::uint32_t batch = 0; batch < phi; batch++) {
                    const std::uint32_t exponent =
                        gl_detail::batchExponent<GLP>(batch);
                    const std::uint32_t root_exponent =
                        (p -
                         static_cast<std::uint32_t>(
                             (static_cast<std::uint64_t>(exponent) * w) % p)) %
                        p;
                    coefficient += static_cast<std::complex<long double>>(
                                       matrices(batch, row, column)) *
                                   gl_detail::wRoot<GLP>(root_exponent);
                }
                w_value(row, column, w) =
                    coefficient / static_cast<long double>(p);
            }
        }
    }

    gl_detail::clear<GLP>(plaintext.poly);
    for (std::uint32_t w = 0; w < phi; w++) {
        for (std::uint32_t x = 0; x < n; x++) {
            for (std::uint32_t y = 0; y < n; y++) {
                std::complex<long double> coefficient = 0;
                for (std::uint32_t row = 0; row < n; row++) {
                    const auto x_factor =
                        std::pow(x_roots[row], -static_cast<int>(x));
                    for (std::uint32_t column = 0; column < n; column++) {
                        coefficient +=
                            w_value(row, column, w) * x_factor *
                            std::pow(x_roots[column], -static_cast<int>(y));
                    }
                }
                coefficient *= scale / static_cast<long double>(n * n);
                const auto real =
                    static_cast<__int128_t>(std::round(coefficient.real()));
                const auto imag =
                    static_cast<__int128_t>(std::round(coefficient.imag()));
                plaintext.poly[y][gl_detail::baseIndex<GLP>(0, x, w)] =
                    ckks_detail::signedToLevel<P, LogQ>(real);
                plaintext.poly[y][gl_detail::baseIndex<GLP>(1, x, w)] =
                    ckks_detail::signedToLevel<P, LogQ>(imag);
            }
        }
    }
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLDecode(GLMatrixBatch<GLP> &matrices,
                     const GLPlaintext<GLP, LogQ, LogDelta> &plaintext)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t n = GLP::matrix_dimension;
    constexpr std::uint32_t phi = GLP::phi;
    const long double inverse_scale = std::ldexp(1.0L, -int(LogDelta));

    std::array<std::complex<long double>, n> x_roots{};
    for (std::uint32_t j = 0; j < n; j++) x_roots[j] = gl_detail::xRoot<GLP>(j);

    for (std::uint32_t batch = 0; batch < phi; batch++) {
        const std::uint32_t exponent = gl_detail::batchExponent<GLP>(batch);
        const auto eta = gl_detail::wRoot<GLP>(exponent);
        for (std::uint32_t row = 0; row < n; row++) {
            for (std::uint32_t column = 0; column < n; column++) {
                std::complex<long double> value = 0;
                for (std::uint32_t y = 0; y < n; y++) {
                    const auto y_factor = std::pow(x_roots[column], int(y));
                    for (std::uint32_t x = 0; x < n; x++) {
                        const auto xy_factor =
                            std::pow(x_roots[row], int(x)) * y_factor;
                        for (std::uint32_t w = 0; w < phi; w++) {
                            const long double real =
                                ckks_detail::levelToLongDouble<P, LogQ>(
                                    plaintext.poly[y][gl_detail::baseIndex<GLP>(
                                        0, x, w)]);
                            const long double imag =
                                ckks_detail::levelToLongDouble<P, LogQ>(
                                    plaintext.poly[y][gl_detail::baseIndex<GLP>(
                                        1, x, w)]);
                            value += std::complex<long double>(real, imag) *
                                     xy_factor * std::pow(eta, int(w));
                        }
                    }
                }
                matrices(batch, row, column) =
                    static_cast<std::complex<double>>(value * inverse_scale);
            }
        }
    }
}

namespace gl_detail {

template <class GLP, std::uint32_t LogQ>
inline void encryptBaseAtLevel(GLBaseCiphertextData<GLP> &ciphertext,
                               const GLBasePolynomial<GLP> &message,
                               const Key<typename GLP::baseP> &key,
                               const CKKSNoise noise)
{
    using P = typename GLP::baseP;
    const auto secret = keyPolynomial<GLP>(key);
    for (std::uint32_t i = 0; i < P::n; i++)
        ciphertext[1][i] = ckks_detail::uniformAtLevel<P, LogQ>();

    auto mask_phase = std::make_unique<GLBasePolynomial<GLP>>();
    baseMultiply<GLP>(*mask_phase, ciphertext[1], secret);
    for (std::uint32_t i = 0; i < P::n; i++) {
        ciphertext[0][i] = ckks_detail::reduceToLevel<P, LogQ>(
            message[i] + ckks_detail::sampleNoise<P, LogQ>(noise) -
            (*mask_phase)[i]);
    }
}

template <class GLP, std::uint32_t LogQ>
inline void encryptAtLevel(GLCiphertextData<GLP> &ciphertext,
                           const GLPolynomial<GLP> &message,
                           const Key<typename GLP::baseP> &key,
                           const CKKSNoise noise)
{
    using P = typename GLP::baseP;
    const auto secret = keyPolynomial<GLP>(key);
    auto mask_phase = std::make_unique<GLBasePolynomial<GLP>>();
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
        for (std::uint32_t i = 0; i < P::n; i++)
            ciphertext[1][y][i] = ckks_detail::uniformAtLevel<P, LogQ>();
        baseMultiply<GLP>(*mask_phase, ciphertext[1][y], secret);
        for (std::uint32_t i = 0; i < P::n; i++) {
            ciphertext[0][y][i] = ckks_detail::reduceToLevel<P, LogQ>(
                message[y][i] + ckks_detail::sampleNoise<P, LogQ>(noise) -
                (*mask_phase)[i]);
        }
    }
}

template <class GLP, std::uint32_t LogQ>
inline void phaseAtLevel(GLPolynomial<GLP> &phase,
                         const GLCiphertextData<GLP> &ciphertext,
                         const Key<typename GLP::baseP> &key)
{
    using P = typename GLP::baseP;
    const auto secret = keyPolynomial<GLP>(key);
    auto mask_phase = std::make_unique<GLBasePolynomial<GLP>>();
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
        baseMultiply<GLP>(*mask_phase, ciphertext[1][y], secret);
        for (std::uint32_t i = 0; i < P::n; i++)
            phase[y][i] = ckks_detail::reduceToLevel<P, LogQ>(
                ciphertext[0][y][i] + (*mask_phase)[i]);
    }
}

}  // namespace gl_detail

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLEncrypt(GLCiphertext<GLP, LogQ, LogDelta> &ciphertext,
                      const GLPlaintext<GLP, LogQ, LogDelta> &plaintext,
                      const Key<typename GLP::baseP> &key,
                      const CKKSNoise noise)
{
    gl_detail::encryptAtLevel<GLP, LogQ>(ciphertext.ct, plaintext.poly, key,
                                         noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLEncrypt(GLCiphertext<GLP, LogQ, LogDelta> &ciphertext,
                      const GLPlaintext<GLP, LogQ, LogDelta> &plaintext,
                      const Key<typename GLP::baseP> &key)
{
    GLEncrypt(ciphertext, plaintext, key, GLNoiseAtLevel<LogQ>());
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLEncrypt(GLCiphertext<GLP, LogQ, LogDelta> &ciphertext,
                      const GLMatrixBatch<GLP> &matrices,
                      const Key<typename GLP::baseP> &key,
                      const CKKSNoise noise = GLNoiseAtLevel<LogQ>())
{
    GLPlaintext<GLP, LogQ, LogDelta> plaintext;
    GLEncode(plaintext, matrices);
    GLEncrypt(ciphertext, plaintext, key, noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLDecrypt(GLPlaintext<GLP, LogQ, LogDelta> &plaintext,
                      const GLCiphertext<GLP, LogQ, LogDelta> &ciphertext,
                      const Key<typename GLP::baseP> &key)
{
    gl_detail::phaseAtLevel<GLP, LogQ>(plaintext.poly, ciphertext.ct, key);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLDecrypt(GLMatrixBatch<GLP> &matrices,
                      const GLCiphertext<GLP, LogQ, LogDelta> &ciphertext,
                      const Key<typename GLP::baseP> &key)
{
    GLPlaintext<GLP, LogQ, LogDelta> plaintext;
    GLDecrypt(plaintext, ciphertext, key);
    GLDecode(matrices, plaintext);
}

namespace gl_detail {

template <class GLP, std::uint32_t LogQ, std::uint32_t BaseBit>
inline std::vector<GLPolynomial<GLP>> activeDecompose(
    const GLPolynomial<GLP> &input)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t rows = (LogQ + BaseBit - 1) / BaseBit;
    static_assert(BaseBit > 0);
    static_assert(rows * BaseBit <= std::numeric_limits<typename P::T>::digits,
                  "GL decomposition needs storage for the top signed digit");

    std::vector<GLPolynomial<GLP>> result(rows);
    auto slice_rows = std::make_unique<std::array<Polynomial<P>, rows>>();
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
        ckks_detail::activeBaseDecomposePolynomialRows<P, LogQ, BaseBit, rows>(
            *slice_rows, input[y]);
        for (std::uint32_t row = 0; row < rows; row++)
            result[row][y] = (*slice_rows)[row];
    }
    return result;
}

template <class GLP, std::uint32_t LogQ, std::uint32_t BaseBit>
inline std::vector<GLBasePolynomial<GLP>> activeDecompose(
    const GLBasePolynomial<GLP> &input)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t rows = (LogQ + BaseBit - 1) / BaseBit;
    static_assert(BaseBit > 0);
    static_assert(rows * BaseBit <= std::numeric_limits<typename P::T>::digits,
                  "GL decomposition needs storage for the top signed digit");

    auto result = std::make_unique<std::array<Polynomial<P>, rows>>();
    ckks_detail::activeBaseDecomposePolynomialRows<P, LogQ, BaseBit, rows>(
        *result, input);
    return std::vector<GLBasePolynomial<GLP>>(result->begin(), result->end());
}

template <class GLP, std::uint32_t LogQ, std::uint32_t BaseBit>
inline void activeRecombine(GLPolynomial<GLP> &result,
                            const std::vector<GLPolynomial<GLP>> &rows)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t row_count = (LogQ + BaseBit - 1) / BaseBit;
    if (rows.size() != row_count)
        throw std::invalid_argument("invalid GL decomposition row count");

    clear<GLP>(result);
    for (std::uint32_t row = 0; row < row_count; row++) {
        const std::uint32_t shift = (row_count - row - 1) * BaseBit;
        for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++)
            for (std::uint32_t i = 0; i < P::n; i++)
                result[y][i] += rows[row][y][i] << shift;
    }
    reduce<GLP, LogQ>(result);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t BaseBit>
inline void activeRecombine(GLBasePolynomial<GLP> &result,
                            const std::vector<GLBasePolynomial<GLP>> &rows)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t row_count = (LogQ + BaseBit - 1) / BaseBit;
    if (rows.size() != row_count)
        throw std::invalid_argument("invalid GL decomposition row count");

    clear<GLP>(result);
    for (std::uint32_t row = 0; row < row_count; row++) {
        const std::uint32_t shift = (row_count - row - 1) * BaseBit;
        for (std::uint32_t i = 0; i < P::n; i++)
            result[i] += rows[row][i] << shift;
    }
    reduce<GLP, LogQ>(result);
}

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t RhsLogQ,
          std::uint32_t LogScale, std::uint32_t BbarBit, bool UseMatrixTrace,
          bool NormalizeMatrixProduct>
inline void productRescaleDD(GLPolynomial<GLP> &result,
                             const GLPolynomial<GLP> &lhs,
                             const GLPolynomial<GLP> &rhs)
{
    using P = typename GLP::baseP;
    using T = typename P::T;
    constexpr std::uint32_t lhs_rows = (LhsLogQ + BbarBit - 1) / BbarBit;
    constexpr std::uint32_t rhs_rows = (RhsLogQ + BbarBit - 1) / BbarBit;
    constexpr std::uint32_t base_log_q = LhsLogQ < RhsLogQ ? LhsLogQ : RhsLogQ;
    constexpr std::uint32_t out_log_q = base_log_q - LogScale;
    static_assert(LogScale < base_log_q);
    static_assert(lhs_rows * BbarBit <= std::numeric_limits<T>::digits);
    static_assert(rhs_rows * BbarBit <= std::numeric_limits<T>::digits);
    static_assert(2 * BbarBit + 2 < 63,
                  "GL DD digit products must fit in int64_t");

    auto lhs_digits = activeDecompose<GLP, LhsLogQ, BbarBit>(lhs);
    auto rhs_digits = activeDecompose<GLP, RhsLogQ, BbarBit>(rhs);
    const std::size_t coefficient_count =
        static_cast<std::size_t>(GLP::matrix_dimension) * P::n;
    auto digit_product = std::make_unique<GLPolynomial<GLP>>();

    auto accumulate_products = [&](auto &accumulators, auto add_shifted) {
        for (std::uint32_t i = 0; i < lhs_rows; i++) {
            const int lhs_shift =
                static_cast<int>((lhs_rows - i - 1) * BbarBit);
            for (std::uint32_t j = 0; j < rhs_rows; j++) {
                const int rhs_shift =
                    static_cast<int>((rhs_rows - j - 1) * BbarBit);
                if constexpr (UseMatrixTrace)
                    traceProduct<GLP>(*digit_product, lhs_digits[i],
                                      rhs_digits[j]);
                else
                    polynomialMultiply<GLP>(*digit_product, lhs_digits[i],
                                            rhs_digits[j]);
                for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
                    for (std::uint32_t coefficient = 0; coefficient < P::n;
                         coefficient++) {
                        std::int64_t value = torusToSignedSmall<P>(
                            (*digit_product)[y][coefficient]);
                        if constexpr (NormalizeMatrixProduct) {
                            const __int128_t scaled =
                                static_cast<__int128_t>(value) *
                                GLP::matrix_dimension;
                            if (scaled <
                                    std::numeric_limits<std::int64_t>::min() ||
                                scaled >
                                    std::numeric_limits<std::int64_t>::max())
                                throw std::overflow_error(
                                    "normalized GL digit product exceeds "
                                    "int64_t");
                            value = static_cast<std::int64_t>(scaled);
                        }
                        const std::size_t index =
                            static_cast<std::size_t>(y) * P::n + coefficient;
                        add_shifted(accumulators[index], value,
                                    lhs_shift + rhs_shift);
                    }
                }
            }
        }
    };

    if constexpr (std::is_same_v<T, __uint128_t>) {
        std::vector<Wide384> accumulators(coefficient_count);
        accumulate_products(
            accumulators, [](Wide384 &accumulator, const std::int64_t value,
                             const int shift) {
                accumulator.add_shifted(static_cast<__int128_t>(value), shift);
            });
        for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
            for (std::uint32_t coefficient = 0; coefficient < P::n;
                 coefficient++) {
                const std::size_t index =
                    static_cast<std::size_t>(y) * P::n + coefficient;
                result[y][coefficient] =
                    ckks_detail::reduceToLevel<P, out_log_q>(
                        ckks_detail::rescaleAccumulator<P, LogScale>(
                            accumulators[index]));
            }
        }
    }
    else {
        static_assert(is_multilimb_uint_v<T>);
        using Accumulator = WideSignedLimbAccumulator<2 * T::limbs + 2>;
        std::vector<Accumulator> accumulators(coefficient_count);
        accumulate_products(
            accumulators,
            [](Accumulator &accumulator, const std::int64_t value,
               const int shift) { accumulator.add_shifted_i64(value, shift); });
        const T divisor = LogScale == 0 ? T{1} : (T{1} << LogScale);
        for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
            for (std::uint32_t coefficient = 0; coefficient < P::n;
                 coefficient++) {
                const std::size_t index =
                    static_cast<std::size_t>(y) * P::n + coefficient;
                if constexpr (LogScale > 0) {
                    accumulators[index].add_shifted_i64(
                        accumulators[index].is_negative() ? -1 : 1,
                        static_cast<int>(LogScale) - 1);
                }
                result[y][coefficient] =
                    ckks_detail::reduceToLevel<P, out_log_q>(
                        accumulators[index].template div_to_torus<T::limbs>(
                            divisor));
            }
        }
    }
}

}  // namespace gl_detail

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
struct GLDDBigKeySwitchKey {
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t primary_bit = PrimaryBit;
    static constexpr std::uint32_t bbar_bit = BbarBit;
    static constexpr std::uint32_t auxiliary_log_q = GLP::auxiliary_log_q;
    static constexpr std::uint32_t key_log_q = LogQ + auxiliary_log_q;
    static constexpr std::uint32_t primary_rows =
        (LogQ + PrimaryBit - 1) / PrimaryBit;
    static constexpr std::uint32_t bbar_rows =
        (key_log_q + BbarBit - 1) / BbarBit;

    static_assert(key_log_q <= std::numeric_limits<typename GLP::T>::digits,
                  "GL q*q0 key-switch modulus exceeds torus storage");
    static_assert(primary_rows * PrimaryBit <=
                  std::numeric_limits<typename GLP::T>::digits);
    static_assert(bbar_rows * BbarBit <=
                  std::numeric_limits<typename GLP::T>::digits);

    std::vector<GLCiphertextData<GLP>> data{};

    void allocate() { data.resize(primary_rows * bbar_rows); }

    GLCiphertextData<GLP> &at(const std::uint32_t primary,
                              const std::uint32_t bbar)
    {
        return data.at(static_cast<std::size_t>(primary) * bbar_rows + bbar);
    }
    const GLCiphertextData<GLP> &at(const std::uint32_t primary,
                                    const std::uint32_t bbar) const
    {
        return data.at(static_cast<std::size_t>(primary) * bbar_rows + bbar);
    }
};

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
struct GLDDSmallKeySwitchKey {
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t primary_bit = PrimaryBit;
    static constexpr std::uint32_t bbar_bit = BbarBit;
    static constexpr std::uint32_t auxiliary_log_q = GLP::auxiliary_log_q;
    static constexpr std::uint32_t key_log_q = LogQ + auxiliary_log_q;
    static constexpr std::uint32_t primary_rows =
        (LogQ + PrimaryBit - 1) / PrimaryBit;
    static constexpr std::uint32_t bbar_rows =
        (key_log_q + BbarBit - 1) / BbarBit;

    static_assert(key_log_q <= std::numeric_limits<typename GLP::T>::digits,
                  "GL q*q0 key-switch modulus exceeds torus storage");
    static_assert(primary_rows * PrimaryBit <=
                  std::numeric_limits<typename GLP::T>::digits);
    static_assert(bbar_rows * BbarBit <=
                  std::numeric_limits<typename GLP::T>::digits);

    std::vector<GLBaseCiphertextData<GLP>> data{};

    void allocate() { data.resize(primary_rows * bbar_rows); }

    GLBaseCiphertextData<GLP> &at(const std::uint32_t primary,
                                  const std::uint32_t bbar)
    {
        return data.at(static_cast<std::size_t>(primary) * bbar_rows + bbar);
    }
    const GLBaseCiphertextData<GLP> &at(const std::uint32_t primary,
                                        const std::uint32_t bbar) const
    {
        return data.at(static_cast<std::size_t>(primary) * bbar_rows + bbar);
    }
};

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDBigKeySwitchKeyGen(
    GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key,
    const GLPolynomial<GLP> &source_secret,
    const Key<typename GLP::baseP> &destination_key,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    using P = typename GLP::baseP;
    using SwitchKey = GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;

    switch_key.allocate();
    auto gadget = std::make_unique<GLPolynomial<GLP>>();
    auto ordinary = std::make_unique<GLCiphertextData<GLP>>();
    for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
         primary++) {
        const std::uint32_t shift =
            (SwitchKey::primary_rows - primary - 1) * PrimaryBit +
            SwitchKey::auxiliary_log_q;
        for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++)
            for (std::uint32_t i = 0; i < P::n; i++)
                (*gadget)[y][i] =
                    ckks_detail::reduceToLevel<P, SwitchKey::key_log_q>(
                        source_secret[y][i] << shift);
        gl_detail::encryptAtLevel<GLP, SwitchKey::key_log_q>(
            *ordinary, *gadget, destination_key, noise);
        for (std::size_t component = 0; component < 2; component++) {
            auto rows =
                gl_detail::activeDecompose<GLP, SwitchKey::key_log_q, BbarBit>(
                    (*ordinary)[component]);
            for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++)
                switch_key.at(primary, bbar)[component] = std::move(rows[bbar]);
        }
    }
}

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDSmallKeySwitchKeyGen(
    GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key,
    const GLBasePolynomial<GLP> &source_secret,
    const Key<typename GLP::baseP> &destination_key,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    using P = typename GLP::baseP;
    using SwitchKey = GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;
    switch_key.allocate();
    GLBasePolynomial<GLP> gadget{};
    GLBaseCiphertextData<GLP> ordinary{};
    for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
         primary++) {
        const std::uint32_t shift =
            (SwitchKey::primary_rows - primary - 1) * PrimaryBit +
            SwitchKey::auxiliary_log_q;
        for (std::uint32_t i = 0; i < P::n; i++)
            gadget[i] = ckks_detail::reduceToLevel<P, SwitchKey::key_log_q>(
                source_secret[i] << shift);
        gl_detail::encryptBaseAtLevel<GLP, SwitchKey::key_log_q>(
            ordinary, gadget, destination_key, noise);
        for (std::size_t component = 0; component < 2; component++) {
            auto rows =
                gl_detail::activeDecompose<GLP, SwitchKey::key_log_q, BbarBit>(
                    ordinary[component]);
            for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++)
                switch_key.at(primary, bbar)[component] = std::move(rows[bbar]);
        }
    }
}

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDBigKeySwitch(
    GLCiphertextData<GLP> &result, const GLPolynomial<GLP> &input,
    const GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key)
{
    using SwitchKey = GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;
    if (switch_key.data.size() !=
        SwitchKey::primary_rows * SwitchKey::bbar_rows)
        throw std::invalid_argument("uninitialized GL big key-switch key");

    auto input_digits =
        gl_detail::activeDecompose<GLP, LogQ, PrimaryBit>(input);
    std::array<std::vector<GLPolynomial<GLP>>, 2> digit_rows{
        std::vector<GLPolynomial<GLP>>(SwitchKey::bbar_rows),
        std::vector<GLPolynomial<GLP>>(SwitchKey::bbar_rows)};
    auto product = std::make_unique<GLPolynomial<GLP>>();
    for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
         primary++) {
        for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++) {
            for (std::size_t component = 0; component < 2; component++) {
                gl_detail::polynomialMultiply<GLP>(
                    *product, input_digits[primary],
                    switch_key.at(primary, bbar)[component]);
                gl_detail::addInPlace<GLP>(digit_rows[component][bbar],
                                           *product);
            }
        }
    }
    for (std::size_t component = 0; component < 2; component++)
        gl_detail::activeRecombine<GLP, SwitchKey::key_log_q, BbarBit>(
            result[component], digit_rows[component]);
    gl_detail::divideRoundLevel<GLP, SwitchKey::key_log_q,
                                SwitchKey::auxiliary_log_q>(result[0]);
    gl_detail::divideRoundLevel<GLP, SwitchKey::key_log_q,
                                SwitchKey::auxiliary_log_q>(result[1]);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDSmallKeySwitch(
    GLCiphertextData<GLP> &result, const GLPolynomial<GLP> &input,
    const GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key)
{
    using SwitchKey = GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;
    if (switch_key.data.size() !=
        SwitchKey::primary_rows * SwitchKey::bbar_rows)
        throw std::invalid_argument("uninitialized GL small key-switch key");

    auto input_digits =
        gl_detail::activeDecompose<GLP, LogQ, PrimaryBit>(input);
    std::array<std::vector<GLPolynomial<GLP>>, 2> digit_rows{
        std::vector<GLPolynomial<GLP>>(SwitchKey::bbar_rows),
        std::vector<GLPolynomial<GLP>>(SwitchKey::bbar_rows)};
    auto product = std::make_unique<GLBasePolynomial<GLP>>();
    for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
         primary++) {
        for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++) {
            for (std::size_t component = 0; component < 2; component++) {
                for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
                    gl_detail::baseMultiply<GLP>(
                        *product, input_digits[primary][y],
                        switch_key.at(primary, bbar)[component]);
                    gl_detail::addInPlace<GLP>(digit_rows[component][bbar][y],
                                               *product);
                }
            }
        }
    }
    for (std::size_t component = 0; component < 2; component++)
        gl_detail::activeRecombine<GLP, SwitchKey::key_log_q, BbarBit>(
            result[component], digit_rows[component]);
    gl_detail::divideRoundLevel<GLP, SwitchKey::key_log_q,
                                SwitchKey::auxiliary_log_q>(result[0]);
    gl_detail::divideRoundLevel<GLP, SwitchKey::key_log_q,
                                SwitchKey::auxiliary_log_q>(result[1]);
}

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
struct GLMatrixRelinKey {
    static constexpr std::uint32_t log_q = LogQ;
    GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> conjugate_key{};
    GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> product_key{};
};

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLMatrixRelinKeyGen(
    GLMatrixRelinKey<GLP, LogQ, PrimaryBit, BbarBit> &relin_key,
    const Key<typename GLP::baseP> &key,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    auto base_secret = gl_detail::keyPolynomial<GLP>(key);
    auto conjugate_secret = std::make_unique<GLPolynomial<GLP>>();
    auto lifted_secret = std::make_unique<GLPolynomial<GLP>>();
    auto product_secret = std::make_unique<GLPolynomial<GLP>>();
    gl_detail::conjugateKeyToY<GLP>(*conjugate_secret, base_secret);
    gl_detail::liftBase<GLP>(*lifted_secret, base_secret);
    gl_detail::polynomialMultiply<GLP>(*product_secret, *lifted_secret,
                                       *conjugate_secret);

    GLDDBigKeySwitchKeyGen(relin_key.conjugate_key, *conjugate_secret, key,
                           noise);
    GLDDBigKeySwitchKeyGen(relin_key.product_key, *product_secret, key, noise);
}

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
struct GLMatrixMultTraits {
    static constexpr std::uint32_t log_scale =
        LhsLogDelta > RhsLogDelta ? LhsLogDelta : RhsLogDelta;
    static constexpr std::uint32_t base_log_q =
        LhsLogQ < RhsLogQ ? LhsLogQ : RhsLogQ;
    static_assert(base_log_q > log_scale);

    static constexpr std::uint32_t log_q = base_log_q - log_scale;
    static constexpr std::uint32_t log_delta =
        LhsLogDelta < RhsLogDelta ? LhsLogDelta : RhsLogDelta;
    using Ciphertext = GLCiphertext<GLP, log_q, log_delta>;
};

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
using GLMatrixMultResult =
    typename GLMatrixMultTraits<GLP, LhsLogQ, LhsLogDelta, RhsLogQ,
                                RhsLogDelta>::Ciphertext;

// Computes lhs * rhs^* slotwise across the matrix batch.  As in the paper,
// the ciphertext must be the left operand so that its secret can be taken out
// of the trace.  NormalizeMatrixProduct cancels the trace's factor 1/n.
template <bool NormalizeMatrixProduct = true, class GLP, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta, std::uint32_t BbarBit = GLP::baseP::B̅gbit>
inline void GLPlaintextMatrixMultiply(
    GLMatrixMultResult<GLP, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta> &result,
    const GLCiphertext<GLP, LhsLogQ, LhsLogDelta> &lhs,
    const GLPlaintext<GLP, RhsLogQ, RhsLogDelta> &rhs)
{
    using Traits =
        GLMatrixMultTraits<GLP, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta>;
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, true, NormalizeMatrixProduct>(
        result[0], lhs[0], rhs.poly);
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, true, NormalizeMatrixProduct>(
        result[1], lhs[1], rhs.poly);
}

template <bool NormalizeMatrixProduct = true, class GLP, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLMatrixMultiply(
    GLMatrixMultResult<GLP, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta> &result,
    const GLCiphertext<GLP, LhsLogQ, LhsLogDelta> &lhs,
    const GLCiphertext<GLP, RhsLogQ, RhsLogDelta> &rhs,
    const GLMatrixRelinKey<GLP,
                           GLMatrixMultTraits<GLP, LhsLogQ, LhsLogDelta,
                                              RhsLogQ, RhsLogDelta>::log_q,
                           PrimaryBit, BbarBit> &relin_key)
{
    using Traits =
        GLMatrixMultTraits<GLP, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta>;
    std::array<GLPolynomial<GLP>, 4> tensor{};
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, true, NormalizeMatrixProduct>(
        tensor[0], lhs[0], rhs[0]);
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, true, NormalizeMatrixProduct>(
        tensor[1], lhs[1], rhs[0]);
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, true, NormalizeMatrixProduct>(
        tensor[2], lhs[0], rhs[1]);
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, true, NormalizeMatrixProduct>(
        tensor[3], lhs[1], rhs[1]);

    auto conjugate_term = std::make_unique<GLCiphertextData<GLP>>();
    auto product_term = std::make_unique<GLCiphertextData<GLP>>();
    GLDDBigKeySwitch(*conjugate_term, tensor[2], relin_key.conjugate_key);
    GLDDBigKeySwitch(*product_term, tensor[3], relin_key.product_key);

    result[0] = std::move(tensor[0]);
    result[1] = std::move(tensor[1]);
    gl_detail::addInPlace<GLP>(result[0], (*conjugate_term)[0]);
    gl_detail::addInPlace<GLP>(result[1], (*conjugate_term)[1]);
    gl_detail::addInPlace<GLP>(result[0], (*product_term)[0]);
    gl_detail::addInPlace<GLP>(result[1], (*product_term)[1]);
    gl_detail::reduce<GLP, Traits::log_q>(result[0]);
    gl_detail::reduce<GLP, Traits::log_q>(result[1]);
}

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
using GLHadamardRelinKey =
    GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLHadamardRelinKeyGen(
    GLHadamardRelinKey<GLP, LogQ, PrimaryBit, BbarBit> &relin_key,
    const Key<typename GLP::baseP> &key,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    const auto secret = gl_detail::keyPolynomial<GLP>(key);
    GLBasePolynomial<GLP> secret_square{};
    gl_detail::baseMultiply<GLP>(secret_square, secret, secret);
    GLDDSmallKeySwitchKeyGen(relin_key, secret_square, key, noise);
}

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta,
          std::uint32_t PrimaryBit, std::uint32_t BbarBit>
inline void GLHadamardMultiply(
    GLMatrixMultResult<GLP, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta> &result,
    const GLCiphertext<GLP, LhsLogQ, LhsLogDelta> &lhs,
    const GLCiphertext<GLP, RhsLogQ, RhsLogDelta> &rhs,
    const GLHadamardRelinKey<GLP,
                             GLMatrixMultTraits<GLP, LhsLogQ, LhsLogDelta,
                                                RhsLogQ, RhsLogDelta>::log_q,
                             PrimaryBit, BbarBit> &relin_key)
{
    using Traits =
        GLMatrixMultTraits<GLP, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta>;
    std::array<GLPolynomial<GLP>, 4> tensor{};
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, false, false>(tensor[0], lhs[0],
                                                       rhs[0]);
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, false, false>(tensor[1], lhs[0],
                                                       rhs[1]);
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, false, false>(tensor[2], lhs[1],
                                                       rhs[0]);
    gl_detail::productRescaleDD<GLP, LhsLogQ, RhsLogQ, Traits::log_scale,
                                BbarBit, false, false>(tensor[3], lhs[1],
                                                       rhs[1]);
    gl_detail::addInPlace<GLP>(tensor[1], tensor[2]);
    gl_detail::reduce<GLP, Traits::log_q>(tensor[1]);

    auto square_term = std::make_unique<GLCiphertextData<GLP>>();
    GLDDSmallKeySwitch(*square_term, tensor[3], relin_key);
    result[0] = std::move(tensor[0]);
    result[1] = std::move(tensor[1]);
    gl_detail::addInPlace<GLP>(result[0], (*square_term)[0]);
    gl_detail::addInPlace<GLP>(result[1], (*square_term)[1]);
    gl_detail::reduce<GLP, Traits::log_q>(result[0]);
    gl_detail::reduce<GLP, Traits::log_q>(result[1]);
}

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
using GLConjugationKey = GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLConjugationKeyGen(
    GLConjugationKey<GLP, LogQ, PrimaryBit, BbarBit> &conjugation_key,
    const Key<typename GLP::baseP> &key,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    constexpr std::uint32_t inverse_x = 4 * GLP::matrix_dimension - 1;
    constexpr std::uint32_t inverse_w = GLP::cyclotomic_order - 1;
    const auto secret = gl_detail::keyPolynomial<GLP>(key);
    GLBasePolynomial<GLP> source_secret{};
    gl_detail::baseAutomorphism<GLP>(source_secret, secret, inverse_x,
                                     inverse_w);
    GLDDSmallKeySwitchKeyGen(conjugation_key, source_secret, key, noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PrimaryBit, std::uint32_t BbarBit>
inline void GLConjugate(
    GLCiphertext<GLP, LogQ, LogDelta> &result,
    const GLCiphertext<GLP, LogQ, LogDelta> &input,
    const GLConjugationKey<GLP, LogQ, PrimaryBit, BbarBit> &conjugation_key)
{
    constexpr std::uint32_t inverse_x = 4 * GLP::matrix_dimension - 1;
    constexpr std::uint32_t inverse_w = GLP::cyclotomic_order - 1;
    GLPolynomial<GLP> body;
    GLPolynomial<GLP> mask;
    gl_detail::polynomialAutomorphism<GLP>(body, input[0], 3, inverse_x,
                                           inverse_x, inverse_w);
    gl_detail::polynomialAutomorphism<GLP>(mask, input[1], 3, inverse_x,
                                           inverse_x, inverse_w);
    GLDDSmallKeySwitch(result.ct, mask, conjugation_key);
    gl_detail::addInPlace<GLP>(result[0], body);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
struct GLRowRotationKey {
    std::uint32_t amount = 0;
    std::uint32_t multiplier = 1;
    GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> switch_key{};
};

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLRowRotationKeyGen(
    GLRowRotationKey<GLP, LogQ, PrimaryBit, BbarBit> &rotation_key,
    const Key<typename GLP::baseP> &key, const std::uint32_t amount,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    rotation_key.amount = amount % GLP::matrix_dimension;
    rotation_key.multiplier =
        gl_detail::powMod(5, rotation_key.amount, 4 * GLP::matrix_dimension);
    const auto secret = gl_detail::keyPolynomial<GLP>(key);
    GLBasePolynomial<GLP> source_secret{};
    gl_detail::baseAutomorphism<GLP>(source_secret, secret,
                                     rotation_key.multiplier, 1);
    GLDDSmallKeySwitchKeyGen(rotation_key.switch_key, source_secret, key,
                             noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PrimaryBit, std::uint32_t BbarBit>
inline void GLRotateRows(
    GLCiphertext<GLP, LogQ, LogDelta> &result,
    const GLCiphertext<GLP, LogQ, LogDelta> &input,
    const GLRowRotationKey<GLP, LogQ, PrimaryBit, BbarBit> &rotation_key)
{
    GLPolynomial<GLP> body;
    GLPolynomial<GLP> mask;
    gl_detail::polynomialAutomorphism<GLP>(body, input[0], 1,
                                           rotation_key.multiplier, 1, 1);
    gl_detail::polynomialAutomorphism<GLP>(mask, input[1], 1,
                                           rotation_key.multiplier, 1, 1);
    GLDDSmallKeySwitch(result.ct, mask, rotation_key.switch_key);
    gl_detail::addInPlace<GLP>(result[0], body);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLRotateColumns(GLCiphertext<GLP, LogQ, LogDelta> &result,
                            const GLCiphertext<GLP, LogQ, LogDelta> &input,
                            const std::uint32_t amount)
{
    const std::uint32_t multiplier = gl_detail::powMod(
        5, amount % GLP::matrix_dimension, 4 * GLP::matrix_dimension);
    gl_detail::polynomialAutomorphism<GLP>(result[0], input[0], 1, 1,
                                           multiplier, 1);
    gl_detail::polynomialAutomorphism<GLP>(result[1], input[1], 1, 1,
                                           multiplier, 1);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
struct GLBatchRotationKey {
    std::uint32_t amount = 0;
    std::uint32_t multiplier = 1;
    GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> switch_key{};
};

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLBatchRotationKeyGen(
    GLBatchRotationKey<GLP, LogQ, PrimaryBit, BbarBit> &rotation_key,
    const Key<typename GLP::baseP> &key, const std::uint32_t amount,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    rotation_key.amount = amount % GLP::phi;
    rotation_key.multiplier = gl_detail::powMod(
        GLP::primitive_root, rotation_key.amount, GLP::cyclotomic_order);
    const auto secret = gl_detail::keyPolynomial<GLP>(key);
    GLBasePolynomial<GLP> source_secret{};
    gl_detail::baseAutomorphism<GLP>(source_secret, secret, 1,
                                     rotation_key.multiplier);
    GLDDSmallKeySwitchKeyGen(rotation_key.switch_key, source_secret, key,
                             noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PrimaryBit, std::uint32_t BbarBit>
inline void GLRotateBatches(
    GLCiphertext<GLP, LogQ, LogDelta> &result,
    const GLCiphertext<GLP, LogQ, LogDelta> &input,
    const GLBatchRotationKey<GLP, LogQ, PrimaryBit, BbarBit> &rotation_key)
{
    GLPolynomial<GLP> body;
    GLPolynomial<GLP> mask;
    gl_detail::polynomialAutomorphism<GLP>(body, input[0], 1, 1, 1,
                                           rotation_key.multiplier);
    gl_detail::polynomialAutomorphism<GLP>(mask, input[1], 1, 1, 1,
                                           rotation_key.multiplier);
    GLDDSmallKeySwitch(result.ct, mask, rotation_key.switch_key);
    gl_detail::addInPlace<GLP>(result[0], body);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
using GLTransposeKey = GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLTransposeKeyGen(
    GLTransposeKey<GLP, LogQ, PrimaryBit, BbarBit> &transpose_key,
    const Key<typename GLP::baseP> &key,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    const auto secret = gl_detail::keyPolynomial<GLP>(key);
    GLPolynomial<GLP> lifted;
    GLPolynomial<GLP> source_secret;
    gl_detail::liftBase<GLP>(lifted, secret);
    gl_detail::polynomialAutomorphism<GLP>(source_secret, lifted, 1, 1, 1, 1,
                                           true);
    GLDDBigKeySwitchKeyGen(transpose_key, source_secret, key, noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PrimaryBit, std::uint32_t BbarBit>
inline void GLTranspose(
    GLCiphertext<GLP, LogQ, LogDelta> &result,
    const GLCiphertext<GLP, LogQ, LogDelta> &input,
    const GLTransposeKey<GLP, LogQ, PrimaryBit, BbarBit> &transpose_key)
{
    GLPolynomial<GLP> body;
    GLPolynomial<GLP> mask;
    gl_detail::polynomialAutomorphism<GLP>(body, input[0], 1, 1, 1, 1, true);
    gl_detail::polynomialAutomorphism<GLP>(mask, input[1], 1, 1, 1, 1, true);
    GLDDBigKeySwitch(result.ct, mask, transpose_key);
    gl_detail::addInPlace<GLP>(result[0], body);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
using GLConjugateTransposeKey =
    GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLConjugateTransposeKeyGen(
    GLConjugateTransposeKey<GLP, LogQ, PrimaryBit, BbarBit> &transpose_key,
    const Key<typename GLP::baseP> &key,
    const CKKSNoise noise = GLNoiseAtLevel<LogQ + GLP::auxiliary_log_q>())
{
    constexpr std::uint32_t inverse_x = 4 * GLP::matrix_dimension - 1;
    constexpr std::uint32_t inverse_w = GLP::cyclotomic_order - 1;
    const auto secret = gl_detail::keyPolynomial<GLP>(key);
    GLPolynomial<GLP> lifted;
    GLPolynomial<GLP> source_secret;
    gl_detail::liftBase<GLP>(lifted, secret);
    gl_detail::polynomialAutomorphism<GLP>(source_secret, lifted, 3, inverse_x,
                                           inverse_x, inverse_w, true);
    GLDDBigKeySwitchKeyGen(transpose_key, source_secret, key, noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PrimaryBit, std::uint32_t BbarBit>
inline void GLConjugateTranspose(
    GLCiphertext<GLP, LogQ, LogDelta> &result,
    const GLCiphertext<GLP, LogQ, LogDelta> &input,
    const GLConjugateTransposeKey<GLP, LogQ, PrimaryBit, BbarBit>
        &transpose_key)
{
    constexpr std::uint32_t inverse_x = 4 * GLP::matrix_dimension - 1;
    constexpr std::uint32_t inverse_w = GLP::cyclotomic_order - 1;
    GLPolynomial<GLP> body;
    GLPolynomial<GLP> mask;
    gl_detail::polynomialAutomorphism<GLP>(body, input[0], 3, inverse_x,
                                           inverse_x, inverse_w, true);
    gl_detail::polynomialAutomorphism<GLP>(mask, input[1], 3, inverse_x,
                                           inverse_x, inverse_w, true);
    GLDDBigKeySwitch(result.ct, mask, transpose_key);
    gl_detail::addInPlace<GLP>(result[0], body);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLAdd(GLCiphertext<GLP, LogQ, LogDelta> &result,
                  const GLCiphertext<GLP, LogQ, LogDelta> &lhs,
                  const GLCiphertext<GLP, LogQ, LogDelta> &rhs)
{
    result = lhs;
    gl_detail::addInPlace<GLP>(result[0], rhs[0]);
    gl_detail::addInPlace<GLP>(result[1], rhs[1]);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

}  // namespace TFHEpp

#endif  // TFHEPP_DEFAULT_128BIT_PARAMS
