#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ckks/ckks.hpp"
#include "modular_ntt.hpp"
#include "params.hpp"
#include "tfhe/key.hpp"

// The approximate GL implementation uses the active-level and multi-limb
// arithmetic introduced by the CKKS implementation.  Those facilities are
// available with TFHEpp's default 128-bit parameter family.
#ifdef TFHEPP_DEFAULT_128BIT_PARAMS

namespace TFHEpp {

namespace gl_detail {

inline bool inOpenMPParallelRegion()
{
#ifdef _OPENMP
    return omp_in_parallel() != 0;
#else
    return false;
#endif
}

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

// Multi-limb storage profiles for the SHIP parameter sets reported in
// ePrint 2026/811.  In particular, n256p17 needs log(PQ)=214 bits, so the
// 128-bit lvl4param torus is not large enough even though it has the right
// polynomial degree.  The four-limb replacement retains degree 8192 and gives
// enough headroom for the 180-bit ciphertext chain plus the 34-bit switching
// modulus.  lvl5param and lvl6param already have the required degrees and
// storage widths for n512p17 and n1024p17.
struct GL256p17BaseParameter {
    static constexpr std::int32_t key_value_max = 1;
    static constexpr std::int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 13;
    static constexpr std::uint32_t n = 1U << nbit;
    static constexpr std::uint32_t k = 1;
    using T = MultiLimbUInt<4>;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t B̅gbit = 8;
};

using GL256p17Parameter = GLParameter<GL256p17BaseParameter, 256, 17, 34>;
using GL512p17Parameter = GLParameter<lvl5param, 512, 17, 92>;
using GL1024p17Parameter = GLParameter<lvl6param, 1024, 17, 220>;

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

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(coefficients_);
    }

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

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(data);
    }
};

template <class GLP>
using GLBaseCiphertextData = std::array<GLBasePolynomial<GLP>, 2>;

// DD evaluation-key rows are centered auxiliary-base digits.  Keeping those
// digits in the full torus type multiplies their payload by as much as 112x
// for the GL paper profiles.  Store the already-decomposed rows in the
// smallest signed native type and expand only the row currently consumed by
// a key switch.  This is an unseeded representation: both RLWE components are
// retained, and no evaluation-key ciphertext is decomposed online.
template <std::uint32_t Bits>
struct GLSignedDigitStorage {
    static_assert(Bits > 0 && Bits <= 30,
                  "packed GL DD digits must fit a signed 32-bit word");
    using type = std::conditional_t<
        (Bits <= 8), std::int8_t,
        std::conditional_t<(Bits <= 16), std::int16_t, std::int32_t>>;
};

template <std::uint32_t Bits>
using GLSignedDigit = typename GLSignedDigitStorage<Bits>::type;

template <class GLP, std::uint32_t Bits>
using GLPackedBasePolynomial = std::array<GLSignedDigit<Bits>, GLP::baseP::n>;

template <class GLP, std::uint32_t Bits>
class GLPackedPolynomial {
public:
    using BasePolynomial = GLPackedBasePolynomial<GLP, Bits>;

    GLPackedPolynomial() : coefficients_(GLP::matrix_dimension) {}

    BasePolynomial &operator[](const std::size_t y) { return coefficients_[y]; }
    const BasePolynomial &operator[](const std::size_t y) const
    {
        return coefficients_[y];
    }

    constexpr std::size_t size() const { return GLP::matrix_dimension; }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(coefficients_);
    }

private:
    std::vector<BasePolynomial> coefficients_;
};

template <class GLP, std::uint32_t Bits>
struct GLPackedCiphertextData {
    std::array<GLPackedPolynomial<GLP, Bits>, 2> data{};

    GLPackedPolynomial<GLP, Bits> &operator[](const std::size_t component)
    {
        return data[component];
    }
    const GLPackedPolynomial<GLP, Bits> &operator[](
        const std::size_t component) const
    {
        return data[component];
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(data);
    }
};

template <class GLP, std::uint32_t Bits>
using GLPackedBaseCiphertextData =
    std::array<GLPackedBasePolynomial<GLP, Bits>, 2>;

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

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(ct);
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
inline void baseMultiplyReference(GLBasePolynomial<GLP> &result,
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

template <class T>
inline std::uint32_t unsignedBitWidth(const T &value)
{
    if constexpr (is_multilimb_uint_v<T>) {
        for (std::size_t limb = T::limbs; limb-- > 0;)
            if (value.limb[limb] != 0)
                return static_cast<std::uint32_t>(
                    64 * limb + 64 - std::countl_zero(value.limb[limb]));
        return 0;
    }
    else {
        static_assert(std::is_same_v<T, __uint128_t>);
        const std::uint64_t high = static_cast<std::uint64_t>(value >> 64);
        if (high != 0)
            return 128U - static_cast<std::uint32_t>(std::countl_zero(high));
        const std::uint64_t low = static_cast<std::uint64_t>(value);
        return low == 0
                   ? 0
                   : 64U - static_cast<std::uint32_t>(std::countl_zero(low));
    }
}

template <class T>
inline bool fullTorusIsNegative(const T &value)
{
    if constexpr (is_multilimb_uint_v<T>)
        return (value.limb[T::limbs - 1] >> 63) != 0;
    else {
        static_assert(std::is_same_v<T, __uint128_t>);
        return (value >> 127) != 0;
    }
}

template <class T>
inline std::uint32_t signedTorusBitWidth(const T &value)
{
    return unsignedBitWidth(fullTorusIsNegative(value) ? T{0} - value : value);
}

template <class GLP>
inline std::uint32_t maxSignedTorusBitWidth(
    const GLBasePolynomial<GLP> &polynomial)
{
    std::uint32_t bits = 0;
    for (const auto &coefficient : polynomial)
        bits = std::max(bits, signedTorusBitWidth(coefficient));
    return bits;
}

template <class GLP>
inline std::uint32_t maxUnsignedTorusBitWidth(
    const GLBasePolynomial<GLP> &polynomial)
{
    std::uint32_t bits = 0;
    for (const auto &coefficient : polynomial)
        bits = std::max(bits, unsignedBitWidth(coefficient));
    return bits;
}

template <class GLP>
inline std::uint32_t maxSignedTorusBitWidth(const GLPolynomial<GLP> &polynomial)
{
    std::uint32_t bits = 0;
    for (const auto &slice : polynomial)
        bits = std::max(bits, maxSignedTorusBitWidth<GLP>(slice));
    return bits;
}

template <class GLP>
inline std::uint32_t maxUnsignedTorusBitWidth(
    const GLPolynomial<GLP> &polynomial)
{
    std::uint32_t bits = 0;
    for (const auto &slice : polynomial)
        bits = std::max(bits, maxUnsignedTorusBitWidth<GLP>(slice));
    return bits;
}

template <class T>
inline std::uint64_t signedTorusResidue(const T &value,
                                        const std::uint64_t modulus)
{
    const bool negative = fullTorusIsNegative(value);
    const T magnitude = negative ? T{0} - value : value;
    std::uint64_t residue;
    if constexpr (is_multilimb_uint_v<T>)
        residue = magnitude.mod_u64(modulus);
    else {
        static_assert(std::is_same_v<T, __uint128_t>);
        residue = static_cast<std::uint64_t>(magnitude % modulus);
    }
    return negative && residue != 0 ? modulus - residue : residue;
}

template <class T>
inline T signedI128ToTorus(const __int128 value)
{
    if constexpr (is_multilimb_uint_v<T>) {
        const bool negative = value < 0;
        const unsigned __int128 magnitude =
            negative ? static_cast<unsigned __int128>(-(value + 1)) + 1
                     : static_cast<unsigned __int128>(value);
        T result{};
        result.limb[0] = static_cast<std::uint64_t>(magnitude);
        if constexpr (T::limbs > 1)
            result.limb[1] = static_cast<std::uint64_t>(magnitude >> 64);
        return negative ? T{0} - result : result;
    }
    else {
        static_assert(std::is_same_v<T, __uint128_t>);
        return static_cast<__uint128_t>(value);
    }
}

template <class T>
inline __uint128_t torusLowU128(const T &value)
{
    if constexpr (is_multilimb_uint_v<T>) {
        __uint128_t result = value.limb[0];
        if constexpr (T::limbs > 1)
            result |= static_cast<__uint128_t>(value.limb[1]) << 64;
        return result;
    }
    else {
        static_assert(std::is_same_v<T, __uint128_t>);
        return value;
    }
}

template <class T>
inline T torusFromLowU128(const __uint128_t value)
{
    if constexpr (is_multilimb_uint_v<T>) {
        T result{};
        result.limb[0] = static_cast<std::uint64_t>(value);
        if constexpr (T::limbs > 1)
            result.limb[1] = static_cast<std::uint64_t>(value >> 64);
        return result;
    }
    else {
        static_assert(std::is_same_v<T, __uint128_t>);
        return value;
    }
}

inline __uint128_t multiplyU128BySignedSmall(const __uint128_t value,
                                             const std::int64_t scalar)
{
    const bool negative = scalar < 0;
    const std::uint64_t magnitude =
        negative ? static_cast<std::uint64_t>(-(scalar + 1)) + 1
                 : static_cast<std::uint64_t>(scalar);
    const __uint128_t product = value * magnitude;
    return negative ? __uint128_t{0} - product : product;
}

constexpr std::uint32_t ceilLog2(const std::size_t value)
{
    if (value <= 1) return 0;
    return static_cast<std::uint32_t>(std::numeric_limits<std::size_t>::digits -
                                      std::countl_zero(value - 1));
}

template <class GLP>
inline constexpr bool supportsWidePrimeNTT = [] {
    constexpr std::uint64_t root_order = 4ULL * GLP::matrix_dimension;
    for (const auto prime : modular_ntt::wide_primes)
        if ((prime.value - 1) % root_order != 0 ||
            (prime.value - 1) % GLP::cyclotomic_order != 0)
            return false;
    return true;
}();

template <class GLP>
class GLBaseNTTPlan {
public:
    static constexpr std::size_t z_dimension = 2 * GLP::matrix_dimension;
    static constexpr std::size_t w_dimension = GLP::phi;
    static constexpr std::size_t coefficient_count = z_dimension * w_dimension;

    explicit GLBaseNTTPlan(const modular_ntt::PrimeModulus prime)
        : prime_(prime), w_plan_(prime), z_plan_(z_dimension, prime)
    {
        static_assert(coefficient_count == GLP::baseP::n);
    }

    std::uint64_t modulus() const { return prime_.value; }

    void forward(std::vector<std::uint64_t> &spectrum,
                 const GLBasePolynomial<GLP> &input) const
    {
        spectrum.assign(coefficient_count, 0);
        forward(std::span<std::uint64_t>(spectrum), input);
    }

    void forward(const std::span<std::uint64_t> spectrum,
                 const GLBasePolynomial<GLP> &input) const
    {
        if (spectrum.size() != coefficient_count)
            throw std::invalid_argument("GL NTT spectrum has the wrong size");
        if (inOpenMPParallelRegion()) {
            std::array<std::uint64_t, w_dimension> w_line{};
            for (std::size_t z = 0; z < z_dimension; z++) {
                const std::uint32_t gaussian =
                    static_cast<std::uint32_t>(z / GLP::matrix_dimension);
                const std::uint32_t x =
                    static_cast<std::uint32_t>(z % GLP::matrix_dimension);
                for (std::uint32_t w = 0; w < w_dimension; w++)
                    w_line[w] = signedTorusResidue(
                        input[baseIndex<GLP>(gaussian, x, w)], prime_.value);
                w_plan_.forward(w_line);
                for (std::size_t w = 0; w < w_dimension; w++)
                    spectrum[w * z_dimension + z] = w_line[w];
            }
            for (std::size_t w = 0; w < w_dimension; w++)
                z_plan_.forwardInBackendOrder(std::span<std::uint64_t>(
                    spectrum.data() + w * z_dimension, z_dimension));
            return;
        }
#pragma omp parallel
        {
            std::array<std::uint64_t, w_dimension> w_line{};
#pragma omp for schedule(static)
            for (std::size_t z = 0; z < z_dimension; z++) {
                const std::uint32_t gaussian =
                    static_cast<std::uint32_t>(z / GLP::matrix_dimension);
                const std::uint32_t x =
                    static_cast<std::uint32_t>(z % GLP::matrix_dimension);
                for (std::uint32_t w = 0; w < w_dimension; w++)
                    w_line[w] = signedTorusResidue(
                        input[baseIndex<GLP>(gaussian, x, w)], prime_.value);
                w_plan_.forward(w_line);
                for (std::size_t w = 0; w < w_dimension; w++)
                    spectrum[w * z_dimension + z] = w_line[w];
            }
        }
#pragma omp parallel for schedule(static)
        for (std::size_t w = 0; w < w_dimension; w++)
            z_plan_.forwardInBackendOrder(std::span<std::uint64_t>(
                spectrum.data() + w * z_dimension, z_dimension));
    }

    void forwardResidues(const std::span<std::uint64_t> spectrum,
                         const std::span<const std::uint64_t> input) const
    {
        if (spectrum.size() != coefficient_count ||
            input.size() != coefficient_count)
            throw std::invalid_argument(
                "GL NTT residue input has the wrong size");
        if (inOpenMPParallelRegion()) {
            std::array<std::uint64_t, w_dimension> w_line{};
            for (std::size_t z = 0; z < z_dimension; z++) {
                const std::uint32_t gaussian =
                    static_cast<std::uint32_t>(z / GLP::matrix_dimension);
                const std::uint32_t x =
                    static_cast<std::uint32_t>(z % GLP::matrix_dimension);
                for (std::uint32_t w = 0; w < w_dimension; w++)
                    w_line[w] = input[baseIndex<GLP>(gaussian, x, w)];
                w_plan_.forward(w_line);
                for (std::size_t w = 0; w < w_dimension; w++)
                    spectrum[w * z_dimension + z] = w_line[w];
            }
            for (std::size_t w = 0; w < w_dimension; w++)
                z_plan_.forwardInBackendOrder(std::span<std::uint64_t>(
                    spectrum.data() + w * z_dimension, z_dimension));
            return;
        }
#pragma omp parallel
        {
            std::array<std::uint64_t, w_dimension> w_line{};
#pragma omp for schedule(static)
            for (std::size_t z = 0; z < z_dimension; z++) {
                const std::uint32_t gaussian =
                    static_cast<std::uint32_t>(z / GLP::matrix_dimension);
                const std::uint32_t x =
                    static_cast<std::uint32_t>(z % GLP::matrix_dimension);
                for (std::uint32_t w = 0; w < w_dimension; w++)
                    w_line[w] = input[baseIndex<GLP>(gaussian, x, w)];
                w_plan_.forward(w_line);
                for (std::size_t w = 0; w < w_dimension; w++)
                    spectrum[w * z_dimension + z] = w_line[w];
            }
        }
#pragma omp parallel for schedule(static)
        for (std::size_t w = 0; w < w_dimension; w++)
            z_plan_.forwardInBackendOrder(std::span<std::uint64_t>(
                spectrum.data() + w * z_dimension, z_dimension));
    }

    void inverse(std::vector<std::uint64_t> &coefficients,
                 std::vector<std::uint64_t> &spectrum) const
    {
        coefficients.assign(coefficient_count, 0);
        inverse(std::span<std::uint64_t>(coefficients),
                std::span<std::uint64_t>(spectrum));
    }

    void inverse(const std::span<std::uint64_t> coefficients,
                 const std::span<std::uint64_t> spectrum) const
    {
        if (spectrum.size() != coefficient_count)
            throw std::invalid_argument("GL NTT spectrum has the wrong size");
        if (coefficients.size() != coefficient_count)
            throw std::invalid_argument(
                "GL NTT coefficient output has the wrong size");
        if (inOpenMPParallelRegion()) {
            for (std::size_t w = 0; w < w_dimension; w++)
                z_plan_.inverseInBackendOrder(std::span<std::uint64_t>(
                    spectrum.data() + w * z_dimension, z_dimension));
            std::array<std::uint64_t, w_dimension> w_line{};
            for (std::size_t z = 0; z < z_dimension; z++) {
                for (std::size_t w = 0; w < w_dimension; w++)
                    w_line[w] = spectrum[w * z_dimension + z];
                w_plan_.inverse(w_line);
                const std::uint32_t gaussian =
                    static_cast<std::uint32_t>(z / GLP::matrix_dimension);
                const std::uint32_t x =
                    static_cast<std::uint32_t>(z % GLP::matrix_dimension);
                for (std::uint32_t w = 0; w < w_dimension; w++)
                    coefficients[baseIndex<GLP>(gaussian, x, w)] = w_line[w];
            }
            return;
        }
#pragma omp parallel for schedule(static)
        for (std::size_t w = 0; w < w_dimension; w++)
            z_plan_.inverseInBackendOrder(std::span<std::uint64_t>(
                spectrum.data() + w * z_dimension, z_dimension));

#pragma omp parallel
        {
            std::array<std::uint64_t, w_dimension> w_line{};
#pragma omp for schedule(static)
            for (std::size_t z = 0; z < z_dimension; z++) {
                for (std::size_t w = 0; w < w_dimension; w++)
                    w_line[w] = spectrum[w * z_dimension + z];
                w_plan_.inverse(w_line);
                const std::uint32_t gaussian =
                    static_cast<std::uint32_t>(z / GLP::matrix_dimension);
                const std::uint32_t x =
                    static_cast<std::uint32_t>(z % GLP::matrix_dimension);
                for (std::uint32_t w = 0; w < w_dimension; w++)
                    coefficients[baseIndex<GLP>(gaussian, x, w)] = w_line[w];
            }
        }
    }

private:
    modular_ntt::PrimeModulus prime_;
    modular_ntt::PrimeCyclotomicNTTPlan<GLP::cyclotomic_order> w_plan_;
    modular_ntt::NegacyclicNTTPlan z_plan_;
};

template <class GLP, std::size_t PrimeIndex>
inline const GLBaseNTTPlan<GLP> &baseNTTPlan()
{
    static_assert(PrimeIndex < modular_ntt::wide_primes.size());
    static const GLBaseNTTPlan<GLP> plan(modular_ntt::wide_primes[PrimeIndex]);
    return plan;
}

template <class GLP>
using GLBaseXAutomorphismSpectrumMap =
    std::array<std::uint32_t, GLBaseNTTPlan<GLP>::z_dimension>;

template <class GLP>
inline GLBaseXAutomorphismSpectrumMap<GLP> baseXAutomorphismSpectrumMap(
    const std::uint32_t x_multiplier)
{
    constexpr std::size_t z_dimension = GLBaseNTTPlan<GLP>::z_dimension;
    constexpr std::uint32_t four_n = 4 * GLP::matrix_dimension;
    if ((x_multiplier & 1U) == 0)
        throw std::invalid_argument("GL X automorphism multiplier must be odd");
    GLBaseXAutomorphismSpectrumMap<GLP> result{};
    for (std::size_t z = 0; z < z_dimension; z++) {
        const std::uint32_t odd_root = static_cast<std::uint32_t>(2 * z + 1);
        const std::uint32_t mapped_root = static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(x_multiplier) * odd_root) % four_n);
        const std::size_t destination =
            modular_ntt::NegacyclicNTTPlan::backendIndex(z_dimension, z);
        result[destination] = static_cast<std::uint32_t>(
            modular_ntt::NegacyclicNTTPlan::backendIndex(
                z_dimension, (mapped_root - 1) / 2));
    }
    return result;
}

template <class GLP>
inline void applyBaseXAutomorphismSpectrum(
    const std::span<std::uint64_t> destination,
    const std::span<const std::uint64_t> source,
    const GLBaseXAutomorphismSpectrumMap<GLP> &z_map)
{
    constexpr std::size_t coefficient_count =
        GLBaseNTTPlan<GLP>::coefficient_count;
    constexpr std::size_t z_dimension = GLBaseNTTPlan<GLP>::z_dimension;
    constexpr std::size_t w_dimension = GLBaseNTTPlan<GLP>::w_dimension;
    if (destination.size() != coefficient_count ||
        source.size() != coefficient_count)
        throw std::invalid_argument("GL NTT spectrum has the wrong size");
    if (destination.data() == source.data())
        throw std::invalid_argument(
            "in-place GL spectrum automorphism is unsupported");
    for (std::size_t w = 0; w < w_dimension; w++)
        for (std::size_t z = 0; z < z_dimension; z++)
            destination[w * z_dimension + z] =
                source[w * z_dimension + z_map[z]];
}

// Set U=Y/X.  Since X^n=Y^n=I, U^n=1, so the full GL ring is a cyclic
// length-n extension of the already transformed (I,X,W) base ring.  A Y slice
// is first multiplied by X^y, transformed in the base ring, and then
// transformed along U.  Spectra use base-frequency-major layout so every
// length-n U transform is contiguous.
template <class GLP>
class GLPolynomialNTTPlan {
public:
    static constexpr std::size_t y_dimension = GLP::matrix_dimension;
    static constexpr std::size_t base_coefficient_count =
        GLBaseNTTPlan<GLP>::coefficient_count;
    static constexpr std::size_t coefficient_count =
        y_dimension * base_coefficient_count;

    explicit GLPolynomialNTTPlan(const modular_ntt::PrimeModulus prime)
        : prime_(prime), base_plan_(prime), y_plan_(y_dimension, prime)
    {
    }

    std::uint64_t modulus() const { return prime_.value; }

    void forward(std::vector<std::uint64_t> &spectrum,
                 const GLPolynomial<GLP> &input) const
    {
        forwardWith(
            spectrum, [&](const std::size_t y, const std::size_t coefficient) {
                return signedTorusResidue(input[y][coefficient], prime_.value);
            });
    }

    template <std::uint32_t Bits>
    void forwardPacked(std::vector<std::uint64_t> &spectrum,
                       const GLPackedPolynomial<GLP, Bits> &input) const
    {
        forwardWith(spectrum,
                    [&](const std::size_t y, const std::size_t coefficient) {
                        return signedSmallResidue(input[y][coefficient]);
                    });
    }

    template <std::uint32_t LogQ, std::uint32_t BaseBit>
    void forwardActiveRows(std::vector<std::uint64_t> &spectra,
                           const GLPolynomial<GLP> &input) const
    {
        using P = typename GLP::baseP;
        constexpr std::uint32_t rows = (LogQ + BaseBit - 1) / BaseBit;
        static_assert(rows * BaseBit <=
                      std::numeric_limits<typename P::T>::digits);
        spectra.assign(static_cast<std::size_t>(rows) * coefficient_count, 0);
#pragma omp parallel
        {
            auto slice_rows =
                std::make_unique<std::array<GLBasePolynomial<GLP>, rows>>();
            std::vector<std::uint64_t> twisted(base_coefficient_count);
            std::vector<std::uint64_t> base_spectrum(base_coefficient_count);
#pragma omp for schedule(static)
            for (std::size_t y = 0; y < y_dimension; y++) {
                ckks_detail::activeBaseDecomposePolynomialRows<P, LogQ, BaseBit,
                                                               rows>(
                    *slice_rows, input[y]);
                for (std::uint32_t row = 0; row < rows; row++) {
                    fillTwistedResidues(
                        twisted, y, [&](const std::size_t coefficient) {
                            return signedTorusResidue(
                                (*slice_rows)[row][coefficient], prime_.value);
                        });
                    base_plan_.forwardResidues(base_spectrum, twisted);
                    std::uint64_t *row_spectrum =
                        spectra.data() +
                        static_cast<std::size_t>(row) * coefficient_count;
                    for (std::size_t base = 0; base < base_coefficient_count;
                         base++)
                        row_spectrum[base * y_dimension + y] =
                            base_spectrum[base];
                }
            }
        }
        for (std::uint32_t row = 0; row < rows; row++)
            finishForward(std::span<std::uint64_t>(
                spectra.data() +
                    static_cast<std::size_t>(row) * coefficient_count,
                coefficient_count));
    }

    // Inverse transforms spectrum and returns canonical coefficients in
    // [Y-slice][base coefficient] order.  spectrum is mutable scratch.
    void inverse(std::vector<std::uint64_t> &coefficients,
                 std::vector<std::uint64_t> &spectrum) const
    {
        if (spectrum.size() != coefficient_count)
            throw std::invalid_argument(
                "GL polynomial NTT spectrum has the wrong size");
#pragma omp parallel for schedule(static)
        for (std::size_t base = 0; base < base_coefficient_count; base++)
            y_plan_.inverseInBackendOrder(std::span<std::uint64_t>(
                spectrum.data() + base * y_dimension, y_dimension));

        coefficients.assign(coefficient_count, 0);
#pragma omp parallel
        {
            std::vector<std::uint64_t> base_spectrum(base_coefficient_count);
            std::vector<std::uint64_t> twisted(base_coefficient_count);
#pragma omp for schedule(static)
            for (std::size_t y = 0; y < y_dimension; y++) {
                for (std::size_t base = 0; base < base_coefficient_count;
                     base++)
                    base_spectrum[base] = spectrum[base * y_dimension + y];
                base_plan_.inverse(twisted, base_spectrum);
                multiplyBaseByXPowerResidues(
                    std::span<std::uint64_t>(
                        coefficients.data() + y * base_coefficient_count,
                        base_coefficient_count),
                    twisted, (4 * y_dimension - y) % (4 * y_dimension));
            }
        }
    }

private:
    template <class Getter>
    void forwardWith(std::vector<std::uint64_t> &spectrum, Getter &&get) const
    {
        spectrum.assign(coefficient_count, 0);
#pragma omp parallel
        {
            std::vector<std::uint64_t> twisted(base_coefficient_count);
            std::vector<std::uint64_t> base_spectrum(base_coefficient_count);
#pragma omp for schedule(static)
            for (std::size_t y = 0; y < y_dimension; y++) {
                fillTwistedResidues(twisted, y,
                                    [&](const std::size_t coefficient) {
                                        return get(y, coefficient);
                                    });
                base_plan_.forwardResidues(base_spectrum, twisted);
                for (std::size_t base = 0; base < base_coefficient_count;
                     base++)
                    spectrum[base * y_dimension + y] = base_spectrum[base];
            }
        }
        finishForward(spectrum);
    }

    void finishForward(const std::span<std::uint64_t> spectrum) const
    {
        if (spectrum.size() != coefficient_count)
            throw std::invalid_argument(
                "GL polynomial NTT spectrum has the wrong size");
#pragma omp parallel for schedule(static)
        for (std::size_t base = 0; base < base_coefficient_count; base++)
            y_plan_.forwardInBackendOrder(std::span<std::uint64_t>(
                spectrum.data() + base * y_dimension, y_dimension));
    }

    template <class Getter>
    void fillTwistedResidues(std::vector<std::uint64_t> &destination,
                             const std::size_t y, Getter &&get) const
    {
        constexpr std::size_t n = GLP::matrix_dimension;
        for (std::uint32_t w = 0; w < GLP::phi; w++) {
            for (std::uint32_t x = 0; x < n; x++) {
                const std::size_t total = x + y;
                const std::uint32_t output_x =
                    static_cast<std::uint32_t>(total % n);
                const std::uint32_t carry =
                    static_cast<std::uint32_t>(total / n);
                for (std::uint32_t gaussian = 0; gaussian < 2; gaussian++) {
                    const std::uint32_t phase = gaussian + carry;
                    std::uint64_t residue = get(baseIndex<GLP>(gaussian, x, w));
                    if ((phase & 2U) != 0)
                        residue = modular_ntt::negate(residue, prime_.value);
                    destination[baseIndex<GLP>(phase & 1U, output_x, w)] =
                        residue;
                }
            }
        }
    }

    void multiplyBaseByXPowerResidues(
        const std::span<std::uint64_t> destination,
        const std::span<const std::uint64_t> source,
        const std::size_t exponent) const
    {
        constexpr std::size_t n = GLP::matrix_dimension;
        for (std::uint32_t w = 0; w < GLP::phi; w++) {
            for (std::uint32_t x = 0; x < n; x++) {
                const std::size_t total = x + exponent;
                const std::uint32_t output_x =
                    static_cast<std::uint32_t>(total % n);
                const std::uint32_t blocks =
                    static_cast<std::uint32_t>(total / n);
                for (std::uint32_t gaussian = 0; gaussian < 2; gaussian++) {
                    const std::uint32_t phase = (gaussian + blocks) & 3U;
                    std::uint64_t residue =
                        source[baseIndex<GLP>(gaussian, x, w)];
                    if ((phase & 2U) != 0)
                        residue = modular_ntt::negate(residue, prime_.value);
                    destination[baseIndex<GLP>(phase & 1U, output_x, w)] =
                        residue;
                }
            }
        }
    }

    template <class Digit>
    std::uint64_t signedSmallResidue(const Digit value) const
    {
        const std::int64_t signed_value = value;
        if (signed_value >= 0)
            return static_cast<std::uint64_t>(signed_value) % prime_.value;
        const std::uint64_t magnitude =
            static_cast<std::uint64_t>(-(signed_value + 1)) + 1;
        const std::uint64_t residue = magnitude % prime_.value;
        return residue == 0 ? 0 : prime_.value - residue;
    }

    modular_ntt::PrimeModulus prime_;
    GLBaseNTTPlan<GLP> base_plan_;
    modular_ntt::Radix2NTTPlan y_plan_;
};

template <class GLP, std::size_t PrimeIndex>
inline const GLPolynomialNTTPlan<GLP> &polynomialNTTPlan()
{
    static_assert(PrimeIndex < modular_ntt::wide_primes.size());
    static const GLPolynomialNTTPlan<GLP> plan(
        modular_ntt::wide_primes[PrimeIndex]);
    return plan;
}

template <class GLP, std::size_t PrimeIndex>
inline void baseMultiplyNTTResidues(std::vector<std::uint64_t> &coefficients,
                                    const GLBasePolynomial<GLP> &lhs,
                                    const GLBasePolynomial<GLP> &rhs)
{
    const auto &plan = baseNTTPlan<GLP, PrimeIndex>();
    std::vector<std::uint64_t> lhs_spectrum;
    std::vector<std::uint64_t> rhs_spectrum;
    plan.forward(lhs_spectrum, lhs);
    plan.forward(rhs_spectrum, rhs);
    for (std::size_t i = 0; i < lhs_spectrum.size(); i++)
        lhs_spectrum[i] = modular_ntt::multiply(
            lhs_spectrum[i], rhs_spectrum[i], plan.modulus());
    plan.inverse(coefficients, lhs_spectrum);
}

template <class GLP>
inline std::uint32_t baseMultiplyNTTPrimeCount(const GLBasePolynomial<GLP> &lhs,
                                               const GLBasePolynomial<GLP> &rhs)
{
    if constexpr (!supportsWidePrimeNTT<GLP>) return 0;

    constexpr std::size_t maximum_terms =
        2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
    const std::uint32_t lhs_bits = maxSignedTorusBitWidth<GLP>(lhs);
    const std::uint32_t rhs_bits = maxSignedTorusBitWidth<GLP>(rhs);
    if (lhs_bits == 0 || rhs_bits == 0) return 1;
    const std::uint32_t required_bits =
        lhs_bits + rhs_bits + ceilLog2(maximum_terms);

    // q/2 is slightly below 2^61 and q0*q1/2 is slightly below 2^123.
    // Keeping one full bit below each threshold makes centered reconstruction
    // valid for every coefficient satisfying the detected signed bound.
    if (required_bits <= 60) return 1;
    if (required_bits <= 122) return 2;
    return 0;
}

template <class GLP>
inline void baseMultiplyNTTDirect(GLBasePolynomial<GLP> &result,
                                  const GLBasePolynomial<GLP> &lhs,
                                  const GLBasePolynomial<GLP> &rhs,
                                  const std::uint32_t prime_count)
{
    using T = typename GLP::T;
    assert(prime_count == 1 || prime_count == 2);

    std::vector<std::uint64_t> first_residues;
    baseMultiplyNTTResidues<GLP, 0>(first_residues, lhs, rhs);
    if (prime_count == 1) {
        for (std::size_t i = 0; i < result.size(); i++)
            result[i] = signedI128ToTorus<T>(modular_ntt::centeredResidue(
                first_residues[i], modular_ntt::wide_primes[0].value));
        return;
    }

    std::vector<std::uint64_t> second_residues;
    baseMultiplyNTTResidues<GLP, 1>(second_residues, lhs, rhs);
    static const modular_ntt::TwoPrimeCRT crt(modular_ntt::wide_primes[0],
                                              modular_ntt::wide_primes[1]);
    for (std::size_t i = 0; i < result.size(); i++)
        result[i] = signedI128ToTorus<T>(
            crt.reconstructSigned(first_residues[i], second_residues[i]));
}

template <class T>
inline T unsignedChunk(const T &value, const std::uint32_t shift,
                       const std::uint32_t bits)
{
    constexpr std::uint32_t width = std::numeric_limits<T>::digits;
    if (shift >= width || bits == 0) return T{0};
    const std::uint32_t available = std::min(bits, width - shift);
    const T mask = available == width ? std::numeric_limits<T>::max()
                                      : (T{1} << available) - T{1};
    return (value >> shift) & mask;
}

// If neither operand is small enough for centered CRT reconstruction, split
// both of their unsigned power-of-two representatives.  Transform every chunk
// once, accumulate products with the same output shift in the NTT domain, and
// perform one inverse transform per prime and output diagonal.  The exact CRT
// bound includes every chunk pair summed on that diagonal.  Unsigned
// decomposition is exact modulo the requested power-of-two level, including
// for negative values represented by their low LogQ bits.
template <class GLP, std::uint32_t OutputBits =
                         std::numeric_limits<typename GLP::T>::digits>
inline bool baseMultiplyNTTDoubleChunk(GLBasePolynomial<GLP> &result,
                                       const GLBasePolynomial<GLP> &lhs,
                                       const GLBasePolynomial<GLP> &rhs)
{
    if constexpr (!supportsWidePrimeNTT<GLP>) return false;

    using T = typename GLP::T;
    constexpr std::size_t maximum_terms =
        2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
    constexpr std::uint32_t convolution_bits = ceilLog2(maximum_terms);
    constexpr std::uint32_t two_prime_safe_bits = 122;
    constexpr std::uint32_t torus_width = OutputBits;
    static_assert(OutputBits > 0 &&
                  OutputBits <= std::numeric_limits<T>::digits);
    constexpr std::uint32_t chunk_bits = [] {
        for (std::uint32_t bits = (two_prime_safe_bits - convolution_bits) / 2;
             bits > 0; bits--) {
            const std::size_t maximum_chunks = (torus_width + bits - 1) / bits;
            if (2 * bits + convolution_bits + ceilLog2(maximum_chunks) <=
                two_prime_safe_bits)
                return bits;
        }
        return 0U;
    }();
    static_assert(chunk_bits > 0);
    const std::uint32_t lhs_width = maxUnsignedTorusBitWidth<GLP>(lhs);
    const std::uint32_t rhs_width = maxUnsignedTorusBitWidth<GLP>(rhs);
    if (lhs_width == 0 || rhs_width == 0) {
        clear<GLP>(result);
        return true;
    }
    const std::size_t lhs_chunks = (lhs_width + chunk_bits - 1) / chunk_bits;
    const std::size_t rhs_chunks = (rhs_width + chunk_bits - 1) / chunk_bits;
    constexpr std::size_t coefficient_count =
        GLBaseNTTPlan<GLP>::coefficient_count;
    const std::array<const GLBaseNTTPlan<GLP> *, 2> plans{
        &baseNTTPlan<GLP, 0>(), &baseNTTPlan<GLP, 1>()};
    std::array<std::vector<std::vector<std::uint64_t>>, 2> lhs_spectra;
    std::array<std::vector<std::vector<std::uint64_t>>, 2> rhs_spectra;
    for (std::size_t prime = 0; prime < 2; prime++) {
        lhs_spectra[prime].resize(
            lhs_chunks, std::vector<std::uint64_t>(coefficient_count));
        rhs_spectra[prime].resize(
            rhs_chunks, std::vector<std::uint64_t>(coefficient_count));
    }
#pragma omp parallel
    {
        auto chunk = std::make_unique<GLBasePolynomial<GLP>>();
#pragma omp for collapse(2) schedule(static)
        for (std::size_t prime = 0; prime < 2; prime++) {
            for (std::size_t row = 0; row < lhs_chunks; row++) {
                const std::uint32_t shift =
                    static_cast<std::uint32_t>(row * chunk_bits);
                for (std::size_t i = 0; i < coefficient_count; i++)
                    (*chunk)[i] = unsignedChunk(lhs[i], shift, chunk_bits);
                plans[prime]->forward(
                    std::span<std::uint64_t>(lhs_spectra[prime][row]), *chunk);
            }
        }
    }
#pragma omp parallel
    {
        auto chunk = std::make_unique<GLBasePolynomial<GLP>>();
#pragma omp for collapse(2) schedule(static)
        for (std::size_t prime = 0; prime < 2; prime++) {
            for (std::size_t row = 0; row < rhs_chunks; row++) {
                const std::uint32_t shift =
                    static_cast<std::uint32_t>(row * chunk_bits);
                for (std::size_t i = 0; i < coefficient_count; i++)
                    (*chunk)[i] = unsignedChunk(rhs[i], shift, chunk_bits);
                plans[prime]->forward(
                    std::span<std::uint64_t>(rhs_spectra[prime][row]), *chunk);
            }
        }
    }

    clear<GLP>(result);
    const std::size_t diagonal_count =
        std::min<std::size_t>(lhs_chunks + rhs_chunks - 1,
                              (torus_width + chunk_bits - 1) / chunk_bits);
    std::array<std::vector<std::uint64_t>, 2> diagonal_coefficients{
        std::vector<std::uint64_t>(diagonal_count * coefficient_count),
        std::vector<std::uint64_t>(diagonal_count * coefficient_count)};
#pragma omp parallel
    {
        std::array<std::vector<std::uint64_t>, 2> accumulators{
            std::vector<std::uint64_t>(coefficient_count),
            std::vector<std::uint64_t>(coefficient_count)};
        std::array<std::vector<std::uint64_t>, 2> coefficients{
            std::vector<std::uint64_t>(coefficient_count),
            std::vector<std::uint64_t>(coefficient_count)};
#pragma omp for schedule(dynamic)
        for (std::size_t diagonal = 0; diagonal < diagonal_count; diagonal++) {
            const std::size_t lhs_begin =
                diagonal >= rhs_chunks ? diagonal - rhs_chunks + 1 : 0;
            const std::size_t lhs_end = std::min(diagonal + 1, lhs_chunks);
            for (std::size_t prime = 0; prime < 2; prime++) {
                auto &accumulator = accumulators[prime];
                std::fill(accumulator.begin(), accumulator.end(), 0);
                const std::uint64_t modulus = plans[prime]->modulus();
                for (std::size_t lhs_row = lhs_begin; lhs_row < lhs_end;
                     lhs_row++) {
                    const std::size_t rhs_row = diagonal - lhs_row;
                    for (std::size_t i = 0; i < coefficient_count; i++) {
                        const std::uint64_t product = modular_ntt::multiply(
                            lhs_spectra[prime][lhs_row][i],
                            rhs_spectra[prime][rhs_row][i], modulus);
                        accumulator[i] =
                            modular_ntt::add(accumulator[i], product, modulus);
                    }
                }
                plans[prime]->inverse(
                    std::span<std::uint64_t>(coefficients[prime]),
                    std::span<std::uint64_t>(accumulator));
                std::copy(coefficients[prime].begin(),
                          coefficients[prime].end(),
                          diagonal_coefficients[prime].begin() +
                              diagonal * coefficient_count);
            }
        }
    }
    static const modular_ntt::TwoPrimeCRT crt(modular_ntt::wide_primes[0],
                                              modular_ntt::wide_primes[1]);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < coefficient_count; i++) {
        T value{};
        for (std::size_t diagonal = 0; diagonal < diagonal_count; diagonal++)
            value +=
                signedI128ToTorus<T>(crt.reconstructSigned(
                    diagonal_coefficients[0][diagonal * coefficient_count + i],
                    diagonal_coefficients[1][diagonal * coefficient_count + i]))
                << (diagonal * chunk_bits);
        result[i] = value;
    }
    return true;
}

template <class GLP>
struct BaseCiphertextTensorNTTWorkspace {
    std::array<std::vector<std::uint64_t>, 2> source_spectra{};
    std::array<std::vector<std::uint64_t>, 2> diagonal_coefficients{};
    std::array<std::vector<std::uint64_t>, 2> accumulators{};
    std::array<std::vector<std::uint64_t>, 2> coefficients{};
    std::vector<unsigned __int128> wide{};
    std::unique_ptr<GLBasePolynomial<GLP>> chunk{};
};

// Compute (a0*b0, a0*b1+a1*b0, a1*b1) together. Chunk spectra for each of
// the four source components are shared, and the two cross products are
// accumulated before each inverse transform. The extra cross-term factor is
// included in the exact two-prime bound.
template <class GLP, std::uint32_t OutputBits>
inline bool baseCiphertextTensorMultiplyNTT(
    std::array<GLBasePolynomial<GLP>, 3> &result,
    const GLBaseCiphertextData<GLP> &lhs, const GLBaseCiphertextData<GLP> &rhs,
    BaseCiphertextTensorNTTWorkspace<GLP> *provided_workspace = nullptr)
{
    if constexpr (!supportsWidePrimeNTT<GLP>) return false;

    using T = typename GLP::T;
    constexpr std::size_t maximum_terms =
        2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
    constexpr std::uint32_t convolution_bits = ceilLog2(maximum_terms);
    constexpr std::uint32_t two_prime_safe_bits = 122;
    static_assert(OutputBits > 0 &&
                  OutputBits <= std::numeric_limits<T>::digits);
    constexpr std::uint32_t chunk_bits = [] {
        for (std::uint32_t bits = (two_prime_safe_bits - convolution_bits) / 2;
             bits > 0; bits--) {
            const std::size_t chunks = (OutputBits + bits - 1) / bits;
            const std::size_t maximum_cross_pairs = 2 * chunks;
            if (2 * bits + convolution_bits + ceilLog2(maximum_cross_pairs) <=
                two_prime_safe_bits)
                return bits;
        }
        return 0U;
    }();
    static_assert(chunk_bits > 0);
    constexpr std::size_t chunk_count =
        (OutputBits + chunk_bits - 1) / chunk_bits;
    constexpr std::size_t diagonal_count = chunk_count;
    constexpr std::size_t coefficient_count =
        GLBaseNTTPlan<GLP>::coefficient_count;
    const std::array<const GLBasePolynomial<GLP> *, 4> sources{
        &lhs[0], &lhs[1], &rhs[0], &rhs[1]};
    const std::array<const GLBaseNTTPlan<GLP> *, 2> plans{
        &baseNTTPlan<GLP, 0>(), &baseNTTPlan<GLP, 1>()};
    static const modular_ntt::TwoPrimeCRT crt(modular_ntt::wide_primes[0],
                                              modular_ntt::wide_primes[1]);

    BaseCiphertextTensorNTTWorkspace<GLP> local_workspace;
    auto &workspace =
        provided_workspace == nullptr ? local_workspace : *provided_workspace;
    const std::size_t source_spectrum_size =
        4 * chunk_count * coefficient_count;
    for (auto &spectra : workspace.source_spectra)
        spectra.resize(source_spectrum_size);
    for (auto &values : workspace.diagonal_coefficients)
        values.resize(diagonal_count * coefficient_count);
    for (auto &values : workspace.accumulators)
        values.resize(coefficient_count);
    for (auto &values : workspace.coefficients)
        values.resize(coefficient_count);
    workspace.wide.resize(coefficient_count);
    if (!workspace.chunk)
        workspace.chunk = std::make_unique<GLBasePolynomial<GLP>>();

    for (std::size_t source = 0; source < sources.size(); source++) {
        for (std::size_t chunk = 0; chunk < chunk_count; chunk++) {
            const std::uint32_t shift =
                static_cast<std::uint32_t>(chunk * chunk_bits);
            for (std::size_t i = 0; i < coefficient_count; i++)
                (*workspace.chunk)[i] =
                    unsignedChunk((*sources[source])[i], shift, chunk_bits);
            for (std::size_t prime = 0; prime < 2; prime++) {
                const std::size_t spectrum_row = source * chunk_count + chunk;
                plans[prime]->forward(
                    std::span<std::uint64_t>(
                        workspace.source_spectra[prime].data() +
                            spectrum_row * coefficient_count,
                        coefficient_count),
                    *workspace.chunk);
            }
        }
    }

    constexpr std::array<std::size_t, 3> pair_counts{1, 2, 1};
    constexpr std::array<std::array<std::size_t, 2>, 3> lhs_sources{{
        {{0, 0}},
        {{0, 1}},
        {{1, 1}},
    }};
    constexpr std::array<std::array<std::size_t, 2>, 3> rhs_sources{{
        {{2, 2}},
        {{3, 2}},
        {{3, 3}},
    }};
    for (std::size_t output = 0; output < result.size(); output++) {
        for (std::size_t diagonal = 0; diagonal < diagonal_count; diagonal++) {
            for (std::size_t prime = 0; prime < 2; prime++) {
                std::fill(workspace.wide.begin(), workspace.wide.end(), 0);
                for (std::size_t pair = 0; pair < pair_counts[output]; pair++) {
                    for (std::size_t lhs_chunk = 0; lhs_chunk <= diagonal;
                         lhs_chunk++) {
                        const std::size_t rhs_chunk = diagonal - lhs_chunk;
                        if (lhs_chunk >= chunk_count ||
                            rhs_chunk >= chunk_count)
                            continue;
                        const std::size_t lhs_row =
                            lhs_sources[output][pair] * chunk_count + lhs_chunk;
                        const std::size_t rhs_row =
                            rhs_sources[output][pair] * chunk_count + rhs_chunk;
                        const std::uint64_t *lhs_spectrum =
                            workspace.source_spectra[prime].data() +
                            lhs_row * coefficient_count;
                        const std::uint64_t *rhs_spectrum =
                            workspace.source_spectra[prime].data() +
                            rhs_row * coefficient_count;
                        for (std::size_t i = 0; i < coefficient_count; i++)
                            workspace.wide[i] += static_cast<unsigned __int128>(
                                                     lhs_spectrum[i]) *
                                                 rhs_spectrum[i];
                    }
                }
                const std::uint64_t modulus = plans[prime]->modulus();
                for (std::size_t i = 0; i < coefficient_count; i++)
                    workspace.accumulators[prime][i] =
                        modular_ntt::reduceWide(workspace.wide[i], modulus);
                plans[prime]->inverse(
                    std::span<std::uint64_t>(workspace.coefficients[prime]),
                    std::span<std::uint64_t>(workspace.accumulators[prime]));
                std::copy(workspace.coefficients[prime].begin(),
                          workspace.coefficients[prime].end(),
                          workspace.diagonal_coefficients[prime].begin() +
                              diagonal * coefficient_count);
            }
        }

        for (std::size_t i = 0; i < coefficient_count; i++) {
            T value{};
            for (std::size_t diagonal = 0; diagonal < diagonal_count;
                 diagonal++)
                value += signedI128ToTorus<T>(crt.reconstructSigned(
                             workspace.diagonal_coefficients
                                 [0][diagonal * coefficient_count + i],
                             workspace.diagonal_coefficients
                                 [1][diagonal * coefficient_count + i]))
                         << (diagonal * chunk_bits);
            result[output][i] = value;
        }
        reduce<GLP, OutputBits>(result[output]);
    }
    return true;
}

template <class GLP>
inline bool baseMultiplyNTT(GLBasePolynomial<GLP> &result,
                            const GLBasePolynomial<GLP> &lhs,
                            const GLBasePolynomial<GLP> &rhs)
{
    if constexpr (!supportsWidePrimeNTT<GLP>) return false;

    const std::uint32_t direct_primes =
        baseMultiplyNTTPrimeCount<GLP>(lhs, rhs);
    if (direct_primes != 0) {
        baseMultiplyNTTDirect<GLP>(result, lhs, rhs, direct_primes);
        return true;
    }

    // A wide torus polynomial multiplied by a bounded signed polynomial is
    // common during encryption and evaluation-key generation.  Expand only
    // the wide operand into unsigned radix-2 chunks; each chunk product then
    // fits the same exact two-prime CRT.  Recombination is modulo the native
    // power-of-two torus, so this remains exact even for a negative wide
    // coefficient represented in two's complement.
    constexpr std::size_t maximum_terms =
        2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
    constexpr std::uint32_t convolution_bits = ceilLog2(maximum_terms);
    constexpr std::uint32_t crt_safe_bits = 122;
    constexpr std::uint32_t torus_width =
        std::numeric_limits<typename GLP::T>::digits;
    const std::uint32_t lhs_bits = maxSignedTorusBitWidth<GLP>(lhs);
    const std::uint32_t rhs_bits = maxSignedTorusBitWidth<GLP>(rhs);
    const bool split_lhs = lhs_bits >= rhs_bits;
    const std::uint32_t bounded_bits = split_lhs ? rhs_bits : lhs_bits;
    if (bounded_bits == 0) {
        clear<GLP>(result);
        return true;
    }
    if (bounded_bits + convolution_bits >= crt_safe_bits)
        return baseMultiplyNTTDoubleChunk<GLP>(result, lhs, rhs);

    const std::uint32_t chunk_bits =
        crt_safe_bits - bounded_bits - convolution_bits;
    const auto &wide = split_lhs ? lhs : rhs;
    const auto &bounded = split_lhs ? rhs : lhs;
    GLBasePolynomial<GLP> chunk{};
    GLBasePolynomial<GLP> chunk_product{};
    clear<GLP>(result);
    for (std::uint32_t shift = 0; shift < torus_width; shift += chunk_bits) {
        bool nonzero = false;
        for (std::size_t i = 0; i < chunk.size(); i++) {
            chunk[i] = unsignedChunk(wide[i], shift, chunk_bits);
            nonzero = nonzero || chunk[i] != typename GLP::T{0};
        }
        if (!nonzero) continue;
        const std::uint32_t chunk_primes =
            baseMultiplyNTTPrimeCount<GLP>(chunk, bounded);
        if (chunk_primes == 0)
            return baseMultiplyNTTDoubleChunk<GLP>(result, lhs, rhs);
        baseMultiplyNTTDirect<GLP>(chunk_product, chunk, bounded, chunk_primes);
        for (std::size_t i = 0; i < result.size(); i++)
            result[i] += chunk_product[i] << shift;
    }
    return true;
}

template <class GLP>
inline void baseMultiply(GLBasePolynomial<GLP> &result,
                         const GLBasePolynomial<GLP> &lhs,
                         const GLBasePolynomial<GLP> &rhs)
{
    if (!baseMultiplyNTT<GLP>(result, lhs, rhs))
        baseMultiplyReference<GLP>(result, lhs, rhs);
}

template <class GLP, std::uint32_t LogQ>
inline void baseMultiplyAtLevel(GLBasePolynomial<GLP> &result,
                                const GLBasePolynomial<GLP> &lhs,
                                const GLBasePolynomial<GLP> &rhs)
{
    static_assert(LogQ > 0 &&
                  LogQ <= std::numeric_limits<typename GLP::T>::digits);
    if constexpr (supportsWidePrimeNTT<GLP>) {
        const std::uint32_t direct_primes =
            baseMultiplyNTTPrimeCount<GLP>(lhs, rhs);
        if (direct_primes != 0) {
            baseMultiplyNTTDirect<GLP>(result, lhs, rhs, direct_primes);
            reduce<GLP, LogQ>(result);
            return;
        }
        if (baseMultiplyNTTDoubleChunk<GLP, LogQ>(result, lhs, rhs)) {
            reduce<GLP, LogQ>(result);
            return;
        }
    }
    baseMultiplyReference<GLP>(result, lhs, rhs);
    reduce<GLP, LogQ>(result);
}

template <class GLP>
inline void polynomialMultiplyReference(GLPolynomial<GLP> &result,
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
inline std::uint32_t polynomialMultiplyNTTPrimeCount(
    const GLPolynomial<GLP> &lhs, const GLPolynomial<GLP> &rhs)
{
    if constexpr (!supportsWidePrimeNTT<GLP>) return 0;
    constexpr std::size_t maximum_terms =
        static_cast<std::size_t>(GLP::matrix_dimension) * 2 *
        GLP::matrix_dimension * (2 * GLP::phi - 1);
    const std::uint32_t lhs_bits = maxSignedTorusBitWidth<GLP>(lhs);
    const std::uint32_t rhs_bits = maxSignedTorusBitWidth<GLP>(rhs);
    if (lhs_bits == 0 || rhs_bits == 0) return 1;
    const std::uint32_t required_bits =
        lhs_bits + rhs_bits + ceilLog2(maximum_terms);
    if (required_bits <= 60) return 1;
    if (required_bits <= 122) return 2;
    if (required_bits <= 126) return 3;
    return 0;
}

template <class GLP, std::size_t PrimeIndex>
inline void polynomialMultiplyNTTResidues(
    std::vector<std::uint64_t> &coefficients, const GLPolynomial<GLP> &lhs,
    const GLPolynomial<GLP> &rhs)
{
    const auto &plan = polynomialNTTPlan<GLP, PrimeIndex>();
    std::vector<std::uint64_t> lhs_spectrum;
    std::vector<std::uint64_t> rhs_spectrum;
    plan.forward(lhs_spectrum, lhs);
    plan.forward(rhs_spectrum, rhs);
    for (std::size_t i = 0; i < lhs_spectrum.size(); i++)
        lhs_spectrum[i] = modular_ntt::multiply(
            lhs_spectrum[i], rhs_spectrum[i], plan.modulus());
    plan.inverse(coefficients, lhs_spectrum);
}

template <class GLP>
inline bool polynomialMultiplyNTT(GLPolynomial<GLP> &result,
                                  const GLPolynomial<GLP> &lhs,
                                  const GLPolynomial<GLP> &rhs)
{
    if constexpr (!supportsWidePrimeNTT<GLP>) return false;

    using T = typename GLP::T;
    const std::uint32_t prime_count =
        polynomialMultiplyNTTPrimeCount<GLP>(lhs, rhs);
    if (prime_count == 0) {
        constexpr std::size_t maximum_terms =
            static_cast<std::size_t>(GLP::matrix_dimension) * 2 *
            GLP::matrix_dimension * (2 * GLP::phi - 1);
        constexpr std::uint32_t convolution_bits = ceilLog2(maximum_terms);
        constexpr std::uint32_t crt_safe_bits = 126;
        const std::uint32_t lhs_bits = maxSignedTorusBitWidth<GLP>(lhs);
        const std::uint32_t rhs_bits = maxSignedTorusBitWidth<GLP>(rhs);
        const bool split_lhs = lhs_bits >= rhs_bits;
        const std::uint32_t bounded_bits = split_lhs ? rhs_bits : lhs_bits;
        if (bounded_bits == 0) {
            clear<GLP>(result);
            return true;
        }
        if (bounded_bits + convolution_bits >= crt_safe_bits) return false;
        const std::uint32_t chunk_bits =
            crt_safe_bits - bounded_bits - convolution_bits;
        const auto &wide = split_lhs ? lhs : rhs;
        const auto &bounded = split_lhs ? rhs : lhs;
        const std::uint32_t wide_bits = maxUnsignedTorusBitWidth<GLP>(wide);
        auto chunk = std::make_unique<GLPolynomial<GLP>>();
        auto chunk_product = std::make_unique<GLPolynomial<GLP>>();
        clear<GLP>(result);
        for (std::uint32_t shift = 0; shift < wide_bits; shift += chunk_bits) {
            bool nonzero = false;
            for (std::size_t y = 0; y < GLP::matrix_dimension; y++)
                for (std::size_t i = 0; i < GLP::baseP::n; i++) {
                    (*chunk)[y][i] =
                        unsignedChunk(wide[y][i], shift, chunk_bits);
                    nonzero = nonzero || (*chunk)[y][i] != typename GLP::T{0};
                }
            if (!nonzero) continue;
            if (!polynomialMultiplyNTT<GLP>(*chunk_product, *chunk, bounded))
                return false;
            for (std::size_t y = 0; y < GLP::matrix_dimension; y++)
                for (std::size_t i = 0; i < GLP::baseP::n; i++)
                    result[y][i] += (*chunk_product)[y][i] << shift;
        }
        return true;
    }

    std::vector<std::uint64_t> first_residues;
    polynomialMultiplyNTTResidues<GLP, 0>(first_residues, lhs, rhs);
    if (prime_count == 1) {
        for (std::size_t y = 0; y < GLP::matrix_dimension; y++)
            for (std::size_t i = 0; i < GLP::baseP::n; i++)
                result[y][i] =
                    signedI128ToTorus<T>(modular_ntt::centeredResidue(
                        first_residues[y * GLP::baseP::n + i],
                        modular_ntt::wide_primes[0].value));
        return true;
    }

    std::vector<std::uint64_t> second_residues;
    polynomialMultiplyNTTResidues<GLP, 1>(second_residues, lhs, rhs);
    if (prime_count == 2) {
        static const modular_ntt::TwoPrimeCRT crt(modular_ntt::wide_primes[0],
                                                  modular_ntt::wide_primes[1]);
        for (std::size_t y = 0; y < GLP::matrix_dimension; y++)
            for (std::size_t i = 0; i < GLP::baseP::n; i++) {
                const std::size_t index = y * GLP::baseP::n + i;
                result[y][i] = signedI128ToTorus<T>(crt.reconstructSigned(
                    first_residues[index], second_residues[index]));
            }
        return true;
    }

    std::vector<std::uint64_t> third_residues;
    polynomialMultiplyNTTResidues<GLP, 2>(third_residues, lhs, rhs);
    static const modular_ntt::ThreePrimeCRT crt(modular_ntt::wide_primes[0],
                                                modular_ntt::wide_primes[1],
                                                modular_ntt::wide_primes[2]);
    for (std::size_t y = 0; y < GLP::matrix_dimension; y++)
        for (std::size_t i = 0; i < GLP::baseP::n; i++) {
            const std::size_t index = y * GLP::baseP::n + i;
            result[y][i] = signedI128ToTorus<T>(crt.reconstructSignedBounded(
                first_residues[index], second_residues[index],
                third_residues[index]));
        }
    return true;
}

template <class GLP>
inline void polynomialMultiply(GLPolynomial<GLP> &result,
                               const GLPolynomial<GLP> &lhs,
                               const GLPolynomial<GLP> &rhs)
{
    if (!polynomialMultiplyNTT<GLP>(result, lhs, rhs))
        polynomialMultiplyReference<GLP>(result, lhs, rhs);
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

// A compact complex scalar used by the W^0-only trace kernel below.  The
// third entry saves one addition in the three-multiply complex product.
struct TraceSmallComplexMultiplier {
    std::int64_t real;
    std::int64_t imag;
    std::int64_t real_plus_imag;
};

// Exact active-level trace by a pre-transformed, signed-small complex matrix.
// The generic traceProduct has to scan both dense GL polynomials and performs
// a full-width torus multiplication for every pair.  StC's X plaintext is
// supported only at W^0 and its coefficients are at most a few bits wider
// than its plaintext scale.  Compact the active ciphertext level to 128 bits,
// retain the input rows in cache across every output row, and use Gauss's
// three-multiply complex product.  All arithmetic wraps modulo 2^128, which
// is exact modulo the requested (strictly smaller) active level.
template <class GLP, std::uint32_t LogQ>
inline void traceProductSmallComplex(
    GLCiphertextData<GLP> &result, const GLCiphertextData<GLP> &lhs,
    const std::span<const TraceSmallComplexMultiplier> multipliers)
{
    using T = typename GLP::T;
    constexpr std::size_t n = GLP::matrix_dimension;
    constexpr std::size_t phi = GLP::phi;
    static_assert(LogQ > 0 && LogQ < 128);
    static_assert(LogQ <= std::numeric_limits<T>::digits);
    if (multipliers.size() != n * n)
        throw std::invalid_argument(
            "GL small-complex trace matrix has the wrong size");

    constexpr __uint128_t active_mask =
        (__uint128_t{1} << LogQ) - __uint128_t{1};
    constexpr std::size_t component_stride = phi * n * n * 2;
    constexpr std::size_t w_stride = n * n * 2;
    constexpr std::size_t x_stride = n * 2;
    std::vector<__uint128_t> active_lhs(2 * component_stride);

#pragma omp parallel for collapse(4) schedule(static)
    for (std::size_t component = 0; component < 2; component++)
        for (std::size_t w = 0; w < phi; w++)
            for (std::size_t x = 0; x < n; x++)
                for (std::size_t z = 0; z < n; z++) {
                    const std::size_t destination =
                        component * component_stride + w * w_stride +
                        x * x_stride + z * 2;
                    active_lhs[destination] =
                        torusLowU128(
                            lhs[component][z][baseIndex<GLP>(0, x, w)]) &
                        active_mask;
                    active_lhs[destination + 1] =
                        torusLowU128(
                            lhs[component][z][baseIndex<GLP>(1, x, w)]) &
                        active_mask;
                }

                // Each (W,X) input line is only 32 KiB for the n512 ciphertext
                // and stays hot while all n output rows consume it.
#pragma omp parallel for collapse(2) schedule(static)
    for (std::size_t w = 0; w < phi; w++) {
        for (std::size_t x = 0; x < n; x++) {
            const __uint128_t *input_lines[2] = {
                active_lhs.data() + w * w_stride + x * x_stride,
                active_lhs.data() + component_stride + w * w_stride +
                    x * x_stride};
            for (std::size_t y = 0; y < n; y++) {
                const auto *matrix_line = multipliers.data() + y * n;
                __uint128_t real_sum[2]{};
                __uint128_t imag_sum[2]{};
                for (std::size_t z = 0; z < n; z++) {
                    const auto &scalar = matrix_line[z];
                    for (std::size_t component = 0; component < 2;
                         component++) {
                        const __uint128_t real = input_lines[component][2 * z];
                        const __uint128_t imag =
                            input_lines[component][2 * z + 1];
                        const __uint128_t real_product =
                            multiplyU128BySignedSmall(real, scalar.real);
                        const __uint128_t imag_product =
                            multiplyU128BySignedSmall(imag, scalar.imag);
                        const __uint128_t sum_product =
                            multiplyU128BySignedSmall(real + imag,
                                                      scalar.real_plus_imag);
                        real_sum[component] += real_product - imag_product;
                        imag_sum[component] +=
                            sum_product - real_product - imag_product;
                    }
                }
                for (std::size_t component = 0; component < 2; component++) {
                    result[component][y][baseIndex<GLP>(0, x, w)] =
                        torusFromLowU128<T>(real_sum[component] & active_mask);
                    result[component][y][baseIndex<GLP>(1, x, w)] =
                        torusFromLowU128<T>(imag_sum[component] & active_mask);
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
    constexpr std::uint32_t n = GLP::matrix_dimension;
    constexpr std::uint32_t four_n = 4 * n;
    constexpr std::uint32_t p = GLP::cyclotomic_order;
    if (x_multiplier % four_n == 1 && w_multiplier % p == 1) {
        result = input;
        return;
    }
    clear<GLP>(result);

    std::array<std::uint32_t, GLP::phi> mapped_ws{};
    for (std::uint32_t w = 0; w < GLP::phi; w++)
        mapped_ws[w] = static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(w_multiplier) * w) % p);
    std::array<std::uint32_t, n> x_powers{};
    std::array<std::uint32_t, n> output_xs{};
    for (std::uint32_t x = 0; x < n; x++) {
        const std::uint32_t mapped_x = static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(x_multiplier) * x) % four_n);
        x_powers[x] = mapped_x / n;
        output_xs[x] = mapped_x % n;
    }
    for (std::uint32_t w = 0; w < GLP::phi; w++) {
        const std::uint32_t mapped_w = mapped_ws[w];
        for (std::uint32_t x = 0; x < n; x++) {
            const std::uint32_t x_power = x_powers[x];
            const std::uint32_t output_x = output_xs[x];
            for (std::uint32_t gaussian = 0; gaussian < 2; gaussian++) {
                const auto value = input[baseIndex<GLP>(gaussian, x, w)];
                if (value == typename GLP::T{0}) continue;
                const std::uint32_t gaussian_power =
                    ((x_multiplier & 3U) * gaussian + x_power) & 3U;
                const std::uint32_t output_gaussian = gaussian_power & 1U;
                const bool negative = gaussian_power >= 2;
                if (mapped_w < GLP::phi) {
                    addSigned(result[baseIndex<GLP>(output_gaussian, output_x,
                                                    mapped_w)],
                              value, negative);
                }
                else {
                    // W^(p-1) = -(1 + W + ... + W^(p-2)).
                    for (std::uint32_t reduced_w = 0; reduced_w < GLP::phi;
                         reduced_w++)
                        addSigned(result[baseIndex<GLP>(output_gaussian,
                                                        output_x, reduced_w)],
                                  value, !negative);
                }
            }
        }
    }
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

// Invert one X canonical embedding.  The GL slots use roots
// zeta^(5^j), where zeta is primitive of order 4n.  Dividing by zeta maps
// them to every n-th root in the permutation (5^j-1)/4, so one ordinary
// radix-2 inverse DFT plus a coefficient twist evaluates the same sum as the
// reference formula in O(n log n).
template <class GLP>
class GLInverseXEmbeddingPlan {
public:
    static constexpr std::size_t dimension = GLP::matrix_dimension;
    using Complex = std::complex<long double>;

    GLInverseXEmbeddingPlan()
    {
        constexpr long double pi = 3.141592653589793238462643383279502884L;
        for (std::size_t slot = 0; slot < dimension; slot++) {
            const std::uint32_t exponent =
                powMod(5, static_cast<std::uint32_t>(slot), 4 * dimension);
            frequency_[slot] = (exponent - 1) / 4;
        }
        for (std::size_t coefficient = 0; coefficient < dimension;
             coefficient++) {
            const long double angle =
                -2.0L * pi * coefficient / (4 * dimension);
            inverse_twist_[coefficient] = {std::cos(angle), std::sin(angle)};
        }
        for (std::size_t length = 2; length <= dimension; length <<= 1) {
            const long double angle = -2.0L * pi / length;
            stage_roots_.emplace_back(std::cos(angle), std::sin(angle));
        }
    }

    void apply(const std::span<Complex> output,
               const std::span<const Complex> input,
               std::vector<Complex> &work) const
    {
        if (output.size() != dimension || input.size() != dimension)
            throw std::invalid_argument(
                "GL X embedding transform has the wrong size");
        work.resize(dimension);
        for (std::size_t slot = 0; slot < dimension; slot++)
            work[frequency_[slot]] = input[slot];

        for (std::size_t i = 1, j = 0; i < dimension; i++) {
            std::size_t bit = dimension >> 1;
            for (; (j & bit) != 0; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap(work[i], work[j]);
        }
        std::size_t stage = 0;
        for (std::size_t length = 2; length <= dimension;
             length <<= 1, stage++) {
            const std::size_t half = length >> 1;
            for (std::size_t block = 0; block < dimension; block += length) {
                Complex twiddle = 1;
                for (std::size_t j = 0; j < half; j++) {
                    const Complex even = work[block + j];
                    const Complex odd = work[block + j + half] * twiddle;
                    work[block + j] = even + odd;
                    work[block + j + half] = even - odd;
                    twiddle *= stage_roots_[stage];
                }
            }
        }
        const long double inverse_dimension =
            1.0L / static_cast<long double>(dimension);
        for (std::size_t coefficient = 0; coefficient < dimension;
             coefficient++)
            output[coefficient] = work[coefficient] *
                                  inverse_twist_[coefficient] *
                                  inverse_dimension;
    }

private:
    std::array<std::uint32_t, dimension> frequency_{};
    std::array<Complex, dimension> inverse_twist_{};
    std::vector<Complex> stage_roots_{};
};

template <class GLP>
inline const GLInverseXEmbeddingPlan<GLP> &inverseXEmbeddingPlan()
{
    static const GLInverseXEmbeddingPlan<GLP> plan;
    return plan;
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

    std::vector<std::complex<long double>> w_values(
        static_cast<std::size_t>(n) * n * phi);
    auto w_value = [&](const std::uint32_t row, const std::uint32_t column,
                       const std::uint32_t w) -> std::complex<long double> & {
        return w_values[(static_cast<std::size_t>(row) * n + column) * phi + w];
    };

    // Invert the W canonical embedding.  Values are supplied at all nonzero
    // p-th roots.  Choosing the degree-(p-2) representative determines the
    // missing value at W=1.  Precompute the tiny phi-by-phi transform so the
    // n^2 matrix entries do not repeatedly evaluate trigonometric functions.
    std::array<std::array<std::complex<long double>, phi>, phi> w_weights{};
    for (std::uint32_t batch = 0; batch < phi; batch++) {
        const std::uint32_t exponent = gl_detail::batchExponent<GLP>(batch);
        const auto at_one_weight = gl_detail::wRoot<GLP>(exponent);
        for (std::uint32_t w = 0; w < phi; w++) {
            const std::uint32_t root_exponent =
                (p - static_cast<std::uint32_t>(
                         (static_cast<std::uint64_t>(exponent) * w) % p)) %
                p;
            w_weights[batch][w] =
                (gl_detail::wRoot<GLP>(root_exponent) - at_one_weight) /
                static_cast<long double>(p);
        }
    }
#pragma omp parallel for collapse(2) schedule(static)
    for (std::uint32_t row = 0; row < n; row++) {
        for (std::uint32_t column = 0; column < n; column++) {
            for (std::uint32_t w = 0; w < phi; w++) {
                std::complex<long double> coefficient = 0;
                for (std::uint32_t batch = 0; batch < phi; batch++)
                    coefficient += static_cast<std::complex<long double>>(
                                       matrices(batch, row, column)) *
                                   w_weights[batch][w];
                w_value(row, column, w) = coefficient;
            }
        }
    }

    gl_detail::clear<GLP>(plaintext.poly);
    const auto &x_plan = gl_detail::inverseXEmbeddingPlan<GLP>();
#pragma omp parallel
    {
        std::vector<std::complex<long double>> after_columns(
            static_cast<std::size_t>(n) * n);
        std::vector<std::complex<long double>> input_line(n);
        std::vector<std::complex<long double>> output_line(n);
        std::vector<std::complex<long double>> work;
#pragma omp for schedule(dynamic)
        for (std::uint32_t w = 0; w < phi; w++) {
            for (std::uint32_t row = 0; row < n; row++) {
                for (std::uint32_t column = 0; column < n; column++)
                    input_line[column] = w_value(row, column, w);
                x_plan.apply(output_line, input_line, work);
                for (std::uint32_t y = 0; y < n; y++)
                    after_columns[static_cast<std::size_t>(row) * n + y] =
                        output_line[y];
            }
            for (std::uint32_t y = 0; y < n; y++) {
                for (std::uint32_t row = 0; row < n; row++)
                    input_line[row] =
                        after_columns[static_cast<std::size_t>(row) * n + y];
                x_plan.apply(output_line, input_line, work);
                for (std::uint32_t x = 0; x < n; x++) {
                    const auto coefficient = output_line[x] * scale;
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

template <class GLP, std::uint32_t Bits>
inline void packDigitPolynomial(GLPackedBasePolynomial<GLP, Bits> &result,
                                const GLBasePolynomial<GLP> &input)
{
    using Digit = GLSignedDigit<Bits>;
    using P = typename GLP::baseP;
    constexpr std::int64_t minimum =
        static_cast<std::int64_t>(std::numeric_limits<Digit>::min());
    constexpr std::int64_t maximum =
        static_cast<std::int64_t>(std::numeric_limits<Digit>::max());
    for (std::uint32_t coefficient = 0; coefficient < P::n; coefficient++) {
        const std::int64_t value = torusToSignedSmall<P>(input[coefficient]);
        if (value < minimum || value > maximum)
            throw std::overflow_error("GL DD digit exceeds packed storage");
        result[coefficient] = static_cast<Digit>(value);
    }
}

template <class GLP, std::uint32_t Bits>
inline void packDigitPolynomial(GLPackedPolynomial<GLP, Bits> &result,
                                const GLPolynomial<GLP> &input)
{
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++)
        packDigitPolynomial<GLP, Bits>(result[y], input[y]);
}

template <class GLP, std::uint32_t Bits>
inline void unpackDigitPolynomial(
    GLBasePolynomial<GLP> &result,
    const GLPackedBasePolynomial<GLP, Bits> &input)
{
    using P = typename GLP::baseP;
    using T = typename P::T;
    for (std::uint32_t coefficient = 0; coefficient < P::n; coefficient++) {
        const std::int64_t value = input[coefficient];
        if constexpr (is_multilimb_uint_v<T>)
            result[coefficient] = T::from_signed_i64(value);
        else {
            static_assert(std::is_same_v<T, __uint128_t>);
            result[coefficient] =
                static_cast<__uint128_t>(static_cast<__int128_t>(value));
        }
    }
}

template <class GLP, std::uint32_t Bits>
inline void unpackDigitPolynomial(GLPolynomial<GLP> &result,
                                  const GLPackedPolynomial<GLP, Bits> &input)
{
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++)
        unpackDigitPolynomial<GLP, Bits>(result[y], input[y]);
}

// A small DD switch multiplies every primary input digit by every packed key
// row.  All of those products are accumulated before the Bbar recomposition,
// so the primary-row sum can stay in the NTT domain.  The bound includes that
// extra sum; this is stricter than the bound for one baseMultiply call.
template <class GLP, class SwitchKey>
inline constexpr std::uint32_t smallKeySwitchAccumulationNTTPrimeCount = [] {
    if constexpr (!supportsWidePrimeNTT<GLP>) return 0U;

    constexpr std::size_t maximum_terms =
        2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
    constexpr std::uint32_t required_bits =
        SwitchKey::primary_bit + SwitchKey::bbar_bit + ceilLog2(maximum_terms) +
        ceilLog2(SwitchKey::primary_rows);
    if constexpr (required_bits <= 60)
        return 1U;
    else if constexpr (required_bits <= 122)
        return 2U;
    else
        return 0U;
}();

template <class GLP, class SwitchKey>
inline constexpr std::uint64_t smallKeySwitchNTTKeyCacheBytes =
    static_cast<std::uint64_t>(
        smallKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>) *
    SwitchKey::primary_rows * SwitchKey::bbar_rows * 2 *
    GLBaseNTTPlan<GLP>::coefficient_count * sizeof(std::uint64_t);

template <class GLP>
inline constexpr std::size_t smallKeySwitchNTTCoefficientBlockSize =
    GLBaseNTTPlan<GLP>::coefficient_count < 4096
        ? GLBaseNTTPlan<GLP>::coefficient_count
        : 4096;

template <class GLP, class SwitchKey>
inline constexpr std::size_t smallKeySwitchBlockedSpectrumOffset(
    const std::uint32_t primary, const std::uint32_t bbar,
    const std::size_t component, const std::size_t block)
{
    constexpr std::size_t coefficient_count =
        GLBaseNTTPlan<GLP>::coefficient_count;
    constexpr std::size_t block_size =
        smallKeySwitchNTTCoefficientBlockSize<GLP>;
    static_assert(coefficient_count % block_size == 0);
    constexpr std::size_t block_count = coefficient_count / block_size;
    return ((((static_cast<std::size_t>(bbar) * 2 + component) * block_count +
               block) *
                  SwitchKey::primary_rows +
              primary) *
             block_size);
}

template <bool CoefficientBlocked, class GLP, class SwitchKey>
inline std::shared_ptr<typename SwitchKey::TransientNTTCache>
prepareSmallKeySwitchNTTCacheLayout(const SwitchKey &switch_key)
{
    constexpr std::uint32_t prime_count =
        smallKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    static_assert(prime_count != 0,
                  "the requested GL DD switch has no exact NTT path");
    constexpr std::size_t coefficient_count =
        GLBaseNTTPlan<GLP>::coefficient_count;
    constexpr std::size_t key_row_count =
        static_cast<std::size_t>(SwitchKey::primary_rows) *
        SwitchKey::bbar_rows * 2;
    const auto key_cache = switch_key.transient_ntt_cache;
    auto &initialize_once = [&]() -> std::once_flag & {
        if constexpr (CoefficientBlocked)
            return key_cache->blocked_initialize_once;
        else
            return key_cache->initialize_once;
    }();
    auto &key_spectra = [&]()
        -> std::array<std::vector<std::uint64_t>, 2> & {
        if constexpr (CoefficientBlocked)
            return key_cache->blocked_spectra;
        else
            return key_cache->spectra;
    }();
    std::call_once(initialize_once, [&] {
        key_spectra[0].resize(key_row_count * coefficient_count);
        if constexpr (prime_count == 2)
            key_spectra[1].resize(key_row_count * coefficient_count);

        const auto &first_plan = baseNTTPlan<GLP, 0>();
        const GLBaseNTTPlan<GLP> *second_plan = nullptr;
        if constexpr (prime_count == 2) second_plan = &baseNTTPlan<GLP, 1>();
#pragma omp parallel
        {
            GLBasePolynomial<GLP> key_row{};
            std::array<std::vector<std::uint64_t>, 2> transformed{};
            if constexpr (CoefficientBlocked) {
                transformed[0].resize(coefficient_count);
                if constexpr (prime_count == 2)
                    transformed[1].resize(coefficient_count);
            }
#pragma omp for collapse(3) schedule(static)
            for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
                 primary++) {
                for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows;
                     bbar++) {
                    for (std::size_t component = 0; component < 2;
                         component++) {
                        unpackDigitPolynomial<GLP, SwitchKey::bbar_bit>(
                            key_row, switch_key.at(primary, bbar)[component]);
                        const std::size_t row =
                            (static_cast<std::size_t>(primary) *
                                 SwitchKey::bbar_rows +
                             bbar) *
                                2 +
                            component;
                        if constexpr (CoefficientBlocked) {
                            constexpr std::size_t block_size =
                                smallKeySwitchNTTCoefficientBlockSize<GLP>;
                            constexpr std::size_t block_count =
                                coefficient_count / block_size;
                            first_plan.forward(
                                std::span<std::uint64_t>(transformed[0]),
                                key_row);
                            if constexpr (prime_count == 2)
                                second_plan->forward(
                                    std::span<std::uint64_t>(transformed[1]),
                                    key_row);
                            for (std::size_t block = 0; block < block_count;
                                 block++) {
                                const std::size_t destination =
                                    smallKeySwitchBlockedSpectrumOffset<
                                        GLP, SwitchKey>(primary, bbar,
                                                        component, block);
                                std::copy_n(
                                    transformed[0].data() + block * block_size,
                                    block_size,
                                    key_spectra[0].data() + destination);
                                if constexpr (prime_count == 2)
                                    std::copy_n(
                                        transformed[1].data() +
                                            block * block_size,
                                        block_size,
                                        key_spectra[1].data() + destination);
                            }
                        }
                        else {
                            first_plan.forward(
                                std::span<std::uint64_t>(
                                    key_spectra[0].data() +
                                        row * coefficient_count,
                                    coefficient_count),
                                key_row);
                            if constexpr (prime_count == 2)
                                second_plan->forward(
                                    std::span<std::uint64_t>(
                                        key_spectra[1].data() +
                                            row * coefficient_count,
                                        coefficient_count),
                                    key_row);
                        }
                    }
                }
            }
        }
    });
    return key_cache;
}

template <class GLP, class SwitchKey>
inline std::shared_ptr<typename SwitchKey::TransientNTTCache>
prepareSmallKeySwitchNTTCache(const SwitchKey &switch_key)
{
    return prepareSmallKeySwitchNTTCacheLayout<false, GLP, SwitchKey>(
        switch_key);
}

template <class GLP, class SwitchKey>
inline std::shared_ptr<typename SwitchKey::TransientNTTCache>
prepareSmallKeySwitchBlockedNTTCache(const SwitchKey &switch_key)
{
    return prepareSmallKeySwitchNTTCacheLayout<true, GLP, SwitchKey>(
        switch_key);
}

template <class SwitchKey>
inline void releaseSmallKeySwitchNTTCache(const SwitchKey &switch_key)
{
    switch_key.transient_ntt_cache =
        std::make_shared<typename SwitchKey::TransientNTTCache>();
}

// PrepareSlice(slice) materializes its primary decomposition, GetInput returns
// one of those base polynomials, and StoreOutput consumes one reconstructed
// coefficient.  Key spectra are transient and never serialized.  Processing
// one Y slice at a time also avoids caching all input spectra for a
// production-size GL polynomial.
template <class GLP, class SwitchKey, class PrepareSlice, class GetInput,
          class StoreOutput>
inline bool accumulateSmallKeySwitchProductsNTT(const std::size_t slice_count,
                                                const SwitchKey &switch_key,
                                                PrepareSlice &&prepare_slice,
                                                GetInput &&get_input,
                                                StoreOutput &&store_output)
{
    using T = typename GLP::T;
    constexpr std::uint32_t prime_count =
        smallKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    if constexpr (prime_count == 0) {
        return false;
    }
    else {
        constexpr std::size_t coefficient_count =
            GLBaseNTTPlan<GLP>::coefficient_count;
        constexpr std::size_t input_row_count = SwitchKey::primary_rows;

        const auto &first_plan = baseNTTPlan<GLP, 0>();
        const GLBaseNTTPlan<GLP> *second_plan = nullptr;
        if constexpr (prime_count == 2) second_plan = &baseNTTPlan<GLP, 1>();
        const std::uint64_t first_modulus = first_plan.modulus();
        const std::uint64_t second_modulus =
            prime_count == 2 ? second_plan->modulus() : 0;
        const auto key_cache =
            prepareSmallKeySwitchNTTCache<GLP, SwitchKey>(switch_key);
        const auto &key_spectra = key_cache->spectra;

        std::array<std::vector<std::uint64_t>, 2> input_spectra;
        std::array<std::vector<std::uint64_t>, 2> accumulators;
        std::array<std::vector<std::uint64_t>, 2> coefficients;
        input_spectra[0].resize(input_row_count * coefficient_count);
        accumulators[0].resize(coefficient_count);
        coefficients[0].resize(coefficient_count);
        if constexpr (prime_count == 2) {
            input_spectra[1].resize(input_row_count * coefficient_count);
            accumulators[1].resize(coefficient_count);
            coefficients[1].resize(coefficient_count);
        }

        for (std::size_t slice = 0; slice < slice_count; slice++) {
            prepare_slice(slice);
            for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
                 primary++) {
                first_plan.forward(std::span<std::uint64_t>(
                                       input_spectra[0].data() +
                                           static_cast<std::size_t>(primary) *
                                               coefficient_count,
                                       coefficient_count),
                                   get_input(primary, slice));
                if constexpr (prime_count == 2) {
                    second_plan->forward(
                        std::span<std::uint64_t>(
                            input_spectra[1].data() +
                                static_cast<std::size_t>(primary) *
                                    coefficient_count,
                            coefficient_count),
                        get_input(primary, slice));
                }
            }

            for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++) {
                for (std::size_t component = 0; component < 2; component++) {
                    std::fill(accumulators[0].begin(), accumulators[0].end(),
                              0);
                    if constexpr (prime_count == 2)
                        std::fill(accumulators[1].begin(),
                                  accumulators[1].end(), 0);

                    for (std::uint32_t primary = 0;
                         primary < SwitchKey::primary_rows; primary++) {
                        const std::size_t row =
                            (static_cast<std::size_t>(primary) *
                                 SwitchKey::bbar_rows +
                             bbar) *
                                2 +
                            component;
                        const std::uint64_t *first_input =
                            input_spectra[0].data() +
                            static_cast<std::size_t>(primary) *
                                coefficient_count;
                        const std::uint64_t *first_key =
                            key_spectra[0].data() + row * coefficient_count;
                        for (std::size_t i = 0; i < coefficient_count; i++) {
                            const std::uint64_t product = modular_ntt::multiply(
                                first_input[i], first_key[i], first_modulus);
                            accumulators[0][i] = modular_ntt::add(
                                accumulators[0][i], product, first_modulus);
                        }

                        if constexpr (prime_count == 2) {
                            const std::uint64_t *second_input =
                                input_spectra[1].data() +
                                static_cast<std::size_t>(primary) *
                                    coefficient_count;
                            const std::uint64_t *second_key =
                                key_spectra[1].data() + row * coefficient_count;
                            for (std::size_t i = 0; i < coefficient_count;
                                 i++) {
                                const std::uint64_t product =
                                    modular_ntt::multiply(second_input[i],
                                                          second_key[i],
                                                          second_modulus);
                                accumulators[1][i] =
                                    modular_ntt::add(accumulators[1][i],
                                                     product, second_modulus);
                            }
                        }
                    }

                    first_plan.inverse(
                        std::span<std::uint64_t>(coefficients[0]),
                        std::span<std::uint64_t>(accumulators[0]));
                    if constexpr (prime_count == 1) {
                        for (std::size_t i = 0; i < coefficient_count; i++)
                            store_output(
                                component, bbar, slice, i,
                                signedI128ToTorus<T>(
                                    modular_ntt::centeredResidue(
                                        coefficients[0][i], first_modulus)));
                    }
                    else {
                        second_plan->inverse(
                            std::span<std::uint64_t>(coefficients[1]),
                            std::span<std::uint64_t>(accumulators[1]));
                        static const modular_ntt::TwoPrimeCRT crt(
                            modular_ntt::wide_primes[0],
                            modular_ntt::wide_primes[1]);
                        for (std::size_t i = 0; i < coefficient_count; i++)
                            store_output(
                                component, bbar, slice, i,
                                signedI128ToTorus<T>(crt.reconstructSigned(
                                    coefficients[0][i], coefficients[1][i])));
                    }
                }
            }
        }
        return true;
    }
}

// Sum several raw base-ring DD switches before reconstruction.  HMux uses
// this to combine all body/mask radix branches under P*Q and therefore pays
// for one set of inverse transforms instead of one set per branch.  Count is
// part of the exact CRT bound; no probabilistic wraparound assumption is used.
//
// The usual half-open balanced digit interval is not exactly equivariant
// under coefficient negation at the -B/2 boundary.  HMux automorphisms need
// that property in order to hoist decomposition across signed coefficient
// permutations.  This variant uses the symmetric tie convention +B/2 for a
// positive value and -B/2 for its negative.  It has the same magnitude bound
// and reconstructs the same active-level torus value.
template <class P, std::uint32_t LogQ, std::uint32_t BaseBit,
          std::size_t RowCount>
inline void activeSymmetricBaseDecomposePolynomialRows(
    std::array<Polynomial<P>, RowCount> &result, const Polynomial<P> &input)
{
    static_assert(LogQ > 0 && LogQ <= ckks_detail::torus_width_v<P>);
    static_assert(BaseBit > 0);
    static_assert(BaseBit < ckks_detail::torus_width_v<P>);
    static_assert(RowCount == (LogQ + BaseBit - 1) / BaseBit);
    static_assert(RowCount * BaseBit <= ckks_detail::torus_width_v<P>);

    using T = typename P::T;
    constexpr T half = T{1} << (BaseBit - 1);
    constexpr T mask = (T{1} << BaseBit) - T{1};
    for (std::uint32_t coefficient = 0; coefficient < P::n; coefficient++) {
        const T centered =
            ckks_detail::centeredLevelToTorus<P, LogQ>(input[coefficient]);
        const bool negative =
            (centered & (T{1} << (ckks_detail::torus_width_v<P> - 1))) != T{0};
        T magnitude = negative ? T{0} - centered : centered;
        for (std::size_t reverse = 0; reverse < RowCount; reverse++) {
            const std::size_t row = RowCount - reverse - 1;
            const T raw = magnitude & mask;
            magnitude >>= BaseBit;
            const bool carry = raw > half;
            T digit_magnitude = carry ? (T{1} << BaseBit) - raw : raw;
            if (carry) magnitude += T{1};
            const bool digit_negative =
                digit_magnitude != T{0} && (negative != carry);
            result[row][coefficient] =
                digit_negative ? T{0} - digit_magnitude : digit_magnitude;
        }
    }
}

template <class GLP, class SwitchKey, std::size_t Count>
struct SmallKeySwitchSumNTTWorkspace {
    using P = typename GLP::baseP;
    static constexpr std::size_t input_row_count = SwitchKey::primary_rows;

    std::array<std::vector<std::uint64_t>, 2> input_spectra{};
    std::array<std::vector<std::uint64_t>, 2> source_spectra{};
    std::array<std::vector<std::uint64_t>, 2> accumulators{};
    std::array<std::vector<std::uint64_t>, 2> coefficients{};
#ifdef USE_HEXL
    std::vector<std::uint64_t> product{};
    // The AVX-512DQ batch kernel wins with the physical-core HMux profile but
    // can become compute-bound with a larger SMT team, so callers opt in.
    bool use_batched_vector_mac = false;
    // Coefficient blocks keep the primary key rows needed by one HMux MAC
    // close together.  This layout is useful only with the batched kernel.
    bool use_coefficient_blocked_key_layout = false;
#else
    std::array<std::vector<unsigned __int128>, 2> wide_accumulators{};
#endif
    std::unique_ptr<std::array<Polynomial<P>, input_row_count>> input_digits{};
};

template <class GLP, class SwitchKey, std::size_t Count>
inline void accumulateSmallKeySwitchPreparedSumNTT(
    GLBaseCiphertextData<GLP> &result,
    const std::array<std::shared_ptr<typename SwitchKey::TransientNTTCache>,
                     Count> &key_caches,
    SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey, Count> &workspace)
{
    using T = typename GLP::T;
    constexpr std::uint32_t prime_count =
        smallKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    constexpr std::size_t coefficient_count =
        GLBaseNTTPlan<GLP>::coefficient_count;
    constexpr std::size_t input_row_count = SwitchKey::primary_rows;
    const auto &first_plan = baseNTTPlan<GLP, 0>();
    const GLBaseNTTPlan<GLP> *second_plan = nullptr;
    if constexpr (prime_count == 2) second_plan = &baseNTTPlan<GLP, 1>();
    const std::uint64_t first_modulus = first_plan.modulus();
    const std::uint64_t second_modulus =
        prime_count == 2 ? second_plan->modulus() : 0;
    const auto &input_spectra = workspace.input_spectra;

    clear<GLP>(result[0]);
    clear<GLP>(result[1]);
    auto &accumulators = workspace.accumulators;
    auto &coefficients = workspace.coefficients;
#ifdef USE_HEXL
    workspace.product.resize(coefficient_count);
#else
    auto &wide_accumulators = workspace.wide_accumulators;
#endif
    accumulators[0].resize(coefficient_count);
    coefficients[0].resize(coefficient_count);
#ifndef USE_HEXL
    wide_accumulators[0].resize(coefficient_count);
#endif
    if constexpr (prime_count == 2) {
        accumulators[1].resize(coefficient_count);
        coefficients[1].resize(coefficient_count);
#ifndef USE_HEXL
        wide_accumulators[1].resize(coefficient_count);
#endif
    }

    for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++) {
        const std::uint32_t shift =
            (SwitchKey::bbar_rows - bbar - 1) * SwitchKey::bbar_bit;
        for (std::size_t component = 0; component < 2; component++) {
            const auto accumulate_prime = [&](const std::size_t prime,
                                              const std::uint64_t modulus) {
                auto &accumulator = accumulators[prime];
                std::fill(accumulator.begin(), accumulator.end(), 0);
#ifdef USE_HEXL
                if (workspace.use_coefficient_blocked_key_layout) {
                    if (!workspace.use_batched_vector_mac ||
                        !modular_ntt::hasFastVectorMultiplyAddBatch)
                        throw std::logic_error(
                            "blocked GL DD spectra require the fast batched "
                            "modular MAC");
                    constexpr std::size_t block_size =
                        smallKeySwitchNTTCoefficientBlockSize<GLP>;
                    constexpr std::size_t block_count =
                        coefficient_count / block_size;
                    constexpr std::size_t product_count =
                        Count * input_row_count;
                    std::array<const std::uint64_t *, product_count>
                        input_pointers{};
                    std::array<const std::uint64_t *, product_count>
                        key_pointers{};
                    for (std::size_t block = 0; block < block_count; block++) {
                        std::size_t product_index = 0;
                        for (std::size_t term = 0; term < Count; term++) {
                            const auto &key_spectra =
                                key_caches[term]->blocked_spectra[prime];
                            for (std::uint32_t primary = 0;
                                 primary < SwitchKey::primary_rows;
                                 primary++) {
                                const std::size_t input_row =
                                    term * input_row_count + primary;
                                input_pointers[product_index] =
                                    input_spectra[prime].data() +
                                    input_row * coefficient_count +
                                    block * block_size;
                                key_pointers[product_index] =
                                    key_spectra.data() +
                                    smallKeySwitchBlockedSpectrumOffset<
                                        GLP, SwitchKey>(primary, bbar,
                                                        component, block);
                                product_index++;
                            }
                        }
                        modular_ntt::multiplyAddVectorBatch(
                            accumulator.data() + block * block_size,
                            input_pointers.data(), key_pointers.data(),
                            product_count, block_size, modulus);
                    }
                    return;
                }
#endif
#ifndef USE_HEXL
                auto &wide = wide_accumulators[prime];
                std::fill(wide.begin(), wide.end(), 0);
                std::size_t batch_count = 0;
                const auto flush = [&] {
                    for (std::size_t i = 0; i < coefficient_count; i++) {
                        accumulator[i] = modular_ntt::add(
                            accumulator[i],
                            modular_ntt::reduceWide(wide[i], modulus), modulus);
                        wide[i] = 0;
                    }
                    batch_count = 0;
                };
#else
                constexpr std::size_t product_count =
                    Count * input_row_count;
                std::array<const std::uint64_t *, product_count>
                    input_pointers{};
                std::array<const std::uint64_t *, product_count>
                    key_pointers{};
                std::size_t product_index = 0;
#endif
                for (std::size_t term = 0; term < Count; term++) {
                    const auto &key_spectra = key_caches[term]->spectra[prime];
                    for (std::uint32_t primary = 0;
                         primary < SwitchKey::primary_rows; primary++) {
                        const std::size_t input_row =
                            term * input_row_count + primary;
                        const std::size_t key_row =
                            (static_cast<std::size_t>(primary) *
                                 SwitchKey::bbar_rows +
                             bbar) *
                                2 +
                            component;
                        const std::uint64_t *input_spectrum =
                            input_spectra[prime].data() +
                            input_row * coefficient_count;
                        const std::uint64_t *key_spectrum =
                            key_spectra.data() + key_row * coefficient_count;
#ifdef USE_HEXL
                        input_pointers[product_index] = input_spectrum;
                        key_pointers[product_index] = key_spectrum;
                        product_index++;
#else
                        for (std::size_t i = 0; i < coefficient_count; i++)
                            wide[i] += static_cast<unsigned __int128>(
                                           input_spectrum[i]) *
                                       key_spectrum[i];
                        batch_count++;
                        if (batch_count == 16) flush();
#endif
                    }
                }
#ifdef USE_HEXL
                if (workspace.use_batched_vector_mac &&
                    modular_ntt::hasFastVectorMultiplyAddBatch) {
                    modular_ntt::multiplyAddVectorBatch(
                        accumulator.data(), input_pointers.data(),
                        key_pointers.data(), product_count, coefficient_count,
                        modulus);
                }
                else {
                    for (std::size_t product = 0; product < product_count;
                         product++) {
                        intel::hexl::EltwiseMultMod(
                            workspace.product.data(), input_pointers[product],
                            key_pointers[product], coefficient_count, modulus,
                            1);
                        intel::hexl::EltwiseAddMod(
                            accumulator.data(), accumulator.data(),
                            workspace.product.data(), coefficient_count,
                            modulus);
                    }
                }
#else
                if (batch_count != 0) flush();
#endif
            };
            accumulate_prime(0, first_modulus);
            if constexpr (prime_count == 2) accumulate_prime(1, second_modulus);

            first_plan.inverse(std::span<std::uint64_t>(coefficients[0]),
                               std::span<std::uint64_t>(accumulators[0]));
            if constexpr (prime_count == 1) {
                for (std::size_t i = 0; i < coefficient_count; i++)
                    result[component][i] +=
                        signedI128ToTorus<T>(modular_ntt::centeredResidue(
                            coefficients[0][i], first_modulus))
                        << shift;
            }
            else {
                second_plan->inverse(std::span<std::uint64_t>(coefficients[1]),
                                     std::span<std::uint64_t>(accumulators[1]));
                static const modular_ntt::TwoPrimeCRT crt(
                    modular_ntt::wide_primes[0], modular_ntt::wide_primes[1]);
                for (std::size_t i = 0; i < coefficient_count; i++)
                    result[component][i] +=
                        signedI128ToTorus<T>(crt.reconstructSigned(
                            coefficients[0][i], coefficients[1][i]))
                        << shift;
            }
        }
    }
    reduce<GLP, SwitchKey::key_log_q>(result[0]);
    reduce<GLP, SwitchKey::key_log_q>(result[1]);
}

template <class GLP, class SwitchKey, std::size_t Count>
inline bool accumulateSmallKeySwitchSumNTT(
    GLBaseCiphertextData<GLP> &result,
    const std::array<const GLBasePolynomial<GLP> *, Count> &inputs,
    const std::array<const SwitchKey *, Count> &switch_keys,
    SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey, Count> *provided_workspace =
        nullptr)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t prime_count =
        smallKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    constexpr std::size_t maximum_terms =
        2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
    constexpr std::uint32_t required_bits =
        SwitchKey::primary_bit + SwitchKey::bbar_bit + ceilLog2(maximum_terms) +
        ceilLog2(SwitchKey::primary_rows) + ceilLog2(Count);
    constexpr bool exact_bound = (prime_count == 1 && required_bits <= 60) ||
                                 (prime_count == 2 && required_bits <= 122);
    if constexpr (prime_count == 0 || !exact_bound) {
        return false;
    }
    else {
        static_assert(Count > 0);
        constexpr std::size_t coefficient_count =
            GLBaseNTTPlan<GLP>::coefficient_count;
        constexpr std::size_t input_row_count = SwitchKey::primary_rows;
        const auto &first_plan = baseNTTPlan<GLP, 0>();
        const GLBaseNTTPlan<GLP> *second_plan = nullptr;
        if constexpr (prime_count == 2) second_plan = &baseNTTPlan<GLP, 1>();

        std::array<std::shared_ptr<typename SwitchKey::TransientNTTCache>,
                   Count>
            key_caches;
        SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey, Count> local_workspace;
        auto &workspace = provided_workspace == nullptr ? local_workspace
                                                        : *provided_workspace;
        auto &input_spectra = workspace.input_spectra;
        const std::size_t input_spectrum_count =
            Count * input_row_count * coefficient_count;
        input_spectra[0].resize(input_spectrum_count);
        if constexpr (prime_count == 2)
            input_spectra[1].resize(input_spectrum_count);

        if (!workspace.input_digits)
            workspace.input_digits =
                std::make_unique<std::array<Polynomial<P>, input_row_count>>();
        auto &input_digits = *workspace.input_digits;
        for (std::size_t term = 0; term < Count; term++) {
            if (inputs[term] == nullptr || switch_keys[term] == nullptr)
                throw std::invalid_argument("null GL DD switch-sum operand");
            if (switch_keys[term]->data.size() !=
                SwitchKey::primary_rows * SwitchKey::bbar_rows)
                throw std::invalid_argument(
                    "uninitialized GL DD switch-sum key");
#ifdef USE_HEXL
            if (workspace.use_coefficient_blocked_key_layout)
                key_caches[term] =
                    prepareSmallKeySwitchBlockedNTTCache<GLP, SwitchKey>(
                        *switch_keys[term]);
            else
#endif
                key_caches[term] =
                    prepareSmallKeySwitchNTTCache<GLP, SwitchKey>(
                        *switch_keys[term]);
            ckks_detail::activeBaseDecomposePolynomialRows<
                P, SwitchKey::log_q, SwitchKey::primary_bit,
                SwitchKey::primary_rows>(input_digits, *inputs[term]);
            for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
                 primary++) {
                const std::size_t input_row = term * input_row_count + primary;
                first_plan.forward(
                    std::span<std::uint64_t>(
                        input_spectra[0].data() + input_row * coefficient_count,
                        coefficient_count),
                    input_digits[primary]);
                if constexpr (prime_count == 2)
                    second_plan->forward(std::span<std::uint64_t>(
                                             input_spectra[1].data() +
                                                 input_row * coefficient_count,
                                             coefficient_count),
                                         input_digits[primary]);
            }
        }

        accumulateSmallKeySwitchPreparedSumNTT<GLP, SwitchKey, Count>(
            result, key_caches, workspace);
        return true;
    }
}

// Hoisted form used by HMux. Every switch input is an automorphism of one of
// a small number of source polynomials. Symmetric decomposition commutes with
// the automorphism's signed coefficient permutation, and evaluation at an
// automorphed X root is just a permutation of the already-computed spectrum.
template <class GLP, class SwitchKey, std::size_t Count,
          std::size_t SourceCount>
inline bool accumulateSmallKeySwitchAutomorphismSumNTT(
    GLBaseCiphertextData<GLP> &result,
    const std::array<const GLBasePolynomial<GLP> *, SourceCount> &sources,
    const std::array<std::size_t, Count> &source_indices,
    const std::array<std::uint32_t, Count> &x_multipliers,
    const std::array<const SwitchKey *, Count> &switch_keys,
    SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey, Count> *provided_workspace =
        nullptr)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t prime_count =
        smallKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    constexpr std::size_t maximum_terms =
        2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
    constexpr std::uint32_t required_bits =
        SwitchKey::primary_bit + SwitchKey::bbar_bit + ceilLog2(maximum_terms) +
        ceilLog2(SwitchKey::primary_rows) + ceilLog2(Count);
    constexpr bool exact_bound = (prime_count == 1 && required_bits <= 60) ||
                                 (prime_count == 2 && required_bits <= 122);
    if constexpr (prime_count == 0 || !exact_bound) {
        return false;
    }
    else {
        static_assert(Count > 0 && SourceCount > 0);
        constexpr std::size_t coefficient_count =
            GLBaseNTTPlan<GLP>::coefficient_count;
        constexpr std::size_t input_row_count = SwitchKey::primary_rows;
        const auto &first_plan = baseNTTPlan<GLP, 0>();
        const GLBaseNTTPlan<GLP> *second_plan = nullptr;
        if constexpr (prime_count == 2) second_plan = &baseNTTPlan<GLP, 1>();

        SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey, Count> local_workspace;
        auto &workspace = provided_workspace == nullptr ? local_workspace
                                                        : *provided_workspace;
        const std::size_t source_spectrum_count =
            SourceCount * input_row_count * coefficient_count;
        workspace.source_spectra[0].resize(source_spectrum_count);
        if constexpr (prime_count == 2)
            workspace.source_spectra[1].resize(source_spectrum_count);
        const std::size_t input_spectrum_count =
            Count * input_row_count * coefficient_count;
        workspace.input_spectra[0].resize(input_spectrum_count);
        if constexpr (prime_count == 2)
            workspace.input_spectra[1].resize(input_spectrum_count);
        if (!workspace.input_digits)
            workspace.input_digits =
                std::make_unique<std::array<Polynomial<P>, input_row_count>>();
        auto &input_digits = *workspace.input_digits;

        for (std::size_t source = 0; source < SourceCount; source++) {
            if (sources[source] == nullptr)
                throw std::invalid_argument("null GL DD hoisted source");
            activeSymmetricBaseDecomposePolynomialRows<P, SwitchKey::log_q,
                                                       SwitchKey::primary_bit,
                                                       SwitchKey::primary_rows>(
                input_digits, *sources[source]);
            for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
                 primary++) {
                const std::size_t source_row =
                    source * input_row_count + primary;
                first_plan.forward(std::span<std::uint64_t>(
                                       workspace.source_spectra[0].data() +
                                           source_row * coefficient_count,
                                       coefficient_count),
                                   input_digits[primary]);
                if constexpr (prime_count == 2)
                    second_plan->forward(
                        std::span<std::uint64_t>(
                            workspace.source_spectra[1].data() +
                                source_row * coefficient_count,
                            coefficient_count),
                        input_digits[primary]);
            }
        }

        std::array<GLBaseXAutomorphismSpectrumMap<GLP>, Count> z_maps{};
        std::array<std::shared_ptr<typename SwitchKey::TransientNTTCache>,
                   Count>
            key_caches;
        for (std::size_t term = 0; term < Count; term++) {
            if (source_indices[term] >= SourceCount ||
                (x_multipliers[term] & 1U) == 0 || switch_keys[term] == nullptr)
                throw std::invalid_argument(
                    "invalid GL DD hoisted automorphism operand");
            if (switch_keys[term]->data.size() !=
                SwitchKey::primary_rows * SwitchKey::bbar_rows)
                throw std::invalid_argument(
                    "uninitialized GL DD hoisted switch key");
#ifdef USE_HEXL
            if (workspace.use_coefficient_blocked_key_layout)
                key_caches[term] =
                    prepareSmallKeySwitchBlockedNTTCache<GLP, SwitchKey>(
                        *switch_keys[term]);
            else
#endif
                key_caches[term] =
                    prepareSmallKeySwitchNTTCache<GLP, SwitchKey>(
                        *switch_keys[term]);
            z_maps[term] =
                baseXAutomorphismSpectrumMap<GLP>(x_multipliers[term]);
        }

        for (std::size_t prime = 0; prime < prime_count; prime++) {
            for (std::size_t term = 0; term < Count; term++) {
                for (std::uint32_t primary = 0;
                     primary < SwitchKey::primary_rows; primary++) {
                    const std::size_t source_row =
                        source_indices[term] * input_row_count + primary;
                    const std::size_t input_row =
                        term * input_row_count + primary;
                    const std::uint64_t *source_spectrum =
                        workspace.source_spectra[prime].data() +
                        source_row * coefficient_count;
                    std::uint64_t *input_spectrum =
                        workspace.input_spectra[prime].data() +
                        input_row * coefficient_count;
                    applyBaseXAutomorphismSpectrum<GLP>(
                        std::span<std::uint64_t>(input_spectrum,
                                                 coefficient_count),
                        std::span<const std::uint64_t>(source_spectrum,
                                                       coefficient_count),
                        z_maps[term]);
                }
            }
        }

        accumulateSmallKeySwitchPreparedSumNTT<GLP, SwitchKey, Count>(
            result, key_caches, workspace);
        return true;
    }
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
    using PackedCiphertext = GLPackedCiphertextData<GLP, BbarBit>;
    static constexpr std::uint64_t packed_ciphertext_payload_bytes =
        static_cast<std::uint64_t>(2) * GLP::matrix_dimension * GLP::baseP::n *
        sizeof(GLSignedDigit<BbarBit>);
    static constexpr std::uint64_t packed_payload_bytes =
        static_cast<std::uint64_t>(primary_rows) * bbar_rows * 2 *
        GLP::matrix_dimension * GLP::baseP::n * sizeof(GLSignedDigit<BbarBit>);

    static_assert(key_log_q <= std::numeric_limits<typename GLP::T>::digits,
                  "GL q*q0 key-switch modulus exceeds torus storage");
    static_assert(primary_rows * PrimaryBit <=
                  std::numeric_limits<typename GLP::T>::digits);
    static_assert(bbar_rows * BbarBit <=
                  std::numeric_limits<typename GLP::T>::digits);

    std::vector<PackedCiphertext> data{};

    void allocate() { data.resize(primary_rows * bbar_rows); }

    PackedCiphertext &at(const std::uint32_t primary, const std::uint32_t bbar)
    {
        return data.at(static_cast<std::size_t>(primary) * bbar_rows + bbar);
    }
    const PackedCiphertext &at(const std::uint32_t primary,
                               const std::uint32_t bbar) const
    {
        return data.at(static_cast<std::size_t>(primary) * bbar_rows + bbar);
    }

    std::uint64_t packedPayloadBytes() const
    {
        return static_cast<std::uint64_t>(data.size()) *
               packed_ciphertext_payload_bytes;
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(data);
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
    using PackedCiphertext = GLPackedBaseCiphertextData<GLP, BbarBit>;
    static constexpr std::uint64_t packed_ciphertext_payload_bytes =
        static_cast<std::uint64_t>(2) * GLP::baseP::n *
        sizeof(GLSignedDigit<BbarBit>);
    static constexpr std::uint64_t packed_payload_bytes =
        static_cast<std::uint64_t>(primary_rows) * bbar_rows * 2 *
        GLP::baseP::n * sizeof(GLSignedDigit<BbarBit>);

    static_assert(key_log_q <= std::numeric_limits<typename GLP::T>::digits,
                  "GL q*q0 key-switch modulus exceeds torus storage");
    static_assert(primary_rows * PrimaryBit <=
                  std::numeric_limits<typename GLP::T>::digits);
    static_assert(bbar_rows * BbarBit <=
                  std::numeric_limits<typename GLP::T>::digits);

    struct TransientNTTCache {
        std::once_flag initialize_once{};
        std::array<std::vector<std::uint64_t>, 2> spectra{};
        std::once_flag blocked_initialize_once{};
        std::array<std::vector<std::uint64_t>, 2> blocked_spectra{};
    };

    std::vector<PackedCiphertext> data{};
    // This exact transform cache is intentionally omitted from serialization.
    // Copies of an immutable evaluation key may share it safely.
    mutable std::shared_ptr<TransientNTTCache> transient_ntt_cache =
        std::make_shared<TransientNTTCache>();

    void allocate()
    {
        data.resize(primary_rows * bbar_rows);
        transient_ntt_cache = std::make_shared<TransientNTTCache>();
    }

    PackedCiphertext &at(const std::uint32_t primary, const std::uint32_t bbar)
    {
        return data.at(static_cast<std::size_t>(primary) * bbar_rows + bbar);
    }
    const PackedCiphertext &at(const std::uint32_t primary,
                               const std::uint32_t bbar) const
    {
        return data.at(static_cast<std::size_t>(primary) * bbar_rows + bbar);
    }

    std::uint64_t packedPayloadBytes() const
    {
        return static_cast<std::uint64_t>(data.size()) *
               packed_ciphertext_payload_bytes;
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(data);
    }
};

namespace gl_detail {

template <class GLP, class SwitchKey>
inline constexpr std::uint32_t bigKeySwitchAccumulationNTTPrimeCount = [] {
    if constexpr (!supportsWidePrimeNTT<GLP>) return 0U;
    constexpr std::size_t maximum_terms =
        static_cast<std::size_t>(GLP::matrix_dimension) * 2 *
        GLP::matrix_dimension * (2 * GLP::phi - 1);
    constexpr std::uint32_t required_bits =
        SwitchKey::primary_bit + SwitchKey::bbar_bit + ceilLog2(maximum_terms) +
        ceilLog2(SwitchKey::primary_rows);
    if constexpr (required_bits <= 60)
        return 1U;
    else if constexpr (required_bits <= 122)
        return 2U;
    else if constexpr (required_bits <= 126)
        return 3U;
    else
        return 0U;
}();

template <class GLP, class SwitchKey>
struct BigKeySwitchNTTCache {
    static constexpr std::uint32_t prime_count =
        bigKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    static constexpr std::size_t row_count =
        static_cast<std::size_t>(SwitchKey::primary_rows) *
        SwitchKey::bbar_rows * 2;
    std::array<std::vector<std::vector<std::uint64_t>>, 3> spectra{};
};

template <class GLP, class SwitchKey>
inline constexpr std::uint64_t bigKeySwitchNTTCacheBytes =
    static_cast<std::uint64_t>(
        bigKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>) *
    SwitchKey::primary_rows * SwitchKey::bbar_rows * 2 *
    GLPolynomialNTTPlan<GLP>::coefficient_count * sizeof(std::uint64_t);

template <class GLP, class SwitchKey>
inline void prepareBigKeySwitchNTTCache(
    BigKeySwitchNTTCache<GLP, SwitchKey> &cache, const SwitchKey &switch_key)
{
    constexpr std::uint32_t prime_count =
        bigKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    static_assert(prime_count != 0,
                  "the requested GL big DD switch has no exact NTT path");
    if (switch_key.data.size() !=
        SwitchKey::primary_rows * SwitchKey::bbar_rows)
        throw std::invalid_argument("uninitialized GL big key-switch key");

    using Plan = GLPolynomialNTTPlan<GLP>;
    std::array<const Plan *, 3> plans{&polynomialNTTPlan<GLP, 0>(), nullptr,
                                      nullptr};
    if constexpr (prime_count >= 2) plans[1] = &polynomialNTTPlan<GLP, 1>();
    if constexpr (prime_count == 3) plans[2] = &polynomialNTTPlan<GLP, 2>();
    constexpr std::size_t row_count =
        BigKeySwitchNTTCache<GLP, SwitchKey>::row_count;
    for (std::uint32_t prime = 0; prime < prime_count; prime++)
        cache.spectra[prime].resize(row_count);

#pragma omp parallel for collapse(2) schedule(static)
    for (std::uint32_t prime = 0; prime < prime_count; prime++) {
        for (std::size_t row = 0; row < row_count; row++) {
            const std::uint32_t primary =
                static_cast<std::uint32_t>(row / (SwitchKey::bbar_rows * 2));
            const std::size_t remainder =
                row % (static_cast<std::size_t>(SwitchKey::bbar_rows) * 2);
            const std::uint32_t bbar =
                static_cast<std::uint32_t>(remainder / 2);
            const std::size_t component = remainder % 2;
            plans[prime]->template forwardPacked<SwitchKey::bbar_bit>(
                cache.spectra[prime][row],
                switch_key.at(primary, bbar)[component]);
        }
    }
}

template <class GLP, class SwitchKey>
inline bool validBigKeySwitchNTTCache(
    const BigKeySwitchNTTCache<GLP, SwitchKey> &cache)
{
    constexpr std::uint32_t prime_count =
        bigKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    constexpr std::size_t row_count =
        BigKeySwitchNTTCache<GLP, SwitchKey>::row_count;
    constexpr std::size_t coefficient_count =
        GLPolynomialNTTPlan<GLP>::coefficient_count;
    for (std::uint32_t prime = 0; prime < prime_count; prime++) {
        if (cache.spectra[prime].size() != row_count) return false;
        for (const auto &row : cache.spectra[prime])
            if (row.size() != coefficient_count) return false;
    }
    return true;
}

}  // namespace gl_detail

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
                gl_detail::packDigitPolynomial<GLP, BbarBit>(
                    switch_key.at(primary, bbar)[component], rows[bbar]);
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
                gl_detail::packDigitPolynomial<GLP, BbarBit>(
                    switch_key.at(primary, bbar)[component], rows[bbar]);
        }
    }
}

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDBigKeySwitch(
    GLCiphertextData<GLP> &result, const GLPolynomial<GLP> &input,
    const GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key,
    const gl_detail::BigKeySwitchNTTCache<
        GLP, GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>> *ntt_cache =
        nullptr)
{
    using SwitchKey = GLDDBigKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;
    if (switch_key.data.size() !=
        SwitchKey::primary_rows * SwitchKey::bbar_rows)
        throw std::invalid_argument("uninitialized GL big key-switch key");

    constexpr std::uint32_t ntt_prime_count =
        gl_detail::bigKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;

    if constexpr (ntt_prime_count != 0) {
        using T = typename GLP::T;
        using Plan = gl_detail::GLPolynomialNTTPlan<GLP>;
        constexpr std::size_t coefficient_count = Plan::coefficient_count;
        constexpr std::size_t output_row_count =
            static_cast<std::size_t>(SwitchKey::bbar_rows) * 2;
        if (ntt_cache != nullptr &&
            !gl_detail::validBigKeySwitchNTTCache<GLP, SwitchKey>(*ntt_cache))
            throw std::invalid_argument("invalid GL big key-switch NTT cache");

        std::array<const Plan *, 3> plans{
            &gl_detail::polynomialNTTPlan<GLP, 0>(), nullptr, nullptr};
        if constexpr (ntt_prime_count >= 2)
            plans[1] = &gl_detail::polynomialNTTPlan<GLP, 1>();
        if constexpr (ntt_prime_count == 3)
            plans[2] = &gl_detail::polynomialNTTPlan<GLP, 2>();

        std::array<std::vector<std::uint64_t>, 3> input_spectra;
        plans[0]->template forwardActiveRows<LogQ, PrimaryBit>(input_spectra[0],
                                                               input);
        if constexpr (ntt_prime_count >= 2)
            plans[1]->template forwardActiveRows<LogQ, PrimaryBit>(
                input_spectra[1], input);
        if constexpr (ntt_prime_count == 3)
            plans[2]->template forwardActiveRows<LogQ, PrimaryBit>(
                input_spectra[2], input);

        // Every Bbar/component/prime inverse is independent.  Retaining its
        // 62-bit residue row lets the outer OpenMP team process all rows in
        // parallel and postpones the torus recombination to one cache-friendly
        // pass over the output.
        std::array<std::vector<std::vector<std::uint64_t>>, 3> residues;
        for (std::uint32_t prime = 0; prime < ntt_prime_count; prime++)
            residues[prime].resize(output_row_count);

#pragma omp parallel for collapse(2) schedule(static)
        for (std::uint32_t prime = 0; prime < ntt_prime_count; prime++) {
            for (std::size_t output_row = 0; output_row < output_row_count;
                 output_row++) {
                const std::uint32_t bbar =
                    static_cast<std::uint32_t>(output_row / 2);
                const std::size_t component = output_row % 2;
                const auto &plan = *plans[prime];
                const std::uint64_t modulus = plan.modulus();
                std::vector<std::uint64_t> accumulator(coefficient_count, 0);
                std::vector<std::uint64_t> transient_key_spectrum;
                for (std::uint32_t primary = 0;
                     primary < SwitchKey::primary_rows; primary++) {
                    const std::size_t key_row =
                        (static_cast<std::size_t>(primary) *
                             SwitchKey::bbar_rows +
                         bbar) *
                            2 +
                        component;
                    const std::uint64_t *key_spectrum;
                    if (ntt_cache == nullptr) {
                        plan.template forwardPacked<BbarBit>(
                            transient_key_spectrum,
                            switch_key.at(primary, bbar)[component]);
                        key_spectrum = transient_key_spectrum.data();
                    }
                    else {
                        key_spectrum =
                            ntt_cache->spectra[prime][key_row].data();
                    }
                    const std::uint64_t *input_row =
                        input_spectra[prime].data() +
                        static_cast<std::size_t>(primary) * coefficient_count;
                    for (std::size_t i = 0; i < coefficient_count; i++) {
                        const std::uint64_t product = modular_ntt::multiply(
                            input_row[i], key_spectrum[i], modulus);
                        accumulator[i] =
                            modular_ntt::add(accumulator[i], product, modulus);
                    }
                }
                plan.inverse(residues[prime][output_row], accumulator);
            }
        }

#pragma omp parallel for schedule(static)
        for (std::size_t index = 0; index < coefficient_count; index++) {
            const std::size_t y = index / GLP::baseP::n;
            const std::size_t coefficient = index % GLP::baseP::n;
            for (std::size_t component = 0; component < 2; component++) {
                T accumulated{};
                for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows;
                     bbar++) {
                    const std::uint32_t shift =
                        (SwitchKey::bbar_rows - bbar - 1) * BbarBit;
                    const std::size_t row =
                        static_cast<std::size_t>(bbar) * 2 + component;
                    if constexpr (ntt_prime_count == 1) {
                        accumulated += gl_detail::signedI128ToTorus<T>(
                                           modular_ntt::centeredResidue(
                                               residues[0][row][index],
                                               plans[0]->modulus()))
                                       << shift;
                    }
                    else if constexpr (ntt_prime_count == 2) {
                        static const modular_ntt::TwoPrimeCRT crt(
                            modular_ntt::wide_primes[0],
                            modular_ntt::wide_primes[1]);
                        accumulated +=
                            gl_detail::signedI128ToTorus<T>(
                                crt.reconstructSigned(residues[0][row][index],
                                                      residues[1][row][index]))
                            << shift;
                    }
                    else {
                        static const modular_ntt::ThreePrimeCRT crt(
                            modular_ntt::wide_primes[0],
                            modular_ntt::wide_primes[1],
                            modular_ntt::wide_primes[2]);
                        accumulated += gl_detail::signedI128ToTorus<T>(
                                           crt.reconstructSignedBounded(
                                               residues[0][row][index],
                                               residues[1][row][index],
                                               residues[2][row][index]))
                                       << shift;
                    }
                }
                result[component][y][coefficient] = accumulated;
            }
        }
        gl_detail::reduce<GLP, SwitchKey::key_log_q>(result[0]);
        gl_detail::reduce<GLP, SwitchKey::key_log_q>(result[1]);
        gl_detail::divideRoundLevel<GLP, SwitchKey::key_log_q,
                                    SwitchKey::auxiliary_log_q>(result[0]);
        gl_detail::divideRoundLevel<GLP, SwitchKey::key_log_q,
                                    SwitchKey::auxiliary_log_q>(result[1]);
        return;
    }

    auto input_digits =
        gl_detail::activeDecompose<GLP, LogQ, PrimaryBit>(input);
    std::array<std::vector<GLPolynomial<GLP>>, 2> digit_rows{
        std::vector<GLPolynomial<GLP>>(SwitchKey::bbar_rows),
        std::vector<GLPolynomial<GLP>>(SwitchKey::bbar_rows)};
    auto product = std::make_unique<GLPolynomial<GLP>>();
    auto key_row = std::make_unique<GLPolynomial<GLP>>();
    for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
         primary++) {
        for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++) {
            for (std::size_t component = 0; component < 2; component++) {
                gl_detail::unpackDigitPolynomial<GLP, BbarBit>(
                    *key_row, switch_key.at(primary, bbar)[component]);
                gl_detail::polynomialMultiply<GLP>(
                    *product, input_digits[primary], *key_row);
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

    if constexpr (gl_detail::smallKeySwitchAccumulationNTTPrimeCount<
                      GLP, SwitchKey> != 0) {
        using P = typename GLP::baseP;
        // Build the immutable key spectrum before entering the Y-slice team.
        // Otherwise the first worker performs this parallelizable work inside
        // an already-active OpenMP region while every other worker waits on
        // the once flag.
        const auto prepared_cache =
            gl_detail::prepareSmallKeySwitchNTTCache<GLP, SwitchKey>(
                switch_key);
        (void)prepared_cache;
        gl_detail::clear<GLP>(result[0]);
        gl_detail::clear<GLP>(result[1]);
#pragma omp parallel for schedule(dynamic)
        for (std::size_t slice = 0; slice < GLP::matrix_dimension; slice++) {
            auto input_digits = std::make_unique<
                std::array<Polynomial<P>, SwitchKey::primary_rows>>();
            const bool used_ntt =
                gl_detail::accumulateSmallKeySwitchProductsNTT<GLP, SwitchKey>(
                    1, switch_key,
                    [&](const std::size_t) {
                        ckks_detail::activeBaseDecomposePolynomialRows<
                            P, LogQ, PrimaryBit, SwitchKey::primary_rows>(
                            *input_digits, input[slice]);
                    },
                    [&](const std::uint32_t primary,
                        const std::size_t) -> const GLBasePolynomial<GLP> & {
                        return (*input_digits)[primary];
                    },
                    [&](const std::size_t component, const std::uint32_t bbar,
                        const std::size_t, const std::size_t coefficient,
                        const typename GLP::T value) {
                        const std::uint32_t shift =
                            (SwitchKey::bbar_rows - bbar - 1) * BbarBit;
                        result[component][slice][coefficient] += value << shift;
                    });
            if (!used_ntt)
                throw std::logic_error(
                    "GL small key switch lost its exact NTT path");
        }
        gl_detail::reduce<GLP, SwitchKey::key_log_q>(result[0]);
        gl_detail::reduce<GLP, SwitchKey::key_log_q>(result[1]);
    }
    else {
        auto input_digits =
            gl_detail::activeDecompose<GLP, LogQ, PrimaryBit>(input);
        std::array<std::vector<GLPolynomial<GLP>>, 2> digit_rows{
            std::vector<GLPolynomial<GLP>>(SwitchKey::bbar_rows),
            std::vector<GLPolynomial<GLP>>(SwitchKey::bbar_rows)};
        auto product = std::make_unique<GLBasePolynomial<GLP>>();
        auto key_row = std::make_unique<GLBasePolynomial<GLP>>();
        for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
             primary++) {
            for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++) {
                for (std::size_t component = 0; component < 2; component++) {
                    gl_detail::unpackDigitPolynomial<GLP, BbarBit>(
                        *key_row, switch_key.at(primary, bbar)[component]);
                    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
                        gl_detail::baseMultiply<GLP>(
                            *product, input_digits[primary][y], *key_row);
                        gl_detail::addInPlace<GLP>(
                            digit_rows[component][bbar][y], *product);
                    }
                }
            }
        }
        for (std::size_t component = 0; component < 2; component++)
            gl_detail::activeRecombine<GLP, SwitchKey::key_log_q, BbarBit>(
                result[component], digit_rows[component]);
    }
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

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(conjugate_key, product_key);
    }
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

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(amount, multiplier, switch_key);
    }
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

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(amount, multiplier, switch_key);
    }
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
        &transpose_key,
    const gl_detail::BigKeySwitchNTTCache<
        GLP, GLConjugateTransposeKey<GLP, LogQ, PrimaryBit, BbarBit>>
        *ntt_cache = nullptr)
{
    constexpr std::uint32_t inverse_x = 4 * GLP::matrix_dimension - 1;
    constexpr std::uint32_t inverse_w = GLP::cyclotomic_order - 1;
    GLPolynomial<GLP> body;
    GLPolynomial<GLP> mask;
    gl_detail::polynomialAutomorphism<GLP>(body, input[0], 3, inverse_x,
                                           inverse_x, inverse_w, true);
    gl_detail::polynomialAutomorphism<GLP>(mask, input[1], 3, inverse_x,
                                           inverse_x, inverse_w, true);
    GLDDBigKeySwitch(result.ct, mask, transpose_key, ntt_cache);
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
