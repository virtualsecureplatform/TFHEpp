#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <type_traits>
#include <tuple>
#include <utility>
#include <vector>

#include "bfv++.hpp"
#include "bfv-slots.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
struct CKKSCiphertext {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);
    static_assert(LogDelta < LogQ);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t log_budget = LogQ - LogDelta;

    TRLWE<P> ct{};
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
struct CKKSPlaintext {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);
    static_assert(LogDelta < LogQ);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;

    Polynomial<P> poly{};
};

template <class P, std::uint32_t LogQ>
struct CKKSRelinKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    std::array<TRLWE<P>, P::l̅> data{};

    TRLWE<P> &operator[](std::size_t i) { return data[i]; }
    const TRLWE<P> &operator[](std::size_t i) const { return data[i]; }
};

template <class P, std::uint32_t LogQ>
struct CKKSAutoKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    std::array<std::array<TRLWE<P>, P::l̅>, P::k> data{};

    std::array<TRLWE<P>, P::l̅> &operator[](std::size_t i) { return data[i]; }
    const std::array<TRLWE<P>, P::l̅> &operator[](std::size_t i) const
    {
        return data[i];
    }
};

template <class P, std::uint32_t LogQ>
using CKKSGaloisKey = std::array<CKKSAutoKey<P, LogQ>, P::nbit + 1>;

struct CKKSNoise {
    double modular_stdev = 0.0;
    std::uint32_t uniform_bits = 0;
};

namespace ckks_detail {

template <std::uint32_t A, std::uint32_t B>
inline constexpr std::uint32_t static_min_v = A < B ? A : B;

template <std::uint32_t A, std::uint32_t B>
inline constexpr std::uint32_t static_max_v = A < B ? B : A;

template <class P>
inline constexpr std::uint32_t torus_width_v =
    std::numeric_limits<typename P::T>::digits;

template <class P>
inline constexpr bool supported_torus_v =
    std::is_same_v<typename P::T, __uint128_t> ||
    is_multilimb_uint_v<typename P::T>;

template <class P>
inline __int128_t unsignedToI128(typename P::T value)
{
    using T = typename P::T;
    if constexpr (std::is_same_v<T, __uint128_t>) {
        return static_cast<__int128_t>(value);
    }
    else {
        static_assert(is_multilimb_uint_v<T>);
        __uint128_t low = value.limb[0];
        if constexpr (T::limbs > 1)
            low |= static_cast<__uint128_t>(value.limb[1]) << 64;
        return static_cast<__int128_t>(low);
    }
}

template <class P>
inline long double unsignedToLongDouble(typename P::T value)
{
    using T = typename P::T;
    if constexpr (std::is_same_v<T, __uint128_t>) {
        constexpr long double limb_scale = 18446744073709551616.0L;
        const auto lo = static_cast<std::uint64_t>(value);
        const auto hi = static_cast<std::uint64_t>(value >> 64);
        return static_cast<long double>(hi) * limb_scale +
               static_cast<long double>(lo);
    }
    else {
        static_assert(is_multilimb_uint_v<T>);
        long double out = 0.0L;
        for (std::size_t i = T::limbs; i-- > 0;)
            out = std::ldexp(out, 64) +
                  static_cast<long double>(value.limb[i]);
        return out;
    }
}

template <class P>
inline __int128_t torusToSigned(typename P::T value)
{
    static_assert(supported_torus_v<P>);
    if constexpr (std::is_same_v<typename P::T, __uint128_t>) {
        return static_cast<__int128_t>(value);
    }
    else {
        static_assert(torus_width_v<P> <= 128,
                      "Use levelToSigned only for <=128-bit multi-limb levels");
        return unsignedToI128<P>(value);
    }
}

template <class P>
inline typename P::T signedToTorus(__int128_t value)
{
    static_assert(supported_torus_v<P>);
    return static_cast<typename P::T>(value);
}

template <class P, std::uint32_t LogQ>
inline constexpr typename P::T levelMask()
{
    static_assert(supported_torus_v<P>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= torus_width_v<P>);
    if constexpr (LogQ == torus_width_v<P>)
        return std::numeric_limits<typename P::T>::max();
    else
        return (typename P::T{1} << LogQ) - typename P::T{1};
}

template <class P, std::uint32_t LogQ>
inline typename P::T reduceToLevel(typename P::T value)
{
    return value & levelMask<P, LogQ>();
}

template <class P, std::uint32_t LogQ>
inline typename P::T signedToLevel(__int128_t value)
{
    return reduceToLevel<P, LogQ>(signedToTorus<P>(value));
}

template <class P, std::uint32_t LogQ>
inline __int128_t levelToSigned(typename P::T value)
{
    value = reduceToLevel<P, LogQ>(value);
    if constexpr (std::is_same_v<typename P::T, __uint128_t> &&
                  LogQ == torus_width_v<P>) {
        return torusToSigned<P>(value);
    }
    else {
        static_assert(LogQ < 128,
                      "levelToSigned returns __int128_t and cannot represent "
                      "larger active CKKS levels");
        const typename P::T sign = typename P::T{1} << (LogQ - 1);
        if ((value & sign) == typename P::T{0}) return unsignedToI128<P>(value);

        const typename P::T modulus = typename P::T{1} << LogQ;
        const typename P::T magnitude = modulus - value;
        return -unsignedToI128<P>(magnitude);
    }
}

template <class P, std::uint32_t LogQ>
inline long double levelToLongDouble(typename P::T value)
{
    using T = typename P::T;
    static_assert(LogQ > 0);
    static_assert(LogQ <= torus_width_v<P>);

    value = reduceToLevel<P, LogQ>(value);
    const T sign = T{1} << (LogQ - 1);
    if ((value & sign) == T{0}) return unsignedToLongDouble<P>(value);

    T magnitude;
    if constexpr (LogQ == torus_width_v<P>)
        magnitude = T{0} - value;
    else
        magnitude = (T{1} << LogQ) - value;
    return -unsignedToLongDouble<P>(magnitude);
}

template <class P, std::uint32_t LogQ>
inline typename P::T centeredLevelToTorus(typename P::T value)
{
    value = reduceToLevel<P, LogQ>(value);
    if constexpr (LogQ == torus_width_v<P>) {
        return value;
    }
    else {
        const typename P::T sign = typename P::T{1} << (LogQ - 1);
        if ((value & sign) == typename P::T{0}) return value;
        return value | ~levelMask<P, LogQ>();
    }
}

template <class P, std::uint32_t LogQ>
inline typename P::T uniformAtLevel()
{
    if constexpr (LogQ == torus_width_v<P>) {
        return UniformTorusRandom<P>();
    }
    else {
        std::uniform_int_distribution<std::uint64_t> dist64(
            0, std::numeric_limits<std::uint64_t>::max());
        typename P::T value = 0;
        std::uint32_t generated = 0;
        while (generated < LogQ) {
            value |= static_cast<typename P::T>(dist64(generator)) << generated;
            generated += 64;
        }
        return reduceToLevel<P, LogQ>(value);
    }
}

inline __int128_t randomSignedBits(std::uint32_t bits)
{
    assert(bits > 0);
    assert(bits < 127);

    std::uniform_int_distribution<std::uint64_t> dist64(
        0, std::numeric_limits<std::uint64_t>::max());
    __uint128_t value = 0;
    std::uint32_t generated = 0;
    while (generated < bits) {
        value |= static_cast<__uint128_t>(dist64(generator)) << generated;
        generated += 64;
    }
    value &= (static_cast<__uint128_t>(1) << bits) - 1;

    const __uint128_t sign = static_cast<__uint128_t>(1) << (bits - 1);
    if ((value & sign) == 0) return static_cast<__int128_t>(value);
    return static_cast<__int128_t>(value) -
           static_cast<__int128_t>(static_cast<__uint128_t>(1) << bits);
}

inline void fftInPlace(std::vector<std::complex<long double>> &a, bool inverse)
{
    const std::size_t n = a.size();
    assert(n != 0 && (n & (n - 1)) == 0);

    for (std::size_t i = 1, j = 0; i < n; i++) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const long double angle =
            (inverse ? 2.0L : -2.0L) * pi / static_cast<long double>(len);
        const std::complex<long double> wlen(std::cos(angle),
                                             std::sin(angle));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<long double> w = 1.0L;
            for (std::size_t j = 0; j < len / 2; j++) {
                const std::complex<long double> u = a[i + j];
                const std::complex<long double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (inverse)
        for (auto &x : a) x /= static_cast<long double>(n);
}

template <class P, std::uint32_t LogQ>
inline typename P::T sampleNoise(const CKKSNoise &noise)
{
    typename P::T value = 0;
    if (noise.modular_stdev != 0.0) {
        std::normal_distribution<long double> distribution(0.0L, 1.0L);
        const long double scaled =
            std::ldexp(static_cast<long double>(noise.modular_stdev) *
                           distribution(generator),
                       LogQ);
        value += signedToLevel<P, LogQ>(longDoubleToI128(scaled));
    }
    if (noise.uniform_bits > 0) {
        assert(noise.uniform_bits < LogQ);
        value += signedToLevel<P, LogQ>(randomSignedBits(noise.uniform_bits));
    }
    return reduceToLevel<P, LogQ>(value);
}

template <class P, std::uint32_t LogQ>
inline void reducePolynomialToLevel(Polynomial<P> &poly)
{
    for (std::uint32_t i = 0; i < P::n; i++)
        poly[i] = reduceToLevel<P, LogQ>(poly[i]);
}

template <class P, std::uint32_t LogQ>
inline void reduceTRLWEToLevel(TRLWE<P> &ct)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        reducePolynomialToLevel<P, LogQ>(ct[c]);
}

template <class P, std::uint32_t LogQ>
inline void reduceTRLWE3ToLevel(TRLWE3<P> &ct)
{
    for (int c = 0; c < 3; c++) reducePolynomialToLevel<P, LogQ>(ct[c]);
}

template <class P>
inline bool fitsU64(typename P::T value)
{
    if constexpr (std::is_same_v<typename P::T, __uint128_t>) {
        return value <= static_cast<__uint128_t>(
                            std::numeric_limits<std::uint64_t>::max());
    }
    else {
        using T = typename P::T;
        static_assert(is_multilimb_uint_v<T>);
        for (std::size_t i = 1; i < T::limbs; i++)
            if (value.limb[i] != 0) return false;
        return true;
    }
}

template <class P>
inline std::uint64_t lowU64(typename P::T value)
{
    if constexpr (std::is_same_v<typename P::T, __uint128_t>)
        return static_cast<std::uint64_t>(value);
    else {
        static_assert(is_multilimb_uint_v<typename P::T>);
        return value.limb[0];
    }
}

template <class P, std::uint32_t LogQ>
inline std::pair<bool, std::uint64_t> smallSignedMagnitude(typename P::T value)
{
    using T = typename P::T;
    static_assert(LogQ > 0);
    static_assert(LogQ <= torus_width_v<P>);

    value = reduceToLevel<P, LogQ>(value);
    const T sign = T{1} << (LogQ - 1);
    const bool negative = (value & sign) != T{0};
    T magnitude = value;
    if (negative) {
        if constexpr (LogQ == torus_width_v<P>)
            magnitude = T{0} - value;
        else
            magnitude = (T{1} << LogQ) - value;
    }

    assert(fitsU64<P>(magnitude));
    return {negative, lowU64<P>(magnitude)};
}

template <class P, std::uint32_t LogQ, std::uint32_t DropBits>
inline typename P::T roundedLevelRightShift(typename P::T value)
{
    using T = typename P::T;
    static_assert(DropBits > 0);
    static_assert(DropBits < LogQ);
    static_assert(LogQ <= torus_width_v<P>);
    constexpr std::uint32_t out_log_q = LogQ - DropBits;

    value = reduceToLevel<P, LogQ>(value);
    const T sign = T{1} << (LogQ - 1);
    const bool negative = (value & sign) != T{0};
    const T half = T{1} << (DropBits - 1);

    if (!negative)
        return reduceToLevel<P, out_log_q>((value + half) >> DropBits);

    T magnitude;
    if constexpr (LogQ == torus_width_v<P>)
        magnitude = T{0} - value;
    else
        magnitude = (T{1} << LogQ) - value;
    const T rounded = (magnitude + half) >> DropBits;
    return reduceToLevel<P, out_log_q>(T{0} - rounded);
}

template <class P, std::uint32_t PlainLogQ>
inline void polyMulTorusBySmallSigned(Polynomial<P> &res,
                                      const Polynomial<P> &torus,
                                      const Polynomial<P> &plain)
{
    using T = typename P::T;
    static_assert(PlainLogQ > 0);
    static_assert(PlainLogQ <= torus_width_v<P>);

    if constexpr (is_multilimb_uint_v<T>) {
        auto pos = std::make_unique<std::array<std::uint64_t, P::n>>();
        auto neg = std::make_unique<std::array<std::uint64_t, P::n>>();
        bool has_pos = false;
        bool has_neg = false;
        for (std::uint32_t i = 0; i < P::n; i++) {
            const auto [negative, magnitude] =
                smallSignedMagnitude<P, PlainLogQ>(plain[i]);
            (*pos)[i] = negative ? 0 : magnitude;
            (*neg)[i] = negative ? magnitude : 0;
            has_pos = has_pos || (!negative && magnitude != 0);
            has_neg = has_neg || (negative && magnitude != 0);
        }

        if (has_pos)
            PolyMulTorusByUnsignedBits<P, 64>(res, torus, *pos);
        else
            for (std::uint32_t i = 0; i < P::n; i++) res[i] = 0;

        if (has_neg) {
            auto tmp = std::make_unique<Polynomial<P>>();
            PolyMulTorusByUnsignedBits<P, 64>(*tmp, torus, *neg);
            for (std::uint32_t i = 0; i < P::n; i++) res[i] -= (*tmp)[i];
        }
    }
    else {
        static_assert(std::is_same_v<T, __uint128_t>);
        for (int i = 0; i < static_cast<int>(P::n); i++) {
            T ri = 0;
            for (int j = 0; j <= i; j++)
                ri += torus[j] * plain[i - j];
            for (int j = i + 1; j < static_cast<int>(P::n); j++)
                ri -= torus[j] * plain[P::n + i - j];
            res[i] = reduceToLevel<P, PlainLogQ>(ri);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void centeredTRLWEAtLevel(TRLWE<P> &out, const TRLWE<P> &in)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++)
            out[c][n] = centeredLevelToTorus<P, LogQ>(in[c][n]);
}

template <class P, std::uint32_t LogScale>
inline typename P::T rescaleAccumulator(Wide384 acc)
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    static_assert(LogScale < torus_width_v<P>);
    if constexpr (LogScale > 0) {
        const bool negative = (acc.w[5] >> 63) != 0;
        acc.add_shifted(negative ? -static_cast<__int128_t>(1)
                                 : static_cast<__int128_t>(1),
                        static_cast<int>(LogScale) - 1);
    }
    const typename P::T divisor = LogScale == 0
                                      ? static_cast<typename P::T>(1)
                                      : (static_cast<typename P::T>(1)
                                         << LogScale);
    return acc.div128(divisor);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void encodeRealPolynomial(Polynomial<P> &poly,
                                 const std::array<double, P::n> &coeffs)
{
    static_assert(LogDelta < LogQ);
    for (std::uint32_t i = 0; i < P::n; i++) {
        const long double scaled =
            std::ldexp(static_cast<long double>(coeffs[i]), LogDelta);
        const auto rounded = static_cast<__int128_t>(std::round(scaled));
        poly[i] = signedToLevel<P, LogQ>(rounded);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void decodeRealPolynomial(std::array<double, P::n> &coeffs,
                                 const Polynomial<P> &poly)
{
    static_assert(LogDelta < LogQ);
    for (std::uint32_t i = 0; i < P::n; i++) {
        const long double value =
            std::ldexp(levelToLongDouble<P, LogQ>(poly[i]),
                       -static_cast<int>(LogDelta));
        coeffs[i] = static_cast<double>(value);
    }
}

template <class P>
inline const std::array<std::uint32_t, P::n / 2> &ckksSlotToEvalIndex()
{
    static const auto table = [] {
        std::array<std::uint32_t, P::n / 2> tbl{};
        std::uint64_t g = 1;
        constexpr std::uint64_t twon = 2 * P::n;
        for (std::uint32_t i = 0; i < P::n / 2; i++) {
            tbl[i] = static_cast<std::uint32_t>((g - 1) / 2);
            g = (g * 5) % twon;
        }
        return tbl;
    }();
    return table;
}

inline std::uint32_t bitReverse(std::uint32_t value, std::uint32_t bits)
{
    std::uint32_t reversed = 0;
    for (std::uint32_t i = 0; i < bits; i++) {
        reversed = (reversed << 1) | (value & 1U);
        value >>= 1;
    }
    return reversed;
}

template <std::uint32_t LogDelta, class P>
inline void fillConjugateSymmetricEvaluations(
    std::array<std::complex<double>, P::n> &evals,
    const std::array<std::complex<double>, P::n / 2> &slots)
{
    const double scale = std::ldexp(1.0, LogDelta);
    const auto &slot_to_eval = ckksSlotToEvalIndex<P>();
    for (std::uint32_t i = 0; i < P::n / 2; i++) {
        const std::uint32_t eval_index = slot_to_eval[i];
        evals[eval_index] = slots[i] * scale;
        evals[P::n - 1 - eval_index] = std::conj(evals[eval_index]);
    }
}

template <class P, std::uint32_t LogQ>
inline void encryptPolynomialAtLevel(TRLWE<P> &ct, const Polynomial<P> &poly,
                                     const Key<P> &key,
                                     CKKSNoise noise = {P::α, 0})
{
    static_assert(supported_torus_v<P>);

    for (std::uint32_t i = 0; i < P::n; i++)
        ct[P::k][i] =
            reduceToLevel<P, LogQ>(poly[i] + sampleNoise<P, LogQ>(noise));

    for (int k = 0; k < static_cast<int>(P::k); k++) {
        for (std::uint32_t i = 0; i < P::n; i++)
            ct[k][i] = uniformAtLevel<P, LogQ>();

        Polynomial<P> partkey{};
        for (std::uint32_t i = 0; i < P::n; i++)
            partkey[i] = key[k * P::n + i];

        Polynomial<P> mask_phase{};
        PolyMul<P>(mask_phase, ct[k], partkey);
        for (std::uint32_t i = 0; i < P::n; i++)
            ct[P::k][i] =
                reduceToLevel<P, LogQ>(ct[P::k][i] + mask_phase[i]);
    }
}

template <class P>
inline void polyMulTorusByBbarDigit(Polynomial<P> &res,
                                    const Polynomial<P> &torus,
                                    const Polynomial<P> &digit)
{
    if constexpr (is_multilimb_uint_v<typename P::T>) {
        PolyMulTorusByDigit<P>(res, torus, digit);
    }
    else {
        static_assert(std::is_same_v<typename P::T, __uint128_t>);
        static_assert(P::l̅ * P::B̅gbit ==
                          std::numeric_limits<typename P::T>::digits,
                      "CKKS Bbar digits must cover the torus");

        using T = typename P::T;
        constexpr int width = std::numeric_limits<T>::digits;
        constexpr T half = static_cast<T>(1) << (P::B̅gbit - 1);
        constexpr T mask = (static_cast<T>(1) << P::B̅gbit) - 1;
        constexpr T offset = [] {
            constexpr int local_width = std::numeric_limits<T>::digits;
            constexpr T local_half = static_cast<T>(1) << (P::B̅gbit - 1);
            T value = 0;
            for (int j = 0; j < static_cast<int>(P::l̅); j++)
                value += local_half << (local_width - (j + 1) * P::B̅gbit);
            return value;
        }();

        for (std::uint32_t n = 0; n < P::n; n++) res[n] = 0;

        auto torus_digit = std::make_unique<Polynomial<P>>();
        auto product = std::make_unique<Polynomial<P>>();
        auto digit_fft = std::make_unique<PolynomialInFD<P>>();
        auto torus_digit_fft = std::make_unique<PolynomialInFD<P>>();
        auto product_fft = std::make_unique<PolynomialInFD<P>>();
        TwistIFFT<P>(*digit_fft, digit);

        for (int j = 0; j < static_cast<int>(P::l̅); j++) {
            const int shift = width - (j + 1) * static_cast<int>(P::B̅gbit);
            for (std::uint32_t n = 0; n < P::n; n++)
                (*torus_digit)[n] =
                    (((torus[n] + offset) >> shift) & mask) - half;

            TwistIFFT<P>(*torus_digit_fft, *torus_digit);
            MulInFD<P::n>(*product_fft, *torus_digit_fft, *digit_fft);
            TwistFFT<P>(*product, *product_fft);

            for (std::uint32_t n = 0; n < P::n; n++)
                res[n] += (*product)[n] << shift;
        }
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline typename P::T ckksEncodeCoeff(double value)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta < LogQ);
    const long double scaled =
        std::ldexp(static_cast<long double>(value), LogDelta);
    return ckks_detail::signedToLevel<P, LogQ>(
        static_cast<__int128_t>(std::round(scaled)));
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline double ckksDecodeCoeff(typename P::T value)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta < LogQ);
    const long double decoded =
        std::ldexp(ckks_detail::levelToLongDouble<P, LogQ>(value),
                   -static_cast<int>(LogDelta));
    return static_cast<double>(decoded);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksEncodePolynomial(Polynomial<P> &poly,
                                 const std::array<double, P::n> &coeffs)
{
    ckks_detail::encodeRealPolynomial<P, LogQ, LogDelta>(poly, coeffs);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksDecodePolynomial(std::array<double, P::n> &coeffs,
                                 const Polynomial<P> &poly)
{
    ckks_detail::decodeRealPolynomial<P, LogQ, LogDelta>(coeffs, poly);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotEncode(
    Polynomial<P> &poly,
    const std::array<std::complex<double>, P::n / 2> &slots)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta < LogQ);
    std::array<std::complex<double>, P::n> evals{};
    ckks_detail::fillConjugateSymmetricEvaluations<LogDelta, P>(evals, slots);

    constexpr double pi = 3.141592653589793238462643383279502884;
    std::vector<std::complex<long double>> spectrum(P::n);
    for (std::uint32_t k = 0; k < P::n; k++)
        spectrum[k] = static_cast<std::complex<long double>>(evals[k]);
    ckks_detail::fftInPlace(spectrum, false);

    for (std::uint32_t j = 0; j < P::n; j++) {
        const long double angle =
            -pi * static_cast<long double>(j) / static_cast<long double>(P::n);
        const auto coeff =
            spectrum[j] *
            std::complex<long double>(std::cos(angle), std::sin(angle)) /
            static_cast<long double>(P::n);
        poly[j] = ckks_detail::signedToLevel<P, LogQ>(
            static_cast<__int128_t>(std::round(coeff.real())));
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotDecode(
    std::array<std::complex<double>, P::n / 2> &slots,
    const Polynomial<P> &poly)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta < LogQ);
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double inv_scale = std::ldexp(1.0, -static_cast<int>(LogDelta));
    const auto &slot_to_eval = ckks_detail::ckksSlotToEvalIndex<P>();
    std::vector<std::complex<long double>> values(P::n);
    for (std::uint32_t j = 0; j < P::n; j++) {
        const long double angle =
            pi * static_cast<long double>(j) / static_cast<long double>(P::n);
        const long double coeff =
            ckks_detail::levelToLongDouble<P, LogQ>(poly[j]);
        values[j] =
            coeff * std::complex<long double>(std::cos(angle), std::sin(angle));
    }
    ckks_detail::fftInPlace(values, true);

    for (std::uint32_t k = 0; k < P::n / 2; k++) {
        const std::uint32_t eval_index = slot_to_eval[k];
        const auto value = values[eval_index] * static_cast<long double>(P::n);
        slots[k] = static_cast<std::complex<double>>(value) * inv_scale;
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksEncrypt(CKKSCiphertext<P, LogQ, LogDelta> &out,
                        const std::array<double, P::n> &coeffs,
                        const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    Polynomial<P> poly;
    ckksEncodePolynomial<P, LogQ, LogDelta>(poly, coeffs);
    ckks_detail::encryptPolynomialAtLevel<P, LogQ>(out.ct, poly, key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksEncryptPolynomial(CKKSCiphertext<P, LogQ, LogDelta> &out,
                                  const Polynomial<P> &poly,
                                  const Key<P> &key,
                                  CKKSNoise noise = {P::α, 0})
{
    ckks_detail::encryptPolynomialAtLevel<P, LogQ>(out.ct, poly, key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksDecrypt(std::array<double, P::n> &out,
                        const CKKSCiphertext<P, LogQ, LogDelta> &ct,
                        const Key<P> &key)
{
    Polynomial<P> phase = trlwePhase<P>(ct.ct, key);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(phase);
    ckksDecodePolynomial<P, LogQ, LogDelta>(out, phase);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotEncrypt(
    CKKSCiphertext<P, LogQ, LogDelta> &out,
    const std::array<std::complex<double>, P::n / 2> &slots,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    Polynomial<P> poly;
    ckksSlotEncode<P, LogQ, LogDelta>(poly, slots);
    ckksEncryptPolynomial<P, LogQ, LogDelta>(out, poly, key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotDecrypt(
    std::array<std::complex<double>, P::n / 2> &slots,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct, const Key<P> &key)
{
    Polynomial<P> phase = trlwePhase<P>(ct.ct, key);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(phase);
    ckksSlotDecode<P, LogQ, LogDelta>(slots, phase);
}

template <class P>
using CKKSSlotVector = std::array<std::complex<double>, P::n / 2>;

template <class P>
inline void rotateCKKSSlotVector(CKKSSlotVector<P> &out,
                                 const CKKSSlotVector<P> &in, int steps)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    for (int i = 0; i < half; i++) out[i] = in[(i + steps) % half];
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
struct CKKSPlainMulTraits {
    static_assert(PlainLogDelta < LogQ);
    static constexpr std::uint32_t log_q = LogQ - PlainLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
using CKKSPlainMulResult =
    typename CKKSPlainMulTraits<P, LogQ, LogDelta,
                                PlainLogDelta>::Ciphertext;

template <class P, std::uint32_t InLogQ, std::uint32_t OutLogQ,
          std::uint32_t LogDelta>
inline void CKKSModRaise(CKKSCiphertext<P, OutLogQ, LogDelta> &res,
                         const CKKSCiphertext<P, InLogQ, LogDelta> &ct)
{
    static_assert(InLogQ < OutLogQ);
    static_assert(OutLogQ <= ckks_detail::torus_width_v<P>);
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            res.ct[c][i] = ckks_detail::reduceToLevel<P, InLogQ>(ct.ct[c][i]);
}

template <class P, std::uint32_t LogQ, std::uint32_t PlainLogDelta>
inline void CKKSPlainMulRescaleTRLWE(TRLWE<P> &res, const TRLWE<P> &ct,
                                     const Polynomial<P> &plain)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    auto product = std::make_unique<Polynomial<P>>();
    for (int c = 0; c <= static_cast<int>(P::k); c++) {
        ckks_detail::polyMulTorusBySmallSigned<P, LogQ>(*product, ct[c], plain);
        for (std::uint32_t i = 0; i < P::n; i++)
            res[c][i] =
                ckks_detail::roundedLevelRightShift<P, LogQ, PlainLogDelta>(
                    (*product)[i]);
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSPlainMulRescale(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSPlaintext<P, LogQ, PlainLogDelta> &plain)
{
    CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(res.ct, ct.ct,
                                                     plain.poly);
}

template <class P, std::uint32_t LogQ>
inline void CKKSAddTRLWEInPlace(TRLWE<P> &acc, const TRLWE<P> &term)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            acc[c][i] = ckks_detail::reduceToLevel<P, LogQ>(acc[c][i] +
                                                            term[c][i]);
}

template <class P, std::uint32_t LogQ>
inline void CKKSSubTRLWEInPlace(TRLWE<P> &acc, const TRLWE<P> &term)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            acc[c][i] = ckks_detail::reduceToLevel<P, LogQ>(acc[c][i] -
                                                            term[c][i]);
}

template <class P, std::uint32_t InLogQ, std::uint32_t OutLogQ,
          std::uint32_t LogDelta>
inline void CKKSLevelReduce(CKKSCiphertext<P, OutLogQ, LogDelta> &res,
                            const CKKSCiphertext<P, InLogQ, LogDelta> &ct)
{
    static_assert(OutLogQ <= InLogQ);
    static_assert(LogDelta < OutLogQ);
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            res.ct[c][i] = ckks_detail::reduceToLevel<P, OutLogQ>(ct.ct[c][i]);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSAddInPlace(CKKSCiphertext<P, LogQ, LogDelta> &acc,
                           const CKKSCiphertext<P, LogQ, LogDelta> &term)
{
    CKKSAddTRLWEInPlace<P, LogQ>(acc.ct, term.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSAdd(CKKSCiphertext<P, LogQ, LogDelta> &res,
                    const CKKSCiphertext<P, LogQ, LogDelta> &lhs,
                    const CKKSCiphertext<P, LogQ, LogDelta> &rhs)
{
    res = lhs;
    CKKSAddInPlace<P, LogQ, LogDelta>(res, rhs);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSubInPlace(CKKSCiphertext<P, LogQ, LogDelta> &acc,
                           const CKKSCiphertext<P, LogQ, LogDelta> &term)
{
    CKKSSubTRLWEInPlace<P, LogQ>(acc.ct, term.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSub(CKKSCiphertext<P, LogQ, LogDelta> &res,
                    const CKKSCiphertext<P, LogQ, LogDelta> &lhs,
                    const CKKSCiphertext<P, LogQ, LogDelta> &rhs)
{
    res = lhs;
    CKKSSubInPlace<P, LogQ, LogDelta>(res, rhs);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSPlainMulByReal(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct, double scalar)
{
    CKKSSlotVector<P> slots{};
    slots.fill({scalar, 0.0});
    CKKSPlaintext<P, LogQ, PlainLogDelta> plain;
    ckksSlotEncode<P, LogQ, PlainLogDelta>(plain.poly, slots);
    CKKSPlainMulRescale<P, LogQ, LogDelta, PlainLogDelta>(res, ct, plain);
}

template <class P, std::uint32_t LogQ, class Rows>
inline void CKKSKeySwitchRows(TRLWE<P> &res, const Polynomial<P> &poly,
                              const Rows &rows)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(P::l̅ * P::B̅gbit ==
                      std::numeric_limits<typename P::T>::digits,
                  "CKKS key switching requires full Bbar decomposition");

    TRLWE<P> poly_as_trlwe{};
    poly_as_trlwe[0] = poly;
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(poly_as_trlwe);

    auto decomposed = std::make_unique<std::array<TRLWE<P>, P::l̅>>();
    TRLWEBaseBbarDecompose<P>(*decomposed, poly_as_trlwe);

    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++) res[c][n] = 0;

    auto product = std::make_unique<Polynomial<P>>();
    for (int j = 0; j < static_cast<int>(P::l̅); j++) {
        const Polynomial<P> &digit = (*decomposed)[j][0];
        for (int c = 0; c <= static_cast<int>(P::k); c++) {
            ckks_detail::polyMulTorusByBbarDigit<P>(*product, rows[j][c],
                                                    digit);
            for (std::uint32_t n = 0; n < P::n; n++)
                res[c][n] =
                    ckks_detail::reduceToLevel<P, LogQ>(res[c][n] +
                                                        (*product)[n]);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSAutoKeyGen(CKKSAutoKey<P, LogQ> &autokey, const uint d,
                           const Key<P> &key,
                           CKKSNoise noise = {P::α, 0})
{
    static_assert(ckks_detail::supported_torus_v<P>);
    constexpr int width = std::numeric_limits<typename P::T>::digits;

    for (int k = 0; k < static_cast<int>(P::k); k++) {
        Polynomial<P> partkey{};
        for (std::uint32_t i = 0; i < P::n; i++)
            partkey[i] = key[k * P::n + i];

        Polynomial<P> autokey_poly{};
        Automorphism<P>(autokey_poly, partkey, d);

        for (int j = 0; j < static_cast<int>(P::l̅); j++) {
            const int shift = width - (j + 1) * static_cast<int>(P::B̅gbit);
            Polynomial<P> gadget{};
            for (std::uint32_t n = 0; n < P::n; n++)
                gadget[n] = ckks_detail::reduceToLevel<P, LogQ>(
                    autokey_poly[n] << shift);
            ckks_detail::encryptPolynomialAtLevel<P, LogQ>(
                autokey[k][j], gadget, key, noise);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSGaloisKeyGen(CKKSGaloisKey<P, LogQ> &gk, const Key<P> &key,
                             CKKSNoise noise = {P::α, 0})
{
    std::uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit); i++) {
        CKKSAutoKeyGen<P, LogQ>(gk[i], static_cast<uint>(d), key, noise);
        d = d * d % (2 * P::n);
    }
    CKKSAutoKeyGen<P, LogQ>(gk[P::nbit], 2 * P::n - 1, key, noise);
}

template <class P, std::uint32_t LogQ>
inline void CKKSEvalAuto(TRLWE<P> &res, const TRLWE<P> &trlwe, const int d,
                         const CKKSAutoKey<P, LogQ> &autokey)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++) res[c][n] = 0;
    Automorphism<P>(res[P::k], trlwe[P::k], d);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(res[P::k]);

    auto temppoly = std::make_unique<Polynomial<P>>();
    auto switched = std::make_unique<TRLWE<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        Automorphism<P>(*temppoly, trlwe[k], d);
        CKKSKeySwitchRows<P, LogQ>(*switched, *temppoly, autokey[k]);
        CKKSSubTRLWEInPlace<P, LogQ>(res, *switched);
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSGaloisKey<P, LogQ> &gk)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    if (steps == 0) {
        res = ct;
        ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
        return;
    }

    std::uint64_t d = 5;
    auto cur = std::make_unique<TRLWE<P>>(ct);
    auto tmp = std::make_unique<TRLWE<P>>();
    for (int i = 0; i < static_cast<int>(P::nbit) - 1; i++) {
        if ((steps >> i) & 1) {
            CKKSEvalAuto<P, LogQ>(*tmp, *cur, static_cast<int>(d), gk[i]);
            *cur = *tmp;
        }
        d = d * d % (2 * P::n);
    }
    res = *cur;
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ>
inline void CKKSConjugateSlots(TRLWE<P> &res, const TRLWE<P> &ct,
                               const CKKSGaloisKey<P, LogQ> &gk)
{
    CKKSEvalAuto<P, LogQ>(res, ct, 2 * P::n - 1, gk[P::nbit]);
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ, std::uint32_t PlainLogDelta>
struct CKKSLinearTransformTerm {
    int baby_step = 0;
    CKKSPlaintext<P, LogQ, PlainLogDelta> plain{};
};

template <class P, std::uint32_t LogQ, std::uint32_t PlainLogDelta>
struct CKKSLinearTransformGiantStep {
    int giant_step = 0;
    std::vector<CKKSLinearTransformTerm<P, LogQ, PlainLogDelta>> terms{};
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
struct CKKSLinearTransformPlan {
    static_assert(PlainLogDelta < LogQ);
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t plain_log_delta = PlainLogDelta;
    static constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;

    int k_step = 0;
    std::vector<CKKSLinearTransformGiantStep<P, LogQ, PlainLogDelta>> groups{};
};

template <class P>
struct CKKSLinearTransformStage {
    std::vector<CKKSSlotVector<P>> diagonals{};
    std::vector<int> rotation_offsets{};
};

template <class P>
using CKKSLinearTransformStages = std::vector<CKKSLinearTransformStage<P>>;

template <class P>
inline void CKKSScaleLinearTransformStages(CKKSLinearTransformStages<P> &stages,
                                           std::complex<double> scalar)
{
    for (auto &stage : stages) {
        if (stage.diagonals.empty()) continue;
        for (auto &diagonal : stage.diagonals)
            for (auto &entry : diagonal) entry *= scalar;
        return;
    }
}

namespace ckks_detail {

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          class IndexSequence>
struct CKKSGaloisKeyChainTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t... Is>
struct CKKSGaloisKeyChainTuple<P, StartLogQ, PlainLogDelta,
                               std::index_sequence<Is...>> {
    using type = std::tuple<CKKSGaloisKey<P, StartLogQ - Is * PlainLogDelta>...>;
};

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
struct CKKSGaloisKeyChain {
    static_assert(StartLogQ > StageCount * PlainLogDelta);

    using Tuple = typename ckks_detail::CKKSGaloisKeyChainTuple<
        P, StartLogQ, PlainLogDelta,
        std::make_index_sequence<StageCount + 1>>::type;

    Tuple keys{};

    template <std::size_t I>
    auto &get()
    {
        return std::get<I>(keys);
    }

    template <std::size_t I>
    const auto &get() const
    {
        return std::get<I>(keys);
    }
};

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
inline void CKKSGaloisKeyChainGenImpl(
    CKKSGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount> &chain,
    const Key<P> &key, CKKSNoise noise)
{
    if constexpr (I <= StageCount) {
        constexpr std::uint32_t log_q = StartLogQ - I * PlainLogDelta;
        CKKSGaloisKeyGen<P, log_q>(chain.template get<I>(), key, noise);
        CKKSGaloisKeyChainGenImpl<I + 1, P, StartLogQ, PlainLogDelta,
                                  StageCount>(chain, key, noise);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
inline void CKKSGaloisKeyChainGen(
    CKKSGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount> &chain,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    ckks_detail::CKKSGaloisKeyChainGenImpl<0, P, StartLogQ, PlainLogDelta,
                                           StageCount>(chain, key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildLinearTransformBSGSPlan(
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const std::vector<CKKSSlotVector<P>> &diagonals,
    const std::vector<int> &rotation_offsets, int k_step)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(diagonals.size() == rotation_offsets.size());
    assert(!diagonals.empty());
    assert(k_step > 0 && k_step <= half);

    int max_j2 = 0;
    for (const int offset : rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        max_j2 = std::max(max_j2, r / k_step);
    }

    std::vector<std::vector<CKKSSlotVector<P>>> accumulated(
        max_j2 + 1, std::vector<CKKSSlotVector<P>>(k_step));
    std::vector<std::vector<bool>> used(max_j2 + 1,
                                        std::vector<bool>(k_step, false));
    for (auto &group : accumulated)
        for (auto &diag : group) diag.fill({0.0, 0.0});

    auto rotated = std::make_unique<CKKSSlotVector<P>>();
    for (std::size_t i = 0; i < rotation_offsets.size(); i++) {
        const int r = ((rotation_offsets[i] % half) + half) % half;
        const int j2 = r / k_step;
        const int j1 = r % k_step;
        rotateCKKSSlotVector<P>(*rotated, diagonals[i], -k_step * j2);
        for (int s = 0; s < half; s++) accumulated[j2][j1][s] += (*rotated)[s];
        used[j2][j1] = true;
    }

    plan.k_step = k_step;
    plan.groups.clear();
    plan.groups.reserve(max_j2 + 1);
    for (int j2 = 0; j2 <= max_j2; j2++) {
        CKKSLinearTransformGiantStep<P, LogQ, PlainLogDelta> group;
        group.giant_step = j2;
        for (int j1 = 0; j1 < k_step; j1++) {
            if (!used[j2][j1]) continue;
            CKKSLinearTransformTerm<P, LogQ, PlainLogDelta> term;
            term.baby_step = j1;
            ckksSlotEncode<P, LogQ, PlainLogDelta>(term.plain.poly,
                                                   accumulated[j2][j1]);
            group.terms.push_back(std::move(term));
        }
        if (!group.terms.empty()) plan.groups.push_back(std::move(group));
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildLinearTransformBSGSPlan(
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, stage.diagonals, stage.rotation_offsets, k_step);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSLinearTransformBSGS(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const CKKSGaloisKey<P, LogQ> &input_gk,
    const CKKSGaloisKey<P, LogQ - PlainLogDelta> &output_gk)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    assert(plan.k_step > 0 && plan.k_step <= static_cast<int>(P::n / 2));
    assert(!plan.groups.empty());

    auto baby = std::make_unique<std::vector<TRLWE<P>>>(plan.k_step);
    (*baby)[0] = ct.ct;
    for (int j1 = 1; j1 < plan.k_step; j1++)
        CKKSRotateSlots<P, LogQ>((*baby)[j1], ct.ct, j1, input_gk);

    bool res_initialized = false;
    auto term = std::make_unique<TRLWE<P>>();
    auto inner = std::make_unique<TRLWE<P>>();
    auto giant_rotated = std::make_unique<TRLWE<P>>();

    for (const auto &group : plan.groups) {
        bool inner_initialized = false;
        for (const auto &entry : group.terms) {
            CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(
                *term, (*baby)[entry.baby_step], entry.plain.poly);

            if (!inner_initialized) {
                *inner = *term;
                inner_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(*inner, *term);
            }
        }

        const TRLWE<P> *to_add = inner.get();
        if (group.giant_step != 0) {
            CKKSRotateSlots<P, out_log_q>(
                *giant_rotated, *inner, plan.k_step * group.giant_step,
                output_gk);
            to_add = giant_rotated.get();
        }

        if (!res_initialized) {
            res.ct = *to_add;
            res_initialized = true;
        }
        else {
            CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, *to_add);
        }
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSLinearTransformBSGS(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const std::vector<CKKSSlotVector<P>> &diagonals,
    const std::vector<int> &rotation_offsets, int k_step,
    const CKKSGaloisKey<P, LogQ> &input_gk,
    const CKKSGaloisKey<P, LogQ - PlainLogDelta> &output_gk)
{
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> plan;
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, diagonals, rotation_offsets, k_step);
    CKKSLinearTransformBSGS<P, LogQ, LogDelta, PlainLogDelta>(
        res, ct, plan, input_gk, output_gk);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
struct CKKSStagedPlainMulTraits {
    static_assert(LogQ > StageCount * PlainLogDelta);
    static constexpr std::uint32_t log_q = LogQ - StageCount * PlainLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
using CKKSStagedPlainMulResult =
    typename CKKSStagedPlainMulTraits<P, LogQ, LogDelta, PlainLogDelta,
                                      StageCount>::Ciphertext;

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
inline void CKKSLinearTransformStagesBSGSImpl(
    CKKSStagedPlainMulResult<P, StartLogQ, LogDelta, PlainLogDelta,
                             StageCount> &res,
    const CKKSCiphertext<P, StartLogQ - I * PlainLogDelta, LogDelta> &ct,
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    int k_step,
    const CKKSGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount>
        &gk_chain)
{
    if constexpr (I == StageCount) {
        res = ct;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * PlainLogDelta;
        CKKSLinearTransformPlan<P, log_q, LogDelta, PlainLogDelta> plan;
        CKKSBuildLinearTransformBSGSPlan<P, log_q, LogDelta, PlainLogDelta>(
            plan, stages[first_stage + I], k_step);

        CKKSPlainMulResult<P, log_q, LogDelta, PlainLogDelta> next;
        CKKSLinearTransformBSGS<P, log_q, LogDelta, PlainLogDelta>(
            next, ct, plan, gk_chain.template get<I>(),
            gk_chain.template get<I + 1>());
        CKKSLinearTransformStagesBSGSImpl<I + 1, P, StartLogQ, LogDelta,
                                          PlainLogDelta, StageCount>(
            res, next, stages, first_stage, k_step, gk_chain);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
inline void CKKSLinearTransformStagesBSGS(
    CKKSStagedPlainMulResult<P, LogQ, LogDelta, PlainLogDelta, StageCount>
        &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    int k_step,
    const CKKSGaloisKeyChain<P, LogQ, PlainLogDelta, StageCount> &gk_chain)
{
    assert(first_stage + StageCount <= stages.size());
    ckks_detail::CKKSLinearTransformStagesBSGSImpl<0, P, LogQ, LogDelta,
                                                   PlainLogDelta, StageCount>(
        res, ct, stages, first_stage, k_step, gk_chain);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
struct CKKSRealLinearTransformPlan {
    static_assert(PlainLogDelta < LogQ);
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t plain_log_delta = PlainLogDelta;
    static constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;

    bool has_direct = false;
    bool has_conjugate = false;
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> direct{};
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> conjugate{};
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildRealLinearTransformPlan(
    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const std::vector<CKKSSlotVector<P>> &direct_diagonals,
    const std::vector<int> &direct_offsets,
    const std::vector<CKKSSlotVector<P>> &conjugate_diagonals,
    const std::vector<int> &conjugate_offsets, int k_step)
{
    assert(!direct_diagonals.empty() || !conjugate_diagonals.empty());
    plan.has_direct = !direct_diagonals.empty();
    plan.has_conjugate = !conjugate_diagonals.empty();
    if (plan.has_direct)
        CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
            plan.direct, direct_diagonals, direct_offsets, k_step);
    if (plan.has_conjugate)
        CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
            plan.conjugate, conjugate_diagonals, conjugate_offsets, k_step);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSRealLinearTransform(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const CKKSGaloisKey<P, LogQ> &input_gk,
    const CKKSGaloisKey<P, LogQ - PlainLogDelta> &output_gk)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    assert(plan.has_direct || plan.has_conjugate);

    bool initialized = false;
    auto term = std::make_unique<CKKSPlainMulResult<P, LogQ, LogDelta,
                                                   PlainLogDelta>>();
    if (plan.has_direct) {
        CKKSLinearTransformBSGS<P, LogQ, LogDelta, PlainLogDelta>(
            *term, ct, plan.direct, input_gk, output_gk);
        res = *term;
        initialized = true;
    }

    if (plan.has_conjugate) {
        auto conjugated = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
        CKKSConjugateSlots<P, LogQ>(conjugated->ct, ct.ct, input_gk);
        CKKSLinearTransformBSGS<P, LogQ, LogDelta, PlainLogDelta>(
            *term, *conjugated, plan.conjugate, input_gk, output_gk);
        if (!initialized) {
            res = *term;
            initialized = true;
        }
        else {
            CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, term->ct);
        }
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
}

namespace ckks_detail {

template <class P>
inline void addCKKSStageDiagonal(CKKSLinearTransformStage<P> &stage, int offset,
                                 const CKKSSlotVector<P> &diag)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    offset = ((offset % half) + half) % half;
    for (std::size_t i = 0; i < stage.rotation_offsets.size(); i++) {
        if (stage.rotation_offsets[i] != offset) continue;
        for (int j = 0; j < half; j++) stage.diagonals[i][j] += diag[j];
        return;
    }
    stage.rotation_offsets.push_back(offset);
    stage.diagonals.push_back(diag);
}

template <class P>
inline std::vector<std::complex<long double>> ckksPowerOfFiveTable()
{
    constexpr int half = static_cast<int>(P::n) / 2;
    std::vector<std::complex<long double>> roots(4 * half + 1);
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    for (int i = 0; i <= 4 * half; i++) {
        const long double angle =
            2.0L * pi * static_cast<long double>(i) /
            static_cast<long double>(4 * half);
        roots[i] = {std::cos(angle), std::sin(angle)};
    }
    return roots;
}

template <class P>
inline std::vector<int> ckksPowerOfFiveExponents()
{
    constexpr int half = static_cast<int>(P::n) / 2;
    std::vector<int> pow5(2 * half + 1);
    pow5[0] = 1;
    for (int i = 1; i <= 2 * half; i++)
        pow5[i] = (pow5[i - 1] * 5) & (4 * half - 1);
    return pow5;
}

}  // namespace ckks_detail

template <class P>
inline std::uint32_t CKKSBitReverseSlotIndex(std::uint32_t index)
{
    static_assert(P::nbit > 0);
    assert(index < P::n / 2);
    return ckks_detail::bitReverse(index, P::nbit - 1);
}

// Maps CKKS slots to packed coefficients in bit-reversed order:
// out[bitreverse(i)] = coeff[i] + i*coeff[i + N/2].
template <class P>
inline void CKKSBuildCoeffToPackedSlotStages(
    CKKSLinearTransformStages<P> &stages)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    constexpr int log_slots = static_cast<int>(P::nbit) - 1;
    static_assert(log_slots > 0);

    stages.clear();
    stages.reserve(log_slots);
    const auto roots = ckks_detail::ckksPowerOfFiveTable<P>();
    const auto pow5 = ckks_detail::ckksPowerOfFiveExponents<P>();
    const long double stage_scale = 0.5L;

    for (int m = half; m >= 2; m >>= 1) {
        CKKSLinearTransformStage<P> stage;
        CKKSSlotVector<P> a{};
        CKKSSlotVector<P> b{};
        CKKSSlotVector<P> c{};
        a.fill({0.0, 0.0});
        b.fill({0.0, 0.0});
        c.fill({0.0, 0.0});

        const int tt = m >> 1;
        const int gap = half / m;
        const int mask = (m << 2) - 1;
        for (int i = 0; i < half; i += m) {
            for (int j = 0; j < tt; j++) {
                const int k = ((m << 2) - (pow5[j] & mask)) * gap;
                const int idx1 = i + j;
                const int idx2 = idx1 + tt;
                a[idx1] = {static_cast<double>(stage_scale), 0.0};
                a[idx2] = static_cast<std::complex<double>>(
                    -stage_scale * roots[k]);
                b[idx1] = {static_cast<double>(stage_scale), 0.0};
                c[idx2] = static_cast<std::complex<double>>(
                    stage_scale * roots[k]);
            }
        }

        ckks_detail::addCKKSStageDiagonal<P>(stage, 0, a);
        ckks_detail::addCKKSStageDiagonal<P>(stage, tt, b);
        ckks_detail::addCKKSStageDiagonal<P>(stage, half - tt, c);
        stages.push_back(std::move(stage));
    }
}

// Inverse of CKKSBuildCoeffToPackedSlotStages.  The input is the same
// bit-reversed packed coefficient layout emitted by CoeffToPackedSlot.
template <class P>
inline void CKKSBuildPackedSlotToCoeffStages(
    CKKSLinearTransformStages<P> &stages)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    constexpr int log_slots = static_cast<int>(P::nbit) - 1;
    static_assert(log_slots > 0);

    stages.clear();
    stages.reserve(log_slots);
    const auto roots = ckks_detail::ckksPowerOfFiveTable<P>();
    const auto pow5 = ckks_detail::ckksPowerOfFiveExponents<P>();

    for (int m = 2; m <= half; m <<= 1) {
        CKKSLinearTransformStage<P> stage;
        CKKSSlotVector<P> a{};
        CKKSSlotVector<P> b{};
        CKKSSlotVector<P> c{};
        a.fill({0.0, 0.0});
        b.fill({0.0, 0.0});
        c.fill({0.0, 0.0});

        const int tt = m >> 1;
        const int gap = half / m;
        const int mask = (m << 2) - 1;
        for (int i = 0; i < half; i += m) {
            for (int j = 0; j < tt; j++) {
                const int k = (pow5[j] & mask) * gap;
                const int idx1 = i + j;
                const int idx2 = idx1 + tt;
                a[idx1] = {1.0, 0.0};
                a[idx2] =
                    static_cast<std::complex<double>>(-roots[k]);
                b[idx1] = static_cast<std::complex<double>>(roots[k]);
                c[idx2] = {1.0, 0.0};
            }
        }

        ckks_detail::addCKKSStageDiagonal<P>(stage, 0, a);
        ckks_detail::addCKKSStageDiagonal<P>(stage, tt, b);
        ckks_detail::addCKKSStageDiagonal<P>(stage, half - tt, c);
        stages.push_back(std::move(stage));
    }
}

template <class P>
inline void CKKSBuildCoeffToPackedSlotDiagonals(
    std::vector<CKKSSlotVector<P>> &direct_diagonals,
    std::vector<int> &direct_offsets,
    std::vector<CKKSSlotVector<P>> &conjugate_diagonals,
    std::vector<int> &conjugate_offsets)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const auto &slot_to_eval = ckks_detail::ckksSlotToEvalIndex<P>();

    direct_diagonals.assign(half, CKKSSlotVector<P>{});
    conjugate_diagonals.assign(half, CKKSSlotVector<P>{});
    direct_offsets.resize(half);
    conjugate_offsets.resize(half);
    for (int d = 0; d < half; d++) {
        direct_offsets[d] = d;
        conjugate_offsets[d] = d;
        direct_diagonals[d].fill({0.0, 0.0});
        conjugate_diagonals[d].fill({0.0, 0.0});
    }

    for (int r = 0; r < half; r++) {
        for (int d = 0; d < half; d++) {
            const int c = (r + d) % half;
            const long double h =
                static_cast<long double>(2 * slot_to_eval[c] + 1);
            const long double theta = pi * h / static_cast<long double>(P::n);
            const auto e0 = std::complex<long double>(
                std::cos(theta * r), -std::sin(theta * r));
            const auto e1 = std::complex<long double>(
                std::cos(theta * (r + half)),
                -std::sin(theta * (r + half)));
            const auto f0 = std::conj(e0);
            const auto f1 = std::conj(e1);
            direct_diagonals[d][r] = static_cast<std::complex<double>>(
                (e0 + std::complex<long double>(0.0L, 1.0L) * e1) /
                static_cast<long double>(P::n));
            conjugate_diagonals[d][r] = static_cast<std::complex<double>>(
                (f0 + std::complex<long double>(0.0L, 1.0L) * f1) /
                static_cast<long double>(P::n));
        }
    }
}

template <class P>
inline void CKKSBuildPackedSlotToCoeffDiagonals(
    std::vector<CKKSSlotVector<P>> &direct_diagonals,
    std::vector<int> &direct_offsets,
    std::vector<CKKSSlotVector<P>> &conjugate_diagonals,
    std::vector<int> &conjugate_offsets)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const auto &slot_to_eval = ckks_detail::ckksSlotToEvalIndex<P>();

    direct_diagonals.assign(half, CKKSSlotVector<P>{});
    conjugate_diagonals.assign(half, CKKSSlotVector<P>{});
    direct_offsets.resize(half);
    conjugate_offsets.resize(half);
    for (int d = 0; d < half; d++) {
        direct_offsets[d] = d;
        conjugate_offsets[d] = d;
        direct_diagonals[d].fill({0.0, 0.0});
        conjugate_diagonals[d].fill({0.0, 0.0});
    }

    for (int i = 0; i < half; i++) {
        const long double h =
            static_cast<long double>(2 * slot_to_eval[i] + 1);
        const long double theta = pi * h / static_cast<long double>(P::n);
        for (int d = 0; d < half; d++) {
            const int c = (i + d) % half;
            const auto e0 = std::complex<long double>(
                std::cos(theta * c), std::sin(theta * c));
            const auto e1 = std::complex<long double>(
                std::cos(theta * (c + half)),
                std::sin(theta * (c + half)));
            direct_diagonals[d][i] = static_cast<std::complex<double>>(
                0.5L *
                (e0 - std::complex<long double>(0.0L, 1.0L) * e1));
            conjugate_diagonals[d][i] = static_cast<std::complex<double>>(
                0.5L *
                (e0 + std::complex<long double>(0.0L, 1.0L) * e1));
        }
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildCoeffToPackedSlotPlan(
    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    int k_step)
{
    std::vector<CKKSSlotVector<P>> direct_diagonals;
    std::vector<CKKSSlotVector<P>> conjugate_diagonals;
    std::vector<int> direct_offsets;
    std::vector<int> conjugate_offsets;
    CKKSBuildCoeffToPackedSlotDiagonals<P>(
        direct_diagonals, direct_offsets, conjugate_diagonals,
        conjugate_offsets);
    CKKSBuildRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, direct_diagonals, direct_offsets, conjugate_diagonals,
        conjugate_offsets, k_step);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildPackedSlotToCoeffPlan(
    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    int k_step)
{
    std::vector<CKKSSlotVector<P>> direct_diagonals;
    std::vector<CKKSSlotVector<P>> conjugate_diagonals;
    std::vector<int> direct_offsets;
    std::vector<int> conjugate_offsets;
    CKKSBuildPackedSlotToCoeffDiagonals<P>(
        direct_diagonals, direct_offsets, conjugate_diagonals,
        conjugate_offsets);
    CKKSBuildRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, direct_diagonals, direct_offsets, conjugate_diagonals,
        conjugate_offsets, k_step);
}

template <class P, std::uint32_t LogQ>
inline std::unique_ptr<CKKSRelinKey<P, LogQ>> makeCKKSRelinKey(
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    static_assert(ckks_detail::supported_torus_v<P>);

    auto relinkey = std::make_unique<CKKSRelinKey<P, LogQ>>();
    auto keysquare = std::make_unique<Polynomial<P>>();
    Polynomial<P> partkey{};
    for (std::uint32_t i = 0; i < P::n; i++) partkey[i] = key[i];
    PolyMul<P>(*keysquare, partkey, partkey);

    constexpr int width = std::numeric_limits<typename P::T>::digits;
    for (int j = 0; j < static_cast<int>(P::l̅); j++) {
        const int shift = width - (j + 1) * static_cast<int>(P::B̅gbit);
        Polynomial<P> gadget{};
        for (std::uint32_t n = 0; n < P::n; n++)
            gadget[n] = ckks_detail::reduceToLevel<P, LogQ>((*keysquare)[n]
                                                            << shift);
        ckks_detail::encryptPolynomialAtLevel<P, LogQ>((*relinkey)[j], gadget,
                                                       key, noise);
    }

    return relinkey;
}

template <class P, std::uint32_t LogQ>
inline void CKKSRelinKeySwitch(TRLWE<P> &res, const Polynomial<P> &poly,
                               const CKKSRelinKey<P, LogQ> &relinkey)
{
    CKKSKeySwitchRows<P, LogQ>(res, poly, relinkey.data);
}

template <class P, std::uint32_t LogQ>
inline void CKKSRelinearization(TRLWE<P> &res, const TRLWE3<P> &mult,
                                const CKKSRelinKey<P, LogQ> &relinkey)
{
    TRLWE<P> squareterm;
    CKKSRelinKeySwitch<P, LogQ>(squareterm, mult[2], relinkey);
    for (std::uint32_t i = 0; i < P::n; i++)
        res[0][i] = ckks_detail::reduceToLevel<P, LogQ>(mult[0][i] +
                                                        squareterm[0][i]);
    for (std::uint32_t i = 0; i < P::n; i++)
        res[1][i] = ckks_detail::reduceToLevel<P, LogQ>(mult[1][i] +
                                                        squareterm[1][i]);
}

template <class P, std::uint32_t LhsLogQ, std::uint32_t RhsLogQ,
          std::uint32_t LogScale>
inline void CKKSTensorProductRescale(TRLWE3<P> &res, const TRLWE<P> &a,
                                     const TRLWE<P> &b)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    constexpr int width = std::numeric_limits<typename P::T>::digits;
    static_assert(P::l̅ * P::B̅gbit == width,
                  "CKKS FullDD multiply requires Bbar digits to cover T");
    static_assert(LogScale < static_cast<std::uint32_t>(width));
    constexpr std::uint32_t base_log_q =
        ckks_detail::static_min_v<LhsLogQ, RhsLogQ>;
    static_assert(base_log_q > LogScale);
    constexpr std::uint32_t out_log_q = base_log_q - LogScale;

    auto a_centered = std::make_unique<TRLWE<P>>();
    auto b_centered = std::make_unique<TRLWE<P>>();
    ckks_detail::centeredTRLWEAtLevel<P, LhsLogQ>(*a_centered, a);
    ckks_detail::centeredTRLWEAtLevel<P, RhsLogQ>(*b_centered, b);

    auto a_dec = std::make_unique<std::array<TRLWE<P>, P::l̅>>();
    auto b_dec = std::make_unique<std::array<TRLWE<P>, P::l̅>>();
    TRLWEBaseBbarDecompose<P>(*a_dec, *a_centered);
    TRLWEBaseBbarDecompose<P>(*b_dec, *b_centered);

    if constexpr (std::is_same_v<typename P::T, __uint128_t>) {
        constexpr bool use_fft_digits =
            2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3 <
            std::numeric_limits<double>::digits;

        auto acc = std::make_unique<std::array<std::array<Wide384, P::n>, 3>>();

        if constexpr (use_fft_digits) {
            auto a_fd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
            auto b_fd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    TwistIFFT<P>((*a_fd)[i][c], (*a_dec)[i][c]);
                    TwistIFFT<P>((*b_fd)[i][c], (*b_dec)[i][c]);
                }
            }

            alignas(64) PolynomialInFD<P> prod_fd;
            Polynomial<P> prod;
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                    const int digit_sum = i + j;
                    const int shift =
                        2 * width -
                        (digit_sum + 2) * static_cast<int>(P::B̅gbit);
                    if (shift <= -width) continue;

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][0], (*b_fd)[j][0]);
                    TwistFFT<P>(prod, prod_fd);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[2][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][1], (*b_fd)[j][1]);
                    TwistFFT<P>(prod, prod_fd);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[1][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][0], (*b_fd)[j][1]);
                    FMAInFD<P::n>(prod_fd, (*a_fd)[i][1], (*b_fd)[j][0]);
                    TwistFFT<P>(prod, prod_fd);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                }
            }
        }
        else {
            Polynomial<P> prod;
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                    const int digit_sum = i + j;
                    const int shift =
                        2 * width -
                        (digit_sum + 2) * static_cast<int>(P::B̅gbit);
                    if (shift <= -width) continue;

                    PolyMulNaive<P>(prod, (*a_dec)[i][0], (*b_dec)[j][0]);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[2][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][1], (*b_dec)[j][1]);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[1][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][0], (*b_dec)[j][1]);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][1], (*b_dec)[j][0]);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                }
            }
        }

        for (int c = 0; c < 3; c++)
            for (std::uint32_t n = 0; n < P::n; n++)
                res[c][n] = ckks_detail::reduceToLevel<P, out_log_q>(
                    ckks_detail::rescaleAccumulator<P, LogScale>((*acc)[c][n]));
    }
    else {
        using Torus = typename P::T;
        using Acc = WideSignedLimbAccumulator<2 * Torus::limbs + 2>;
        auto acc = std::make_unique<std::array<std::array<Acc, P::n>, 3>>();

        auto accumulate = [&](int component, const Polynomial<P> &prod,
                              int shift) {
            for (std::uint32_t n = 0; n < P::n; n++)
                (*acc)[component][n].add_shifted_i64(
                    multilimb_to_signed_i64(prod[n]), shift);
        };

        constexpr bool use_fft_digits =
            use_multilimb_digit_fft_v<P> &&
            2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3 <
                std::numeric_limits<double>::digits;

        if constexpr (use_fft_digits) {
            constexpr int digit_product_bits =
                2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3;
            constexpr int fd_batch_slack =
                std::numeric_limits<double>::digits - digit_product_bits;
            static_assert(fd_batch_slack > 0);
            constexpr int fd_batch_size =
                std::min(static_cast<int>(P::l̅), 1 << fd_batch_slack);

            auto a_fd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
            auto b_fd = std::make_unique<std::array<TRLWEInFD<P>, P::l̅>>();
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    TwistIFFTDigit<P>((*a_fd)[i][c], (*a_dec)[i][c]);
                    TwistIFFTDigit<P>((*b_fd)[i][c], (*b_dec)[i][c]);
                }
            }

            auto sum_fd = std::make_unique<std::array<PolynomialInFD<P>, 3>>();
            auto prod = std::make_unique<Polynomial<P>>();
            for (int digit_sum = 0;
                 digit_sum <= 2 * static_cast<int>(P::l̅) - 2; digit_sum++) {
                const int i_begin =
                    std::max(0, digit_sum - static_cast<int>(P::l̅) + 1);
                const int i_end =
                    std::min(static_cast<int>(P::l̅) - 1, digit_sum);
                const int shift =
                    2 * width -
                    (digit_sum + 2) * static_cast<int>(P::B̅gbit);

                for (int batch_begin = i_begin; batch_begin <= i_end;
                     batch_begin += fd_batch_size) {
                    for (int c = 0; c < 3; c++) (*sum_fd)[c].fill(0.0);
                    const int batch_end =
                        std::min(i_end, batch_begin + fd_batch_size - 1);
                    for (int i = batch_begin; i <= batch_end; i++) {
                        const int j = digit_sum - i;
                        FMAInFD<P::n>((*sum_fd)[2], (*a_fd)[i][0],
                                      (*b_fd)[j][0]);
                        FMAInFD<P::n>((*sum_fd)[1], (*a_fd)[i][1],
                                      (*b_fd)[j][1]);
                        FMAInFD<P::n>((*sum_fd)[0], (*a_fd)[i][0],
                                      (*b_fd)[j][1]);
                        FMAInFD<P::n>((*sum_fd)[0], (*a_fd)[i][1],
                                      (*b_fd)[j][0]);
                    }
                    for (int c = 0; c < 3; c++) {
                        TwistFFTDigitProduct<P>(*prod, (*sum_fd)[c]);
                        accumulate(c, *prod, shift);
                    }
                }
            }
        }
        else {
            auto prod = std::make_unique<Polynomial<P>>();
            for (int i = 0; i < static_cast<int>(P::l̅); i++) {
                for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                    const int digit_sum = i + j;
                    const int shift =
                        2 * width -
                        (digit_sum + 2) * static_cast<int>(P::B̅gbit);

                    PolyMulDigit<P>(*prod, (*a_dec)[i][0], (*b_dec)[j][0]);
                    accumulate(2, *prod, shift);
                    PolyMulDigit<P>(*prod, (*a_dec)[i][1], (*b_dec)[j][1]);
                    accumulate(1, *prod, shift);
                    PolyMulDigit<P>(*prod, (*a_dec)[i][0], (*b_dec)[j][1]);
                    accumulate(0, *prod, shift);
                    PolyMulDigit<P>(*prod, (*a_dec)[i][1], (*b_dec)[j][0]);
                    accumulate(0, *prod, shift);
                }
            }
        }

        const Torus divisor =
            LogScale == 0 ? Torus{1} : (Torus{1} << LogScale);
        for (int c = 0; c < 3; c++) {
            for (std::uint32_t n = 0; n < P::n; n++) {
                if constexpr (LogScale > 0) {
                    (*acc)[c][n].add_shifted_i64(
                        (*acc)[c][n].is_negative() ? -1 : 1,
                        static_cast<int>(LogScale) - 1);
                }
                res[c][n] = ckks_detail::reduceToLevel<P, out_log_q>(
                    (*acc)[c][n].template div_to_torus<Torus::limbs>(divisor));
            }
        }
    }
}

template <class P, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
struct CKKSMultTraits {
    static constexpr std::uint32_t log_scale =
        ckks_detail::static_max_v<LhsLogDelta, RhsLogDelta>;
    static constexpr std::uint32_t base_log_q =
        ckks_detail::static_min_v<LhsLogQ, RhsLogQ>;
    static_assert(base_log_q > log_scale);

    static constexpr std::uint32_t log_q = base_log_q - log_scale;
    static constexpr std::uint32_t log_delta =
        ckks_detail::static_min_v<LhsLogDelta, RhsLogDelta>;

    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
using CKKSMultResult =
    typename CKKSMultTraits<P, LhsLogQ, LhsLogDelta, RhsLogQ,
                            RhsLogDelta>::Ciphertext;

template <class P, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
inline void CKKSMult(
    CKKSMultResult<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta> &res,
    const CKKSCiphertext<P, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<P, RhsLogQ, RhsLogDelta> &rhs,
    const CKKSRelinKey<P, CKKSMultTraits<P, LhsLogQ, LhsLogDelta, RhsLogQ,
                                        RhsLogDelta>::log_q> &relinkey)
{
    using Traits =
        CKKSMultTraits<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta>;

    TRLWE3<P> mult;
    CKKSTensorProductRescale<P, LhsLogQ, RhsLogQ, Traits::log_scale>(
        mult, lhs.ct, rhs.ct);
    CKKSRelinearization<P, Traits::log_q>(res.ct, mult, relinkey);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSquare(
    CKKSMultResult<P, LogQ, LogDelta, LogQ, LogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSRelinKey<P,
                       CKKSMultTraits<P, LogQ, LogDelta, LogQ,
                                      LogDelta>::log_q> &relinkey)
{
    CKKSMult<P>(res, ct, ct, relinkey);
}

inline double CKKSPlainEvalModSineDegree5(double x)
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    const double x2 = x * x;
    return x * (1.0 - (two_pi * two_pi / 6.0) * x2 +
                (two_pi * two_pi * two_pi * two_pi / 120.0) * x2 * x2);
}

inline double CKKSPlainEvalModSineDegree3(double x)
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    return x * (1.0 - (two_pi * two_pi / 6.0) * x * x);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
struct CKKSEvalModSineDegree3Traits {
    static_assert(LogQ > 3 * LogDelta + CoeffLogDelta);
    static constexpr std::uint32_t x2_log_q = LogQ - LogDelta;
    static constexpr std::uint32_t x3_log_q = LogQ - 2 * LogDelta;
    static constexpr std::uint32_t term_input_log_q = x3_log_q;
    static constexpr std::uint32_t log_q = x3_log_q - CoeffLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
using CKKSEvalModSineDegree3Result =
    typename CKKSEvalModSineDegree3Traits<P, LogQ, LogDelta,
                                          CoeffLogDelta>::Ciphertext;

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
struct CKKSEvalModSineDegree3RelinKeys {
    using Traits =
        CKKSEvalModSineDegree3Traits<P, LogQ, LogDelta, CoeffLogDelta>;

    CKKSRelinKey<P, Traits::x2_log_q> x2{};
    CKKSRelinKey<P, Traits::x3_log_q> x3{};
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
inline void CKKSEvalModSineDegree3KeyGen(
    CKKSEvalModSineDegree3RelinKeys<P, LogQ, LogDelta, CoeffLogDelta> &keys,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    using Traits =
        CKKSEvalModSineDegree3Traits<P, LogQ, LogDelta, CoeffLogDelta>;
    keys.x2 = *makeCKKSRelinKey<P, Traits::x2_log_q>(key, noise);
    keys.x3 = *makeCKKSRelinKey<P, Traits::x3_log_q>(key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
inline void CKKSEvalModSineDegree3(
    CKKSEvalModSineDegree3Result<P, LogQ, LogDelta, CoeffLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSEvalModSineDegree3RelinKeys<P, LogQ, LogDelta, CoeffLogDelta>
        &keys)
{
    using Traits =
        CKKSEvalModSineDegree3Traits<P, LogQ, LogDelta, CoeffLogDelta>;
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    constexpr double c1 = -(two_pi * two_pi / 6.0);

    CKKSMultResult<P, LogQ, LogDelta, LogQ, LogDelta> x2;
    CKKSMult<P>(x2, ct, ct, keys.x2);

    CKKSMultResult<P, Traits::x2_log_q, LogDelta, LogQ, LogDelta> x3;
    CKKSMult<P>(x3, x2, ct, keys.x3);

    CKKSCiphertext<P, Traits::term_input_log_q, LogDelta> x_term_input;
    CKKSLevelReduce<P, LogQ, Traits::term_input_log_q, LogDelta>(
        x_term_input, ct);

    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x_term;
    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x3_term;
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x_term, x_term_input, 1.0);
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x3_term, x3, c1);

    CKKSAdd<P, Traits::log_q, LogDelta>(res, x_term, x3_term);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
struct CKKSEvalModSineDegree5Traits {
    static_assert(LogQ > 4 * LogDelta + CoeffLogDelta);
    static constexpr std::uint32_t x2_log_q = LogQ - LogDelta;
    static constexpr std::uint32_t x3_log_q = LogQ - 2 * LogDelta;
    static constexpr std::uint32_t x5_log_q = LogQ - 3 * LogDelta;
    static constexpr std::uint32_t term_input_log_q = x5_log_q;
    static constexpr std::uint32_t log_q = x5_log_q - CoeffLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
using CKKSEvalModSineDegree5Result =
    typename CKKSEvalModSineDegree5Traits<P, LogQ, LogDelta,
                                          CoeffLogDelta>::Ciphertext;

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
struct CKKSEvalModSineDegree5RelinKeys {
    using Traits =
        CKKSEvalModSineDegree5Traits<P, LogQ, LogDelta, CoeffLogDelta>;

    CKKSRelinKey<P, Traits::x2_log_q> x2{};
    CKKSRelinKey<P, Traits::x3_log_q> x3{};
    CKKSRelinKey<P, Traits::x5_log_q> x5{};
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
inline void CKKSEvalModSineDegree5KeyGen(
    CKKSEvalModSineDegree5RelinKeys<P, LogQ, LogDelta, CoeffLogDelta> &keys,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    using Traits =
        CKKSEvalModSineDegree5Traits<P, LogQ, LogDelta, CoeffLogDelta>;
    keys.x2 = *makeCKKSRelinKey<P, Traits::x2_log_q>(key, noise);
    keys.x3 = *makeCKKSRelinKey<P, Traits::x3_log_q>(key, noise);
    keys.x5 = *makeCKKSRelinKey<P, Traits::x5_log_q>(key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
inline void CKKSEvalModSineDegree5(
    CKKSEvalModSineDegree5Result<P, LogQ, LogDelta, CoeffLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSEvalModSineDegree5RelinKeys<P, LogQ, LogDelta, CoeffLogDelta>
        &keys)
{
    using Traits =
        CKKSEvalModSineDegree5Traits<P, LogQ, LogDelta, CoeffLogDelta>;
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    constexpr double c1 = -(two_pi * two_pi / 6.0);
    constexpr double c2 = two_pi * two_pi * two_pi * two_pi / 120.0;

    CKKSMultResult<P, LogQ, LogDelta, LogQ, LogDelta> x2;
    CKKSMult<P>(x2, ct, ct, keys.x2);

    CKKSMultResult<P, Traits::x2_log_q, LogDelta, LogQ, LogDelta> x3;
    CKKSMult<P>(x3, x2, ct, keys.x3);

    CKKSMultResult<P, Traits::x3_log_q, LogDelta, Traits::x2_log_q, LogDelta>
        x5;
    CKKSMult<P>(x5, x3, x2, keys.x5);

    CKKSCiphertext<P, Traits::term_input_log_q, LogDelta> x_term_input;
    CKKSCiphertext<P, Traits::term_input_log_q, LogDelta> x3_term_input;
    CKKSLevelReduce<P, LogQ, Traits::term_input_log_q, LogDelta>(
        x_term_input, ct);
    CKKSLevelReduce<P, Traits::x3_log_q, Traits::term_input_log_q, LogDelta>(
        x3_term_input, x3);

    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x_term;
    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x3_term;
    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x5_term;
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x_term, x_term_input, 1.0);
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x3_term, x3_term_input, c1);
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x5_term, x5, c2);

    CKKSAdd<P, Traits::log_q, LogDelta>(res, x_term, x3_term);
    CKKSAddInPlace<P, Traits::log_q, LogDelta>(res, x5_term);
}

}  // namespace TFHEpp
