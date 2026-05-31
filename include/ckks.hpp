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

#include "bfv++.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
struct CKKSCiphertext {
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
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
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);
    static_assert(LogDelta < LogQ);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;

    Polynomial<P> poly{};
};

template <class P, std::uint32_t LogQ>
struct CKKSRelinKey {
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    std::array<TRLWE<P>, P::l̅> data{};

    TRLWE<P> &operator[](std::size_t i) { return data[i]; }
    const TRLWE<P> &operator[](std::size_t i) const { return data[i]; }
};

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
inline __int128_t torusToSigned(typename P::T value)
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    return static_cast<__int128_t>(value);
}

template <class P>
inline typename P::T signedToTorus(__int128_t value)
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    return static_cast<typename P::T>(value);
}

template <class P, std::uint32_t LogQ>
inline constexpr typename P::T levelMask()
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    static_assert(LogQ > 0);
    static_assert(LogQ <= torus_width_v<P>);
    if constexpr (LogQ == torus_width_v<P>)
        return std::numeric_limits<typename P::T>::max();
    else
        return (static_cast<typename P::T>(1) << LogQ) - 1;
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
    if constexpr (LogQ == torus_width_v<P>) {
        return torusToSigned<P>(value);
    }
    else {
        const typename P::T sign = static_cast<typename P::T>(1) << (LogQ - 1);
        if ((value & sign) == 0) return static_cast<__int128_t>(value);

        const typename P::T modulus = static_cast<typename P::T>(1) << LogQ;
        const typename P::T magnitude = modulus - value;
        return -static_cast<__int128_t>(magnitude);
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

template <class P, std::uint32_t LogQ>
inline void centeredTRLWEAtLevel(TRLWE<P> &out, const TRLWE<P> &in)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++)
            out[c][n] = signedToTorus<P>(levelToSigned<P, LogQ>(in[c][n]));
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
            std::ldexp(static_cast<long double>(levelToSigned<P, LogQ>(poly[i])),
                       -static_cast<int>(LogDelta));
        coeffs[i] = static_cast<double>(value);
    }
}

template <std::uint32_t LogDelta, class P>
inline void fillConjugateSymmetricEvaluations(
    std::array<std::complex<double>, P::n> &evals,
    const std::array<std::complex<double>, P::n / 2> &slots)
{
    const double scale = std::ldexp(1.0, LogDelta);
    for (std::uint32_t i = 0; i < P::n / 2; i++) {
        evals[i] = slots[i] * scale;
        evals[P::n - 1 - i] = std::conj(evals[i]);
    }
}

template <class P, std::uint32_t LogQ>
inline void encryptPolynomialAtLevel(TRLWE<P> &ct, const Polynomial<P> &poly,
                                     const Key<P> &key,
                                     CKKSNoise noise = {P::α, 0})
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");

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
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
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
            (*torus_digit)[n] = (((torus[n] + offset) >> shift) & mask) - half;

        TwistIFFT<P>(*torus_digit_fft, *torus_digit);
        MulInFD<P::n>(*product_fft, *torus_digit_fft, *digit_fft);
        TwistFFT<P>(*product, *product_fft);

        for (std::uint32_t n = 0; n < P::n; n++)
            res[n] += (*product)[n] << shift;
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline typename P::T ckksEncodeCoeff(double value)
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    static_assert(LogDelta < LogQ);
    const long double scaled =
        std::ldexp(static_cast<long double>(value), LogDelta);
    return ckks_detail::signedToLevel<P, LogQ>(
        static_cast<__int128_t>(std::round(scaled)));
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline double ckksDecodeCoeff(typename P::T value)
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    static_assert(LogDelta < LogQ);
    const long double decoded =
        std::ldexp(static_cast<long double>(
                       ckks_detail::levelToSigned<P, LogQ>(value)),
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
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    static_assert(LogDelta < LogQ);
    std::array<std::complex<double>, P::n> evals{};
    ckks_detail::fillConjugateSymmetricEvaluations<LogDelta, P>(evals, slots);

    constexpr double pi = 3.141592653589793238462643383279502884;
    std::array<double, P::n> coeffs{};
    for (std::uint32_t j = 0; j < P::n; j++) {
        std::complex<long double> sum = 0.0L;
        for (std::uint32_t k = 0; k < P::n; k++) {
            const long double angle =
                -pi * static_cast<long double>((2 * k + 1) * j) /
                static_cast<long double>(P::n);
            sum += static_cast<std::complex<long double>>(evals[k]) *
                   std::complex<long double>(std::cos(angle),
                                             std::sin(angle));
        }
        coeffs[j] = static_cast<double>(sum.real() /
                                        static_cast<long double>(P::n));
    }

    for (std::uint32_t i = 0; i < P::n; i++)
        poly[i] = ckks_detail::signedToLevel<P, LogQ>(
            static_cast<__int128_t>(std::round(coeffs[i])));
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotDecode(
    std::array<std::complex<double>, P::n / 2> &slots,
    const Polynomial<P> &poly)
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    static_assert(LogDelta < LogQ);
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double inv_scale = std::ldexp(1.0, -static_cast<int>(LogDelta));
    for (std::uint32_t k = 0; k < P::n / 2; k++) {
        std::complex<long double> value = 0.0L;
        for (std::uint32_t j = 0; j < P::n; j++) {
            const long double angle =
                pi * static_cast<long double>((2 * k + 1) * j) /
                static_cast<long double>(P::n);
            const long double coeff = static_cast<long double>(
                ckks_detail::levelToSigned<P, LogQ>(poly[j]));
            value += coeff *
                     std::complex<long double>(std::cos(angle),
                                               std::sin(angle));
        }
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

template <class P, std::uint32_t LogQ>
inline std::unique_ptr<CKKSRelinKey<P, LogQ>> makeCKKSRelinKey(
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");

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
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    static_assert(P::l̅ * P::B̅gbit ==
                      std::numeric_limits<typename P::T>::digits,
                  "CKKS relinearization requires full Bbar decomposition");

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
            ckks_detail::polyMulTorusByBbarDigit<P>(
                *product, relinkey[j][c], digit);
            for (std::uint32_t n = 0; n < P::n; n++)
                res[c][n] =
                    ckks_detail::reduceToLevel<P, LogQ>(res[c][n] +
                                                        (*product)[n]);
        }
    }
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
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
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

}  // namespace TFHEpp
