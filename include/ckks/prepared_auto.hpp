#pragma once

#include "ckks.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

template <class P, std::uint32_t LogQ>
constexpr std::uint32_t CKKSPreparedAutoSecondaryBits()
{
    constexpr auto rows = CKKSKeySwitchRowCountForLevel<P, LogQ>();
    if constexpr (P::B̅gbit < 16 && std::numeric_limits<typename P::T>::digits % 16 == 0 &&
        P::B̅gbit + 16 + P::nbit + std::bit_width(rows) + 3 < std::numeric_limits<double>::digits)
        return 16;
    else return P::B̅gbit;
}

// Default caches retain the raw primary gadget. An explicit wider primary
// radix selects a subset of existing rows; it changes switching noise and
// needs an accuracy check for the caller's scale/auxiliary-lift schedule.
template <class P, std::uint32_t LogQ, std::uint32_t SecondaryBits = CKKSPreparedAutoSecondaryBits<P, LogQ>(),
          std::uint32_t PrimaryBits = P::B̅gbit>
struct CKKSPreparedAutoKeyFFT {
    static_assert(LogQ > 0 && LogQ <= std::numeric_limits<typename P::T>::digits);
    static_assert(P::l̅ * P::B̅gbit == std::numeric_limits<typename P::T>::digits);
    static_assert(SecondaryBits > 0 && SecondaryBits <= 32);
    static_assert(std::numeric_limits<typename P::T>::digits % SecondaryBits == 0);
    static_assert(PrimaryBits > 0 && PrimaryBits <= 32);
    static_assert(PrimaryBits % P::B̅gbit == 0);
    static_assert(std::numeric_limits<typename P::T>::digits % PrimaryBits == 0);
    static constexpr std::uint32_t primary_bits = PrimaryBits;
    static constexpr std::size_t rows = (LogQ + PrimaryBits - 1) / PrimaryBits;
    static constexpr std::size_t first = std::numeric_limits<typename P::T>::digits / PrimaryBits - rows;
    static constexpr std::size_t secondary_first =
        (std::numeric_limits<typename P::T>::digits - LogQ) / SecondaryBits;
    static constexpr std::size_t secondary_rows =
        std::numeric_limits<typename P::T>::digits / SecondaryBits - secondary_first;
    static constexpr std::size_t spectrum_count = P::k * rows * secondary_rows * (P::k + 1);
    std::vector<PolynomialInFD<P>> spectra;

    const auto &at(std::size_t source, std::size_t primary,
                   std::size_t secondary, std::size_t component) const
    {
        return spectra[(((source * rows + primary) * secondary_rows + secondary) * (P::k + 1)) + component];
    }
};

template <class P, std::uint32_t LogQ, std::uint32_t SecondaryBits = CKKSPreparedAutoSecondaryBits<P, LogQ>(),
          std::uint32_t PrimaryBits = P::B̅gbit>
inline std::unique_ptr<CKKSPreparedAutoKeyFFT<P, LogQ, SecondaryBits, PrimaryBits>> makeCKKSPreparedAutoKeyFFT(
    const CKKSAutoKey<P, LogQ> &raw)
{
    static_assert(use_multilimb_digit_fft_v<P>);
    static_assert(P::l̅ * P::B̅gbit == std::numeric_limits<typename P::T>::digits);
    using Cache = CKKSPreparedAutoKeyFFT<P, LogQ, SecondaryBits, PrimaryBits>;
    static_assert(PrimaryBits + SecondaryBits + P::nbit + std::bit_width(Cache::rows) + 3 <
                  std::numeric_limits<double>::digits);
    using T = typename P::T;
    constexpr int width = std::numeric_limits<T>::digits;
    constexpr T half = T{1} << (SecondaryBits - 1);
    constexpr T mask = (T{1} << SecondaryBits) - T{1};
    constexpr T offset = [] {
        T value{};
        for (std::size_t j = 0; j < width / SecondaryBits; j++)
            value += (T{1} << (SecondaryBits - 1)) << (width - (j + 1) * SecondaryBits);
        return value;
    }();
    auto result = std::make_unique<Cache>();
    result->spectra.resize(Cache::spectrum_count);
#pragma omp parallel
    {
        auto digit = std::make_unique<Polynomial<P>>();
#pragma omp for collapse(3) schedule(static)
        for (std::size_t source = 0; source < P::k; source++)
            for (std::size_t row = 0; row < Cache::rows; row++)
                for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++) {
                    constexpr auto raw_first = CKKSKeySwitchFirstRowForLevel<P, LogQ>();
                    const auto raw_row = (Cache::first + row + 1) * PrimaryBits / P::B̅gbit - 1 - raw_first;
                    const auto shift = width - (Cache::secondary_first + secondary + 1) * SecondaryBits;
                    for (std::size_t c = 0; c <= P::k; c++) {
                        for (std::size_t i = 0; i < P::n; i++)
                            (*digit)[i] = (((raw[source][raw_row][c][i] + offset) >> shift) & mask) - half;
                        auto &spectrum = result->spectra[
                            (((source * Cache::rows + row) * Cache::secondary_rows + secondary) * (P::k + 1)) + c];
                        TwistIFFTDigit<P>(spectrum, *digit);
                    }
                }
    }
    return result;
}

namespace ckks_detail {
template <class P, std::uint32_t Bits>
struct AutoPrimaryParam : P {
    static constexpr std::uint32_t B̅gbit = Bits;
    static constexpr std::uint32_t l̅ = std::numeric_limits<typename P::T>::digits / Bits;
};
} // namespace ckks_detail

template <class P, std::uint32_t LogQ, std::uint32_t SecondaryBits, std::uint32_t PrimaryBits>
inline void CKKSPreparedKeySwitchRows(TRLWE<P> &result, const Polynomial<P> &poly,
    const CKKSPreparedAutoKeyFFT<P, LogQ, SecondaryBits, PrimaryBits> &key, std::size_t source)
{
    static_assert(use_multilimb_digit_fft_v<P>);
    using Cache = CKKSPreparedAutoKeyFFT<P, LogQ, SecondaryBits, PrimaryBits>;
    // Sum all primary rows before rounding the inverse FFT. Integer products
    // have ample mantissa headroom even after accumulation across those rows.
    static_assert(PrimaryBits + SecondaryBits + P::nbit + std::bit_width(Cache::rows) + 3 <
                  std::numeric_limits<double>::digits);
    if (source >= P::k || key.spectra.size() != Cache::spectrum_count)
        throw std::invalid_argument("invalid prepared auto key");
    auto centered = std::make_unique<Polynomial<P>>();
    ckks_detail::centeredPolynomialAtLevel<P, LogQ>(*centered, poly);
    auto digits = std::make_unique<std::array<Polynomial<P>, Cache::rows>>();
    ckks_detail::baseBbarDecomposePolynomialRows<ckks_detail::AutoPrimaryParam<P, PrimaryBits>, Cache::first>(*digits, *centered);
    std::vector<PolynomialInFD<P>> transformed(Cache::rows);
#pragma omp parallel for schedule(static) if(P::n >= 8192)
    for (std::size_t row = 0; row < Cache::rows; row++)
        TwistIFFTDigit<P>(transformed[row], (*digits)[row]);
    // FFT products fit signed64; expand only when recombining torus coefficients.
    std::vector<std::array<std::array<std::uint64_t, P::n>, P::k + 1>> products(Cache::secondary_rows);
#pragma omp parallel
    {
        auto accumulated = std::make_unique<PolynomialInFD<P>>();
#pragma omp for collapse(2) schedule(static)
        for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++)
            for (std::size_t c = 0; c <= P::k; c++) {
                accumulated->fill(0.0);
                for (std::size_t row = 0; row < Cache::rows; row++)
                    FMAInFD<P::n>(*accumulated, transformed[row], key.at(source, row, secondary, c));
                TwistFFTDigitProduct64<P>(products[secondary][c], *accumulated);
            }
    }
    if constexpr (P::n >= 8192) {
        // Each worker owns complete coefficients. Keep the secondary order,
        // but avoid streaming the full torus output once per secondary row.
#pragma omp parallel for collapse(2) schedule(static)
        for (std::size_t c = 0; c <= P::k; c++)
            for (std::size_t i = 0; i < P::n; i++) {
                typename P::T coefficient{};
                for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++) {
                    const auto shift = std::numeric_limits<typename P::T>::digits -
                        (Cache::secondary_first + secondary + 1) * SecondaryBits;
                    coefficient += static_cast<typename P::T>(static_cast<std::int64_t>(products[secondary][c][i])) << shift;
                }
                result[c][i] = ckks_detail::reduceToLevel<P, LogQ>(coefficient);
            }
    }
    else {
        for (auto &component : result) component.fill(typename P::T{});
        for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++) {
            const auto shift = std::numeric_limits<typename P::T>::digits -
                (Cache::secondary_first + secondary + 1) * SecondaryBits;
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t i = 0; i < P::n; i++)
                    result[c][i] += static_cast<typename P::T>(static_cast<std::int64_t>(products[secondary][c][i])) << shift;
        }
        ckks_detail::reduceTRLWEToLevel<P, LogQ>(result);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t SecondaryBits, std::uint32_t PrimaryBits>
inline void CKKSEvalAuto(TRLWE<P> &result, const TRLWE<P> &input, int exponent,
    const CKKSPreparedAutoKeyFFT<P, LogQ, SecondaryBits, PrimaryBits> &key)
{
    for (auto &component : result) component.fill(typename P::T{});
    Automorphism<P>(result[P::k], input[P::k], exponent);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(result[P::k]);
    auto mask = std::make_unique<Polynomial<P>>();
    auto switched = std::make_unique<TRLWE<P>>();
    for (std::size_t source = 0; source < P::k; source++) {
        Automorphism<P>(*mask, input[source], exponent);
        CKKSPreparedKeySwitchRows<P, LogQ>(*switched, *mask, key, source);
        CKKSSubTRLWEInPlace<P, LogQ>(result, *switched);
    }
}

// Apply the same two-axis decomposition as the prepared cache, but transform
// raw key digits on demand. Accumulating primary rows before the inverse FFT
// avoids one inverse transform and torus recomposition per primary row without
// retaining the quadratic-size key spectrum cache.
template <class P, std::uint32_t LogQ,
          std::uint32_t SecondaryBits = CKKSPreparedAutoSecondaryBits<P, LogQ>()>
inline void CKKSEvalAutoStreamedFFT(TRLWE<P> &result, const TRLWE<P> &input,
                                  int exponent, const CKKSAutoKey<P, LogQ> &key)
{
    static_assert(use_multilimb_digit_fft_v<P>);
    using Cache = CKKSPreparedAutoKeyFFT<P, LogQ, SecondaryBits>;
    using T = typename P::T;
    constexpr int width = std::numeric_limits<T>::digits;
    static_assert(P::l̅ * P::B̅gbit == width);
    static_assert(P::B̅gbit + SecondaryBits + P::nbit + std::bit_width(Cache::rows) + 3 <
                  std::numeric_limits<double>::digits);
    constexpr T half = T{1} << (SecondaryBits - 1);
    constexpr T mask = (T{1} << SecondaryBits) - T{1};
    constexpr T offset = [] {
        T value{};
        for (std::size_t j = 0; j < width / SecondaryBits; j++)
            value += (T{1} << (SecondaryBits - 1)) << (width - (j + 1) * SecondaryBits);
        return value;
    }();
    for (auto &component : result) component.fill(T{});
    Automorphism<P>(result[P::k], input[P::k], exponent);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(result[P::k]);
    auto centered = std::make_unique<Polynomial<P>>();
    auto digits = std::make_unique<std::array<Polynomial<P>, Cache::rows>>();
    std::vector<PolynomialInFD<P>> transformed(Cache::rows);
    std::vector<std::array<Polynomial<P>, P::k + 1>> products(Cache::secondary_rows);
    for (std::size_t source = 0; source < P::k; source++) {
        Automorphism<P>(*centered, input[source], exponent);
        ckks_detail::centeredPolynomialAtLevel<P, LogQ>(*centered, *centered);
        ckks_detail::baseBbarDecomposePolynomialRows<P, Cache::first>(*digits, *centered);
        for (std::size_t row = 0; row < Cache::rows; row++)
            TwistIFFTDigit<P>(transformed[row], (*digits)[row]);
#pragma omp parallel
        {
            auto digit = std::make_unique<Polynomial<P>>();
            auto spectrum = std::make_unique<PolynomialInFD<P>>();
            auto accumulated = std::make_unique<PolynomialInFD<P>>();
#pragma omp for collapse(2) schedule(static)
            for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++)
                for (std::size_t c = 0; c <= P::k; c++) {
                    accumulated->fill(0.0);
                    const auto shift = width - (Cache::secondary_first + secondary + 1) * SecondaryBits;
                    for (std::size_t row = 0; row < Cache::rows; row++) {
                        for (std::size_t i = 0; i < P::n; i++)
                            (*digit)[i] = (((key[source][row][c][i] + offset) >> shift) & mask) - half;
                        TwistIFFTDigit<P>(*spectrum, *digit);
                        FMAInFD<P::n>(*accumulated, transformed[row], *spectrum);
                    }
                    TwistFFTDigitProduct<P>(products[secondary][c], *accumulated);
                }
        }
        for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++) {
            const auto shift = width - (Cache::secondary_first + secondary + 1) * SecondaryBits;
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t i = 0; i < P::n; i++)
                    result[c][i] -= products[secondary][c][i] << shift;
        }
    }
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(result);
}

namespace ckks_detail {
template <class P, std::uint32_t LogQ, std::uint32_t AuxiliaryLogScale, class Operation>
inline void applyWithAuxiliaryScale(TRLWE<P> &result, const TRLWE<P> &input, Operation operation)
{
    static_assert(AuxiliaryLogScale > 0 && LogQ + AuxiliaryLogScale <= std::numeric_limits<typename P::T>::digits);
    auto lifted = std::make_unique<TRLWE<P>>();
    auto evaluated = std::make_unique<TRLWE<P>>();
    for (std::size_t c = 0; c <= P::k; c++)
        for (std::size_t i = 0; i < P::n; i++)
            (*lifted)[c][i] = reduceToLevel<P, LogQ>(input[c][i]) << AuxiliaryLogScale;
    operation(*evaluated, *lifted);
    for (std::size_t c = 0; c <= P::k; c++)
        for (std::size_t i = 0; i < P::n; i++)
            result[c][i] = roundedLevelRightShift<P, LogQ + AuxiliaryLogScale, AuxiliaryLogScale>((*evaluated)[c][i]);
}
}  // namespace ckks_detail

// Lift by a public power of two, switch at the enlarged modulus, then
// return to the original scale/level. Evaluation keys must be generated at
// LogQ+AuxiliaryLogScale; include that actual key modulus in security tuning.
// Input encryption error is preserved, while switch-key error is divided by
// 2^AuxiliaryLogScale. Final two-component rounding is still present.
template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t AuxiliaryLogScale, class AutoKey>
inline void CKKSEvalAutoWithAuxiliaryScale(
    CKKSCiphertext<P, LogQ, LogDelta> &result,
    const CKKSCiphertext<P, LogQ, LogDelta> &input,
    int exponent, const AutoKey &key)
{
    static_assert(AuxiliaryLogScale > 0);
    static_assert(LogDelta < LogQ);
    static_assert(LogQ + AuxiliaryLogScale <= std::numeric_limits<typename P::T>::digits);
    ckks_detail::applyWithAuxiliaryScale<P, LogQ, AuxiliaryLogScale>(result.ct, input.ct,
        [&](auto &out, const auto &in) { CKKSEvalAuto<P, LogQ + AuxiliaryLogScale>(out, in, exponent, key); });
}

template <class P, std::uint32_t LogQ>
struct CKKSPreparedRotationKeyFFT {
    static_assert(P::nbit > 1 && P::nbit <= 32);
    std::array<std::unique_ptr<CKKSPreparedAutoKeyFFT<P, LogQ>>, P::nbit - 1> bits;
    std::size_t byteSize() const
    {
        std::size_t bytes = 0;
        for (const auto &key : bits)
            if (key) bytes += key->spectra.size() * sizeof(PolynomialInFD<P>);
        return bytes;
    }
};

// Select only binary rotation keys needed by the caller's public workload.
template <class P, std::uint32_t LogQ>
inline auto makeCKKSPreparedRotationKeyFFT(const CKKSGaloisKey<P, LogQ> &raw,
                                         std::uint32_t bit_mask)
{
    using Cache = CKKSPreparedRotationKeyFFT<P, LogQ>;
    constexpr auto allowed = (std::uint32_t{1} << (P::nbit - 1)) - 1;
    if (bit_mask & ~allowed) throw std::invalid_argument("invalid prepared rotation bit mask");
    auto result = std::make_unique<Cache>();
    for (std::size_t bit = 0; bit < result->bits.size(); bit++)
        if ((bit_mask >> bit) & 1)
            result->bits[bit] = makeCKKSPreparedAutoKeyFFT<P, LogQ>(raw[bit]);
    return result;
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &result, const TRLWE<P> &input, int steps,
                            const CKKSPreparedRotationKeyFFT<P, LogQ> &key)
{
    constexpr int slots = P::n / 2;
    const auto normalized = ((steps % slots) + slots) % slots;
    for (std::size_t bit = 0; bit < key.bits.size(); bit++)
        if (((normalized >> bit) & 1) && !key.bits[bit])
            throw std::invalid_argument("missing prepared binary rotation key");
    ckks_detail::CKKSRotateSlotsBinaryChain<P, LogQ>(result, input, steps,
        [&](std::size_t bit) -> const CKKSPreparedAutoKeyFFT<P, LogQ> & {
            return *key.bits[bit];
        });
}

// Apply a batch of public rotations with one derived binary key resident.
// ciphertext_at(i) returns a stable TRLWE reference. Bits are applied in the
// same ascending order as CKKSRotateSlots; an auxiliary lift is divided only
// after the entire chain. Return the maximum resident key payload in bytes.
template <class P, std::uint32_t LogQ, std::uint32_t AuxiliaryLogScale = 0,
          class CiphertextProvider>
inline std::size_t CKKSRotateSlotsBatchWithBoundedFFT(
    std::span<const int> steps, CiphertextProvider &&ciphertext_at,
    const CKKSGaloisKey<P, LogQ + AuxiliaryLogScale> &raw)
{
    static_assert(use_multilimb_digit_fft_v<P>);
    static_assert(LogQ + AuxiliaryLogScale <= std::numeric_limits<typename P::T>::digits);
    constexpr int slots = P::n / 2;
    constexpr auto key_q = LogQ + AuxiliaryLogScale;
    std::vector<int> normalized(steps.size());
    std::uint32_t used = 0;
    for (std::size_t i = 0; i < steps.size(); i++) {
        normalized[i] = ((steps[i] % slots) + slots) % slots;
        used |= normalized[i];
    }
    if (steps.empty()) return 0;
    for (std::size_t i = 0; i < steps.size(); i++) {
        auto &ct = ciphertext_at(i);
        for (auto &component : ct)
            for (auto &value : component)
                value = ckks_detail::reduceToLevel<P, LogQ>(value) << AuxiliaryLogScale;
    }
    auto evaluated = std::make_unique<TRLWE<P>>();
    std::size_t peak_bytes = 0;
    std::uint64_t exponent = 5;
    for (std::size_t bit = 0; bit < P::nbit - 1; bit++) {
        if ((used >> bit) & 1) {
            auto cache = makeCKKSPreparedAutoKeyFFT<P, key_q>(raw[bit]);
            peak_bytes = std::max(peak_bytes, cache->spectra.size() * sizeof(PolynomialInFD<P>));
            for (std::size_t i = 0; i < steps.size(); i++)
                if ((normalized[i] >> bit) & 1) {
                    auto &ct = ciphertext_at(i);
                    CKKSEvalAuto<P, key_q>(*evaluated, ct, static_cast<int>(exponent), *cache);
                    ct = *evaluated;
                }
        }
        exponent = exponent * exponent % (2 * P::n);
    }
    for (std::size_t i = 0; i < steps.size(); i++) {
        auto &ct = ciphertext_at(i);
        if constexpr (AuxiliaryLogScale != 0) {
            for (auto &component : ct)
                for (auto &value : component)
                    value = ckks_detail::roundedLevelRightShift<P, key_q, AuxiliaryLogScale>(value);
        }
        else ckks_detail::reduceTRLWEToLevel<P, LogQ>(ct);
    }
    return peak_bytes;
}


// Evaluation-owned derived public cache. The coordinator performs rotations
// serially; each switch may use its existing internal OpenMP workers. Do not
// share this mutable cache between concurrent evaluations.
template <class P, std::uint32_t LogQ, std::uint32_t PrimaryBits=P::B̅gbit>
struct CKKSLazyPreparedGaloisKeyFFT {
    using Auto = CKKSPreparedAutoKeyFFT<P, LogQ,CKKSPreparedAutoSecondaryBits<P,LogQ>(),PrimaryBits>;
    const CKKSGaloisKey<P, LogQ> &raw;
    mutable std::array<std::unique_ptr<Auto>, P::nbit + 1> keys;
    mutable std::size_t prepared_keys = 0, peak_spectrum_bytes = 0;

    explicit CKKSLazyPreparedGaloisKeyFFT(const CKKSGaloisKey<P, LogQ> &source)
        : raw(source) {}
    std::size_t spectrum_bytes() const
    {
        std::size_t bytes = 0;
        for (const auto &key : keys)
            if (key) bytes += key->spectra.size() * sizeof(PolynomialInFD<P>);
        return bytes;
    }
    const Auto &get(std::size_t index) const
    {
        if (!(index < P::nbit - 1 || index == P::nbit))
            throw std::invalid_argument("invalid S2C prepared Galois key index");
        if (!keys[index]) {
            keys[index] = makeCKKSPreparedAutoKeyFFT<P,LogQ,CKKSPreparedAutoSecondaryBits<P,LogQ>(),PrimaryBits>(raw[index]);
            ++prepared_keys;
            peak_spectrum_bytes = std::max(peak_spectrum_bytes, spectrum_bytes());
        }
        return *keys[index];
    }
    // The BSGS evaluator already calls this hook after baby-step generation
    // and at output stage boundaries. Retain the raw serialized key only.
    void clear_cached_keys() const
    {
        for (auto &key : keys) key.reset();
    }
};

template <class P, std::uint32_t LogQ, std::uint32_t PrimaryBits=P::B̅gbit>
inline void CKKSRotateSlots(TRLWE<P> &result, const TRLWE<P> &input, int steps,
                           const CKKSLazyPreparedGaloisKeyFFT<P, LogQ,PrimaryBits> &key)
{
    ckks_detail::CKKSRotateSlotsBinaryChain<P, LogQ>(result, input, steps,
        [&](std::size_t bit) -> const typename CKKSLazyPreparedGaloisKeyFFT<P,LogQ,PrimaryBits>::Auto & {
            return key.get(bit);
        });
}

template <class P, std::uint32_t LogQ, std::uint32_t PrimaryBits=P::B̅gbit>
inline void CKKSConjugateSlots(TRLWE<P> &result, const TRLWE<P> &input,
                             const CKKSLazyPreparedGaloisKeyFFT<P, LogQ,PrimaryBits> &key)
{
    CKKSEvalAuto<P, LogQ>(result, input, 2 * P::n - 1, key.get(P::nbit));
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(result);
}

}  // namespace TFHEpp
#endif
