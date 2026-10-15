#pragma once

#include "compact_dd.hpp"
#include <modular_ntt.hpp>
#include <span>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Unlike the FFT cache, this backend can multiply signed 64-bit primary
// digits by signed 32-bit secondary digits without mantissa rounding.
// A wider primary gadget still needs adequate ciphertext scale or auxiliary
// headroom; changing the arithmetic backend does not remove switching noise.
template <class P, std::uint32_t Q, std::uint32_t Bg = 64,
          std::uint32_t Bbar = 32>
struct CKKSPreparedDDRelinKeyNTT {
    static_assert(P::k == 1 && is_multilimb_uint_v<typename P::T>);
    static_assert(P::T::limbs >= 2 && Bg > 0 && Bg <= 64);
    static_assert(Bbar > 0 && Bbar <= 32);
    using Key = CKKSCompactDDRelinKey<P, Q, Bg, Bbar>;
    static constexpr auto log_q = Q, relin_bgbit = Bg, relin_bbarbit = Bbar;
    static constexpr std::size_t primary_rows = Key::primary_rows;
    static constexpr std::size_t secondary_rows = Key::bbar_rows;
    static_assert(primary_rows * Bg <= P::T::digits);
    static_assert(secondary_rows * Bbar <= P::T::digits);

    static constexpr bool fitsSinglePrime()
    {
        const auto limit = modular_ntt::degree65536_primes[0].value / 2;
        return (static_cast<unsigned __int128>(1) << (Bg-1)) <=
            (limit-1)/(static_cast<unsigned __int128>(1) << (Bbar-1))/P::n/primary_rows;
    }
    static constexpr std::size_t prime_count = fitsSinglePrime() ? 1 : 2;
    std::array<std::unique_ptr<modular_ntt::NegacyclicNTTPlan>, prime_count> plans;
    std::vector<modular_ntt::ShoupMultiplier> spectra;

    static constexpr bool fitsCRT()
    {
        const unsigned __int128 limit = static_cast<unsigned __int128>(
            modular_ntt::degree65536_primes[0].value) *
            modular_ntt::degree65536_primes[1].value / 2;
        return (static_cast<unsigned __int128>(1) << (Bg - 1)) <=
            (limit - 1) / (static_cast<unsigned __int128>(1) << (Bbar - 1)) /
            P::n / primary_rows;
    }
    static_assert(fitsCRT(), "DD digit convolution exceeds exact CRT range");

    static std::size_t index(std::size_t prime, std::size_t row,
        std::size_t secondary, std::size_t component, std::size_t coefficient)
    {
        return ((((prime * primary_rows + row) * secondary_rows + secondary) *
            2 + component) * P::n + coefficient);
    }
    explicit CKKSPreparedDDRelinKeyNTT(const Key &key)
        : spectra(prime_count * primary_rows * secondary_rows * 2 * P::n)
    {
        const auto &primes = modular_ntt::degree65536_primes;
        for (std::size_t prime = 0; prime < prime_count; ++prime)
            plans[prime] = std::make_unique<modular_ntt::NegacyclicNTTPlan>(
                P::n, primes[prime]);
#pragma omp parallel
        {
            std::vector<std::uint64_t> values(P::n);
#pragma omp for collapse(4) schedule(static)
            for (std::size_t prime = 0; prime < prime_count; ++prime)
                for (std::size_t row = 0; row < primary_rows; ++row)
                    for (std::size_t secondary = 0; secondary < secondary_rows; ++secondary)
                        for (std::size_t component = 0; component < 2; ++component) {
                            const auto modulus = primes[prime].value;
                            for (std::size_t i = 0; i < P::n; ++i) {
                                const std::int64_t value = key.data[row][secondary][component][i];
                                values[i] = value < 0 ? modulus - std::uint64_t(-value)
                                                     : std::uint64_t(value);
                            }
                            plans[prime]->forward(values);
                            for (std::size_t i = 0; i < P::n; ++i)
                                spectra[index(prime, row, secondary, component, i)] =
                                    modular_ntt::ShoupMultiplier(values[i], modulus);
                        }
        }
    }
    std::size_t byteSize() const { return spectra.size() * sizeof(spectra.front()); }

    static typename P::T signedTorus(__int128 value)
    {
        const bool negative = value < 0;
        const unsigned __int128 magnitude = negative
            ? static_cast<unsigned __int128>(0) - static_cast<unsigned __int128>(value)
            : static_cast<unsigned __int128>(value);
        typename P::T result{};
        result.limb[0] = std::uint64_t(magnitude);
        result.limb[1] = std::uint64_t(magnitude >> 64);
        return negative ? typename P::T{} - result : result;
    }
    void apply(TRLWE<P> &result, const Polynomial<P> &input) const
    {
        if (spectra.size() != prime_count * primary_rows * secondary_rows * 2 * P::n ||
            !plans[0] || (prime_count == 2 && !plans[1]))
            throw std::invalid_argument("invalid NTT DD cache");
        auto digits = std::make_unique<std::array<Polynomial<P>, primary_rows>>();
        ckks_detail::activeBaseDecomposePolynomialRows<P, Q, Bg, primary_rows>(
            *digits, input);
        std::vector<std::uint64_t> transformed(prime_count * primary_rows * P::n);
#pragma omp parallel for collapse(2) schedule(static)
        for (std::size_t prime = 0; prime < prime_count; ++prime)
            for (std::size_t row = 0; row < primary_rows; ++row) {
                const auto modulus = plans[prime]->modulus();
                auto values = std::span<std::uint64_t>(transformed.data() +
                    (prime * primary_rows + row) * P::n, P::n);
                for (std::size_t i = 0; i < P::n; ++i) {
                    const auto value = std::bit_cast<std::int64_t>((*digits)[row][i].limb[0]);
                    const auto magnitude = value < 0 ? std::uint64_t(0) - std::uint64_t(value)
                                                    : std::uint64_t(value);
                    const auto reduced = magnitude % modulus;
                    values[i] = value < 0 && reduced ? modulus - reduced : reduced;
                }
                plans[prime]->forward(values);
            }
        std::vector<std::uint64_t> products(prime_count * secondary_rows * 2 * P::n);
        const auto position = [](std::size_t prime, std::size_t secondary,
                                  std::size_t component, std::size_t i) {
            return ((prime * secondary_rows + secondary) * 2 + component) * P::n + i;
        };
#pragma omp parallel for collapse(3) schedule(static)
        for (std::size_t prime = 0; prime < prime_count; ++prime)
            for (std::size_t secondary = 0; secondary < secondary_rows; ++secondary)
                for (std::size_t component = 0; component < 2; ++component) {
                    const auto modulus = plans[prime]->modulus();
                    auto values = std::span<std::uint64_t>(products.data() +
                        position(prime, secondary, component, 0), P::n);
                    // Each lazy Shoup product is below 2p. Starting with a
                    // canonical sum, three products stay below 7p < 2^64.
                    // Reduce once per three rows, and traverse short contiguous
                    // coefficient tiles so the key streams stay cache friendly.
                    if (modulus < (std::uint64_t{1} << 61)) {
                        for (std::size_t begin = 0; begin < P::n; begin += 256)
                            for (std::size_t row = 0; row < primary_rows; row += 3)
                                for (std::size_t i = begin; i < std::min(begin+256,std::size_t(P::n)); ++i) {
                                    auto sum = values[i];
                                    for (std::size_t j = 0; j < 3 && row+j < primary_rows; ++j)
                                        sum += spectra[index(prime,row+j,secondary,component,i)].applyLazy(
                                            transformed[(prime*primary_rows+row+j)*P::n+i],modulus);
                                    if (sum >= 4*modulus) sum -= 4*modulus;
                                    if (sum >= 2*modulus) sum -= 2*modulus;
                                    if (sum >= modulus) sum -= modulus;
                                    values[i] = sum;
                                }
                    }
                    else {
                        for (std::size_t i = 0; i < P::n; ++i) {
                            std::uint64_t sum = 0;
                            for (std::size_t row = 0; row < primary_rows; ++row)
                                sum = modular_ntt::add(sum,spectra[index(prime,row,
                                    secondary,component,i)].apply(transformed[
                                        (prime*primary_rows+row)*P::n+i],modulus),modulus);
                            values[i] = sum;
                        }
                    }
                    plans[prime]->inverse(values);
                }
        const modular_ntt::TwoPrimeCRT crt(modular_ntt::degree65536_primes[0],
                                           modular_ntt::degree65536_primes[1]);
#pragma omp parallel for collapse(2) schedule(static)
        for (std::size_t component = 0; component < 2; ++component)
            for (std::size_t i = 0; i < P::n; ++i) {
                typename P::T sum{};
                for (std::size_t secondary = 0; secondary < secondary_rows; ++secondary)
                    {
                        __int128 exact;
                        const auto first = products[position(0,secondary,component,i)];
                        if constexpr (prime_count == 1) {
                            const auto modulus = modular_ntt::degree65536_primes[0].value;
                            exact = first <= modulus/2 ? static_cast<__int128>(first)
                                : -static_cast<__int128>(modulus-first);
                        }
                        else exact = crt.reconstructSigned(first,products[position(1,secondary,component,i)]);
                        sum += signedTorus(exact) << ((secondary_rows-secondary-1)*Bbar);
                    }
                result[component][i] = ckks_detail::reduceToLevel<P, Q>(sum);
            }
    }
};

namespace ckks_detail {
template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
struct is_ckks_dd_relin_key<CKKSPreparedDDRelinKeyNTT<P, Q, Bg, Bbar>>
    : std::true_type {};
}

template <std::uint32_t NewBbar = 32, class P, std::uint32_t Q,
          std::uint32_t Bg, std::uint32_t OldBbar>
inline auto makeCKKSDDRelinKeyNTT(const CKKSCompactDDRelinKey<P, Q, Bg, OldBbar> &key)
{
    using Cache = CKKSPreparedDDRelinKeyNTT<P, Q, Bg, NewBbar>;
    if constexpr (NewBbar == OldBbar) return std::make_shared<Cache>(key);
    else {
        auto converted = makeCKKSCompactDDRelinKeyWithSecondaryBase<NewBbar>(key);
        return std::make_shared<Cache>(*converted);
    }
}

template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
inline void CKKSRelinKeySwitch(TRLWE<P> &result, const Polynomial<P> &input,
    const CKKSPreparedDDRelinKeyNTT<P, Q, Bg, Bbar> &key)
{ key.apply(result, input); }

// Generation does not multiply primary and secondary digits. The existing
// generator's FFT headroom restriction remains intact; this entry point
// validates the integer NTT backend's bound instead.
template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
inline void CKKSCompactDDRelinKeyGenForNTT(CKKSCompactDDRelinKey<P, Q, Bg, Bbar> &result,
    const Key<P> &secret, CKKSNoise noise = {P::α, 0})
{
    using Compact = CKKSCompactDDRelinKey<P, Q, Bg, Bbar>;
    static_assert(CKKSPreparedDDRelinKeyNTT<P, Q, Bg, Bbar>::fitsCRT());
    auto square = std::make_unique<Polynomial<P>>();
    PolyMulDigit<P>(*square, secret, secret);
    auto gadget = std::make_unique<Polynomial<P>>();
    auto ordinary = std::make_unique<TRLWE<P>>();
    auto digits = std::make_unique<std::array<TRLWE<P>, Compact::bbar_rows>>();
    for (std::size_t row = 0; row < Compact::primary_rows; ++row) {
        const auto shift = (Compact::primary_rows - row - 1) * Bg;
        for (std::size_t i = 0; i < P::n; ++i)
            (*gadget)[i] = ckks_detail::reduceToLevel<P, Q>((*square)[i] << shift);
        ckks_detail::encryptPolynomialAtLevel<P, Q>(*ordinary, *gadget, secret, noise);
        ckks_detail::activeBaseDecomposeTRLWERows<P, Q, Bbar, Compact::bbar_rows>(
            *digits, *ordinary);
        for (std::size_t secondary = 0; secondary < Compact::bbar_rows; ++secondary)
            for (std::size_t component = 0; component < P::k + 1; ++component)
                for (std::size_t i = 0; i < P::n; ++i)
                    result.data[row][secondary][component][i] = Compact::pack(
                        (*digits)[secondary][component][i]);
    }
}

template <class P, std::uint32_t DataQ, std::uint32_t AuxiliaryLogScale,
          std::uint32_t Bg = 64, std::uint32_t Bbar = 32>
struct CKKSAuxiliaryDDRelinKeyNTT {
    static_assert(AuxiliaryLogScale > 0);
    using Cache = CKKSPreparedDDRelinKeyNTT<P, DataQ + AuxiliaryLogScale, Bg, Bbar>;
    static constexpr auto log_q = DataQ;
    std::shared_ptr<const Cache> cache;
    explicit CKKSAuxiliaryDDRelinKeyNTT(std::shared_ptr<const Cache> value)
        : cache(std::move(value))
    { if (!cache) throw std::invalid_argument("missing auxiliary NTT DD cache"); }
};

namespace ckks_detail {
template <class P, std::uint32_t Q, std::uint32_t Auxiliary,
          std::uint32_t Bg, std::uint32_t Bbar>
struct is_ckks_dd_relin_key<CKKSAuxiliaryDDRelinKeyNTT<P, Q, Auxiliary, Bg, Bbar>>
    : std::true_type {};
}

// Lift and rescale the nonlinear switched term, then let the standard
// relinearization add the tensor's unchanged linear terms at DataQ.
template <class P, std::uint32_t DataQ, std::uint32_t AuxiliaryLogScale,
          std::uint32_t Bg, std::uint32_t Bbar>
inline void CKKSRelinKeySwitch(TRLWE<P> &result, const Polynomial<P> &input,
    const CKKSAuxiliaryDDRelinKeyNTT<P, DataQ, AuxiliaryLogScale, Bg, Bbar> &key)
{
    constexpr auto key_q = DataQ + AuxiliaryLogScale;
    auto lifted = std::make_unique<Polynomial<P>>();
    for (std::size_t i = 0; i < P::n; ++i)
        (*lifted)[i] = ckks_detail::reduceToLevel<P, key_q>(input[i] << AuxiliaryLogScale);
    auto switched = std::make_unique<TRLWE<P>>();
    key.cache->apply(*switched, *lifted);
    for (std::size_t component = 0; component < 2; ++component)
        for (std::size_t i = 0; i < P::n; ++i)
            result[component][i] = ckks_detail::roundedLevelRightShift<P, key_q,
                AuxiliaryLogScale>((*switched)[component][i]);
    ckks_detail::reduceTRLWEToLevel<P, DataQ>(result);
}

} // namespace TFHEpp
#endif
