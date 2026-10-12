#pragma once

#include <bit>
#include "prepared_dd.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Secondary DD rows are balanced signed digits, not ordinary torus
// ciphertext coefficients. Store each digit at its actual integer width.
// This explicit format is separate from CKKSDDRelinKey's legacy archive.
template <class P, std::uint32_t LogQ, std::uint32_t Bgbit = P::Bgbit,
          std::uint32_t Bbarbit = P::B̅gbit>
struct CKKSCompactDDRelinKey {
    static_assert(is_multilimb_uint_v<typename P::T>);
    static_assert(Bbarbit > 0 && Bbarbit <= 32);
    using Raw = CKKSDDRelinKey<P, LogQ, Bgbit, Bbarbit>;
    using Digit = std::conditional_t<(Bbarbit <= 8), std::int8_t,
                  std::conditional_t<(Bbarbit <= 16), std::int16_t, std::int32_t>>;
    static constexpr auto log_q = LogQ;
    static constexpr auto relin_bgbit = Bgbit, relin_bbarbit = Bbarbit;
    static constexpr auto primary_rows = Raw::primary_rows, bbar_rows = Raw::bbar_rows;
    using Polynomial = std::array<Digit, P::n>;
    using Row = std::array<Polynomial, P::k + 1>;
    std::array<std::array<Row, bbar_rows>, primary_rows> data{};

    static Digit pack(const typename P::T &value)
    {
        const auto digit = std::bit_cast<std::int64_t>(value.limb[0]);
        constexpr auto half = std::int64_t{1} << (Bbarbit - 1);
        if (digit < -half || digit >= half || P::T::from_signed_i64(digit) != value)
            throw std::invalid_argument("non-digit coefficient in compact DD key");
        return static_cast<Digit>(digit);
    }
    template <class Archive>
    void serialize(Archive &archive)
    {
        const std::array<std::uint64_t, 7> expected{
            0x434b4b5343444431ULL, P::n, P::k, LogQ, Bgbit, Bbarbit, sizeof(Digit)};
        auto metadata = expected;
        archive(metadata);
        if (metadata != expected) throw std::runtime_error("compact DD key metadata mismatch");
        archive(data);
    }
};

namespace ckks_detail {
template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
struct is_ckks_dd_relin_key<CKKSCompactDDRelinKey<P, Q, Bg, Bbar>> : std::true_type {};
}

template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
inline auto makeCKKSCompactDDRelinKey(const CKKSDDRelinKey<P, Q, Bg, Bbar> &raw)
{
    using Compact = CKKSCompactDDRelinKey<P, Q, Bg, Bbar>;
    auto compact = std::make_unique<Compact>();
    for (std::size_t i = 0; i < Compact::primary_rows; i++)
        for (std::size_t j = 0; j < Compact::bbar_rows; j++)
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t n = 0; n < P::n; n++)
                    compact->data[i][j][c][n] = Compact::pack(raw.data[i][j][c][n]);
    return compact;
}

// Expand only when explicitly requested; generation and switching below
// never allocate this full-width quadratic matrix.
template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
inline auto expandCKKSCompactDDRelinKey(const CKKSCompactDDRelinKey<P, Q, Bg, Bbar> &compact)
{
    using Compact = CKKSCompactDDRelinKey<P, Q, Bg, Bbar>;
    auto raw = std::make_unique<typename Compact::Raw>();
    for (std::size_t i = 0; i < Compact::primary_rows; i++)
        for (std::size_t j = 0; j < Compact::bbar_rows; j++)
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t n = 0; n < P::n; n++)
                    raw->data[i][j][c][n] = P::T::from_signed_i64(compact.data[i][j][c][n]);
    return raw;
}

// Select the existing encrypted rows whose gadget exponents match a wider
// primary base. This reuses coefficient keys without secret access or key
// generation. Wider primary digits change switching noise; callers must
// qualify accuracy before adopting the derived key.
template <std::uint32_t NewBg, class P, std::uint32_t Q,
          std::uint32_t OldBg, std::uint32_t Bbar>
inline auto makeCKKSCompactDDRelinKeyWithPrimaryBase(
    const CKKSCompactDDRelinKey<P, Q, OldBg, Bbar> &source)
{
    static_assert(NewBg >= OldBg && NewBg % OldBg == 0);
    using Source = CKKSCompactDDRelinKey<P, Q, OldBg, Bbar>;
    using Result = CKKSCompactDDRelinKey<P, Q, NewBg, Bbar>;
    static_assert(Result::primary_rows * NewBg <= std::numeric_limits<typename P::T>::digits);
    auto result = std::make_unique<Result>();
    for (std::size_t primary = 0; primary < Result::primary_rows; primary++) {
        const auto exponent = (Result::primary_rows - primary - 1) * NewBg;
        const auto source_row = Source::primary_rows - 1 - exponent / OldBg;
        result->data[primary] = source.data[source_row];
    }
    return result;
}

// Change only the public secondary representation. Reconstruct each encrypted
// primary row modulo Q and decompose it again; no secret or fresh encryption is
// needed, and the primary gadget, ciphertext modulus and encryption noise stay
// unchanged. The original compact key remains the serialized source of truth.
template <std::uint32_t NewBbar, class P, std::uint32_t Q,
          std::uint32_t Bg, std::uint32_t OldBbar>
inline auto makeCKKSCompactDDRelinKeyWithSecondaryBase(
    const CKKSCompactDDRelinKey<P, Q, Bg, OldBbar> &source)
{
    using Source = CKKSCompactDDRelinKey<P, Q, Bg, OldBbar>;
    using Result = CKKSCompactDDRelinKey<P, Q, Bg, NewBbar>;
    static_assert(Result::primary_rows == Source::primary_rows);
    static_assert(Result::bbar_rows * NewBbar <= std::numeric_limits<typename P::T>::digits);
    constexpr auto half = std::int64_t{1} << (OldBbar - 1);
    auto result = std::make_unique<Result>();
    auto ordinary = std::make_unique<TRLWE<P>>();
    auto digits = std::make_unique<std::array<TRLWE<P>, Result::bbar_rows>>();
    for (std::size_t primary = 0; primary < Source::primary_rows; primary++) {
        for (auto &component : *ordinary) component.fill(typename P::T{});
        for (std::size_t secondary = 0; secondary < Source::bbar_rows; secondary++) {
            const auto shift = (Source::bbar_rows - secondary - 1) * OldBbar;
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t i = 0; i < P::n; i++) {
                    const std::int64_t digit = source.data[primary][secondary][c][i];
                    if (digit < -half || digit >= half)
                        throw std::invalid_argument("invalid source secondary DD digit");
                    (*ordinary)[c][i] += P::T::from_signed_i64(digit) << shift;
                }
        }
        ckks_detail::reduceTRLWEToLevel<P, Q>(*ordinary);
        ckks_detail::activeBaseDecomposeTRLWERows<P, Q, NewBbar, Result::bbar_rows>(*digits, *ordinary);
        for (std::size_t secondary = 0; secondary < Result::bbar_rows; secondary++)
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t i = 0; i < P::n; i++)
                    result->data[primary][secondary][c][i] = Result::pack((*digits)[secondary][c][i]);
    }
    return result;
}

template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
inline void CKKSCompactDDRelinKeyGen(CKKSCompactDDRelinKey<P, Q, Bg, Bbar> &result,
    const Key<P> &secret, CKKSNoise noise = {P::α, 0})
{
    using Compact = CKKSCompactDDRelinKey<P, Q, Bg, Bbar>;
    static_assert(Bg + Bbar + P::nbit + 3 < std::numeric_limits<double>::digits);
    static_assert(Compact::primary_rows * Bg <= std::numeric_limits<typename P::T>::digits);
    static_assert(Compact::bbar_rows * Bbar <= std::numeric_limits<typename P::T>::digits);
    auto square = std::make_unique<Polynomial<P>>();
    auto part = std::make_unique<Polynomial<P>>();
    for (std::size_t n = 0; n < P::n; n++) (*part)[n] = secret[n];
    PolyMulDigit<P>(*square, *part, *part);
    auto gadget = std::make_unique<Polynomial<P>>();
    auto ordinary = std::make_unique<TRLWE<P>>();
    auto digits = std::make_unique<std::array<TRLWE<P>, Compact::bbar_rows>>();
    for (std::size_t i = 0; i < Compact::primary_rows; i++) {
        const auto shift = (Compact::primary_rows - i - 1) * Bg;
        for (std::size_t n = 0; n < P::n; n++)
            (*gadget)[n] = ckks_detail::reduceToLevel<P, Q>((*square)[n] << shift);
        ckks_detail::encryptPolynomialAtLevel<P, Q>(*ordinary, *gadget, secret, noise);
        ckks_detail::activeBaseDecomposeTRLWERows<P, Q, Bbar, Compact::bbar_rows>(*digits, *ordinary);
        for (std::size_t j = 0; j < Compact::bbar_rows; j++)
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t n = 0; n < P::n; n++)
                    result.data[i][j][c][n] = Compact::pack((*digits)[j][c][n]);
    }
}

template <class P, std::uint32_t Q, std::uint32_t Bg = P::Bgbit,
          std::uint32_t Bbar = P::B̅gbit>
inline auto makeCKKSCompactDDRelinKey(const Key<P> &secret, CKKSNoise noise = {P::α, 0})
{
    auto result = std::make_unique<CKKSCompactDDRelinKey<P, Q, Bg, Bbar>>();
    CKKSCompactDDRelinKeyGen(*result, secret, noise);
    return result;
}

template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
inline auto makeCKKSDDRelinKeyFFT(const CKKSCompactDDRelinKey<P, Q, Bg, Bbar> &compact)
{
    using Compact = CKKSCompactDDRelinKey<P, Q, Bg, Bbar>;
    static_assert(use_multilimb_digit_fft_v<P>);
    static_assert(Bg + Bbar + P::nbit + std::bit_width(Compact::primary_rows) + 3 <
                  std::numeric_limits<double>::digits,
                  "CKKS DD relin FFT accumulation needs mantissa headroom");
    auto result = std::make_unique<CKKSDDRelinKeyFFT<P, Q, Bg, Bbar>>();
    auto expanded = std::make_unique<Polynomial<P>>();
    for (std::size_t i = 0; i < Compact::primary_rows; i++)
        for (std::size_t j = 0; j < Compact::bbar_rows; j++)
            for (std::size_t c = 0; c <= P::k; c++) {
                for (std::size_t n = 0; n < P::n; n++)
                    (*expanded)[n] = P::T::from_signed_i64(compact.data[i][j][c][n]);
                TwistIFFTDigit<P>(result->data[i][j][c], *expanded);
            }
    return result;
}

template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
inline void CKKSRelinKeySwitch(TRLWE<P> &result, const Polynomial<P> &input,
                              const CKKSCompactDDRelinKey<P, Q, Bg, Bbar> &key)
{
    using Compact = CKKSCompactDDRelinKey<P, Q, Bg, Bbar>;
    static_assert(Bg + Bbar + P::nbit + 3 < std::numeric_limits<double>::digits);
    auto primary = std::make_unique<std::array<Polynomial<P>, Compact::primary_rows>>();
    ckks_detail::activeBaseDecomposePolynomialRows<P, Q, Bg, Compact::primary_rows>(*primary, input);
    auto accum = std::make_unique<std::array<TRLWE<P>, Compact::bbar_rows>>();
    auto expanded = std::make_unique<Polynomial<P>>();
    auto product = std::make_unique<Polynomial<P>>();
    for (std::size_t i = 0; i < Compact::primary_rows; i++)
        for (std::size_t j = 0; j < Compact::bbar_rows; j++)
            for (std::size_t c = 0; c <= P::k; c++) {
                for (std::size_t n = 0; n < P::n; n++)
                    (*expanded)[n] = P::T::from_signed_i64(key.data[i][j][c][n]);
                PolyMulDigit<P>(*product, *expanded, (*primary)[i]);
                for (std::size_t n = 0; n < P::n; n++) (*accum)[j][c][n] += (*product)[n];
            }
    ckks_detail::activeBaseRecombineTRLWERows<P, Q, Bbar, Compact::bbar_rows>(result, *accum);
}

namespace ckks_detail {
template <class P, std::uint32_t StartQ, std::uint32_t Delta, std::uint32_t Bg,
          std::uint32_t Bbar, class Seq>
struct CompactDDChainTuple;
template <class P, std::uint32_t StartQ, std::uint32_t Delta, std::uint32_t Bg,
          std::uint32_t Bbar, std::size_t... Is>
struct CompactDDChainTuple<P, StartQ, Delta, Bg, Bbar, std::index_sequence<Is...>> {
    using type = std::tuple<CKKSCompactDDRelinKey<P, StartQ - (Is + 1) * Delta, Bg, Bbar>...>;
};
}

template <class P, std::uint32_t StartQ, std::uint32_t Delta, std::size_t Depth,
          std::uint32_t Bg = P::Bgbit, std::uint32_t Bbar = P::B̅gbit>
struct CKKSCompactDDRelinKeyChain {
    static_assert(Depth > 0 && StartQ > Depth * Delta);
    using Tuple = typename ckks_detail::CompactDDChainTuple<P, StartQ, Delta, Bg, Bbar,
        std::make_index_sequence<Depth>>::type;
    Tuple keys;
    template <std::size_t I> auto &get() { return std::get<I>(keys); }
    template <std::size_t I> const auto &get() const { return std::get<I>(keys); }
    template <class Archive> void serialize(Archive &archive) { archive(keys); }
};

template <std::uint32_t NewBg, class P, std::uint32_t StartQ,
          std::uint32_t Delta, std::size_t Depth, std::uint32_t OldBg,
          std::uint32_t Bbar>
inline auto makeCKKSCompactDDRelinKeyChainWithPrimaryBase(
    const CKKSCompactDDRelinKeyChain<P, StartQ, Delta, Depth, OldBg, Bbar> &source)
{
    auto result = std::make_unique<CKKSCompactDDRelinKeyChain<P, StartQ, Delta, Depth, NewBg, Bbar>>();
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        ((result->template get<Is>() = *makeCKKSCompactDDRelinKeyWithPrimaryBase<NewBg>(
            source.template get<Is>())), ...);
    }(std::make_index_sequence<Depth>{});
    return result;
}

template <class P, std::uint32_t StartQ, std::uint32_t Delta, std::size_t Depth,
          std::uint32_t Bg, std::uint32_t Bbar>
inline void CKKSDDRelinKeyChainGen(CKKSCompactDDRelinKeyChain<P, StartQ, Delta, Depth, Bg, Bbar> &chain,
    const Key<P> &secret, CKKSNoise noise = {P::α, 0})
{
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        (CKKSCompactDDRelinKeyGen(chain.template get<Is>(), secret, noise), ...);
    }(std::make_index_sequence<Depth>{});
}

template <class P, std::uint32_t StartQ, std::uint32_t Delta, std::size_t Depth,
          std::uint32_t Bg, std::uint32_t Bbar>
inline auto makeCKKSPreparedDDRelinKeyChainFFT(
    const CKKSCompactDDRelinKeyChain<P, StartQ, Delta, Depth, Bg, Bbar> &compact)
{
    static_assert(Bg + Bbar + P::nbit + std::bit_width(
        CKKSDDRelinPrimaryRowCountForLevel<P, StartQ - Delta, Bg>()) + 3 <
        std::numeric_limits<double>::digits);
    auto result = std::make_unique<CKKSPreparedDDRelinKeyChainFFT<P, StartQ, Delta, Depth, Bg, Bbar>>();
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        ((std::get<Is>(result->keys) = makeCKKSDDRelinKeyFFT(compact.template get<Is>())), ...);
    }(std::make_index_sequence<Depth>{});
    return result;
}

// Prepare a derived FFT chain at a different secondary base while retaining
// the original serialized compact keys. Only one converted coefficient key
// is live at a time; all primary encrypted rows remain identical modulo Q.
template <std::uint32_t NewBbar, class P, std::uint32_t StartQ,
          std::uint32_t Delta, std::size_t Depth, std::uint32_t Bg,
          std::uint32_t OldBbar>
inline auto makeCKKSPreparedDDRelinKeyChainFFTWithSecondaryBase(
    const CKKSCompactDDRelinKeyChain<P, StartQ, Delta, Depth, Bg, OldBbar> &source)
{
    using Cache = CKKSPreparedDDRelinKeyChainFFT<P, StartQ, Delta, Depth, Bg, NewBbar>;
    auto result = std::make_unique<Cache>();
    auto prepare = [&]<std::size_t I>() {
        auto converted = makeCKKSCompactDDRelinKeyWithSecondaryBase<NewBbar>(source.template get<I>());
        std::get<I>(result->keys) = makeCKKSDDRelinKeyFFT(*converted);
    };
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        (prepare.template operator()<Is>(), ...);
    }(std::make_index_sequence<Depth>{});
    return result;
}
} // namespace TFHEpp
#endif
