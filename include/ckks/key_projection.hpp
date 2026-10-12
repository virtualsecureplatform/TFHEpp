#pragma once

#include "ckks.hpp"
#include "compact_dd.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Public modulus reduction of a retained automorphism/secret-switch key.
// Gadget powers are anchored to the storage width, so the lower-level key
// uses a suffix of the same rows. Reduce coefficients, without rescaling:
// dividing them would also divide the encrypted gadget and invalidate the key.
// The encryption secret, automorphism and primary radix remain the same.
template <class P, std::uint32_t TargetLogQ, std::uint32_t SourceLogQ>
inline void CKKSReduceAutoKey(
    CKKSAutoKey<P, TargetLogQ> &result,
    const CKKSAutoKey<P, SourceLogQ> &source)
{
    static_assert(TargetLogQ <= SourceLogQ,
                  "retained evaluation keys can only be reduced to a lower modulus");
    static_assert(P::l̅ * P::B̅gbit == std::numeric_limits<typename P::T>::digits);
    constexpr auto offset = CKKSKeySwitchFirstRowForLevel<P, TargetLogQ>() -
                            CKKSKeySwitchFirstRowForLevel<P, SourceLogQ>();
    constexpr auto rows = CKKSAutoKey<P, TargetLogQ>::key_switch_rows;
    static_assert(offset + rows <= CKKSAutoKey<P, SourceLogQ>::key_switch_rows);
    for (std::size_t k = 0; k < P::k; ++k)
        for (std::size_t r = 0; r < rows; ++r)
            for (std::size_t c = 0; c <= P::k; ++c)
                for (std::size_t i = 0; i < P::n; ++i)
                    result[k][r][c][i] = ckks_detail::reduceToLevel<P, TargetLogQ>(
                        source[k][offset + r][c][i]);
}

template <class P, std::uint32_t TargetLogQ, std::uint32_t SourceLogQ>
inline void CKKSReduceGaloisKey(
    CKKSGaloisKey<P, TargetLogQ> &result,
    const CKKSGaloisKey<P, SourceLogQ> &source)
{
    for (std::size_t i = 0; i < source.size(); ++i)
        CKKSReduceAutoKey<P, TargetLogQ, SourceLogQ>(result[i], source[i]);
}

// Lower a retained compact DD key without fresh encryption. Primary gadget
// powers match a suffix of the source rows. Secondary digits are a signed
// representation of encrypted coefficients: recompose, reduce modulo TargetQ,
// and decompose again, rather than simply truncating their row array.
template <std::uint32_t TargetQ, class P, std::uint32_t SourceQ,
          std::uint32_t Bg, std::uint32_t Bbar>
inline auto makeCKKSCompactDDRelinKeyAtLevel(
    const CKKSCompactDDRelinKey<P, SourceQ, Bg, Bbar> &source)
{
    static_assert(TargetQ <= SourceQ,
                  "retained DD keys can only be reduced to a lower modulus");
    using Source = CKKSCompactDDRelinKey<P, SourceQ, Bg, Bbar>;
    using Result = CKKSCompactDDRelinKey<P, TargetQ, Bg, Bbar>;
    constexpr auto offset = Source::primary_rows - Result::primary_rows;
    constexpr auto half = std::int64_t{1} << (Bbar - 1);
    auto result = std::make_unique<Result>();
    auto ordinary = std::make_unique<TRLWE<P>>();
    auto digits = std::make_unique<std::array<TRLWE<P>, Result::bbar_rows>>();
    for (std::size_t primary = 0; primary < Result::primary_rows; ++primary) {
        for (auto &component : *ordinary) component.fill(typename P::T{});
        for (std::size_t secondary = 0; secondary < Source::bbar_rows; ++secondary) {
            const auto shift = (Source::bbar_rows - secondary - 1) * Bbar;
            for (std::size_t c = 0; c <= P::k; ++c)
                for (std::size_t i = 0; i < P::n; ++i) {
                    const std::int64_t digit = source.data[offset + primary][secondary][c][i];
                    if (digit < -half || digit >= half)
                        throw std::invalid_argument("invalid retained secondary DD digit");
                    (*ordinary)[c][i] += P::T::from_signed_i64(digit) << shift;
                }
        }
        ckks_detail::reduceTRLWEToLevel<P, TargetQ>(*ordinary);
        ckks_detail::activeBaseDecomposeTRLWERows<P, TargetQ, Bbar,
            Result::bbar_rows>(*digits, *ordinary);
        for (std::size_t secondary = 0; secondary < Result::bbar_rows; ++secondary)
            for (std::size_t c = 0; c <= P::k; ++c)
                for (std::size_t i = 0; i < P::n; ++i)
                    result->data[primary][secondary][c][i] =
                        Result::pack((*digits)[secondary][c][i]);
    }
    return result;
}

}  // namespace TFHEpp
#endif
