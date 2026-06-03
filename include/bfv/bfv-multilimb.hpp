#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace TFHEpp {

template <std::size_t Limbs>
struct MultiLimbUInt {
    static_assert(Limbs > 0);

    std::array<uint64_t, Limbs> limb{};

    static constexpr std::size_t limbs = Limbs;
    static constexpr std::size_t digits = 64 * Limbs;

    constexpr MultiLimbUInt() = default;

    template <class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
    constexpr MultiLimbUInt(I value)
    {
        if constexpr (std::is_signed_v<I>) {
            if (value < 0) limb.fill(std::numeric_limits<uint64_t>::max());
        }
        limb[0] = static_cast<uint64_t>(value);
    }

    static constexpr MultiLimbUInt max_value()
    {
        MultiLimbUInt res;
        res.limb.fill(std::numeric_limits<uint64_t>::max());
        return res;
    }

    static constexpr MultiLimbUInt from_signed_i64(int64_t value)
    {
        return MultiLimbUInt(value);
    }

    constexpr explicit operator uint64_t() const { return limb[0]; }

    constexpr bool bit(std::size_t i) const
    {
        return i < digits && ((limb[i / 64] >> (i % 64)) & 1U) != 0;
    }

    constexpr MultiLimbUInt &operator+=(const MultiLimbUInt &rhs)
    {
        uint64_t carry = 0;
        for (std::size_t i = 0; i < Limbs; i++) {
            const uint64_t a = limb[i];
            const uint64_t s1 = a + rhs.limb[i];
            const uint64_t c1 = s1 < a;
            const uint64_t s2 = s1 + carry;
            const uint64_t c2 = s2 < s1;
            limb[i] = s2;
            carry = c1 | c2;
        }
        return *this;
    }

    constexpr MultiLimbUInt &operator-=(const MultiLimbUInt &rhs)
    {
        uint64_t borrow = 0;
        for (std::size_t i = 0; i < Limbs; i++) {
            const uint64_t a = limb[i];
            const uint64_t b = rhs.limb[i];
            const uint64_t d1 = a - b;
            const uint64_t b1 = a < b;
            const uint64_t d2 = d1 - borrow;
            const uint64_t b2 = d1 < borrow;
            limb[i] = d2;
            borrow = b1 | b2;
        }
        return *this;
    }

    constexpr MultiLimbUInt &operator&=(const MultiLimbUInt &rhs)
    {
        for (std::size_t i = 0; i < Limbs; i++) limb[i] &= rhs.limb[i];
        return *this;
    }

    constexpr MultiLimbUInt &operator|=(const MultiLimbUInt &rhs)
    {
        for (std::size_t i = 0; i < Limbs; i++) limb[i] |= rhs.limb[i];
        return *this;
    }

    constexpr MultiLimbUInt &operator^=(const MultiLimbUInt &rhs)
    {
        for (std::size_t i = 0; i < Limbs; i++) limb[i] ^= rhs.limb[i];
        return *this;
    }

    constexpr MultiLimbUInt operator~() const
    {
        MultiLimbUInt res;
        for (std::size_t i = 0; i < Limbs; i++) res.limb[i] = ~limb[i];
        return res;
    }

    constexpr MultiLimbUInt operator-() const
    {
        MultiLimbUInt res = ~*this;
        res += MultiLimbUInt(1);
        return res;
    }

    constexpr MultiLimbUInt &operator<<=(std::size_t shift)
    {
        if (shift >= digits) {
            *this = {};
            return *this;
        }
        const std::size_t limb_shift = shift / 64;
        const std::size_t bit_shift = shift % 64;
        std::array<uint64_t, Limbs> out{};
        for (std::size_t i = Limbs; i-- > 0;) {
            if (i < limb_shift) continue;
            uint64_t v = limb[i - limb_shift] << bit_shift;
            if (bit_shift != 0 && i > limb_shift)
                v |= limb[i - limb_shift - 1] >> (64 - bit_shift);
            out[i] = v;
        }
        limb = out;
        return *this;
    }

    constexpr MultiLimbUInt &operator>>=(std::size_t shift)
    {
        if (shift >= digits) {
            *this = {};
            return *this;
        }
        const std::size_t limb_shift = shift / 64;
        const std::size_t bit_shift = shift % 64;
        std::array<uint64_t, Limbs> out{};
        for (std::size_t i = 0; i < Limbs; i++) {
            if (i + limb_shift >= Limbs) continue;
            uint64_t v = limb[i + limb_shift] >> bit_shift;
            if (bit_shift != 0 && i + limb_shift + 1 < Limbs)
                v |= limb[i + limb_shift + 1] << (64 - bit_shift);
            out[i] = v;
        }
        limb = out;
        return *this;
    }

    constexpr void add_limb_at(std::size_t pos, uint64_t value)
    {
        uint64_t carry = value;
        for (std::size_t i = pos; i < Limbs && carry != 0; i++) {
            const uint64_t old = limb[i];
            limb[i] += carry;
            carry = limb[i] < old;
        }
    }

    constexpr MultiLimbUInt &operator*=(const MultiLimbUInt &rhs)
    {
        MultiLimbUInt res;
        for (std::size_t i = 0; i < Limbs; i++) {
            for (std::size_t j = 0; j + i < Limbs; j++) {
                const unsigned __int128 prod =
                    static_cast<unsigned __int128>(limb[i]) * rhs.limb[j];
                res.add_limb_at(i + j, static_cast<uint64_t>(prod));
                if (i + j + 1 < Limbs)
                    res.add_limb_at(i + j + 1, static_cast<uint64_t>(prod >> 64));
            }
        }
        *this = res;
        return *this;
    }

    constexpr MultiLimbUInt div_u64(uint64_t divisor) const
    {
        MultiLimbUInt q;
        uint64_t rem = 0;
        for (std::size_t i = Limbs; i-- > 0;) {
            const unsigned __int128 cur =
                (static_cast<unsigned __int128>(rem) << 64) | limb[i];
            q.limb[i] = static_cast<uint64_t>(cur / divisor);
            rem = static_cast<uint64_t>(cur % divisor);
        }
        return q;
    }

    constexpr uint64_t mod_u64(uint64_t divisor) const
    {
        uint64_t rem = 0;
        for (std::size_t i = Limbs; i-- > 0;) {
            const unsigned __int128 cur =
                (static_cast<unsigned __int128>(rem) << 64) | limb[i];
            rem = static_cast<uint64_t>(cur % divisor);
        }
        return rem;
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(limb);
    }
};

template <class T>
struct is_multilimb_uint : std::false_type {};

template <std::size_t L>
struct is_multilimb_uint<MultiLimbUInt<L>> : std::true_type {};

template <class T>
inline constexpr bool is_multilimb_uint_v = is_multilimb_uint<T>::value;

namespace detail {

template <std::size_t Limbs>
constexpr bool limb_array_bit(const std::array<uint64_t, Limbs> &x,
                              std::size_t bit)
{
    return bit < 64 * Limbs && ((x[bit / 64] >> (bit % 64)) & 1U) != 0;
}

template <std::size_t Limbs>
constexpr void add_limb_at(std::array<uint64_t, Limbs> &x, std::size_t pos,
                           uint64_t value)
{
    uint64_t carry = value;
    for (std::size_t i = pos; i < Limbs && carry != 0; i++) {
        const uint64_t old = x[i];
        x[i] += carry;
        carry = x[i] < old;
    }
}

template <std::size_t Limbs>
constexpr void sub_limb_at(std::array<uint64_t, Limbs> &x, std::size_t pos,
                           uint64_t value)
{
    uint64_t borrow = value != 0 ? 1 : 0;
    uint64_t sub = value;
    for (std::size_t i = pos; i < Limbs && borrow != 0; i++) {
        const uint64_t old = x[i];
        x[i] -= sub;
        borrow = old < sub ? 1 : 0;
        sub = borrow;
    }
}

template <std::size_t Limbs>
constexpr void twos_complement(std::array<uint64_t, Limbs> &x)
{
    uint64_t carry = 1;
    for (std::size_t i = 0; i < Limbs; i++) {
        const uint64_t v = ~x[i];
        x[i] = v + carry;
        carry = carry != 0 && x[i] == 0 ? 1 : 0;
    }
}

template <std::size_t Limbs>
constexpr void shift_left_one(std::array<uint64_t, Limbs> &x)
{
    uint64_t carry = 0;
    for (std::size_t i = 0; i < Limbs; i++) {
        const uint64_t next = x[i] >> 63;
        x[i] = (x[i] << 1) | carry;
        carry = next;
    }
}

template <std::size_t OutLimbs>
constexpr bool rem_ge_divisor(
    const std::array<uint64_t, OutLimbs + 1> &rem,
    const MultiLimbUInt<OutLimbs> &divisor)
{
    if (rem[OutLimbs] != 0) return true;
    for (std::size_t i = OutLimbs; i-- > 0;) {
        if (rem[i] != divisor.limb[i]) return rem[i] > divisor.limb[i];
    }
    return true;
}

template <std::size_t OutLimbs>
constexpr void rem_sub_divisor(std::array<uint64_t, OutLimbs + 1> &rem,
                               const MultiLimbUInt<OutLimbs> &divisor)
{
    uint64_t borrow = 0;
    for (std::size_t i = 0; i < OutLimbs; i++) {
        const uint64_t a = rem[i];
        const uint64_t b = divisor.limb[i];
        const uint64_t d1 = a - b;
        const uint64_t b1 = a < b;
        const uint64_t d2 = d1 - borrow;
        const uint64_t b2 = d1 < borrow;
        rem[i] = d2;
        borrow = b1 | b2;
    }
    rem[OutLimbs] -= borrow;
}

}  // namespace detail

template <std::size_t Limbs>
struct WideSignedLimbAccumulator {
    static_assert(Limbs > 0);

    std::array<uint64_t, Limbs> limb{};

    static constexpr std::size_t limbs = Limbs;
    static constexpr std::size_t digits = 64 * Limbs;

    constexpr bool is_negative() const { return (limb[Limbs - 1] >> 63) != 0; }

    constexpr void add(const WideSignedLimbAccumulator &rhs)
    {
        for (std::size_t i = 0; i < Limbs; i++)
            detail::add_limb_at(limb, i, rhs.limb[i]);
    }

    constexpr void add_shifted_u64(uint64_t value, int shift)
    {
        if (value == 0) return;
        if (shift < 0) {
            const int rshift = -shift;
            if (rshift >= 64) return;
            value >>= rshift;
            shift = 0;
            if (value == 0) return;
        }
        if (shift >= static_cast<int>(digits)) return;

        const std::size_t limb_offset = static_cast<std::size_t>(shift) / 64;
        const int bit_offset = shift % 64;
        detail::add_limb_at(limb, limb_offset, value << bit_offset);
        if (bit_offset != 0 && limb_offset + 1 < Limbs)
            detail::add_limb_at(limb, limb_offset + 1,
                                value >> (64 - bit_offset));
    }

    constexpr void sub_shifted_u64(uint64_t value, int shift)
    {
        if (value == 0) return;
        if (shift < 0) {
            const int rshift = -shift;
            if (rshift >= 64) return;
            value >>= rshift;
            shift = 0;
            if (value == 0) return;
        }
        if (shift >= static_cast<int>(digits)) return;

        const std::size_t limb_offset = static_cast<std::size_t>(shift) / 64;
        const int bit_offset = shift % 64;
        detail::sub_limb_at(limb, limb_offset, value << bit_offset);
        if (bit_offset != 0 && limb_offset + 1 < Limbs)
            detail::sub_limb_at(limb, limb_offset + 1,
                                value >> (64 - bit_offset));
    }

    constexpr void add_shifted_i64(int64_t value, int shift)
    {
        if (value == 0) return;
        if (shift < 0) {
            const int rshift = -shift;
            if (rshift >= 63)
                value = value < 0 ? -1 : 0;
            else
                value >>= rshift;
            shift = 0;
            if (value == 0) return;
        }

        const bool negative = value < 0;
        const uint64_t magnitude =
            negative ? uint64_t{0} - static_cast<uint64_t>(value)
                     : static_cast<uint64_t>(value);
        if (negative)
            sub_shifted_u64(magnitude, shift);
        else
            add_shifted_u64(magnitude, shift);
    }

    template <std::size_t OutLimbs>
    constexpr MultiLimbUInt<OutLimbs> div_unsigned_to_torus(
        const MultiLimbUInt<OutLimbs> &divisor) const
    {
        MultiLimbUInt<OutLimbs> quotient;
        std::array<uint64_t, OutLimbs + 1> rem{};

        for (std::size_t bit = digits; bit-- > 0;) {
            detail::shift_left_one(rem);
            if (detail::limb_array_bit(limb, bit)) rem[0] |= 1;
            if (detail::rem_ge_divisor(rem, divisor)) {
                detail::rem_sub_divisor(rem, divisor);
                if (bit < MultiLimbUInt<OutLimbs>::digits)
                    quotient.limb[bit / 64] |= uint64_t{1} << (bit % 64);
            }
        }
        return quotient;
    }

    template <std::size_t OutLimbs>
    constexpr MultiLimbUInt<OutLimbs> div_to_torus(
        const MultiLimbUInt<OutLimbs> &divisor) const
    {
        if (!is_negative()) return div_unsigned_to_torus(divisor);

        WideSignedLimbAccumulator magnitude{limb};
        detail::twos_complement(magnitude.limb);
        return -magnitude.template div_unsigned_to_torus<OutLimbs>(divisor);
    }
};

template <std::size_t L>
constexpr uint64_t round_mul_u64_div_pow2(const MultiLimbUInt<L> &x,
                                          uint64_t multiplier)
{
    std::array<uint64_t, L + 1> product{};
    uint64_t carry = 0;
    for (std::size_t i = 0; i < L; i++) {
        const unsigned __int128 cur =
            static_cast<unsigned __int128>(x.limb[i]) * multiplier + carry;
        product[i] = static_cast<uint64_t>(cur);
        carry = static_cast<uint64_t>(cur >> 64);
    }
    product[L] = carry;

    const uint64_t old = product[L - 1];
    product[L - 1] += uint64_t{1} << 63;
    if (product[L - 1] < old) product[L]++;
    return product[L];
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator+(MultiLimbUInt<L> lhs,
                                     const MultiLimbUInt<L> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <std::size_t L, class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator+(MultiLimbUInt<L> lhs, I rhs)
{
    lhs += MultiLimbUInt<L>(rhs);
    return lhs;
}

template <class I, std::size_t L, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator+(I lhs, MultiLimbUInt<L> rhs)
{
    rhs += MultiLimbUInt<L>(lhs);
    return rhs;
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator-(MultiLimbUInt<L> lhs,
                                     const MultiLimbUInt<L> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template <std::size_t L, class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator-(MultiLimbUInt<L> lhs, I rhs)
{
    lhs -= MultiLimbUInt<L>(rhs);
    return lhs;
}

template <class I, std::size_t L, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator-(I lhs, const MultiLimbUInt<L> &rhs)
{
    MultiLimbUInt<L> res(lhs);
    res -= rhs;
    return res;
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator*(MultiLimbUInt<L> lhs,
                                     const MultiLimbUInt<L> &rhs)
{
    lhs *= rhs;
    return lhs;
}

template <std::size_t L, class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator*(MultiLimbUInt<L> lhs, I rhs)
{
    lhs *= MultiLimbUInt<L>(rhs);
    return lhs;
}

template <class I, std::size_t L, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator*(I lhs, MultiLimbUInt<L> rhs)
{
    rhs *= MultiLimbUInt<L>(lhs);
    return rhs;
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator&(MultiLimbUInt<L> lhs,
                                     const MultiLimbUInt<L> &rhs)
{
    lhs &= rhs;
    return lhs;
}

template <std::size_t L, class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator&(MultiLimbUInt<L> lhs, I rhs)
{
    lhs &= MultiLimbUInt<L>(rhs);
    return lhs;
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator|(MultiLimbUInt<L> lhs,
                                     const MultiLimbUInt<L> &rhs)
{
    lhs |= rhs;
    return lhs;
}

template <std::size_t L, class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator|(MultiLimbUInt<L> lhs, I rhs)
{
    lhs |= MultiLimbUInt<L>(rhs);
    return lhs;
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator^(MultiLimbUInt<L> lhs,
                                     const MultiLimbUInt<L> &rhs)
{
    lhs ^= rhs;
    return lhs;
}

template <std::size_t L, class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr MultiLimbUInt<L> operator^(MultiLimbUInt<L> lhs, I rhs)
{
    lhs ^= MultiLimbUInt<L>(rhs);
    return lhs;
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator<<(MultiLimbUInt<L> lhs, std::size_t shift)
{
    lhs <<= shift;
    return lhs;
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator>>(MultiLimbUInt<L> lhs, std::size_t shift)
{
    lhs >>= shift;
    return lhs;
}

template <std::size_t L>
constexpr MultiLimbUInt<L> operator/(const MultiLimbUInt<L> &lhs,
                                     uint64_t rhs)
{
    return lhs.div_u64(rhs);
}

template <std::size_t L>
constexpr uint64_t operator%(const MultiLimbUInt<L> &lhs, uint64_t rhs)
{
    return lhs.mod_u64(rhs);
}

template <std::size_t L>
constexpr bool operator==(const MultiLimbUInt<L> &a,
                          const MultiLimbUInt<L> &b)
{
    return a.limb == b.limb;
}

template <std::size_t L, class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr bool operator==(const MultiLimbUInt<L> &a, I b)
{
    return a == MultiLimbUInt<L>(b);
}

template <class I, std::size_t L, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr bool operator==(I a, const MultiLimbUInt<L> &b)
{
    return MultiLimbUInt<L>(a) == b;
}

template <std::size_t L>
constexpr bool operator!=(const MultiLimbUInt<L> &a,
                          const MultiLimbUInt<L> &b)
{
    return !(a == b);
}

template <std::size_t L, class I, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr bool operator!=(const MultiLimbUInt<L> &a, I b)
{
    return !(a == b);
}

template <class I, std::size_t L, std::enable_if_t<std::is_integral_v<I>, int> = 0>
constexpr bool operator!=(I a, const MultiLimbUInt<L> &b)
{
    return !(a == b);
}

template <std::size_t L>
constexpr bool operator<(const MultiLimbUInt<L> &a, const MultiLimbUInt<L> &b)
{
    for (std::size_t i = L; i-- > 0;) {
        if (a.limb[i] != b.limb[i]) return a.limb[i] < b.limb[i];
    }
    return false;
}

template <std::size_t L>
constexpr bool operator>(const MultiLimbUInt<L> &a, const MultiLimbUInt<L> &b)
{
    return b < a;
}

template <std::size_t L>
constexpr bool operator<=(const MultiLimbUInt<L> &a, const MultiLimbUInt<L> &b)
{
    return !(b < a);
}

template <std::size_t L>
constexpr bool operator>=(const MultiLimbUInt<L> &a, const MultiLimbUInt<L> &b)
{
    return !(a < b);
}

template <std::size_t L>
constexpr int64_t multilimb_to_signed_i64(const MultiLimbUInt<L> &x)
{
    return static_cast<int64_t>(x.limb[0]);
}

}  // namespace TFHEpp

namespace std {

template <std::size_t L>
class numeric_limits<TFHEpp::MultiLimbUInt<L>> {
public:
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = false;
    static constexpr bool is_integer = true;
    static constexpr bool is_exact = true;
    static constexpr bool has_infinity = false;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = true;
    static constexpr int digits = static_cast<int>(TFHEpp::MultiLimbUInt<L>::digits);
    static constexpr int digits10 = 0;
    static constexpr int max_digits10 = 0;
    static constexpr int radix = 2;

    static constexpr TFHEpp::MultiLimbUInt<L> min() noexcept { return {}; }
    static constexpr TFHEpp::MultiLimbUInt<L> lowest() noexcept { return {}; }
    static constexpr TFHEpp::MultiLimbUInt<L> max() noexcept
    {
        return TFHEpp::MultiLimbUInt<L>::max_value();
    }
};

}  // namespace std
