#pragma once

#include <array>
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
};

template <class T>
struct is_multilimb_uint : std::false_type {};

template <std::size_t L>
struct is_multilimb_uint<MultiLimbUInt<L>> : std::true_type {};

template <class T>
inline constexpr bool is_multilimb_uint_v = is_multilimb_uint<T>::value;

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
