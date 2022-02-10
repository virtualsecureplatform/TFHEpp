#pragma once

#include <array>
#include <cloudkey.hpp>
#include <cstddef>
#include <cstdint>
#include <gate.hpp>
#include <limits>
#include <memory>
#include <type_traits>

#if __cplusplus == 202002L
#include <bit>
#endif

namespace TFHEpp {

namespace detail::integer {
template <typename P>
void full_adder(TLWE<P>& c1, TLWE<P>& s, const TLWE<P>& a, const TLWE<P>& b,
                const TLWE<P>& c0, const EvalKey& ek);
}

template <size_t N, typename TLWEPARAM>
class tfhe_uintN_t {
    static_assert(std::is_same_v<TLWEPARAM, lvl0param> ||
                  std::is_same_v<TLWEPARAM, lvlMparam>);

public:
    using TLWETYPE = TLWE<TLWEPARAM>;

    std::array<TLWETYPE, N> data;
    const EvalKey* ek;

private:
    template <typename T>
    static constexpr size_t clz(T n)  // count left zeros
    {
#if __cplusplus == 202002L
        return std::countl_zero(n);
#else
        return (n == 0) ? std::numeric_limits<T>::digits
                        : (!(n & static_cast<T>(1)
                                     << (std::numeric_limits<T>::digits - 1))
                               ? (1 + clz(static_cast<T>(n << 1)))
                               : 0);
#endif
    }

    template <typename T>
    static constexpr size_t msb(T n)
    {
        return std::numeric_limits<T>::digits - clz(n);
    }

    static TLWETYPE encrypt(const uint8_t x, const SecretKey* const sk)
    {
        typename TLWEPARAM::T mu = x ? TLWEPARAM::μ : -TLWEPARAM::μ;
        if (sk)
            return tlweSymEncrypt<TLWEPARAM>(mu, TLWEPARAM::α,
                                             sk->key.get<TLWEPARAM>());
        else {
            TLWETYPE r = {};
            r[TLWEPARAM::n] = mu;
            return r;
        }
    }

    static bool decrypt(const TLWETYPE& x, const SecretKey* const sk)
    {
        return tlweSymDecrypt<TLWEPARAM>(x, sk->key.get<TLWEPARAM>());
    }

public:
    tfhe_uintN_t() : ek() {}

    explicit tfhe_uintN_t(const EvalKey* const ek)
        : tfhe_uintN_t(0, nullptr, ek)
    {
    }

    ~tfhe_uintN_t() {}

    tfhe_uintN_t(const tfhe_uintN_t& x) : data(x.data), ek(x.ek) {}

    tfhe_uintN_t<N, TLWEPARAM>& operator=(const tfhe_uintN_t<N, TLWEPARAM>& x)
    {
        this->ek = x.ek;
        this->data = x.data;
        return *this;
    }

    tfhe_uintN_t(const uint64_t x, const SecretKey* const sk,
                 const EvalKey* const ek)
        : ek(ek)
    {
        uint64_t mask = 1;
        for (size_t i = 0; i < N; mask <<= 1, ++i) {
            data[i] = encrypt((x & mask) ? 1 : 0, sk);
        }
    }

    template <typename T>
    T decrypt(const SecretKey* const sk) const
    {
        static_assert(std::is_integral_v<T>);
        static_assert(!std::is_signed_v<T>);
        constexpr size_t d = std::numeric_limits<T>::digits;
        constexpr size_t width = N < d ? N : d;
        T res = 0, mask = 1;
        for (size_t i = 0; i < width; mask <<= 1, ++i) {
            if (decrypt(data[i], sk)) res |= mask;
        }
        return res;
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        for (auto& x : data) {
            archive(x);
        }
    }

    template <size_t M>
    void add(const tfhe_uintN_t<M, TLWEPARAM>& x,
             tfhe_uintN_t<N, TLWEPARAM>& res,
             tfhe_uintN_t<1, TLWEPARAM>& carry) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res_2;
        res.ek = carry.ek = this->ek;

        TLWETYPE zero;
        HomCONSTANTZERO<TLWEPARAM>(carry.data[0]);
        HomCONSTANTZERO<TLWEPARAM>(zero);
        for (size_t i = 0; i < N; ++i) {
            detail::integer::full_adder<TLWEPARAM>(
                carry.data[0], res_2.data[i], this->data[i],
                (i < M) ? x.data[i] : zero, carry.data[0], *ek);
        }
 
#pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            res.data[i] = res_2.data[i];
        }
    }

    void add(const uint64_t x, tfhe_uintN_t<N, TLWEPARAM>& res,
             tfhe_uintN_t<1, TLWEPARAM>& carry) const
    {
        res.ek = carry.ek = this->ek;

        TLWETYPE zero, one;
        HomCONSTANTZERO<TLWEPARAM>(carry.data[0]);
        HomCONSTANTZERO<TLWEPARAM>(zero);
        HomCONSTANTONE<TLWEPARAM>(one);

        bool skip = true;
        for (size_t i = 0; i < N; ++i) {
            const bool is_one = (i < 64) && (x & (1ull << i));
            if (skip && is_one) skip = false;

            if (skip)
                HomCOPY(res.data[i], this->data[i]);
            else {
                detail::integer::full_adder<TLWEPARAM>(
                    carry.data[0], res.data[i], is_one ? one : zero,
                    this->data[i], carry.data[0], *ek);
            }
        }
    }

    void sub(const tfhe_uintN_t<N, TLWEPARAM>& x,
             tfhe_uintN_t<N, TLWEPARAM>& res,
             tfhe_uintN_t<1, TLWEPARAM>& borrow) const
    {
        res.ek = borrow.ek = this->ek;

        HomCONSTANTONE<TLWEPARAM>(borrow.data[0]);
        TLWETYPE not_x;
        for (size_t i = 0; i < N; ++i) {
            HomNOT(not_x, x.data[i]);
            detail::integer::full_adder<TLWEPARAM>(borrow.data[0], res.data[i],
                                                   not_x, this->data[i],
                                                   borrow.data[0], *ek);
        }

        HomNOT(borrow.data[0], borrow.data[0]);
    }

    void sub(const uint64_t x, tfhe_uintN_t<N, TLWEPARAM>& res,
             tfhe_uintN_t<1, TLWEPARAM>& borrow) const
    {
        res.ek = borrow.ek = this->ek;

        TLWETYPE zero, one;
        HomCONSTANTONE<TLWEPARAM>(borrow.data[0]);
        HomCONSTANTZERO<TLWEPARAM>(zero);
        HomCONSTANTONE<TLWEPARAM>(one);

        bool skip = true;
        for (size_t i = 0; i < N; ++i) {
            bool is_zero = (i < 64) && (x & (1ull << i));
            if (skip && is_zero) skip = false;

            if (skip)
                HomCOPY(res.data[i], this->data[i]);
            else {
                detail::integer::full_adder<TLWEPARAM>(
                    borrow.data[0], res.data[i], is_zero ? zero : one,
                    this->data[i], borrow.data[0], *ek);
            }
        }

        HomNOT(borrow.data[0], borrow.data[0]);
    }

    // res = (*this) * x (res: M bits, *this: N bits, x: N bits)
    template <size_t M>
    void mul(const tfhe_uintN_t<N, TLWEPARAM>& x,
             tfhe_uintN_t<M, TLWEPARAM>& res) const
    {
        res.ek = this->ek;

        tfhe_uintN_t<M, TLWEPARAM> res_2(this->ek);
        tfhe_uintN_t<N, TLWEPARAM> carries(this->ek);
        constexpr size_t L = N < M ? M : N;
#pragma omp parallel for
        for (size_t i = 0; i < L; ++i) {
            if (i < M) HomCONSTANTZERO<TLWEPARAM>(res_2.data[i]);
            if (i < N) HomCONSTANTZERO<TLWEPARAM>(carries.data[i]);
        }

        constexpr size_t j_max = N <= M ? N : M;
        for (size_t j = 0; j < j_max; ++j) {
            const size_t i_max = N <= (M - j) ? N : M - j;
#pragma omp parallel for
            for (size_t i = 0; i < i_max; ++i) {
                TLWETYPE and_ij;
                // Generating Partial Product
                HomAND(and_ij, this->data[i], x.data[j], *ek);
                // Carry Saved Adder
                detail::integer::full_adder<TLWEPARAM>(
                    carries.data[i], res_2.data[i + j], and_ij,
                    res_2.data[i + j], carries.data[i], *ek);
            }
        }

        if (N < M) {
            TLWETYPE c;
            HomCONSTANTZERO<TLWEPARAM>(c);
            for (size_t i = 0; i < M - N; ++i) {
                detail::integer::full_adder<TLWEPARAM>(
                    c, res_2.data[N + i], carries.data[i], res_2.data[N + i], c,
                    *ek);
            }
        }

#pragma omp parallel for
        for (size_t i = 0; i < M; ++i) HomCOPY(res.data[i], res_2.data[i]);
    }

    template <size_t M>
    void mul(const uint64_t x, tfhe_uintN_t<M, TLWEPARAM>& res) const
    {
        res.ek = this->ek;

        int no_skip = 0;
        TLWETYPE zero, one;
        tfhe_uintN_t<M, TLWEPARAM> res_2(this->ek);
        tfhe_uintN_t<N, TLWEPARAM> carries(this->ek);
        constexpr size_t L = N < M ? M : N;
#pragma omp parallel for
        for (size_t i = 0; i < L; ++i) {
            if (i < M) HomCONSTANTZERO<TLWEPARAM>(res_2.data[i]);
            if (i < N) HomCONSTANTZERO<TLWEPARAM>(carries.data[i]);
        }
        HomCONSTANTZERO<TLWEPARAM>(zero);
        HomCONSTANTONE<TLWEPARAM>(one);

        constexpr size_t j_max = N <= M ? N : M;
        for (size_t j = 0; j < j_max; ++j) {
            const uint64_t j_mask = 1ull << j;
            const uint64_t j_masked_x = x & j_mask;
            if (!no_skip && j < 64 && j_masked_x) no_skip = 1;

            const size_t i_max = N <= M - j ? N : M - j;
#pragma omp parallel for
            for (size_t i = 0; i < i_max; ++i) {
                TLWETYPE and_ij;

                if (no_skip) {
                    if (no_skip == 1)
                        res_2.data[i + j] = this->data[i];
                    else {
                        if (64 <= j)
                            HomCONSTANTZERO<TLWEPARAM>(and_ij);
                        else {
                            if (j_masked_x)
                                HomAND(and_ij, this->data[i], one, *ek);
                            else
                                HomAND(and_ij, this->data[i], zero, *ek);
                        }
                        detail::integer::full_adder<TLWEPARAM>(
                            carries.data[i], res_2.data[i + j], and_ij,
                            res_2.data[i + j], carries.data[i], *ek);
                    }
                }
            }

            if (no_skip == 1) no_skip = 2;
        }

        if (N < M) {
            TLWETYPE c;
            HomCONSTANTZERO<TLWEPARAM>(c);
            for (size_t i = 0; i < M - N; ++i) {
                detail::integer::full_adder<TLWEPARAM>(
                    c, res_2.data[N + i], carries.data[i], res_2.data[N + i], c,
                    *ek);
            }
        }

#pragma omp parallel for
        for (size_t i = 0; i < M; ++i) HomCOPY(res.data[i], res_2.data[i]);
    }

    void div(const tfhe_uintN_t<N, TLWEPARAM>& x, tfhe_uintN_t<N, TLWEPARAM>& q,
             tfhe_uintN_t<N, TLWEPARAM>& r) const
    {
        q.ek = r.ek = this->ek;

        tfhe_uintN_t<N, TLWEPARAM> r_2;
        tfhe_uintN_t<N, TLWEPARAM> subtrahend;
#pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            HomNOT(subtrahend.data[i], x.data[i]);
            HomCONSTANTZERO<TLWEPARAM>(r_2.data[i]);
        }

        for (size_t j = 0; j < N; ++j) {
            for (size_t i = 0; i < N - 1; ++i) {
                HomCOPY(r_2.data[N - i - 1], r_2.data[N - i - 2]);
            }
            r_2.data[0] = this->data[N - j - 1];  // r is minuend

            tfhe_uintN_t<N, TLWEPARAM> tmp(ek);
            TLWETYPE ib;
            HomCONSTANTONE<TLWEPARAM>(ib);
            for (size_t i = 0; i < N; ++i) {
                detail::integer::full_adder<TLWEPARAM>(
                    ib, tmp.data[i], r_2.data[i], subtrahend.data[i], ib, *ek);
            }
#pragma omp parallel for
            for (size_t i = 0; i < N; ++i) {
                HomMUX(r_2.data[i], ib, tmp.data[i], r_2.data[i], *ek);
            }
            q.data[N - j - 1] = ib;
        }

#pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            HomCOPY(r.data[i], r_2.data[i]);
        }
    }

    void div(const uint64_t x, tfhe_uintN_t<N, TLWEPARAM>& q,
             tfhe_uintN_t<N, TLWEPARAM>& r) const
    {
        q.ek = r.ek = this->ek;

        tfhe_uintN_t<N, TLWEPARAM> r_2;
        tfhe_uintN_t<N, TLWEPARAM> subtrahend;
#pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            if (64 <= N || !(x & (1ull << i)))
                HomCONSTANTONE<TLWEPARAM>(subtrahend.data[i]);
            else
                HomCONSTANTZERO<TLWEPARAM>(subtrahend.data[i]);
            HomCONSTANTZERO<TLWEPARAM>(r_2.data[i]);
        }

        for (size_t j = 0; j < N; ++j) {
            for (size_t i = 0; i < N - 1; ++i) {
                HomCOPY(r_2.data[N - i - 1], r_2.data[N - i - 2]);
            }
            r_2.data[0] = this->data[N - j - 1];  // r is minuend

            tfhe_uintN_t<N, TLWEPARAM> tmp(ek);
            TLWETYPE ib;
            HomCONSTANTONE<TLWEPARAM>(ib);
            for (size_t i = 0; i < N; ++i) {
                detail::integer::full_adder<TLWEPARAM>(
                    ib, tmp.data[i], r_2.data[i], subtrahend.data[i], ib, *ek);
            }
#pragma omp parallel for
            for (size_t i = 0; i < N; ++i) {
                HomMUX(r_2.data[i], ib, tmp.data[i], r_2.data[i], *ek);
            }
            q.data[N - j - 1] = ib;
        }

#pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            HomCOPY(r.data[i], r_2.data[i]);
        }
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator+(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res;
        tfhe_uintN_t<1, TLWEPARAM> carry;
        this->add(x, res, carry);
        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator+(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res;
        tfhe_uintN_t<1, TLWEPARAM> carry;
        this->add(x, res, carry);
        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator-(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res;
        tfhe_uintN_t<1, TLWEPARAM> borrow;
        this->sub(x, res, borrow);
        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator-(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res;
        tfhe_uintN_t<1, TLWEPARAM> borrow;
        this->sub(x, res, borrow);
        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator*(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res;
        this->mul(x, res);
        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator*(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res;
        this->mul(x, res);
        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator/(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> q, r;
        this->div(x, q, r);
        return q;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator/(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> q, r;
        this->div(x, q, r);
        return q;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator%(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> q, r;
        this->div(x, q, r);
        return r;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator%(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> q, r;
        this->div(x, q, r);
        return r;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator~() const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(this->ek);
#pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            HomNOT(res.data[i], this->data[i]);
        }
        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator&(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(ek);

#pragma omp parallel for
        for (size_t i = 0; i < N; i++) {
            HomAND(res.data[i], this->data[i], x.data[i], *ek);
        }

        return res;
    }

    template <typename T = const tfhe_uintN_t<N, TLWEPARAM>>
    typename std::enable_if<N != 1, T>::type operator&(
        const tfhe_uintN_t<1, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(ek);

#pragma omp parallel for
        for (size_t i = 0; i < N; i++) {
            HomAND(res.data[i], this->data[i], x.data[0], *ek);
        }

        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator&(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(ek);

#pragma omp parallel for
        for (size_t i = 0; i < N; i++) {
            const bool is_one = (i < 64) && (x & (1ull << i));
            if (is_one)
                HomCOPY(res.data[i], this->data[i]);
            else
                HomCONSTANTZERO<TLWEPARAM>(res.data[i]);
        }

        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator|(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(ek);

#pragma omp parallel for
        for (size_t i = 0; i < N; i++) {
            HomOR(res.data[i], this->data[i], x.data[i], *ek);
        }

        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator|(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(ek);

#pragma omp parallel for
        for (size_t i = 0; i < N; i++) {
            const bool is_one = (i < 64) && (x & (1ull << i));
            if (is_one)
                HomCONSTANTONE<TLWEPARAM>(res.data[i]);
            else
                HomCOPY(res.data[i], this->data[i]);
        }

        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator^(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(ek);

#pragma omp parallel for
        for (size_t i = 0; i < N; i++) {
            HomXOR(res.data[i], this->data[i], x.data[i], *ek);
        }

        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator^(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(ek);

#pragma omp parallel for
        for (size_t i = 0; i < N; i++) {
            const bool is_one = (i < 64) && (x & (1ull << i));
            if (is_one)
                HomNOT(res.data[i], this->data[i]);
            else
                HomCOPY(res.data[i], this->data[i]);
        }

        return res;
    }

    template <size_t M>
    const tfhe_uintN_t<N, TLWEPARAM> operator<<(
        const tfhe_uintN_t<M, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(*this);

        TLWETYPE zero;
        HomCONSTANTZERO<TLWEPARAM>(zero);

        constexpr size_t n_cnt_bits = msb(N - 1) < M ? msb(N - 1) : M;
        tfhe_uintN_t<n_cnt_bits, TLWEPARAM> cnt;
#pragma omp parallel for
        for (size_t i = 0; i < n_cnt_bits; ++i) {
            cnt.data[i] = x.data[i];
        }
        size_t n = 1;
        constexpr size_t j_max = msb(N - 1) < M ? msb(N - 1) : M;
        for (size_t j = 0; j < j_max; n <<= 1, ++j) {
            for (size_t i = 0; i < N; ++i) {
                const size_t idx = N - i - 1;
                HomMUX(res.data[idx], cnt.data[j],
                       (idx < n) ? zero : res.data[idx - n], res.data[idx],
                       *ek);
            }
        }

        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator<<(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(*this);

        const size_t shift_cnt = (x < N) ? x : N - 1;
#pragma omp for
        for (size_t i = 0; i < N - shift_cnt; ++i) {
            res.data[N - i - 1] = res.data[N - i - 1 - shift_cnt];
        }
#pragma omp for
        for (size_t i = 0; i < shift_cnt; ++i) {
            HomCONSTANTZERO<TLWEPARAM>(res.data[i]);
        }

        return res;
    }

    template <size_t M>
    const tfhe_uintN_t<N, TLWEPARAM> operator>>(
        const tfhe_uintN_t<M, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(*this);

        TLWETYPE zero;
        HomCONSTANTZERO<TLWEPARAM>(zero);

        constexpr size_t n_cnt_bits = msb(N - 1) < M ? msb(N - 1) : M;
        tfhe_uintN_t<n_cnt_bits, TLWEPARAM> cnt;
#pragma omp parallel for
        for (size_t i = 0; i < n_cnt_bits; ++i) {
            cnt.data[i] = x.data[i];
        }

        size_t n = 1;
        for (size_t j = 0; j < n_cnt_bits; n <<= 1, ++j) {
            for (size_t i = 0; i < N; ++i) {
                HomMUX(res.data[i], cnt.data[j],
                       (N <= i + n) ? zero : res.data[i + n], res.data[i], *ek);
            }
        }

        return res;
    }

    const tfhe_uintN_t<N, TLWEPARAM> operator>>(const uint64_t x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res(*this);

        const size_t shift_cnt = (x < N) ? x : N - 1;
#pragma omp for
        for (size_t i = 0; i < N - shift_cnt; ++i) {
            res.data[i] = res.data[i + shift_cnt];
        }
#pragma omp for
        for (size_t i = 0; i < shift_cnt; ++i) {
            HomCONSTANTZERO<TLWEPARAM>(res.data[N - i - 1]);
        }

        return res;
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator==(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> xnor;

#pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            HomXNOR(xnor.data[i], this->data[i], x.data[i], *ek);
        }

        for (size_t i = 1; i < N; ++i) {
            HomAND(xnor.data[0], xnor.data[0], xnor.data[i], *ek);
        }

        tfhe_uintN_t<1, TLWEPARAM> res(ek);
        res.data[0] = xnor.data[0];
        return res;
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator==(const uint64_t x) const
    {
        tfhe_uintN_t<1, TLWEPARAM> res(ek);

        if (N < 64) {
            const uint64_t mask = 0xffff'ffff'ffff'ffff - ((1ull << N) - 1);
            if (x & mask) {
                HomCONSTANTZERO<TLWEPARAM>(res.data[0]);
                return res;
            }
        }

        tfhe_uintN_t<N, TLWEPARAM> xnor;
        TLWETYPE zero, one;
        HomCONSTANTZERO<TLWEPARAM>(zero);
        HomCONSTANTONE<TLWEPARAM>(one);

#pragma omp parallel for
        for (size_t i = 0; i < N; ++i) {
            const bool is_one = (i < 64) && (x & (1ull << i));
            HomXNOR(xnor.data[i], this->data[i], is_one ? one : zero, *ek);
        }

        for (size_t i = 1; i < N; ++i) {
            HomAND(xnor.data[0], xnor.data[0], xnor.data[i], *ek);
        }

        res.data[0] = xnor.data[0];
        return res;
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator!=(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        return ~(*this == x);
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator!=(const uint64_t x) const
    {
        return ~(*this == x);
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator<(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        tfhe_uintN_t<N, TLWEPARAM> res;
        tfhe_uintN_t<1, TLWEPARAM> borrow;
        this->sub(x, res, borrow);
        return borrow;
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator<(const uint64_t x) const
    {
        tfhe_uintN_t<1, TLWEPARAM> res(ek);

        if (N < 64) {
            const uint64_t mask = 0xffff'ffff'ffff'ffff - ((1ull << N) - 1);
            if (x & mask) {
                HomCONSTANTONE<TLWEPARAM>(res.data[0]);
                return res;
            }
        }

        tfhe_uintN_t<N, TLWEPARAM> tmp;
        this->sub(x, tmp, res);
        return res;
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator>=(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        return ~(*this < x);
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator>=(const uint64_t x) const
    {
        return ~(*this < x);
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator>(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        return (x < *this);
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator>(const uint64_t x) const
    {
        tfhe_uintN_t<1, TLWEPARAM> res(ek);
        if (N < 64) {
            const uint64_t mask = 0xffff'ffff'ffff'ffff - ((1ull << N) - 1);
            if (x & mask) {
                HomCONSTANTZERO<TLWEPARAM>(res.data[0]);
                return res;
            }
        }

        tfhe_uintN_t<N, TLWEPARAM> tmp;
        this->sub(x, tmp, res);
        for (size_t i = 1; i < N; ++i) {
            HomOR(tmp.data[0], tmp.data[0], tmp.data[i], *ek);
        }
        HomANDNY(res.data[0], res.data[0], tmp.data[0], *ek);
        return res;
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator<=(
        const tfhe_uintN_t<N, TLWEPARAM>& x) const
    {
        return ~(x < *this);
    }

    const tfhe_uintN_t<1, TLWEPARAM> operator<=(const uint64_t x) const
    {
        return ~(*this > x);
    }
};

template <typename T>
using tfhe_uint64_t = tfhe_uintN_t<64, T>;
template <typename T>
using tfhe_uint32_t = tfhe_uintN_t<32, T>;
template <typename T>
using tfhe_uint16_t = tfhe_uintN_t<16, T>;
template <typename T>
using tfhe_uint8_t = tfhe_uintN_t<8, T>;

}  // namespace TFHEpp
