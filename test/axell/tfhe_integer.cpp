#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>
#include <type_traits>
#include <vector>

// using TP = TFHEpp::lvl0param;
using TP = TFHEpp::lvlMparam;

template <typename T>
T mask_bits(const T x, const size_t n)
{
    T ret = 0, mask = 1;
    for (size_t i = 0; i < n; ++i) {
        ret |= (x & mask);
        mask <<= 1;
    }
    return ret;
}

template <typename T>
T make_mask(const size_t n)
{
    T ret = 0;
    for (size_t i = 0; i < n; ++i) {
        ret |= (1UL << i);
    }
    return ret;
}

template <typename T, typename U>
void init_right_val(T& t, const U& u, const TFHEpp::SecretKey* const sk,
                    const TFHEpp::EvalKey* const ek)
{
    if constexpr (std::is_unsigned_v<T>)
        t = static_cast<T>(u);
    else
        t = T(u, sk, ek);
}

template <typename T, typename U>
void init_right_val(T& t, const U& u, const size_t bits,
                    const TFHEpp::SecretKey* const sk, const TFHEpp::EvalKey* const ek)
{
    if constexpr (std::is_unsigned_v<T>)
        t = static_cast<T>(u) & make_mask<T>(bits);
    else
        t = T(u, sk, ek);
}

template <typename T, typename U>
size_t check_mismatch(std::vector<T> results, std::vector<T> expects,
                      std::vector<U> input, const size_t n,
                      const size_t result_offset)
{
    size_t n_mismatch = 0;
    for (size_t i = 0; i < n; ++i) {
        const size_t ridx = result_offset + i;
        if (results[ridx] != expects[i]) {
            const size_t width = std::numeric_limits<T>::digits / 4;
            std::cout << std::hex << std::setfill('0') << "result mismatch: 0x"
                      << std::setw(width) << input[i] << " -> 0x"
                      << std::setw(width) << results[ridx] << " (expected: 0x"
                      << std::setw(width) << expects[i] << ")" << std::endl;
            ++n_mismatch;
        }
    }
    return n_mismatch;
}

template <typename T, typename U0, typename U1>
size_t check_mismatch(std::vector<T> results, std::vector<T> expects,
                      std::vector<U0> input_x, std::vector<U1> input_y,
                      const size_t n, const size_t result_offset)
{
    size_t n_mismatch = 0;
    for (size_t i = 0; i < n; ++i) {
        const size_t ridx = result_offset + i;
        if (results[ridx] != expects[i]) {
            const size_t width = std::numeric_limits<T>::digits / 4;
            std::cout << std::hex << std::setfill('0') << "result mismatch: (0x"
                      << std::setw(width) << input_x[i] << ", 0x"
                      << std::setw(width) << input_y[i] << ") -> 0x"
                      << std::setw(width) << results[ridx] << " (expected: 0x"
                      << std::setw(width) << expects[i] << ")" << std::endl;
            ++n_mismatch;
        }
    }
    return n_mismatch;
}

template <typename T, size_t left_bits, typename Tr, size_t right_bits,
          typename Trnd, typename F, typename Fc>
void test(const size_t num_test, F func, Fc func_tfhe,
          const bool shift_test = false)
{
    constexpr bool x_y_same_type =
        std::is_same_v<TFHEpp::tfhe_uintN_t<left_bits, TP>, Tr>;
    constexpr size_t NUM_CASE = x_y_same_type ? 5 : 2;

    std::random_device seed_gen;
    std::uniform_int_distribution<Trnd> rnd(0,
                                            std::numeric_limits<Trnd>::max());

    const Trnd l_mask = make_mask<Trnd>(left_bits);
    const Trnd r_mask = make_mask<Trnd>(right_bits);
    std::vector<Trnd> y_orig(num_test);
    std::vector<T> x(num_test), y(num_test), z(num_test), w(num_test),
        ans(NUM_CASE * num_test);
    std::vector<bool> z_cb(num_test, 0), w_cb(num_test, 0),
        ans_cb(NUM_CASE * num_test);  // carry, borrow
    for (size_t i = 0; i < num_test; ++i) {
        y_orig[i] = rnd(seed_gen);
        x[i] = static_cast<T>(rnd(seed_gen)) & l_mask;
        y[i] = static_cast<T>(y_orig[i]) & r_mask;

        uint8_t tmp = 0;
        func(x[i], y[i], z[i], tmp);
        z_cb[i] = (tmp != 0);

        tmp = 0;
        T x_shift = shift_test ? mask_bits(x[i], right_bits) : x[i];
        func(x_shift, x_shift, w[i], tmp);
        w_cb[i] = (tmp != 0);

        // std::cout << std::hex << "x=" << x[i] << ", y=" << y[i] << "(" <<
        // y_orig[i] << "), z=" << z[i] << ", z_cb=" << (int)z_cb[i] <<
        // ", w[i]=" << w[i] << ", w_cb=" << (int)w_cb[i]  << std::dec <<
        // std::endl;
    }

    std::chrono::system_clock::time_point start, end;
    double elapsed = 0.0;

    for (size_t i = 0; i < num_test; ++i) {
        std::cout << '.' << std::flush;
        std::unique_ptr<TFHEpp::SecretKey> sk_uptr(new TFHEpp::SecretKey());
         const TFHEpp::SecretKey* const sk = sk_uptr.get();
        std::unique_ptr<TFHEpp::EvalKey> ek_uptr(new TFHEpp::EvalKey());
        TFHEpp::EvalKey* const ek = ek_uptr.get();
        ek->emplacebkfft<TFHEpp::lvl01param>(*sk);
        ek->emplaceiksk<TFHEpp::lvl10param>(*sk);

        size_t ans_idx;
        TFHEpp::tfhe_uintN_t<left_bits, TP> cx(x[i], sk, ek), cans(0, sk, ek),
            cx2, cy2;
        TFHEpp::tfhe_uintN_t<1, TP> ccb(0, sk, ek);
        Tr cy;
        init_right_val(cy, y_orig[i], right_bits, sk, ek);

        start = std::chrono::system_clock::now();
        func_tfhe(cx, cy, cans, ccb);
        end = std::chrono::system_clock::now();
        elapsed +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        ans_idx = i + 0 * num_test;
        ans[ans_idx] = cans.template decrypt<T>(sk);
        ans_cb[ans_idx] = ccb.template decrypt<uint8_t>(sk) != 0;

        std::cout << '1' << std::flush;
        cx2 = cx;
        func_tfhe(cx2, cy, cx2, ccb);
        ans_idx = i + 1 * num_test;
        ans[ans_idx] = cx2.template decrypt<T>(sk);
        ans_cb[ans_idx] = ccb.template decrypt<uint8_t>(sk) != 0;

        if constexpr (2 < NUM_CASE) {
            TFHEpp::tfhe_uintN_t<left_bits, TP> cx_shift(
                mask_bits(x[i], right_bits), sk, ek);

            std::cout << '2' << std::flush;
            cx2 = shift_test ? cx_shift : cx;
            func_tfhe(cx2, cx2, cans, ccb);
            ans_idx = i + 2 * num_test;
            ans[ans_idx] = cans.template decrypt<T>(sk);
            ans_cb[ans_idx] = ccb.template decrypt<uint8_t>(sk) != 0;

            std::cout << '3' << std::flush;
            cx2 = shift_test ? cx_shift : cx;
            func_tfhe(cx2, cx2, cx2, ccb);
            ans_idx = i + 3 * num_test;
            ans[ans_idx] = cx2.template decrypt<T>(sk);
            ans_cb[ans_idx] = ccb.template decrypt<uint8_t>(sk) != 0;

            std::cout << '4' << std::flush;
            cy2 = cy;
            func_tfhe(cx, cy2, cy2, ccb);
            ans_idx = i + 4 * num_test;
            ans[ans_idx] = cy2.template decrypt<T>(sk);
            ans_cb[ans_idx] = ccb.template decrypt<uint8_t>(sk) != 0;
        }
    }

    std::cout << std::endl
              << (elapsed / num_test) << "ms (" << elapsed << "ms)"
              << std::endl;

    size_t n_mismatch_ans = 0, n_mismatch_cb = 0;
    for (size_t i = 0; i < NUM_CASE; ++i) {
        const bool ans_is_z = (i != 2 && i != 3);
        const size_t off = i * num_test;
        const size_t n_mans = check_mismatch(ans, ans_is_z ? z : w, x,
                                             ans_is_z ? y : x, num_test, off);
        n_mismatch_ans += n_mans;
        if (n_mans != 0) {
            std::cout << std::dec << n_mismatch_ans << " mismatch(es) found ("
                      << i << ")" << std::endl;
        }
        const size_t n_mcb = check_mismatch(ans_cb, ans_is_z ? z_cb : w_cb, x,
                                            ans_is_z ? y : x, num_test, off);
        n_mismatch_cb += n_mcb;
        if (n_mcb != 0) {
            std::cout << std::dec << n_mismatch_cb
                      << " mismatch(es) found (carry/borrow, " << i << ")"
                      << std::endl;
        }
    }

    if (n_mismatch_ans == 0 && n_mismatch_cb == 0)
        std::cout << "Passed" << std::endl;
}

template <typename T, typename Tr, size_t in_bits, typename Tout,
          size_t out_bits, typename Trnd>
void test_mul(const size_t num_test)
{
    constexpr bool can_copy_x_y =
        std::is_same_v<TFHEpp::tfhe_uintN_t<in_bits, TP>, Tr> &&
        (in_bits == out_bits);
    constexpr size_t NUM_CASE = can_copy_x_y ? 4 : 1;

    std::random_device seed_gen;
    std::uniform_int_distribution<Trnd> rnd(0,
                                            std::numeric_limits<Trnd>::max());

    std::vector<Trnd> y_orig(num_test);
    std::vector<T> x(num_test), y(num_test);
    std::vector<Tout> z(num_test), w(num_test), ans(NUM_CASE * num_test);
    for (size_t i = 0; i < num_test; ++i) {
        y_orig[i] = rnd(seed_gen);
        x[i] = mask_bits(rnd(seed_gen), in_bits),
        y[i] = mask_bits(y_orig[i], in_bits),
        z[i] = mask_bits(static_cast<Tout>(static_cast<uint64_t>(x[i]) *
                                           static_cast<uint64_t>(y[i])),
                         out_bits);
        w[i] = mask_bits(static_cast<Tout>(static_cast<uint64_t>(x[i]) *
                                           static_cast<uint64_t>(x[i])),
                         out_bits);
        // std::cout << std::hex << "x=" << x[i] << ", y=" << y[i] << ", z=" <<
        // z[i] << ", w[i]=" << w[i] << std::dec << std::endl;
    }

    std::chrono::system_clock::time_point start, end;
    double elapsed = 0.0;

    for (size_t i = 0; i < num_test; ++i) {
        std::cout << '.' << std::flush;
        std::unique_ptr<TFHEpp::SecretKey> sk_uptr(new TFHEpp::SecretKey());
         const TFHEpp::SecretKey* const sk = sk_uptr.get();
        std::unique_ptr<TFHEpp::EvalKey> ek_uptr(new TFHEpp::EvalKey());
        TFHEpp::EvalKey* const ek = ek_uptr.get();
        ek->emplacebkfft<TFHEpp::lvl01param>(*sk);
        ek->emplaceiksk<TFHEpp::lvl10param>(*sk);

        TFHEpp::tfhe_uintN_t<in_bits, TP> cx(x[i], sk, ek), cx2, cy2;
        TFHEpp::tfhe_uintN_t<out_bits, TP> cans(0, sk, ek);
        Tr cy;
        init_right_val(cy, y[i], in_bits, sk, ek);

        start = std::chrono::system_clock::now();
        cx.mul(cy, cans);
        end = std::chrono::system_clock::now();
        elapsed +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        ans[i + 0 * num_test] = cans.template decrypt<Tout>(sk);

        if constexpr (can_copy_x_y) {
            std::cout << '1' << std::flush;
            cx2 = cx;
            cx2.mul(cy, cx2);
            ans[i + 1 * num_test] = cx2.template decrypt<Tout>(sk);

            std::cout << '2' << std::flush;
            cy2 = cy;
            cx.mul(cy2, cy2);
            ans[i + 2 * num_test] = cy2.template decrypt<Tout>(sk);

            std::cout << '3' << std::flush;
            cx2 = cx;
            cx2.mul(cx2, cx2);
            ans[i + 3 * num_test] = cx2.template decrypt<Tout>(sk);
        }
    }

    std::cout << std::endl
              << (elapsed / num_test) << "ms (" << elapsed << "ms)"
              << std::endl;

    size_t n_mismatch = 0;
    for (size_t i = 0; i < NUM_CASE; ++i) {
        const bool ans_is_z = (i != 3);
        const size_t n_m = check_mismatch(
            ans, ans_is_z ? z : w, x, ans_is_z ? y : x, num_test, i * num_test);
        if (n_m != 0) {
            std::cout << std::dec << n_mismatch << " mismatch(es) found (" << i
                      << ")" << std::endl;
        }
        n_mismatch += n_m;
    }

    if (n_mismatch == 0) std::cout << "Passed" << std::endl;
}

template <typename T, typename Tc, typename Tr, size_t divisor_bits,
          typename Trnd, typename F, typename Fc>
void test_div(const size_t num_test, F func, Fc func_tfhe)
{
    constexpr bool Tr_is_int = std::is_integral_v<Tr>;
    constexpr size_t NUM_CASE = Tr_is_int ? 3 : 9;

    std::random_device seed_gen;
    std::uniform_int_distribution<Trnd> rnd(0,
                                            std::numeric_limits<Trnd>::max());

    std::vector<Trnd> y_orig(num_test);
    std::vector<T> x(num_test), y(num_test);
    std::vector<T> q(num_test), r(num_test);
    std::vector<T> s(num_test), t(num_test);
    std::vector<T> ans_q(NUM_CASE * num_test), ans_r(NUM_CASE * num_test);
    for (size_t i = 0; i < num_test; ++i) {
        y_orig[i] = rnd(seed_gen);
        x[i] = mask_bits(rnd(seed_gen), std::numeric_limits<T>::digits);
        y[i] = mask_bits(y_orig[i], divisor_bits);
        func(x[i], y[i], q[i], r[i]);
        func(x[i], x[i], s[i], t[i]);
        // std::cout << std::hex << "x=" << x[i] << ", y=" << y[i] << ", q=" <<
        // q[i] << ", r=" << r[i] << ", s=" << s[i] << ", t=" << t[i] <<
        // std::dec << std::endl;
    }

    std::chrono::system_clock::time_point start, end;
    double elapsed = 0.0;

    for (size_t i = 0; i < num_test; ++i) {
        std::cout << '.' << std::flush;
        std::unique_ptr<TFHEpp::SecretKey> sk_uptr(new TFHEpp::SecretKey());
        const TFHEpp::SecretKey* const sk = sk_uptr.get();
        std::unique_ptr<TFHEpp::EvalKey> ek_uptr(new TFHEpp::EvalKey());
        TFHEpp::EvalKey* const ek = ek_uptr.get();
        ek->emplacebkfft<TFHEpp::lvl01param>(*sk);
        ek->emplaceiksk<TFHEpp::lvl10param>(*sk);

        Tc cx(x[i], sk, ek), cq(0, sk, ek), cr(0, sk, ek), cx2, cy2;
        Tr cy;
        init_right_val(cy, y[i], divisor_bits, sk, ek);

        start = std::chrono::system_clock::now();
        func_tfhe(cx, cy, cq, cr);
        end = std::chrono::system_clock::now();
        elapsed +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        ans_q[i + 0 * num_test] = cq.template decrypt<T>(sk);
        ans_r[i + 0 * num_test] = cr.template decrypt<T>(sk);

        std::cout << '1' << std::flush;
        cx2 = cx;
        func_tfhe(cx2, cy, cx2, cr);
        ans_q[i + 1 * num_test] = cx2.template decrypt<T>(sk);
        ans_r[i + 1 * num_test] = cr.template decrypt<T>(sk);

        std::cout << '2' << std::flush;
        cx2 = cx;
        func_tfhe(cx2, cy, cq, cx2);
        ans_q[i + 2 * num_test] = cq.template decrypt<T>(sk);
        ans_r[i + 2 * num_test] = cx2.template decrypt<T>(sk);

        if constexpr (!Tr_is_int) {
            std::cout << '3' << std::flush;
            cy2 = cy;
            func_tfhe(cx, cy2, cy2, cr);
            ans_q[i + 3 * num_test] = cy2.template decrypt<T>(sk);
            ans_r[i + 3 * num_test] = cr.template decrypt<T>(sk);

            std::cout << '4' << std::flush;
            cy2 = cy;
            func_tfhe(cx, cy2, cq, cy2);
            ans_q[i + 4 * num_test] = cq.template decrypt<T>(sk);
            ans_r[i + 4 * num_test] = cy2.template decrypt<T>(sk);

            std::cout << '5' << std::flush;
            cx2 = cx;
            cy2 = cy;
            func_tfhe(cx2, cy2, cx2, cy2);
            ans_q[i + 5 * num_test] = cx2.template decrypt<T>(sk);
            ans_r[i + 5 * num_test] = cy2.template decrypt<T>(sk);

            std::cout << '6' << std::flush;
            cx2 = cx;
            cy2 = cy;
            func_tfhe(cx2, cy2, cy2, cx2);
            ans_q[i + 6 * num_test] = cy2.template decrypt<T>(sk);
            ans_r[i + 6 * num_test] = cx2.template decrypt<T>(sk);

            std::cout << '7' << std::flush;
            cx2 = cx;
            func_tfhe(cx2, cx2, cx2, cr);
            ans_q[i + 7 * num_test] = cx2.template decrypt<T>(sk);
            ans_r[i + 7 * num_test] = cr.template decrypt<T>(sk);

            std::cout << '8' << std::flush;
            cx2 = cx;
            func_tfhe(cx2, cx2, cq, cx2);
            ans_q[i + 8 * num_test] = cq.template decrypt<T>(sk);
            ans_r[i + 8 * num_test] = cx2.template decrypt<T>(sk);
        }
    }

    std::cout << std::endl
              << (elapsed / num_test) << "ms (" << elapsed << "ms)"
              << std::endl;

    size_t n_mismatch_q = 0, n_mismatch_r = 0;
    for (size_t i = 0; i < NUM_CASE; ++i) {
        const bool ans_is_qr = (i != 7) && (i != 8);
        const size_t n_mq =
            check_mismatch(ans_q, ans_is_qr ? q : s, x, ans_is_qr ? y : x,
                           num_test, i * num_test);
        n_mismatch_q += n_mq;
        if (n_mq != 0) {
            std::cout << std::dec << n_mq << " mismatch(es) found (q, " << i
                      << ")" << std::endl;
        }
        const size_t n_mr =
            check_mismatch(ans_r, ans_is_qr ? r : t, x, ans_is_qr ? y : x,
                           num_test, i * num_test);
        n_mismatch_r += n_mr;
        if (n_mr != 0) {
            std::cout << std::dec << n_mr << " mismatch(es) found (r, " << i
                      << ")" << std::endl;
        }
    }

    if (n_mismatch_q == 0 && n_mismatch_r == 0)
        std::cout << "Passed" << std::endl;
}

template <typename Tin, typename Tcin, typename Tout, typename Tcout,
          typename F>
void test_unary_operator(const size_t num_test, F func)
{
    std::random_device seed_gen;
    std::uniform_int_distribution<Tin> rnd(0, std::numeric_limits<Tin>::max());

    std::vector<Tin> x(num_test);
    std::vector<Tout> z(num_test), ans(num_test);
    for (size_t i = 0; i < num_test; ++i) {
        x[i] = rnd(seed_gen);
        func(x[i], z[i]);
        // std::cout << std::hex << "x=" << x[i] << ", z=" << z[i] << std::dec
        // << std::endl;
    }

    std::chrono::system_clock::time_point start, end;
    double elapsed = 0.0;

    for (size_t i = 0; i < num_test; ++i) {
        std::cout << '.' << std::flush;
        std::unique_ptr<TFHEpp::SecretKey> sk_uptr(new TFHEpp::SecretKey());
         const TFHEpp::SecretKey* const sk = sk_uptr.get();
        std::unique_ptr<TFHEpp::EvalKey> ek_uptr(new TFHEpp::EvalKey());
        TFHEpp::EvalKey* const ek = ek_uptr.get();
        ek->emplacebkfft<TFHEpp::lvl01param>(*sk);
        ek->emplaceiksk<TFHEpp::lvl10param>(*sk);

        Tcin cx(x[i], sk, ek);
        Tcout cans(0, sk, ek);

        start = std::chrono::system_clock::now();
        func(cx, cans);
        end = std::chrono::system_clock::now();

        elapsed +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        ans[i] = cans.template decrypt<Tout>(sk);
    }

    std::cout << std::endl
              << (elapsed / num_test) << "ms (" << elapsed << "ms)"
              << std::endl;

    const size_t n_mismatch = check_mismatch(ans, z, x, num_test, 0);
    if (n_mismatch != 0)
        std::cout << std::dec << n_mismatch << " mismatch(es) found"
                  << std::endl;

    if (n_mismatch == 0) std::cout << "Passed" << std::endl;
}

template <typename T, size_t left_bits, typename Tr, typename Tcr,
          size_t right_bits, typename Trnd, typename F>
void test_comparison(const size_t num_test, F func)
{
    constexpr size_t NUM_CASE = 2;

    std::random_device seed_gen;
    std::uniform_int_distribution<Trnd> rnd(0,
                                            std::numeric_limits<Trnd>::max());

    const Trnd l_mask = make_mask<Trnd>(left_bits);
    const Trnd r_mask = make_mask<Trnd>(right_bits);
    std::vector<Trnd> y_orig(num_test);
    std::vector<T> x(num_test);
    std::vector<Tr> y(num_test);
    std::vector<bool> z(num_test), w(num_test), ans(num_test * NUM_CASE);
    for (size_t i = 0; i < num_test; ++i) {
        // x and y are equal in num_test/4 cases
        if (i < num_test / 4) {
            x[i] = rnd(seed_gen) & l_mask;
            y[i] = y_orig[i] = x[i];
        }
        else {
            x[i] = rnd(seed_gen) & l_mask;
            y[i] = (y_orig[i] = rnd(seed_gen)) & r_mask;
        }
        func(x[i], y[i], z[i]);
        func(x[i], x[i], w[i]);
        // std::cout << std::hex << "x=" << x[i] << ", y=" << y[i] <<
        // "(orig=" << y_orig[i] << "), z=" << z[i] << ", w=" << w[i] <<
        // std::dec << std::endl;
    }

    std::chrono::system_clock::time_point start, end;
    double elapsed = 0.0;

    for (size_t i = 0; i < num_test; ++i) {
        std::cout << '.' << std::flush;
        std::unique_ptr<TFHEpp::SecretKey> sk_uptr(new TFHEpp::SecretKey());
         const TFHEpp::SecretKey* const sk = sk_uptr.get();
        std::unique_ptr<TFHEpp::EvalKey> ek_uptr(new TFHEpp::EvalKey());
        TFHEpp::EvalKey* const ek = ek_uptr.get();
        ek->emplacebkfft<TFHEpp::lvl01param>(*sk);
        ek->emplaceiksk<TFHEpp::lvl10param>(*sk);

        TFHEpp::tfhe_uintN_t<left_bits, TP> cx(x[i], sk, ek), cx2;
        TFHEpp::tfhe_uintN_t<1, TP> cans(0, sk, ek);
        Tcr cy;
        init_right_val(cy, y_orig[i], right_bits, sk, ek);

        start = std::chrono::system_clock::now();
        func(cx, cy, cans);
        end = std::chrono::system_clock::now();
        elapsed +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count();
        ans[i + 0 * num_test] = (cans.template decrypt<uint8_t>(sk) != 0);

        std::cout << '1' << std::flush;
        cx2 = cx;
        func(cx2, cx2, cans);
        ans[i + 1 * num_test] = (cans.template decrypt<uint8_t>(sk) != 0);
    }

    std::cout << std::endl
              << (elapsed / num_test) << "ms (" << elapsed << "ms)"
              << std::endl;

    size_t n_mismatch = 0;
    for (size_t i = 0; i < NUM_CASE; ++i) {
        const bool ans_is_z = (i == 0);
        const size_t n =
            ans_is_z ? check_mismatch(ans, z, x, y, num_test, i * num_test)
                     : check_mismatch(ans, w, x, x, num_test, i * num_test);
        n_mismatch += n;
        if (n != 0) {
            std::cout << std::dec << n_mismatch << " mismatch(es) found (" << i
                      << ")" << std::endl;
        }
    }

    if (n_mismatch == 0) std::cout << "Passed" << std::endl;
}

template <typename T, typename U, typename V>
void add_func(const T& x, const U& y, T& z, [[maybe_unused]] V& _)
{
    z = x + y;
}
auto add = [](auto&&... args) { add_func(args...); };

template <typename T, typename U, typename V>
void sub_func(const T& x, const U& y, T& z, [[maybe_unused]] V& _)
{
    z = x - y;
}
auto sub = [](auto&&... args) { sub_func(args...); };

template <typename T, typename U, typename V>
void mul_func(const T& x, const U& y, T& z, [[maybe_unused]] V& _)
{
    z = x * y;
}
auto mul = [](auto&&... args) { mul_func(args...); };

template <typename T, typename U>
void div_func(const T& x, const U& y, T& q, T& r)
{
    // &q and &r may be point to x or y.
    T xx, x_save;
    xx = x_save = x;
    U yy, y_save;
    yy = y_save = y;
    q = xx / yy;
    xx = x_save;
    yy = y_save;
    r = xx % yy;
}
auto div_ = [](auto&&... args) { div_func(args...); };

template <typename T, typename U, typename V>
void and_func(const T& x, const U& y, T& z, [[maybe_unused]] V& _)
{
    z = x & y;
}
auto and_ = [](auto&&... args) { and_func(args...); };

template <typename T, typename U, typename V>
void or_func(const T& x, const U& y, T& z, [[maybe_unused]] V& _)
{
    z = x | y;
}
auto or_ = [](auto&&... args) { or_func(args...); };

template <typename T, typename U, typename V>
void xor_func(const T& x, const U& y, T& z, [[maybe_unused]] V& _)
{
    z = x ^ y;
}
auto xor_ = [](auto&&... args) { xor_func(args...); };

template <typename T>
void not_func(const T& x, T& z)
{
    z = ~x;
}
auto not_ = [](auto&&... args) { not_func(args...); };

template <typename T, typename U, typename V>
void l_shift_func(const T& x, const U& y, T& z, [[maybe_unused]] V& _)
{
    z = x << y;
}
auto l_shift = [](auto&&... args) { l_shift_func(args...); };

template <typename T, typename U, typename V>
void r_shift_func(const T& x, const U& y, T& z, [[maybe_unused]] V& _)
{
    z = x >> y;
}
auto r_shift = [](auto&&... args) { r_shift_func(args...); };

template <typename T, typename U, typename V>
void equal_func(const T& x, const U& y, V& z)
{
    z = (x == y);
}
auto equal = [](auto&&... args) { equal_func(args...); };

template <typename T, typename U, typename V>
void not_equal_func(const T& x, const U& y, V& z)
{
    z = (x != y);
}
auto not_equal = [](auto&&... args) { not_equal_func(args...); };

template <typename T, typename U, typename V>
void below_func(const T& x, const U& y, V& z)
{
    z = (x < y);
}
auto below = [](auto&&... args) { below_func(args...); };

template <typename T, typename U, typename V>
void below_or_equal_func(const T& x, const U& y, V& z)
{
    z = (x <= y);
}
auto below_or_equal = [](auto&&... args) { below_or_equal_func(args...); };

template <typename T, typename U, typename V>
void above_func(const T& x, const U& y, V& z)
{
    z = (x > y);
}
auto above = [](auto&&... args) { above_func(args...); };

template <typename T, typename U, typename V>
void above_or_equal_func(const T& x, const U& y, V& z)
{
    z = (x >= y);
}
auto above_or_equal = [](auto&&... args) { above_or_equal_func(args...); };

using u64 = uint64_t;
using t64 = TFHEpp::tfhe_uint64_t<TP>;
using u32 = uint32_t;
using t32 = TFHEpp::tfhe_uint32_t<TP>;
using u16 = uint16_t;
using t16 = TFHEpp::tfhe_uint16_t<TP>;
using t8 = TFHEpp::tfhe_uint8_t<TP>;
using t1 = TFHEpp::tfhe_uintN_t<1, TP>;

#define COMMA ,
#define COMPILE_TIME_FOR(FUNC, SUFFIX)                   \
    template <size_t from, size_t to>                    \
    struct compile_time_for_##SUFFIX {                   \
        void operator()()                                \
        {                                                \
            FUNC;                                        \
            compile_time_for_##SUFFIX<from + 1, to>()(); \
        }                                                \
    };                                                   \
    template <size_t to>                                 \
    struct compile_time_for_##SUFFIX<to, to> {           \
        void operator()() {}                             \
    };
COMPILE_TIME_FOR(
    test_mul<u16 COMMA t16 COMMA 16 COMMA u32 COMMA from COMMA u16>(5),
    test_mul1);
COMPILE_TIME_FOR(
    test_mul<u16 COMMA u64 COMMA 16 COMMA u32 COMMA from COMMA u64>(5),
    test_mul2);
#undef COMPILE_TIME_FOR
#undef COMMA

int main()
{
    std::cout << "------ Test of operator+() ------" << std::endl;
    test<u32, 32, t32, 32, u32>(10, add, add);
    test<u32, 32, u64, 64, u64>(10, add, add);
    test<u64, 64, u64, 64, u64>(10, add, add);

    std::cout << "------ Test of add() ------" << std::endl;
    test<u32, 32, t32, 32, u32>(
        10,
        [](const u32 x, const u32 y, u32& z, uint8_t& c) {
            const u32 d = std::numeric_limits<u32>::max() - y;
            c = d <= x;
            z = x + y;
        },
        [](const t32& x, const t32& y, t32& res, t1& carry) {
            x.add(y, res, carry);
        });
    test<u16, 16, t8, 8, u16>(
        10,
        [](const u16 x, const u16 y, u16& z, uint8_t& c) {
            const u16 d = std::numeric_limits<u16>::max() - y;
            c = d <= x;
            z = x + y;
        },
        [](const t16& x, const t8& y, t16& res, t1& carry) {
            x.add(y, res, carry);
        });
    test<u32, 16, TFHEpp::tfhe_uintN_t<24, TP>, 24, u32>(
        10,
        [](const u32 x, const u32 y, u32& z, uint8_t& c) {
            const u32 mask = make_mask<u32>(16);
            const u32 new_y = y & mask;
            const u32 d = std::numeric_limits<u16>::max() - new_y;
            c = d <= x;
            z = (x + new_y) & mask;
        },
        [](const t16& x, const TFHEpp::tfhe_uintN_t<24, TP>& y, t16& res,
           t1& carry) { x.add(y, res, carry); });
    test<u32, 32, u64, 64, u64>(
        10,
        [](const u32 x, const u32 y, u32& z, uint8_t& c) {
            const u32 d = std::numeric_limits<u16>::max() - y;
            c = d <= x;
            z = x + y;
        },
        [](const t32& x, const u64 y, t32& res, t1& carry) {
            x.add(y, res, carry);
        });
    test<u64, 64, u64, 64, u64>(
        10,
        [](const u64 x, const u64 y, u64& z, uint8_t& c) {
            const u64 d = std::numeric_limits<u64>::max() - y;
            c = d <= x;
            z = x + y;
        },
        [](const t64& x, const u64 y, t64& res, t1& carry) {
            x.add(y, res, carry);
        });

    std::cout << "------ Test of operator-() ------" << std::endl;
    test<u32, 32, t32, 32, u32>(10, sub, sub);
    test<u32, 32, u64, 64, u64>(10, sub, sub);
    test<u64, 64, u64, 64, u64>(10, sub, sub);

    std::cout << "------ Test of sub() ------" << std::endl;
    test<u32, 32, t32, 32, u32>(
        10,
        [](const u32 x, const u32 y, u32& z, uint8_t& b) {
            b = x < y;
            z = x - y;
        },
        [](const t32& x, const t32& y, t32& res, t1& borrow) {
            x.sub(y, res, borrow);
        });
    test<u32, 32, u64, 64, u64>(
        10,
        [](const u32 x, const u32 y, u32& z, uint8_t& b) {
            b = x < y;
            z = x - y;
        },
        [](const t32& x, const u64 y, t32& res, t1& borrow) {
            x.sub(y, res, borrow);
        });
    test<u64, 64, u64, 64, u64>(
        10,
        [](const u64 x, const u64 y, u64& z, uint8_t& b) {
            b = x < y;
            z = x - y;
        },
        [](const t64& x, const u64 y, t64& res, t1& borrow) {
            x.sub(y, res, borrow);
        });

    std::cout << "------ Test of operator*() ------" << std::endl;
    test<u16, 16, t16, 16, u16>(5, mul, mul);
    test<u16, 16, u32, 32, u32>(5, mul, mul);
    test<u16, 16, u16, 16, u16>(5, mul, mul);

    std::cout << "------ Test of mul() ------" << std::endl;
    compile_time_for_test_mul1<1, 32 + 1>()();
    compile_time_for_test_mul2<1, 32 + 1>()();
    // test_mul<u16, t16, 16, u32, 16, u16>(5);
    // test_mul<u16, t16, 16, u32, 32, u16>(5);
    // test_mul<u16, u64, 16, u32, 32, u64>(5);

    std::cout << "------ Test of operator/(), "
                 "operator%() ------"
              << std::endl;
    test_div<u32, t32, t32, 16, u32>(5, div_, div_);
    test_div<u32, t32, u32, 16, u32>(5, div_, div_);

    std::cout << "------ Test of div() ------" << std::endl;
    test_div<u32, t32, t32, 16, u32>(
        5, div_,
        [](const t32& x, const t32& y, t32& q, t32& r) { x.div(y, q, r); });
    test_div<u32, t32, u32, 16, u32>(
        5, div_,
        [](const t32& x, const u32& y, t32& q, t32& r) { x.div(y, q, r); });

    std::cout << "------ Test of operator&() ------" << std::endl;
    test<u32, 32, t32, 32, u32>(10, and_, and_);
    test<u32, 32, u64, 64, u64>(10, and_, and_);
    test<u64, 64, u64, 64, u64>(10, and_, and_);
    test<u32, 32, t1, 1, u32>(
        20,
        [](const u32 x, const u32 y, u32& z, [[maybe_unused]] uint8_t& _) {
            if (y == 0)
                z = 0;
            else
                z = x;
        },
        [](const t32& x, const t1& y, t32& z, [[maybe_unused]] t1& _) { z = x & y; });

    std::cout << "------ Test of operator|() ------" << std::endl;
    test<u32, 32, t32, 32, u32>(10, or_, or_);
    test<u32, 32, u64, 64, u64>(10, or_, or_);
    test<u64, 64, u64, 64, u64>(10, or_, or_);

    std::cout << "------ Test of operator^() ------" << std::endl;
    test<u32, 32, t32, 32, u32>(10, xor_, xor_);
    test<u32, 32, u64, 64, u64>(10, xor_, xor_);
    test<u64, 64, u64, 64, u64>(10, xor_, xor_);

    std::cout << "------ Test of operator<<() ------" << std::endl;
    test<u32, 32, t32, 5, u32>(10, l_shift, l_shift, true);
    test<u32, 32, u64, 5, u64>(10, l_shift, l_shift, true);
    test<u64, 64, u64, 6, u64>(10, l_shift, l_shift, true);

    std::cout << "------ Test of operator>>() ------" << std::endl;
    test<u32, 32, t32, 5, u32>(10, r_shift, r_shift, true);
    test<u32, 32, u64, 5, u64>(10, r_shift, r_shift, true);
    test<u64, 64, u64, 6, u64>(10, r_shift, r_shift, true);

    std::cout << "------ Test of operator~() ------" << std::endl;
    test_unary_operator<u32, t32, u32, t32>(20, not_);

    std::cout << "------ Test of operator==() ------" << std::endl;
    test_comparison<u32, 32, u32, t32, 32, u32>(10, equal);
    test_comparison<u32, 32, u64, u64, 64, u64>(10, equal);
    test_comparison<u64, 64, u64, u64, 64, u64>(10, equal);

    std::cout << "------ Test of operator!=() ------" << std::endl;
    test_comparison<u32, 32, u32, t32, 32, u32>(10, not_equal);
    test_comparison<u32, 32, u64, u64, 64, u64>(10, not_equal);
    test_comparison<u64, 64, u64, u64, 64, u64>(10, not_equal);

    std::cout << "------ Test of operator<() ------" << std::endl;
    test_comparison<u32, 32, u32, t32, 32, u32>(10, below);
    test_comparison<u32, 32, u64, u64, 64, u64>(10, below);
    test_comparison<u64, 64, u64, u64, 64, u64>(10, below);

    std::cout << "------ Test of operator>=() ------" << std::endl;
    test_comparison<u32, 32, u32, t32, 32, u32>(10, above_or_equal);
    test_comparison<u32, 32, u64, u64, 64, u64>(10, above_or_equal);
    test_comparison<u64, 64, u64, u64, 64, u64>(10, above_or_equal);

    std::cout << "------ Test of operator>() ------" << std::endl;
    test_comparison<u32, 32, u32, t32, 32, u32>(10, above);
    test_comparison<u32, 32, u64, u64, 64, u64>(10, above);
    test_comparison<u64, 64, u64, u64, 64, u64>(10, above);

    std::cout << "------ Test of operator<=() ------" << std::endl;
    test_comparison<u32, 32, u32, t32, 32, u32>(10, below_or_equal);
    test_comparison<u32, 32, u64, u64, 64, u64>(10, below_or_equal);
    test_comparison<u64, 64, u64, u64, 64, u64>(10, below_or_equal);

    std::cout << "Finished." << std::endl;

    return 0;
}
