// bfv_prime_power_slots.cpp -- checks the p^2 SIMD parameter adapter used by
// the first BFV bootstrapping stages.

#include <cstdint>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl3simdparam;
    using BP = TFHEpp::bfvboot::PrimePower2Param<P>;
    using BK = TFHEpp::bfvboot::BootstrapKey<P>;
    using T = typename BP::T;

    static_assert(P::bfv_bootstrap_digit_error_bound == 15);
    static_assert(P::bfv_bootstrap_linear_bsgs_step == 64);
    static_assert(BP::bfv_bootstrap_digit_error_bound ==
                  P::bfv_bootstrap_digit_error_bound);
    static_assert(BK::digit_error_bound == P::bfv_bootstrap_digit_error_bound);
    static_assert(BK::default_linear_bsgs_step ==
                  P::bfv_bootstrap_linear_bsgs_step);

    [[maybe_unused]] BK (*make_bk)(const TFHEpp::Key<P> &, bool) =
        &TFHEpp::bfvboot::MakeBootstrapKey<P>;
    [[maybe_unused]] void (*bootstrap_noisy_slots)(TFHEpp::TRLWE<BP> &,
                                                   const TFHEpp::TRLWE<P> &,
                                                   const BK &) =
        &TFHEpp::bfvboot::BootstrapNoisySlots<P>;

    constexpr uint64_t p = BP::base_plain_modulus;
    constexpr uint64_t q = static_cast<uint64_t>(BP::plain_modulus);
    constexpr int n = static_cast<int>(BP::n);
    constexpr uint64_t B = P::bfv_bootstrap_digit_error_bound;
    static_assert(2 * B + 1 <= p);

    std::cout << "BFV p^2 SIMD parameter test" << std::endl;
    std::cout << "  n=" << n << "  p=" << p << "  p^2=" << q << std::endl;

    const uint64_t psi_order = TFHEpp::powmod64(BP::simd_psi, 2 * BP::n, q);
    const uint64_t psi_half = TFHEpp::powmod64(BP::simd_psi, BP::n, q);
    if (psi_order != 1 || psi_half != q - 1) {
        std::cerr << "  FAIL root lift: psi^(2n)=" << psi_order
                  << " psi^n=" << psi_half << std::endl;
        return 1;
    }

    const auto removal =
        TFHEpp::digitext::GetLowestDigitRemovalPolynomialOverRange(p, B);
    const std::array<uint64_t, 3> sample_messages{0, 1, p - 1};
    for (uint64_t m : sample_messages) {
        for (int64_t e = -static_cast<int64_t>(B);
             e <= static_cast<int64_t>(B); e++) {
            const uint64_t centered_e =
                e < 0 ? q - static_cast<uint64_t>(-e)
                      : static_cast<uint64_t>(e);
            const uint64_t x =
                (static_cast<unsigned __int128>(p) * m +
                 centered_e) %
                q;
            const uint64_t got =
                TFHEpp::digitext::plainEvalMod(removal, x, q);
            const uint64_t want =
                static_cast<uint64_t>(
                    (static_cast<unsigned __int128>(p) * m) % q);
            if (got != want) {
                std::cerr << "  FAIL default bounded removal: m=" << m
                          << " e=" << e << " expected=" << want
                          << " got=" << got << std::endl;
                return 1;
            }
        }
    }

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, q - 1);

    for (int test = 0; test < 5; test++) {
        std::array<uint64_t, BP::n> slots_in{}, slots_out{};
        TFHEpp::Polynomial<BP> poly{};

        for (auto &v : slots_in) v = slot_dist(engine);
        TFHEpp::SlotEncode<BP>(poly, slots_in);
        TFHEpp::SlotDecode<BP>(slots_out, poly);

        int mismatches = 0;
        int first_bad = -1;
        for (int i = 0; i < n; i++) {
            if (slots_out[i] != slots_in[i]) {
                if (first_bad < 0) first_bad = i;
                mismatches++;
            }
        }
        if (mismatches != 0) {
            std::cerr << "  FAIL encode/decode test=" << test
                      << " mismatches=" << mismatches
                      << " first_bad=" << first_bad
                      << " expected=" << slots_in[first_bad]
                      << " got=" << slots_out[first_bad] << std::endl;
            return 1;
        }
        std::cout << "." << std::flush;
    }
    std::cout << " PASS" << std::endl;

    std::uniform_int_distribution<int> binary(0, 1);
    TFHEpp::Key<BP> key{};
    for (auto &v : key)
        v = binary(engine) ? static_cast<T>(1) : static_cast<T>(-1);

    std::array<uint64_t, BP::n> poly_in{}, poly_out{};
    for (auto &v : poly_in) v = slot_dist(engine);

    TFHEpp::TRLWE<BP> ct{};
    TFHEpp::bfvboot::BfvPolyEncrypt<BP>(ct, poly_in, key);
    TFHEpp::bfvboot::BfvPolyDecrypt<BP>(poly_out, ct, key);

    int mismatches = 0;
    int first_bad = -1;
    for (int i = 0; i < n; i++) {
        if (poly_out[i] != poly_in[i]) {
            if (first_bad < 0) first_bad = i;
            mismatches++;
        }
    }
    if (mismatches != 0) {
        std::cerr << "  FAIL p^2 poly encrypt/decrypt mismatches="
                  << mismatches << " first_bad=" << first_bad
                  << " expected=" << poly_in[first_bad]
                  << " got=" << poly_out[first_bad] << std::endl;
        return 1;
    }

    TFHEpp::TRLWE<P> zero_in{};
    TFHEpp::TRLWE<BP> zero_enc_sk{}, zero_out{};
    TFHEpp::bfvboot::NoisyDecrypt<P, BP>(zero_out, zero_in, zero_enc_sk);
    for (int c = 0; c <= static_cast<int>(BP::k); c++) {
        for (int i = 0; i < n; i++) {
            if (zero_out[c][i] != 0) {
                std::cerr << "  FAIL zero NoisyDecrypt at component=" << c
                          << " index=" << i << std::endl;
                return 1;
            }
        }
    }

    std::cout << "PASS (BFV p^2 SIMD parameter test)" << std::endl;
    return 0;
}
