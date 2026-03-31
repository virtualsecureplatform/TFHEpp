// simd_ops.cpp — Test SIMD slot-wise addition and multiplication
//
// Encrypts two slot vectors, performs:
//   (a) coefficient-wise TRLWE addition  → expects element-wise add mod t
//   (b) TRLWEMultFullDD                  → expects element-wise mul mod t
//
// Uses lvl3simdparam (n=4096, t=114689, DD l̅=8) so multiplication uses the
// Double Decomposition path automatically (P::l̅ > 1).

#include <cassert>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{


    using P = TFHEpp::lvl3simdparam;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);

    std::cout << "SIMD slot operations test" << std::endl;
    std::cout << "  n=" << P::n << "  t=" << t
              << "  l̅=" << P::l̅ << "  B̅gbit=" << P::B̅gbit << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, t - 1);
    std::uniform_int_distribution<int> binary(0, 1);

    // Generate secret key
    TFHEpp::Key<P> key;
    for (auto &v : key)
        v = binary(engine) ? static_cast<typename P::T>(1)
                           : static_cast<typename P::T>(-1);

    // Generate relinearization key on heap (avoids ~6MB stack allocation in
    // relinKeyFFTgen for lvl3simdparam — uses makeRelinKeyFFT from bfv-slots.hpp)
    auto relinkeyfft = TFHEpp::makeRelinKeyFFT<P>(key);

    // ------------------------------------------------------------------
    // Test 1: SIMD addition
    // ------------------------------------------------------------------
    std::cout << "  Testing SIMD addition... " << std::flush;
    {
        constexpr int num_test = 3;
        int failures = 0;
        for (int test = 0; test < num_test; test++) {
            std::array<uint64_t, P::n> slots_a, slots_b, slots_sum, expected;
            for (auto &v : slots_a) v = slot_dist(engine);
            for (auto &v : slots_b) v = slot_dist(engine);
            for (int i = 0; i < static_cast<int>(P::n); i++)
                expected[i] = (slots_a[i] + slots_b[i]) % t;

            TFHEpp::TRLWE<P> ct_a, ct_b, ct_sum;
            TFHEpp::trlweSlotEncrypt<P>(ct_a, slots_a, key);
            TFHEpp::trlweSlotEncrypt<P>(ct_b, slots_b, key);

            // Addition: coefficient-wise
            for (int k = 0; k <= static_cast<int>(P::k); k++)
                for (int i = 0; i < static_cast<int>(P::n); i++)
                    ct_sum[k][i] = ct_a[k][i] + ct_b[k][i];

            TFHEpp::trlweSlotDecrypt<P>(slots_sum, ct_sum, key);

            for (int i = 0; i < static_cast<int>(P::n); i++) {
                if (slots_sum[i] != expected[i]) {
                    std::cerr << "\n    FAIL test=" << test << " slot=" << i
                              << " expected=" << expected[i]
                              << " got=" << slots_sum[i] << std::endl;
                    failures++;
                    break;
                }
            }
        }
        if (failures == 0)
            std::cout << "PASS" << std::endl;
        else {
            std::cout << "FAIL (" << failures << " failures)" << std::endl;
            return 1;
        }
    }

    // ------------------------------------------------------------------
    // Test 2: SIMD multiplication (TRLWEMultFullDD)
    // ------------------------------------------------------------------
    std::cout << "  Testing SIMD multiplication (TRLWEMultFullDD)... " << std::flush;
    {
        constexpr int num_test = 2;
        int failures = 0;
        for (int test = 0; test < num_test; test++) {
            std::array<uint64_t, P::n> slots_a, slots_b, slots_mul, expected;
            // Use small values to keep noise budget comfortable
            std::uniform_int_distribution<uint64_t> small_dist(0, 100);
            for (auto &v : slots_a) v = small_dist(engine);
            for (auto &v : slots_b) v = small_dist(engine);
            for (int i = 0; i < static_cast<int>(P::n); i++)
                expected[i] = (slots_a[i] * slots_b[i]) % t;

            TFHEpp::TRLWE<P> ct_a, ct_b, ct_mul;
            TFHEpp::trlweSlotEncrypt<P>(ct_a, slots_a, key);
            TFHEpp::trlweSlotEncrypt<P>(ct_b, slots_b, key);

            TFHEpp::TRLWEMultFullDD<P>(ct_mul, ct_a, ct_b, *relinkeyfft);

            TFHEpp::trlweSlotDecrypt<P>(slots_mul, ct_mul, key);

            for (int i = 0; i < static_cast<int>(P::n); i++) {
                if (slots_mul[i] != expected[i]) {
                    std::cerr << "\n    FAIL test=" << test << " slot=" << i
                              << " expected=" << expected[i]
                              << " got=" << slots_mul[i] << std::endl;
                    failures++;
                    break;
                }
            }
        }
        if (failures == 0)
            std::cout << "PASS" << std::endl;
        else {
            std::cout << "FAIL (" << failures << " failures)" << std::endl;
            return 1;
        }
    }

    std::cout << "PASS (all SIMD ops tests)" << std::endl;
    return 0;
}
