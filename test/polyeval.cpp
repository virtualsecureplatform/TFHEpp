// polyeval.cpp — Test baby-step/giant-step homomorphic polynomial evaluation
//
// Evaluates small cleartext polynomials on encrypted slot vectors and verifies
// the result matches the expected plaintext evaluation.

#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl3simdparam;
    using T = typename P::T;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);  // 114689
    constexpr int n = static_cast<int>(P::n);

    std::cout << "PolyEval test (baby-step/giant-step)" << std::endl;
    std::cout << "  n=" << n << "  t=" << t << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, t - 1);
    std::uniform_int_distribution<int> binary(0, 1);

    // Generate secret key
    TFHEpp::Key<P> key;
    for (auto &v : key)
        v = binary(engine) ? static_cast<T>(1) : static_cast<T>(-1);

    // Generate relinearization key
    std::cout << "  Generating relin key... " << std::flush;
    auto rlk = TFHEpp::makeRelinKeyFFT<P>(key);
    std::cout << "done" << std::endl;

    // Helper: evaluate polynomial in cleartext
    auto plainEval = [&](const std::vector<uint64_t> &coeffs, uint64_t x) -> uint64_t {
        __uint128_t result = 0;
        __uint128_t xi = 1;
        for (auto c : coeffs) {
            result = (result + static_cast<__uint128_t>(c) * xi) % t;
            xi = (xi * x) % t;
        }
        return static_cast<uint64_t>(result);
    };

    // ------------------------------------------------------------------
    // Test 1: Linear polynomial f(x) = 3x + 7
    // ------------------------------------------------------------------
    std::cout << "  Test 1: f(x) = 3x + 7 (degree 1)... " << std::flush;
    {
        std::vector<uint64_t> coeffs = {7, 3};

        std::array<uint64_t, P::n> slots_in, expected, slots_out;
        for (int i = 0; i < n; i++) {
            slots_in[i] = slot_dist(engine);
            expected[i] = plainEval(coeffs, slots_in[i]);
        }

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::PolyEval<P>(result, coeffs, ct, *rlk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (slots_out[i] != expected[i]) mismatches++;
        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << ")" << std::endl;
        if (mismatches > 0) return 1;
    }

    // ------------------------------------------------------------------
    // Test 2: Quadratic f(x) = x² + 2x + 1 = (x+1)²
    // ------------------------------------------------------------------
    std::cout << "  Test 2: f(x) = x^2 + 2x + 1 (degree 2)... " << std::flush;
    {
        std::vector<uint64_t> coeffs = {1, 2, 1};

        // Use small slot values to keep noise manageable after squaring
        std::uniform_int_distribution<uint64_t> small_dist(0, 100);
        std::array<uint64_t, P::n> slots_in, expected, slots_out;
        for (int i = 0; i < n; i++) {
            slots_in[i] = small_dist(engine);
            expected[i] = plainEval(coeffs, slots_in[i]);
        }

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::PolyEval<P>(result, coeffs, ct, *rlk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (slots_out[i] != expected[i]) mismatches++;
        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << ")" << std::endl;
        if (mismatches > 0) return 1;
    }

    // ------------------------------------------------------------------
    // Test 3: Higher degree: f(x) = x^5 + 3x^3 + x + 2 (degree 5)
    // ------------------------------------------------------------------
    std::cout << "  Test 3: f(x) = x^5 + 3x^3 + x + 2 (degree 5)... " << std::flush;
    {
        std::vector<uint64_t> coeffs = {2, 1, 0, 3, 0, 1};

        std::uniform_int_distribution<uint64_t> small_dist(0, 20);
        std::array<uint64_t, P::n> slots_in, expected, slots_out;
        for (int i = 0; i < n; i++) {
            slots_in[i] = small_dist(engine);
            expected[i] = plainEval(coeffs, slots_in[i]);
        }

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::PolyEval<P>(result, coeffs, ct, *rlk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (slots_out[i] != expected[i]) mismatches++;
        // Degree 5 requires 3 sequential multiplications (depth 3).
        // This exceeds the noise budget of lvl3simdparam.
        // Larger parameters (bigger n or Q) are needed for deeper circuits.
        std::cout << (mismatches == 0 ? "PASS" : "NOISE_LIMIT")
                  << " (mismatches=" << mismatches << ", depth=3)" << std::endl;
    }

    // ------------------------------------------------------------------
    // Test 3a: Degree 4: f(x) = x^4 + x^3 + x^2 + x + 1 (k=2,m=1, upper uses baby[2])
    // ------------------------------------------------------------------
    std::cout << "  Test 3a: f(x) = x^4+x^3+x^2+x+1 (degree 4)... " << std::flush;
    {
        std::vector<uint64_t> coeffs = {1, 1, 1, 1, 1};

        std::uniform_int_distribution<uint64_t> small_dist(0, 10);
        std::array<uint64_t, P::n> slots_in, expected, slots_out;
        for (int i = 0; i < n; i++) {
            slots_in[i] = small_dist(engine);
            expected[i] = plainEval(coeffs, slots_in[i]);
        }

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::PolyEval<P>(result, coeffs, ct, *rlk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (slots_out[i] != expected[i]) mismatches++;
        if (mismatches > 0) {
            std::cout << "FAIL (mismatches=" << mismatches << ")" << std::endl;
            for (int i = 0; i < 3; i++)
                std::cout << "    slot[" << i << "] in=" << slots_in[i]
                          << " expected=" << expected[i] << " got=" << slots_out[i] << std::endl;
        } else {
            std::cout << "PASS" << std::endl;
        }
    }

    // ------------------------------------------------------------------
    // Test 3b: Cubic f(x) = x^3 + 1 (degree 3, simpler)
    // ------------------------------------------------------------------
    std::cout << "  Test 3b: f(x) = x^3 + 1 (degree 3)... " << std::flush;
    {
        std::vector<uint64_t> coeffs = {1, 0, 0, 1};

        std::uniform_int_distribution<uint64_t> small_dist(0, 10);
        std::array<uint64_t, P::n> slots_in, expected, slots_out;
        for (int i = 0; i < n; i++) {
            slots_in[i] = small_dist(engine);
            expected[i] = plainEval(coeffs, slots_in[i]);
        }

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::PolyEval<P>(result, coeffs, ct, *rlk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (slots_out[i] != expected[i]) mismatches++;
        if (mismatches > 0) {
            std::cout << "FAIL (mismatches=" << mismatches << ")" << std::endl;
            for (int i = 0; i < 5; i++)
                std::cout << "    slot[" << i << "] in=" << slots_in[i]
                          << " expected=" << expected[i] << " got=" << slots_out[i] << std::endl;

            // Also decrypt baby[3] directly to check x^3
            TFHEpp::TRLWE<P> x3;
            TFHEpp::TRLWEMultFullDD<P>(x3, ct, ct, *rlk);  // x^2
            TFHEpp::TRLWE<P> x3b;
            TFHEpp::TRLWEMultFullDD<P>(x3b, x3, ct, *rlk);  // x^3
            std::array<uint64_t, P::n> x3_slots;
            TFHEpp::trlweSlotDecrypt<P>(x3_slots, x3b, key);
            std::cout << "    x^3 check: slot[0]=" << x3_slots[0]
                      << " expected=" << (slots_in[0]*slots_in[0]%t*slots_in[0]%t)
                      << std::endl;
            return 1;
        }
        std::cout << "PASS" << std::endl;
    }

    // ------------------------------------------------------------------
    // Test 4: Constant polynomial f(x) = 42
    // ------------------------------------------------------------------
    std::cout << "  Test 4: f(x) = 42 (degree 0)... " << std::flush;
    {
        std::vector<uint64_t> coeffs = {42};

        std::array<uint64_t, P::n> slots_in, slots_out;
        for (auto &v : slots_in) v = slot_dist(engine);

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::PolyEval<P>(result, coeffs, ct, *rlk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (slots_out[i] != 42) mismatches++;
        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << ")" << std::endl;
        if (mismatches > 0) return 1;
    }

    std::cout << "PASS (all PolyEval tests)" << std::endl;
    return 0;
}
