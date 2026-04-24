// linear_transform.cpp — Test SlotPtxtMul and LinearTransform primitives
//
// Verifies:
//   1. SlotPtxtMul: pointwise slot multiplication with a plaintext vector
//      equals what we'd get by multiplying slot-by-slot.
//   2. LinearTransform: Σ_i d_i ∘ rotate(ct, r_i) gives the expected
//      slot-wise linear combination under the Galois intra-row rotation.

#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl3simdparam;
    using T = typename P::T;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int n = static_cast<int>(P::n);
    constexpr int half = n / 2;

    std::cout << "LinearTransform primitives test" << std::endl;
    std::cout << "  n=" << n << "  t=" << t << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, t - 1);
    std::uniform_int_distribution<uint64_t> small_dist(0, 100);
    std::uniform_int_distribution<int> binary(0, 1);

    TFHEpp::Key<P> key;
    for (auto &v : key)
        v = binary(engine) ? static_cast<T>(1) : static_cast<T>(-1);

    std::cout << "  Generating GaloisKey (" << (P::nbit + 1) << " keys)... " << std::flush;
    auto gk = std::make_unique<TFHEpp::GaloisKey<P>>();
    TFHEpp::GaloisKeyGen<P>(*gk, key);
    std::cout << "done" << std::endl;

    int failures = 0;

    // -----------------------------------------------------------------------
    // Test 1: SlotPtxtMul — pointwise multiplication
    // -----------------------------------------------------------------------
    std::cout << "  Test 1: SlotPtxtMul (pointwise)... " << std::flush;
    {
        std::array<uint64_t, P::n> slots_in{}, ptxt_vec{}, slots_out{}, expected{};
        // Keep values small so s_i · v_i stays within t.
        for (int i = 0; i < n; i++) {
            slots_in[i] = small_dist(engine);
            ptxt_vec[i] = small_dist(engine);
            expected[i] = (slots_in[i] * ptxt_vec[i]) % t;
        }

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::SlotPtxtMul<P>(result, ct, ptxt_vec);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (slots_out[i] != expected[i]) mismatches++;

        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;
        if (mismatches > 0) {
            for (int i = 0; i < 3; i++)
                std::cout << "    slot[" << i << "] s=" << slots_in[i]
                          << " v=" << ptxt_vec[i]
                          << " expected=" << expected[i]
                          << " got=" << slots_out[i] << std::endl;
            failures++;
        }
    }

    // -----------------------------------------------------------------------
    // Test 2: LinearTransform with a single identity diagonal
    //         (d_0 = all-ones, r_0 = 0).  Should be identity.
    // -----------------------------------------------------------------------
    std::cout << "  Test 2: LinearTransform identity... " << std::flush;
    {
        std::array<uint64_t, P::n> slots_in{}, slots_out{};
        for (auto &v : slots_in) v = slot_dist(engine);

        std::array<uint64_t, P::n> ones{};
        for (auto &v : ones) v = 1;

        std::vector<std::array<uint64_t, P::n>> diagonals = {ones};
        std::vector<int> offsets = {0};

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::LinearTransform<P>(result, ct, diagonals, offsets, *gk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (slots_out[i] != slots_in[i]) mismatches++;

        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;
        if (mismatches > 0) failures++;
    }

    // -----------------------------------------------------------------------
    // Test 3: LinearTransform as a rotation
    //         d_0 = all-ones, r_0 = 1.  Should rotate each row by 1.
    // -----------------------------------------------------------------------
    std::cout << "  Test 3: LinearTransform single rotation... " << std::flush;
    {
        std::array<uint64_t, P::n> slots_in{}, slots_out{};
        for (auto &v : slots_in) v = slot_dist(engine);

        std::array<uint64_t, P::n> ones{};
        for (auto &v : ones) v = 1;

        std::vector<std::array<uint64_t, P::n>> diagonals = {ones};
        std::vector<int> offsets = {1};

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::LinearTransform<P>(result, ct, diagonals, offsets, *gk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++) {
            uint64_t expected = (i < half)
                ? slots_in[(i + 1) % half]
                : slots_in[half + (i - half + 1) % half];
            if (slots_out[i] != expected) mismatches++;
        }

        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;
        if (mismatches > 0) failures++;
    }

    // -----------------------------------------------------------------------
    // Test 4: Two-term LinearTransform — compute s_i + s_{i+1} per row
    //         diagonals = [(1,...,1), (1,...,1)], offsets = [0, 1]
    // -----------------------------------------------------------------------
    std::cout << "  Test 4: LinearTransform two-term (s + rot(s,1))... " << std::flush;
    {
        std::array<uint64_t, P::n> slots_in{}, slots_out{};
        for (auto &v : slots_in) v = small_dist(engine);  // keep small to avoid t wrap

        std::array<uint64_t, P::n> ones{};
        for (auto &v : ones) v = 1;

        std::vector<std::array<uint64_t, P::n>> diagonals = {ones, ones};
        std::vector<int> offsets = {0, 1};

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::LinearTransform<P>(result, ct, diagonals, offsets, *gk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++) {
            uint64_t expected;
            if (i < half)
                expected = (slots_in[i] + slots_in[(i + 1) % half]) % t;
            else
                expected = (slots_in[i] + slots_in[half + (i - half + 1) % half]) % t;
            if (slots_out[i] != expected) mismatches++;
        }

        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;
        if (mismatches > 0) {
            for (int i = 0; i < 3; i++) {
                uint64_t expected;
                if (i < half)
                    expected = (slots_in[i] + slots_in[(i + 1) % half]) % t;
                else
                    expected = (slots_in[i] + slots_in[half + (i - half + 1) % half]) % t;
                std::cout << "    slot[" << i << "] expected=" << expected
                          << " got=" << slots_out[i] << std::endl;
            }
            failures++;
        }
    }

    // -----------------------------------------------------------------------
    // Test 5: LinearTransform with distinct diagonals
    //         diagonal[0] = a, diagonal[1] = b, offsets = [0, 1]
    //         output[i] = a[i]·s[i] + b[i]·s[(i+1) mod half]
    // -----------------------------------------------------------------------
    std::cout << "  Test 5: LinearTransform distinct diagonals... " << std::flush;
    {
        std::array<uint64_t, P::n> slots_in{}, slots_out{};
        for (auto &v : slots_in) v = small_dist(engine);

        std::array<uint64_t, P::n> a{}, b{};
        for (int i = 0; i < n; i++) {
            a[i] = small_dist(engine);
            b[i] = small_dist(engine);
        }

        std::vector<std::array<uint64_t, P::n>> diagonals = {a, b};
        std::vector<int> offsets = {0, 1};

        TFHEpp::TRLWE<P> ct, result;
        TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
        TFHEpp::LinearTransform<P>(result, ct, diagonals, offsets, *gk);
        TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

        int mismatches = 0;
        for (int i = 0; i < n; i++) {
            uint64_t rot;
            if (i < half) rot = slots_in[(i + 1) % half];
            else rot = slots_in[half + (i - half + 1) % half];
            uint64_t expected = (a[i] * slots_in[i] + b[i] * rot) % t;
            if (slots_out[i] != expected) mismatches++;
        }

        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;
        if (mismatches > 0) {
            for (int i = 0; i < 3; i++) {
                uint64_t rot;
                if (i < half) rot = slots_in[(i + 1) % half];
                else rot = slots_in[half + (i - half + 1) % half];
                uint64_t expected = (a[i] * slots_in[i] + b[i] * rot) % t;
                std::cout << "    slot[" << i << "] s=" << slots_in[i]
                          << " rot_s=" << rot
                          << " a=" << a[i] << " b=" << b[i]
                          << " expected=" << expected
                          << " got=" << slots_out[i] << std::endl;
            }
            failures++;
        }
    }

    if (failures == 0) {
        std::cout << "PASS (all LinearTransform tests)" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL (" << failures << " test(s))" << std::endl;
        return 1;
    }
}
