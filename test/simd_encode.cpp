// simd_encode.cpp — Test SlotEncode / SlotDecode round-trip for lvl3simdparam
//
// Verifies that encoding a slot vector to a polynomial and decoding it back
// recovers the original values exactly, exercising NTTmod / INTTmod correctness.

#include <cassert>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl3simdparam;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);

    constexpr int num_test = 100;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, t - 1);

    std::cout << "SlotEncode / SlotDecode round-trip test" << std::endl;
    std::cout << "  n=" << P::n << "  t=" << t << std::endl;

    int failures = 0;
    for (int test = 0; test < num_test; test++) {
        std::array<uint64_t, P::n> slots_in, slots_out;
        for (uint64_t &v : slots_in) v = slot_dist(engine);

        TFHEpp::Polynomial<P> poly;
        TFHEpp::SlotEncode<P>(poly, slots_in);
        TFHEpp::SlotDecode<P>(slots_out, poly);

        for (int i = 0; i < static_cast<int>(P::n); i++) {
            if (slots_in[i] != slots_out[i]) {
                std::cerr << "  FAIL test=" << test << " slot=" << i
                          << " expected=" << slots_in[i]
                          << " got=" << slots_out[i] << std::endl;
                failures++;
                break;
            }
        }
    }

    // Test with all-zero slots
    {
        std::array<uint64_t, P::n> zeros{}, out;
        TFHEpp::Polynomial<P> poly;
        TFHEpp::SlotEncode<P>(poly, zeros);
        TFHEpp::SlotDecode<P>(out, poly);
        for (int i = 0; i < static_cast<int>(P::n); i++)
            assert(out[i] == 0);
    }

    // Test with all-one slots
    {
        std::array<uint64_t, P::n> ones, out;
        ones.fill(1);
        TFHEpp::Polynomial<P> poly;
        TFHEpp::SlotEncode<P>(poly, ones);
        TFHEpp::SlotDecode<P>(out, poly);
        for (int i = 0; i < static_cast<int>(P::n); i++)
            assert(out[i] == 1);
    }

    if (failures == 0) {
        std::cout << "PASS (" << num_test << " random + 2 edge-case tests)" << std::endl;
        return 0;
    }
    else {
        std::cout << "FAIL (" << failures << " failures)" << std::endl;
        return 1;
    }
}
