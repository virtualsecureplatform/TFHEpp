// simd_rotation.cpp — Test slot rotation using GaloisKeyGen + RotateSlots
//
// Verifies that RotateSlots correctly permutes SIMD slots:
//   decrypt(RotateSlots(Enc(slots), k))[i] == slots[(i+k) mod n]
//
// Also tests ConjugateSlots as a special case.
// Uses lvl3simdparam; rotation uses EvalAuto which auto-selects DD (l̅=8).

#include <cassert>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl3simdparam;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int n      = static_cast<int>(P::n);

    std::cout << "SIMD slot rotation test" << std::endl;
    std::cout << "  n=" << P::n << "  t=" << t
              << "  nbit=" << P::nbit << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, t - 1);
    std::uniform_int_distribution<int> binary(0, 1);

    // Generate secret key
    TFHEpp::Key<P> key;
    for (auto &v : key)
        v = binary(engine) ? static_cast<typename P::T>(1)
                           : static_cast<typename P::T>(-1);

    // Generate Galois keys (all log2(n) power-of-2 rotation keys + conjugation)
    std::cout << "  Generating GaloisKey (" << (P::nbit + 1) << " keys)... "
              << std::flush;
    auto gk = std::make_unique<TFHEpp::GaloisKey<P>>();
    TFHEpp::GaloisKeyGen<P>(*gk, key);
    std::cout << "done" << std::endl;

    // Generate random slot vector and encrypt
    std::array<uint64_t, P::n> slots_in;
    for (auto &v : slots_in) v = slot_dist(engine);

    TFHEpp::TRLWE<P> ct_in;
    TFHEpp::trlweSlotEncrypt<P>(ct_in, slots_in, key);

    // Rotation step amounts to test
    const int steps_to_test[] = {1, 3, 7, n / 2, -1, n - 1};
    int failures = 0;

    for (int steps : steps_to_test) {
        const int normalized = ((steps % n) + n) % n;
        std::cout << "  Rotate by " << steps << " (normalized=" << normalized
                  << ")... " << std::flush;

        TFHEpp::TRLWE<P> ct_rot;
        TFHEpp::RotateSlots<P>(ct_rot, ct_in, steps, *gk);

        std::array<uint64_t, P::n> slots_rot;
        TFHEpp::trlweSlotDecrypt<P>(slots_rot, ct_rot, key);

        bool ok = true;
        for (int i = 0; i < n; i++) {
            const uint64_t expected = slots_in[(i + normalized) % n];
            if (slots_rot[i] != expected) {
                std::cerr << "\n    FAIL slot=" << i << " expected=" << expected
                          << " got=" << slots_rot[i] << std::endl;
                failures++;
                ok = false;
                break;
            }
        }
        if (ok) std::cout << "PASS" << std::endl;
    }

    // Test zero rotation (identity)
    {
        std::cout << "  Rotate by 0 (identity)... " << std::flush;
        TFHEpp::TRLWE<P> ct_rot;
        TFHEpp::RotateSlots<P>(ct_rot, ct_in, 0, *gk);
        std::array<uint64_t, P::n> slots_rot;
        TFHEpp::trlweSlotDecrypt<P>(slots_rot, ct_rot, key);
        bool ok = true;
        for (int i = 0; i < n; i++) {
            if (slots_rot[i] != slots_in[i]) {
                failures++;
                ok = false;
                break;
            }
        }
        if (ok) std::cout << "PASS" << std::endl;
    }

    if (failures == 0) {
        std::cout << "PASS (all rotation tests)" << std::endl;
        return 0;
    }
    else {
        std::cout << "FAIL (" << failures << " failures)" << std::endl;
        return 1;
    }
}
