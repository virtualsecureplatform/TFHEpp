// modulusswitch.cpp — Test BFV-style ModulusSwitch (Q → Qp)
//
// Encrypts a SIMD slot vector at the default BFV scaling (Δ = Q/t),
// mod-switches to a smaller modulus Qp, and verifies that decryption
// with scaling (Qp/t) still recovers the original slots.

#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <tfhe++.hpp>

template <class P>
static int RunOne(uint64_t Qp, const char *label, std::default_random_engine &engine,
                  std::uniform_int_distribution<uint64_t> &slot_dist,
                  const TFHEpp::Key<P> &key)
{
    constexpr int n = static_cast<int>(P::n);

    std::array<uint64_t, P::n> slots_in, slots_out;
    for (auto &v : slots_in) v = slot_dist(engine);

    TFHEpp::TRLWE<P> ct, ct_switched;
    TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);
    TFHEpp::ModulusSwitch<P>(ct_switched, ct, Qp);

    // Verify that all mask/body coefficients are actually reduced into [0, Qp).
    for (int k = 0; k <= static_cast<int>(P::k); k++) {
        for (int i = 0; i < n; i++) {
            if (ct_switched[k][i] >= static_cast<typename P::T>(Qp)) {
                std::cerr << "    FAIL (" << label << "): coefficient out of range"
                          << " k=" << k << " i=" << i << std::endl;
                return 1;
            }
        }
    }

    TFHEpp::trlweSlotDecryptModSwitched<P>(slots_out, ct_switched, key, Qp);

    int mismatches = 0;
    for (int i = 0; i < n; i++)
        if (slots_out[i] != slots_in[i]) mismatches++;

    std::cout << "    " << label << " (Qp=" << Qp << "): "
              << (mismatches == 0 ? "PASS" : "FAIL")
              << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;

    if (mismatches > 0) {
        for (int i = 0; i < 3; i++)
            std::cout << "      slot[" << i << "] in=" << slots_in[i]
                      << " out=" << slots_out[i] << std::endl;
    }
    return mismatches == 0 ? 0 : 1;
}

int main()
{
    using P = TFHEpp::lvl3simdparam;
    using T = typename P::T;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int n = static_cast<int>(P::n);

    std::cout << "ModulusSwitch test (BFV Q → Qp)" << std::endl;
    std::cout << "  n=" << n << "  t=" << t << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, t - 1);
    std::uniform_int_distribution<int> binary(0, 1);

    TFHEpp::Key<P> key;
    for (auto &v : key)
        v = binary(engine) ? static_cast<T>(1) : static_cast<T>(-1);

    int failures = 0;

    // Large Qp — plenty of noise budget.
    failures += RunOne<P>(static_cast<uint64_t>(1) << 60, "Qp=2^60", engine, slot_dist, key);

    // Power of two close to Q_mod_t threshold.
    failures += RunOne<P>(static_cast<uint64_t>(1) << 50, "Qp=2^50", engine, slot_dist, key);

    // Odd prime-like Qp: 2^50 + 33 (odd, coprime to 2^128).
    failures += RunOne<P>((static_cast<uint64_t>(1) << 50) + 33, "Qp=2^50+33", engine, slot_dist, key);

    // Small Qp near the threshold: must be large enough that fresh noise
    // (~ α·Q·√n ≈ 2^29 for lvl3simdparam) rescaled by Qp/Q plus rounding
    // error (~n) stays below Qp/(2t). Need Qp > 2·t·n ≈ 2^30.
    failures += RunOne<P>(static_cast<uint64_t>(1) << 40, "Qp=2^40", engine, slot_dist, key);

    // t·2^16: Qp/t = 2^16 of headroom, about the minimum that works.
    failures += RunOne<P>(static_cast<uint64_t>(t) << 16, "Qp=t·2^16", engine, slot_dist, key);

    // Tighter: t·2^13, near the noise boundary for lvl3simdparam.
    failures += RunOne<P>(static_cast<uint64_t>(t) << 13, "Qp=t·2^13", engine, slot_dist, key);

    if (failures == 0)
        std::cout << "PASS (all ModulusSwitch tests)" << std::endl;
    else
        std::cout << "FAIL (" << failures << " cases failed)" << std::endl;
    return failures == 0 ? 0 : 1;
}
