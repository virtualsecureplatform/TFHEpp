// coeff_to_slot.cpp — Homomorphic CoeffToSlot end-to-end test.
//
// Encrypts a polynomial P directly (BFV-encoded, phase ≈ Δ·P), applies
// CoeffToSlot via LinearTransformBSGS + ConjugateSlots, and verifies that
// the resulting ciphertext decrypts to slots = (p_0, p_1, ..., p_{n-1}).

#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <tfhe++.hpp>

// Encrypt a polynomial P directly (no slot permutation) with BFV scaling.
// Phase after encryption ≈ Δ · P + noise, so trlweSlotDecrypt will recover
// SlotDecode(P) = NTT of P's coefficients under our GaloisPermutation.
template <class P_>
static void BfvPolyEncrypt(TFHEpp::TRLWE<P_> &ct,
                           const std::array<uint64_t, P_::n> &poly_coefs,
                           const TFHEpp::Key<P_> &key)
{
    constexpr typename P_::T delta = P_::delta_int;
    constexpr uint64_t r = P_::Q_mod_t;
    constexpr uint64_t t_val = static_cast<uint64_t>(P_::plain_modulus);

    TFHEpp::Polynomial<P_> scaled;
    for (uint32_t i = 0; i < P_::n; i++) {
        uint64_t m = poly_coefs[i] % t_val;
        scaled[i] = static_cast<typename P_::T>(m) * delta
                  + static_cast<typename P_::T>(m * r / t_val);
    }
    TFHEpp::trlweSymEncrypt<P_>(ct, scaled, key);
}

int main()
{
    using P = TFHEpp::lvl3simdparam;
    using T = typename P::T;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int n = static_cast<int>(P::n);
    constexpr int half = n / 2;

    std::cout << "CoeffToSlot (homomorphic) test" << std::endl;
    std::cout << "  n=" << n << "  t=" << t << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> small_dist(0, 16);
    std::uniform_int_distribution<int> binary(0, 1);

    TFHEpp::Key<P> key;
    for (auto &v : key)
        v = binary(engine) ? static_cast<T>(1) : static_cast<T>(-1);

    std::cout << "  Generating GaloisKey (" << (P::nbit + 1) << " keys)... " << std::flush;
    auto gk = std::make_unique<TFHEpp::GaloisKey<P>>();
    TFHEpp::GaloisKeyGen<P>(*gk, key);
    std::cout << "done" << std::endl;

    std::cout << "  Building CoeffToSlot diagonals (~128 MB matrix)... " << std::flush;
    auto D_same  = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    auto D_cross = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    TFHEpp::c2s::BuildCoeffToSlotKeys<P>(*D_same, *D_cross);
    std::cout << "done" << std::endl;

    // Random polynomial with small coefficients (keeps SlotPtxtMul noise low).
    std::array<uint64_t, P::n> poly_coefs{}, slots_out{};
    for (auto &c : poly_coefs) c = small_dist(engine);

    std::cout << "  Encrypting polynomial... " << std::flush;
    TFHEpp::TRLWE<P> ct_in, ct_out;
    BfvPolyEncrypt<P>(ct_in, poly_coefs, key);
    std::cout << "done" << std::endl;

    std::cout << "  Running CoeffToSlot (BSGS k=64)... " << std::flush;
    TFHEpp::c2s::CoeffToSlot<P>(ct_out, ct_in, *D_same, *D_cross, 64, *gk);
    std::cout << "done" << std::endl;

    std::cout << "  Decrypting output slots... " << std::flush;
    TFHEpp::trlweSlotDecrypt<P>(slots_out, ct_out, key);
    std::cout << "done" << std::endl;

    int mismatches = 0;
    int first_bad = -1;
    for (int i = 0; i < n; i++) {
        if (slots_out[i] != poly_coefs[i]) {
            if (first_bad < 0) first_bad = i;
            mismatches++;
        }
    }

    std::cout << "  Result: " << (mismatches == 0 ? "PASS" : "FAIL")
              << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;
    if (mismatches > 0) {
        std::cout << "  First bad slot i=" << first_bad
                  << " expected=" << poly_coefs[first_bad]
                  << " got=" << slots_out[first_bad] << std::endl;
        for (int i = 0; i < 5; i++)
            std::cout << "    slot[" << i << "] expected=" << poly_coefs[i]
                      << " got=" << slots_out[i] << std::endl;
        return 1;
    }
    std::cout << "PASS (CoeffToSlot homomorphic test)" << std::endl;
    return 0;
}
