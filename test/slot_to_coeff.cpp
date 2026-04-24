// slot_to_coeff.cpp — Homomorphic SlotToCoeff end-to-end test.
//
// Encrypts a random slot vector with trlweSlotEncrypt, applies
// SlotToCoeff (which applies the SlotDecode matrix S^{-1} in slot domain),
// and verifies that the output ciphertext's polynomial has BFV plaintext
// coefficients equal to the original slot vector.

#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <tfhe++.hpp>

// Decrypt a TRLWE as a BFV-encoded polynomial (no slot permutation).
// Returns the round(phase · t / Q) per coefficient, i.e. the "plaintext
// polynomial" before SlotDecode.
template <class P_>
static void BfvPolyDecrypt(std::array<uint64_t, P_::n> &poly,
                           const TFHEpp::TRLWE<P_> &ct,
                           const TFHEpp::Key<P_> &key)
{
    TFHEpp::Polynomial<P_> phase = TFHEpp::trlwePhase<P_>(ct, key);
    for (uint32_t i = 0; i < P_::n; i++)
        poly[i] = TFHEpp::bfvDecodeCoeff<P_>(phase[i]);
}

int main()
{
    using P = TFHEpp::lvl3simdparam;
    using T = typename P::T;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int n = static_cast<int>(P::n);

    std::cout << "SlotToCoeff (homomorphic) test" << std::endl;
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

    std::cout << "  Building SlotToCoeff diagonals (~128 MB matrix)... " << std::flush;
    auto D_same  = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    auto D_cross = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    TFHEpp::c2s::BuildSlotToCoeffKeys<P>(*D_same, *D_cross);
    std::cout << "done" << std::endl;

    // Random slot vector with small values (keeps SlotPtxtMul noise low).
    std::array<uint64_t, P::n> slots_in{}, poly_out{};
    for (auto &v : slots_in) v = small_dist(engine);

    std::cout << "  Encrypting slots... " << std::flush;
    TFHEpp::TRLWE<P> ct_in, ct_out;
    TFHEpp::trlweSlotEncrypt<P>(ct_in, slots_in, key);
    std::cout << "done" << std::endl;

    std::cout << "  Running SlotToCoeff (BSGS k=64)... " << std::flush;
    TFHEpp::c2s::SlotToCoeff<P>(ct_out, ct_in, *D_same, *D_cross, 64, *gk);
    std::cout << "done" << std::endl;

    std::cout << "  Decrypting output as polynomial... " << std::flush;
    BfvPolyDecrypt<P>(poly_out, ct_out, key);
    std::cout << "done" << std::endl;

    int mismatches = 0;
    int first_bad = -1;
    for (int i = 0; i < n; i++) {
        if (poly_out[i] != slots_in[i]) {
            if (first_bad < 0) first_bad = i;
            mismatches++;
        }
    }

    std::cout << "  Result: " << (mismatches == 0 ? "PASS" : "FAIL")
              << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;
    if (mismatches > 0) {
        std::cout << "  First bad coef i=" << first_bad
                  << " expected=" << slots_in[first_bad]
                  << " got=" << poly_out[first_bad] << std::endl;
        for (int i = 0; i < 5; i++)
            std::cout << "    coef[" << i << "] expected=" << slots_in[i]
                      << " got=" << poly_out[i] << std::endl;
        return 1;
    }
    std::cout << "PASS (SlotToCoeff homomorphic test)" << std::endl;
    return 0;
}
