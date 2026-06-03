// c2s_factorization.cpp — Plaintext verification of the CoeffToSlot
// decomposition.
//
// Computes the full SlotEncode matrix S (slot-domain target for a
// homomorphic CoeffToSlot), splits it into same-row and cross-row parts,
// and checks that
//
//    S · v  ==  PlainLinearTransform(v; D_same)
//               + PlainLinearTransform(Conj(v); D_cross)
//
// holds exactly for several random slot vectors v.  If this passes, the
// homomorphic version (replacing the plaintext linear transform with
// LinearTransformBSGS and Conj with ConjugateSlots) will produce the right
// ciphertext up to the usual BFV noise.

#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl3simdparam;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int n = static_cast<int>(P::n);

    std::cout << "CoeffToSlot factorization (plaintext verification)" << std::endl;
    std::cout << "  n=" << n << "  t=" << t << std::endl;

    std::cout << "  Building SlotEncode matrix S (" << n << "×" << n
              << " uint64, ~" << (static_cast<uint64_t>(n) * n * 8 / (1 << 20))
              << " MB)... " << std::flush;
    auto S = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    TFHEpp::c2s::ComputeSlotEncodeMatrix<P>(*S);
    std::cout << "done" << std::endl;

    std::cout << "  Extracting same/cross-row diagonals... " << std::flush;
    auto D_same  = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    auto D_cross = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    TFHEpp::c2s::ExtractSameRowDiagonals<P>(*D_same,  *S);
    TFHEpp::c2s::ExtractCrossRowDiagonals<P>(*D_cross, *S);
    std::cout << "done"
              << " (same=" << D_same->size()
              << ", cross=" << D_cross->size() << ")" << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, t - 1);

    constexpr int num_trials = 4;
    int failures = 0;
    for (int trial = 0; trial < num_trials; trial++) {
        std::array<uint64_t, P::n> v{}, y_ref{}, y_decomp{};
        for (auto &x : v) x = slot_dist(engine);

        // Reference: apply S directly via matrix-vector multiplication.
        TFHEpp::c2s::MatVecMulMod<P>(y_ref, *S, v);

        // Decomposition: same-row + cross-row via conjugation.
        TFHEpp::c2s::PlainApplyS<P>(y_decomp, v, *D_same, *D_cross);

        int mismatches = 0;
        int first_bad = -1;
        for (int i = 0; i < n; i++) {
            if (y_ref[i] != y_decomp[i]) {
                if (first_bad < 0) first_bad = i;
                mismatches++;
            }
        }
        std::cout << "  Trial " << trial << ": "
                  << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << "/" << n << ")";
        if (mismatches > 0 && first_bad >= 0)
            std::cout << " first_bad i=" << first_bad
                      << " ref=" << y_ref[first_bad]
                      << " decomp=" << y_decomp[first_bad];
        std::cout << std::endl;
        if (mismatches > 0) failures++;
    }

    // Sanity check: SlotEncode(v) should equal S · v for any v.
    std::cout << "  SlotEncode(v) vs S·v self-consistency... " << std::flush;
    {
        std::array<uint64_t, P::n> v{}, y_matvec{};
        for (auto &x : v) x = slot_dist(engine);
        TFHEpp::c2s::MatVecMulMod<P>(y_matvec, *S, v);

        TFHEpp::Polynomial<P> poly;
        TFHEpp::SlotEncode<P>(poly, v);

        int mismatches = 0;
        for (int i = 0; i < n; i++)
            if (static_cast<uint64_t>(poly[i]) != y_matvec[i]) mismatches++;
        std::cout << (mismatches == 0 ? "PASS" : "FAIL")
                  << " (mismatches=" << mismatches << "/" << n << ")" << std::endl;
        if (mismatches > 0) failures++;
    }

    if (failures == 0) {
        std::cout << "PASS (CoeffToSlot factorization verified in plaintext)" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL (" << failures << " trial(s) mismatched)" << std::endl;
        return 1;
    }
}
