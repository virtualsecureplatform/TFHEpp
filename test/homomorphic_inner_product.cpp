// homomorphic_inner_product.cpp — Test HomomorphicInnerProduct for BFV bootstrapping
//
// Encrypts a slot vector, runs HomomorphicInnerProduct (using HalfTRGSW boot key),
// verifies the result decrypts to the original plaintext.

#include <cassert>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    using P = TFHEpp::lvl3simdparam;
    using T = typename P::T;
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int n = static_cast<int>(P::n);

    std::cout << "HomomorphicInnerProduct test (HalfTRGSW boot key)" << std::endl;
    std::cout << "  n=" << P::n << "  t=" << t << std::endl;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint64_t> slot_dist(0, t - 1);
    std::uniform_int_distribution<int> binary(0, 1);

    TFHEpp::Key<P> key;
    for (auto &v : key)
        v = binary(engine) ? static_cast<T>(1) : static_cast<T>(-1);

    std::cout << "  Generating boot key (HalfTRGSWFFT)... " << std::flush;
    auto bootkeyfft = TFHEpp::makeBootKeyFFT<P>(key);
    std::cout << "done" << std::endl;

    std::cout << "  Testing HomomorphicInnerProduct... " << std::flush;
    {
        constexpr int num_test = 3;
        int failures = 0;
        for (int test = 0; test < num_test; test++) {
            std::array<uint64_t, P::n> slots_in, slots_out;
            for (auto &v : slots_in) v = slot_dist(engine);

            TFHEpp::TRLWE<P> ct;
            TFHEpp::trlweSlotEncrypt<P>(ct, slots_in, key);

            TFHEpp::TRLWE<P> result;
            TFHEpp::HomomorphicInnerProduct<P>(result, ct, *bootkeyfft);

            TFHEpp::trlweSlotDecrypt<P>(slots_out, result, key);

            bool ok = true;
            for (int i = 0; i < n; i++) {
                if (slots_out[i] != slots_in[i]) {
                    std::cerr << "\n    FAIL test=" << test << " slot=" << i
                              << " expected=" << slots_in[i]
                              << " got=" << slots_out[i] << std::endl;
                    ok = false;
                    failures++;
                    break;
                }
            }
            if (ok) std::cout << "." << std::flush;
        }
        if (failures == 0)
            std::cout << " PASS" << std::endl;
        else {
            std::cout << " FAIL (" << failures << "/" << num_test << " failed)" << std::endl;
            return 1;
        }
    }

    std::cout << "PASS (HomomorphicInnerProduct test)" << std::endl;
    return 0;
}
