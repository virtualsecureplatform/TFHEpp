// pbs_tfhers.cpp — TFHEpp gate-bootstrapping benchmark for comparison with
// tfhe-rs DEFAULT_PARAMETERS.
//
// Measures NAND gate throughput using the tfhe-rs-compatible parameter set
// (lvl01param for PBS, lvl10param for key switching):
//   n=805, k=3, N=512, pbs_level=2, pbs_base_log=10, ks_level=5, ks_base_log=3
//
// Build with -DUSE_TFHE_RS=ON.
// Matching Rust benchmark: bench-tfhers/src/main.rs

#include <cassert>
#include <chrono>
#include <iostream>
#include <random>

#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;

    // Gate bootstrapping uses lvl01param (lvl0→lvl1) + lvl10param (lvl1→lvl0)
    // which is the default for HomNAND.
    using brP = TFHEpp::lvl01param;
    using iksP = TFHEpp::lvl10param;
    using LWEin = TFHEpp::TLWE<typename brP::domainP>;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);
    ek.emplaceiksk<iksP>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    std::vector<LWEin> ca(num_test), cb(num_test), res(num_test);
    std::array<bool, num_test> pa, pb;

    for (uint32_t i = 0; i < num_test; i++) {
        pa[i] = binary(engine) > 0;
        pb[i] = binary(engine) > 0;
        TFHEpp::tlweSymEncrypt<typename brP::domainP>(
            ca[i], pa[i] ? brP::domainP::μ : -brP::domainP::μ,
            sk.key.get<typename brP::domainP>());
        TFHEpp::tlweSymEncrypt<typename brP::domainP>(
            cb[i], pb[i] ? brP::domainP::μ : -brP::domainP::μ,
            sk.key.get<typename brP::domainP>());
    }

    // Warmup
    TFHEpp::HomNAND<brP, brP::targetP::μ, iksP>(res[0], ca[0], cb[0], ek);

    auto start = std::chrono::system_clock::now();
    for (uint32_t i = 0; i < num_test; i++)
        TFHEpp::HomNAND<brP, brP::targetP::μ, iksP>(res[i], ca[i], cb[i], ek);
    auto end = std::chrono::system_clock::now();

    double elapsed_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();

    std::cout << elapsed_ms / num_test << "ms" << std::endl;

    // Correctness check
    for (uint32_t i = 0; i < num_test; i++) {
        assert(TFHEpp::tlweSymDecrypt<typename iksP::targetP>(
                   res[i], sk.key.get<typename iksP::targetP>()) ==
               !(pa[i] && pb[i]));
    }
    std::cout << "Passed" << std::endl;
}
