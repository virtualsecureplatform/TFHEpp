// pbs_tfhers.cpp — TFHEpp gate-bootstrapping benchmark for comparison with
// tfhe-rs DEFAULT_PARAMETERS.
//
// Measures BlindRotate (BR), IdentityKeySwitch (IKS), and full NAND gate
// throughput using the tfhe-rs-compatible parameter set
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
    using lvl0param = typename brP::domainP;
    using lvl1param = typename brP::targetP;
    using LWEin = TFHEpp::TLWE<lvl0param>;
    using LWEbig = TFHEpp::TLWE<lvl1param>;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);
    ek.emplaceiksk<iksP>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    std::vector<LWEin> ca(num_test), cb(num_test);
    std::array<bool, num_test> pa, pb;

    for (uint32_t i = 0; i < num_test; i++) {
        pa[i] = binary(engine) > 0;
        pb[i] = binary(engine) > 0;
        TFHEpp::tlweSymEncrypt<lvl0param>(
            ca[i], pa[i] ? lvl0param::μ : -lvl0param::μ,
            sk.key.get<lvl0param>());
        TFHEpp::tlweSymEncrypt<lvl0param>(
            cb[i], pb[i] ? lvl0param::μ : -lvl0param::μ,
            sk.key.get<lvl0param>());
    }

    // Prepare NAND inputs: -ca - cb + μ (replicating HomGate with casign=-1,
    // cbsign=-1, offset=μ)
    std::vector<LWEin> prepared(num_test);
    for (uint32_t i = 0; i < num_test; i++) {
        for (int j = 0; j <= lvl0param::k * lvl0param::n; j++)
            prepared[i][j] = -ca[i][j] - cb[i][j];
        prepared[i][lvl0param::k * lvl0param::n] += lvl0param::μ;
    }

    const auto testvector =
        TFHEpp::μpolygen<lvl1param, lvl1param::μ>();

    // Warmup
    {
        LWEbig tmp;
        TFHEpp::GateBootstrappingTLWE2TLWE<brP>(
            tmp, prepared[0], ek.getbkfft<brP>(), testvector);
        LWEin tmp2;
        TFHEpp::IdentityKeySwitch<iksP>(tmp2, tmp, ek.getiksk<iksP>());
    }

    // --- BR benchmark ---
    std::vector<LWEbig> tlwelvl1(num_test);
    auto br_start = std::chrono::system_clock::now();
    for (uint32_t i = 0; i < num_test; i++)
        TFHEpp::GateBootstrappingTLWE2TLWE<brP>(
            tlwelvl1[i], prepared[i], ek.getbkfft<brP>(), testvector);
    auto br_end = std::chrono::system_clock::now();

    // --- IKS benchmark ---
    std::vector<LWEin> res(num_test);
    auto iks_start = std::chrono::system_clock::now();
    for (uint32_t i = 0; i < num_test; i++)
        TFHEpp::IdentityKeySwitch<iksP>(
            res[i], tlwelvl1[i], ek.getiksk<iksP>());
    auto iks_end = std::chrono::system_clock::now();

    // --- Full NAND benchmark ---
    std::vector<LWEin> nand_res(num_test);
    auto nand_start = std::chrono::system_clock::now();
    for (uint32_t i = 0; i < num_test; i++)
        TFHEpp::HomNAND<brP, brP::targetP::μ, iksP>(
            nand_res[i], ca[i], cb[i], ek);
    auto nand_end = std::chrono::system_clock::now();

    double br_ms =
        std::chrono::duration<double, std::milli>(br_end - br_start).count();
    double iks_ms =
        std::chrono::duration<double, std::milli>(iks_end - iks_start).count();
    double nand_ms =
        std::chrono::duration<double, std::milli>(nand_end - nand_start)
            .count();

    std::cout << "BR: " << br_ms / num_test << "ms" << std::endl;
    std::cout << "IKS: " << iks_ms / num_test << "ms" << std::endl;
    std::cout << "NAND: " << nand_ms / num_test << "ms" << std::endl;

    // Correctness check (split BR+IKS)
    for (uint32_t i = 0; i < num_test; i++) {
        assert(TFHEpp::tlweSymDecrypt<typename iksP::targetP>(
                   res[i], sk.key.get<typename iksP::targetP>()) ==
               !(pa[i] && pb[i]));
    }
    // Correctness check (full NAND)
    for (uint32_t i = 0; i < num_test; i++) {
        assert(TFHEpp::tlweSymDecrypt<typename iksP::targetP>(
                   nand_res[i], sk.key.get<typename iksP::targetP>()) ==
               !(pa[i] && pb[i]));
    }
    std::cout << "Passed" << std::endl;
}
