#ifdef USE_PERF
#include <gperftools/profiler.h>
#endif

#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 10;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl21param;

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    std::vector<std::array<uint8_t, privksP::targetP::n>> pa(num_test);
    std::vector<std::array<typename privksP::targetP::T, privksP::targetP::n>>
        pmu(num_test);
    std::vector<uint8_t> pones(num_test);
    std::array<bool, privksP::targetP::n> pres;
    for (std::array<uint8_t, privksP::targetP::n> &i : pa)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < privksP::targetP::n; j++)
            pmu[i][j] = pa[i][j] ? privksP::targetP::μ : -privksP::targetP::μ;
    for (int i = 0; i < num_test; i++) pones[i] = true;
    alignas(64) std::vector<TFHEpp::TRLWE<typename privksP::targetP>> ca(
        num_test);
    alignas(64) std::vector<TFHEpp::TLWE<typename iksP::domainP>> cones(
        num_test);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> bootedTGSW(
        num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = TFHEpp::trlweSymEncrypt<typename privksP::targetP>(
            pmu[i], sk->key.get<typename privksP::targetP>());
    cones = TFHEpp::bootsSymEncrypt<typename iksP::domainP>(pones, *sk);

    std::chrono::system_clock::time_point start, end;
#ifdef USE_PERF
    ProfilerStart("cb.prof");
#endif
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW[test],
                                                            cones[test], ek);
    }
    end = std::chrono::system_clock::now();
#ifdef USE_PERF
    ProfilerStop();
#endif
    for (int test = 0; test < num_test; test++) {
        TFHEpp::trgswfftExternalProduct<typename privksP::targetP>(
            ca[test], ca[test], bootedTGSW[test]);
        pres = TFHEpp::trlweSymDecrypt<typename privksP::targetP>(
            ca[test], sk->key.get<typename privksP::targetP>());
        for (int i = 0; i < privksP::targetP::n; i++)
            assert(pres[i] == pa[test][i]);
    }
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}
