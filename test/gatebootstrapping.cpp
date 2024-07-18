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
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    using bkP = TFHEpp::lvl02param;
    using iksP = TFHEpp::lvl20param;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<bkP>(sk);
    ek.emplaceiksk<iksP>(sk);
    std::vector<TFHEpp::TLWE<typename iksP::domainP>> tlwe(num_test),
        bootedtlwe(num_test);
    std::array<bool, num_test> p;
    for (int i = 0; i < num_test; i++) p[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++)
        tlwe[i] = TFHEpp::tlweSymEncrypt<typename iksP::domainP>(
            p[i] ? iksP::domainP::μ : -iksP::domainP::μ,
            sk.key.get<typename iksP::domainP>());

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
#ifdef USE_PERF
    ProfilerStart("gb.prof");
#endif
    for (int test = 0; test < num_test; test++) {
        TFHEpp::GateBootstrapping<iksP, bkP, bkP::targetP::μ>(bootedtlwe[test],
                                                              tlwe[test], ek);
    }
#ifdef USE_PERF
    ProfilerStop();
#endif

    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
    for (int i = 0; i < num_test; i++) {
        bool p2 = TFHEpp::tlweSymDecrypt<typename iksP::domainP>(
            bootedtlwe[i], sk.key.get<typename iksP::domainP>());
        assert(p[i] == p2);
    }
    std::cout << "Passed" << std::endl;
}