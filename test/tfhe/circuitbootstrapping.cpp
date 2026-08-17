#ifdef USE_PERF
#include <gperftools/profiler.h>
#endif

#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <sstream>
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
    std::stringstream serialized_secret{std::ios::binary | std::ios::out |
                                        std::ios::in};
    {
        cereal::PortableBinaryOutputArchive archive{serialized_secret};
        archive(*sk);
    }
    {
        cereal::PortableBinaryInputArchive archive{serialized_secret};
        archive(*sk);
    }
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebk<TFHEpp::lvl01param>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    std::stringstream serialized_key{std::ios::binary | std::ios::out |
                                     std::ios::in};
    {
        cereal::PortableBinaryOutputArchive archive{serialized_key};
        archive(ek);
    }
    TFHEpp::EvalKey serialized_ek;
    {
        cereal::PortableBinaryInputArchive archive{serialized_key};
        archive(serialized_ek);
    }

    std::vector<std::array<uint8_t, privksP::targetP::n>> pa(num_test);
    std::vector<std::array<typename privksP::targetP::T, privksP::targetP::n>>
        pmu(num_test);
    std::vector<uint8_t> pones(num_test);
    std::vector<uint8_t> pzeros(num_test, 0);
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
    alignas(64) std::vector<TFHEpp::TLWE<typename bkP::domainP>> directcones(
        num_test);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>,
                TFHEpp::AlignedAllocator<
                    TFHEpp::TRGSWFFT<typename privksP::targetP>, 64>>
        bootedTGSW(num_test);

    for (int i = 0; i < num_test; i++)
        TFHEpp::trlweSymEncrypt<typename privksP::targetP>(
            ca[i], pmu[i], sk->key.get<typename privksP::targetP>());
    TFHEpp::bootsSymEncrypt<typename iksP::domainP>(cones, pones, *sk);
    TFHEpp::bootsSymEncrypt<typename bkP::domainP>(directcones, pones, *sk);

    std::chrono::system_clock::time_point start, end;
#ifdef USE_PERF
    ProfilerStart("cb.prof");
#endif
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        TFHEpp::CircuitBootstrapping<iksP, bkP, privksP>(bootedTGSW[test],
                                                         cones[test], serialized_ek);
        TFHEpp::CircuitBootstrapping<bkP, privksP>(bootedTGSW[test],
                                                    directcones[test], serialized_ek);
    }
    end = std::chrono::system_clock::now();
#ifdef USE_PERF
    ProfilerStop();
#endif
    for (int test = 0; test < num_test; test++) {
        TFHEpp::ExternalProduct<typename privksP::targetP>(
            ca[test], ca[test], bootedTGSW[test]);
        pres = TFHEpp::trlweSymDecrypt<typename privksP::targetP>(
            ca[test], sk->key.get<typename privksP::targetP>());
        for (int i = 0; i < privksP::targetP::n; i++)
            assert(pres[i] == pa[test][i]);

        TFHEpp::TLWE<TFHEpp::lvl1param> extracted;
        TFHEpp::TLWE<TFHEpp::lvl0param> switched;
        TFHEpp::SampleExtractIndex<TFHEpp::lvl1param>(extracted, ca[test], 0);
        TFHEpp::IdentityKeySwitch<iksP>(
            switched, extracted, serialized_ek.getiksk<iksP>());
        assert(TFHEpp::tlweSymDecrypt<TFHEpp::lvl0param>(switched, *sk) ==
               pa[test][0]);
    }

    std::vector<TFHEpp::TLWE<typename bkP::domainP>> directzeros;
    TFHEpp::bootsSymEncrypt<typename bkP::domainP>(directzeros, pzeros, *sk);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>,
                TFHEpp::AlignedAllocator<
                    TFHEpp::TRGSWFFT<typename privksP::targetP>, 64>>
        directZeroTGSW(num_test), directZeroInvTGSW(num_test);
    for (int test = 0; test < num_test; test++) {
        TFHEpp::CircuitBootstrappingWithInv<bkP, privksP>(
            directZeroTGSW[test], directZeroInvTGSW[test], directzeros[test],
            serialized_ek);

        TFHEpp::TRLWE<typename privksP::targetP> normal, inverted;
        TFHEpp::trlweSymEncrypt<typename privksP::targetP>(
            normal, pmu[test], sk->key.get<typename privksP::targetP>());
        inverted = normal;
        TFHEpp::ExternalProduct<typename privksP::targetP>(
            normal, normal, directZeroTGSW[test]);
        TFHEpp::ExternalProduct<typename privksP::targetP>(
            inverted, inverted, directZeroInvTGSW[test]);
        pres = TFHEpp::trlweSymDecrypt<typename privksP::targetP>(
            normal, sk->key.get<typename privksP::targetP>());
        for (int i = 0; i < privksP::targetP::n; i++) assert(!pres[i]);
        pres = TFHEpp::trlweSymDecrypt<typename privksP::targetP>(
            inverted, sk->key.get<typename privksP::targetP>());
        for (int i = 0; i < privksP::targetP::n; i++)
            assert(pres[i] == pa[test][i]);
    }

    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>,
                TFHEpp::AlignedAllocator<
                    TFHEpp::TRGSWFFT<typename privksP::targetP>, 64>>
        directTGSW(num_test), directInvTGSW(num_test);
    for (int test = 0; test < num_test; test++) {
        TFHEpp::CircuitBootstrappingWithInv<bkP, privksP>(
            directTGSW[test], directInvTGSW[test], directcones[test],
            serialized_ek);

        TFHEpp::ExternalProduct<typename privksP::targetP>(
            ca[test], ca[test], directTGSW[test]);
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
