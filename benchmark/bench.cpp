#include <memory>
#include <random>

#include "../include/tfhe++.hpp"
#include "google-benchmark/include/benchmark/benchmark.h"

void BM_TRGSWenc(benchmark::State& state)
{
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    TFHEpp::TRGSWFFT<TFHEpp::lvl1param> res;
    for (auto _ : state)
        res = TFHEpp::trgswfftSymEncrypt<TFHEpp::lvl1param>(
            {}, TFHEpp::lvl1param::α, sk->key.lvl1);
    ;
}

void BM_HomGate(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    TFHEpp::TLWE<TFHEpp::lvl0param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
            binary(engine), TFHEpp::lvl0param::α, sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> cb =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
            binary(engine), TFHEpp::lvl0param::α, sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> res;
    for (auto _ : state) TFHEpp::HomNAND<TFHEpp::lvl0param>(res, ca, cb, ek);
}

void BM_HomMUX(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    TFHEpp::TLWE<TFHEpp::lvl0param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
            binary(engine), TFHEpp::lvl0param::α, sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> cb =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
            binary(engine), TFHEpp::lvl0param::α, sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> cs =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
            binary(engine), TFHEpp::lvl0param::α, sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> res;
    for (auto _ : state) TFHEpp::HomMUX<TFHEpp::lvl0param>(res, cs, ca, cb, ek);
}

void BM_TLWE2TRLWE(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    TFHEpp::TLWE<TFHEpp::lvl0param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
            binary(engine), TFHEpp::lvl0param::α, sk->key.lvl0);
    TFHEpp::TRLWE<TFHEpp::lvl1param> res;
    for (auto _ : state)
        TFHEpp::BlindRotate<TFHEpp::lvl01param>(
            res, ca, *ek.bkfftlvl01,
            TFHEpp::μpolygen<TFHEpp::lvl1param, TFHEpp::lvl1param::μ>());
}

void BM_IKS(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    TFHEpp::TLWE<TFHEpp::lvl1param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine), TFHEpp::lvl1param::α, sk->key.lvl1);
    TFHEpp::TLWE<TFHEpp::lvl0param> res;
    for (auto _ : state)
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl10param>(res, ca, *ek.iksklvl10);
}

void BM_SEI(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    std::array<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::n> pmu;
    for (int j = 0; j < TFHEpp::lvl1param::n; j++)
        pmu[j] = binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;
    TFHEpp::TRLWE<TFHEpp::lvl1param> ca =
        TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(pmu, TFHEpp::lvl1param::α,
                                                   sk->key.lvl1);
    TFHEpp::TLWE<TFHEpp::lvl1param> res;
    for (auto _ : state)
        TFHEpp::SampleExtractIndex<TFHEpp::lvl1param>(res, ca, 0);
}

void BM_CMUX(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    std::array<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::n> pmu1, pmu0;
    for (int j = 0; j < TFHEpp::lvl1param::n; j++)
        pmu1[j] = binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;
    for (int j = 0; j < TFHEpp::lvl1param::n; j++)
        pmu0[j] = binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;
    TFHEpp::TRLWE<TFHEpp::lvl1param> c0 =
        TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(pmu0, TFHEpp::lvl1param::α,
                                                   sk->key.lvl1);
    TFHEpp::TRLWE<TFHEpp::lvl1param> c1 =
        TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(pmu1, TFHEpp::lvl1param::α,
                                                   sk->key.lvl1);
    const TFHEpp::Polynomial<TFHEpp::lvl1param> plainpoly = {binary(engine)};
    TFHEpp::TRGSWFFT<TFHEpp::lvl1param> cs =
        TFHEpp::trgswfftSymEncrypt<TFHEpp::lvl1param>(
            plainpoly, TFHEpp::lvl1param::α, sk->key.lvl1);
    TFHEpp::TRLWE<TFHEpp::lvl1param> res;
    for (auto _ : state) TFHEpp::CMUXFFT<TFHEpp::lvl1param>(res, cs, c1, c0);
}

void BM_ExternalProduct(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    std::array<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::n> pmu1, pmu0;
    for (int j = 0; j < TFHEpp::lvl1param::n; j++)
        pmu1[j] = binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;
    TFHEpp::TRLWE<TFHEpp::lvl1param> c0 =
        TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(pmu0, TFHEpp::lvl1param::α,
                                                   sk->key.lvl1);
    const TFHEpp::Polynomial<TFHEpp::lvl1param> plainpoly = {binary(engine)};
    TFHEpp::TRGSWFFT<TFHEpp::lvl1param> cs =
        TFHEpp::trgswfftSymEncrypt<TFHEpp::lvl1param>(
            plainpoly, TFHEpp::lvl1param::α, sk->key.lvl1);
    TFHEpp::TRLWE<TFHEpp::lvl1param> res;
    for (auto _ : state)
        TFHEpp::trgswfftExternalProduct<TFHEpp::lvl1param>(res, c0, cs);
}

void BM_CB(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    TFHEpp::EvalKey ek;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl21param;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);
    TFHEpp::TLWE<TFHEpp::lvl1param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine), TFHEpp::lvl1param::α, sk->key.lvl1);
    TFHEpp::TRGSWFFT<TFHEpp::lvl1param> res;
    for (auto _ : state)
        TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(res, ca, ek);
}

BENCHMARK(BM_TRGSWenc)
    ->Iterations(1)
    ->Repetitions(100)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_HomGate)
    ->Iterations(1)
    ->Repetitions(100)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_HomMUX)->Iterations(1)->Repetitions(10)->DisplayAggregatesOnly(
    true);
BENCHMARK(BM_TLWE2TRLWE)
    ->Iterations(1)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_IKS)->Iterations(1)->Repetitions(10)->DisplayAggregatesOnly(true);
BENCHMARK(BM_SEI)->Iterations(1)->Repetitions(10)->DisplayAggregatesOnly(true);
BENCHMARK(BM_CMUX)->Iterations(1)->Repetitions(10)->DisplayAggregatesOnly(true);
BENCHMARK(BM_ExternalProduct)
    ->Iterations(1)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true);
// BENCHMARK(BM_CB)->Iterations(1)->Repetitions(10)->DisplayAggregatesOnly(true);
BENCHMARK_MAIN();
