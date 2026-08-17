#include <memory>
#include <random>

#include "../include/tfhe++.hpp"
#include "google-benchmark/include/benchmark/benchmark.h"

void BM_TLWE2TRLWE(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    TFHEpp::TLWE<TFHEpp::lvl0param> ca;
    TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(
        ca, static_cast<TFHEpp::lvl0param::T>(binary(engine)),
        TFHEpp::lvl0param::α, sk->key.getSubset<TFHEpp::lvl0param>());
    TFHEpp::TRLWE<TFHEpp::lvl1param> res;
    for (auto _ : state)
        TFHEpp::BlindRotate<TFHEpp::lvl01param>(
            res, ca, ek.getbkfft<TFHEpp::lvl01param>(),
            TFHEpp::μpolygen<TFHEpp::lvl1param, TFHEpp::lvl1param::μ>());
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
    TFHEpp::TRLWE<TFHEpp::lvl1param> c0, c1;
    TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(c0, pmu0,
                                               sk->key.getSubset<TFHEpp::lvl1param>());
    TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(c1, pmu1,
                                               sk->key.getSubset<TFHEpp::lvl1param>());
    const TFHEpp::Polynomial<TFHEpp::lvl1param> plainpoly = {binary(engine)};
    TFHEpp::TRGSWFFT<TFHEpp::lvl1param> cs;
    TFHEpp::trgswSymEncrypt<TFHEpp::lvl1param>(cs, plainpoly,
                                               TFHEpp::lvl1param::α,
                                               sk->key.getSubset<TFHEpp::lvl1param>());
    TFHEpp::TRLWE<TFHEpp::lvl1param> res;
    for (auto _ : state) TFHEpp::CMUXFFT<TFHEpp::lvl1param>(res, cs, c1, c0);
}

void BM_ExternalProduct(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    std::array<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::n> pmu0;
    for (int j = 0; j < TFHEpp::lvl1param::n; j++)
        pmu0[j] = binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ;
    TFHEpp::TRLWE<TFHEpp::lvl1param> c0;
    TFHEpp::trlweSymEncrypt<TFHEpp::lvl1param>(c0, pmu0,
                                               sk->key.getSubset<TFHEpp::lvl1param>());
    const TFHEpp::Polynomial<TFHEpp::lvl1param> plainpoly = {binary(engine)};
    TFHEpp::TRGSWFFT<TFHEpp::lvl1param> cs;
    TFHEpp::trgswSymEncrypt<TFHEpp::lvl1param>(cs, plainpoly,
                                               TFHEpp::lvl1param::α,
                                               sk->key.getSubset<TFHEpp::lvl1param>());
    TFHEpp::TRLWE<TFHEpp::lvl1param> res;
    for (auto _ : state)
        TFHEpp::ExternalProduct<TFHEpp::lvl1param>(res, c0, cs);
}

BENCHMARK(BM_TLWE2TRLWE)
    ->Iterations(1)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_CMUX)->Iterations(1)->Repetitions(10)->DisplayAggregatesOnly(true);
BENCHMARK(BM_ExternalProduct)
    ->Iterations(1)
    ->Repetitions(10)
    ->DisplayAggregatesOnly(true);
BENCHMARK_MAIN();
