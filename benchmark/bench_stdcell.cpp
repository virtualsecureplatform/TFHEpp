#include <memory>
#include <random>

#include "../include/tfhe++.hpp"
#include "google-benchmark/include/benchmark/benchmark.h"

void BM_HalfAdder(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    const std::unique_ptr<TFHEpp::GateKey> gk(new TFHEpp::GateKey(*sk));
    const TFHEpp::TLWE<TFHEpp::lvl1param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    const TFHEpp::TLWE<TFHEpp::lvl1param> cb =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    TFHEpp::TLWE<TFHEpp::lvl1param> carry, sum;
    for (auto _ : state) TFHEpp::HomHalfAdder(carry, sum, ca, cb, *gk);
}

void BM_2BRFullAdder(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    const std::unique_ptr<TFHEpp::GateKey> gk(new TFHEpp::GateKey(*sk));
    const TFHEpp::TLWE<TFHEpp::lvl1param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    const TFHEpp::TLWE<TFHEpp::lvl1param> cb =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    const TFHEpp::TLWE<TFHEpp::lvl1param> cc =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::lvl1param::μ : -TFHEpp::lvl1param::μ,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    TFHEpp::TLWE<TFHEpp::lvl1param> carry, sum;
    for (auto _ : state) TFHEpp::Hom2BRFullAdder(carry, sum, ca, cb, cc, *gk);
}

void BM_FullAdder(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    const std::unique_ptr<TFHEpp::GateKey> gk(new TFHEpp::GateKey(*sk));
    const TFHEpp::TLWE<TFHEpp::lvl1param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::DEF_oneover12 : -TFHEpp::DEF_oneover12,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    const TFHEpp::TLWE<TFHEpp::lvl1param> cb =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::DEF_oneover12 : -TFHEpp::DEF_oneover12,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    const TFHEpp::TLWE<TFHEpp::lvl1param> cc =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::DEF_oneover12 : -TFHEpp::DEF_oneover12,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    TFHEpp::TLWE<TFHEpp::lvl1param> carry, sum;
    for (auto _ : state) TFHEpp::HomFullAdder(carry, sum, ca, cb, cc, *gk);
}

void BM_AOI3(benchmark::State& state)
{
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    const std::unique_ptr<TFHEpp::GateKey> gk(new TFHEpp::GateKey(*sk));
    const TFHEpp::TLWE<TFHEpp::lvl1param> ca =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::DEF_oneover12 : -TFHEpp::DEF_oneover12,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    const TFHEpp::TLWE<TFHEpp::lvl1param> cb =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::DEF_oneover12 : -TFHEpp::DEF_oneover12,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    const TFHEpp::TLWE<TFHEpp::lvl1param> cc =
        TFHEpp::tlweSymEncrypt<TFHEpp::lvl1param>(
            binary(engine) ? TFHEpp::DEF_oneover12 : -TFHEpp::DEF_oneover12,
            TFHEpp::lvl1param::α, sk->key.lvl1);
    TFHEpp::TLWE<TFHEpp::lvl1param> res;
    for (auto _ : state) TFHEpp::HomAOI3(res, ca, cb, cc, *gk);
}

BENCHMARK(BM_HalfAdder)
    ->Iterations(1)
    ->Repetitions(100)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_2BRFullAdder)
    ->Iterations(1)
    ->Repetitions(100)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_FullAdder)
    ->Iterations(1)
    ->Repetitions(100)
    ->DisplayAggregatesOnly(true);
BENCHMARK(BM_AOI3)->Iterations(1)->Repetitions(100)->DisplayAggregatesOnly(
    true);
BENCHMARK_MAIN();