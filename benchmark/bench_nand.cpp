#include "google-benchmark/include/benchmark/benchmark.h"
#include "../include/tfhe++.hpp"
#include <random>

void BM_HomGate(benchmark::State& state){
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const TFHEpp::SecretKey* sk = new TFHEpp::SecretKey();
    const TFHEpp::GateKey* gk = new TFHEpp::GateKey(*sk);
    TFHEpp::TLWE<TFHEpp::lvl0param> ca = TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(binary(engine),TFHEpp::lvl0param::α,sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> cb = TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(binary(engine),TFHEpp::lvl0param::α,sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> res;
    for (auto _ : state) TFHEpp::HomNAND(res, ca, cb, *gk);
}

void BM_HomMUX(benchmark::State& state){
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
    const TFHEpp::SecretKey* sk = new TFHEpp::SecretKey();
    const TFHEpp::GateKey* gk = new TFHEpp::GateKey(*sk);
    TFHEpp::TLWE<TFHEpp::lvl0param> ca = TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(binary(engine),TFHEpp::lvl0param::α,sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> cb = TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(binary(engine),TFHEpp::lvl0param::α,sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> cs = TFHEpp::tlweSymEncrypt<TFHEpp::lvl0param>(binary(engine),TFHEpp::lvl0param::α,sk->key.lvl0);
    TFHEpp::TLWE<TFHEpp::lvl0param> res;
    for (auto _ : state) TFHEpp::HomMUX(res, cs, ca, cb, *gk);
}

BENCHMARK(BM_HomGate)->Iterations(1)->Repetitions(10);
BENCHMARK(BM_HomMUX)->Iterations(1)->Repetitions(10);
BENCHMARK_MAIN();
