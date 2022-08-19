#ifdef USE_PERF
#include <gperftools/profiler.h>
#endif

#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    std::cout << lvl1param::k << std::endl;
    constexpr uint32_t num_test = 1000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey* sk = new SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    vector<uint8_t> pa(num_test);
    vector<uint8_t> pb(num_test);
    vector<uint8_t> pres(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pb[i] = binary(engine) > 0;
    vector<TLWE<TFHEpp::lvl1param>> ca(num_test);
    vector<TLWE<TFHEpp::lvl1param>> cb(num_test);
    vector<TLWE<TFHEpp::lvl1param>> cres(num_test);

    ca = bootsSymEncrypt(pa, *sk);
    cb = bootsSymEncrypt(pb, *sk);

    chrono::system_clock::time_point start, end;
#ifdef USE_PERF
    ProfilerStart("nand.prof");
#endif
    start = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        HomNAND(cres[test], ca[test], cb[test], ek);
    }

    end = chrono::system_clock::now();
#ifdef USE_PERF
    ProfilerStop();
#endif
    pres = bootsSymDecrypt(cres, *sk);
    for (int i = 0; i < num_test; i++) assert(pres[i] == !(pa[i] & pb[i]));
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}