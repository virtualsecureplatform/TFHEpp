#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

#include <gperftools/profiler.h>

using namespace std;
using namespace TFHEpp;

int main()
{
    const uint32_t num_test = 1000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey* sk = new SecretKey();
    GateKey* gk = new GateKey(*sk);
    vector<uint8_t> pa(num_test);
    vector<uint8_t> pb(num_test);
    vector<uint8_t> pres(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pb[i] = binary(engine) > 0;
    vector<TLWElvl0> ca(num_test);
    vector<TLWElvl0> cb(num_test);
    vector<TLWElvl0> cres(num_test);

    ca = bootsSymEncrypt(pa, *sk);
    cb = bootsSymEncrypt(pb, *sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    // ProfilerStart("nand.prof");
    for (int test = 0; test < num_test; test++) {
        HomNAND(cres[test], ca[test], cb[test], *gk);
    }
    // ProfilerStop();

    end = chrono::system_clock::now();
    pres = bootsSymDecrypt(cres, *sk);
    for (int i = 0; i < num_test; i++) assert(pres[i] == !(pa[i] & pb[i]));
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}