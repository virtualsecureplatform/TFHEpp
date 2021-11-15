#include <cassert>
#include <chrono>
#include <iostream>
#include <limits>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    const uint32_t num_test = 1000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey* sk = new SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    vector<uint8_t> pa(num_test);
    vector<uint8_t> pb(num_test);
    vector<uint8_t> pc(num_test);
    vector<uint8_t> anssum(num_test);
    vector<uint8_t> anscarry(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pb[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pc[i] = binary(engine) > 0;
    vector<TLWE<lvl1param>> ca(num_test);
    vector<TLWE<lvl1param>> cb(num_test);
    vector<TLWE<lvl1param>> cc(num_test);
    vector<TLWE<lvl1param>> carry(num_test);
    vector<TLWE<lvl1param>> sum(num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = tlweSymEncrypt<lvl1param>(
            pa[i]
                ? (1ULL << std::numeric_limits<typename lvl1param::T>::digits) /
                      12
                : -(
                      (1ULL
                       << std::numeric_limits<typename lvl1param::T>::digits) /
                      12),
            lvl1param::α, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        cb[i] = tlweSymEncrypt<lvl1param>(
            pb[i]
                ? (1ULL << std::numeric_limits<typename lvl1param::T>::digits) /
                      12
                : -(
                      (1ULL
                       << std::numeric_limits<typename lvl1param::T>::digits) /
                      12),
            lvl1param::α, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        cc[i] = tlweSymEncrypt<lvl1param>(
            pc[i]
                ? (1ULL << std::numeric_limits<typename lvl1param::T>::digits) /
                      12
                : -(
                      (1ULL
                       << std::numeric_limits<typename lvl1param::T>::digits) /
                      12),
            lvl1param::α, sk->key.lvl1);
    ;

    chrono::system_clock::time_point start, end;

    start = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        HomFullAdder(carry[test], sum[test], ca[test], cb[test], cc[test], ek);
    }

    end = chrono::system_clock::now();
    for (int i = 0; i < num_test; i++)
        anssum[i] = tlweSymDecrypt<lvl1param>(sum[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        anscarry[i] = tlweSymDecrypt<lvl1param>(carry[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        assert(anssum[i] == ((pa[i] + pb[i] + pc[i]) & 1));
    for (int i = 0; i < num_test; i++)
        assert(anscarry[i] == (((pa[i] + pb[i] + pc[i]) & 2) >> 1));
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}