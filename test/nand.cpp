#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    const uint32_t num_test = 10;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey sk;
    CloudKey ck(sk);
    vector<bool> pa(num_test);
    vector<bool> pb(num_test);
    vector<bool> pres(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pb[i] = binary(engine) > 0;
    vector<array<uint32_t, DEF_n + 1>> ca(num_test);
    vector<array<uint32_t, DEF_n + 1>> cb(num_test);
    vector<array<uint32_t, DEF_n + 1>> cres(num_test);

    ca = bootsSymEncrypt(pa, sk);
    cb = bootsSymEncrypt(pb, sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        HomNAND(cres[test], ca[test], cb[test], ck);
    }
    end = chrono::system_clock::now();
    pres = bootsSymDecrypt(cres, sk);
    for (int i = 0; i < num_test; i++) assert(pres[i] == !(pa[i] & pb[i]));
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}