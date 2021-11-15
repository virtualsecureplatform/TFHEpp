#include <cassert>
#include <chrono>
#include <iostream>
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

    ca = bootsSymEncrypt(pa, *sk);
    cb = bootsSymEncrypt(pb, *sk);
    cc = bootsSymEncrypt(pc, *sk);

    chrono::system_clock::time_point start, end;

    start = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        TLWE<lvl1param> hcarry, hsum, csum;
        HomXOR(hsum, ca[test], cb[test], ek);
        HomAND(hcarry, ca[test], cb[test], ek);
        HomXOR(sum[test], hsum, cc[test], ek);
        HomAND(csum, hsum, cc[test], ek);
        HomOR(carry[test], csum, hcarry, ek);
    }

    end = chrono::system_clock::now();
    anssum = bootsSymDecrypt(sum, *sk);
    anscarry = bootsSymDecrypt(carry, *sk);
    for (int i = 0; i < num_test; i++)
        assert(anssum[i] == ((pa[i] + pb[i] + pc[i]) & 1));
    for (int i = 0; i < num_test; i++)
        assert(anscarry[i] == (((pa[i] + pb[i] + pc[i]) & 2) >> 1));
    cout << "Passed Normal" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;

    // start = chrono::system_clock::now();

    // for (int test = 0; test < num_test; test++) {
    //     TLWE<lvl0param> hcarry, hsum, csum;
    //     HomHalfAdder(hcarry, hsum, ca[test], cb[test], *gk);
    //     HomHalfAdder(csum, sum[test], hsum, cc[test], *gk);
    //     HomOR(carry[test], csum, hcarry, *gk);
    // }

    // end = chrono::system_clock::now();
    // anssum = bootsSymDecrypt(sum, *sk);
    // anscarry = bootsSymDecrypt(carry, *sk);
    // for (int i = 0; i < num_test; i++)
    //     assert(anssum[i] == ((pa[i] + pb[i] + pc[i]) & 1));
    // for (int i = 0; i < num_test; i++)
    //     assert(anscarry[i] == (((pa[i] + pb[i] + pc[i]) & 2) >> 1));
    // cout << "Passed 3in1" << endl;
    // elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
    // start)
    //               .count();
    // cout << elapsed / num_test << "ms" << endl;
}
