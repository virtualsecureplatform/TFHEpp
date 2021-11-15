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
    vector<uint8_t> pxor(num_test);
    vector<uint8_t> pnand(num_test);
    vector<uint8_t> pnor(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pb[i] = binary(engine) > 0;
    vector<TLWE<lvl1param>> ca(num_test);
    vector<TLWE<lvl1param>> cb(num_test);
    vector<TLWE<lvl1param>> cxor(num_test);
    vector<TLWE<lvl1param>> cnand(num_test);
    vector<TLWE<lvl1param>> cnor(num_test);

    ca = bootsSymEncrypt(pa, *sk);
    cb = bootsSymEncrypt(pb, *sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        HomXORNANDNOR(cxor[test], cnand[test], cnor[test], ca[test],
                           cb[test], ek);
    }

    end = chrono::system_clock::now();
    for (int i = 0; i < num_test; i++)
        pxor[i] = tlweSymDecrypt<lvl1param>(cxor[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        pnand[i] = tlweSymDecrypt<lvl1param>(cnand[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        pnor[i] = tlweSymDecrypt<lvl1param>(cnor[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++) assert(pxor[i] == (pa[i] ^ pb[i]));
    for (int i = 0; i < num_test; i++) assert(pnand[i] == !(pa[i] & pb[i]));
    for (int i = 0; i < num_test; i++) assert(pnor[i] == !(pa[i] | pb[i]));
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}
