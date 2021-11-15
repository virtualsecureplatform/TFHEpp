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
    vector<uint8_t> pd(num_test);
    vector<uint8_t> pres(num_test);
    for (int i = 0; i < num_test; i++) pa[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pb[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pc[i] = binary(engine) > 0;
    for (int i = 0; i < num_test; i++) pd[i] = binary(engine) > 0;
    vector<TLWE<lvl1param>> ca(num_test);
    vector<TLWE<lvl1param>> cb(num_test);
    vector<TLWE<lvl1param>> cc(num_test);
    vector<TLWE<lvl1param>> cd(num_test);
    vector<TLWE<lvl1param>> cres(num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] =
            tlweSymEncrypt<lvl1param>(pa[i] ? DEF_oneover12 : -(DEF_oneover12),
                                      lvl1param::α, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        cb[i] =
            tlweSymEncrypt<lvl1param>(pb[i] ? DEF_oneover12 : -(DEF_oneover12),
                                      lvl1param::α, sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        cc[i] =
            tlweSymEncrypt<lvl1param>(pc[i] ? DEF_oneover12 : -(DEF_oneover12),
                                      lvl1param::α, sk->key.lvl1);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        HomOAI3(cres[test], ca[test], cb[test], cc[test], ek);
    }

    end = chrono::system_clock::now();
    for (int i = 0; i < num_test; i++)
        pres[i] = tlweSymDecrypt<lvl1param>(cres[i], sk->key.lvl1);
    for (int i = 0; i < num_test; i++)
        assert(pres[i] == (!((pa[i] | pb[i]) & pc[i])));
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}
