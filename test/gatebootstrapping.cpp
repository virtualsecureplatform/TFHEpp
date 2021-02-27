#include <cassert>
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

    SecretKey sk;
    GateKey* gk = new GateKey(sk);
    for (int test = 0; test < num_test; test++) {
        bool p = binary(engine) > 0;
        TLWE<lvl0param> tlwe = tlweSymEncryptlvl0(
            p ? lvl0param::μ : -lvl0param::μ, lvl0param::α, sk.key.lvl0);
        TLWE<lvl0param> bootedtlwe;
        GateBootstrapping(bootedtlwe, tlwe, *gk);
        bool p2 = tlweSymDecryptlvl0(bootedtlwe, sk.key.lvl0);
        assert(p == p2);
    }
    cout << "Passed" << endl;
}