#include<tfhe++.hpp>
#include<cassert>
#include<random>
#include<iostream>

using namespace std;
using namespace TFHEpp;

int main(){
    const uint32_t num_test = 1000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey sk;
    CloudKey ck(sk);
    for(int test = 0;test<num_test;test++){
        bool p = binary(engine)>0;
        array<uint32_t,DEF_n+1> tlwe = tlweSymEncryptlvl1(p?DEF_MU:-DEF_MU,DEF_α,sk.key.lvl0);
        array<uint32_t,DEF_n+1> bootedtlwe;
        GateBootstrapping(bootedtlwe,tlwe,ck);
        bool p2 = tlweSymDecryptlvl1(bootedtlwe,sk.key.lvl0);
        assert(p == p2);
    }
    cout<<"Passed"<<endl;
}