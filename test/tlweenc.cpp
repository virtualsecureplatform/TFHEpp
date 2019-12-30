#include<tfhe++.hpp>
#include<cassert>
#include<vector>
#include<random>

using namespace TFHEpp;

int main(){
    const uint32_t num_test = 10000;

    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    lweKey sk;

    vector<bool> p(num_test);
    for(bool i:p) i = binary(engine)>0;
    vector<array<uint32_t,DEF_n+1>> c(num_test);
    c = bootsSymEncrypt(p,sk);
    vector<bool> p2(num_test);
    p2 = bootsSymDecrypt(c,sk);
    for(int i = 0; i<num_test;i++) assert(p[i] == p2[i]);
}