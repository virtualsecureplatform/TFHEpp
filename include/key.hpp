#pragma once

#include<array>
#include<random>
#include<algorithm>

#include<params.hpp>

namespace TFHEpp{
    using namespace std;
    struct lweKey{
        array<uint32_t,DEF_n> lvl0;
        array<uint32_t,DEF_N> lvl1;
        array<uint64_t,DEF_N> lvl2;
        lweKey(){
            random_device engine;
            uniform_int_distribution<uint32_t> binary(0,1);
            for(uint32_t &i:lvl0) i = binary(engine);
            for(uint32_t &i:lvl1) i = binary(engine);
            for(uint64_t &i:lvl2) i = binary(engine);
        }
    };
    

    struct SecretKey{
        lweKey key;
        lweParams params;
        SecretKey(){
            lweKey lwekey;
            key = lwekey;
            lweParams lweparams;
            params = lweparams;
        }
    };
}