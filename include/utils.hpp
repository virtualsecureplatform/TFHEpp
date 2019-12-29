#pragma once
#include<array>
#include<algorithm>
#include<cmath>

#include<params.hpp>

namespace TFHEpp{
    using namespace std;

    inline double ttod32(const uint32_t torus){
        return static_cast<double>(static_cast<int32_t>(torus));
    }

    template<uint32_t N = DEF_N>
    inline array<double,N> ttod32_array(const array<uint32_t,DEF_N> torus_array){
        array<double,N> res;
        transform(torus_array.cbegin(),torus_array.cend(),res.begin(),ttod32);
        return res;
    }
}