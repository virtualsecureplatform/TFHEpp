#pragma once

#include<array>
#include<cstdint>

#include<params.hpp>

namespace TFHEpp{
    using namespace std;

    array<array<array<double,DEF_N>,2>,2*DEF_l> trgswfftSymEncryptlvl1(int32_t p, double α, array<uint32_t,DEF_N> &key);
    void trgswfftExternalProductlvl1(array<array<uint32_t,DEF_N>,2> &restrlwe, const array<array<array<double,DEF_N>,2>,2*DEF_l> trgswfft, const array<array<uint32_t,DEF_N>,2> &trlwe);
}