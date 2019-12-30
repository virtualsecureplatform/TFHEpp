#pragma once 

#include<array>

namespace TFHEpp{
    using namespace std;
    array<array<uint32_t,DEF_N>,2>trlweSymEncryptlvl1(const array<uint32_t,DEF_N> p, const double α, const array<uint32_t,DEF_N> key);
    array<bool,DEF_N> trlweSymDecryptlvl1(array<array<uint32_t,DEF_N>,2> c, array<uint32_t,DEF_N> key);
}