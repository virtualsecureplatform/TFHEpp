#pragma once

#include<array>
#include<vector>
#include<cstdint>

#include<params.hpp>
#include<key.hpp>

namespace TFHEpp{
    using namespace std;

    vector<array<uint32_t,DEF_n+1>> bootsSymEncrypt(const vector<bool> p,const lweKey sk);
    vector<bool> bootsSymDecrypt(const vector<array<uint32_t,DEF_n+1>> c, const lweKey sk);
}