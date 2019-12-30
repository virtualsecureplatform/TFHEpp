#pragma once

#include<array>
#include<vector>
#include<cstdint>

#include<params.hpp>
#include<key.hpp>

namespace TFHEpp{
    using namespace std;

    array<uint32_t,DEF_n+1> tlweSymEncryptlvl1(const uint32_t p, const double α, const array<uint32_t,DEF_n> &key);
    vector<array<uint32_t,DEF_n+1>> bootsSymEncrypt(const vector<bool> &p,const SecretKey &sk);
    vector<bool> bootsSymDecrypt(const vector<array<uint32_t,DEF_n+1>> &c, const SecretKey &sk);
}