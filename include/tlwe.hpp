#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include <key.hpp>
#include <params.hpp>

namespace TFHEpp {
using namespace std;

array<uint32_t, DEF_n + 1> tlweSymEncryptlvl1(
    const uint32_t p, const double α, const array<uint32_t, DEF_n> &key);
bool tlweSymDecryptlvl1(const array<uint32_t, DEF_n + 1> &c,
                        const array<uint32_t, DEF_n> &key);
vector<array<uint32_t, DEF_n + 1>> bootsSymEncrypt(const vector<bool> &p,
                                                   const SecretKey &sk);
vector<bool> bootsSymDecrypt(const vector<array<uint32_t, DEF_n + 1>> &c,
                             const SecretKey &sk);
}  // namespace TFHEpp