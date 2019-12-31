#pragma once

#include <array>
#include <cstdint>
#include <key.hpp>
#include <params.hpp>
#include <vector>

namespace TFHEpp {
using namespace std;

TLWElvl0 tlweSymEncryptlvl0(const uint32_t p, const double α,
                            const Keylvl0 &key);
bool tlweSymDecryptlvl0(const TLWElvl0 &c, const Keylvl0 &key);
vector<TLWElvl0> bootsSymEncrypt(const vector<bool> &p, const SecretKey &sk);
vector<bool> bootsSymDecrypt(const vector<TLWElvl0> &c, const SecretKey &sk);
}  // namespace TFHEpp