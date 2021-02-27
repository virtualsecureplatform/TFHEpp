#pragma once

#include <array>
#include <cstdint>
#include <key.hpp>
#include <params.hpp>
#include <vector>

namespace TFHEpp {
using namespace std;

TLWE<lvl0param> tlweSymEncryptlvl0(const lvl0param::T p, const double α,
                                   const Key<lvl0param> &key);
bool tlweSymDecryptlvl0(const TLWE<lvl0param> &c, const Key<lvl0param> &key);
vector<TLWE<lvl0param>> bootsSymEncrypt(const vector<uint8_t> &p,
                                        const SecretKey &sk);
vector<uint8_t> bootsSymDecrypt(const vector<TLWE<lvl0param>> &c,
                                const SecretKey &sk);
}  // namespace TFHEpp