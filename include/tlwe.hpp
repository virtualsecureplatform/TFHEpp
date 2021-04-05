#pragma once

#include <array>
#include <cstdint>
#include <key.hpp>
#include <params.hpp>
#include <utils.hpp>
#include <vector>

namespace TFHEpp {
using namespace std;

template <class P>
array<typename P::T, P::n + 1> tlweSymEncrypt(
    const typename P::T p, const double α,
    const array<typename P::T, P::n> &key);

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key);

vector<TLWE<lvl0param>> bootsSymEncrypt(const vector<uint8_t> &p,
                                        const SecretKey &sk);
vector<uint8_t> bootsSymDecrypt(const vector<TLWE<lvl0param>> &c,
                                const SecretKey &sk);
}  // namespace TFHEpp