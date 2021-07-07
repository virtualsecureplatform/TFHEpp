#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "./key.hpp"
#include "./params.hpp"
#include "./utils.hpp"

namespace TFHEpp {
using namespace std;

template <class P>
array<typename P::T, P::n + 1> tlweSymEncrypt(
    const typename P::T p, const double α,
    const array<typename P::T, P::n> &key);

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key);

template <class P = lvl0param>
vector<TLWE<P>> bootsSymEncrypt(const vector<uint8_t> &p,
                                        const SecretKey &sk);
template <class P = lvl0param>
vector<uint8_t> bootsSymDecrypt(const vector<TLWE<P>> &c,
                                const SecretKey &sk);
}  // namespace TFHEpp