#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "key.hpp"
#include "params.hpp"

namespace TFHEpp {
using namespace std;

template <class P>
TLWE<P> tlweSymEncrypt(const typename P::T p, const double α,
                       const Key<P> &key);

template <class P>
TLWE<P> tlweSymIntEncrypt(const typename P::T p, const double α,
                          const Key<P> &key);

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key);

template <class P>
typename P::T tlweSymIntDecrypt(const TLWE<P> &c, const Key<P> &key);

template <class P = lvl1param>
vector<TLWE<P>> bootsSymEncrypt(const vector<uint8_t> &p, const SecretKey &sk);
template <class P = lvl1param>
vector<uint8_t> bootsSymDecrypt(const vector<TLWE<P>> &c, const SecretKey &sk);
}  // namespace TFHEpp