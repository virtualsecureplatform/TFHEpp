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
inline array<typename P::T, P::n + 1> tlweSymEncrypt(
    const typename P::T p, const double α,
    const array<typename P::T, P::n> &key)
{
    uniform_int_distribution<typename P::T> Torusdist(
        0, numeric_limits<typename P::T>::max());
    array<typename P::T, P::n + 1> res = {};
    res[P::n] = ModularGaussian<P>(p, α);
    for (int i = 0; i < P::n; i++) {
        res[i] = Torusdist(generator);
        res[P::n] += res[i] * key[i];
    }
    return res;
}

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = c[P::n];
    for (int i = 0; i < P::n; i++) phase -= c[i] * key[i];
    bool res =
        static_cast<typename make_signed<typename P::T>::type>(phase) > 0;
    return res;
}

vector<TLWE<lvl0param>> bootsSymEncrypt(const vector<uint8_t> &p,
                                        const SecretKey &sk);
vector<uint8_t> bootsSymDecrypt(const vector<TLWE<lvl0param>> &c,
                                const SecretKey &sk);
}  // namespace TFHEpp