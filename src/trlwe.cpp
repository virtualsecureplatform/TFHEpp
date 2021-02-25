#include <randen.h>

#include <array>
#include <mulfft.hpp>
#include <params.hpp>
#include <random>
#include <utils.hpp>
#include <trlwe.hpp>

namespace TFHEpp {
using namespace std;

TRLWE<lvl1param> trlweSymEncryptZerolvl1(const double α, const Key<lvl1param> &key)
{
    return trlweSymEncryptZero<lvl1param>(α,key);
}

TRLWE<lvl2param> trlweSymEncryptZerolvl2(const double α, const Key<lvl2param> &key)
{
    return trlweSymEncryptZero<lvl2param>(α,key);
}

template<class P>
TRLWE<P> trlweSymEncrypt(const array<typename P::T, P::n> &p, const double α,
                              const Key<P> &key)
{
    TRLWE<P> c;
    c = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::n; i++) c[1][i] += p[i];
    return c;
}

TRLWE<lvl1param> trlweSymEncryptlvl1(const array<lvl1param::T, lvl1param::n> &p, const double α,
                              const Key<lvl1param> &key)
{
    return trlweSymEncrypt<lvl1param>(p,α,key);
}

TRLWE<lvl2param> trlweSymEncryptlvl2(const array<lvl2param::T, lvl2param::n> &p,
                              const double α, const Key<lvl2param> &key)
{
    return trlweSymEncrypt<lvl2param>(p,α,key);
}

template<class P>
array<bool, P::n> trlweSymDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> mulres;
    PolyMul<P>(mulres, c[0], key);
    Polynomial<P> phase = c[1];
    for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];

    array<bool, P::n> p;
    for (int i = 0; i < P::n; i++) p[i] = static_cast<make_signed<P::T>>(phase[i]) > 0;
    return p;
}

array<bool, lvl1param::n> trlweSymDecryptlvl1(const TRLWE<lvl1param> &c, const Key<lvl1param> &key)
{
    return trlweSymDecrypt<lvl1param>(c,key);
}

array<bool, lvl2param::n> trlweSymDecryptlvl2(const TRLWE<lvl2param> &c,
                                          const Key<lvl2param> &key)
{
    return trlweSymDecrypt<lvl2param>(c,key);
}

void SampleExtractIndexlvl1(TLWE<lvl1param> &tlwe, const TRLWE<lvl1param> &trlwe,
                            const int index)
{
    SampleExtractIndex<lvl1param>(tlwe, trlwe, index);
}

void SampleExtractIndexlvl2(TLWE<lvl2param> &tlwe, const TRLWE<lvl2param> &trlwe,
                            const int index)
{
    SampleExtractIndex<lvl2param>(tlwe, trlwe, index);
}
}  // namespace TFHEpp