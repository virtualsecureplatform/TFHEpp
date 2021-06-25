#include <randen.h>

#include <array>
#include <cstdint>
#include <key.hpp>
#include <limits>
#include <params.hpp>
#include <random>
#include <tlwe.hpp>
#include <type_traits>
#include <utils.hpp>
#include <vector>

namespace TFHEpp {
using namespace std;

template <class P>
array<typename P::T, P::n + 1> tlweSymEncrypt(
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
#define INST(P)                                                \
    template array<typename P::T, P::n + 1> tlweSymEncrypt<P>( \
        const typename P::T p, const double α,                 \
        const array<typename P::T, P::n> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = c[P::n];
    for (int i = 0; i < P::n; i++) phase -= c[i] * key[i];
    bool res =
        static_cast<typename make_signed<typename P::T>::type>(phase) > 0;
    return res;
}
#define INST(P) \
    template bool tlweSymDecrypt<P>(const TLWE<P> &c, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

vector<TLWE<lvl0param>> bootsSymEncrypt(const vector<uint8_t> &p,
                                        const SecretKey &sk)
{
    vector<TLWE<lvl0param>> c(p.size());
    for (int i = 0; i < p.size(); i++)
        c[i] = tlweSymEncrypt<lvl0param>(p[i] ? lvl0param::μ : -lvl0param::μ,
                                         lvl0param::α, sk.key.lvl0);
    return c;
}

vector<uint8_t> bootsSymDecrypt(const vector<TLWE<lvl0param>> &c,
                                const SecretKey &sk)
{
    vector<uint8_t> p(c.size());
    for (int i = 0; i < p.size(); i++)
        p[i] = tlweSymDecrypt<lvl0param>(c[i], sk.key.lvl0);
    return p;
}
}  // namespace TFHEpp