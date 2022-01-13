#include "trlwe.hpp"

#include "mulfft.hpp"

namespace TFHEpp {
using namespace std;

template <class P>
TRLWE<P> trlweSymEncryptZero(const double α, const Key<P> &key)
{
    uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
    TRLWE<P> c;
    for (typename P::T &i : c[P::k]) i = ModularGaussian<P>(0, α);
    for (int k = 0; k < P::k; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator);
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < P::n; i++) c[P::k][i] += temp[i];
    }
    return c;
}
#define INST(P) \
    template TRLWE<P> trlweSymEncryptZero<P>(const double α, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
TRLWE<P> trlweSymEncrypt(const array<typename P::T, P::n> &p, const double α,
                         const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
    return c;
}
#define INST(P)                                                               \
    template TRLWE<P> trlweSymEncrypt<P>(const array<typename P::T, P::n> &p, \
                                         const double α, const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
TRLWE<P> trlweSymIntEncrypt(const array<typename P::T, P::n> &p, const double α,
                            const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::Δ * p[i]);
    return c;
}
#define INST(P)                                              \
    template TRLWE<P> trlweSymIntEncrypt<P>(                 \
        const array<typename P::T, P::n> &p, const double α, \
        const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
array<bool, P::n> trlweSymDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }

    array<bool, P::n> p;
    for (int i = 0; i < P::n; i++)
        p[i] = static_cast<typename make_signed<typename P::T>::type>(
                   phase[i]) > 0;
    return p;
}
#define INST(P)                                                      \
    template array<bool, P::n> trlweSymDecrypt<P>(const TRLWE<P> &c, \
                                                  const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
Polynomial<P> trlweSymIntDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }

    Polynomial<P> p;
    for (int i = 0; i < P::n; i++)
        p[i] = static_cast<typename P::T>(std::round(phase[i] / P::Δ)) %
               P::plain_modulus;
    return p;
}
#define INST(P)                                                     \
    template Polynomial<P> trlweSymIntDecrypt<P>(const TRLWE<P> &c, \
                                                 const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

template <class P>
void SampleExtractIndex(TLWE<P> &tlwe, const TRLWE<P> &trlwe, const int index)
{
    for (int k = 0; k < P::k; k++) {
        for (int i = 0; i <= index; i++)
            tlwe[k * P::n + i] = trlwe[k][index - i];
        for (int i = index + 1; i < P::n; i++)
            tlwe[k * P::n + i] = -trlwe[k][P::n + index - i];
    }
    tlwe[P::k * P::n] = trlwe[P::k][index];
}
#define INST(P)                                                                \
    template void SampleExtractIndex<P>(TLWE<P> & tlwe, const TRLWE<P> &trlwe, \
                                        const int index)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST
}  // namespace TFHEpp
