#include <randen.h>

#include <array>
#include <cstdint>
#include <key.hpp>
#include <limits>
#include <params.hpp>
#include <random>
#include <type_traits>
#include <utils.hpp>
#include <vector>

namespace TFHEpp {
using namespace std;
static randen::Randen<uint64_t> engine;

template <typename T, uint32_t n>
inline array<T, n + 1> tlweSymEncrypt(const T p, const double α,
                                      const array<uint32_t, n> &key)
{
    uniform_int_distribution<T> Torusdist(0, numeric_limits<T>::max());
    array<uint32_t, n + 1> res = {};
    res[n] = gaussian32(p, α);
    for (int i = 0; i < n; i++) {
        res[i] = Torusdist(engine);
        res[n] += res[i] * key[i];
    }
    return res;
}

TLWE<lvl0param> tlweSymEncryptlvl0(const lvl0param::T p, const double α,
                            const Key<lvl0param> &key)
{
    return tlweSymEncrypt<lvl0param::T, lvl0param::n>(p, α, key);
}

template <typename T, uint32_t n>
bool tlweSymDecrypt(const array<T, n + 1> &c, const array<T, n> &key)
{
    T phase = c[n];
    for (int i = 0; i < n; i++) phase -= c[i] * key[i];
    bool res = static_cast<typename make_signed<T>::type>(phase) > 0;
    return res;
}

bool tlweSymDecryptlvl0(const TLWE<lvl0param> &c, const Key<lvl0param> &key)
{
    return tlweSymDecrypt<lvl0param::T, lvl0param::n>(c, key);
}

vector<TLWE<lvl0param>> bootsSymEncrypt(const vector<uint8_t> &p, const SecretKey &sk)
{
    vector<TLWE<lvl0param>> c(p.size());
    for (int i = 0; i < p.size(); i++)
        c[i] = tlweSymEncryptlvl0(p[i] ? lvl0param::μ : -lvl0param::μ, lvl0param::α, sk.key.lvl0);
    return c;
}

vector<uint8_t> bootsSymDecrypt(const vector<TLWE<lvl0param>> &c, const SecretKey &sk)
{
    vector<uint8_t> p(c.size());
    for (int i = 0; i < p.size(); i++)
        p[i] = tlweSymDecryptlvl0(c[i], sk.key.lvl0);
    return p;
}
}  // namespace TFHEpp