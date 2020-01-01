#include <array>
#include <cstdint>
#include <vector>
#include <limits>
#include <random>

#include <randen.h>

#include <key.hpp>
#include <params.hpp>
#include <type_traits>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;
static randen::Randen<uint64_t> engine;

template <typename T = uint32_t, uint32_t n = DEF_n>
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

TLWElvl0 tlweSymEncryptlvl0(const uint32_t p, const double α,
                            const Keylvl0 &key)
{
    return tlweSymEncrypt<uint32_t, DEF_n>(p, α, key);
}

template <typename T = uint32_t, uint32_t n = DEF_n>
bool tlweSymDecrypt(const array<T, n + 1> &c, const array<T, n> &key)
{
    T phase = c[n];
    for (int i = 0; i < n; i++) phase -= c[i] * key[i];
    bool res = static_cast<typename make_signed<T>::type>(phase) > 0;
    return res;
}

bool tlweSymDecryptlvl0(const TLWElvl0 &c, const Keylvl0 &key)
{
    return tlweSymDecrypt<uint32_t, DEF_n>(c, key);
}

vector<TLWElvl0> bootsSymEncrypt(const vector<uint8_t> &p, const SecretKey &sk)
{
    vector<TLWElvl0> c(p.size());
    for (int i = 0; i < p.size(); i++)
        c[i] = tlweSymEncryptlvl0(p[i] ? DEF_μ : -DEF_μ, DEF_α, sk.key.lvl0);
    return c;
}

vector<uint8_t> bootsSymDecrypt(const vector<TLWElvl0> &c, const SecretKey &sk)
{
    vector<uint8_t> p(c.size());
    for (int i = 0; i < p.size(); i++)
        p[i] = tlweSymDecryptlvl0(c[i], sk.key.lvl0);
    return p;
}
}  // namespace TFHEpp