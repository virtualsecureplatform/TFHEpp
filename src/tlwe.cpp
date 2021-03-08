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