#include <array>
#include <mulfft.hpp>
#include <params.hpp>
#include <random>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;
static random_device engine;

TRLWElvl1 trlweSymEncryptZerolvl1(const double α,
                                                         const Keylvl1 &key)
{
    uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);
<<<<<<< HEAD
    TRLWElvl1 c;
=======
    array<array<uint32_t, DEF_N>, 2> c;
>>>>>>> master
    for (uint32_t &i : c[0]) i = Torus32dist(engine);
    PolyMullvl1(c[1], c[0], key);
    for (uint32_t &i : c[1]) i += gaussian32(0, α);
    return c;
}

TRLWElvl2 trlweSymEncryptZerolvl2(const double α,
                                                         const Keylvl2 &key)
{
    uniform_int_distribution<uint64_t> Torus64dist(0, UINT64_MAX);
    TRLWElvl2 c;
    for (uint64_t &i : c[0]) i = Torus64dist(engine);
    PolyMullvl2(c[1], c[0], key);
    for (uint64_t &i : c[1]) i += gaussian64(0, α);
    return c;
}

TRLWElvl1 trlweSymEncryptlvl1(
    const array<uint32_t, DEF_N> &p, const double α, const Keylvl1 &key)
{
    TRLWElvl1 c;
    c = trlweSymEncryptZerolvl1(α, key);
    for (int i = 0; i < DEF_N; i++) c[1][i] += p[i];
    return c;
}

TRLWElvl2 trlweSymEncryptlvl2(
    const array<uint64_t, DEF_nbar> &p, const double α, const Keylvl2 &key)
{
    TRLWElvl2 c;
    c = trlweSymEncryptZerolvl2(α, key);
    for (int i = 0; i < DEF_nbar; i++) c[1][i] += p[i];
    return c;
}

array<bool, DEF_N> trlweSymDecryptlvl1(const TRLWElvl1 &c, const Keylvl1 &key)
{
    Polynomiallvl1 mulres;
    PolyMullvl1(mulres, c[0], key);
    Polynomiallvl1 phase = c[1];
    for (int i = 0; i < DEF_N; i++) phase[i] -= mulres[i];

    array<bool, DEF_N> p;
    for (int i = 0; i < DEF_N; i++) p[i] = static_cast<int32_t>(phase[i]) > 0;
    return p;
}

template <typename T = uint32_t, uint32_t N = DEF_N>
inline void SampleExtractIndex(array<T, N + 1> &tlwe,
                               const array<array<T, N>, 2> &trlwe,
                               const int index)
{
    for (int i = 0; i <= index; i++) tlwe[i] = trlwe[0][index - i];
    for (int i = index + 1; i < N; i++) tlwe[i] = -trlwe[0][N + index - i];
    tlwe[N] = trlwe[1][index];
}

void SampleExtractIndexlvl1(TLWElvl1 &tlwe, const TRLWElvl1 &trlwe,
                            const int index)
{
    SampleExtractIndex<uint32_t, DEF_N>(tlwe, trlwe, index);
}

void SampleExtractIndexlvl2(TLWElvl2 &tlwe, const TRLWElvl2 &trlwe,
                            const int index)
{
    SampleExtractIndex<uint64_t, DEF_nbar>(tlwe, trlwe, index);
}
}  // namespace TFHEpp