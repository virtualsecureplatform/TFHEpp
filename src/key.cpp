#include <key.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;

lweKey::lweKey()
{
    uniform_int_distribution<int32_t> lvl0gen(lvl0param::key_value_min,
                                              lvl0param::key_value_max);
    uniform_int_distribution<uint32_t> binary(0, 1);
    for (typename lvl0param::T &i : lvl0) i = lvl0gen(generator);
    // for (typename lvl0param::T &i : lvl0) i = binary(generator);
    for (typename lvl1param::T &i : lvl1) i = binary(generator);
    for (typename lvl2param::T &i : lvl2) i = binary(generator);
}

template <class P>
Key<P> lweKey::get() const
{
    if constexpr (std::is_same_v<P, lvl0param>)
        return lvl0;
    else if constexpr (std::is_same_v<P, lvl1param>)
        return lvl1;
    else if constexpr (std::is_same_v<P, lvl2param>)
        return lvl2;
}
#define INST(P) template Key<P> lweKey::get<P>() const;
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

}  // namespace TFHEpp
