#include <key.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;

lweKey::lweKey()
{
    uniform_int_distribution<int32_t> lvl0gen(lvl0param::key_value_min,
                                              lvl0param::key_value_max);
    uniform_int_distribution<int32_t> lvl1gen(lvl1param::key_value_min,
                                              lvl1param::key_value_max);
    uniform_int_distribution<int32_t> lvl2gen(lvl2param::key_value_min,
                                              lvl2param::key_value_max);
    for (typename lvl0param::T &i : lvl0) i = lvl0gen(generator);
    for (typename lvl1param::T &i : lvl1) i = lvl1gen(generator);
    for (typename lvl2param::T &i : lvl2) i = lvl2gen(generator);
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
    else if constexpr (std::is_same_v<P, lvlMparam>)
        return lvl1;
    else
        static_assert(false_v<P>, "Undefined Secret Key!");
}
#define INST(P) template Key<P> lweKey::get<P>() const;
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

}  // namespace TFHEpp
