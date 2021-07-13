#include <key.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;

lweKey::lweKey()
{
    randen::Randen<uint64_t> engine;
    uniform_int_distribution<uint32_t> binary(0, 1);
    for (typename lvl0param::T &i : lvl0) i = binary(engine);
    for (typename lvl1param::T &i : lvl1) i = binary(engine);
    for (typename lvl2param::T &i : lvl2) i = binary(engine);
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
    else
        static_assert(false_v<P>, "Undefined Secret Key!");
}
#define INST(P) template Key<P> lweKey::get<P>() const;
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST

}  // namespace TFHEpp
