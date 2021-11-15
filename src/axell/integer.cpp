#include <axell/gate.hpp>
#include <axell/integer.hpp>
#include <axell/mpparam.hpp>
#include <gate.hpp>
#include <params.hpp>

namespace TFHEpp::detail::integer {

template <>
void full_adder<lvlMparam>(TLWE<lvlMparam>& c1, TLWE<lvlMparam>& s,
                           const TLWE<lvlMparam>& a, const TLWE<lvlMparam>& b,
                           const TLWE<lvlMparam>& c0, const EvalKey& ek)
{
    HomFullAdder(c1, s, a, b, c0, ek);
}

template <>
void full_adder<lvl1param>(TLWE<lvl1param>& c1, TLWE<lvl1param>& s,
                           const TLWE<lvl1param>& a, const TLWE<lvl1param>& b,
                           const TLWE<lvl1param>& c0, const EvalKey& ek)
{
    TLWE<lvl1param> tmp;
    HomXOR(tmp, a, b, ek);
    HomXOR(s, tmp, c0, ek);
    HomMUX(c1, tmp, c0, a, ek);
}

}  // namespace TFHEpp::detail::integer
