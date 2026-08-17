#pragma once
#include "../key.hpp"

namespace TFHEpp {
#define INST(P)                                          \
    extern template Key<P> lweKey::getSubset<P>() const; \
    extern template Key<P> lweKey::getIndependent<P>() const;
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST
}  // namespace TFHEpp
