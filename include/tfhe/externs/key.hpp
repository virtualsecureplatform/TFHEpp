#pragma once
#include "../key.hpp"

namespace TFHEpp{
#define INST(P) extern template Key<P> lweKey::get<P>() const;
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST
}