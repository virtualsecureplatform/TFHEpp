#include "key.hpp"
namespace TFHEpp {
#define INST(P) template Key<P> lweKey::get<P>() const;
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST
}  // namespace TFHEpp