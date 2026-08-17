#include "tfhe/key.hpp"
namespace TFHEpp {
#define INST(P)                                   \
    template Key<P> lweKey::getSubset<P>() const; \
    template Key<P> lweKey::getIndependent<P>() const;
TFHEPP_EXPLICIT_INSTANTIATION_TLWE(INST)
#undef INST
}  // namespace TFHEpp
