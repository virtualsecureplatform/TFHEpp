#include <cloudkey.hpp>
#include <detwfa.hpp>
#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>
#include <mulfft.hpp>
#include <params.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>
#include <utils.hpp>

namespace TFHEpp {

#define INST(P)                                     \
    template void GateBootstrappingTLWE2TLWEFFT<P>( \
        TLWE<typename P::targetP> & res,            \
        const TLWE<typename P::domainP> &tlwe,      \
        const BootstrappingKeyFFT<P> &bkfft,        \
        const Polynomial<typename P::targetP> &testvector)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

}  // namespace TFHEpp
