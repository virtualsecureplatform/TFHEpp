#include "cloudkey.hpp"

namespace TFHEpp {

GateKeywoFFT::GateKeywoFFT(const SecretKey &sk)
{
    // Generete bkfft
    bkgen<lvl01param>(bklvl01, sk);

    // Generete ksk
    ikskgen<lvl10param>(ksk, sk);
}

GateKey::GateKey(const SecretKey &sk)
{
    // Generete bkfft
    bkfftgen<lvl01param>(bkfftlvl01, sk);

    // Generete ksk
    ikskgen<lvl10param>(ksk, sk);
}

GateKey::GateKey(const GateKeywoFFT &gkwofft)
{
    for (int i = 0; i < lvl01param::domainP::n; i++)
        bkfftlvl01[i] = ApplyFFT2trgsw<lvl1param>(gkwofft.bklvl01[i]);
    ksk = gkwofft.ksk;
}

GateKeyNTT::GateKeyNTT(const GateKeywoFFT &gkwofft)
{
    for (int i = 0; i < lvl01param::domainP::n; i++)
        bknttlvl01[i] = ApplyNTT2trgsw<lvl1param>(gkwofft.bklvl01[i]);
    ksk = gkwofft.ksk;
}

template <class bsP, class privksP>
CircuitKey<bsP, privksP>::CircuitKey(const SecretKey &sk)
{
    // Generete bkfft
    bkfftgen<bsP>(bkfft, sk);

    // Generate privksk
    TFHEpp::Polynomial<typename privksP::targetP> poly = {1};
    privkskgen<privksP>(privksk[1], poly, sk);
    for(int i = 0; i<privksP::targetP::n; i++) poly[i] = -sk.key.get<typename privksP::targetP>()[i];
    privkskgen<privksP>(privksk[0], poly, sk);
}
#define INST(bsP, privksP) \
    template CircuitKey<bsP, privksP>::CircuitKey(const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

}  // namespace TFHEpp
