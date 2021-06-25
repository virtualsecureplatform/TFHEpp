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
    array<typename privksP::domainP::T, privksP::domainP::n + 1> key;
    for (int i = 0; i < privksP::domainP::n; i++) key[i] = sk.key.lvl2[i];
    key[privksP::domainP::n] = -1;
#pragma omp parallel for collapse(4)
    for (int z = 0; z < 2; z++)
        for (int i = 0; i <= privksP::domainP::n; i++)
            for (int j = 0; j < privksP::t; j++)
                for (typename privksP::targetP::T u = 0;
                     u < (1 << privksP::basebit) - 1; u++) {
                    TRLWE<typename privksP::targetP> c =
                        trlweSymEncryptZero<typename privksP::targetP>(
                            privksP::α,
                            sk.key.get<typename privksP::targetP>());
                    c[z][0] += (u + 1) * key[i]
                               << (numeric_limits<
                                       typename privksP::targetP::T>::digits -
                                   (j + 1) * privksP::basebit);
                    privksk[z][i][j][u] = c;
                }
}
#define INST(bsP, privksP) \
    template CircuitKey<bsP, privksP>::CircuitKey(const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

}  // namespace TFHEpp
