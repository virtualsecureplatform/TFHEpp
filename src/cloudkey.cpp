#include "cloudkey.hpp"

namespace TFHEpp {

template <class P>
void bkgen(BootstrappingKey<P> &bk, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bk[i] = trgswSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P) \
    template void bkgen<P>(BootstrappingKey<P> & bk, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void bkfftgen(BootstrappingKeyFFT<P> &bkfft, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bkfft[i] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P)                                               \
    template void bkfftgen<P>(BootstrappingKeyFFT<P> & bkfft, \
                              const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void bknttgen(BootstrappingKeyNTT<P> &bkntt, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++) {
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bkntt[i] = trgswnttSymEncrypt<typename P::targetP>(
            plainpoly, P::targetP::α, sk.key.get<typename P::targetP>());
    }
}
#define INST(P)                                               \
    template void bknttgen<P>(BootstrappingKeyNTT<P> & bkntt, \
                              const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

template <class P>
void ikskgen(KeySwitchingKey<P> &ksk, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (uint32_t k = 0; k < (1 << P::basebit) - 1; k++)
                ksk[i][j][k] = tlweSymEncrypt<typename P::targetP>(
                    sk.key.get<typename P::domainP>()[i] * (k + 1) *
                        (1ULL
                         << (numeric_limits<typename P::targetP::T>::digits -
                             (j + 1) * P::basebit)),
                    P::α, sk.key.get<typename P::targetP>());
}
#define INST(P) \
    template void ikskgen<P>(KeySwitchingKey<P> & ksk, const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH(INST)
#undef INST

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
    for (int i = 0; i < privksP::targetP::n; i++)
        poly[i] = -sk.key.get<typename privksP::targetP>()[i];
    privkskgen<privksP>(privksk[0], poly, sk);
}
#define INST(bsP, privksP) \
    template CircuitKey<bsP, privksP>::CircuitKey(const SecretKey &sk)
TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(INST)
#undef INST

}  // namespace TFHEpp
