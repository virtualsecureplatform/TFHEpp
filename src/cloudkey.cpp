#include "cloudkey.hpp"

namespace TFHEpp {
template <class P>
void bkgen(BootStrappingKey<P> &bkfft, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++)
        bkfft[i] = trgswSymEncrypt<typename P::targetP>(
            static_cast<typename make_signed<typename P::targetP::T>::type>(
                sk.key.get<typename P::domainP>()[i]),
            P::targetP::α, sk.key.get<typename P::targetP>());
}

template <class P>
inline void bkfftgen(BootStrappingKeyFFT<P> &bkfft, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++)
        bkfft[i] = trgswfftSymEncrypt<typename P::targetP>(
            static_cast<typename make_signed<typename P::targetP::T>::type>(
                sk.key.get<typename P::domainP>()[i]),
            P::targetP::α, sk.key.get<typename P::targetP>());
}

template <class P>
inline void ikskgen(KeySwitchingKey<P> &ksk, const SecretKey &sk)
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
TFHEPP_EXPLICIT_INSTANTIATION_LVL0221_0222(INST)
#undef INST

}  // namespace TFHEpp
