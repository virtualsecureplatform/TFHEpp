#pragma once

#include "../thirdparties/cereal/include/cereal/archives/portable_binary.hpp"
#include "../thirdparties/cereal/include/cereal/types/array.hpp"
#include "params.hpp"
#include "tlwe.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"

namespace TFHEpp {

template <class P>
void bkgen(BootstrappingKey<P> &bkfft, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++)
        bkfft[i] = trgswSymEncrypt<typename P::targetP>(
            static_cast<typename make_signed<typename P::targetP::T>::type>(
                sk.key.get<typename P::domainP>()[i]),
            P::targetP::α, sk.key.get<typename P::targetP>());
}

template <class P>
inline void bkfftgen(BootstrappingKeyFFT<P> &bkfft, const SecretKey &sk)
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

struct GateKeywoFFT {
    BootstrappingKey<lvl01param> bklvl01;
    KeySwitchingKey<lvl10param> ksk;
    GateKeywoFFT(const SecretKey &sk);
    GateKeywoFFT() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(ksk, bklvl01);
    }
};

struct GateKey {
    BootstrappingKeyFFT<lvl01param> bkfftlvl01;
    KeySwitchingKey<lvl10param> ksk;
    GateKey(const SecretKey &sk);
    GateKey(const GateKeywoFFT &gkwofft);
    GateKey() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(ksk, bkfftlvl01);
    }
};

template <class bsP, class privksP>
struct CircuitKey {
    BootstrappingKeyFFT<bsP> bkfft;
    std::array<PrivKeySwitchKey<privksP>, 2> privksk;
    CircuitKey(const SecretKey &sk);
    CircuitKey() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(privksk, bkfft);
    }
};

template <class CBbsP, class CBprivksP, class CMksP>
struct CloudKey {
    GateKey gk;
    CircuitKey<CBbsP, CBprivksP> ck;
    KeySwitchingKey<CMksP> ksk;
    lweParams params;
    CloudKey(SecretKey sk) : gk(sk), ck(sk) { ikskgen<CMksP>(ksk, sk); }
    CloudKey() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(gk.ksk, gk.bkfftlvl01, ck.privksk, ck.bkfft, params);
    }
};
}  // namespace TFHEpp