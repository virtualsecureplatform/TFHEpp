#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>

#include "params.hpp"
#include "tlwe.hpp"
#include "trgsw.hpp"
#include "trlwe.hpp"

namespace TFHEpp {

template <class P>
inline void bkgen(BootstrappingKey<P> &bk, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++){
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bk[i] = trgswSymEncrypt<typename P::targetP>(
            plainpoly,
            P::targetP::α, sk.key.get<typename P::targetP>());
    }
}

template <class P>
inline void bkfftgen(BootstrappingKeyFFT<P> &bkfft, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++){
        Polynomial<typename P::targetP> plainpoly = {};
        plainpoly[0] = sk.key.get<typename P::domainP>()[i];
        bkfft[i] = trgswfftSymEncrypt<typename P::targetP>(
            plainpoly,
            P::targetP::α, sk.key.get<typename P::targetP>());
    }
}

template <class P>
inline void tlwe2trlweikskkgen(TLWE2TRLWEIKSKey<P> &iksk, const SecretKey &sk)
{
    for (int i = 0; i < P::domainP::n; i++)
        for (int j = 0; j < P::t; j++)
            for (uint32_t k = 0; k < (1 << P::basebit) - 1; k++){
                Polynomial<typename P::targetP> p = {};
                p[0] = sk.key.get<typename P::domainP>()[i] * (k + 1) *
                        (1ULL
                         << (numeric_limits<typename P::targetP::T>::digits -
                             (j + 1) * P::basebit));
                iksk[i][j][k] = trlweSymEncrypt<typename P::targetP>(
                    p,
                    P::α, sk.key.get<typename P::targetP>());
            }
}

template <class P>
inline void annihilatekeyegen(AnnihilateKey<P> &ahk, const SecretKey &sk){
    for(int i = 0; i < P::nbit; i++){
        Polynomial<P> autokey;
        Automorphism<P>(autokey, sk.key.get<P>(), (1<<(P::nbit-i))+1);
        ahk[i] = trgswfftSymEncrypt<P>(autokey, P::α, sk.key.get<P>());
    }
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

struct GateKeyNTT {
    BootstrappingKeyNTT<lvl01param> bknttlvl01;
    KeySwitchingKey<lvl10param> ksk;
    // GateKey(const SecretKey &sk);
    GateKeyNTT(const GateKeywoFFT &gkwofft);
    GateKeyNTT() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(ksk, bknttlvl01);
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