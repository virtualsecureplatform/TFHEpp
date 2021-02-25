#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <params.hpp>
#include <tlwe.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>

namespace TFHEpp {

template<class P>
inline BootStrappingKeyFFT<P> bkfftgen(SecretKey &sk){
    BootStrappingKeyFFT<P> bkfft;
    for (int i = 0; i < P::domainP::n; i++)
            bkfft[i] = trgswfftSymEncrypt<typename P::targetP>(
                static_cast<typename make_signed<typename P::targetP::T>::type>(sk.key.get<typename P::domainP>()[i]), P::targetP::α, sk.key.get<typename P::targetP>());
    return bkfft;
}

struct GateKey {
    BootStrappingKeyFFT<lvl01param> bkfftlvl01;
    KeySwitchingKey<lvl10param> ksk;
    GateKey(SecretKey sk){
        //Generete bkfft
        bkfftlvl01 = bkfftgen<lvl01param>(sk);
        
        //Generete ksk
        for (int i = 0; i < lvl1param::n; i++)
            for (int j = 0; j < lvl10param::t; j++)
                for (uint32_t k = 0; k < (1 << lvl10param::basebit) - 1; k++)
                    ksk[i][j][k] =
                        tlweSymEncryptlvl0(sk.key.get<lvl1param>()[i] * (k + 1) *
                                            (1ULL << (numeric_limits<typename lvl0param::T>::digits - (j + 1) * lvl10param::basebit)),
                                        lvl0param::α, sk.key.get<lvl0param>());
    }
    GateKey() {}
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(ksk, bkfftlvl01);
    }
};

template<class bsP,class privksP>
struct CircuitKey {
    BootStrappingKeyFFT<bsP> bkfft;
    PrivKeySwitchKey<privksP> privksk;
    CircuitKey(SecretKey sk){
        //Generete bkfft
        bkfft = bkfftgen<bsP>(sk);

        //Generate privksk
        array<typename privksP::domainP::T, privksP::domainP::n + 1> key;
        for (int i = 0; i < privksP::domainP::n; i++) key[i] = sk.key.lvl2[i];
        key[privksP::domainP::n] = -1;
        #pragma omp parallel for collapse(4)
        for (int z = 0; z < 2; z++)
            for (int i = 0; i <= privksP::domainP::n; i++)
                for (int j = 0; j < privksP::t; j++)
                    for (typename privksP::targetP::T u = 0; u < (1 << privksP::basebit) - 1; u++) {
                        TRLWE<typename privksP::targetP> c =
                            trlweSymEncryptZero<typename privksP::targetP>(privksP::targetP::α, sk.key.get<typename privksP::targetP>());
                        c[z][0] += (u + 1) * key[i]
                                << (numeric_limits<typename privksP::targetP::T>::digits - (j + 1) * privksP::basebit);
                        privksk[z][i][j][u] = c;
                    }
    };
    CircuitKey() {}
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(privksk, bkfft);
    }
};

using CircuitKeylvl01 = CircuitKey<lvl02param,lvl21param>;
using CircuitKeylvl02 = CircuitKey<lvl02param,lvl22param>;

struct CloudKey {
    GateKey gk;
    CircuitKeylvl01 ck;
    lweParams params;
    CloudKey(SecretKey sk) : gk(sk), ck(sk) {}
    CloudKey() {}
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(gk.ksk, gk.bkfftlvl01, ck.privksk, ck.bkfft, params);
    }
};
}  // namespace TFHEpp