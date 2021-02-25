#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <params.hpp>
#include <tlwe.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>

namespace TFHEpp {
struct GateKey {
    BootStrappingKeyFFT<lvl01param> bkfftlvl01;
    KeySwitchingKey<lvl10param> ksk;
    GateKey(SecretKey sk);
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
        for (int i = 0; i < bsP::domainP::n; i++)
            bkfft[i] = trgswfftSymEncrypt<bsP::targetP>(
                static_cast<make_signed<typename bsP::targetP::T>>(sk.key.get<typename bsP::domainP>()[i]), bsP::targetP::α, sk.key.get<typename bsP::targetP>());

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