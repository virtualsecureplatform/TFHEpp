#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <params.hpp>
#include <tlwe.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>

namespace TFHEpp {

struct GateKeywoFFT {
    BootStrappingKey<lvl01param> bklvl01;
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
    BootStrappingKeyFFT<lvl01param> bkfftlvl01;
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
    BootStrappingKeyFFT<bsP> bkfft;
    std::array<PrivKeySwitchKey<privksP>, 2> privksk;
    CircuitKey(const SecretKey &sk);
    CircuitKey() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(privksk, bkfft);
    }
};

using CircuitKeylvl01 = CircuitKey<lvl02param, lvl21param>;
using CircuitKeylvl02 = CircuitKey<lvl02param, lvl22param>;

struct CloudKey {
    GateKey gk;
    CircuitKeylvl01 ck;
    lweParams params;
    CloudKey(SecretKey sk) : gk(sk), ck(sk) {}
    CloudKey() {}
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(gk.ksk, gk.bkfftlvl01, ck.privksk, ck.bkfft, params);
    }
};
}  // namespace TFHEpp