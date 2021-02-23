#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <params.hpp>
#include <tlwe.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>

namespace TFHEpp {
struct GateKey {
    KeySwitchingKey<lvl10param> ksk;
    BootStrappingKeyFFT<lvl01param> bkfftlvl01;
    GateKey(SecretKey sk);
    GateKey() {}
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(ksk, bkfftlvl01);
    }
};

struct CircuitKey {
    PrivKeySwitchKey<lvl21param> privksk;
    BootStrappingKeyFFT<lvl02param> bkfftlvl02;
    CircuitKey(SecretKey sk);
    CircuitKey() {}
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(privksk, bkfftlvl02);
    }
};

struct CircuitKeylvl22 {
    PrivKeySwitchKey<lvl22param> privksk;
    BootStrappingKeyFFT<lvl02param> bkfftlvl02;
    CircuitKeylvl22(SecretKey sk);
    CircuitKeylvl22() {}
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(privksk, bkfftlvl02);
    }
};

struct CloudKey {
    GateKey gk;
    CircuitKey ck;
    lweParams params;
    CloudKey(SecretKey sk) : gk(sk), ck(sk) {}
    CloudKey() {}
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(gk.ksk, gk.bkfftlvl01, ck.privksk, ck.bkfftlvl02, params);
    }
};
}  // namespace TFHEpp