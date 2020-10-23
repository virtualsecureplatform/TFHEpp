#pragma once

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <params.hpp>
#include <tlwe.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>

namespace TFHEpp {
struct GateKey {
    KeySwitchingKey ksk;
    BootStrappingKeyFFTlvl01 bkfftlvl01;
    GateKey(SecretKey sk);
    GateKey() {}
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(ksk, bkfftlvl01);
    }
};

struct CircuitKey {
    PrivKeySwitchKey privksk;
    BootStrappingKeyFFTlvl02 bkfftlvl02;
    CircuitKey(SecretKey sk);
    CircuitKey() {}
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
        archive(gk.ksk, gk.bkfftlvl01, ck.privksk, ck.bkfftlvl02, params.n,
                params.α, params.Nbit, params.N, params.l, params.Bgbit,
                params.Bg, params.αbk, params.t, params.basebit, params.αks,
                params.μ, params.nbarbit, params.nbar, params.lbar,
                params.Bgbitbar, params.Bgbit, params.αbklvl02, params.tbar,
                params.basebitlvl21, params.αprivks, params.μbar);
    }
};
}  // namespace TFHEpp