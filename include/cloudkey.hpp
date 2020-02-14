#pragma once

#include <params.hpp>
#include <tlwe.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>

namespace TFHEpp {
struct GateKey {
    KeySwitchingKey ksk;
    BootStrappingKeylvl01 bklvl01;
    GateKey(SecretKey sk);
    GateKey() {}
};

struct CircuitKey {
    PrivKeySwitchKey privksk;
    BootStrappingKeylvl02 bklvl02;
    CircuitKey(SecretKey sk);
    CircuitKey() {}
};

struct CloudKey {
    GateKey gk;
    CircuitKey ck;
    lweParams params;
    CloudKey(SecretKey sk) : gk(sk), ck(sk) {}
};
}  // namespace TFHEpp