#pragma once

#include <params.hpp>
#include <tlwe.hpp>
#include <trgsw.hpp>
#include <trlwe.hpp>

#include <iostream>

namespace TFHEpp {
struct GateKey {
    KeySwitchingKey ksk;
    BootStrappingKeyFFTlvl01 bkfftlvl01;
    lweParams params;
    GateKey(SecretKey sk)
    {
        for (int i = 0; i < DEF_N; i++)
            for (int j = 0; j < DEF_t; j++)
                for (uint32_t k = 0; k < (1 << DEF_basebit) - 1; k++)
                    ksk[i][j][k] = tlweSymEncryptlvl0(
                        sk.key.lvl1[i] * (k + 1) *
                            (1U << (32 - (j + 1) * DEF_basebit)),
                        DEF_αks, sk.key.lvl0);
        for (int i = 0; i < DEF_n; i++)
            bkfftlvl01[i] = trgswfftSymEncryptlvl1(
                static_cast<int32_t>(sk.key.lvl0[i]), DEF_αbk, sk.key.lvl1);
    }
};

struct CircuitKey {
    PrivKeySwitchKey privksk;
    BootStrappingKeyFFTlvl02 bkfftlvl02;
    lweParams params;
    CircuitKey(SecretKey sk)
    {
        array<uint32_t, DEF_nbar + 1> key;
        for (int i = 0; i < DEF_nbar; i++) key[i] = sk.key.lvl2[i];
        key[DEF_nbar] = -1;
        for (int z = 0; z < 2; z++)
            for (int i = 0; i <= DEF_nbar; i++)
                for (int j = 0; j < DEF_tbar; j++)
                    for (int u = 0; u < (1 << DEF_basebitlvl21) - 1; u++) {
                        cout << z << ":" << i << ":" << j << ":" << u << endl;
                        TRLWElvl1 c =
                            trlweSymEncryptZerolvl1(DEF_αprivks, sk.key.lvl1);
                        c[z][0] += (u + 1) * key[i]
                                   << (32 - (j + 1) * DEF_basebitlvl21);
                        privksk[z][i][j][u] = c;
                    }
        for (int i = 0; i < DEF_n; i++)
            bkfftlvl02[i] =
                trgswfftSymEncryptlvl2(static_cast<int32_t>(sk.key.lvl0[i]),
                                       DEF_αbklvl02, sk.key.lvl2);
    }
};
}  // namespace TFHEpp