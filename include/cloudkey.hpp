#pragma once

#include <params.hpp>
#include <tlwe.hpp>
#include <trgsw.hpp>

namespace TFHEpp {
struct CloudKey {
    KeySwitchingKey ksk;
    array<TRGSWFFTlvl1, DEF_n> bkfftlvl01;
    lweParams params;
    CloudKey(SecretKey sk)
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
}  // namespace TFHEpp