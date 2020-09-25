#pragma once

#include <randen.h>

#include <algorithm>
#include <array>
#include <params.hpp>
#include <random>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>

namespace TFHEpp {
using namespace std;
struct lweKey {
    Keylvl0 lvl0;
    Keylvl1 lvl1;
    Keylvl2 lvl2;
    lweKey()
    {
        randen::Randen<uint64_t> engine;
        uniform_int_distribution<uint32_t> binary(0, 1);
        for (uint32_t &i : lvl0) i = binary(engine);
        for (uint32_t &i : lvl1) i = binary(engine);
        for (uint64_t &i : lvl2) i = binary(engine);
    }
};

struct SecretKey {
    lweKey key;
    lweParams params;
    SecretKey()
    {
        lweKey lwekey;
        key = lwekey;
    }
    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(key.lvl0,key.lvl1,key.lvl2, 
        params.n,params.α,params.Nbit,params.N,params.l,params.Bgbit,params.Bg,params.αbk,params.t,params.basebit,params.αks,params.μ,
        params.nbarbit,params.nbar,params.lbar,params.Bgbitbar,params.Bgbit,params.αbklvl02,params.tbar,params.basebitlvl21,params.αprivks,params.μbar);
    }
};
}  // namespace TFHEpp