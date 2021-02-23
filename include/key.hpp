#pragma once

#include <randen.h>

#include <algorithm>
#include <array>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <params.hpp>
#include <random>

namespace TFHEpp {
using namespace std;
struct lweKey {
    Key<lvl0param> lvl0;
    Key<lvl1param> lvl1;
    Key<lvl2param> lvl2;
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
    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(key.lvl0, key.lvl1, key.lvl2, params);
    }
};
}  // namespace TFHEpp