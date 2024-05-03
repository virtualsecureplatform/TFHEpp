#pragma once

#ifdef USE_RANDEN
#include <randen.h>
#endif

#include <algorithm>
#include <array>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>

#include "lweParams.hpp"
#include "params.hpp"
#include "random"

namespace TFHEpp {
using namespace std;
struct lweKey {
    Key<lvl0param> lvl0;
    Key<lvlhalfparam> lvlhalf;
    Key<lvl1param> lvl1;
    Key<lvl2param> lvl2;
    Key<lvl3param> lvl3;
    lweKey();
    template <class P>
    Key<P> get() const;
};

struct SecretKey {
    lweKey key;
    lweParams params;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(key.lvl0, key.lvlhalf, key.lvl1, key.lvl2, key.lvl3, params);
    }
};
}  // namespace TFHEpp
