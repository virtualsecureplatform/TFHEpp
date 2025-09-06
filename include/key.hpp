#pragma once

#include <algorithm>
#include <array>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>

#include "lweParams.hpp"
#include "params.hpp"
#include "utils.hpp"
#include "random"

namespace TFHEpp {
using namespace std;
struct lweKey {
    Key<lvl0param> lvl0;
    Key<lvlhalfparam> lvlhalf;
    Key<lvl1param> lvl1;
    Key<lvl2param> lvl2;
    Key<lvl3param> lvl3;
    lweKey(){
        std::uniform_int_distribution<int32_t> lvl0gen(lvl0param::key_value_min,
                                                    lvl0param::key_value_max);
        std::uniform_int_distribution<int32_t> lvlhalfgen(
            lvlhalfparam::key_value_min, lvlhalfparam::key_value_max);
        std::uniform_int_distribution<int32_t> lvl1gen(lvl1param::key_value_min,
                                                    lvl1param::key_value_max);
        std::uniform_int_distribution<int32_t> lvl2gen(lvl2param::key_value_min,
                                                    lvl2param::key_value_max);
        for (typename lvl0param::T &i : lvl0) i = lvl0gen(generator);
        for (typename lvlhalfparam::T &i : lvlhalf) i = lvlhalfgen(generator);
        for (typename lvl1param::T &i : lvl1) i = lvl1gen(generator);
        for (typename lvl2param::T &i : lvl2) i = lvl2gen(generator);
    #ifdef USE_SUBSET_KEY
        for (int i = 0; i < lvl1param::k * lvl1param::n; i++) lvl2[i] = lvl1[i];
    #endif
    }
    template <class P>
    Key<P> get() const
    {
        if constexpr (std::is_same_v<P, lvl0param>)
            return lvl0;
        else if constexpr (std::is_same_v<P, lvlhalfparam>)
            return lvlhalf;
        else if constexpr (std::is_same_v<P, lvl1param> || std::is_same_v<P, AHlvl1param>)
            return lvl1;
        else if constexpr (std::is_same_v<P, lvl2param> || std::is_same_v<P, AHlvl2param>)
            return lvl2;
        else if constexpr (std::is_same_v<P, lvl3param>)
            return lvl3;
        else
            static_assert(false_v<typename P::T>, "not supported parameter");
    }
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
