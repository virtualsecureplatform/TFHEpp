#pragma once

#import <algorithm>
#import <array>
#import <cereal/archives/portable_binary.hpp>
#import <cereal/types/array.hpp>

#import "lweParams.hpp"
#import "params.hpp"
#import "random"

namespace TFHEpp {
using namespace std;
struct lweKey {
    Key<lvl0param> lvl0;
    Key<lvlhalfparam> lvlhalf;
    Key<lvl1param> lvl1;
    Key<lvl2param> lvl2;
    Key<lvl3param> lvl3;
    lweKey(); // Constructor definition can remain in key.cpp
    template <class P>
    Key<P> get() const { // Definition moved here
        if constexpr (std::is_same_v<P, lvl0param>)
            return lvl0;
        else if constexpr (std::is_same_v<P, lvlhalfparam>)
            return lvlhalf;
        else if constexpr (std::is_same_v<P, lvl1param>)
            return lvl1;
        else if constexpr (std::is_same_v<P, lvl2param>)
            return lvl2;
        else if constexpr (std::is_same_v<P, lvl3param>)
            return lvl3;
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
