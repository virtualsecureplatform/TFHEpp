#pragma once

#include <randen.h>

#include <algorithm>
#include <array>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <params.hpp>
#include <random>
#include <utils.hpp>

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
        for (typename lvl0param::T &i : lvl0) i = binary(engine);
        for (typename lvl1param::T &i : lvl1) i = binary(engine);
        for (typename lvl2param::T &i : lvl2) i = binary(engine);
    }
    template<class P>
    inline Key<P> get(){
        if constexpr(std::is_same_v<P,lvl0param>) return lvl0;
        else if constexpr(std::is_same_v<P,lvl1param>) return lvl1;
        else if constexpr(std::is_same_v<P,lvl2param>) return lvl2;
        else static_assert(false_v<P>, "Undefined Secret Key!");
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