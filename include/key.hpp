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

struct lweParams {
    lvl0param lvl0;
    lvl1param lvl1;
    lvl2param lvl2;
    lvl01param lvl01;
    lvl02param lvl02;
    lvl10param lvl10;
    lvl21param lvl21;
    lvl22param lvl22;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(lvl0, lvl1, lvl2, lvl01, lvl02, lvl10, lvl21, lvl22);
    }

// https://cpprefjp.github.io/lang/cpp20/consistent_comparison.html
    auto operator<=>(const lweParams&) const = default;
};

struct lweKey {
    Key<lvl0param> lvl0;
    Key<lvl1param> lvl1;
    Key<lvl2param> lvl2;
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
        archive(key.lvl0, key.lvl1, key.lvl2, params);
    }
};
}  // namespace TFHEpp