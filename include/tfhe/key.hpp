#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <type_traits>

#include "lweParams.hpp"
#include "params.hpp"
#include "random"
#include "utils.hpp"

namespace TFHEpp {
using namespace std;

template <class P, class = void>
struct hasell : false_type {};

template <class P>
struct hasell<P, void_t<decltype(P::ell)>> : true_type {};

template <class P>
void keyGen(Key<P>& key)
{
    if constexpr (hasell<P>::value) {
        static_assert(P::k == 1, "block-binary key generation assumes k = 1");
        static_assert(P::ell > 0, "block-binary ell must be positive");
        static_assert(P::n % P::ell == 0,
                      "block-binary dimension must be divisible by ell");
        static_assert(P::key_value_min == 0 && P::key_value_max == 1,
                      "block-binary key generation assumes binary key values");

        fill(key.begin(), key.end(), typename P::T{0});
        uniform_int_distribution<uint32_t> blockgen(0, P::ell);
        for (uint32_t block = 0; block < P::n / P::ell; block++) {
            const uint32_t offset = blockgen(generator);
            if (offset < P::ell)
                key[block * P::ell + offset] = typename P::T{1};
        }
    }
    else {
        uniform_int_distribution<int32_t> gen(P::key_value_min,
                                              P::key_value_max);
        for (typename P::T& i : key) i = gen(generator);
    }
}

struct lweKey {
    std::tuple<Key<lvl0param>, Key<lvlhalfparam>, Key<lvl1param>,
               Key<lvl2param>, Key<lvl3param>>
        keys;
    lweKey()
    {
        keyGen<lvl0param>(std::get<Key<lvl0param>>(keys));
        keyGen<lvlhalfparam>(std::get<Key<lvlhalfparam>>(keys));
        keyGen<lvl1param>(std::get<Key<lvl1param>>(keys));
        keyGen<lvl2param>(std::get<Key<lvl2param>>(keys));
#ifdef USE_SUBSET_KEY
        for (int i = 0; i < lvl1param::k * lvl1param::n; i++)
            std::get<Key<lvl2param>>(keys)[i] =
                std::get<Key<lvl1param>>(keys)[i];
#endif
    }
    template <class P>
    Key<P> get() const
    {
        return std::get<Key<P>>(keys);
    }
};

struct SecretKey {
    lweKey key;
    lweParams params;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(key.keys, params);
    }
};
}  // namespace TFHEpp
