#pragma once

#include <algorithm>
#include <array>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <cstdint>
#include <numeric>
#include <type_traits>
#include <vector>

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

template <class P, class = void>
struct has_bfv_key_hamming_weight : false_type {};

template <class P>
struct has_bfv_key_hamming_weight<P,
                                  void_t<decltype(P::bfv_key_hamming_weight)>>
    : true_type {};

// Sample the exact secret law used by the sparse-ternary security proof:
// choose a support uniformly from all subsets of size h, then choose every
// nonzero sign independently and uniformly.  The partial Fisher-Yates pass
// makes the ordered support uniform using only h swaps.
template <class P, class URBG>
void fixedWeightTernaryKeyGen(Key<P>& key, URBG& engine)
{
    static_assert(has_bfv_key_hamming_weight<P>::value,
                  "fixed-weight key generation requires a Hamming weight");
    static_assert(P::k == 1, "fixed-weight key generation assumes k = 1");
    static_assert(P::key_value_min == -1 && P::key_value_max == 1,
                  "fixed-weight key generation assumes ternary key values");
    static_assert(P::bfv_key_hamming_weight <= P::n,
                  "fixed-weight key Hamming weight exceeds the dimension");

    using T = typename P::T;
    fill(key.begin(), key.end(), T{0});

    vector<uint32_t> positions(P::n);
    iota(positions.begin(), positions.end(), uint32_t{0});
    for (uint32_t selected = 0; selected < P::bfv_key_hamming_weight;
         selected++) {
        uniform_int_distribution<uint32_t> supportgen(selected, P::n - 1);
        const uint32_t replacement = supportgen(engine);
        swap(positions[selected], positions[replacement]);
    }

    uniform_int_distribution<uint32_t> signgen(0, 1);
    const T minus_one = T{0} - T{1};
    for (uint32_t selected = 0; selected < P::bfv_key_hamming_weight;
         selected++)
        key[positions[selected]] = signgen(engine) == 0 ? minus_one : T{1};
}

// Check both the support size and the ternary alphabet.  Bootstrap-key
// construction uses this for parameter types carrying a fixed-weight marker,
// so externally supplied keys cannot silently violate the proof assumption.
template <class P>
bool isFixedWeightTernaryKey(const Key<P>& key, const uint32_t expected_weight)
{
    using T = typename P::T;
    const T zero{0};
    const T one{1};
    const T minus_one = zero - one;
    uint32_t weight = 0;
    for (const T& coefficient : key) {
        if (coefficient == zero) continue;
        if (coefficient != one && coefficient != minus_one) return false;
        weight++;
    }
    return weight == expected_weight;
}

template <class P>
void keyGen(Key<P>& key)
{
    if constexpr (has_bfv_key_hamming_weight<P>::value) {
        fixedWeightTernaryKeyGen<P>(key, generator);
    }
    else if constexpr (hasell<P>::value) {
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

template <class K>
struct KeyPair {
    K subset;
    K independent;

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(subset, independent);
    }
};

struct lweKey {
#ifdef USE_BLOCK_BINARY
    struct blockbinaryaesKey {
        KeyPair<Key<blockbinaryaeslvl2param>> value;

        template <class Archive>
        void serialize(Archive& archive)
        {
            archive(value);
        }
    };
#endif
    std::tuple<Key<lvl0param>, KeyPair<Key<lvlhalfparam>>,
               KeyPair<Key<lvl1param>>, KeyPair<Key<lvl2param>>,
               KeyPair<Key<lvl3param>>
#ifdef USE_BLOCK_BINARY
               ,
               blockbinaryaesKey
#endif
               >
        keys;
    lweKey()
    {
        keyGen<lvl0param>(std::get<Key<lvl0param>>(keys));
        auto& halfkeys = std::get<KeyPair<Key<lvlhalfparam>>>(keys);
        auto& lvl1keys = std::get<KeyPair<Key<lvl1param>>>(keys);
        auto& lvl2keys = std::get<KeyPair<Key<lvl2param>>>(keys);
        auto& lvl3keys = std::get<KeyPair<Key<lvl3param>>>(keys);
        keyGen<lvlhalfparam>(halfkeys.subset);
        keyGen<lvlhalfparam>(halfkeys.independent);
        keyGen<lvl1param>(lvl1keys.subset);
        keyGen<lvl1param>(lvl1keys.independent);
        keyGen<lvl2param>(lvl2keys.subset);
        keyGen<lvl2param>(lvl2keys.independent);
        keyGen<lvl3param>(lvl3keys.subset);
        keyGen<lvl3param>(lvl3keys.independent);
#ifdef USE_BLOCK_BINARY
        auto& blockaeskeys = std::get<blockbinaryaesKey>(keys).value;
        keyGen<blockbinaryaeslvl2param>(blockaeskeys.subset);
        keyGen<blockbinaryaeslvl2param>(blockaeskeys.independent);
#endif
        if constexpr (lvlhalfparam::k * lvlhalfparam::n >=
                      lvl0param::k * lvl0param::n)
            for (int i = 0; i < lvl0param::k * lvl0param::n; i++)
                halfkeys.subset[i] = static_cast<typename lvlhalfparam::T>(
                    std::get<Key<lvl0param>>(keys)[i]);
        static_assert(lvl1param::k * lvl1param::n >=
                      lvl0param::k * lvl0param::n);
        for (int i = 0; i < lvl0param::k * lvl0param::n; i++)
            lvl1keys.subset[i] = static_cast<typename lvl1param::T>(
                std::get<Key<lvl0param>>(keys)[i]);
        static_assert(lvl2param::k * lvl2param::n >=
                      lvl1param::k * lvl1param::n);
        for (int i = 0; i < lvl1param::k * lvl1param::n; i++)
            lvl2keys.subset[i] = static_cast<lvl2param::T>(
                static_cast<std::make_signed_t<lvl1param::T>>(
                    lvl1keys.subset[i]));
        static_assert(lvl3param::k * lvl3param::n >=
                      lvl2param::k * lvl2param::n);
        for (int i = 0; i < lvl2param::k * lvl2param::n; i++)
            lvl3keys.subset[i] = static_cast<lvl3param::T>(
                static_cast<std::make_signed_t<lvl2param::T>>(
                    lvl2keys.subset[i]));
#ifdef USE_BLOCK_BINARY
        static_assert(blockbinaryaeslvl2param::k * blockbinaryaeslvl2param::n >=
                      lvl1param::k * lvl1param::n);
        for (int i = 0; i < lvl1param::k * lvl1param::n; i++)
            blockaeskeys.subset[i] = static_cast<blockbinaryaeslvl2param::T>(
                static_cast<std::make_signed_t<lvl1param::T>>(
                    lvl1keys.subset[i]));
#endif
    }
    template <class P>
    Key<P> getSubset() const
    {
#ifdef USE_BLOCK_BINARY
        if constexpr (std::is_same_v<P, blockbinaryaeslvl2param> ||
                      std::is_same_v<P, blockbinaryaesAHlvl2param>)
            return std::get<blockbinaryaesKey>(keys).value.subset;
        else
#endif
            if constexpr (std::is_same_v<P, lvl0param>)
            return std::get<Key<lvl0param>>(keys);
        else
            return std::get<KeyPair<Key<P>>>(keys).subset;
    }

    template <class P>
    Key<P> getIndependent() const
    {
#ifdef USE_BLOCK_BINARY
        if constexpr (std::is_same_v<P, blockbinaryaeslvl2param> ||
                      std::is_same_v<P, blockbinaryaesAHlvl2param>)
            return std::get<blockbinaryaesKey>(keys).value.independent;
        else
#endif
            if constexpr (std::is_same_v<P, lvl0param>)
            return std::get<Key<lvl0param>>(keys);
        else
            return std::get<KeyPair<Key<P>>>(keys).independent;
    }
};

struct SecretKey {
    lweKey key;
    lweParams params;

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(key.keys, params);
    }
};
}  // namespace TFHEpp
