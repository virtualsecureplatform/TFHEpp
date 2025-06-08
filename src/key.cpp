#import <key.hpp>
#import <utils.hpp>

namespace TFHEpp {
lweKey::lweKey()
{
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
    for(int i = 0; i < lvl1param::k*lvl1param::n; i++) lvl2[i] = lvl1[i];
    #endif
}

// Definition of lweKey::get<P>() const removed from here.
// The #define INST(P) and #undef INST are also removed.

}  // namespace TFHEpp
