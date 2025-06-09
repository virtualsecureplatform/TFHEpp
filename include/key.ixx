export module tfhepp:key;
import std;

// Cereal imports - noted for special handling with modules
#import <cereal/archives/portable_binary.hpp>
#import <cereal/types/array.hpp>

export import tfhepp:lweParams;
export import tfhepp:params;

export namespace TFHEpp {
// using namespace std; // Removed: Bad practice in module interface
export struct lweKey {
    Key<lvl0param> lvl0;
    Key<lvlhalfparam> lvlhalf;
    Key<lvl1param> lvl1;
    Key<lvl2param> lvl2;
    Key<lvl3param> lvl3;
    lweKey(); // Constructor definition can remain in key.cpp or be moved to an impl unit
    export template <class P>
    Key<P> get() const { // Definition remains in interface as it's a template
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
        // else static_assert(false, "Unsupported param type"); // Optional: or allow compilation to fail
    }
};

export struct SecretKey {
    lweKey key;
    lweParams params;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(key.lvl0, key.lvlhalf, key.lvl1, key.lvl2, key.lvl3, params);
    }
};
}  // namespace TFHEpp
