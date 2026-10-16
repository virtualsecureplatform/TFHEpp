#pragma once
#include <array>
#include <bit>
#include <cstdint>

namespace TFHEpp {
// Candidate for the CPU N16384/q440-selector/q433-DD/q435-native schedule.
// Noise is the continuous sampler's integer-scale sigma, before truncation.
// Its security/correctness qualification is separate from this parameter type.
struct CKKSNarrowLL14RetunedNoise {
    static constexpr double selector_sigma = 32768.0;
    static constexpr double dd_sigma = 4096.0;
    static constexpr double native_auto_sigma = 4096.0;
    static constexpr std::uint32_t native_auto_primary_bits = 32;
    // Distinguish archived low-noise selectors from the retuned profile.
    static constexpr std::uint64_t selector_magic = 0x4e4152524f573033ULL;
    static constexpr auto dd_metadata() {
        return std::array<std::uint64_t, 6>{0x4e44444e30303032ULL, 16384,
            433, 20, 20, std::bit_cast<std::uint64_t>(dd_sigma)};
    }
};
}
