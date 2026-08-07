#pragma once

#include <cstdint>
#include <limits>

#include "../../params.hpp"

// The CLPX scheme-switching parameters assume the modular-Gaussian lvl2param
// noise members of the default 128-bit parameter set.
#ifdef TFHEPP_DEFAULT_128BIT_PARAMS

namespace TFHEpp {

// Parameters for the TFHE-to-CLPX scheme-switching setting of Nagai et al.
// The CLPX radix/base is fixed to 2.  Three base-2^13 decomposition levels
// retain the noise margin of the original four base-2^9 levels while reducing
// the cost and size of the scheme-switch bootstrapping key by one quarter.
struct SS2CLPXlvl2param : lvl2param {
    static constexpr std::uint32_t l = 3;
    static constexpr std::uint32_t lₐ = l;
    static constexpr std::uint32_t Bgbit = 13;
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr std::uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(static_cast<T>(1)
                            << (std::numeric_limits<T>::digits - 1));
};

struct SS2CLPXlvl02param {
    using domainP = lvl0param;
    using targetP = SS2CLPXlvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr std::uint32_t Addends = 2;
#else
    static constexpr std::uint32_t Addends = 1;
#endif
};

struct SS2CLPXlvlh2param {
    using domainP = lvlhalfparam;
    using targetP = SS2CLPXlvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr std::uint32_t Addends = 2;
#else
    static constexpr std::uint32_t Addends = 1;
#endif
};

struct SS2CLPXlvl22param {
    static constexpr std::uint32_t t = lvl22param::t;
    static constexpr std::uint32_t basebit = lvl22param::basebit;
    static constexpr ErrorDistribution errordist = lvl22param::errordist;
    static const inline double α = SS2CLPXlvl2param::α;
    using domainP = SS2CLPXlvl2param;
    using targetP = SS2CLPXlvl2param;
};

}  // namespace TFHEpp

#endif  // TFHEPP_DEFAULT_128BIT_PARAMS
