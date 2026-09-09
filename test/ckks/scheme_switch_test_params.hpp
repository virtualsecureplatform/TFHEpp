#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <tfhe++.hpp>

namespace scheme_switch_test {
// Precision-oriented functional profile: REJECTED by the 128-bit calibration
// (111.87-bit CKKS screen). The accompanying compact TFHE-760 switch also
// fails independently. See docs/SchemeSwitchingCalibration.md. The extra
// limbs leave a 58-bit scale through the complete reverse circuit.
struct FullCKKS {
    static constexpr std::int32_t key_value_min = -1, key_value_max = 1;
    static constexpr std::uint32_t nbit = 15, n = 1 << nbit, k = 1;
    using T = TFHEpp::MultiLimbUInt<16>;
    static constexpr std::uint32_t l = 5, lₐ = 5, Bgbit = 16, Bgₐbit = 16;
    static constexpr T Bg = T{1} << Bgbit, Bgₐ = Bg;
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr std::uint32_t l̅ = 128, l̅ₐ = 128;
    static constexpr std::uint32_t B̅gbit = 8, B̅gₐbit = 8;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static inline const double α = std::ldexp(3.2, -1000);
    static constexpr double ckks_min_noise_stddev = 3.2;
};
}
