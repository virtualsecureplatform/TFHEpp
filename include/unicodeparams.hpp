#pragma once
#include "params.hpp"

namespace TFHEpp {
#ifdef USE_80BIT_SECURITY
constexpr double DEF_α = 2.44e-5;
constexpr double DEF_αbk = 3.73e-9;
constexpr double DEF_αks = 2.44e-5;
constexpr uint32_t DEF_μ = 1U << 29;

const double DEF_αbklvl02 = std::pow(2.0, -44);
const double DEF_αprivks = std::pow(2, -31);
constexpr uint64_t DEF_μbar = 1UL << 61;
#else
const double DEF_α = std::pow(2.0, -15);
const double DEF_αbk = std::pow(2.0, -25);
const double DEF_αks = DEF_α;
constexpr uint32_t DEF_μ = 1U << 29;

const double DEF_αbklvl02 = std::pow(2.0, -44);
const double DEF_αprivks = std::pow(2, -31);
constexpr uint64_t DEF_μbar = 1UL << 61;
#endif

struct lweParams {
    uint32_t n = DEF_N;
    double α = DEF_α;
    uint32_t Nbit = DEF_Nbit;
    uint32_t N = DEF_N;
    uint32_t l = DEF_l;
    uint32_t Bgbit = DEF_Bgbit;
    uint32_t Bg = DEF_Bg;
    double αbk = DEF_αbk;
    uint32_t t = DEF_t;
    uint32_t basebit = DEF_basebit;
    double αks = DEF_α;
    uint32_t μ = DEF_μ;

    uint32_t nbarbit = DEF_nbarbit;
    uint32_t nbar = DEF_nbar;
    uint32_t lbar = DEF_lbar;
    uint32_t Bgbitbar = DEF_Bgbitbar;
    uint32_t Bgbar = DEF_Bgbar;
    double αbklvl02 = DEF_αbklvl02;
    uint32_t tbar = DEF_tbar;
    uint32_t basebitlvl21 = DEF_basebitlvl21;
    double αprivks = DEF_αprivks;
    uint64_t μbar = DEF_μbar;
};
}  // namespace TFHEpp