#pragma once

#include <cstdint>
#include <cmath>

constexpr uint32_t DEF_n = 634;
constexpr uint32_t DEF_Nbit = 10;
constexpr uint32_t DEF_N = 1 << DEF_Nbit;
constexpr uint32_t DEF_l = 3;
constexpr uint32_t DEF_Bgbit = 6;
constexpr uint32_t DEF_Bg = 1 << DEF_Bgbit;
constexpr uint32_t DEF_t = 8;
constexpr uint32_t DEF_basebit = 2;
const double DEF_α = std::pow(2.0, -15);
const double DEF_αbk = std::pow(2.0, -25);
const double DEF_αks = DEF_α;

constexpr uint32_t DEF_nbarbit = 11;
constexpr uint32_t DEF_nbar = 1 << DEF_nbarbit;
constexpr uint32_t DEF_lbar = 4;
constexpr uint32_t DEF_Bgbitbar = 9;
constexpr uint32_t DEF_Bgbar = 1 << DEF_Bgbitbar;
constexpr uint32_t DEF_tbar = 10;
constexpr uint32_t DEF_basebitlvl21 = 3;
const double DEF_αbklvl02 = std::pow(2.0, -44);
const double DEF_αprivks = std::pow(2, -25);

constexpr uint32_t DEF_tlvl22 = 10;
constexpr uint32_t DEF_basebitlvl22 = 3;