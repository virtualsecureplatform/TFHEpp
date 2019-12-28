#pragma once
#include <cstdint>
#include <cmath>

const uint32_t DEF_n = 500;
const double DEF_α = 2.44e-5;
const uint32_t DEF_Nbit = 10;
const uint32_t DEF_N = 1<<DEF_Nbit;
const uint32_t DEF_l = 2;
const uint32_t DEF_Bgbit = 10;
const uint32_t DEF_Bg = 1<<DEF_Bgbit;
const double DEF_αbk = 3.73e-9;
const uint32_t DEF_t = 8;
const uint32_t DEF_basebit = 2;
const double DEF_αks = 2.44e-5;

const uint32_t DEF_nbarbit = 11;
const uint32_t DEF_nbar = 1<<DEF_nbarbit;
const uint32_t DEF_lbar = 4;
const uint32_t DEF_Bgibitbar = 9;
const uint32_t DEF_Bgbar = 1<<DEF_Bgibitbar;
const double DEF_αbklvl02 = std::pow(2.0,-44);
const u_int32_t DEF_tbar = 10;
const uint32_t DEF_basebitlvl21 = 3;
const double DEF_αprivks = std::pow(2,-31);
