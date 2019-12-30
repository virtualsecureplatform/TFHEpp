#pragma once
#include <cstdint>
#include <cmath>

namespace TFHEpp{
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
    const uint32_t DEF_MU = 1<<29;

    const uint32_t DEF_nbarbit = 11;
    const uint32_t DEF_nbar = 1<<DEF_nbarbit;
    const uint32_t DEF_lbar = 4;
    const uint32_t DEF_Bgbitbar = 9;
    const uint32_t DEF_Bgbar = 1<<DEF_Bgbitbar;
    const double DEF_αbklvl02 = std::pow(2.0,-44);
    const uint32_t DEF_tbar = 10;
    const uint32_t DEF_basebitlvl21 = 3;
    const double DEF_αprivks = std::pow(2,-31);

    struct lweParams{
        uint32_t n;
        double α;
        uint32_t Nbit;
        uint32_t N;
        uint32_t l;
        uint32_t Bgbit;
        uint32_t Bg;
        double αbk;
        uint32_t t;
        uint32_t basebit;
        double αks;
        uint32_t MU;

        uint32_t nbarbit;
        uint32_t nbar;
        uint32_t lbar;
        uint32_t Bgbitbar;
        uint32_t Bgbar;
        double αbklvl02;
        uint32_t tbar;
        uint32_t basebitlvl21;
        double αprivks;
        lweParams(){
            n = DEF_n;
            α = DEF_α;
            Nbit = DEF_Nbit;
            N = DEF_N;
            l = DEF_l;
            Bgbit = DEF_Bgbit;
            Bg = DEF_Bg;
            αbk = DEF_αbk;
            t = DEF_t;
            basebit = DEF_basebit;
            αks = DEF_αks;
            MU = DEF_MU;
            nbarbit = DEF_nbarbit;
            nbar = DEF_nbar;
            lbar = DEF_lbar;
            Bgbitbar = DEF_Bgbitbar;
            Bgbar = DEF_Bgbar;
            αbklvl02 = DEF_αbklvl02;
            tbar = DEF_tbar;
            basebitlvl21 = DEF_basebitlvl21;
            αprivks = DEF_αprivks;
        }
    };
}