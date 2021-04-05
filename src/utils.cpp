#include "utils.hpp"

namespace TFHEpp {
template <class P>
typename P::T ModularGaussian(typename P::T center, double stdev)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>) {
        // 32bit fixed-point number version
        std::normal_distribution<double> distribution(0., stdev);
        double err = distribution(generator);
        return center + dtot32(err);
    }
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        // 64bit fixed-point number version
        static const double _2p64 = std::pow(2., 64);
        std::normal_distribution<double> distribution(0., 1.0);
        const double val = stdev * distribution(generator) * _2p64;
        const uint64_t ival = static_cast<typename P::T>(val);
        return ival + center;
    }
    else
        static_assert(false_v<typename P::T>, "Undefined Modular Gaussian!");
}
#define INST(P)                                                     \
    template typename P::T ModularGaussian<P>(typename P::T center, \
                                              double stdev)
TFHEPP_EXPLICIT_INST_WRT_LVL0_1_2(INST)
#undef INST

template <class P>
typename P::T modSwitchFromTorus(uint32_t phase)
{
    constexpr uint32_t Mbit = P::nbit + 1;
    constexpr uint32_t Msize = 1U << Mbit;
    if constexpr (std::is_same_v<typename P::T, uint32_t>)
        return (phase + (1U << (31 - Mbit))) >> (32 - Mbit);
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        typename P::T interv =
            ((1ULL << 63) / Msize) * 2;  // width of each intervall
        typename P::T half_interval =
            interv / 2;  // begin of the first intervall

        // Mod Switching (as in modSwitchFromTorus32)
        typename P::T temp = (static_cast<typename P::T>(phase) << 32) +
                             half_interval;  // RIVEDI
        return temp / interv;
    }
    else
        static_assert(false_v<typename P::T>, "Undefined modSwitchFromTorus!");
}
#define INST(P) template typename P::T modSwitchFromTorus<P>(uint32_t phase)
TFHEPP_EXPLICIT_INST_WRT_LVL1_2(INST)
#undef INST

template <class P>
void PolynomialMulByXai(Polynomial<P> &res, const Polynomial<P> &poly,
                        const typename P::T a)
{
    if (a == 0)
        return;
    else if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a];
    }
    else {
        const typename P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa];
    }
}
#define INST(P)                          \
    template void PolynomialMulByXai<P>( \
        Polynomial<P> & res, const Polynomial<P> &poly, const typename P::T a)
TFHEPP_EXPLICIT_INST_WRT_LVL0_1_2(INST)
#undef INST

template <class P>
void PolynomialMulByXaiMinusOne(Polynomial<P> &res, const Polynomial<P> &poly,
                                const typename P::T a)
{
    if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n] - poly[i];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a] - poly[i];
    }
    else {
        const typename P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n] - poly[i];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa] - poly[i];
    }
}
#define INST(P)                                  \
    template void PolynomialMulByXaiMinusOne<P>( \
        Polynomial<P> & res, const Polynomial<P> &poly, const typename P::T a)
TFHEPP_EXPLICIT_INST_WRT_LVL0_1_2(INST)
#undef INST

}  // namespace TFHEpp
