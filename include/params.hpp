#pragma once
#include <array>
#include <cmath>
#include <cstdint>

namespace TFHEpp {
using namespace std;

constexpr uint32_t DEF_μ = 1U << 29;
constexpr uint64_t DEF_μbar = 1ULL << 61;

// Use old 80bit security parameters. It is faster, but not recommended.

#ifdef USE_80BIT_SECURITY
#include <params/CGGI16.hpp>
#elif defined(USE_CGGI19)
#include <params/CGGI19.hpp >
#else
#include <params/128bit.hpp>
#endif

using Keylvl0 = array<uint32_t, DEF_n>;
using Keylvl1 = array<uint32_t, DEF_N>;
using Keylvl2 = array<uint64_t, DEF_nbar>;

using TLWElvl0 = array<uint32_t, DEF_n + 1>;
using TLWElvl1 = array<uint32_t, DEF_N + 1>;
using TLWElvl2 = array<uint64_t, DEF_nbar + 1>;

using Polynomiallvl1 = array<uint32_t, DEF_N>;
using Polynomiallvl2 = array<uint64_t, DEF_nbar>;
using PolynomialInFDlvl1 = array<double, DEF_N>;
using PolynomialInFDlvl2 = array<double, DEF_nbar>;

using TRLWElvl1 = array<Polynomiallvl1, 2>;
using TRLWElvl2 = array<Polynomiallvl2, 2>;
using TRLWEInFDlvl1 = array<PolynomialInFDlvl1, 2>;
using TRLWEInFDlvl2 = array<PolynomialInFDlvl2, 2>;
using DecomposedTRLWElvl1 = array<Polynomiallvl1, 2 * DEF_l>;
using DecomposedTRLWElvl2 = array<Polynomiallvl2, 2 * DEF_lbar>;
using DecomposedTRLWEInFDlvl1 = array<PolynomialInFDlvl1, 2 * DEF_l>;
using DecomposedTRLWEInFDlvl2 = array<PolynomialInFDlvl2, 2 * DEF_lbar>;

using TRGSWlvl1 = array<TRLWElvl1, 2 * DEF_l>;
using TRGSWlvl2 = array<TRLWElvl2, 2 * DEF_lbar>;
using TRGSWFFTlvl1 = array<TRLWEInFDlvl1, 2 * DEF_l>;
using TRGSWFFTlvl2 = array<TRLWEInFDlvl2, 2 * DEF_lbar>;

using BootStrappingKeyFFTlvl01 = array<TRGSWFFTlvl1, DEF_n>;
using BootStrappingKeyFFTlvl02 = array<TRGSWFFTlvl2, DEF_n>;

using KeySwitchingKey =
    array<array<array<TLWElvl0, (1 << DEF_basebit) - 1>, DEF_t>, DEF_N>;
using PrivKeySwitchKey =
    array<array<array<array<TRLWElvl1, (1 << DEF_basebitlvl21) - 1>, DEF_tbar>,
                DEF_nbar + 1>,
          2>;
using PrivKeySwitchlvl22Key = array<
    array<array<array<TRLWElvl2, (1 << DEF_basebitlvl22) - 1>, DEF_tlvl22>,
          DEF_nbar + 1>,
    2>;

struct lweParams {
    uint32_t n = DEF_n;
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