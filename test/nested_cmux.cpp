#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using TRGSWLvl1FFT = TFHEpp::TRGSWFFT<TFHEpp::lvl1param>;
using TRLWELvl1 = TFHEpp::TRLWE<TFHEpp::lvl1param>;
using PolyLvl1 = TFHEpp::Polynomial<TFHEpp::lvl1param>;
using SecretKey = TFHEpp::SecretKey;
using Lvl1 = TFHEpp::lvl1param;

PolyLvl1 uint2weight(uint64_t n)
{
    PolyLvl1 w;
    const uint32_t mu = 1u << 31;
    for (size_t i = 0; i < Lvl1::n; i++)
        if (i < 64)
            w[i] = ((n >> i) & 1u) ? mu : 0;
        else
            w[i] = 0;
    return w;
}

TRLWELvl1 trivial_TRLWELvl1(const PolyLvl1 &src)
{
    TRLWELvl1 ret = {};
    ret[TFHEpp::lvl1param::k] = src;
    return ret;
}

void dump_histgram_of_phase_of_TRLWELvl1(std::ostream &os, const PolyLvl1 &src)
{
    std::vector<size_t> hist(10, 0);
    for (size_t i = 0; i < Lvl1::n; i++) {
        size_t v =
            static_cast<size_t>(src[i] / (std::pow<double>(2.0, 32) / 10.0));
        hist.at(v)++;
    }
    for (size_t i = 0; i < 10; i++)
        os << 0.1 * i << ": \t" << hist.at(i) << "\n";
    os << "\n";
}

int main()
{
    const size_t N = 500;
    std::cout << "N = " << N << std::endl;

    SecretKey skey;

    std::vector<TRGSWLvl1FFT, TFHEpp::AlignedAllocator<
                                  TFHEpp::TRGSWFFT<TFHEpp::lvl1param>, 64>>
        guard;
    TFHEpp::Polynomial<TFHEpp::lvl1param> plainpoly = {};
    plainpoly[0] = 1;
    for (size_t i = 0; i < N; i++)
        guard.push_back(
            TFHEpp::trgswfftSymEncrypt<Lvl1>(plainpoly, skey.key.lvl1));

    alignas(64) TRLWELvl1 c1 = trivial_TRLWELvl1(uint2weight(1)),
                          c0 = trivial_TRLWELvl1(uint2weight(0));
    alignas(64) TRLWELvl1 res = c1;
    dump_histgram_of_phase_of_TRLWELvl1(
        std::cout, TFHEpp::trlwePhase<TFHEpp::lvl1param>(res, skey.key.lvl1));
    for (size_t i = 0; i < N; i++) {
        alignas(64) TRLWELvl1 tmp = res;
        TFHEpp::CMUXFFT<Lvl1>(res, guard.at(i), tmp, c0);
    }
    dump_histgram_of_phase_of_TRLWELvl1(
        std::cout, TFHEpp::trlwePhase<TFHEpp::lvl1param>(res, skey.key.lvl1));

    /*
    PolyLvl1 testvec1 = {}, testvec2 = {};
    for (size_t i = 0; i < Lvl1::n; i++) {
        testvec1.at(i) = (1u << 29);
        testvec2.at(i) = (1u << 29);
    }
    TRLWELvl1 c1 = trivial_TRLWELvl1(testvec2),
              c0 = trivial_TRLWELvl1(testvec1);
    TRLWELvl1 res = c1;
    dump_histgram_of_phase_of_TRLWELvl1(std::cout,
                                        phase_of_TRLWELvl1(res, skey));
    for (size_t i = 0; i < N; i++) {
        TRLWELvl1 trlwe0 = res, trlwe1 = {};
        size_t k = rand() % (2 * Lvl1::n);
        TFHEpp::PolynomialMulByXai<Lvl1>(trlwe1[0], trlwe0[0], k);
        TFHEpp::PolynomialMulByXai<Lvl1>(trlwe1[1], trlwe0[1], k);
        TFHEpp::CMUXFFT<Lvl1>(res, guard.at(i), trlwe1, trlwe0);
    }
    dump_histgram_of_phase_of_TRLWELvl1(std::cout,
                                        phase_of_TRLWELvl1(res, skey));
    */
}
