#include <algorithm>
#include <cassert>
#include <iostream>
#include <raintt.hpp>
#include <random>
#include <tfhe++.hpp>

namespace raintt {
template <uint32_t Nbit, uint8_t radixbit>
void NTTdebug(
    std::array<SWord, 1 << Nbit> &res,
    const std::array<SWord, 1 << Nbit> &table,
    std::array<std::array<SWord, 1 << Nbit>,
               Nbit / radixbit + ((Nbit % radixbit) > 0 ? 1 : 0)> &debug)
{
    constexpr uint8_t remainder = ((Nbit - 1) % radixbit) + 1;
    constexpr uint32_t size = 1U << remainder;
    constexpr uint32_t num_block = 1U << (Nbit - remainder);
    for (uint32_t block = 0; block < num_block; block++)
        NTTradixButterfly<remainder>(&res[size * block], size);
    debug[0] = res;
    int stage = 0;
    for (uint8_t sizebit = remainder + radixbit; sizebit <= Nbit;
         sizebit += radixbit) {
        const uint32_t size = 1U << sizebit;
        const uint32_t num_block = 1U << (Nbit - sizebit);
        for (uint32_t block = 0; block < num_block; block++)
            NTTradix<Nbit, radixbit>(&res[size * block], size, num_block,
                                     table);
        stage++;
        debug[stage] = res;
    }
}
}  // namespace raintt

int main()
{
    constexpr uint32_t num_test = 1000;
    constexpr uint radixbit = 3;
    constexpr uint8_t remainder =
        ((TFHEpp::lvl1param::nbit - 1) % radixbit) + 1;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> Bgdist(0, TFHEpp::lvl1param::Bg);
    std::uniform_int_distribution<uint32_t> Torus32dist(
        0, std::numeric_limits<typename TFHEpp::lvl1param::T>::max());
    std::uniform_int_distribution<uint> Pdist(0, raintt::P - 1);
    std::uniform_int_distribution<int> sPdist(-raintt::P, raintt::P - 1);

    std::cout << "wordbits:" << raintt::wordbits << std::endl;

    for (int i = 0; i < num_test; i++) {
        raintt::Word a = Pdist(engine);
        raintt::Word tres =
            (static_cast<raintt::DoubleWord>(a) * raintt::R) % raintt::P;
        raintt::Word res = raintt::MulREDC(a, raintt::R2);
        if (res != tres) {
            std::cout << "REDC:" << static_cast<uint>(tres) << ":"
                      << static_cast<uint>(res) << ":" << static_cast<uint>(a)
                      << std::endl;
            exit(1);
        }
    }
    for (int i = 0; i < num_test; i++) {
        raintt::SWord a = sPdist(engine);
        raintt::SWord tres =
            ((static_cast<raintt::DoubleSWord>(a) * raintt::R) % raintt::P +
             raintt::P) %
            raintt::P;
        raintt::SWord res = raintt::MulSREDC(a, raintt::R2);
        res = res < 0 ? res + raintt::P : res;
        if (res != tres) {
            std::cout << "SREDC:" << i << ":" << static_cast<int>(tres) << ":"
                      << static_cast<int>(res) << ":" << static_cast<int>(a)
                      << std::endl;
            exit(1);
        }
    }

    const std::unique_ptr<const std::array<
        std::array<std::array<raintt::SWord, TFHEpp::lvl1param::n>, 2>, 2>>
        tablelvl1 = raintt::TableGen<TFHEpp::lvl1param::nbit>();
    const std::unique_ptr<
        const std::array<std::array<raintt::SWord, TFHEpp::lvl1param::n>, 2>>
        twistlvl1 = raintt::TwistGen<TFHEpp::lvl1param::nbit, radixbit>();

    for (int i = 0; i < TFHEpp::lvl1param::n; i++)
        assert((raintt::MulREDC((*tablelvl1)[0][0][i],
                                (*tablelvl1)[1][0][i])) == raintt::R);
    for (int i = 1; i < TFHEpp::lvl1param::n; i++)
        assert((*tablelvl1)[0][0][i] != raintt::R);
    // for (int i = 0; i < TFHEpp::lvl1param::n; i++)
    //     assert((raintt::MulREDC((*twistlvl1)[0][i], (*twistlvl1)[1][i])) ==
    //            (*twistlvl1)[0][0]);
    for (int i = 1; i < TFHEpp::lvl1param::n; i++)
        assert((*twistlvl1)[0][i] != 1);
    // assert(((static_cast<raintt::DoubleWord>((*tablelvl1)[1][TFHEpp::lvl1param::n
    // / 2 - 1]) * (*tablelvl1)[1][1])%raintt::P)
    //    == raintt::P - 1);

    std::cout << "Start NTT only test." << std::endl;
    for (int test = 0; test < num_test; test++) {
        // std::array<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::n> a,res;
        std::array<raintt::DoubleSWord, TFHEpp::lvl1param::n> res, temp;
        TFHEpp::Polynomial<TFHEpp::lvl1param> a;
        // for (typename TFHEpp::lvl1param::T &i : a) i = Bgdist(engine);
        // for (typename TFHEpp::lvl1param::T &i : a) i = Pdist(engine);
        for (typename TFHEpp::lvl1param::T &i : a) i = sPdist(engine);
        // for (uint32_t &i : a) i = Bgdist(engine) - TFHEpp::lvl1param::Bg / 2;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            res[i] =
                static_cast<raintt::DoubleSWord>(static_cast<int32_t>(a[i]));
        temp = res;
        raintt::INTT<TFHEpp::lvl1param::nbit, radixbit>(res, (*tablelvl1)[1]);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            if ((i & ((1 << remainder) - 1)) > 1)
                res[i] = raintt::MulSREDC(res[i], raintt::R2);
        // raintt::INTT<TFHEpp::lvl1param::nbit, 1>(temp, (*tablelvl1)[1]);
        // for (int i = 0; i < TFHEpp::lvl1param::n/2+1; i++)
        // if (temp[i] != res[i])
        //  std::cout << i << ":"
        //  <<static_cast<int>(res[i])<<":"<<static_cast<int>(temp[i])<<std::endl;
        raintt::NTT<TFHEpp::lvl1param::nbit, 3>(res, (*tablelvl1)[0]);
        // radix4
        // for(int i = 0; i < TFHEpp::lvl1param::n>>2; i++){
        //     res[i + (TFHEpp::lvl1param::n>>2)] = raintt::MulSREDC(res[i +
        //     (TFHEpp::lvl1param::n>>2)],raintt::R2); res[i + 3 *
        //     (TFHEpp::lvl1param::n>>2)] = raintt::MulSREDC(res[i + 3 *
        //     (TFHEpp::lvl1param::n>>2)],raintt::R2);
        // }
        // radix8
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            if ((i >> (TFHEpp::lvl1param::nbit - 3) & 3) != 0)
                res[i] = raintt::MulSREDC(res[i], raintt::R2);

        const raintt::Word invN =
            (static_cast<raintt::DoubleWord>(
                 raintt::inv_mod<TFHEpp::lvl1param::n>()) *
             raintt::R) %
            raintt::P;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            res[i] = raintt::MulSREDC(res[i], invN);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            res[i] = res[i] < 0 ? res[i] + raintt::P : res[i];
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            a[i] = static_cast<int32_t>(a[i]) < 0 ? a[i] + raintt::P : a[i];
        for (int i = 0; i < TFHEpp::lvl1param::n / 2 + 2; i++)
            if (a[i] != res[i])
                std::cout << i << ":" << static_cast<int>(res[i]) << ":"
                          << static_cast<int32_t>(a[i]) << std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            assert(a[i] == res[i]);
        }
    }
    std::cout << "NTT only test PASS" << std::endl;

    std::cout << "Start LVL1 test." << std::endl;
    for (int test = 0; test < num_test; test++) {
        TFHEpp::Polynomial<TFHEpp::lvl1param> a, res;
        for (typename TFHEpp::lvl1param::T &i : a) i = sPdist(engine);
        // for (typename TFHEpp::lvl1param::T &i : a) i = Pdist(engine);
        // for (uint32_t &i : a) i = Bgdist(engine) - TFHEpp::lvl1param::Bg / 2;
        std::array<raintt::DoubleSWord, TFHEpp::lvl1param::n> resntt;
        raintt::TwistINTT<typename TFHEpp::lvl1param::T,
                          TFHEpp::lvl1param::nbit, false>(
            resntt, a, (*tablelvl1)[1], (*twistlvl1)[1]);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            if ((i & ((1 << remainder) - 1)) > 1)
                resntt[i] = raintt::MulSREDC(resntt[i], raintt::R2);
        raintt::TwistNTT<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::nbit,
                         false>(res, resntt, (*tablelvl1)[0], (*twistlvl1)[0]);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            res[i] =
                static_cast<int32_t>(res[i]) < 0 ? res[i] + raintt::P : res[i];
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            a[i] = static_cast<int32_t>(a[i]) < 0 ? a[i] + raintt::P : a[i];
        for (int i = 0; i < TFHEpp::lvl1param::n / 2; i++)
            if (a[i] != res[i])
                std::cout << i << ":" << static_cast<int32_t>(res[i]) << ":"
                          << static_cast<int32_t>(a[i]) << std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) assert(a[i] == res[i]);
    }
    std::cout << "NTT witout modswitch Passed" << std::endl;
    for (int test = 0; test < num_test; test++) {
        typename TFHEpp::lvl1param::T a = Torus32dist(engine);
        raintt::Word b =
            (((static_cast<uint64_t>(a) * raintt::K) << raintt::shiftamount) +
             a + (1ULL << (32 - 1))) >>
            32;
        typename TFHEpp::lvl1param::T c =
            (static_cast<uint64_t>(b) *
                 ((1ULL << (32 + raintt::wordbits - 1)) / raintt::P) +
             (1ULL << (raintt::wordbits - 1 - 1))) >>
            (raintt::wordbits - 1);
        assert(std::abs(static_cast<int>(a - c)) <=
               (1 << (32 - raintt::wordbits + 1)));
    }
    std::cout << "Modswitch Passed" << std::endl;
    for (int test = 0; test < num_test; test++) {
        // std::array<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::n> a,res;
        TFHEpp::Polynomial<TFHEpp::lvl1param> a, res;
        for (typename TFHEpp::lvl1param::T &i : a) i = Torus32dist(engine);
        std::array<raintt::DoubleSWord, TFHEpp::lvl1param::n> resntt;
        raintt::TwistINTT<typename TFHEpp::lvl1param::T,
                          TFHEpp::lvl1param::nbit, true>(
            resntt, a, (*tablelvl1)[1], (*twistlvl1)[1]);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            if ((i & ((1 << remainder) - 1)) > 1)
                resntt[i] = raintt::MulSREDC(resntt[i], raintt::R2);
        raintt::TwistNTT<typename TFHEpp::lvl1param::T, TFHEpp::lvl1param::nbit,
                         true>(res, resntt, (*tablelvl1)[0], (*twistlvl1)[0]);
        // for (int i = 0; i < TFHEpp::lvl1param::n/2; i++)
        // if(std::abs(static_cast<int>(res[i] - a[i])) >
        // 4)std::cout<<res[i]<<":"<<a[i]<<std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            assert(std::abs(static_cast<int>(res[i] - a[i])) <=
                   (1 << (32 - raintt::wordbits + 1)));
    }
    std::cout << "NTT with modswitch Passed" << std::endl;

    for (int test = 0; test < num_test; test++) {
        TFHEpp::Polynomial<TFHEpp::lvl1param> a, b, polymul;
        for (typename TFHEpp::lvl1param::T &i : a) i = Bgdist(engine);
        // for (uint32_t &i : a) i = Bgdist(engine) - TFHEpp::lvl1param::Bg / 2;
        for (typename TFHEpp::lvl1param::T &i : b) i = Pdist(engine);

        raintt::PolyMullvl1<typename TFHEpp::lvl1param::T,
                            TFHEpp::lvl1param::nbit, false, false>(
            polymul, a, b, (*tablelvl1), (*twistlvl1));

        TFHEpp::Polynomial<TFHEpp::lvl1param> naieve = {};
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            for (int j = 0; j <= i; j++)
                naieve[i] =
                    ((static_cast<int64_t>(naieve[i]) +
                      (static_cast<int64_t>(static_cast<int32_t>(a[j])) *
                       b[i - j]) %
                          raintt::P) +
                     raintt::P) %
                    raintt::P;
            for (int j = i + 1; j < TFHEpp::lvl1param::n; j++)
                naieve[i] = ((static_cast<int64_t>(naieve[i]) -
                              static_cast<int64_t>(static_cast<int32_t>(a[j])) *
                                  b[TFHEpp::lvl1param::n + i - j]) %
                                 raintt::P +
                             raintt::P) %
                            raintt::P;
        }
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            polymul[i] = static_cast<int32_t>(polymul[i]) < 0
                             ? polymul[i] + raintt::P
                             : polymul[i];
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            assert(std::abs(static_cast<int>(naieve[i] - polymul[i])) <= 1);
        }
    }
    std::cout << "PolyMul without modsiwtch Passed" << std::endl;

    for (int test = 0; test < num_test; test++) {
        TFHEpp::Polynomial<TFHEpp::lvl1param> a, b, polymul;
        // for (uint32_t &i : a) i = Bgdist(engine) - TFHEpp::lvl1param::Bg / 2;
        for (typename TFHEpp::lvl1param::T &i : a) i = Bgdist(engine);
        for (typename TFHEpp::lvl1param::T &i : b) i = Torus32dist(engine);

        raintt::PolyMullvl1<typename TFHEpp::lvl1param::T,
                            TFHEpp::lvl1param::nbit, false, true>(
            polymul, a, b, (*tablelvl1), (*twistlvl1));

        TFHEpp::Polynomial<TFHEpp::lvl1param> naieve = {};
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            for (int j = 0; j <= i; j++)
                naieve[i] += static_cast<int32_t>(a[j]) * b[i - j];
            for (int j = i + 1; j < TFHEpp::lvl1param::n; j++)
                naieve[i] -= static_cast<int32_t>(a[j]) *
                             b[TFHEpp::lvl1param::n + i - j];
        }
        // for (int i = 0; i < TFHEpp::lvl1param::n / 2; i++)
        //     std::cout << i << ":" << naieve[i] << ":" << polymul[i]
        //               << std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            // assert(std::abs(static_cast<int>(naieve[i] - polymul[i])) <=
            // (1U<<(TFHEpp::lvl1param::nbit+4)));
            assert(std::abs(static_cast<int>(naieve[i] - polymul[i])) <=
                   (1U << (TFHEpp::lvl1param::nbit +
                           (32 - raintt::wordbits + 1) + 4)));
        }
    }
    std::cout << "PolyMul with modsiwtch Passed" << std::endl;

    // uniform_int_distribution<uint64_t> Bgbardist(0, DEF_Bgbar);
    // uniform_int_distribution<uint64_t> Torus64dist(0, UINT64_MAX);

    // cout << "Start LVL2 test." << endl;
    // for (int test = 0; test < num_test; test++) {
    //     Polynomiallvl2 a;
    //     for (uint64_t &i : a) i = Torus64dist(engine);
    //     PolynomialInFDlvl2 resfft;
    //     TFHEpp::TwistIFFTlvl2(resfft, a);
    //     Polynomiallvl2 res;
    //     TFHEpp::TwistFFTlvl2(res, resfft);
    //     for (int i = 0; i < DEF_N; i++)
    //         assert(abs(static_cast<int64_t>(a[i] - res[i])) <= (1 << 14));
    // }
    // cout << "FFT Passed" << endl;

    return 0;
}
