#include <algorithm>
#include <cassert>
#include <cuhe++.hpp>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

namespace cuHEpp {
template <uint32_t Nbit, uint8_t radixbit>
void NTTdebug(
    std::array<INTorus, 1 << Nbit> &res,
    const std::array<INTorus, 1 << Nbit> &table,
    std::array<std::array<INTorus, 1 << Nbit>,
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
}  // namespace cuHEpp

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> Bgdist(0, TFHEpp::lvl1param::Bg);
    std::uniform_int_distribution<uint32_t> Torus32dist(
        0, std::numeric_limits<typename TFHEpp::lvl1param::T>::max());

    const std::unique_ptr<
        const std::array<std::array<cuHEpp::INTorus, TFHEpp::lvl1param::n>, 2>>
        twistlvl1 = cuHEpp::TwistGen<TFHEpp::lvl1param::nbit>();
    const std::unique_ptr<
        const std::array<std::array<cuHEpp::INTorus, TFHEpp::lvl1param::n>, 2>>
        tablelvl1 = cuHEpp::TableGen<TFHEpp::lvl1param::nbit>();

    for (int i = 0; i < TFHEpp::lvl1param::n; i++)
        assert(((*twistlvl1)[0][i] * (*twistlvl1)[1][i]).value == 1);
    for (int i = 1; i < TFHEpp::lvl1param::n; i++)
        assert((*twistlvl1)[0][i].value != 1);
    for (int i = 0; i < TFHEpp::lvl1param::n; i++)
        assert(((*tablelvl1)[0][i] * (*tablelvl1)[1][i]).value == 1);
    for (int i = 1; i < TFHEpp::lvl1param::n; i++)
        assert((*tablelvl1)[0][i].value != 1);
    assert(((*tablelvl1)[1][TFHEpp::lvl1param::n / 2 - 1] * (*tablelvl1)[1][1])
               .value == cuHEpp::P - 1);

    // std::cout<<"Bit Reverse"<<std::endl;
    // for(int i = 0; i<4; i++)
    // std::cout<<i<<":"<<cuHEpp::BitReverse<4>(i)<<std::endl;

    // std::cout << "Start INTT radix match test." << std::endl;
    // for (int test = 0; test < num_test; test++) {
    //     // std::array<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::n>
    //     a,res; std::array<cuHEpp::INTorus,TFHEpp::lvl1param::n> res1,res2;
    //     TFHEpp::Polynomial<TFHEpp::lvl1param> a;
    //     for (uint32_t &i : a) i = Bgdist(engine);
    //     for (int i = 0; i < TFHEpp::lvl1param::n; i++) res1[i] = res2[i] =
    //     cuHEpp::INTorus(a[i]); for (int i = 0; i < TFHEpp::lvl1param::n; i++)
    //         assert(a[i] == res1[i].value);
    //     cuHEpp::INTT<TFHEpp::lvl1param::nbit,5>(res1, tablelvl1[1]);
    //     cuHEpp::INTT<TFHEpp::lvl1param::nbit,6>(res2, tablelvl1[1]);
    //     // for (int i = 0; i < TFHEpp::lvl1param::n; i++)
    //         // std::cout<<res1[i].value<<":"<<res2[i].value<<std::endl;
    //     for (int i = 0; i < TFHEpp::lvl1param::n; i++){
    //         // std::cout<<i<<std::endl;
    //         assert(res1[i].value==res2[i].value);
    //     }
    // }
    // std::cout << "INTT radix match PASS" << std::endl;

    // std::cout<<"Start radixButterfly Test"<<std::endl;
    // for (int test = 0; test < num_test; test++) {
    //     constexpr uint8_t radixbit = 3;
    //     constexpr uint32_t radix = 1U<<radixbit;
    //     std::array<uint32_t,1<<radixbit> a;
    //     std::array<cuHEpp::INTorus,1<<radixbit> res;
    //     for (uint32_t &i : a) i = Bgdist(engine)/32;
    //     const cuHEpp::INTorus shifter =
    //     cuHEpp::INTorus(1)<<(TFHEpp::lvl1param::nbit-radixbit); for (int i =
    //     0; i < radix; i++) res[i] = cuHEpp::INTorus(a[i]);
    //     cuHEpp::INTTradixButterfly<radixbit>(res.data(), 1<<radixbit);
    //     for (int i = 0; i < radix; i++) res[i] *= shifter;
    //     cuHEpp::NTTradixButterfly<radixbit>(res.data(), 1<<radixbit);
    //     const cuHEpp::INTorus invN =
    //     cuHEpp::InvPow2(TFHEpp::lvl1param::nbit);
    //     // for (int i = 0; i < radix; i++)
    //         // std::cout<<res[i].value<<":"<<a[i]<<std::endl;

    //     for (int i = 0; i < radix; i++)
    //         res[i] *= invN;
    //     for (int i = 0; i < (1<<radixbit); i++)
    //         assert(a[i] == res[i].value);
    // }
    // std::cout<<"radixButterfly Test PASS"<<std::endl;

    // std::cout<<"Start radix size NTT test"<<std::endl;
    // for (int test = 0; test < num_test; test++) {
    //     constexpr uint8_t radixbit = 3;
    //     constexpr uint8_t sizebit = 11;
    //     constexpr uint32_t size = 1U<<sizebit;
    //     const std::array<std::array<cuHEpp::INTorus,size>,2> table =
    //     cuHEpp::TableGen<sizebit>(); std::array<uint32_t,size> a;
    //     std::array<cuHEpp::INTorus,size> res;
    //     for (uint32_t &i : a) i = Bgdist(engine);

    //     for (int i = 0; i < size; i++) res[i] = cuHEpp::INTorus(a[i]);
    //     cuHEpp::INTT<sizebit,radixbit>(res, table[1]);
    //     cuHEpp::NTT<sizebit,radixbit>(res, table[0]);
    //     const cuHEpp::INTorus invN = cuHEpp::InvPow2(sizebit);
    //     for (int i = 0; i < size; i++)
    //         res[i] *= invN;
    //     // for (int i = 0; i < radix; i++)
    //         // std::cout<<a[i]<<":"<<res[i].value<<std::endl;
    //     for (int i = 0; i < size; i++)
    //         assert(a[i] == res[i].value);
    // }
    // std::cout<<"radix size NTT PASS"<<std::endl;

    // std::cout<<"Start INTT debug test"<<std::endl;
    // for (int test = 0; test < num_test; test++) {
    //     std::array<cuHEpp::INTorus,TFHEpp::lvl1param::n> res,tres;
    //     std::array<std::array<cuHEpp::INTorus,TFHEpp::lvl1param::n>,4> debug;
    //     std::array<std::array<cuHEpp::INTorus,TFHEpp::lvl1param::n>,TFHEpp::lvl1param::nbit>
    //     tdebug; TFHEpp::Polynomial<TFHEpp::lvl1param> a; for (uint32_t &i :
    //     a) i = Bgdist(engine); for (int i = 0; i < TFHEpp::lvl1param::n; i++)
    //     res[i] = cuHEpp::INTorus(a[i]); for (int i = 0; i <
    //     TFHEpp::lvl1param::n; i++)
    //         assert(a[i] == res[i].value);
    //     cuHEpp::INTT<TFHEpp::lvl1param::nbit,3>(res, tablelvl1[1]);
    //     for (int i = 0; i < TFHEpp::lvl1param::n; i++) tres[i] = res[i];
    //     cuHEpp::NTTdebug<TFHEpp::lvl1param::nbit,1>(tres,
    //     tablelvl1[0],tdebug);
    //     cuHEpp::NTTdebug<TFHEpp::lvl1param::nbit,3>(res, tablelvl1[0],debug);
    //     for(int i = 0;i<TFHEpp::lvl1param::n;i++)
    //     assert(debug[0][i].value==tdebug[0][i].value); for(int stage = 0;
    //     stage<3; stage++){
    //         std::cout<<stage<<std::endl;
    //         std::cout<<debug[2][186].value<<":"<<debug[2][250].value<<std::endl;
    //         for(int i = 0;i<TFHEpp::lvl1param::n;i++)
    //         if(debug[stage+1][i].value!=tdebug[3*(stage+1)][i].value)
    //         std::cout<<"index:"<<i<<",value:"<<debug[stage+1][i].value<<":"<<tdebug[3*(stage+1)][i].value<<std::endl;
    //         for(int i = 0;i<TFHEpp::lvl1param::n;i++)
    //         assert(debug[stage+1][i].value==tdebug[3*(stage+1)][i].value);
    //     }
    // }
    // std::cout<<"INTT debug PASS"<<std::endl;
    // {
    //     constexpr uint8_t radixbit = 3;
    //     constexpr uint32_t radix = 1U<<radixbit;
    //     constexpr uint32_t leading = 58;
    //     const std::array<uint64_t,1<<radixbit> a =
    //     {768,4294967295,13070072148996651454,2891126649239098500,17210682317612006867,17933073005031899459,15032082446739191299,6595323810082472591};
    //     const std::array<uint64_t,1<<radixbit> tres =
    //     {0,0,1024,1024,1024,0,1024,2048}; const
    //     std::array<uint64_t,1<<radixbit> fres =
    //     {1152903912152367104,17293840157262217217,1152903912152368128,17293840157262218241,1152903912152368128,17293840157262217217,1152903912152368128,17293840157262219265};
    //     std::array<cuHEpp::INTorus, radix> res,ires;
    //     for (int i = 0; i < radix; i++) res[i] = cuHEpp::INTorus(a[i],false);
    //     for (uint32_t i = 1; i < radix; i++) res[i] *=
    //     tablelvl1[0][cuHEpp::BitReverse<radixbit>(i)*leading];

    //     for (int i = 0; i < radix; i++) ires[i] =
    //     cuHEpp::INTorus(tres[i]/radix,false);
    //     cuHEpp::INTTradixButterfly<radixbit>(ires.data(), radix);
    //     for(int i = 0; i<radix;i++)
    //     std::cout<<res[i].value<<":"<<ires[i].value<<":"<<tablelvl1[0][cuHEpp::BitReverse<radixbit>(i)*leading].value<<":"<<tablelvl1[1][cuHEpp::BitReverse<radixbit>(i)*leading].value<<":"<<(tablelvl1[0][cuHEpp::BitReverse<radixbit>(i)*leading]*tablelvl1[1][cuHEpp::BitReverse<radixbit>(i)*leading]).value<<std::endl;
    //     for(int i = 0;i<radix;i++) assert(res[i].value == ires[i].value);

    //     cuHEpp::NTTradixButterfly<radixbit>(res.data(), radix);

    //     for(int i = 0; i<radix;i++)
    //     std::cout<<res[i].value<<":"<<tres[i]<<":"<<fres[i]<<std::endl;
    //     for(int i = 0; i<radix;i++) assert(res[i].value==fres[i]);
    //     for(int i = 0; i<radix;i++) assert(res[i].value==tres[i]);
    // }

    std::cout << "Start NTT only test." << std::endl;
    for (int test = 0; test < num_test; test++) {
        // std::array<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::n> a,res;
        std::array<cuHEpp::INTorus, TFHEpp::lvl1param::n> res;
        TFHEpp::Polynomial<TFHEpp::lvl1param> a;
        for (typename TFHEpp::lvl1param::T &i : a) i = Bgdist(engine) / 32;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            res[i] = cuHEpp::INTorus(a[i]);
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            assert(a[i] == res[i].value);
        cuHEpp::INTT<TFHEpp::lvl1param::nbit, 3>(res, (*tablelvl1)[1]);
        // for (int i = 0; i < TFHEpp::lvl1param::n; i++)
        // std::cout<<res[i].value<<":"<<a[i]<<std::endl;
        cuHEpp::NTT<TFHEpp::lvl1param::nbit, 6>(res, (*tablelvl1)[0]);

        const cuHEpp::INTorus invN = cuHEpp::InvPow2(TFHEpp::lvl1param::nbit);
        // std::cout<<"InvN:"<<invN.value<<std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) res[i] *= invN;
        // for (int i = 0; i < TFHEpp::lvl1param::n; i++)
        // std::cout<<res[i].value<<":"<<a[i]<<std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++)
            if (a[i] != res[i].value) std::cout << i << std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            assert(a[i] == res[i].value);
        }
    }
    std::cout << "NTT only test PASS" << std::endl;

    std::cout << "Start LVL1 test." << std::endl;
    for (int test = 0; test < num_test; test++) {
        // std::array<typename TFHEpp::lvl1param::T,TFHEpp::lvl1param::n> a,res;
        TFHEpp::Polynomial<TFHEpp::lvl1param> a, res;
        for (typename TFHEpp::lvl1param::T &i : a) i = Torus32dist(engine);
        std::array<cuHEpp::INTorus, TFHEpp::lvl1param::n> resntt;
        cuHEpp::TwistINTT<typename TFHEpp::lvl1param::T,
                          TFHEpp::lvl1param::nbit>(resntt, a, (*tablelvl1)[1],
                                                   (*twistlvl1)[1]);
        cuHEpp::TwistNTT<typename TFHEpp::lvl1param::T,
                         TFHEpp::lvl1param::nbit>(res, resntt, (*tablelvl1)[0],
                                                  (*twistlvl1)[0]);
        // for (int i = 0; i < TFHEpp::lvl1param::n; i++)
        // std::cout<<res[i]<<":"<<a[i]<<std::endl;
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) assert(a[i] == res[i]);
    }
    std::cout << "NTT Passed" << std::endl;

    for (int test = 0; test < num_test; test++) {
        TFHEpp::Polynomial<TFHEpp::lvl1param> a, b, polymul;
        // for (uint32_t &i : a) i = Bgdist(engine) - TFHEpp::lvl1param::Bg / 2;
        for (typename TFHEpp::lvl1param::T &i : a) i = Bgdist(engine);
        for (typename TFHEpp::lvl1param::T &i : b) i = Torus32dist(engine);

        cuHEpp::PolyMullvl1<typename TFHEpp::lvl1param::T,
                            TFHEpp::lvl1param::nbit>(
            polymul, a, b, (*tablelvl1), (*twistlvl1));

        TFHEpp::Polynomial<TFHEpp::lvl1param> naieve = {};
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            for (int j = 0; j <= i; j++)
                naieve[i] += static_cast<int32_t>(a[j]) * b[i - j];
            for (int j = i + 1; j < TFHEpp::lvl1param::n; j++)
                naieve[i] -= static_cast<int32_t>(a[j]) *
                             b[TFHEpp::lvl1param::n + i - j];
        }
        for (int i = 0; i < TFHEpp::lvl1param::n; i++) {
            assert(std::abs(static_cast<int>(naieve[i] - polymul[i])) <= 1);
        }
    }
    std::cout << "PolyMul Passed" << std::endl;

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
