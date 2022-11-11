#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

template <class P, uint32_t int_width>
using TBSR = std::array<TFHEpp::TRLWE<P>, int_width>;

template <class P, uint32_t int_width>
inline void TBSRgen(TBSR<P, int_width> &tbsr)
{
    tbsr = {};
    for (int i = 0; i < int_width; i++)
        for (int j = 0; j < (1 << int_width); j++)
            tbsr[i][1][j] = ((j & (1U << i)) > 0) ? P::μ : -P::μ;
}

template <class P, uint32_t int_width>
void AddByTBSR(TBSR<P, int_width + 1> &res,
               const std::array<TFHEpp::TRGSWFFT<P>, int_width> &A,
               const std::array<TFHEpp::TRGSWFFT<P>, int_width> &B)
{
    TBSRgen<P, int_width + 1>(res);
    for (int i = 0; i < int_width; i++)
        for (int j = 0; j < int_width + 1; j++)
            TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<TFHEpp::lvl2param>(
                res[j], A[i], 2 * P::n - (1U << i));
    for (int i = 0; i < int_width; i++)
        for (int j = 0; j < int_width + 1; j++)
            TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<TFHEpp::lvl2param>(
                res[j], B[i], 2 * P::n - (1U << i));
}

int main()
{
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    constexpr uint32_t int_width = 8;
    static_assert(TFHEpp::lvl1param::nbit >= int_width);

    constexpr uint32_t num_test = 1;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> randint(0, 1 << int_width);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);
    TFHEpp::KeySwitchingKey<iksP> *iksk = new TFHEpp::KeySwitchingKey<iksP>();
    TFHEpp::ikskgen<iksP>(*iksk, *sk);

    std::array<uint32_t, num_test> a, b;
    std::array<TBSR<TFHEpp::lvl2param, int_width + 1>, num_test> res;
    for (uint32_t &i : a) i = randint(engine);
    for (uint32_t &i : b) i = randint(engine);

    std::array<std::vector<uint8_t>, num_test> mua, mub;
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < int_width; j++) mua[i].push_back(a[i] & (1 << j));
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < int_width; j++) mub[i].push_back(b[i] & (1 << j));

    std::array<std::vector<TFHEpp::TLWE<TFHEpp::lvl1param>>, num_test> tlwea,
        tlweb;
    for (int i = 0; i < num_test; i++)
        tlwea[i] = TFHEpp::bootsSymEncrypt(mua[i], *sk);
    for (int i = 0; i < num_test; i++)
        tlweb[i] = TFHEpp::bootsSymEncrypt(mub[i], *sk);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        std::array<TFHEpp::TRGSWFFT<TFHEpp::lvl2param>, int_width> trgswa,
            trgswb;
        for (int i = 0; i < int_width; i++)
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(
                trgswa[i], tlwea[test][i], ek);
        for (int i = 0; i < int_width; i++)
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(
                trgswb[i], tlweb[test][i], ek);
        AddByTBSR<TFHEpp::lvl2param, int_width>(res[test], trgswa, trgswb);
    }
    end = std::chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        std::array<TFHEpp::TLWE<TFHEpp::lvl2param>, int_width + 1> restlwe;
        for (int i = 0; i < int_width + 1; i++)
            TFHEpp::SampleExtractIndex<TFHEpp::lvl2param>(restlwe[i],
                                                          res[test][i], 0);

        std::array<uint32_t, int_width + 1> resbit;
        for (int i = 0; i < int_width + 1; i++)
            resbit[i] = TFHEpp::tlweSymDecrypt<TFHEpp::lvl2param>(restlwe[i],
                                                                  sk->key.lvl2)
                            ? 1
                            : 0;
        int resint = 0;
        for (int i = 0; i < int_width + 1; i++) resint += resbit[i] << i;
        std::cout << a[test] << ":" << b[test] << ":" << resint << std::endl;
        assert(resint == a[test] + b[test]);
    }
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}