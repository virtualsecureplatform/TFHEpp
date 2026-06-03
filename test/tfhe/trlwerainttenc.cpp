#include <array>
#include <cassert>
#include <iostream>
#include <tfhe++.hpp>

namespace {

template <class P>
bool run_rainttenc_test()
{
    // THIS TEST CAUSES COMPILE ERROR WITH CLANG 19 AND 20
#ifndef __clang__
    constexpr uint32_t num_test = 1000;
    if constexpr (TFHEpp::hasq<P> && TFHEpp::hasqbit<P>) {
        if constexpr (P::qbit == raintt::wordbits) {
            for (int test = 0; test < num_test; test++) {
                std::random_device seed_gen;
                std::default_random_engine engine(seed_gen());
                std::uniform_int_distribution<typename P::T> binary(0, 1);

                TFHEpp::lweKey key;
                std::array<bool, P::n> p;
                for (bool &i : p) i = binary(engine) > 0;
                std::array<typename P::T, P::n> pmu = {};
                for (int i = 0; i < P::n; i++)
                    pmu[i] = p[i] ? P::μ : -P::μ;
                TFHEpp::TRLWERAINTT<P> craintt;
                TFHEpp::trlweSymEncrypt<P>(craintt, pmu, 3, key.get<P>());
                TFHEpp::TRLWE<P> c;
                for (int k = 0; k <= P::k; k++) {
                    raintt::TwistNTT<typename P::T, P::nbit, false>(
                        c[k], craintt[k], (*TFHEpp::raintttable)[0],
                        (*TFHEpp::raintttwist)[0]);
                    for (int i = 0; i < P::n; i++)
                        c[k][i] = static_cast<int32_t>(c[k][i]) < 0
                                      ? c[k][i] + raintt::P
                                      : c[k][i];
                }
                std::array<bool, P::n> p2 =
                    TFHEpp::trlweSymDecrypt<P>(c, key.get<P>());
                // for (int i = 0; i < TFHEpp::lvl1param::n; i++)
                // std::cout<<test<<":"<<i<<":"<<c[TFHEpp::lvl1param::k][i]<<":"<<(p[i]?1:0)<<":"<<(p2[i]?1:0)<<std::endl;
                for (int i = 0; i < P::n; i++) assert(p[i] == p2[i]);
            }
            return true;
        }
    }
#endif
    return false;
}

}  // namespace

int main()
{
    std::cout << "lvl1" << std::endl;
    bool tested = run_rainttenc_test<TFHEpp::lvl1param>();
    if (!tested)
        std::cout << "Nothing to do" << std::endl;
    std::cout << "Passed" << std::endl;
}
