#include <iostream>
#include <tfhe++.hpp>

namespace TFHEpp {
template <class P>
inline Polynomial<P> decryptTRLWE3(const TRLWE3<P> &c, const Key<P> &key)
{
    std::array<typename P::T, P::n> partkey;
    for (int i = 0; i < P::n; i++) partkey[i] = key[0 * P::n + i];
    Polynomial<P> mulres, p, keysquare;
    PolyMul<P>(mulres, c[0], partkey);
    Polynomial<P> phase = c[1];
    for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];

    PolyMulNaive<P>(keysquare, partkey, partkey);
    PolyMulNaive<P>(mulres, c[2], keysquare);
    for (int i = 0; i < P::n; i++) phase[i] += mulres[i];

    for (int i = 0; i < P::n; i++)
        p[i] = static_cast<typename P::T>(std::round(phase[i] / P::Δ)) %
               P::plain_modulus;
    return p;
}
}  // namespace TFHEpp

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = TFHEpp::lvl1param;

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey();

    std::cout << "Add Test" << std::endl;
    for (int test = 0; test < num_test; test++) {
        std::uniform_int_distribution<typename P::T> message(
            0, P::plain_modulus - 1);

        TFHEpp::Polynomial<P> p0, p1, pres;
        for (typename P::T &i : p0) i = message(engine);
        for (typename P::T &i : p1) i = message(engine);

        TFHEpp::TRLWE<P> c0 =
            TFHEpp::trlweSymIntEncrypt<P>(p0, sk->key.get<P>());
        TFHEpp::TRLWE<P> c1 =
            TFHEpp::trlweSymIntEncrypt<P>(p1, sk->key.get<P>());
        TFHEpp::TRLWE<P> cres;
        for (int i = 0; i < 2 * P::n; i++) cres[0][i] = c0[0][i] + c1[0][i];
        pres = TFHEpp::trlweSymIntDecrypt<P>(cres, sk->key.get<P>());
        // for (int i = 0; i < P::n; i++)
        // std::cout<<p0[i]<<":"<<p1[i]<<std::endl;
        for (int i = 0; i < P::n; i++)
            assert(pres[i] == (p0[i] + p1[i]) % P::plain_modulus);
    }
    std::cout << "Passed" << std::endl;

    std::cout << "Mul without Relinerarization Test" << std::endl;
    for (int test = 0; test < num_test; test++) {
        std::uniform_int_distribution<typename P::T> message(
            0, P::plain_modulus - 1);

        TFHEpp::Polynomial<P> p0, p1, pres, ptrue;
        for (typename P::T &i : p0) i = message(engine);
        for (typename P::T &i : p1) i = message(engine);

        TFHEpp::TRLWE<P> c0 =
            TFHEpp::trlweSymIntEncrypt<P>(p0, sk->key.get<P>());
        TFHEpp::TRLWE<P> c1 =
            TFHEpp::trlweSymIntEncrypt<P>(p1, sk->key.get<P>());
        TFHEpp::TRLWE3<P> cres;
        TFHEpp::TRLWEMultWithoutRelinerization<P>(cres, c0, c1);
        pres = TFHEpp::decryptTRLWE3<P>(cres, sk->key.get<P>());

        TFHEpp::PolyMulNaive<P>(ptrue, p0, p1);
        for (int i = 0; i < P::n; i++) ptrue[i] %= P::plain_modulus;

        // for (int i = 0; i < P::n; i++)
        // std::cout<<pres[i]<<":"<<ptrue[i]<<std::endl;
        for (int i = 0; i < P::n; i++) assert(pres[i] == ptrue[i]);
    }
    std::cout << "Passed" << std::endl;

    TFHEpp::relinKeyFFT<P> relinkeyfft =
        TFHEpp::relinKeyFFTgen<P>(sk->key.get<P>());

    std::cout << "Mul Test" << std::endl;
    for (int test = 0; test < num_test; test++) {
        std::uniform_int_distribution<typename P::T> message(
            0, P::plain_modulus - 1);

        TFHEpp::Polynomial<P> p0, p1, pres, ptrue;
        for (typename P::T &i : p0) i = message(engine);
        for (typename P::T &i : p1) i = message(engine);

        TFHEpp::TRLWE<P> c0 =
            TFHEpp::trlweSymIntEncrypt<P>(p0, sk->key.get<P>());
        TFHEpp::TRLWE<P> c1 =
            TFHEpp::trlweSymIntEncrypt<P>(p1, sk->key.get<P>());
        TFHEpp::TRLWE<P> cres;
        TFHEpp::TRLWEMult<P>(cres, c0, c1, relinkeyfft);
        pres = TFHEpp::trlweSymIntDecrypt<P>(cres, sk->key.get<P>());

        TFHEpp::PolyMulNaive<P>(ptrue, p0, p1);
        for (int i = 0; i < P::n; i++) ptrue[i] %= P::plain_modulus;

        // for (int i = 0; i < P::n; i++)
        //     std::cout<<pres[i]<<":"<<ptrue[i]<<std::endl;
        for (int i = 0; i < P::n; i++) assert(pres[i] == ptrue[i]);
    }
    std::cout << "Passed" << std::endl;
}