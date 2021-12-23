#include <bits/stdint-uintn.h>

#include <iostream>
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = TFHEpp::lvl1param;
    using privksP = TFHEpp::lvl11param;

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey();

    TFHEpp::relinKeyFFT<P> relinkeyfft =
        TFHEpp::relinKeyFFTgen<P>(sk->key.get<P>());
    TFHEpp::PrivateKeySwitchingKey<privksP> privksk;
    TFHEpp::privkskgen<privksP>(privksk, {1}, *sk);

    for (int test = 0; test < num_test; test++) {
        std::uniform_int_distribution<typename P::T> message(
            0, P::plain_modulus - 1);

        typename P::T p0, p1, pres, ptrue;
        p0 = message(engine);
        p1 = message(engine);
        ptrue = (p0 * p1) % P::plain_modulus;

        TFHEpp::TLWE<P> c0 =
            TFHEpp::tlweSymIntEncrypt<P>(p0, P::α, sk->key.get<P>());
        TFHEpp::TLWE<P> c1 =
            TFHEpp::tlweSymIntEncrypt<P>(p1, P::α, sk->key.get<P>());
        TFHEpp::TLWE<P> cres;
        TFHEpp::TLWEMult<privksP>(cres, c0, c1, relinkeyfft, privksk);
        pres = TFHEpp::tlweSymIntDecrypt<P>(cres, sk->key.get<P>());

        std::cout << p0 << ":" << p1 << ":" << pres << ":" << ptrue
                  << std::endl;
        assert(pres == ptrue);
    }
    std::cout << "Passed" << std::endl;
}