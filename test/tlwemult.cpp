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

    std::unique_ptr<TFHEpp::relinKeyFFT<P>> relinkeyfft =
        std::make_unique<TFHEpp::relinKeyFFT<P>>();
    *relinkeyfft = TFHEpp::relinKeyFFTgen<P>(sk->key.get<P>());
    std::unique_ptr<TFHEpp::PrivateKeySwitchingKey<privksP>> privksk =
        std::make_unique<TFHEpp::PrivateKeySwitchingKey<privksP>>();
    TFHEpp::privkskgen<privksP>(*privksk, {1}, *sk);

    int count = 0;

    for (int test = 0; test < num_test; test++) {
        std::uniform_int_distribution<typename P::T> message(
            0, P::plain_modulus - 1);

        typename P::T p0, p1, pres, ptrue;
        p0 = message(engine);
        p1 = message(engine);
        ptrue = (p0 * p1) % P::plain_modulus;

        TFHEpp::TLWE<P> c0 = TFHEpp::tlweSymIntEncrypt<P>(p0, sk->key.get<P>());
        TFHEpp::TLWE<P> c1 = TFHEpp::tlweSymIntEncrypt<P>(p1, sk->key.get<P>());
        TFHEpp::TLWE<P> cres;
        TFHEpp::TLWEMult<privksP>(cres, c0, c1, *relinkeyfft, *privksk);
        pres = TFHEpp::tlweSymIntDecrypt<P>(cres, sk->key.get<P>());

        if (pres != ptrue) count++;
        // std::cout << p0 << ":" << p1 << ":" << pres << ":" << ptrue
        //           << std::endl;
        // assert(pres == ptrue);
    }
    std::cout << "Number of erros: " << count << std::endl;
    std::cout << "Currently, we does not support working parameters!"
              << std::endl;
}