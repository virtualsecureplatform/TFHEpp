#include <cassert>
#include <cmath>
#include <iostream>
#include <trgsw.hpp>

int main()
{
    using P = TFHEpp::lvl1param;

    TFHEpp::SecretKey sk;

    TFHEpp::Polynomial<P> plain = {};
    plain[0] = 1;

    TFHEpp::TRGSW<P> trgsw;
    TFHEpp::trgswSymEncrypt<P>(trgsw, plain, sk.key.get<P>());

    const TFHEpp::TRGSWFFT<P> trgswfft = TFHEpp::ApplyFFT2trgsw<P>(trgsw);
    const TFHEpp::TRGSWFNT<P> trgswfnt = TFHEpp::ApplyFNT2trgsw<P>(trgsw);

    TFHEpp::TRLWE<P> input;
    for (auto &poly : input)
        for (auto &coef : poly) coef = TFHEpp::UniformTorusRandom<P>();

    TFHEpp::TRLWE<P> fftres;
    TFHEpp::TRLWE<P> fntres;
    TFHEpp::ExternalProduct<P>(fftres, input, trgswfft);
    TFHEpp::ExternalProduct<P>(fntres, input, trgswfnt);

    for (int k = 0; k < P::k + 1; k++)
        for (int i = 0; i < P::n; i++) {
            const int32_t diff =
                static_cast<int32_t>(fftres[k][i] - fntres[k][i]);
            assert(std::abs(diff) <= 8);
        }

    std::cout << "Passed" << std::endl;
    return 0;
}
