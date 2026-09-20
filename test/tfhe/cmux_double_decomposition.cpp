#include <tfhe++.hpp>
#include <iostream>
#include <stdexcept>
#include <random>

using namespace TFHEpp;

template <class P>
void check_profile()
{
    // Independent sparse integer convolution checks transform size, exact
    // digit rounding, and sign extension, without any encryption noise.
    std::mt19937 rng(20260920);
    Polynomial<P> a{}, b{}, expected_product{};
    for (unsigned i = 0; i < P::n; ++i)
        b[i] = static_cast<typename P::T>(int(rng() & 65535) - 32768);
    for (unsigned i = 0; i < 64; ++i) {
        a[i] = static_cast<typename P::T>(int(rng() & 65535) - 32768);
        for (unsigned j = 0; j < P::n; ++j) {
            if (i + j < P::n) expected_product[i + j] += a[i] * b[j];
            else expected_product[i + j - P::n] -= a[i] * b[j];
        }
    }
    PolynomialInFD<P> af{}, bf{}, productf{};
    TwistIFFT<P>(af, a); TwistIFFT<P>(bf, b);
    MulInFD<P::n>(productf, af, bf);
    Polynomial<P> actual_product{};
    TwistFFT<P>(actual_product, productf);
    if (actual_product != expected_product)
        throw std::runtime_error("integer digit FFT convolution mismatch");
    Key<P> key;
    keyGen<P>(key);
    Polynomial<P> plain{};
    plain[0] = 17 * P::delta_int;
    plain[1] = 23 * P::delta_int;
    plain[P::n - 1] = 5 * P::delta_int;
    auto input = std::make_unique<TRLWE<P>>();
    trlweSymEncrypt<P>(*input, plain, P::α, key);
    auto selector = std::make_unique<TRGSWFFT<P>>();
    for (unsigned bit = 0; bit <= 1; ++bit) {
        Polynomial<P> message{};
        message[0] = bit;
        trgswSymEncrypt<P>(*selector, message, key);
        for (unsigned rotation : {0U, 1U, P::n, P::n + 1}) {
            auto acc = std::make_unique<TRLWE<P>>(*input);
            CMUXFFTwithPolynomialMulByXaiMinusOne<P>(*acc, *selector, rotation);
            Polynomial<P> expected;
            PolynomialMulByXai<P>(expected, plain, bit ? rotation : 0);
            const auto phase = trlwePhase<P>(*acc, key);
            for (unsigned i = 0; i < P::n; ++i) {
                const typename P::T diff = phase[i] - expected[i];
                const typename P::T error = diff >> 127 ? typename P::T{0} - diff : diff;
                if (error >= P::delta_int / 8)
                    throw std::runtime_error("DD CMUX rotation mismatch");
            }
        }
    }
    std::cout << "n=" << P::n << " integer FFT and DD CMUX encrypted selectors/rotations Passed\n";
}

int main() try
{
    check_profile<lvl3simdparam>();
    check_profile<BFVSwitchParamFor<16>>();
}
catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
