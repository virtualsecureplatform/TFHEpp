#include <tfhe++.hpp>
#include <iostream>
#include <stdexcept>

int main() try
{
    using namespace TFHEpp;
    using P = lvl3simdparam;
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
                const P::T diff = phase[i] - expected[i];
                const P::T error = diff >> 127 ? P::T{0} - diff : diff;
                if (error >= P::delta_int / 8)
                    throw std::runtime_error("DD CMUX rotation mismatch");
            }
        }
    }
    std::cout << "DD CMUX encrypted selectors and rotations Passed\n";
}
catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
