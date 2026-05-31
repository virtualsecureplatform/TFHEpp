#include <cmath>
#include <iostream>
#include <limits>
#include <tfhe++.hpp>

namespace {

using P = TFHEpp::lvl3simdparam;

long double abs_i128(__int128_t value)
{
    return value < 0 ? -static_cast<long double>(value)
                     : static_cast<long double>(value);
}

}  // namespace

int main()
{
    constexpr int samples = 512;
    long double sum_square = 0.0L;
    int nonzero = 0;

    for (int i = 0; i < samples; i++) {
        const auto noise = TFHEpp::ModularGaussian<P>(0, P::α);
        const auto centered = static_cast<__int128_t>(noise);
        if (centered != 0) nonzero++;
        const long double magnitude = abs_i128(centered);
        sum_square += magnitude * magnitude;
    }

    const long double rms = std::sqrt(sum_square / samples);
    const long double expected =
        std::ldexp(static_cast<long double>(P::α),
                   std::numeric_limits<typename P::T>::digits);

    std::cout << "nonzero=" << nonzero << "/" << samples
              << " rms=" << static_cast<double>(rms)
              << " expected=" << static_cast<double>(expected) << std::endl;

    if (nonzero < samples * 9 / 10) {
        std::cerr << "128-bit ModularGaussian collapsed to zero" << std::endl;
        return 1;
    }
    if (rms < expected / 4 || rms > expected * 4) {
        std::cerr << "128-bit ModularGaussian has unexpected scale"
                  << std::endl;
        return 1;
    }

    std::cout << "Passed" << std::endl;
}
