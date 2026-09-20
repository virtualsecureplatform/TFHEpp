#include <tfhe++.hpp>
#include <iostream>
#include <stdexcept>

using namespace TFHEpp;

// Isolate the 32-bit reverse switch. Only its three evaluation keys are
// needed; do not regenerate the much larger forward bootstrap key here.
int main() try
{
    using P = BFVSwitchParamFor<16>;
    using Reverse = BFVToTFHEKeySwitchParamFor<16>;
    Key<P> ring_key;
    Key<lvl1param> io_key;
    Key<lvlhalfparam> half_key;
    keyGen<P>(ring_key);
    keyGen<lvl1param>(io_key);
    keyGen<lvlhalfparam>(half_key);
    TFHEBFVSwitchKeyFor<16> key;
    key.to_half = std::make_unique<KeySwitchingKey<lvl1hparam>>();
    ikskgen<lvl1hparam>(*key.to_half, io_key, half_key);
    key.to_tfhe = std::make_unique<KeySwitchingKey<Reverse>>();
    ikskgen<Reverse>(*key.to_tfhe, ring_key, io_key);
    key.boolean_bootstrap = std::make_unique<BootstrappingKeyFFT<lvlh1param>>();
    bkfftgen<lvlh1param>(*key.boolean_bootstrap, half_key, io_key);
    const std::array<uint32_t, 16> values{
        0, 1, 2, 3, 15, 16, 255, 256, 65535, 65536,
        0x3fffffffU, 0x40000000U, 0x80000000U,
        3298317160U, 4294836225U, 0xffffffffU};
    auto ciphertext = std::make_unique<TRLWE<P>>();
    for (uint32_t value : values) {
        for (int sign : {-1, 0, 1}) {
            Polynomial<P> plain{};
            plain[0] = P::T{value} * P::delta_int;
            // 1/16 plaintext unit is below the selected 256-product envelope.
            if (sign < 0) plain[0] -= P::delta_int / 16;
            if (sign > 0) plain[0] += P::delta_int / 16;
            trlweSymEncrypt<P>(*ciphertext, plain, P::α, ring_key);
            std::array<TLWE<lvl1param>, 32> result;
            BFVToTFHE(result, *ciphertext, key);
            uint32_t decoded = 0;
            for (unsigned bit = 0; bit < 32; ++bit)
                decoded |= uint32_t(tlweSymDecrypt<lvl1param>(result[bit], io_key)) << bit;
            if (decoded != value)
                throw std::runtime_error("BFV16 reverse mismatch: " + std::to_string(value)
                    + " error sign " + std::to_string(sign) + " got " + std::to_string(decoded));
        }
    }
    std::cout << "BFV16 reverse: 48 encrypted boundary/error cases Passed\n";
}
catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
}
