#include <tfhe++.hpp>

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string_view>

#ifdef TFHEPP_DEFAULT_128BIT_PARAMS
using namespace TFHEpp;

BFVSwitchInput encrypt_bits(unsigned value, const Key<lvl1param> &key)
{
    BFVSwitchInput bits;
    for (unsigned i = 0; i < bits.size(); ++i)
        tlweSymEncrypt<lvl1param>(bits[i], (value >> i) & 1
            ? lvl1param::μ : -lvl1param::μ, lvl1param::α, key);
    return bits;
}

unsigned decrypt_bits(const BFVSwitchOutput &bits, const Key<lvl1param> &key)
{
    unsigned value = 0;
    for (unsigned i = 0; i < bits.size(); ++i)
        value |= unsigned(tlweSymDecrypt<lvl1param>(bits[i], key)) << i;
    return value;
}

int main(int argc, char **argv) try
{
    const bool full = argc == 2 && std::string_view(argv[1]) == "--full";
    if (argc > 2 || (argc == 2 && !full))
        throw std::invalid_argument("usage: bfv_scheme_switch [--full]");
    Key<lvl1param> io_key;
    Key<lvlhalfparam> half_key;
    Key<BFVSwitchParam> bfv_key;
    keyGen<lvl1param>(io_key);
    keyGen<lvlhalfparam>(half_key);
    keyGen<BFVSwitchParam>(bfv_key);
    TFHEBFVSwitchKey key;
    std::cout << "Generating BFV switching keys" << std::endl;
    TFHEBFVSwitchKeyGen(key, io_key, half_key, bfv_key);
    std::cout << "Keys ready" << std::endl;
    using P = BFVSwitchParam;
    auto a = std::make_unique<TRLWE<P>>();
    auto b = std::make_unique<TRLWE<P>>();
    auto product = std::make_unique<TRLWE<P>>();
    // Isolate PBS phase scaling and the one-coefficient trace before testing
    // the full converter. These checks also diagnose backend regressions.
    for (unsigned bit = 0; bit < 2; ++bit) {
        TLWE<lvlhalfparam> small;
        tlweSymEncrypt<lvlhalfparam>(small, bit ? lvlhalfparam::μ : -lvlhalfparam::μ,
                                    lvlhalfparam::α, half_key);
        TLWE<P> booted;
        Polynomial<P> lut;
        lut.fill(P::delta_int / 2);
        GateBootstrappingTLWE2TLWE<TFHEToBFVBootstrapParam>(booted, small, *key.to_bfv, lut);
        booted[P::n] += P::delta_int / 2;
        const auto phase = tlweSymPhase<P>(booted, bfv_key);
        const unsigned decoded = unsigned((phase + P::delta_int / 2) >> 112);
        std::cout << "PBS bit=" << bit << " decoded=" << decoded
                  << " phase/Delta=" << static_cast<long double>(phase) / static_cast<long double>(P::delta_int)
                  << std::endl;
        if (decoded != bit) throw std::runtime_error("isolated PBS mismatch");
        InvSampleExtractIndex<P>(*a, booted, 0);
        AnnihilateKeySwitching<P>(*a, *a, *key.packing);
        const auto packed = trlwePhase<P>(*a, bfv_key);
        const unsigned packed_value = unsigned((packed[0] + P::delta_int / 2) >> 112);
        std::cout << "Packed bit=" << bit << " decoded=" << packed_value << std::endl;
        if (packed_value != bit) throw std::runtime_error("isolated packing mismatch");
    }
    const unsigned limit = full ? 256 : 4;
    for (unsigned i = 0; i < limit; ++i) {
        const unsigned x = full ? i : std::array<unsigned, 4>{0, 1, 128, 255}[i];
        TFHEToBFV(*a, encrypt_bits(x, io_key), key);
        const auto phase = trlwePhase<P>(*a, bfv_key);
        for (unsigned j = 0; j < P::n; ++j) {
            const unsigned decoded = unsigned((phase[j] + P::delta_int / 2) >> 112);
            if (decoded != (j == 0 ? x : 0))
                throw std::runtime_error("forward mismatch x=" + std::to_string(x)
                    + " coefficient=" + std::to_string(j) + " decoded=" + std::to_string(decoded));
        }
        if (!full || i % 16 == 15)
            std::cout << "Forward checked " << x << std::endl;
    }
    // Fresh BFV inputs isolate reverse conversion from the forward switch.
    const std::array<unsigned, 20> boundaries{
        0, 1, 2, 3, 4, 15, 16, 63, 64, 255, 256, 1023, 1024,
        4095, 4096, 16383, 16384, 32768, 65025, 65535};
    for (const auto x : boundaries) {
      for (const int error_sign : {-1, 0, 1}) {
        Polynomial<P> plain{};
        plain[0] = P::T{x} * P::delta_int;
        // +/-1/512 plaintext unit lies within the checked-in conditional
        // product envelope and exercises signed errors around carry edges.
        if (error_sign < 0) plain[0] -= P::delta_int / 512;
        if (error_sign > 0) plain[0] += P::delta_int / 512;
        trlweSymEncrypt<P>(*a, plain, P::α, bfv_key);
        BFVSwitchOutput out;
        BFVToTFHE(out, *a, key);
        const auto decoded = decrypt_bits(out, io_key);
        if (decoded != x)
            throw std::runtime_error("reverse mismatch " + std::to_string(x)
                + " got " + std::to_string(decoded));
      }
    }
    std::cout << "Reverse boundary cases passed" << std::endl;
    for (const auto [x, y] : std::array<std::pair<unsigned, unsigned>, 6>{
             {{0, 255}, {1, 255}, {255, 255}, {128, 128}, {127, 129}, {254, 255}}}) {
        TFHEToBFV(*a, encrypt_bits(x, io_key), key);
        TFHEToBFV(*b, encrypt_bits(y, io_key), key);
        TRLWEMultFullDD<P>(*product, *a, *b, *key.relin);
        BFVSwitchOutput out;
        BFVToTFHE(out, *product, key);
        const unsigned decoded = decrypt_bits(out, io_key);
        if (decoded != x * y)
            throw std::runtime_error("product mismatch " + std::to_string(x)
                + "*" + std::to_string(y) + " got " + std::to_string(decoded));
        std::cout << x << "*" << y << "=" << decoded << std::endl;
    }
    std::cout << "BFV scheme switching Passed" << std::endl;
}
catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
}
#else
int main() { std::cerr << "Requires default 128-bit parameters\n"; return 1; }
#endif
