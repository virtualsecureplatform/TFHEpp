#include <array>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <tfhe++.hpp>

namespace {

struct FixedWeightTestParam {
    using T = uint32_t;
    static constexpr int32_t key_value_min = -1;
    static constexpr int32_t key_value_max = 1;
    static constexpr uint32_t n = 8;
    static constexpr uint32_t k = 1;
    static constexpr uint32_t bfv_key_hamming_weight = 3;
};

int fail(const char *message)
{
    std::cerr << "fixed-weight ternary key test: " << message << std::endl;
    return 1;
}

}  // namespace

int main()
{
    using P = FixedWeightTestParam;
    using T = typename P::T;
    static_assert(TFHEpp::has_bfv_key_hamming_weight<P>::value);

    std::mt19937_64 engine(0x4657485445524e41ULL);
    std::array<uint32_t, P::n> support_counts{};
    uint32_t positive_count = 0;
    uint32_t negative_count = 0;
    constexpr uint32_t trials = 4096;
    const T minus_one = T{0} - T{1};

    for (uint32_t trial = 0; trial < trials; trial++) {
        TFHEpp::Key<P> key{};
        TFHEpp::fixedWeightTernaryKeyGen<P>(key, engine);
        if (!TFHEpp::isFixedWeightTernaryKey<P>(key, P::bfv_key_hamming_weight))
            return fail("sampler violated the exact support invariant");

        for (uint32_t i = 0; i < P::n; i++) {
            if (key[i] == T{0}) continue;
            support_counts[i]++;
            if (key[i] == T{1})
                positive_count++;
            else if (key[i] == minus_one)
                negative_count++;
            else
                return fail("sampler emitted a non-ternary coefficient");
        }
    }

    for (const uint32_t count : support_counts)
        if (count == 0 || count == trials)
            return fail("support selection did not vary across trials");
    if (positive_count == 0 || negative_count == 0)
        return fail("independent sign selection did not emit both signs");

    TFHEpp::Key<P> generated{};
    TFHEpp::keyGen<P>(generated);
    if (!TFHEpp::isFixedWeightTernaryKey<P>(generated,
                                            P::bfv_key_hamming_weight))
        return fail("generic keyGen did not select the fixed-weight branch");

    generated[0] = T{2};
    if (TFHEpp::isFixedWeightTernaryKey<P>(generated,
                                           P::bfv_key_hamming_weight))
        return fail("validator accepted a non-ternary coefficient");

#ifdef TFHEPP_DEFAULT_128BIT_PARAMS
    using BootP = TFHEpp::lvl5bootparam;
    auto boot_key = std::make_unique<TFHEpp::Key<BootP>>();
    TFHEpp::keyGen<BootP>(*boot_key);
    if (!TFHEpp::isFixedWeightTernaryKey<BootP>(*boot_key,
                                                BootP::bfv_key_hamming_weight))
        return fail("lvl5bootparam key does not match the proof sampler");

    boot_key->fill(BootP::T{0});
    bool rejected = false;
    try {
        [[maybe_unused]] auto invalid =
            TFHEpp::bfvboot::MakeBootstrapKey<BootP>(*boot_key, false);
    }
    catch (const std::invalid_argument &) {
        rejected = true;
    }
    if (!rejected)
        return fail("bootstrap-key construction accepted an invalid secret");
#endif

    std::cout << "fixed-weight ternary key test: PASS" << std::endl;
    return 0;
}
