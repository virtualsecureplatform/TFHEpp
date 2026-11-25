#include <ckks/tutti.hpp>

#include <iostream>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace {

template <std::uint32_t Nbit, std::size_t Stride, std::uint32_t Auxiliary,
          std::uint32_t OutputDelta>
struct TestConfig {
    using Param = TFHEpp::CKKSTuttiParam<Nbit, 4>;
    static constexpr std::uint32_t base_log_q = 36, input_log_delta = 32;
    static constexpr std::uint32_t s2c_plain_log_delta = 40;
    static constexpr std::uint32_t log_delta = 40, indicator_log_delta = 44;
    static constexpr std::uint32_t plain_log_delta = 40;
    static constexpr std::uint32_t output_log_q = 120, output_log_delta = OutputDelta;
    static constexpr std::size_t blocks = 4, subring_stride = Stride;
    static constexpr std::uint32_t auxiliary_log_scale = Auxiliary;
    static constexpr std::uint32_t sparse_digit_bits = 4;
    static constexpr double message_bound = 0.25;
    static constexpr bool minimax_sine = true;
};

template <class Config>
void testBootstrap(std::size_t dense_weight)
{
    using namespace TFHEpp;
    using Schedule = CKKSTuttiSchedule<Config>;
    using P = typename Schedule::Param;
    Key<P> dense, sparse;
    CKKSTuttiDenseKeyGen<P>(dense, dense_weight);
    CKKSTuttiSparseKeyGen<Schedule>(sparse);
    auto key = std::make_unique<CKKSTuttiKey<Schedule>>();
    auto invalid = sparse;
    invalid[0] = typename P::T{0};
    bool rejected = false;
    try { CKKSTuttiKeyGen(*key, dense, invalid); }
    catch (const std::invalid_argument &) { rejected = true; }
    if (!rejected) throw std::runtime_error("unanchored sparse key accepted");
    CKKSTuttiKeyGen(*key, dense, sparse);

    CKKSSlotVector<P> slots, decoded;
    for (std::size_t i = 0; i < slots.size(); ++i)
        slots[i] = (static_cast<int>(i % 16) - 8) / 32.0;
    typename Schedule::Input input;
    ckksSlotEncrypt<P>(input, slots, dense);
    typename Schedule::Output output;
    CKKSTuttiBootstrap(output, input, *key);
    ckksSlotDecrypt<P>(decoded, output, dense);
    double worst = 0;
    for (std::size_t i = 0; i < slots.size(); ++i)
        worst = std::max(worst, std::abs(decoded[i] - slots[i]));
    std::cout << "N = " << P::n << ", refresh error " << worst << '\n';
    if (!(worst < 3e-4))
        throw std::runtime_error("TUTTI refresh error " + std::to_string(worst));

    constexpr auto q = Schedule::output_log_q;
    constexpr auto scale = Schedule::output_log_delta;
    CKKSDDRelinKeyChain<P, q + scale, scale, 1> relin;
    CKKSDDRelinKeyChainGen<P, q + scale, scale, 1>(relin, dense);
    CKKSCiphertext<P, q - scale, scale> squared;
    CKKSMultRelinThenRescale<P>(squared, output, output, relin.template get<0>());
    ckksSlotDecrypt<P>(decoded, squared, dense);
    for (std::size_t i = 0; i < slots.size(); ++i)
        if (!(std::abs(decoded[i] - slots[i] * slots[i]) < 3e-4))
            throw std::runtime_error("TUTTI post-refresh square error");
}

} // namespace
#endif

int main()
{
#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
    // N = 32 uses the reference root encoder, N = 512 the digit FFT.
    testBootstrap<TestConfig<5, 1, 0, 40>>(0);
    testBootstrap<TestConfig<5, 4, 20, 36>>(16);
    testBootstrap<TestConfig<9, 4, 20, 36>>(64);
    std::cout << "TUTTI encrypted refresh and square passed\n";
#endif
}
