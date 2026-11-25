#pragma once

#include "subring_switch.hpp"
#include "tlwe_coefficient_packing.hpp"
#include "tutti.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// TUTTI schedule for refreshing TFHE Boolean TLWEs in bulk: the LL13 ring,
// sparse key and scales, without the opening slot-to-coefficient map and
// without application levels. The slots receive 2 sin(2 pi m / q0).
struct CKKSTuttiBooleanLL13Config {
    using Param = CKKSTuttiParam<13, 4>;
    static constexpr std::uint32_t base_log_q = 24, input_log_delta = 20;
    static constexpr std::uint32_t s2c_plain_log_delta = 14;  // closing map
    static constexpr std::uint32_t log_delta = 22, indicator_log_delta = 24;
    static constexpr std::uint32_t plain_log_delta = 18;
    static constexpr std::uint32_t output_log_q = 46, output_log_delta = 20;
    static constexpr std::size_t blocks = 64, subring_stride = 4;
    static constexpr std::uint32_t auxiliary_log_scale = 20;
    static constexpr std::uint32_t sparse_digit_bits = 1;
    static constexpr std::size_t dense_weight = 512;
    static constexpr double message_bound = 1.0;
    static constexpr bool minimax_sine = false;
    static constexpr bool coefficient_input = true;
    static constexpr double amplitude = 1.0;
};

// Batch refresh of 32-bit Boolean TLWEs (phase +-1/8) through one TUTTI
// bootstrap. L is the TFHE TRLWE parameter whose key the TLWEs use.
//   1. Functional packing places the TLWE phases in consecutive coefficients
//      under the dense CKKS key at 2^48.
//   2. Dividing by 2^23 and reducing to q0 doubles the phase: the two
//      Boolean values land on the peaks of the sine, where its derivative
//      vanishes and input error is squashed.
//   3. The TUTTI bootstrap leaves +-2 in the slots.
//   4. A slot-to-coefficient map returns the values to coefficients; eight
//      more bits of shift give phase +-1/8 at 2^32.
//   5. A subring switch and sample extraction return TLWEs under the TFHE key.
template <class Config, class L>
struct CKKSTuttiBooleanRefreshKey {
    using Schedule = CKKSTuttiSchedule<Config>;
    using P = typename Schedule::Param;
    static_assert(use_multilimb_digit_fft_v<P>);
    static_assert(!Schedule::opening_s2c && Config::amplitude == 1.0);
    static_assert(Config::base_log_q <= 47);
    static constexpr std::size_t batch = P::n / 2;
    // Closing map: slots at the bootstrap output, coefficients at 2^32.
    using Linear = CKKSRefreshLinearSchedule<P, 32, Config::output_log_delta,
        Config::s2c_plain_log_delta, 0, true, Config::auxiliary_log_scale>;
    static_assert(Linear::input_log_q == Schedule::output_log_q);
    static_assert(Config::output_log_delta + 1 <= 29);
    using PreparedPlan = CKKSPreparedCompactRealLinearTransformFFT<P,
        Linear::evaluation_input_log_q, Linear::evaluation_input_log_delta,
        Config::s2c_plain_log_delta>;

    CKKSTuttiKey<Schedule> bootstrap;
    TLWECoefficientPackingKey48<P, L> packing;
    typename Linear::S2CPlan s2c_plan;
    CKKSGaloisKey<P, Linear::evaluation_input_log_q> s2c_input_galois;
    CKKSGaloisKey<P, Linear::evaluation_coefficient_log_q> s2c_output_galois;
    CKKSToTFHESubringSwitchKey32<P, L> to_tfhe;
    // Derived caches; they refer to the members above.
    std::unique_ptr<PreparedTLWECoefficientPacking48<P, L>> packing_ntt;
    std::unique_ptr<PreparedPlan> s2c_fft;
    std::unique_ptr<CKKSRetainedPreparedGaloisKeyFFT<P, Linear::evaluation_input_log_q>>
        s2c_input_fft;
    std::unique_ptr<CKKSRetainedPreparedGaloisKeyFFT<P, Linear::evaluation_coefficient_log_q>>
        s2c_output_fft;
    std::unique_ptr<PreparedCKKSToTFHESubringSwitch32<P, L>> to_tfhe_rns;

    CKKSTuttiBooleanRefreshKey() = default;
    CKKSTuttiBooleanRefreshKey(const CKKSTuttiBooleanRefreshKey &) = delete;
    CKKSTuttiBooleanRefreshKey &operator=(const CKKSTuttiBooleanRefreshKey &) = delete;
};

struct CKKSTuttiBooleanTimings {
    double packing_seconds = 0, bootstrap_seconds = 0, decoding_seconds = 0,
           switching_seconds = 0;
    CKKSTuttiTimings bootstrap;
};

template <class Config, class L>
inline void CKKSTuttiBooleanPrepareKey(CKKSTuttiBooleanRefreshKey<Config, L> &key)
{
    using Key = CKKSTuttiBooleanRefreshKey<Config, L>;
    using P = typename Key::P;
    using Linear = typename Key::Linear;
    key.packing_ntt = std::make_unique<PreparedTLWECoefficientPacking48<P, L>>(key.packing);
    key.s2c_fft = std::make_unique<typename Key::PreparedPlan>(
        key.s2c_plan, std::uint64_t{4} << 30);
    key.s2c_input_fft = std::make_unique<CKKSRetainedPreparedGaloisKeyFFT<
        P, Linear::evaluation_input_log_q>>(key.s2c_input_galois);
    key.s2c_output_fft = std::make_unique<CKKSRetainedPreparedGaloisKeyFFT<
        P, Linear::evaluation_coefficient_log_q>>(key.s2c_output_galois);
    key.to_tfhe_rns = std::make_unique<PreparedCKKSToTFHESubringSwitch32<P, L>>(key.to_tfhe);
}

template <class Config, class L>
inline void CKKSTuttiBooleanKeyGen(CKKSTuttiBooleanRefreshKey<Config, L> &key,
    const Key<L> &tfhe, const Key<typename Config::Param> &dense,
    const Key<typename Config::Param> &sparse)
{
    using RefreshKey = CKKSTuttiBooleanRefreshKey<Config, L>;
    using P = typename RefreshKey::P;
    using Linear = typename RefreshKey::Linear;
    const CKKSNoise noise{P::α, 0};
    CKKSTuttiKeyGen(key.bootstrap, dense, sparse, noise);
    TLWECoefficientPackingKeyGen48(key.packing, tfhe, dense);
    CKKSBuildPackedSlotToCoeffPlan<P, Linear::evaluation_input_log_q,
        Linear::evaluation_input_log_delta, Config::s2c_plain_log_delta>(
        key.s2c_plan, Linear::bsgs_step);
    CKKSGaloisKeyGen<P, Linear::evaluation_input_log_q>(key.s2c_input_galois, dense, noise);
    CKKSGaloisKeyGen<P, Linear::evaluation_coefficient_log_q>(key.s2c_output_galois, dense, noise);
    CKKSToTFHESubringSwitchKeyGen32(key.to_tfhe, dense, tfhe);
    CKKSTuttiBooleanPrepareKey(key);
}

// Refresh up to N/2 TLWEs. The result has the same key and +-1/8 convention.
// `pre_switch`, when given, receives the 2^32 coefficient ciphertext under
// the dense key that feeds the output switch, for an owner-side audit.
template <class Config, class L>
inline void CKKSTuttiBooleanRefresh(std::vector<TLWE<L>> &out,
    const std::vector<TLWE<L>> &in, const CKKSTuttiBooleanRefreshKey<Config, L> &key,
    CKKSTuttiBooleanTimings *timings = nullptr,
    CKKSCiphertext<typename Config::Param, 32, 29> *pre_switch = nullptr)
{
    using RefreshKey = CKKSTuttiBooleanRefreshKey<Config, L>;
    using Schedule = typename RefreshKey::Schedule;
    using P = typename RefreshKey::P;
    using Linear = typename RefreshKey::Linear;
    using Clock = std::chrono::steady_clock;
    constexpr auto base_q = Config::base_log_q;
    constexpr auto aux = Config::auxiliary_log_scale;
    if (in.empty() || in.size() > RefreshKey::batch)
        throw std::invalid_argument("invalid TUTTI Boolean batch size");
    if (!key.packing_ntt) throw std::invalid_argument("TUTTI Boolean key is not prepared");
    auto last = Clock::now();
    const auto record = [&](double &field) {
        const auto now = Clock::now();
        field += std::chrono::duration<double>(now - last).count();
        last = now;
    };
    CKKSTuttiBooleanTimings local;
    auto &times = timings ? *timings : local;

    auto packed = std::make_unique<CKKSCiphertext<P, 48, 44>>();
    key.packing_ntt->apply(*packed, in);
    // phase * 2^48 -> 2 * phase * 2^q0 (mod 2^q0).
    auto lowered = std::make_unique<CKKSCiphertext<P, base_q + 1, 44 - (47 - base_q)>>();
    CKKSScaleDown<P, 48, 44, 47 - base_q>(*lowered, *packed);
    auto coefficients = std::make_unique<typename Schedule::Coefficients>();
    for (std::size_t c = 0; c <= P::k; ++c)
        for (std::size_t i = 0; i < P::n; ++i)
            coefficients->ct[c][i] = ckks_detail::reduceToLevel<P, base_q>(lowered->ct[c][i]);
    record(times.packing_seconds);

    auto slots = std::make_unique<typename Schedule::Output>();
    CKKSTuttiBootstrapCoefficients<Schedule>(*slots, *coefficients, key.bootstrap,
        &times.bootstrap);
    record(times.bootstrap_seconds);

    auto lifted = std::make_unique<typename Linear::EvaluationInput>();
    for (std::size_t c = 0; c <= P::k; ++c)
        for (std::size_t i = 0; i < P::n; ++i)
            lifted->ct[c][i] = ckks_detail::reduceToLevel<P, Schedule::output_log_q>(
                slots->ct[c][i]) << aux;
    auto raised = std::make_unique<typename Linear::EvaluationCoefficients>();
    CKKSRealLinearTransform(*raised, *lifted, *key.s2c_fft, *key.s2c_input_fft,
        *key.s2c_output_fft);
    auto decoded = std::make_unique<typename Linear::Coefficients>();
    CKKSScaleDown<P, Linear::evaluation_coefficient_log_q,
        Linear::evaluation_input_log_delta, aux>(*decoded, *raised);
    // +-2 at scale 2^output_log_delta becomes +-2^29 modulo 2^32.
    auto extracted = std::make_unique<CKKSCiphertext<P, 32, 29>>();
    for (std::size_t c = 0; c <= P::k; ++c)
        for (std::size_t i = 0; i < P::n; ++i)
            extracted->ct[c][i] = ckks_detail::reduceToLevel<P, 32>(
                decoded->ct[c][i] << (28 - Config::output_log_delta));
    record(times.decoding_seconds);

    key.to_tfhe_rns->apply(out, *extracted, in.size());
    record(times.switching_seconds);
    if (pre_switch) *pre_switch = *extracted;
}

} // namespace TFHEpp
#endif
