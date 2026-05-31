#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>
#include <tuple>

namespace {

using P = TFHEpp::ckkslvl3param;

struct SmallCKKSParam {
    static constexpr std::uint32_t nbit = 4;
    static constexpr std::uint32_t n = 1 << nbit;
};

struct SmallMultiLimbCKKSParam {
    static constexpr std::uint32_t nbit = 4;
    static constexpr std::uint32_t n = 1 << nbit;
    using T = TFHEpp::MultiLimbUInt<3>;
};

void fill_key(TFHEpp::Key<P> &key)
{
    for (std::size_t i = 0; i < P::n; i++) {
        const int v = static_cast<int>(i % 3) - 1;
        key[i] = static_cast<typename P::T>(v);
    }
}

template <class Param>
double max_error(const TFHEpp::CKKSSlotVector<Param> &got,
                 const TFHEpp::CKKSSlotVector<Param> &want)
{
    double err = 0.0;
    for (std::size_t i = 0; i < Param::n / 2; i++)
        err = std::max(err, std::abs(got[i] - want[i]));
    return err;
}

double max_error(const TFHEpp::CKKSSlotVector<P> &got,
                 const TFHEpp::CKKSSlotVector<P> &want)
{
    return max_error<P>(got, want);
}

void require_close(const TFHEpp::CKKSSlotVector<P> &got,
                   const TFHEpp::CKKSSlotVector<P> &want, double tol,
                   const char *label)
{
    const double err = max_error(got, want);
    std::cout << label << " max_error=" << err << std::endl;
    if (err > tol) {
        std::cerr << label << " exceeded tolerance " << tol << std::endl;
        std::exit(1);
    }
}

template <class Param>
void apply_complex_linear(
    TFHEpp::CKKSSlotVector<Param> &out,
    const TFHEpp::CKKSSlotVector<Param> &in,
    const std::vector<TFHEpp::CKKSSlotVector<Param>> &diagonals,
    const std::vector<int> &offsets)
{
    constexpr int half = static_cast<int>(Param::n) / 2;
    out.fill({0.0, 0.0});
    for (std::size_t d = 0; d < diagonals.size(); d++) {
        const int offset = ((offsets[d] % half) + half) % half;
        for (int i = 0; i < half; i++)
            out[i] += diagonals[d][i] * in[(i + offset) % half];
    }
}

template <class Param>
void apply_real_linear(
    TFHEpp::CKKSSlotVector<Param> &out,
    const TFHEpp::CKKSSlotVector<Param> &in,
    const std::vector<TFHEpp::CKKSSlotVector<Param>> &direct_diags,
    const std::vector<int> &direct_offsets,
    const std::vector<TFHEpp::CKKSSlotVector<Param>> &conj_diags,
    const std::vector<int> &conj_offsets)
{
    auto direct = std::make_unique<TFHEpp::CKKSSlotVector<Param>>();
    auto conj_part = std::make_unique<TFHEpp::CKKSSlotVector<Param>>();
    auto conj_in = std::make_unique<TFHEpp::CKKSSlotVector<Param>>();
    for (std::size_t i = 0; i < Param::n / 2; i++)
        (*conj_in)[i] = std::conj(in[i]);
    apply_complex_linear<Param>(*direct, in, direct_diags, direct_offsets);
    apply_complex_linear<Param>(*conj_part, *conj_in, conj_diags,
                                conj_offsets);
    for (std::size_t i = 0; i < Param::n / 2; i++)
        out[i] = (*direct)[i] + (*conj_part)[i];
}

template <class Param>
void apply_complex_stages(
    TFHEpp::CKKSSlotVector<Param> &out,
    const TFHEpp::CKKSSlotVector<Param> &in,
    const TFHEpp::CKKSLinearTransformStages<Param> &stages)
{
    auto cur = std::make_unique<TFHEpp::CKKSSlotVector<Param>>(in);
    auto next = std::make_unique<TFHEpp::CKKSSlotVector<Param>>();
    for (const auto &stage : stages) {
        apply_complex_linear<Param>(*next, *cur, stage.diagonals,
                                    stage.rotation_offsets);
        std::swap(cur, next);
    }
    out = *cur;
}

void test_dense_coeff_slot_diagonals()
{
    using S = SmallCKKSParam;
    constexpr int half = static_cast<int>(S::n) / 2;
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const auto &slot_to_eval = TFHEpp::ckks_detail::ckksSlotToEvalIndex<S>();

    std::array<double, S::n> coeffs{};
    for (std::size_t i = 0; i < S::n; i++)
        coeffs[i] = static_cast<double>(static_cast<int>(i % 7) - 3) / 16.0;

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    auto packed = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    for (int i = 0; i < half; i++) {
        const long double h =
            static_cast<long double>(2 * slot_to_eval[i] + 1);
        std::complex<long double> sum = 0.0L;
        for (std::size_t j = 0; j < S::n; j++) {
            const long double angle =
                pi * h * static_cast<long double>(j) /
                static_cast<long double>(S::n);
            sum += static_cast<long double>(coeffs[j]) *
                   std::complex<long double>(std::cos(angle),
                                             std::sin(angle));
        }
        (*slots)[i] = static_cast<std::complex<double>>(sum);
        (*packed)[i] = {coeffs[i], coeffs[i + half]};
    }

    std::vector<TFHEpp::CKKSSlotVector<S>> direct_diags;
    std::vector<TFHEpp::CKKSSlotVector<S>> conj_diags;
    std::vector<int> direct_offsets;
    std::vector<int> conj_offsets;
    TFHEpp::CKKSBuildCoeffToPackedSlotDiagonals<S>(
        direct_diags, direct_offsets, conj_diags, conj_offsets);

    auto got_packed = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    apply_real_linear<S>(*got_packed, *slots, direct_diags, direct_offsets,
                         conj_diags, conj_offsets);
    const double c2s_err = max_error<S>(*got_packed, *packed);
    std::cout << "CKKS dense coeff-to-packed-slot max_error=" << c2s_err
              << std::endl;
    if (c2s_err > 1e-10) std::exit(1);

    TFHEpp::CKKSBuildPackedSlotToCoeffDiagonals<S>(
        direct_diags, direct_offsets, conj_diags, conj_offsets);
    auto got_slots = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    apply_real_linear<S>(*got_slots, *packed, direct_diags, direct_offsets,
                         conj_diags, conj_offsets);
    const double stc_err = max_error<S>(*got_slots, *slots);
    std::cout << "CKKS dense packed-slot-to-coeff max_error=" << stc_err
              << std::endl;
    if (stc_err > 1e-10) std::exit(1);
}

void test_factorized_coeff_slot_stages()
{
    using S = SmallCKKSParam;
    constexpr int half = static_cast<int>(S::n) / 2;
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const auto &slot_to_eval = TFHEpp::ckks_detail::ckksSlotToEvalIndex<S>();

    std::array<double, S::n> coeffs{};
    for (std::size_t i = 0; i < S::n; i++)
        coeffs[i] = static_cast<double>(static_cast<int>(i % 11) - 5) / 32.0;

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    auto packed_bitrev = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    for (int i = 0; i < half; i++) {
        const long double h =
            static_cast<long double>(2 * slot_to_eval[i] + 1);
        std::complex<long double> sum = 0.0L;
        for (std::size_t j = 0; j < S::n; j++) {
            const long double angle =
                pi * h * static_cast<long double>(j) /
                static_cast<long double>(S::n);
            sum += static_cast<long double>(coeffs[j]) *
                   std::complex<long double>(std::cos(angle),
                                             std::sin(angle));
        }
        (*slots)[i] = static_cast<std::complex<double>>(sum);

        const auto br = TFHEpp::CKKSBitReverseSlotIndex<S>(
            static_cast<std::uint32_t>(i));
        (*packed_bitrev)[br] = {coeffs[i], coeffs[i + half]};
    }

    TFHEpp::CKKSLinearTransformStages<S> c2s_stages;
    TFHEpp::CKKSBuildCoeffToPackedSlotStages<S>(c2s_stages);
    if (c2s_stages.size() != S::nbit - 1) std::exit(1);
    for (const auto &stage : c2s_stages)
        if (stage.rotation_offsets.size() > 3) std::exit(1);

    auto got_packed = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    apply_complex_stages<S>(*got_packed, *slots, c2s_stages);
    const double c2s_err = max_error<S>(*got_packed, *packed_bitrev);
    std::cout << "CKKS factorized coeff-to-packed-slot max_error=" << c2s_err
              << std::endl;
    if (c2s_err > 1e-10) std::exit(1);

    TFHEpp::CKKSLinearTransformStages<S> fused_c2s_stages;
    TFHEpp::CKKSFuseLinearTransformStages<S>(fused_c2s_stages, c2s_stages, 2);
    if (fused_c2s_stages.size() != 2) std::exit(1);
    auto fused_got_packed = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    apply_complex_stages<S>(*fused_got_packed, *slots, fused_c2s_stages);
    const double fused_c2s_err = max_error<S>(*fused_got_packed, *packed_bitrev);
    std::cout << "CKKS fused factorized coeff-to-packed-slot max_error="
              << fused_c2s_err << std::endl;
    if (fused_c2s_err > 1e-10) std::exit(1);

    auto scaled_expected = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    auto scaled_got = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    constexpr std::complex<double> stage_scale{0.25, -0.125};
    for (int i = 0; i < half; i++)
        (*scaled_expected)[i] = stage_scale * (*packed_bitrev)[i];
    TFHEpp::CKKSScaleLinearTransformStages<S>(c2s_stages, stage_scale);
    apply_complex_stages<S>(*scaled_got, *slots, c2s_stages);
    const double scaled_err = max_error<S>(*scaled_got, *scaled_expected);
    std::cout << "CKKS scaled factorized stage max_error=" << scaled_err
              << std::endl;
    if (scaled_err > 1e-10) std::exit(1);

    TFHEpp::CKKSLinearTransformStages<S> stc_stages;
    TFHEpp::CKKSBuildPackedSlotToCoeffStages<S>(stc_stages);
    if (stc_stages.size() != S::nbit - 1) std::exit(1);
    for (const auto &stage : stc_stages)
        if (stage.rotation_offsets.size() > 3) std::exit(1);

    auto got_slots = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    apply_complex_stages<S>(*got_slots, *packed_bitrev, stc_stages);
    const double stc_err = max_error<S>(*got_slots, *slots);
    std::cout << "CKKS factorized packed-slot-to-coeff max_error=" << stc_err
              << std::endl;
    if (stc_err > 1e-10) std::exit(1);

    TFHEpp::CKKSLinearTransformStages<S> fused_stc_stages;
    TFHEpp::CKKSFuseLinearTransformStages<S>(fused_stc_stages, stc_stages, 2);
    if (fused_stc_stages.size() != 2) std::exit(1);
    auto fused_got_slots = std::make_unique<TFHEpp::CKKSSlotVector<S>>();
    apply_complex_stages<S>(*fused_got_slots, *packed_bitrev, fused_stc_stages);
    const double fused_stc_err = max_error<S>(*fused_got_slots, *slots);
    std::cout << "CKKS fused factorized packed-slot-to-coeff max_error="
              << fused_stc_err << std::endl;
    if (fused_stc_err > 1e-10) std::exit(1);
}

void test_lvl6_factorized_stage_shape()
{
    using L = TFHEpp::lvl6param;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<L>;
    static_assert(Schedule::input_log_q == 48);
    static_assert(Schedule::boot_log_q == 880);
    static_assert(Schedule::raw_linear_stage_count == 14);
    static_assert(Schedule::coeff_to_slot_level_count == 4);
    static_assert(Schedule::slot_to_coeff_level_count == 4);
    static_assert(Schedule::evalmod_polynomial_depth == 7);
    static_assert(Schedule::evalmod_depth == 10);
    static_assert(Schedule::after_coeff_to_slot_log_q == 800);
    static_assert(Schedule::after_evalmod_log_q == 400);
    static_assert(Schedule::output_log_q == 320);
    static_assert(Schedule::message_ratio == 256.0);
    static_assert(Schedule::coeff_to_slot_scaling_factor == 1.0 / 4096.0);
    static_assert(Schedule::slot_to_coeff_scaling_factor == 256.0);
    static_assert(Schedule::OutputCiphertext::log_q == Schedule::output_log_q);

    TFHEpp::CKKSLinearTransformStages<L> c2s_stages;
    TFHEpp::CKKSBuildCoeffToPackedSlotStages<L>(c2s_stages);
    if (c2s_stages.size() != L::nbit - 1) std::exit(1);

    std::size_t c2s_diag_count = 0;
    for (const auto &stage : c2s_stages) {
        if (stage.rotation_offsets.size() > 3) std::exit(1);
        c2s_diag_count += stage.rotation_offsets.size();
    }

    TFHEpp::CKKSLinearTransformStages<L> stc_stages;
    TFHEpp::CKKSBuildPackedSlotToCoeffStages<L>(stc_stages);
    if (stc_stages.size() != L::nbit - 1) std::exit(1);

    std::size_t stc_diag_count = 0;
    for (const auto &stage : stc_stages) {
        if (stage.rotation_offsets.size() > 3) std::exit(1);
        stc_diag_count += stage.rotation_offsets.size();
    }

    std::cout << "CKKS lvl6 factorized C2S/STC stages="
              << c2s_stages.size() << "/" << stc_stages.size()
              << " diagonals=" << c2s_diag_count << "/" << stc_diag_count
              << std::endl;

    TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    TFHEpp::CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    if (linear_plan.coeff_to_slot_stages.size() !=
        Schedule::coeff_to_slot_level_count)
        std::exit(1);
    if (linear_plan.slot_to_coeff_stages.size() !=
        Schedule::slot_to_coeff_level_count)
        std::exit(1);

    std::size_t fused_c2s_diag_count = 0;
    for (const auto &stage : linear_plan.coeff_to_slot_stages) {
        if (stage.rotation_offsets.size() > 81) std::exit(1);
        fused_c2s_diag_count += stage.rotation_offsets.size();
    }
    std::size_t fused_stc_diag_count = 0;
    for (const auto &stage : linear_plan.slot_to_coeff_stages) {
        if (stage.rotation_offsets.size() > 81) std::exit(1);
        fused_stc_diag_count += stage.rotation_offsets.size();
    }
    std::cout << "CKKS lvl6 dense bootstrap fused C2S/STC levels="
              << linear_plan.coeff_to_slot_stages.size() << "/"
              << linear_plan.slot_to_coeff_stages.size()
              << " diagonals=" << fused_c2s_diag_count << "/"
              << fused_stc_diag_count << " output_log_q="
              << Schedule::output_log_q << std::endl;
}

void test_bounded_cos_evalmod_plain()
{
    using L = TFHEpp::lvl6param;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<L>;
    const auto poly = TFHEpp::CKKSBuildBoundedCosEvalModPolynomial<Schedule>();
    constexpr double pi = 3.141592653589793238462643383279502884;

    if (poly.chebyshev_coeffs.size() != Schedule::evalmod_degree + 1)
        std::exit(1);
    if (poly.power_coeffs.size() != Schedule::evalmod_degree + 1)
        std::exit(1);

    double max_basis_err = 0.0;
    for (int i = -32; i <= 32; i++) {
        const double x = static_cast<double>(i) / 32.0;
        const double cheb =
            TFHEpp::CKKSEvaluateChebyshevUnit(poly.chebyshev_coeffs, x);
        const double power =
            TFHEpp::CKKSEvaluatePowerPolynomial(poly.power_coeffs, x);
        max_basis_err = std::max(max_basis_err, std::abs(cheb - power));
    }

    double max_sine_err = 0.0;
    double max_message_err = 0.0;
    double max_power_err = 0.0;
    for (int mask = -static_cast<int>(Schedule::evalmod_k) + 1;
         mask < static_cast<int>(Schedule::evalmod_k); mask++) {
        for (const double msg : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
            const double masked =
                static_cast<double>(mask) * Schedule::message_ratio + msg;
            const double got =
                TFHEpp::CKKSPlainEvalModBoundedCos(poly, masked);
            const double got_power =
                TFHEpp::CKKSPlainEvalModBoundedCosPower(poly, masked);
            const double want_sine =
                Schedule::message_ratio *
                std::sin(2.0 * pi * masked / Schedule::message_ratio) /
                (2.0 * pi);
            max_power_err = std::max(max_power_err, std::abs(got - got_power));
            max_sine_err = std::max(max_sine_err, std::abs(got - want_sine));
            max_message_err =
                std::max(max_message_err, std::abs(got - msg));
        }
    }

    std::cout << "CKKS bounded cosine EvalMod sine max_error="
              << max_sine_err << " message max_error=" << max_message_err
              << " basis max_error=" << max_basis_err
              << " power max_error=" << max_power_err << std::endl;
    if (max_basis_err > 1e-8) std::exit(1);
    if (max_power_err > 1e-8) std::exit(1);
    if (max_sine_err > 1e-8) std::exit(1);
    if (max_message_err > 2e-4) std::exit(1);
}

void test_multilimb_slot_decode_high_level()
{
    using M = SmallMultiLimbCKKSParam;
    constexpr std::uint32_t log_q = 160;
    constexpr std::uint32_t log_delta = 50;

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    slots->fill({0.0, 0.0});
    for (std::size_t i = 0; i < M::n / 2; i++) {
        const double re =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 16.0;
        const double im =
            static_cast<double>(static_cast<int>(i % 3) - 1) / 32.0;
        (*slots)[i] = {re, im};
    }

    TFHEpp::Polynomial<M> poly;
    TFHEpp::ckksSlotEncode<M, log_q, log_delta>(poly, *slots);
    TFHEpp::ckksSlotDecode<M, log_q, log_delta>(*decoded, poly);
    const double err = max_error<M>(*decoded, *slots);
    std::cout << "CKKS multi-limb high-level slot decode max_error=" << err
              << std::endl;
    if (err > 1e-10) std::exit(1);

    const auto neg = TFHEpp::ckksEncodeCoeff<M, log_q, log_delta>(-0.125);
    const double neg_decoded =
        TFHEpp::ckksDecodeCoeff<M, log_q, log_delta>(neg);
    if (std::abs(neg_decoded + 0.125) > 1e-12) std::exit(1);
}

void test_sine_evalmod(const TFHEpp::Key<P> &key)
{
    constexpr std::uint32_t eval_log_q = 128;
    constexpr std::uint32_t eval_log_delta = 36;
    constexpr std::uint32_t coeff_log_delta = 12;
    using EvalCt = TFHEpp::CKKSCiphertext<P, eval_log_q, eval_log_delta>;
    using EvalOut =
        TFHEpp::CKKSEvalModSineDegree3Result<P, eval_log_q, eval_log_delta,
                                             coeff_log_delta>;
    static_assert(EvalOut::log_q == 44);

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    slots->fill({0.0, 0.0});
    expected->fill({0.0, 0.0});
    for (std::size_t i = 0; i < 64; i++) {
        const double x =
            static_cast<double>(static_cast<int>(i % 17) - 8) / 512.0;
        (*slots)[i] = {x, 0.0};
        (*expected)[i] = {TFHEpp::CKKSPlainEvalModSineDegree3(x), 0.0};
    }

    auto keys = std::make_unique<
        TFHEpp::CKKSEvalModSineDegree3RelinKeys<P, eval_log_q,
                                                eval_log_delta,
                                                coeff_log_delta>>();
    TFHEpp::CKKSEvalModSineDegree3KeyGen<P, eval_log_q, eval_log_delta,
                                         coeff_log_delta>(*keys, key,
                                                          {0.0, 0});
    using EvalRelinChain =
        TFHEpp::CKKSRelinKeyChain<P, eval_log_q, eval_log_delta, 2>;
    static_assert(std::tuple_size<typename EvalRelinChain::Tuple>::value == 2);
    auto relin_chain = std::make_unique<EvalRelinChain>();
    TFHEpp::CKKSRelinKeyChainGen<P, eval_log_q, eval_log_delta, 2>(
        *relin_chain, key, {0.0, 0});

    auto ct = std::make_unique<EvalCt>();
    TFHEpp::ckksSlotEncrypt<P, eval_log_q, eval_log_delta>(*ct, *slots, key,
                                                           {0.0, 0});

    using X2 = TFHEpp::CKKSMultResult<P, eval_log_q, eval_log_delta, eval_log_q,
                                      eval_log_delta>;
    using X3 = TFHEpp::CKKSMultResult<P, X2::log_q, X2::log_delta, eval_log_q,
                                      eval_log_delta>;

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto power_expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, eval_log_q, eval_log_delta>(*decoded, *ct, key);
    require_close(*decoded, *slots, 0.01, "CKKS EvalMod fresh input");

    auto transparent = std::make_unique<EvalCt>();
    TFHEpp::CKKSSetTransparentReal<P, eval_log_q, eval_log_delta>(*transparent,
                                                                  -0.125);
    TFHEpp::ckksSlotDecrypt<P, eval_log_q, eval_log_delta>(*decoded,
                                                           *transparent, key);
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*power_expected)[i] = {-0.125, 0.0};
    require_close(*decoded, *power_expected, 1e-9, "CKKS transparent real");

    auto shifted = std::make_unique<EvalCt>(*ct);
    TFHEpp::CKKSAddPlainRealInPlace<P, eval_log_q, eval_log_delta>(*shifted,
                                                                   0.25);
    TFHEpp::CKKSSubPlainRealInPlace<P, eval_log_q, eval_log_delta>(*shifted,
                                                                   0.125);
    TFHEpp::ckksSlotDecrypt<P, eval_log_q, eval_log_delta>(*decoded, *shifted,
                                                           key);
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*power_expected)[i] = (*slots)[i] + std::complex<double>{0.125, 0.0};
    require_close(*decoded, *power_expected, 0.01, "CKKS add/sub plain real");

    auto plain_ct = std::make_unique<EvalCt>();
    plain_ct->ct[0].fill(0);
    TFHEpp::ckksSlotEncode<P, eval_log_q, eval_log_delta>(plain_ct->ct[1],
                                                          *slots);

    auto x2 = std::make_unique<X2>();
    TFHEpp::CKKSMult<P>(*x2, *plain_ct, *plain_ct, relin_chain->get<0>());
    TFHEpp::ckksSlotDecrypt<P, X2::log_q, X2::log_delta>(*decoded, *x2, key);
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*power_expected)[i] = (*slots)[i] * (*slots)[i];
    require_close(*decoded, *power_expected, 0.01,
                  "CKKS EvalMod zero-mask x^2");

    TFHEpp::CKKSMult<P>(*x2, *ct, *ct, relin_chain->get<0>());
    TFHEpp::ckksSlotDecrypt<P, X2::log_q, X2::log_delta>(*decoded, *x2, key);
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*power_expected)[i] = (*slots)[i] * (*slots)[i];
    require_close(*decoded, *power_expected, 0.01, "CKKS EvalMod x^2");

    auto x3 = std::make_unique<X3>();
    TFHEpp::CKKSMult<P>(*x3, *x2, *ct, relin_chain->get<1>());
    TFHEpp::ckksSlotDecrypt<P, X3::log_q, X3::log_delta>(*decoded, *x3, key);
    for (std::size_t i = 0; i < P::n / 2; i++)
        (*power_expected)[i] = (*slots)[i] * (*slots)[i] * (*slots)[i];
    require_close(*decoded, *power_expected, 0.01, "CKKS EvalMod x^3");

    auto out = std::make_unique<EvalOut>();
    TFHEpp::CKKSEvalModSineDegree3<P, eval_log_q, eval_log_delta,
                                   coeff_log_delta>(*out, *ct, *keys);

    TFHEpp::ckksSlotDecrypt<P, EvalOut::log_q, EvalOut::log_delta>(
        *decoded, *out, key);
    require_close(*decoded, *expected, 0.03, "CKKS sine EvalMod degree-3");
}

void test_power_polynomial_evaluator(const TFHEpp::Key<P> &key)
{
    constexpr std::uint32_t log_q = 128;
    constexpr std::uint32_t log_delta = 36;
    constexpr std::uint32_t coeff_log_delta = 12;
    constexpr std::size_t degree = 3;
    using EvalCt = TFHEpp::CKKSCiphertext<P, log_q, log_delta>;
    using Traits =
        TFHEpp::CKKSPowerPolynomialEvaluatorTraits<P, log_q, log_delta,
                                                   coeff_log_delta, degree>;
    using EvalOut =
        TFHEpp::CKKSPowerPolynomialResult<P, log_q, log_delta, coeff_log_delta,
                                          degree>;
    static_assert(Traits::power_depth == 2);
    static_assert(EvalOut::log_q == log_q - 2 * log_delta - coeff_log_delta);

    std::vector<double> coeffs(degree + 1, 0.0);
    coeffs[0] = 0.125;
    coeffs[1] = 1.0;
    coeffs[2] = -0.5;
    coeffs[3] = 0.25;

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    slots->fill({0.0, 0.0});
    expected->fill({0.0, 0.0});
    for (std::size_t i = 0; i < P::n / 2; i++) {
        const double x =
            static_cast<double>(static_cast<int>(i % 17) - 8) / 64.0;
        (*slots)[i] = {x, 0.0};
        (*expected)[i] = {TFHEpp::CKKSEvaluatePowerPolynomial(coeffs, x), 0.0};
    }

    auto relin_chain = std::make_unique<typename Traits::RelinKeyChain>();
    TFHEpp::CKKSRelinKeyChainGen<P, log_q, log_delta, Traits::power_depth>(
        *relin_chain, key, {0.0, 0});

    auto ct = std::make_unique<EvalCt>();
    TFHEpp::ckksSlotEncrypt<P, log_q, log_delta>(*ct, *slots, key, {0.0, 0});

    auto out = std::make_unique<EvalOut>();
    TFHEpp::CKKSEvalPowerPolynomial<P, log_q, log_delta, coeff_log_delta,
                                    degree>(*out, *ct, coeffs, *relin_chain);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    TFHEpp::ckksSlotDecrypt<P, EvalOut::log_q, EvalOut::log_delta>(
        *decoded, *out, key);
    require_close(*decoded, *expected, 0.05, "CKKS power polynomial degree-3");
}

}  // namespace

int main()
{
    test_dense_coeff_slot_diagonals();
    test_factorized_coeff_slot_stages();
    test_lvl6_factorized_stage_shape();
    test_bounded_cos_evalmod_plain();
    test_multilimb_slot_decode_high_level();

    constexpr std::uint32_t low_log_q = 82;
    constexpr std::uint32_t boot_log_q = 110;
    constexpr std::uint32_t log_delta = 40;
    constexpr std::uint32_t plain_log_delta = 20;
    constexpr double tol = 0.05;

    using LowCt = TFHEpp::CKKSCiphertext<P, low_log_q, log_delta>;
    using BootCt = TFHEpp::CKKSCiphertext<P, boot_log_q, log_delta>;
    using LinearOut =
        TFHEpp::CKKSPlainMulResult<P, boot_log_q, log_delta, plain_log_delta>;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_key(*key);
    test_sine_evalmod(*key);
    test_power_polynomial_evaluator(*key);

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    slots->fill({0.0, 0.0});
    for (std::size_t i = 0; i < 32; i++) {
        const double re = static_cast<double>(static_cast<int>(i % 9) - 4) /
                          128.0;
        const double im = static_cast<double>(static_cast<int>(i % 5) - 2) /
                          128.0;
        (*slots)[i] = {re, im};
    }
    (*slots)[P::n / 2 - 1] = {-0.0625, 0.03125};

    auto low_ct = std::make_unique<LowCt>();
    TFHEpp::ckksSlotEncrypt<P, low_log_q, log_delta>(*low_ct, *slots, *key,
                                                     {0.0, 0});

    auto raised = std::make_unique<BootCt>();
    TFHEpp::CKKSModRaise<P, low_log_q, boot_log_q, log_delta>(*raised,
                                                              *low_ct);

    auto low_phase = std::make_unique<TFHEpp::Polynomial<P>>(
        TFHEpp::trlwePhase<P>(low_ct->ct, *key));
    auto raised_phase = std::make_unique<TFHEpp::Polynomial<P>>(
        TFHEpp::trlwePhase<P>(raised->ct, *key));
    for (std::uint32_t i = 0; i < P::n; i++) {
        const auto low =
            TFHEpp::ckks_detail::reduceToLevel<P, low_log_q>((*low_phase)[i]);
        const auto raised_low = TFHEpp::ckks_detail::reduceToLevel<P, low_log_q>(
            (*raised_phase)[i]);
        if (low != raised_low) {
            std::cerr << "mod raise phase mismatch at coeff " << i << std::endl;
            return 1;
        }
    }
    std::cout << "CKKS mod raise preserves low-level phase" << std::endl;

    auto randomized_raised = std::make_unique<BootCt>();
    TFHEpp::CKKSModRaiseRandomized<P, low_log_q, boot_log_q, log_delta>(
        *randomized_raised, *low_ct);
    auto randomized_phase = std::make_unique<TFHEpp::Polynomial<P>>(
        TFHEpp::trlwePhase<P>(randomized_raised->ct, *key));
    bool randomized_has_high_bits = false;
    for (int c = 0; c <= static_cast<int>(P::k); c++) {
        for (std::uint32_t i = 0; i < P::n; i++) {
            const auto low_coeff =
                TFHEpp::ckks_detail::reduceToLevel<P, low_log_q>(
                    low_ct->ct[c][i]);
            const auto randomized_low =
                TFHEpp::ckks_detail::reduceToLevel<P, low_log_q>(
                    randomized_raised->ct[c][i]);
            if (low_coeff != randomized_low) {
                std::cerr << "randomized mod raise coefficient low bits mismatch"
                          << std::endl;
                return 1;
            }
            if (randomized_raised->ct[c][i] != randomized_low)
                randomized_has_high_bits = true;
        }
    }
    for (std::uint32_t i = 0; i < P::n; i++) {
        const auto low =
            TFHEpp::ckks_detail::reduceToLevel<P, low_log_q>((*low_phase)[i]);
        const auto randomized_low =
            TFHEpp::ckks_detail::reduceToLevel<P, low_log_q>(
                (*randomized_phase)[i]);
        if (low != randomized_low) {
            std::cerr << "randomized mod raise phase mismatch at coeff " << i
                      << std::endl;
            return 1;
        }
    }
    if (!randomized_has_high_bits) {
        std::cerr << "randomized mod raise did not populate high bits"
                  << std::endl;
        return 1;
    }
    std::cout << "CKKS randomized mod raise preserves low-level phase"
              << std::endl;

    auto boot_ct = std::make_unique<BootCt>();
    TFHEpp::ckksSlotEncrypt<P, boot_log_q, log_delta>(*boot_ct, *slots, *key,
                                                      {0.0, 0});

    auto boot_gk = std::make_unique<TFHEpp::CKKSGaloisKey<P, boot_log_q>>();
    TFHEpp::CKKSGaloisKeyGen<P, boot_log_q>(*boot_gk, *key, {0.0, 0});
    auto out_gk =
        std::make_unique<TFHEpp::CKKSGaloisKey<P, LinearOut::log_q>>();
    TFHEpp::CKKSGaloisKeyGen<P, LinearOut::log_q>(*out_gk, *key, {0.0, 0});

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    {
        constexpr std::size_t stage_count = 2;
        constexpr int stage_k_step = 4;
        using StageOut =
            TFHEpp::CKKSStagedPlainMulResult<P, boot_log_q, log_delta,
                                             plain_log_delta, stage_count>;

        TFHEpp::CKKSLinearTransformStages<P> c2s_stages;
        TFHEpp::CKKSBuildCoeffToPackedSlotStages<P>(c2s_stages);

        auto stage_gks = std::make_unique<
            TFHEpp::CKKSGaloisKeyChain<P, boot_log_q, plain_log_delta,
                                       stage_count>>();
        TFHEpp::CKKSGaloisKeyChainGen<P, boot_log_q, plain_log_delta,
                                      stage_count>(*stage_gks, *key,
                                                   {0.0, 0});

        auto staged = std::make_unique<StageOut>();
        TFHEpp::CKKSLinearTransformStagesBSGS<P, boot_log_q, log_delta,
                                              plain_log_delta, stage_count>(
            *staged, *boot_ct, c2s_stages, 0, stage_k_step, *stage_gks);
        TFHEpp::ckksSlotDecrypt<P, StageOut::log_q, StageOut::log_delta>(
            *decoded, *staged, *key);

        TFHEpp::CKKSLinearTransformStages<P> prefix_stages(
            c2s_stages.begin(), c2s_stages.begin() + stage_count);
        apply_complex_stages<P>(*expected, *slots, prefix_stages);
        require_close(*decoded, *expected, tol,
                      "CKKS homomorphic factorized C2S prefix");
    }

    {
        auto conjugated = std::make_unique<BootCt>();
        TFHEpp::CKKSConjugateSlots<P, boot_log_q>(conjugated->ct, boot_ct->ct,
                                                  *boot_gk);
        TFHEpp::ckksSlotDecrypt<P, boot_log_q, log_delta>(*decoded,
                                                          *conjugated, *key);
        for (std::size_t i = 0; i < P::n / 2; i++)
            (*expected)[i] = std::conj((*slots)[i]);
        require_close(*decoded, *expected, tol, "CKKS conjugation");
    }

    std::vector<TFHEpp::CKKSSlotVector<P>> diagonals(4);
    for (auto &d : diagonals) d.fill({0.0, 0.0});
    for (std::size_t i = 0; i < P::n / 2; i++) {
        diagonals[0][i] = {0.25, 0.0};
        diagonals[1][i] = {0.125, 0.0625};
        diagonals[2][i] = {-0.0625, 0.03125};
        diagonals[3][i] = {i % 2 == 0 ? 0.03125 : -0.03125, -0.015625};
    }
    const std::vector<int> offsets{0, 1, 5, -3};
    constexpr int k_step = 4;
    TFHEpp::CKKSLinearTransformPlan<P, boot_log_q, log_delta, plain_log_delta>
        plan;
    TFHEpp::CKKSBuildLinearTransformBSGSPlan<P, boot_log_q, log_delta,
                                             plain_log_delta>(
        plan, diagonals, offsets, k_step);

    auto transformed = std::make_unique<LinearOut>();
    TFHEpp::CKKSLinearTransformBSGS<P, boot_log_q, log_delta,
                                    plain_log_delta>(
        *transformed, *boot_ct, plan, *boot_gk, *out_gk);
    TFHEpp::ckksSlotDecrypt<P, LinearOut::log_q, LinearOut::log_delta>(
        *decoded, *transformed, *key);

    constexpr int half = static_cast<int>(P::n) / 2;
    for (int i = 0; i < half; i++) {
        (*expected)[i] =
            diagonals[0][i] * (*slots)[i] +
            diagonals[1][i] * (*slots)[(i + 1) % half] +
            diagonals[2][i] * (*slots)[(i + 5) % half] +
            diagonals[3][i] * (*slots)[(i - 3 + half) % half];
    }
    require_close(*decoded, *expected, tol, "CKKS precomputed BSGS");

    std::vector<TFHEpp::CKKSSlotVector<P>> direct_real_diags(1);
    std::vector<TFHEpp::CKKSSlotVector<P>> conj_real_diags(1);
    direct_real_diags[0].fill({1.0, 0.0});
    conj_real_diags[0].fill({1.0, 0.0});
    const std::vector<int> zero_offset{0};

    TFHEpp::CKKSRealLinearTransformPlan<P, boot_log_q, log_delta,
                                        plain_log_delta>
        real_plan;
    TFHEpp::CKKSBuildRealLinearTransformPlan<P, boot_log_q, log_delta,
                                             plain_log_delta>(
        real_plan, direct_real_diags, zero_offset, conj_real_diags, zero_offset,
        1);
    auto real_part = std::make_unique<LinearOut>();
    TFHEpp::CKKSRealLinearTransform<P, boot_log_q, log_delta, plain_log_delta>(
        *real_part, *boot_ct, real_plan, *boot_gk, *out_gk);
    TFHEpp::ckksSlotDecrypt<P, LinearOut::log_q, LinearOut::log_delta>(
        *decoded, *real_part, *key);
    for (int i = 0; i < half; i++) (*expected)[i] = {2 * (*slots)[i].real(), 0};
    require_close(*decoded, *expected, tol, "CKKS real-linear transform");

    std::cout << "Passed" << std::endl;
}
