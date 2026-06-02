#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <tfhe++.hpp>
#include <tuple>
#include <type_traits>

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

struct TinyMultiLimbCKKSParam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 4;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 1;
    static constexpr std::uint32_t lₐ = 1;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t Bgₐbit = 16;
    using T = TFHEpp::MultiLimbUInt<5>;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bgₐ = T{1} << Bgₐbit;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = 0.0;
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 20;
    static constexpr std::uint32_t l̅ₐ = 20;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;
};

struct TinyDeepMultiLimbCKKSParam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 4;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 1;
    static constexpr std::uint32_t lₐ = 1;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t Bgₐbit = 16;
    using T = TFHEpp::MultiLimbUInt<9>;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bgₐ = T{1} << Bgₐbit;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = 0.0;
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 36;
    static constexpr std::uint32_t l̅ₐ = 36;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;
};

void fill_key(TFHEpp::Key<P> &key)
{
    for (std::size_t i = 0; i < P::n; i++) {
        const int v = static_cast<int>(i % 3) - 1;
        key[i] = static_cast<typename P::T>(v);
    }
}

template <class Param>
void fill_test_key(TFHEpp::Key<Param> &key)
{
    for (std::size_t i = 0; i < Param::n; i++) {
        const int v = static_cast<int>(i % 3) - 1;
        key[i] = static_cast<typename Param::T>(v);
    }
}

template <class Param>
bool trlwe_equal(const TFHEpp::TRLWE<Param> &lhs,
                 const TFHEpp::TRLWE<Param> &rhs)
{
    for (int c = 0; c <= static_cast<int>(Param::k); c++)
        for (std::uint32_t i = 0; i < Param::n; i++)
            if (lhs[static_cast<std::size_t>(c)][i] !=
                rhs[static_cast<std::size_t>(c)][i])
                return false;
    return true;
}

void test_seeded_key_switch_rows()
{
    using Param = TinyDeepMultiLimbCKKSParam;
    constexpr std::uint32_t log_q = 95;
    constexpr std::uint32_t row_count =
        TFHEpp::CKKSKeySwitchRowCountForLevel<Param, log_q>();
    static_assert(row_count > 1);
    static_assert(TFHEpp::CKKSSeededKeySwitchRowByteSize<Param>() <
                  TFHEpp::CKKSKeySwitchRowByteSize<Param>());

    auto key = std::make_unique<TFHEpp::Key<Param>>();
    fill_test_key<Param>(*key);

    using SeededRow = TFHEpp::CKKSSeededKeySwitchRow<Param, log_q>;
    std::array<SeededRow, row_count> seeded_rows{};
    std::array<TFHEpp::TRLWE<Param>, row_count> expanded_rows{};

    for (std::uint32_t j = 0; j < row_count; j++) {
        TFHEpp::Polynomial<Param> gadget{};
        for (std::uint32_t i = 0; i < Param::n; i++) {
            const auto value = static_cast<__int128_t>((j + 3) * (i + 5)) - 41;
            gadget[i] = TFHEpp::ckks_detail::signedToLevel<Param, log_q>(value);
        }
        TFHEpp::ckks_detail::encryptPolynomialAtLevel<Param, log_q>(
            seeded_rows[j], gadget, *key, {0.0, 0});
        TFHEpp::ckks_detail::expandSeededKeySwitchRow<Param, log_q>(
            expanded_rows[j], seeded_rows[j]);

        TFHEpp::TRLWE<Param> expanded_again{};
        TFHEpp::ckks_detail::expandSeededKeySwitchRow<Param, log_q>(
            expanded_again, seeded_rows[j]);
        if (!trlwe_equal<Param>(expanded_rows[j], expanded_again)) {
            std::cerr << "seeded key-switch row expansion is not deterministic"
                      << std::endl;
            std::exit(1);
        }
    }

    TFHEpp::Polynomial<Param> input{};
    for (std::uint32_t i = 0; i < Param::n; i++) {
        const auto value = static_cast<__int128_t>(
            static_cast<int>(i % 7) - 3);
        input[i] = TFHEpp::ckks_detail::signedToLevel<Param, log_q>(value);
    }

    TFHEpp::TRLWE<Param> expanded_out{};
    TFHEpp::TRLWE<Param> seeded_out{};
    TFHEpp::CKKSKeySwitchRows<Param, log_q>(expanded_out, input,
                                            expanded_rows);
    TFHEpp::CKKSKeySwitchRows<Param, log_q>(seeded_out, input, seeded_rows);
    if (!trlwe_equal<Param>(expanded_out, seeded_out)) {
        std::cerr << "seeded key-switch rows differ from expanded rows"
                  << std::endl;
        std::exit(1);
    }
    std::cout << "CKKS seeded key-switch rows match expanded rows" << std::endl;
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
void require_close_param(const TFHEpp::CKKSSlotVector<Param> &got,
                         const TFHEpp::CKKSSlotVector<Param> &want, double tol,
                         const char *label)
{
    const double err = max_error<Param>(got, want);
    std::cout << label << " max_error=" << err << std::endl;
    if (err > tol) {
        std::cerr << label << " exceeded tolerance " << tol << std::endl;
        std::exit(1);
    }
}

template <class Schedule>
struct CountingDenseBootstrapProvider {
    using EvalModTraits = TFHEpp::CKKSDenseEvalModBoundedCosTraits<Schedule>;

    const TFHEpp::CKKSDenseBootstrapKey<Schedule> &key;
    mutable std::array<std::size_t, Schedule::coeff_to_slot_level_count + 1>
        coeff_to_slot_gets{};
    mutable std::array<std::size_t, Schedule::coeff_to_slot_level_count + 1>
        coeff_to_slot_releases{};
    mutable std::array<std::size_t, Schedule::slot_to_coeff_level_count + 1>
        slot_to_coeff_gets{};
    mutable std::array<std::size_t, Schedule::slot_to_coeff_level_count + 1>
        slot_to_coeff_releases{};
    mutable std::array<std::size_t,
                       EvalModTraits::PolynomialTraits::power_depth>
        polynomial_relin_gets{};
    mutable std::array<std::size_t,
                       EvalModTraits::PolynomialTraits::power_depth>
        polynomial_relin_releases{};
    mutable std::array<std::size_t, Schedule::evalmod_double_angle>
        double_angle_relin_gets{};
    mutable std::array<std::size_t, Schedule::evalmod_double_angle>
        double_angle_relin_releases{};
    mutable std::array<std::size_t,
                       EvalModTraits::InverseTraits::power_depth>
        inverse_relin_gets{};
    mutable std::array<std::size_t,
                       EvalModTraits::InverseTraits::power_depth>
        inverse_relin_releases{};
    mutable std::size_t packed_conjugate_gets = 0;
    mutable std::size_t packed_conjugate_releases = 0;
    mutable std::size_t evalmod_polynomial_gets = 0;

    explicit CountingDenseBootstrapProvider(
        const TFHEpp::CKKSDenseBootstrapKey<Schedule> &bootstrap_key)
        : key(bootstrap_key)
    {}

    const TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan() const
    {
        return key.linear_plan;
    }

    const TFHEpp::CKKSBoundedCosEvalModPolynomial &evalmod_polynomial() const
    {
        evalmod_polynomial_gets++;
        return key.evalmod_polynomial;
    }

    template <std::size_t I>
    const auto &coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        coeff_to_slot_gets[I]++;
        return key.coeff_to_slot_galois.template get<I>();
    }

    template <std::size_t I>
    void release_coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        coeff_to_slot_releases[I]++;
    }

    const auto &packed_conjugate_galois() const
    {
        packed_conjugate_gets++;
        return key.packed_conjugate_galois;
    }

    void release_packed_conjugate_galois() const
    {
        packed_conjugate_releases++;
    }

    template <std::size_t I>
    const auto &polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::power_depth);
        polynomial_relin_gets[I]++;
        return key.evalmod_relin.polynomial.template get<I>();
    }

    template <std::size_t I>
    void release_polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::power_depth);
        polynomial_relin_releases[I]++;
    }

    template <std::size_t I>
    const auto &double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        double_angle_relin_gets[I]++;
        return key.evalmod_relin.double_angle.template get<I>();
    }

    template <std::size_t I>
    void release_double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        double_angle_relin_releases[I]++;
    }

    template <std::size_t I>
    const auto &inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        inverse_relin_gets[I]++;
        return key.evalmod_relin.inverse.template get<I>();
    }

    template <std::size_t I>
    void release_inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        inverse_relin_releases[I]++;
    }

    template <std::size_t I>
    const auto &slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        slot_to_coeff_gets[I]++;
        return key.slot_to_coeff_galois.template get<I>();
    }

    template <std::size_t I>
    void release_slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        slot_to_coeff_releases[I]++;
    }
};

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

template <class Param>
void coeffs_to_slots(TFHEpp::CKKSSlotVector<Param> &slots,
                     const std::array<double, Param::n> &coeffs)
{
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const auto &slot_to_eval = TFHEpp::ckks_detail::ckksSlotToEvalIndex<Param>();
    for (std::size_t i = 0; i < Param::n / 2; i++) {
        const long double h =
            static_cast<long double>(2 * slot_to_eval[i] + 1);
        std::complex<long double> sum = 0.0L;
        for (std::size_t j = 0; j < Param::n; j++) {
            const long double angle =
                pi * h * static_cast<long double>(j) /
                static_cast<long double>(Param::n);
            sum += static_cast<long double>(coeffs[j]) *
                   std::complex<long double>(std::cos(angle),
                                             std::sin(angle));
        }
        slots[i] = static_cast<std::complex<double>>(sum);
    }
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
    static_assert(Schedule::input_log_q == 58);
    static_assert(Schedule::boot_log_q == 880);
    static_assert(Schedule::linear_plain_log_delta == 50);
    static_assert(Schedule::coeff_to_slot_plain_log_delta == 50);
    static_assert(Schedule::component_split_plain_log_delta == 50);
    static_assert(Schedule::slot_to_coeff_plain_log_delta == 35);
    static_assert(Schedule::coeff_to_slot_fuse_radix == 5);
    static_assert(Schedule::slot_to_coeff_fuse_radix == 5);
    static_assert(Schedule::linear_bsgs_step == 128);
    static_assert(Schedule::hybrid_giant_direct_popcount_threshold == 4);
    static_assert(Schedule::raw_linear_stage_count == 14);
    static_assert(Schedule::coeff_to_slot_level_count == 3);
    static_assert(Schedule::slot_to_coeff_level_count == 3);
    static_assert(Schedule::evalmod_polynomial_depth == 6);
    static_assert(Schedule::evalmod_depth == 9);
    static_assert(Schedule::modraise_mask_bound == 0);
    static_assert(Schedule::after_coeff_to_slot_log_q == 730);
    static_assert(Schedule::after_component_split_log_q == 680);
    static_assert(Schedule::after_evalmod_log_q == 230);
    static_assert(Schedule::output_log_q == 125);
    static_assert(Schedule::message_ratio == 256.0);
    static_assert(Schedule::coeff_to_slot_scaling_factor == 1.0 / 4096.0);
    static_assert(Schedule::slot_to_coeff_scaling_factor == 256.0);
    static_assert(Schedule::OutputCiphertext::log_q == Schedule::output_log_q);
    using DenseEvalTraits = TFHEpp::CKKSDenseEvalModBoundedCosTraits<Schedule>;
    using DenseEvalOut = TFHEpp::CKKSDenseEvalModBoundedCosResult<Schedule>;
    static_assert(DenseEvalTraits::polynomial_log_q == 380);
    static_assert(DenseEvalOut::log_q == Schedule::after_evalmod_log_q);

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
    if (linear_plan.slot_to_coeff_imag_stages.size() !=
        Schedule::slot_to_coeff_level_count)
        std::exit(1);

    std::size_t fused_c2s_diag_count = 0;
    for (const auto &stage : linear_plan.coeff_to_slot_stages) {
        if (stage.rotation_offsets.empty()) std::exit(1);
        fused_c2s_diag_count += stage.rotation_offsets.size();
    }
    std::size_t fused_stc_diag_count = 0;
    for (const auto &stage : linear_plan.slot_to_coeff_stages) {
        if (stage.rotation_offsets.empty()) std::exit(1);
        fused_stc_diag_count += stage.rotation_offsets.size();
    }
    std::cout << "CKKS lvl6 dense bootstrap fused C2S/STC levels="
              << linear_plan.coeff_to_slot_stages.size() << "/"
              << linear_plan.slot_to_coeff_stages.size()
              << " diagonals=" << fused_c2s_diag_count << "/"
              << fused_stc_diag_count << " output_log_q="
              << Schedule::output_log_q << std::endl;

    std::size_t total_c2s_baby_rotations = 0;
    std::size_t max_c2s_baby_rotations = 0;
    std::size_t c2s_current_evalautos = 0;
    std::size_t c2s_hybrid_evalautos = 0;
    std::size_t c2s_direct_evalautos = 0;
    for (const auto &stage : linear_plan.coeff_to_slot_stages) {
        const std::size_t rotations =
            TFHEpp::CKKSLinearTransformStageBabyRotationCount<L>(
                stage, Schedule::linear_bsgs_step);
        total_c2s_baby_rotations += rotations;
        max_c2s_baby_rotations = std::max(max_c2s_baby_rotations, rotations);
        c2s_current_evalautos +=
            TFHEpp::CKKSLinearTransformStageRotationEvalAutoCount<L>(
                stage, Schedule::linear_bsgs_step);
        c2s_hybrid_evalautos +=
            TFHEpp::CKKSLinearTransformStageHybridGiantRotationEvalAutoCount<L>(
                stage, Schedule::linear_bsgs_step,
                Schedule::hybrid_giant_direct_popcount_threshold);
        c2s_direct_evalautos +=
            TFHEpp::CKKSLinearTransformStageDirectRotationEvalAutoCount<L>(
                stage, Schedule::linear_bsgs_step);
    }
    if (total_c2s_baby_rotations == 0 ||
        max_c2s_baby_rotations >=
            static_cast<std::size_t>(Schedule::linear_bsgs_step))
        std::exit(1);
    const std::size_t stc_current_evalautos =
        TFHEpp::CKKSLinearTransformStagesDualInputSharedTailRotationEvalAutoCount<
            L>(linear_plan.slot_to_coeff_stages, 0,
               linear_plan.slot_to_coeff_stages.size(),
               Schedule::linear_bsgs_step);
    const std::size_t stc_direct_evalautos =
        TFHEpp::
            CKKSLinearTransformStagesDualInputSharedTailDirectRotationEvalAutoCount<
                L>(linear_plan.slot_to_coeff_stages, 0,
                   linear_plan.slot_to_coeff_stages.size(),
                   Schedule::linear_bsgs_step);
    const std::size_t stc_hybrid_evalautos =
        TFHEpp::
            CKKSLinearTransformStagesDualInputSharedTailHybridGiantRotationEvalAutoCount<
                L>(linear_plan.slot_to_coeff_stages, 0,
                   linear_plan.slot_to_coeff_stages.size(),
                   Schedule::linear_bsgs_step,
                   Schedule::hybrid_giant_direct_popcount_threshold);
    if (c2s_direct_evalautos == 0 ||
        c2s_direct_evalautos >= c2s_current_evalautos)
        std::exit(1);
    if (stc_direct_evalautos == 0 ||
        stc_direct_evalautos >= stc_current_evalautos)
        std::exit(1);
    if (c2s_hybrid_evalautos <= c2s_direct_evalautos ||
        c2s_hybrid_evalautos >= c2s_current_evalautos)
        std::exit(1);
    if (stc_hybrid_evalautos <= stc_direct_evalautos ||
        stc_hybrid_evalautos >= stc_current_evalautos)
        std::exit(1);

    TFHEpp::CKKSDenseBootstrapRotationKeyUsage<Schedule> rotation_usage;
    TFHEpp::CKKSBuildDenseBootstrapRotationKeyUsage<Schedule>(rotation_usage,
                                                              linear_plan);
    const std::size_t planned_key_indices =
        TFHEpp::CKKSDenseBootstrapRotationKeyUsageCount<Schedule>(
            rotation_usage);
    constexpr std::size_t full_key_indices =
        TFHEpp::CKKSDenseBootstrapFullGaloisKeyIndexCount<Schedule>();
    static_assert(Schedule::supports_post_bootstrap_product);
    static_assert(
        TFHEpp::CKKSAutoKeySwitchRowCount<L, Schedule::boot_log_q>() == 55);
    static_assert(
        TFHEpp::CKKSRelinKeySwitchRowCount<L, Schedule::output_log_q>() == 8);
    static_assert(
        TFHEpp::CKKSDenseBootstrapEvalModRelinKeyCount<Schedule>() == 8);
    if (planned_key_indices == 0 || planned_key_indices >= full_key_indices)
        std::exit(1);
    if (planned_key_indices != 57 || full_key_indices != 144)
        std::exit(1);
    if (TFHEpp::CKKSRotationKeyIndexSetCount<L>(
            rotation_usage.packed_conjugate) != 1)
        std::exit(1);
    const std::size_t sparse_galois_rows =
        TFHEpp::CKKSDenseBootstrapSparseGaloisKeySwitchRowCount<Schedule>(
            rotation_usage);
    const std::size_t sparse_key_rows =
        TFHEpp::CKKSDenseBootstrapSparseKeySwitchRowCount<Schedule>(
            rotation_usage);
    const std::size_t sparse_key_bytes =
        TFHEpp::CKKSDenseBootstrapSparseKeyByteEstimate<Schedule>(
            rotation_usage);
    TFHEpp::CKKSDenseBootstrapDirectRotationKeyUsage<Schedule>
        direct_rotation_usage;
    TFHEpp::CKKSBuildDenseBootstrapDirectRotationKeyUsage<Schedule>(
        direct_rotation_usage, linear_plan);
    TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        hybrid_rotation_usage;
    TFHEpp::CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        hybrid_rotation_usage, linear_plan);
    const std::size_t direct_key_indices =
        TFHEpp::CKKSDenseBootstrapDirectRotationKeyUsageCount<Schedule>(
            direct_rotation_usage);
    const std::size_t hybrid_key_indices =
        TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsageCount<Schedule>(
            hybrid_rotation_usage);
    const std::size_t direct_key_rows =
        TFHEpp::CKKSDenseBootstrapDirectKeySwitchRowCount<Schedule>(
            direct_rotation_usage);
    const std::size_t hybrid_key_rows =
        TFHEpp::CKKSDenseBootstrapHybridGiantKeySwitchRowCount<Schedule>(
            hybrid_rotation_usage);
    const std::size_t direct_key_bytes =
        TFHEpp::CKKSDenseBootstrapDirectKeyByteEstimate<Schedule>(
            direct_rotation_usage);
    const std::size_t hybrid_key_bytes =
        TFHEpp::CKKSDenseBootstrapHybridGiantKeyByteEstimate<Schedule>(
            hybrid_rotation_usage);
    const std::size_t direct_streamed_peak_rows =
        TFHEpp::CKKSDenseBootstrapDirectStreamedKeySwitchPeakRowCount<
            Schedule>(direct_rotation_usage);
    const std::size_t hybrid_streamed_peak_rows =
        TFHEpp::CKKSDenseBootstrapHybridGiantStreamedKeySwitchPeakRowCount<
            Schedule>(hybrid_rotation_usage);
    const std::size_t direct_streamed_peak_bytes =
        TFHEpp::CKKSDenseBootstrapDirectStreamedPeakKeyByteEstimate<Schedule>(
            direct_rotation_usage);
    const std::size_t hybrid_streamed_peak_bytes =
        TFHEpp::CKKSDenseBootstrapHybridGiantStreamedPeakKeyByteEstimate<
            Schedule>(hybrid_rotation_usage);
    const std::size_t streamed_peak_rows =
        TFHEpp::CKKSDenseBootstrapStreamedKeySwitchPeakRowCount<Schedule>(
            rotation_usage);
    const std::size_t streamed_peak_bytes =
        TFHEpp::CKKSDenseBootstrapStreamedPeakKeyByteEstimate<Schedule>(
            rotation_usage);
    constexpr std::size_t full_key_rows =
        TFHEpp::CKKSDenseBootstrapFullKeySwitchRowCount<Schedule>();
    constexpr std::size_t full_key_bytes =
        TFHEpp::CKKSDenseBootstrapFullKeyByteEstimate<Schedule>();
    if (sparse_key_rows != sparse_galois_rows +
                               TFHEpp::CKKSDenseBootstrapEvalModKeySwitchRowCount<
                                   Schedule>())
        std::exit(1);
    if (sparse_key_bytes !=
        sparse_key_rows * TFHEpp::CKKSKeySwitchRowByteSize<L>())
        std::exit(1);
    if (streamed_peak_bytes !=
        streamed_peak_rows * TFHEpp::CKKSKeySwitchRowByteSize<L>())
        std::exit(1);
    if (sparse_key_rows >= full_key_rows || sparse_key_bytes >= full_key_bytes)
        std::exit(1);
    if (streamed_peak_rows == 0 || streamed_peak_rows >= sparse_key_rows)
        std::exit(1);
    if (direct_key_indices <= planned_key_indices ||
        direct_key_rows <= sparse_key_rows ||
        direct_key_bytes <= sparse_key_bytes)
        std::exit(1);
    if (direct_streamed_peak_rows == 0 ||
        direct_streamed_peak_bytes !=
            direct_streamed_peak_rows * TFHEpp::CKKSKeySwitchRowByteSize<L>())
        std::exit(1);
    if (hybrid_key_indices == 0 || hybrid_key_indices >= direct_key_indices ||
        hybrid_key_rows == 0 || hybrid_key_rows >= direct_key_rows)
        std::exit(1);
    if (hybrid_key_rows >= sparse_key_rows ||
        hybrid_key_bytes >= sparse_key_bytes)
        std::exit(1);
    if (hybrid_streamed_peak_rows == 0 ||
        hybrid_streamed_peak_bytes !=
            hybrid_streamed_peak_rows * TFHEpp::CKKSKeySwitchRowByteSize<L>())
        std::exit(1);
    std::cout << "CKKS lvl6 dense bootstrap rotation key indices planned/full="
              << planned_key_indices << "/" << full_key_indices
              << " c2s_baby_rotations=" << total_c2s_baby_rotations
              << std::endl;
    std::cout << "CKKS lvl6 dense bootstrap rotation evalautos current/direct "
              << "c2s=" << c2s_current_evalautos << "/"
              << c2s_direct_evalautos << " stc=" << stc_current_evalautos
              << "/" << stc_direct_evalautos << std::endl;
    std::cout << "CKKS lvl6 dense bootstrap rotation evalautos hybrid "
              << "c2s=" << c2s_hybrid_evalautos
              << " stc=" << stc_hybrid_evalautos
              << " direct_popcount_threshold="
              << Schedule::hybrid_giant_direct_popcount_threshold
              << std::endl;
    std::cout << "CKKS lvl6 dense bootstrap key rows sparse/full="
              << sparse_key_rows << "/" << full_key_rows
              << " bytes=" << sparse_key_bytes << "/" << full_key_bytes
              << std::endl;
    std::cout << "CKKS lvl6 dense bootstrap streamed peak rows="
              << streamed_peak_rows << " bytes=" << streamed_peak_bytes
              << std::endl;
    std::cout << "CKKS lvl6 dense bootstrap direct key indices/rows="
              << direct_key_indices << "/" << direct_key_rows
              << " bytes=" << direct_key_bytes
              << " streamed_peak_rows=" << direct_streamed_peak_rows
              << " streamed_peak_bytes=" << direct_streamed_peak_bytes
              << std::endl;
    std::cout << "CKKS lvl6 dense bootstrap hybrid key indices/rows="
              << hybrid_key_indices << "/" << hybrid_key_rows
              << " bytes=" << hybrid_key_bytes
              << " streamed_peak_rows=" << hybrid_streamed_peak_rows
              << " streamed_peak_bytes=" << hybrid_streamed_peak_bytes
              << std::endl;
}

void test_popcount3_direct_rotation()
{
    using M = TinyDeepMultiLimbCKKSParam;
    constexpr std::uint32_t log_q = 520;
    constexpr std::uint32_t log_delta = 40;
    constexpr int steps = 7;

    auto key = std::make_unique<TFHEpp::Key<M>>();
    fill_test_key<M>(*key);

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    for (std::size_t i = 0; i < M::n / 2; i++) {
        (*slots)[i] = {
            static_cast<double>(static_cast<int>(i % 17) - 8) / 512.0,
            static_cast<double>(static_cast<int>((3 * i) % 19) - 9) / 1024.0};
    }

    auto ct = std::make_unique<TFHEpp::CKKSCiphertext<M, log_q, log_delta>>();
    TFHEpp::ckksSlotEncrypt<M, log_q, log_delta>(*ct, *slots, *key, {0.0, 0});

    TFHEpp::CKKSRotationKeyIndexSet<M> binary_indices{};
    TFHEpp::CKKSDirectRotationKeyIndexSet<M> direct_indices{};
    TFHEpp::CKKSClearRotationKeyIndexSet<M>(binary_indices);
    TFHEpp::CKKSClearDirectRotationKeyIndexSet<M>(direct_indices);
    TFHEpp::CKKSMarkRotationPowerKeyIndices<M>(binary_indices, steps);
    TFHEpp::CKKSMarkDirectRotationKeyIndex<M>(direct_indices, steps);

    auto binary_gk =
        std::make_unique<TFHEpp::CKKSSparseGaloisKey<M, log_q>>();
    auto direct_gk =
        std::make_unique<TFHEpp::CKKSDirectSparseGaloisKey<M, log_q>>();
    TFHEpp::CKKSSparseGaloisKeyGen<M, log_q>(*binary_gk, *key,
                                             binary_indices, {M::α, 0});
    TFHEpp::CKKSDirectSparseGaloisKeyGen<M, log_q>(*direct_gk, *key,
                                                   direct_indices, {M::α, 0});

    auto binary_rotated =
        std::make_unique<TFHEpp::CKKSCiphertext<M, log_q, log_delta>>();
    auto direct_rotated =
        std::make_unique<TFHEpp::CKKSCiphertext<M, log_q, log_delta>>();
    TFHEpp::CKKSRotateSlots<M, log_q>(binary_rotated->ct, ct->ct, steps,
                                      *binary_gk);
    TFHEpp::CKKSRotateSlots<M, log_q>(direct_rotated->ct, ct->ct, steps,
                                      *direct_gk);

    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto binary_decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto direct_decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    TFHEpp::rotateCKKSSlotVector<M>(*expected, *slots, steps);
    TFHEpp::ckksSlotDecrypt<M, log_q, log_delta>(*binary_decoded,
                                                 *binary_rotated, *key);
    TFHEpp::ckksSlotDecrypt<M, log_q, log_delta>(*direct_decoded,
                                                 *direct_rotated, *key);

    require_close_param<M>(*binary_decoded, *expected, 1e-6,
                           "CKKS popcount-3 binary rotation");
    require_close_param<M>(*direct_decoded, *expected, 1e-6,
                           "CKKS popcount-3 direct rotation");
    require_close_param<M>(*direct_decoded, *binary_decoded, 1e-6,
                           "CKKS popcount-3 direct/binary rotation");
}

void test_dense_bootstrap_api_shape()
{
    using M = TinyDeepMultiLimbCKKSParam;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        M, 40, 8, 560, 40, 3, 31, 4, 2, 0, 40, 2>;
    static_assert(Schedule::input_log_q == 48);
    static_assert(Schedule::boot_log_q == 560);
    static_assert(Schedule::linear_bsgs_step == 2);
    static_assert(Schedule::coeff_to_slot_level_count == 1);
    static_assert(Schedule::slot_to_coeff_level_count == 1);
    static_assert(Schedule::after_coeff_to_slot_log_q == 520);
    static_assert(Schedule::after_component_split_log_q == 480);
    static_assert(Schedule::evalmod_polynomial_depth == 6);
    static_assert(Schedule::evalmod_depth == 8);
    static_assert(Schedule::modraise_mask_bound == 0);
    static_assert(Schedule::after_evalmod_log_q == 160);
    static_assert(Schedule::output_log_q == 120);

    using FreshProductTraits = TFHEpp::CKKSDenseBootstrapProductTraits<
        Schedule, Schedule::input_log_q + 2 * Schedule::log_delta,
        Schedule::log_delta,
        Schedule::input_log_q + 2 * Schedule::log_delta,
        Schedule::log_delta>;
    using PostBootstrapProductTraits =
        TFHEpp::CKKSDenseBootstrapProductTraits<
            Schedule, Schedule::output_log_q, Schedule::log_delta,
            Schedule::output_log_q, Schedule::log_delta>;
    static_assert(FreshProductTraits::ProductCiphertext::log_q ==
                  Schedule::input_log_q + Schedule::log_delta);
    static_assert(PostBootstrapProductTraits::ProductCiphertext::log_q ==
                  Schedule::post_bootstrap_product_log_q);
    static_assert(FreshProductTraits::product_level_compatible);
    static_assert(PostBootstrapProductTraits::product_level_compatible);

    using Lvl6FastSchedule = TFHEpp::lvl6CKKSDenseBootstrapFastSchedule;
    using Lvl6CompactSchedule = TFHEpp::lvl6CKKSDenseBootstrapCompactSchedule;
    using Lvl6InverseSchedule =
        TFHEpp::lvl6CKKSDenseBootstrapInverseSchedule;
    using Lvl6TunedSchedule = TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule;
    static_assert(std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapSchedule,
                                 Lvl6FastSchedule>);
    static_assert(std::is_same_v<typename Lvl6FastSchedule::Param,
                                 TFHEpp::lvl6param>);
    static_assert(Lvl6FastSchedule::evalmod_k == 18);
    static_assert(Lvl6FastSchedule::evalmod_degree == 34);
    static_assert(Lvl6FastSchedule::boot_log_q == 888);
    static_assert(Lvl6FastSchedule::slot_to_coeff_plain_log_delta == 25);
    static_assert(Lvl6FastSchedule::output_log_q == 113);
    static_assert(Lvl6FastSchedule::supports_post_bootstrap_product);
    static_assert(Lvl6FastSchedule::post_bootstrap_product_slack == 5);
    static_assert(Lvl6FastSchedule::hybrid_giant_direct_popcount_threshold == 3);
    static_assert(Lvl6CompactSchedule::hybrid_giant_direct_popcount_threshold ==
                  4);
    using Lvl6PostBootstrapProductTraits =
        TFHEpp::CKKSDenseBootstrapProductTraits<
            Lvl6FastSchedule, Lvl6FastSchedule::output_log_q,
            Lvl6FastSchedule::log_delta, Lvl6FastSchedule::output_log_q,
            Lvl6FastSchedule::log_delta>;
    static_assert(
        Lvl6PostBootstrapProductTraits::ProductCiphertext::log_q ==
        Lvl6FastSchedule::post_bootstrap_product_log_q);
    static_assert(Lvl6PostBootstrapProductTraits::scale_drop == 0);
    static_assert(Lvl6PostBootstrapProductTraits::product_level_compatible);
    static_assert(std::is_same_v<typename Lvl6InverseSchedule::Param,
                                 TFHEpp::lvl6param>);
    static_assert(Lvl6InverseSchedule::log_delta == 40);
    static_assert(Lvl6InverseSchedule::evalmod_degree == 52);
    static_assert(Lvl6InverseSchedule::evalmod_inv_degree == 3);
    static_assert(Lvl6InverseSchedule::evalmod_log_q_consumption == 560);
    static_assert(Lvl6InverseSchedule::coeff_to_slot_plain_log_delta == 50);
    static_assert(Lvl6InverseSchedule::component_split_plain_log_delta == 50);
    static_assert(Lvl6InverseSchedule::slot_to_coeff_plain_log_delta == 20);
    static_assert(Lvl6InverseSchedule::output_log_q == 60);
    static_assert(std::is_same_v<typename Lvl6TunedSchedule::Param,
                                 TFHEpp::lvl6param>);
    static_assert(Lvl6TunedSchedule::log_delta == 52);
    static_assert(Lvl6TunedSchedule::input_log_q == 60);
    static_assert(Lvl6TunedSchedule::evalmod_degree == 52);
    static_assert(Lvl6TunedSchedule::evalmod_double_angle == 4);
    static_assert(Lvl6TunedSchedule::evalmod_inv_degree == 7);
    static_assert(Lvl6TunedSchedule::evalmod_log_q_consumption == 780);
    static_assert(Lvl6TunedSchedule::coeff_to_slot_plain_log_delta == 52);
    static_assert(Lvl6TunedSchedule::component_split_plain_log_delta == 52);
    static_assert(Lvl6TunedSchedule::slot_to_coeff_plain_log_delta == 30);
    static_assert(Lvl6TunedSchedule::after_evalmod_log_q == 216);
    static_assert(Lvl6TunedSchedule::output_log_q == 156);
    static_assert(Lvl6TunedSchedule::supports_post_bootstrap_product);
    static_assert(Lvl6TunedSchedule::post_bootstrap_product_slack == 44);
    using Lvl6SparseKey = TFHEpp::CKKSDenseBootstrapKey<Lvl6FastSchedule>;
    using Lvl6HybridKey =
        TFHEpp::CKKSDenseBootstrapHybridGiantKey<Lvl6FastSchedule>;
    using Lvl6CompactHybridKey =
        TFHEpp::CKKSDenseBootstrapHybridGiantKey<Lvl6CompactSchedule>;
    using Lvl6InverseHybridKey =
        TFHEpp::CKKSDenseBootstrapHybridGiantKey<Lvl6InverseSchedule>;
    using Lvl6TunedHybridKey =
        TFHEpp::CKKSDenseBootstrapHybridGiantKey<Lvl6TunedSchedule>;
    using Lvl6SparseFilesystemProvider =
        TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Lvl6FastSchedule>;
    using Lvl6HybridFilesystemProvider =
        TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
            Lvl6FastSchedule>;
    using Lvl6CompactHybridFilesystemProvider =
        TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
            Lvl6CompactSchedule>;
    using Lvl6InverseHybridFilesystemProvider =
        TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
            Lvl6InverseSchedule>;
    using Lvl6TunedHybridFilesystemProvider =
        TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
            Lvl6TunedSchedule>;
    static_assert(std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapInput,
                                 typename Lvl6FastSchedule::InputCiphertext>);
    static_assert(std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapOutput,
                                 typename Lvl6FastSchedule::OutputCiphertext>);
    static_assert(std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapSparseKey,
                                 Lvl6SparseKey>);
    static_assert(std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapHybridGiantKey,
                                 TFHEpp::lvl6CKKSDenseBootstrapFastHybridGiantKey>);
    static_assert(std::is_same_v<
                  TFHEpp::lvl6CKKSDenseBootstrapHybridGiantKey, Lvl6HybridKey>);
    static_assert(std::is_same_v<
                  TFHEpp::lvl6CKKSDenseBootstrapCompactHybridGiantKey,
                  Lvl6CompactHybridKey>);
    static_assert(std::is_same_v<
                  TFHEpp::
                      lvl6CKKSDenseBootstrapSparseFilesystemKeyProvider,
                  Lvl6SparseFilesystemProvider>);
    static_assert(std::is_same_v<
                  TFHEpp::
                      lvl6CKKSDenseBootstrapHybridGiantFilesystemKeyProvider,
                  TFHEpp::
                      lvl6CKKSDenseBootstrapFastHybridGiantFilesystemKeyProvider>);
    static_assert(std::is_same_v<
                  TFHEpp::
                      lvl6CKKSDenseBootstrapFastHybridGiantFilesystemKeyProvider,
                  Lvl6HybridFilesystemProvider>);
    static_assert(std::is_same_v<
                  TFHEpp::
                      lvl6CKKSDenseBootstrapCompactHybridGiantFilesystemKeyProvider,
                  Lvl6CompactHybridFilesystemProvider>);
    static_assert(std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapInverseInput,
                                 typename Lvl6InverseSchedule::InputCiphertext>);
    static_assert(std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapInverseOutput,
                                 typename Lvl6InverseSchedule::OutputCiphertext>);
    static_assert(
        std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapInverseHybridGiantKey,
                       Lvl6InverseHybridKey>);
    static_assert(std::is_same_v<
                  TFHEpp::
                      lvl6CKKSDenseBootstrapInverseHybridGiantFilesystemKeyProvider,
                  Lvl6InverseHybridFilesystemProvider>);
    static_assert(
        std::is_same_v<TFHEpp::lvl6CKKSDenseBootstrapTunedHybridGiantKey,
                       Lvl6TunedHybridKey>);
    static_assert(std::is_same_v<
                  TFHEpp::
                      lvl6CKKSDenseBootstrapTunedHybridGiantFilesystemKeyProvider,
                  Lvl6TunedHybridFilesystemProvider>);
    using Lvl6TunedFreshProductTraits =
        TFHEpp::CKKSDenseBootstrapProductTraits<
            Lvl6TunedSchedule,
            Lvl6TunedSchedule::input_log_q + 2 * Lvl6TunedSchedule::log_delta,
            Lvl6TunedSchedule::log_delta,
            Lvl6TunedSchedule::input_log_q + 2 * Lvl6TunedSchedule::log_delta,
            Lvl6TunedSchedule::log_delta>;
    using Lvl6TunedPostBootstrapProductTraits =
        TFHEpp::CKKSDenseBootstrapProductTraits<
            Lvl6TunedSchedule, Lvl6TunedSchedule::output_log_q,
            Lvl6TunedSchedule::log_delta, Lvl6TunedSchedule::output_log_q,
            Lvl6TunedSchedule::log_delta>;
    static_assert(
        Lvl6TunedFreshProductTraits::ProductCiphertext::log_q ==
        Lvl6TunedSchedule::input_log_q + Lvl6TunedSchedule::log_delta);
    static_assert(Lvl6TunedFreshProductTraits::scale_drop == 0);
    static_assert(Lvl6TunedFreshProductTraits::product_level_compatible);
    static_assert(
        Lvl6TunedPostBootstrapProductTraits::ProductCiphertext::log_q ==
        Lvl6TunedSchedule::post_bootstrap_product_log_q);
    static_assert(Lvl6TunedPostBootstrapProductTraits::scale_drop == 0);
    static_assert(
        Lvl6TunedPostBootstrapProductTraits::product_level_compatible);

    using BootstrapKey = TFHEpp::CKKSDenseBootstrapKey<Schedule>;
    using HybridBootstrapKey =
        TFHEpp::CKKSDenseBootstrapHybridGiantKey<Schedule>;
    using InMemoryProvider =
        TFHEpp::CKKSDenseBootstrapInMemoryKeyProvider<Schedule>;
    using HybridFilesystemProvider =
        TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule>;
    static_assert(std::is_same_v<
                  typename BootstrapKey::CoeffToSlotGaloisKeyChain,
                  TFHEpp::CKKSSparseGaloisKeyChain<
                      M, Schedule::boot_log_q,
                      Schedule::coeff_to_slot_plain_log_delta,
                      Schedule::coeff_to_slot_level_count>>);
    static_assert(std::is_same_v<
                  typename BootstrapKey::SlotToCoeffGaloisKeyChain,
                  TFHEpp::CKKSSparseGaloisKeyChain<
                      M, Schedule::after_evalmod_log_q,
                      Schedule::slot_to_coeff_plain_log_delta,
                      Schedule::slot_to_coeff_level_count>>);
    static_assert(std::is_same_v<
                  TFHEpp::CKKSDenseBootstrapCoeffToSlotGaloisKey<Schedule, 0>,
                  TFHEpp::CKKSSparseGaloisKey<M, Schedule::boot_log_q>>);
    static_assert(std::is_same_v<
                  TFHEpp::CKKSDenseBootstrapCoeffToSlotGaloisKey<Schedule, 1>,
                  TFHEpp::CKKSSparseGaloisKey<
                      M, Schedule::boot_log_q -
                             Schedule::coeff_to_slot_plain_log_delta>>);
    static_assert(std::is_same_v<
                  TFHEpp::CKKSDenseBootstrapPackedConjugateGaloisKey<Schedule>,
                  TFHEpp::CKKSSparseGaloisKey<
                      M, Schedule::after_coeff_to_slot_log_q>>);
    static_assert(std::is_same_v<
                  TFHEpp::CKKSDenseBootstrapSlotToCoeffGaloisKey<Schedule, 1>,
                  TFHEpp::CKKSSparseGaloisKey<
                      M, Schedule::after_evalmod_log_q -
                             Schedule::slot_to_coeff_plain_log_delta>>);
    static_assert(std::is_same_v<
                  typename HybridBootstrapKey::CoeffToSlotGaloisKeyChain,
                  TFHEpp::CKKSHybridSparseGaloisKeyChain<
                      M, Schedule::boot_log_q,
                      Schedule::coeff_to_slot_plain_log_delta,
                      Schedule::coeff_to_slot_level_count>>);
    static_assert(std::is_same_v<
                  typename HybridBootstrapKey::SlotToCoeffGaloisKeyChain,
                  TFHEpp::CKKSHybridSparseGaloisKeyChain<
                      M, Schedule::after_evalmod_log_q,
                      Schedule::slot_to_coeff_plain_log_delta,
                      Schedule::slot_to_coeff_level_count>>);
    static_assert(std::is_same_v<
                  TFHEpp::CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKey<
                      Schedule, 0>,
                  TFHEpp::CKKSHybridSparseGaloisKey<M,
                                                    Schedule::boot_log_q>>);
    static_assert(std::is_same_v<
                  TFHEpp::CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKey<
                      Schedule, 1>,
                  TFHEpp::CKKSHybridSparseGaloisKey<
                      M, Schedule::after_evalmod_log_q -
                             Schedule::slot_to_coeff_plain_log_delta>>);
    using KeyGenFn = void (*)(BootstrapKey &, const TFHEpp::Key<M> &,
                              TFHEpp::CKKSNoise);
    using BootstrapFn = void (*)(typename Schedule::OutputCiphertext &,
                                 const typename Schedule::InputCiphertext &,
                                 const BootstrapKey &);
    using TimedBootstrapFn = void (*)(
        typename Schedule::OutputCiphertext &,
        const typename Schedule::InputCiphertext &, const BootstrapKey &,
        TFHEpp::CKKSDenseBootstrapTimings &);
    using HybridKeyGenFn = void (*)(HybridBootstrapKey &,
                                    const TFHEpp::Key<M> &,
                                    TFHEpp::CKKSNoise);
    using HybridBootstrapFn =
        void (*)(typename Schedule::OutputCiphertext &,
                 const typename Schedule::InputCiphertext &,
                 const HybridBootstrapKey &);
    using TimedHybridBootstrapFn = void (*)(
        typename Schedule::OutputCiphertext &,
        const typename Schedule::InputCiphertext &, const HybridBootstrapKey &,
        TFHEpp::CKKSDenseBootstrapTimings &);
    using ProviderBootstrapFn =
        void (*)(typename Schedule::OutputCiphertext &,
                 const typename Schedule::InputCiphertext &,
                 const InMemoryProvider &);
    using HybridFilesystemBootstrapFn =
        void (*)(typename Schedule::OutputCiphertext &,
                 const typename Schedule::InputCiphertext &,
                 const HybridFilesystemProvider &);
    [[maybe_unused]] KeyGenFn keygen =
        &TFHEpp::CKKSDenseBootstrapKeyGen<Schedule>;
    [[maybe_unused]] BootstrapFn bootstrap =
        &TFHEpp::CKKSDenseBootstrap<Schedule>;
    [[maybe_unused]] TimedBootstrapFn timed_bootstrap =
        &TFHEpp::CKKSDenseBootstrapTimed<Schedule>;
    [[maybe_unused]] HybridKeyGenFn hybrid_keygen =
        &TFHEpp::CKKSDenseBootstrapHybridGiantKeyGen<Schedule>;
    [[maybe_unused]] HybridBootstrapFn hybrid_bootstrap =
        &TFHEpp::CKKSDenseBootstrapHybridGiant<Schedule>;
    [[maybe_unused]] TimedHybridBootstrapFn timed_hybrid_bootstrap =
        &TFHEpp::CKKSDenseBootstrapHybridGiantTimed<Schedule>;
    [[maybe_unused]] ProviderBootstrapFn provider_bootstrap =
        &TFHEpp::CKKSDenseBootstrapWithKeyProvider<Schedule,
                                                   InMemoryProvider>;
    [[maybe_unused]] HybridFilesystemBootstrapFn
        hybrid_filesystem_bootstrap =
            &TFHEpp::CKKSDenseBootstrapWithKeyProvider<
                Schedule, HybridFilesystemProvider>;

    TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    TFHEpp::CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    if (linear_plan.coeff_to_slot_stages.size() !=
        Schedule::coeff_to_slot_level_count)
        std::exit(1);
    if (linear_plan.slot_to_coeff_stages.size() !=
        Schedule::slot_to_coeff_level_count)
        std::exit(1);
    if (linear_plan.slot_to_coeff_imag_stages.size() !=
        Schedule::slot_to_coeff_level_count)
        std::exit(1);

    TFHEpp::CKKSDenseBootstrapRotationKeyUsage<Schedule> rotation_usage;
    TFHEpp::CKKSBuildDenseBootstrapRotationKeyUsage<Schedule>(rotation_usage,
                                                              linear_plan);
    const std::size_t planned_key_indices =
        TFHEpp::CKKSDenseBootstrapRotationKeyUsageCount<Schedule>(
            rotation_usage);
    constexpr std::size_t full_key_indices =
        TFHEpp::CKKSDenseBootstrapFullGaloisKeyIndexCount<Schedule>();
    if (planned_key_indices == 0 || planned_key_indices >= full_key_indices)
        std::exit(1);
    TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        hybrid_usage;
    TFHEpp::CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        hybrid_usage, linear_plan);
    const std::size_t hybrid_key_indices =
        TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsageCount<Schedule>(
            hybrid_usage);
    if (hybrid_key_indices == 0 || hybrid_key_indices >= full_key_indices)
        std::exit(1);
}

void test_dense_bootstrap_plain_split_pipeline()
{
    using M = TinyDeepMultiLimbCKKSParam;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        M, 40, 8, 560, 40, 3, 31, 4, 2, 0, 40, 2>;
    const auto poly = TFHEpp::CKKSBuildBoundedCosEvalModPolynomial<Schedule>();

    std::array<double, M::n> coeffs{};
    std::array<double, M::n> masked_coeffs{};
    for (std::size_t i = 0; i < M::n; i++) {
        coeffs[i] =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 64.0;
        const int mask =
            static_cast<int>(i % (2 * Schedule::evalmod_k - 1)) -
            static_cast<int>(Schedule::evalmod_k) + 1;
        masked_coeffs[i] =
            static_cast<double>(mask) * Schedule::message_ratio + coeffs[i];
    }

    auto masked_slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto wanted_slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    coeffs_to_slots<M>(*masked_slots, masked_coeffs);
    coeffs_to_slots<M>(*wanted_slots, coeffs);

    TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    TFHEpp::CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);

    auto packed = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    apply_complex_stages<M>(*packed, *masked_slots,
                            linear_plan.coeff_to_slot_stages);

    auto real_eval = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto imag_eval = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto combined_eval = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    for (std::size_t i = 0; i < M::n / 2; i++) {
        const double real_value =
            TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                poly, packed->at(i).real()) /
            Schedule::message_ratio;
        const double imag_value =
            TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                poly, packed->at(i).imag()) /
            Schedule::message_ratio;
        (*real_eval)[i] = {real_value, 0.0};
        (*imag_eval)[i] = {imag_value, 0.0};
        (*combined_eval)[i] = {real_value, imag_value};
    }

    auto direct_stc = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    apply_complex_stages<M>(*direct_stc, *combined_eval,
                            linear_plan.slot_to_coeff_stages);

    auto real_stc = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto imag_stc = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto split_stc = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    apply_complex_stages<M>(*real_stc, *real_eval,
                            linear_plan.slot_to_coeff_stages);
    apply_complex_stages<M>(*imag_stc, *imag_eval,
                            linear_plan.slot_to_coeff_imag_stages);
    for (std::size_t i = 0; i < M::n / 2; i++)
        (*split_stc)[i] = (*real_stc)[i] + (*imag_stc)[i];

    require_close_param<M>(*split_stc, *direct_stc, 1e-10,
                           "CKKS dense split EvalMod STC mirror");
    require_close_param<M>(*split_stc, *wanted_slots, 0.02,
                           "CKKS dense plaintext bootstrap mirror");
}

void test_dense_bootstrap_encrypted_pipeline()
{
    using M = TinyDeepMultiLimbCKKSParam;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        M, 40, 8, 560, 40, 3, 31, 4, 2, 0, 40, 2>;

    std::array<double, M::n> coeffs{};
    std::array<double, M::n> masked_coeffs{};
    for (std::size_t i = 0; i < M::n; i++) {
        coeffs[i] =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 64.0;
        const int mask =
            static_cast<int>(i % (2 * Schedule::evalmod_k - 1)) -
            static_cast<int>(Schedule::evalmod_k) + 1;
        masked_coeffs[i] =
            static_cast<double>(mask) * Schedule::message_ratio + coeffs[i];
    }

    auto input_slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto wanted_slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    coeffs_to_slots<M>(*input_slots, masked_coeffs);
    coeffs_to_slots<M>(*wanted_slots, coeffs);

    TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    TFHEpp::CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);

    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    apply_complex_stages<M>(*expected, *input_slots,
                            linear_plan.coeff_to_slot_stages);

    auto key = std::make_unique<TFHEpp::Key<M>>();
    fill_test_key<M>(*key);

    auto ct = std::make_unique<typename Schedule::BootstrapCiphertext>();
    TFHEpp::ckksSlotEncrypt<M, Schedule::boot_log_q, Schedule::log_delta>(
        *ct, *input_slots, *key, {0.0, 0});

    auto gks = std::make_unique<typename TFHEpp::CKKSGaloisKeyChain<
        M, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>>();
    TFHEpp::CKKSGaloisKeyChainGen<M, Schedule::boot_log_q,
                                  Schedule::coeff_to_slot_plain_log_delta,
                                  Schedule::coeff_to_slot_level_count>(
        *gks, *key, {0.0, 0});

    auto out = std::make_unique<typename Schedule::CoeffToSlotCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGS<
        M, Schedule::boot_log_q, Schedule::log_delta,
        Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>(
        *out, *ct, linear_plan.coeff_to_slot_stages, 0,
        Schedule::linear_bsgs_step, *gks);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    TFHEpp::ckksSlotDecrypt<M, Schedule::after_coeff_to_slot_log_q,
                            Schedule::log_delta>(*decoded, *out, *key);
    require_close_param<M>(*decoded, *expected, 1e-6,
                           "CKKS dense encrypted C2S");

    std::array<TFHEpp::CKKSDirectRotationKeyIndexSet<M>,
               Schedule::coeff_to_slot_level_count + 1>
        direct_c2s_usage{};
    TFHEpp::CKKSClearDirectRotationKeyIndexSets<M>(direct_c2s_usage);
    for (std::size_t i = 0; i < linear_plan.coeff_to_slot_stages.size(); i++)
        TFHEpp::CKKSCollectLinearTransformStageDirectRotationKeyIndices<M>(
            direct_c2s_usage[i], direct_c2s_usage[i + 1],
            linear_plan.coeff_to_slot_stages[i],
            Schedule::linear_bsgs_step);
    auto direct_c2s_gks =
        std::make_unique<TFHEpp::CKKSDirectSparseGaloisKeyChain<
            M, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
            Schedule::coeff_to_slot_level_count>>();
    TFHEpp::CKKSDirectSparseGaloisKeyChainGen<
        M, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>(*direct_c2s_gks,
                                             direct_c2s_usage, *key,
                                             {0.0, 0});
    auto direct_out =
        std::make_unique<typename Schedule::CoeffToSlotCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGSDirect<
        M, Schedule::boot_log_q, Schedule::log_delta,
        Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>(
        *direct_out, *ct, linear_plan.coeff_to_slot_stages, 0,
        Schedule::linear_bsgs_step, *direct_c2s_gks);
    TFHEpp::ckksSlotDecrypt<M, Schedule::after_coeff_to_slot_log_q,
                            Schedule::log_delta>(*decoded, *direct_out, *key);
    require_close_param<M>(*decoded, *expected, 1e-6,
                           "CKKS dense encrypted direct C2S");

    std::array<TFHEpp::CKKSRotationKeyIndexSet<M>,
               Schedule::coeff_to_slot_level_count + 1>
        hybrid_c2s_binary_usage{};
    std::array<TFHEpp::CKKSDirectRotationKeyIndexSet<M>,
               Schedule::coeff_to_slot_level_count + 1>
        hybrid_c2s_direct_usage{};
    TFHEpp::CKKSClearRotationKeyIndexSets<M>(hybrid_c2s_binary_usage);
    TFHEpp::CKKSClearDirectRotationKeyIndexSets<M>(
        hybrid_c2s_direct_usage);
    for (std::size_t i = 0; i < linear_plan.coeff_to_slot_stages.size(); i++)
        TFHEpp::CKKSCollectLinearTransformStageHybridGiantRotationKeyIndices<
            M>(hybrid_c2s_binary_usage[i],
               hybrid_c2s_binary_usage[i + 1],
               hybrid_c2s_direct_usage[i + 1],
               linear_plan.coeff_to_slot_stages[i],
               Schedule::linear_bsgs_step);
    auto hybrid_c2s_gks =
        std::make_unique<TFHEpp::CKKSHybridSparseGaloisKeyChain<
            M, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
            Schedule::coeff_to_slot_level_count>>();
    TFHEpp::CKKSHybridSparseGaloisKeyChainGen<
        M, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>(
        *hybrid_c2s_gks, hybrid_c2s_binary_usage, hybrid_c2s_direct_usage,
        *key, {0.0, 0});
    auto hybrid_out =
        std::make_unique<typename Schedule::CoeffToSlotCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGS<
        M, Schedule::boot_log_q, Schedule::log_delta,
        Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>(
        *hybrid_out, *ct, linear_plan.coeff_to_slot_stages, 0,
        Schedule::linear_bsgs_step, *hybrid_c2s_gks);
    TFHEpp::ckksSlotDecrypt<M, Schedule::after_coeff_to_slot_log_q,
                            Schedule::log_delta>(*decoded, *hybrid_out, *key);
    require_close_param<M>(*decoded, *expected, 1e-6,
                           "CKKS dense encrypted hybrid C2S");

    auto conjugate_gk =
        std::make_unique<TFHEpp::CKKSGaloisKey<
            M, Schedule::after_coeff_to_slot_log_q>>();
    TFHEpp::CKKSGaloisKeyGen<M, Schedule::after_coeff_to_slot_log_q>(
        *conjugate_gk, *key, {0.0, 0});

    auto real_component =
        std::make_unique<typename Schedule::ComponentCiphertext>();
    TFHEpp::CKKSExtractRealSlots<M, Schedule::after_coeff_to_slot_log_q,
                                 Schedule::log_delta,
                                 Schedule::component_split_plain_log_delta>(
        *real_component, *out, *conjugate_gk);
    auto imag_component =
        std::make_unique<typename Schedule::ComponentCiphertext>();
    TFHEpp::CKKSExtractImagSlots<M, Schedule::after_coeff_to_slot_log_q,
                                 Schedule::log_delta,
                                 Schedule::component_split_plain_log_delta>(
        *imag_component, *out, *conjugate_gk);

    auto expected_real = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto expected_imag = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    for (std::size_t i = 0; i < M::n / 2; i++) {
        (*expected_real)[i] = {expected->at(i).real(), 0.0};
        (*expected_imag)[i] = {expected->at(i).imag(), 0.0};
    }

    TFHEpp::ckksSlotDecrypt<M, Schedule::after_component_split_log_q,
                            Schedule::log_delta>(*decoded, *real_component,
                                                 *key);
    require_close_param<M>(*decoded, *expected_real, 1e-6,
                           "CKKS dense encrypted C2S real split");
    TFHEpp::ckksSlotDecrypt<M, Schedule::after_component_split_log_q,
                            Schedule::log_delta>(*decoded, *imag_component,
                                                 *key);
    require_close_param<M>(*decoded, *expected_imag, 1e-6,
                           "CKKS dense encrypted C2S imag split");

    const auto poly = TFHEpp::CKKSBuildBoundedCosEvalModPolynomial<Schedule>();
    auto evalmod_keys =
        std::make_unique<TFHEpp::CKKSDenseEvalModBoundedCosRelinKeys<
            Schedule>>();
    TFHEpp::CKKSDenseEvalModBoundedCosKeyGen<Schedule>(*evalmod_keys, *key,
                                                       {0.0, 0});

    auto real_eval =
        std::make_unique<TFHEpp::CKKSDenseEvalModBoundedCosResult<Schedule>>();
    TFHEpp::CKKSDenseEvalModBoundedCosNormalized<Schedule>(
        *real_eval, *real_component, poly, *evalmod_keys);
    auto imag_eval =
        std::make_unique<TFHEpp::CKKSDenseEvalModBoundedCosResult<Schedule>>();
    TFHEpp::CKKSDenseEvalModBoundedCosNormalized<Schedule>(
        *imag_eval, *imag_component, poly, *evalmod_keys);

    auto expected_real_eval = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto expected_imag_eval = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    for (std::size_t i = 0; i < M::n / 2; i++) {
        (*expected_real_eval)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, expected_real->at(i).real()) /
                 Schedule::message_ratio,
             0.0};
        (*expected_imag_eval)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosNormalizedPower(
                 poly, expected_imag->at(i).real()) /
                 Schedule::message_ratio,
             0.0};
    }

    TFHEpp::ckksSlotDecrypt<M, Schedule::after_evalmod_log_q,
                            Schedule::log_delta>(*decoded, *real_eval, *key);
    require_close_param<M>(*decoded, *expected_real_eval, 1e-5,
                           "CKKS dense encrypted real EvalMod");
    TFHEpp::ckksSlotDecrypt<M, Schedule::after_evalmod_log_q,
                            Schedule::log_delta>(*decoded, *imag_eval, *key);
    require_close_param<M>(*decoded, *expected_imag_eval, 1e-5,
                           "CKKS dense encrypted imag EvalMod");

    auto slot_to_coeff_gks =
        std::make_unique<typename TFHEpp::CKKSGaloisKeyChain<
            M, Schedule::after_evalmod_log_q,
            Schedule::slot_to_coeff_plain_log_delta,
            Schedule::slot_to_coeff_level_count>>();
    TFHEpp::CKKSGaloisKeyChainGen<M, Schedule::after_evalmod_log_q,
                                  Schedule::slot_to_coeff_plain_log_delta,
                                  Schedule::slot_to_coeff_level_count>(
        *slot_to_coeff_gks, *key, {0.0, 0});

    auto real_out = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGS<
        M, Schedule::after_evalmod_log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *real_out, *real_eval, linear_plan.slot_to_coeff_stages, 0,
        Schedule::linear_bsgs_step, *slot_to_coeff_gks);
    auto imag_out = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGS<
        M, Schedule::after_evalmod_log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *imag_out, *imag_eval, linear_plan.slot_to_coeff_imag_stages, 0,
        Schedule::linear_bsgs_step, *slot_to_coeff_gks);

    auto combined_out =
        std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSAdd<M, Schedule::output_log_q, Schedule::log_delta>(
        *combined_out, *real_out, *imag_out);

    auto fused_out = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGSDualInputSharedTail<
        M, Schedule::after_evalmod_log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *fused_out, *real_eval, *imag_eval, linear_plan.slot_to_coeff_stages,
        linear_plan.slot_to_coeff_imag_stages, 0,
        Schedule::linear_bsgs_step, *slot_to_coeff_gks);

    auto expected_real_out = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto expected_imag_out = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto expected_combined_out =
        std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    apply_complex_stages<M>(*expected_real_out, *expected_real_eval,
                            linear_plan.slot_to_coeff_stages);
    apply_complex_stages<M>(*expected_imag_out, *expected_imag_eval,
                            linear_plan.slot_to_coeff_imag_stages);
    for (std::size_t i = 0; i < M::n / 2; i++)
        (*expected_combined_out)[i] =
            (*expected_real_out)[i] + (*expected_imag_out)[i];

    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *combined_out, *key);
    require_close_param<M>(*decoded, *expected_combined_out, 1e-3,
                           "CKKS dense encrypted STC recombination");

    auto fused_decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *fused_decoded, *fused_out, *key);
    require_close_param<M>(*fused_decoded, *expected_combined_out, 1e-3,
                           "CKKS dense encrypted fused STC recombination");
    require_close_param<M>(*fused_decoded, *decoded, 1e-3,
                           "CKKS dense encrypted fused STC matches split STC");
    require_close_param<M>(*fused_decoded, *wanted_slots, 0.02,
                           "CKKS dense encrypted bootstrap without modraise");
}

void test_dense_bootstrap_fused_stc_shared_tail()
{
    using M = TinyDeepMultiLimbCKKSParam;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        M, 30, 8, 540, 30, 1, 31, 4, 2, 0, 30, 2>;
    static_assert(Schedule::slot_to_coeff_level_count == 3);
    static_assert(Schedule::after_evalmod_log_q == 180);
    static_assert(Schedule::output_log_q == 90);

    TFHEpp::CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    TFHEpp::CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);

    auto real_slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto imag_slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    for (std::size_t i = 0; i < M::n / 2; i++) {
        (*real_slots)[i] =
            {static_cast<double>(static_cast<int>(i) - 3) / 128.0, 0.0};
        (*imag_slots)[i] =
            {static_cast<double>(4 - static_cast<int>(i)) / 160.0, 0.0};
    }

    auto key = std::make_unique<TFHEpp::Key<M>>();
    fill_test_key<M>(*key);
    auto real_ct = std::make_unique<typename Schedule::EvalModCiphertext>();
    auto imag_ct = std::make_unique<typename Schedule::EvalModCiphertext>();
    TFHEpp::ckksSlotEncrypt<M, Schedule::after_evalmod_log_q,
                            Schedule::log_delta>(*real_ct, *real_slots, *key,
                                                 {0.0, 0});
    TFHEpp::ckksSlotEncrypt<M, Schedule::after_evalmod_log_q,
                            Schedule::log_delta>(*imag_ct, *imag_slots, *key,
                                                 {0.0, 0});

    auto slot_to_coeff_gks =
        std::make_unique<typename TFHEpp::CKKSGaloisKeyChain<
            M, Schedule::after_evalmod_log_q,
            Schedule::slot_to_coeff_plain_log_delta,
            Schedule::slot_to_coeff_level_count>>();
    TFHEpp::CKKSGaloisKeyChainGen<M, Schedule::after_evalmod_log_q,
                                  Schedule::slot_to_coeff_plain_log_delta,
                                  Schedule::slot_to_coeff_level_count>(
        *slot_to_coeff_gks, *key, {0.0, 0});

    auto real_out = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGS<
        M, Schedule::after_evalmod_log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *real_out, *real_ct, linear_plan.slot_to_coeff_stages, 0,
        Schedule::linear_bsgs_step, *slot_to_coeff_gks);
    auto imag_out = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGS<
        M, Schedule::after_evalmod_log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *imag_out, *imag_ct, linear_plan.slot_to_coeff_imag_stages, 0,
        Schedule::linear_bsgs_step, *slot_to_coeff_gks);
    auto split_out = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSAdd<M, Schedule::output_log_q, Schedule::log_delta>(
        *split_out, *real_out, *imag_out);

    auto fused_out = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGSDualInputSharedTail<
        M, Schedule::after_evalmod_log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *fused_out, *real_ct, *imag_ct, linear_plan.slot_to_coeff_stages,
        linear_plan.slot_to_coeff_imag_stages, 0,
        Schedule::linear_bsgs_step, *slot_to_coeff_gks);

    std::array<TFHEpp::CKKSDirectRotationKeyIndexSet<M>,
               Schedule::slot_to_coeff_level_count + 1>
        direct_stc_usage{};
    TFHEpp::CKKSClearDirectRotationKeyIndexSets<M>(direct_stc_usage);
    for (std::size_t i = 0; i < linear_plan.slot_to_coeff_stages.size(); i++) {
        TFHEpp::CKKSCollectLinearTransformStageDirectRotationKeyIndices<M>(
            direct_stc_usage[i], direct_stc_usage[i + 1],
            linear_plan.slot_to_coeff_stages[i],
            Schedule::linear_bsgs_step);
        TFHEpp::CKKSCollectLinearTransformStageDirectRotationKeyIndices<M>(
            direct_stc_usage[i], direct_stc_usage[i + 1],
            linear_plan.slot_to_coeff_imag_stages[i],
            Schedule::linear_bsgs_step);
    }
    auto direct_slot_to_coeff_gks =
        std::make_unique<TFHEpp::CKKSDirectSparseGaloisKeyChain<
            M, Schedule::after_evalmod_log_q,
            Schedule::slot_to_coeff_plain_log_delta,
            Schedule::slot_to_coeff_level_count>>();
    TFHEpp::CKKSDirectSparseGaloisKeyChainGen<
        M, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(*direct_slot_to_coeff_gks,
                                             direct_stc_usage, *key,
                                             {0.0, 0});
    auto direct_fused_out =
        std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGSDualInputSharedTailDirect<
        M, Schedule::after_evalmod_log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *direct_fused_out, *real_ct, *imag_ct,
        linear_plan.slot_to_coeff_stages, linear_plan.slot_to_coeff_imag_stages,
        0, Schedule::linear_bsgs_step, *direct_slot_to_coeff_gks);

    std::array<TFHEpp::CKKSRotationKeyIndexSet<M>,
               Schedule::slot_to_coeff_level_count + 1>
        hybrid_stc_binary_usage{};
    std::array<TFHEpp::CKKSDirectRotationKeyIndexSet<M>,
               Schedule::slot_to_coeff_level_count + 1>
        hybrid_stc_direct_usage{};
    TFHEpp::CKKSClearRotationKeyIndexSets<M>(hybrid_stc_binary_usage);
    TFHEpp::CKKSClearDirectRotationKeyIndexSets<M>(
        hybrid_stc_direct_usage);
    for (std::size_t i = 0; i < linear_plan.slot_to_coeff_stages.size(); i++) {
        TFHEpp::CKKSCollectLinearTransformStageHybridGiantRotationKeyIndices<
            M>(hybrid_stc_binary_usage[i],
               hybrid_stc_binary_usage[i + 1],
               hybrid_stc_direct_usage[i + 1],
               linear_plan.slot_to_coeff_stages[i],
               Schedule::linear_bsgs_step);
        TFHEpp::CKKSCollectLinearTransformStageHybridGiantRotationKeyIndices<
            M>(hybrid_stc_binary_usage[i],
               hybrid_stc_binary_usage[i + 1],
               hybrid_stc_direct_usage[i + 1],
               linear_plan.slot_to_coeff_imag_stages[i],
               Schedule::linear_bsgs_step);
    }
    auto hybrid_slot_to_coeff_gks =
        std::make_unique<TFHEpp::CKKSHybridSparseGaloisKeyChain<
            M, Schedule::after_evalmod_log_q,
            Schedule::slot_to_coeff_plain_log_delta,
            Schedule::slot_to_coeff_level_count>>();
    TFHEpp::CKKSHybridSparseGaloisKeyChainGen<
        M, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *hybrid_slot_to_coeff_gks, hybrid_stc_binary_usage,
        hybrid_stc_direct_usage, *key, {0.0, 0});
    auto hybrid_fused_out =
        std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSLinearTransformStagesBSGSDualInputSharedTail<
        M, Schedule::after_evalmod_log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        *hybrid_fused_out, *real_ct, *imag_ct,
        linear_plan.slot_to_coeff_stages, linear_plan.slot_to_coeff_imag_stages,
        0, Schedule::linear_bsgs_step, *hybrid_slot_to_coeff_gks);

    auto expected_real = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto expected_imag = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    apply_complex_stages<M>(*expected_real, *real_slots,
                            linear_plan.slot_to_coeff_stages);
    apply_complex_stages<M>(*expected_imag, *imag_slots,
                            linear_plan.slot_to_coeff_imag_stages);
    for (std::size_t i = 0; i < M::n / 2; i++)
        (*expected)[i] = (*expected_real)[i] + (*expected_imag)[i];

    auto split_decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto fused_decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto direct_fused_decoded =
        std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto hybrid_fused_decoded =
        std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *split_decoded, *split_out, *key);
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *fused_decoded, *fused_out, *key);
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *direct_fused_decoded, *direct_fused_out, *key);
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *hybrid_fused_decoded, *hybrid_fused_out, *key);
    require_close_param<M>(*split_decoded, *expected, 1e-3,
                           "CKKS dense encrypted split STC shared tail");
    require_close_param<M>(*fused_decoded, *expected, 1e-3,
                           "CKKS dense encrypted fused STC shared tail");
    require_close_param<M>(*direct_fused_decoded, *expected, 1e-3,
                           "CKKS dense encrypted direct fused STC shared tail");
    require_close_param<M>(*hybrid_fused_decoded, *expected, 1e-3,
                           "CKKS dense encrypted hybrid fused STC shared tail");
    require_close_param<M>(*fused_decoded, *split_decoded, 1e-3,
                           "CKKS dense encrypted fused/split shared tail");
    require_close_param<M>(*direct_fused_decoded, *fused_decoded, 1e-3,
                           "CKKS dense encrypted direct/current shared tail");
    require_close_param<M>(*hybrid_fused_decoded, *fused_decoded, 1e-3,
                           "CKKS dense encrypted hybrid/current shared tail");
}

void test_dense_bootstrap_e2e_smoke()
{
    using M = TinyDeepMultiLimbCKKSParam;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        M, 40, 8, 560, 40, 3, 31, 4, 2, 0, 40, 2>;

    std::array<double, M::n> coeffs{};
    for (std::size_t i = 0; i < M::n; i++)
        coeffs[i] =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 64.0;

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    coeffs_to_slots<M>(*slots, coeffs);

    auto key = std::make_unique<TFHEpp::Key<M>>();
    fill_test_key<M>(*key);

    auto bootstrap_key =
        std::make_unique<TFHEpp::CKKSDenseBootstrapKey<Schedule>>();
    TFHEpp::CKKSDenseBootstrapKeyGen<Schedule>(*bootstrap_key, *key,
                                               {0.0, 0});
    TFHEpp::CKKSDenseBootstrapRotationKeyUsage<Schedule> expected_usage;
    TFHEpp::CKKSBuildDenseBootstrapRotationKeyUsage<Schedule>(
        expected_usage, bootstrap_key->linear_plan);
    if (bootstrap_key->coeff_to_slot_galois.template get<0>().available !=
        expected_usage.coeff_to_slot[0])
        std::exit(1);
    if (bootstrap_key->coeff_to_slot_galois.template get<1>().available !=
        expected_usage.coeff_to_slot[1])
        std::exit(1);
    if (bootstrap_key->slot_to_coeff_galois.template get<0>().available !=
        expected_usage.slot_to_coeff[0])
        std::exit(1);
    if (bootstrap_key->slot_to_coeff_galois.template get<1>().available !=
        expected_usage.slot_to_coeff[1])
        std::exit(1);
    if (!bootstrap_key->packed_conjugate_galois.has(M::nbit) ||
        TFHEpp::CKKSRotationKeyIndexSetCount<M>(
            bootstrap_key->packed_conjugate_galois.available) != 1)
        std::exit(1);

    auto hybrid_bootstrap_key = std::make_unique<
        TFHEpp::CKKSDenseBootstrapHybridGiantKey<Schedule>>();
    TFHEpp::CKKSDenseBootstrapHybridGiantKeyGen<Schedule>(
        *hybrid_bootstrap_key, *key, {0.0, 0});
    TFHEpp::CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        expected_hybrid_usage;
    TFHEpp::CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        expected_hybrid_usage, hybrid_bootstrap_key->linear_plan);
    if (hybrid_bootstrap_key->coeff_to_slot_galois.template get<0>()
            .binary.available != expected_hybrid_usage.coeff_to_slot_binary[0])
        std::exit(1);
    if (hybrid_bootstrap_key->coeff_to_slot_galois.template get<1>()
            .binary.available != expected_hybrid_usage.coeff_to_slot_binary[1])
        std::exit(1);
    if (hybrid_bootstrap_key->coeff_to_slot_galois.template get<0>()
            .direct.keys.size() !=
        TFHEpp::CKKSDirectRotationKeyIndexSetCount<M>(
            expected_hybrid_usage.coeff_to_slot_direct[0]))
        std::exit(1);
    if (hybrid_bootstrap_key->coeff_to_slot_galois.template get<1>()
            .direct.keys.size() !=
        TFHEpp::CKKSDirectRotationKeyIndexSetCount<M>(
            expected_hybrid_usage.coeff_to_slot_direct[1]))
        std::exit(1);
    if (hybrid_bootstrap_key->slot_to_coeff_galois.template get<0>()
            .binary.available != expected_hybrid_usage.slot_to_coeff_binary[0])
        std::exit(1);
    if (hybrid_bootstrap_key->slot_to_coeff_galois.template get<1>()
            .binary.available != expected_hybrid_usage.slot_to_coeff_binary[1])
        std::exit(1);
    if (hybrid_bootstrap_key->slot_to_coeff_galois.template get<0>()
            .direct.keys.size() !=
        TFHEpp::CKKSDirectRotationKeyIndexSetCount<M>(
            expected_hybrid_usage.slot_to_coeff_direct[0]))
        std::exit(1);
    if (hybrid_bootstrap_key->slot_to_coeff_galois.template get<1>()
            .direct.keys.size() !=
        TFHEpp::CKKSDirectRotationKeyIndexSetCount<M>(
            expected_hybrid_usage.slot_to_coeff_direct[1]))
        std::exit(1);
    if (hybrid_bootstrap_key->packed_conjugate_galois.available !=
            expected_hybrid_usage.packed_conjugate ||
        !hybrid_bootstrap_key->packed_conjugate_galois.has(M::nbit))
        std::exit(1);

    TFHEpp::CKKSDenseBootstrapCoeffToSlotGaloisKey<Schedule, 0>
        streamed_c2s0;
    TFHEpp::CKKSDenseBootstrapCoeffToSlotGaloisKeyGen<Schedule, 0>(
        streamed_c2s0, expected_usage, *key, {0.0, 0});
    if (streamed_c2s0.available != expected_usage.coeff_to_slot[0])
        std::exit(1);
    TFHEpp::CKKSDenseBootstrapPackedConjugateGaloisKey<Schedule>
        streamed_packed;
    TFHEpp::CKKSDenseBootstrapPackedConjugateGaloisKeyGen<Schedule>(
        streamed_packed, expected_usage, *key, {0.0, 0});
    if (streamed_packed.available != expected_usage.packed_conjugate ||
        !streamed_packed.has(M::nbit))
        std::exit(1);
    TFHEpp::CKKSDenseBootstrapSlotToCoeffGaloisKey<Schedule, 1>
        streamed_stc1;
    TFHEpp::CKKSDenseBootstrapSlotToCoeffGaloisKeyGen<Schedule, 1>(
        streamed_stc1, expected_usage, *key, {0.0, 0});
    if (streamed_stc1.available != expected_usage.slot_to_coeff[1])
        std::exit(1);

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    TFHEpp::ckksSlotEncrypt<M, Schedule::input_log_q, Schedule::log_delta>(
        *input, *slots, *key, {0.0, 0});

    CountingDenseBootstrapProvider<Schedule> counting_provider(*bootstrap_key);
    auto provider_output =
        std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapWithKeyProvider<Schedule>(
        *provider_output, *input, counting_provider);
    if (counting_provider.coeff_to_slot_gets[0] != 1 ||
        counting_provider.coeff_to_slot_gets[1] != 1 ||
        counting_provider.slot_to_coeff_gets[0] != 1 ||
        counting_provider.slot_to_coeff_gets[1] != 1 ||
        counting_provider.packed_conjugate_gets != 2 ||
        counting_provider.evalmod_polynomial_gets != 2)
        std::exit(1);
    if (counting_provider.polynomial_relin_gets[0] != 2 ||
        counting_provider.polynomial_relin_gets[1] != 4 ||
        counting_provider.polynomial_relin_gets[2] != 8 ||
        counting_provider.polynomial_relin_gets[3] != 16 ||
        counting_provider.polynomial_relin_gets[4] != 30 ||
        counting_provider.double_angle_relin_gets[0] != 2 ||
        counting_provider.double_angle_relin_gets[1] != 2)
        std::exit(1);
    if (counting_provider.coeff_to_slot_releases[0] != 1 ||
        counting_provider.coeff_to_slot_releases[1] != 1 ||
        counting_provider.slot_to_coeff_releases[0] != 1 ||
        counting_provider.slot_to_coeff_releases[1] != 1 ||
        counting_provider.packed_conjugate_releases != 1)
        std::exit(1);
    if (counting_provider.polynomial_relin_releases[0] != 2 ||
        counting_provider.polynomial_relin_releases[1] != 2 ||
        counting_provider.polynomial_relin_releases[2] != 2 ||
        counting_provider.polynomial_relin_releases[3] != 2 ||
        counting_provider.polynomial_relin_releases[4] != 2 ||
        counting_provider.double_angle_relin_releases[0] != 2 ||
        counting_provider.double_angle_relin_releases[1] != 2)
        std::exit(1);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *provider_output, *key);
    require_close_param<M>(*decoded, *slots, 0.02,
                           "CKKS dense encrypted provider bootstrap e2e");

    auto hybrid_output =
        std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapTimings hybrid_timings;
    TFHEpp::CKKSDenseBootstrapHybridGiantTimed<Schedule>(
        *hybrid_output, *input, *hybrid_bootstrap_key, hybrid_timings);
    if (hybrid_timings.total_ms() <= 0.0 ||
        hybrid_timings.coeff_to_slot_ms <= 0.0 ||
        hybrid_timings.real_evalmod_ms <= 0.0 ||
        hybrid_timings.slot_to_coeff_ms <= 0.0)
        std::exit(1);
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *hybrid_output, *key);
    require_close_param<M>(*decoded, *slots, 0.02,
                           "CKKS dense encrypted hybrid bootstrap e2e");

    const std::filesystem::path hybrid_slice_dir =
        std::filesystem::temp_directory_path() /
        "tfhepp_ckks_dense_bootstrap_hybrid_slices";
    std::filesystem::remove_all(hybrid_slice_dir);
    std::size_t generated_hybrid_slices = 0;
    while (
        TFHEpp::CKKSDenseBootstrapHybridGiantKeyGenNextMissingToDirectory<
            Schedule>(hybrid_slice_dir, *key, {0.0, 0})) {
        generated_hybrid_slices++;
    }
    if (generated_hybrid_slices == 0 ||
        !TFHEpp::CKKSDenseBootstrapHybridGiantKeyDirectoryComplete<Schedule>(
            hybrid_slice_dir) ||
        !TFHEpp::CKKSDenseBootstrapHybridGiantKeyDirectoryManifestMatches<
            Schedule>(hybrid_slice_dir) ||
        TFHEpp::CKKSDenseBootstrapKeyDirectoryManifestMatches<Schedule>(
            hybrid_slice_dir))
        std::exit(1);
    TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions hybrid_resume_options;
    hybrid_resume_options.overwrite_existing = false;
    TFHEpp::CKKSDenseBootstrapHybridGiantKeyGenToDirectory<Schedule>(
        hybrid_slice_dir, *key, {0.0, 0}, hybrid_resume_options);
    TFHEpp::CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule>
        hybrid_filesystem_provider(hybrid_slice_dir);
    auto hybrid_filesystem_output =
        std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapWithKeyProvider<Schedule>(
        *hybrid_filesystem_output, *input, hybrid_filesystem_provider);
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *hybrid_filesystem_output, *key);
    require_close_param<M>(*decoded, *slots, 0.02,
                           "CKKS dense hybrid filesystem bootstrap e2e");
    std::filesystem::remove_all(hybrid_slice_dir);

    const std::filesystem::path slice_dir =
        std::filesystem::temp_directory_path() /
        "tfhepp_ckks_dense_bootstrap_slices";
    std::filesystem::remove_all(slice_dir);
    TFHEpp::CKKSDenseBootstrapKeyGenToDirectory<Schedule>(slice_dir, *key,
                                                          {0.0, 0});
    if (!std::filesystem::exists(slice_dir / "linear_plan.bin") ||
        !std::filesystem::exists(slice_dir / "coeff_to_slot_galois_0.bin") ||
        !std::filesystem::exists(slice_dir / "packed_conjugate_galois.bin") ||
        !std::filesystem::exists(slice_dir / "polynomial_relin_0.bin") ||
        !std::filesystem::exists(slice_dir / "double_angle_relin_0.bin") ||
        !std::filesystem::exists(slice_dir / "slot_to_coeff_galois_1.bin"))
        std::exit(1);
    TFHEpp::CKKSDenseBootstrapFilesystemKeyProvider<Schedule>
        filesystem_provider(slice_dir);
    auto filesystem_output =
        std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapWithKeyProvider<Schedule>(
        *filesystem_output, *input, filesystem_provider);
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *filesystem_output, *key);
    require_close_param<M>(*decoded, *slots, 0.02,
                           "CKKS dense filesystem-slice bootstrap e2e");
    std::filesystem::remove_all(slice_dir);

    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapTimings timings;
    TFHEpp::CKKSDenseBootstrapTimed<Schedule>(*output, *input, *bootstrap_key,
                                              timings);
    if (timings.total_ms() <= 0.0 || timings.coeff_to_slot_ms <= 0.0 ||
        timings.real_evalmod_ms <= 0.0 || timings.slot_to_coeff_ms <= 0.0)
        std::exit(1);

    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *output, *key);
    require_close_param<M>(*decoded, *slots, 0.02,
                           "CKKS dense encrypted bootstrap e2e");
}

void test_dense_bootstrap_inverse_e2e_smoke()
{
    using M = TinyDeepMultiLimbCKKSParam;
    using Schedule = TFHEpp::CKKSDenseBootstrapSchedule<
        M, 30, 8, 560, 30, 3, 31, 4, 2, 5, 30, 2>;
    static_assert(Schedule::evalmod_inv_degree == 5);
    static_assert(Schedule::evalmod_depth == 12);
    static_assert(Schedule::output_log_q == 110);

    std::array<double, M::n> coeffs{};
    for (std::size_t i = 0; i < M::n; i++)
        coeffs[i] =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 64.0;

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    coeffs_to_slots<M>(*slots, coeffs);

    auto key = std::make_unique<TFHEpp::Key<M>>();
    fill_test_key<M>(*key);

    auto bootstrap_key =
        std::make_unique<TFHEpp::CKKSDenseBootstrapKey<Schedule>>();
    TFHEpp::CKKSDenseBootstrapKeyGen<Schedule>(*bootstrap_key, *key,
                                               {0.0, 0});

    auto input = std::make_unique<typename Schedule::InputCiphertext>();
    TFHEpp::ckksSlotEncrypt<M, Schedule::input_log_q, Schedule::log_delta>(
        *input, *slots, *key, {0.0, 0});

    CountingDenseBootstrapProvider<Schedule> counting_provider(*bootstrap_key);
    auto output = std::make_unique<typename Schedule::OutputCiphertext>();
    TFHEpp::CKKSDenseBootstrapWithKeyProvider<Schedule>(
        *output, *input, counting_provider);

    if (counting_provider.coeff_to_slot_gets[0] != 1 ||
        counting_provider.coeff_to_slot_gets[1] != 1 ||
        counting_provider.slot_to_coeff_gets[0] != 1 ||
        counting_provider.slot_to_coeff_gets[1] != 1 ||
        counting_provider.packed_conjugate_gets != 2 ||
        counting_provider.evalmod_polynomial_gets != 2)
        std::exit(1);
    if (counting_provider.polynomial_relin_gets[0] != 2 ||
        counting_provider.polynomial_relin_gets[1] != 4 ||
        counting_provider.polynomial_relin_gets[2] != 8 ||
        counting_provider.polynomial_relin_gets[3] != 16 ||
        counting_provider.polynomial_relin_gets[4] != 30 ||
        counting_provider.double_angle_relin_gets[0] != 2 ||
        counting_provider.double_angle_relin_gets[1] != 2 ||
        counting_provider.inverse_relin_gets[0] != 2 ||
        counting_provider.inverse_relin_gets[1] != 4 ||
        counting_provider.inverse_relin_gets[2] != 2)
        std::exit(1);
    if (counting_provider.coeff_to_slot_releases[0] != 1 ||
        counting_provider.coeff_to_slot_releases[1] != 1 ||
        counting_provider.slot_to_coeff_releases[0] != 1 ||
        counting_provider.slot_to_coeff_releases[1] != 1 ||
        counting_provider.packed_conjugate_releases != 1)
        std::exit(1);
    if (counting_provider.polynomial_relin_releases[0] != 2 ||
        counting_provider.polynomial_relin_releases[1] != 2 ||
        counting_provider.polynomial_relin_releases[2] != 2 ||
        counting_provider.polynomial_relin_releases[3] != 2 ||
        counting_provider.polynomial_relin_releases[4] != 2 ||
        counting_provider.double_angle_relin_releases[0] != 2 ||
        counting_provider.double_angle_relin_releases[1] != 2 ||
        counting_provider.inverse_relin_releases[0] != 2 ||
        counting_provider.inverse_relin_releases[1] != 2 ||
        counting_provider.inverse_relin_releases[2] != 2)
        std::exit(1);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    TFHEpp::ckksSlotDecrypt<M, Schedule::output_log_q, Schedule::log_delta>(
        *decoded, *output, *key);
    require_close_param<M>(*decoded, *slots, 0.05,
                           "CKKS dense inverse encrypted bootstrap e2e");
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
    if (max_sine_err > 1e-6) std::exit(1);
    if (max_message_err > 2e-4) std::exit(1);
}

void test_bounded_cos_evalmod_inverse_plain()
{
    using L = TinyDeepMultiLimbCKKSParam;
    using Schedule =
        TFHEpp::CKKSDenseBootstrapSchedule<L, 25, 8, 500, 25, 5, 34, 18, 3, 5>;
    const auto poly = TFHEpp::CKKSBuildBoundedCosEvalModPolynomial<Schedule>();

    if (poly.sqrt_coeff != 1.0) std::exit(1);
    if (TFHEpp::CKKSBuildEvalModInversePowerCoefficients(
            poly, Schedule::evalmod_inv_degree)
            .size() != Schedule::evalmod_inv_degree + 1)
        std::exit(1);

    double max_message_err = 0.0;
    double max_power_err = 0.0;
    for (int mask = -static_cast<int>(Schedule::evalmod_k) + 1;
         mask < static_cast<int>(Schedule::evalmod_k); mask++) {
        for (const double msg : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
            const double masked =
                static_cast<double>(mask) * Schedule::message_ratio + msg;
            const double got = TFHEpp::CKKSPlainEvalModBoundedCos(
                poly, masked, Schedule::evalmod_inv_degree);
            const double got_power = TFHEpp::CKKSPlainEvalModBoundedCosPower(
                poly, masked, Schedule::evalmod_inv_degree);
            max_power_err = std::max(max_power_err, std::abs(got - got_power));
            max_message_err = std::max(max_message_err, std::abs(got - msg));
        }
    }

    std::cout << "CKKS inverse EvalMod plain message max_error="
              << max_message_err << " power max_error=" << max_power_err
              << std::endl;
    if (max_power_err > 1e-7) std::exit(1);
    if (max_message_err > 1e-7) std::exit(1);
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

void test_bounded_cos_evalmod_homomorphic()
{
    using M = TinyMultiLimbCKKSParam;
    constexpr std::uint32_t log_q = 300;
    constexpr std::uint32_t log_delta = 40;
    constexpr std::uint32_t coeff_log_delta = 20;
    constexpr std::size_t degree = 7;
    constexpr std::uint32_t k = 2;
    constexpr std::uint32_t log_message_ratio = 8;
    constexpr std::uint32_t double_angle = 2;
    using EvalCt = TFHEpp::CKKSCiphertext<M, log_q, log_delta>;
    using Traits =
        TFHEpp::CKKSEvalModBoundedCosTraits<M, log_q, log_delta,
                                            coeff_log_delta, degree,
                                            double_angle>;
    using EvalOut =
        TFHEpp::CKKSEvalModBoundedCosResult<M, log_q, log_delta,
                                            coeff_log_delta, degree,
                                            double_angle>;
    static_assert(Traits::polynomial_log_q ==
                  log_q - Traits::PolynomialTraits::power_depth * log_delta -
                      coeff_log_delta);
    static_assert(EvalOut::log_q ==
                  log_q - Traits::PolynomialTraits::power_depth * log_delta -
                      coeff_log_delta - double_angle * log_delta);

    const auto poly = TFHEpp::CKKSBuildBoundedCosEvalModPolynomial(
        k, degree, log_message_ratio, double_angle);
    const double normalizer =
        static_cast<double>(k) * poly.message_ratio * poly.q_diff;

    auto key = std::make_unique<TFHEpp::Key<M>>();
    fill_test_key<M>(*key);

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    slots->fill({0.0, 0.0});
    expected->fill({0.0, 0.0});
    for (std::size_t i = 0; i < M::n / 2; i++) {
        const int mask = static_cast<int>(i % (2 * k - 1)) -
                         static_cast<int>(k) + 1;
        const double msg =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 4.0;
        const double masked = static_cast<double>(mask) * poly.message_ratio + msg;
        (*slots)[i] = {masked / normalizer, 0.0};
        (*expected)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosPower(poly, masked) /
                 poly.message_ratio,
             0.0};
    }

    auto keys = std::make_unique<
        TFHEpp::CKKSEvalModBoundedCosRelinKeys<M, log_q, log_delta,
                                               coeff_log_delta, degree,
                                               double_angle>>();
    TFHEpp::CKKSEvalModBoundedCosKeyGen<M, log_q, log_delta, coeff_log_delta,
                                        degree, double_angle>(*keys, *key,
                                                              {0.0, 0});

    auto ct = std::make_unique<EvalCt>();
    TFHEpp::ckksSlotEncrypt<M, log_q, log_delta>(*ct, *slots, *key, {0.0, 0});

    auto out = std::make_unique<EvalOut>();
    TFHEpp::CKKSEvalModBoundedCosNormalized<M, log_q, log_delta,
                                            coeff_log_delta, degree,
                                            double_angle>(*out, *ct, poly,
                                                          *keys);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    TFHEpp::ckksSlotDecrypt<M, EvalOut::log_q, EvalOut::log_delta>(
        *decoded, *out, *key);
    require_close_param<M>(*decoded, *expected, 0.03,
                           "CKKS bounded cosine EvalMod homomorphic");
}

void test_bounded_cos_evalmod_inverse_homomorphic()
{
    using M = TinyMultiLimbCKKSParam;
    constexpr std::uint32_t log_q = 300;
    constexpr std::uint32_t log_delta = 25;
    constexpr std::uint32_t coeff_log_delta = 15;
    constexpr std::size_t degree = 7;
    constexpr std::uint32_t k = 2;
    constexpr std::uint32_t log_message_ratio = 8;
    constexpr std::uint32_t double_angle = 2;
    constexpr std::size_t inverse_degree = 5;
    using EvalCt = TFHEpp::CKKSCiphertext<M, log_q, log_delta>;
    using Traits =
        TFHEpp::CKKSEvalModBoundedCosTraits<M, log_q, log_delta,
                                            coeff_log_delta, degree,
                                            double_angle, inverse_degree>;
    using EvalOut =
        TFHEpp::CKKSEvalModBoundedCosResult<M, log_q, log_delta,
                                            coeff_log_delta, degree,
                                            double_angle, inverse_degree>;
    static_assert(Traits::inverse_enabled);
    static_assert(Traits::InverseTraits::power_depth == 3);
    static_assert(EvalOut::log_q ==
                  Traits::after_double_angle_log_q -
                      Traits::InverseTraits::power_depth * log_delta -
                      coeff_log_delta);

    const auto poly = TFHEpp::CKKSBuildBoundedCosEvalModPolynomial(
        k, degree, log_message_ratio, double_angle, 1.0, true);
    const double normalizer =
        static_cast<double>(k) * poly.message_ratio * poly.q_diff;

    auto key = std::make_unique<TFHEpp::Key<M>>();
    fill_test_key<M>(*key);

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    slots->fill({0.0, 0.0});
    expected->fill({0.0, 0.0});
    for (std::size_t i = 0; i < M::n / 2; i++) {
        const int mask =
            static_cast<int>(i % (2 * k - 1)) - static_cast<int>(k) + 1;
        const double msg =
            static_cast<double>(static_cast<int>(i % 5) - 2) / 4.0;
        const double masked =
            static_cast<double>(mask) * poly.message_ratio + msg;
        (*slots)[i] = {masked / normalizer, 0.0};
        (*expected)[i] =
            {TFHEpp::CKKSPlainEvalModBoundedCosPower(poly, masked,
                                                     inverse_degree) /
                 poly.message_ratio,
             0.0};
    }

    auto keys = std::make_unique<
        TFHEpp::CKKSEvalModBoundedCosRelinKeys<M, log_q, log_delta,
                                               coeff_log_delta, degree,
                                               double_angle, inverse_degree>>();
    TFHEpp::CKKSEvalModBoundedCosKeyGen<M, log_q, log_delta, coeff_log_delta,
                                        degree, double_angle, inverse_degree>(
        *keys, *key, {0.0, 0});

    auto ct = std::make_unique<EvalCt>();
    TFHEpp::ckksSlotEncrypt<M, log_q, log_delta>(*ct, *slots, *key, {0.0, 0});

    auto out = std::make_unique<EvalOut>();
    TFHEpp::CKKSEvalModBoundedCosNormalized<M, log_q, log_delta,
                                            coeff_log_delta, degree,
                                            double_angle, inverse_degree>(
        *out, *ct, poly, *keys);

    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<M>>();
    TFHEpp::ckksSlotDecrypt<M, EvalOut::log_q, EvalOut::log_delta>(
        *decoded, *out, *key);
    require_close_param<M>(*decoded, *expected, 0.05,
                           "CKKS bounded cosine inverse EvalMod homomorphic");
}

}  // namespace

int main()
{
    test_seeded_key_switch_rows();
    test_dense_coeff_slot_diagonals();
    test_factorized_coeff_slot_stages();
    test_lvl6_factorized_stage_shape();
    test_popcount3_direct_rotation();
    test_dense_bootstrap_api_shape();
    test_dense_bootstrap_plain_split_pipeline();
    test_dense_bootstrap_encrypted_pipeline();
    test_dense_bootstrap_fused_stc_shared_tail();
    test_dense_bootstrap_e2e_smoke();
    test_dense_bootstrap_inverse_e2e_smoke();
    test_bounded_cos_evalmod_plain();
    test_bounded_cos_evalmod_inverse_plain();
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
    test_bounded_cos_evalmod_homomorphic();
    test_bounded_cos_evalmod_inverse_homomorphic();

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
