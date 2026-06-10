#if defined(USE_80BIT_SECURITY) || defined(USE_COMPRESS) || \
    defined(USE_CGGI19) || defined(USE_CONCRETE) || defined(USE_TFHE_RS) || \
    defined(USE_TERNARY)

#include <iostream>

int main()
{
    std::cout << "CKKS multi-limb torus test skipped for non-default "
                 "parameters"
              << std::endl;
    return 0;
}

#else

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include <tfhe++.hpp>

namespace {

using TinyTorus = TFHEpp::MultiLimbUInt<3>;

struct TinyCKKSMLParam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 3;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 2;
    static constexpr std::uint32_t lₐ = 2;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t Bgₐbit = 16;
    using T = TinyTorus;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bgₐ = T{1} << Bgₐbit;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static constexpr double α = 0.0;
    static constexpr T μ = T{1} << 189;
    static constexpr uint64_t plain_modulus_u64 = 257;
    static constexpr uint32_t plain_modulusbit = 9;
    static constexpr T plain_modulus = T{plain_modulus_u64};
    static constexpr double Δ = 0.0;
    static constexpr T delta_int =
        std::numeric_limits<T>::max() / plain_modulus_u64;
    static constexpr uint64_t Q_mod_t =
        (std::numeric_limits<T>::max() % plain_modulus_u64) + 1;
    static constexpr std::uint32_t l̅ = 12;
    static constexpr std::uint32_t l̅ₐ = 12;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;
};

constexpr uint64_t pow2_mod(unsigned bits, uint64_t mod)
{
    uint64_t r = 1 % mod;
    for (unsigned i = 0; i < bits; i++)
        r = static_cast<uint64_t>((static_cast<unsigned __int128>(r) * 2) %
                                  mod);
    return r;
}

void require(bool ok, const char *message)
{
    if (!ok) {
        std::cerr << "FAIL: " << message << std::endl;
        std::exit(1);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
typename P::T encode_integer(double value)
{
    const auto scaled = static_cast<__int128_t>(
        std::llround(std::ldexp(value, LogDelta)));
    return TFHEpp::ckks_detail::signedToLevel<P, LogQ>(scaled);
}

template <class P>
void fill_diagnostic_key(TFHEpp::Key<P> &key)
{
    for (std::uint32_t i = 0; i < P::k * P::n; i++)
        key[i] = P::T::from_signed_i64(static_cast<int>(i % 3) - 1);
}

template <std::uint32_t LogQ, std::uint32_t RelinBgbit = TFHEpp::lvl6param::Bgbit,
          std::uint32_t RelinBbarbit = TFHEpp::lvl6param::B̅gbit>
void run_lvl6_seeded_dd_relin_diagnostic_at(const char *label)
{
    using P = TFHEpp::lvl6param;
    constexpr std::uint32_t log_delta = 40;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_diagnostic_key<P>(*key);

    auto coeffs = std::make_unique<std::array<double, P::n>>();
    coeffs->fill(0.0);
    (*coeffs)[0] = 0.125;
    (*coeffs)[1] = -0.0625;
    (*coeffs)[17] = 0.03125;
    (*coeffs)[P::n - 2] = -0.015625;

    auto poly = std::make_unique<TFHEpp::Polynomial<P>>();
    TFHEpp::ckksEncodePolynomial<P, LogQ, log_delta>(*poly, *coeffs);

    auto relinkey = TFHEpp::makeCKKSRelinKey<P, LogQ>(*key, {0.0, 0});
    auto seeded_dd_relinkey =
        TFHEpp::makeCKKSSeededDDRelinKey<P, LogQ, RelinBgbit, RelinBbarbit>(
            *key, {0.0, 0});

    auto standard = std::make_unique<TFHEpp::TRLWE<P>>();
    auto seeded_dd = std::make_unique<TFHEpp::TRLWE<P>>();
    TFHEpp::CKKSRelinKeySwitch<P>(*standard, *poly, *relinkey);
    TFHEpp::CKKSRelinKeySwitch<P>(*seeded_dd, *poly, *seeded_dd_relinkey);

    auto standard_ct =
        std::make_unique<TFHEpp::CKKSCiphertext<P, LogQ, log_delta>>();
    auto seeded_dd_ct =
        std::make_unique<TFHEpp::CKKSCiphertext<P, LogQ, log_delta>>();
    standard_ct->ct = *standard;
    seeded_dd_ct->ct = *seeded_dd;

    auto standard_decoded = std::make_unique<std::array<double, P::n>>();
    auto seeded_dd_decoded = std::make_unique<std::array<double, P::n>>();
    TFHEpp::ckksDecrypt<P, LogQ, log_delta>(*standard_decoded, *standard_ct,
                                            *key);
    TFHEpp::ckksDecrypt<P, LogQ, log_delta>(*seeded_dd_decoded, *seeded_dd_ct,
                                            *key);

    double max_delta = 0.0;
    for (std::uint32_t i = 0; i < P::n; i++) {
        max_delta = std::max(
            max_delta,
            std::abs((*standard_decoded)[i] - (*seeded_dd_decoded)[i]));
    }
    std::cout << label << "_log_q=" << LogQ
              << " standard_rows="
              << TFHEpp::CKKSRelinKeySwitchRowCount<P, LogQ>()
              << " dd_bgbit=" << RelinBgbit
              << " dd_bbarbit=" << RelinBbarbit
              << " dd_primary_rows="
              << TFHEpp::CKKSDDRelinPrimaryRowCountForLevel<P, LogQ,
                                                            RelinBgbit>()
              << " dd_bbar_rows="
              << TFHEpp::CKKSDDRelinBbarRowCountForLevel<P, LogQ,
                                                         RelinBbarbit>()
              << " max_delta=" << max_delta << '\n';
    require(max_delta < 0.01, "lvl6 seeded DD relin matches standard relin");
}

template <std::uint32_t LogQ, std::uint32_t RelinBgbit = TFHEpp::lvl6param::Bgbit,
          std::uint32_t RelinBbarbit = TFHEpp::lvl6param::B̅gbit>
void run_lvl6_seeded_dd_mult_diagnostic_at(const char *label,
                                           TFHEpp::CKKSNoise noise)
{
    using P = TFHEpp::lvl6param;
    constexpr std::uint32_t log_delta = 40;
    using Ct = TFHEpp::CKKSCiphertext<P, LogQ, log_delta>;
    using ProductCt =
        TFHEpp::CKKSMultResult<P, LogQ, log_delta, LogQ, log_delta>;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_diagnostic_key<P>(*key);

    auto lhs = std::make_unique<std::array<double, P::n>>();
    auto rhs = std::make_unique<std::array<double, P::n>>();
    auto expected = std::make_unique<std::array<double, P::n>>();
    lhs->fill(0.0);
    rhs->fill(0.0);
    expected->fill(0.0);
    (*lhs)[0] = 0.125;
    (*lhs)[1] = -0.0625;
    (*lhs)[17] = 0.03125;
    (*lhs)[P::n - 2] = -0.015625;
    (*rhs)[0] = -0.1875;
    (*rhs)[2] = 0.078125;
    (*rhs)[19] = -0.0390625;
    (*rhs)[P::n - 3] = 0.0234375;

    const auto add_product = [&](std::uint32_t ai, double av, std::uint32_t bi,
                                 double bv) {
        const std::uint32_t raw = ai + bi;
        if (raw < P::n)
            (*expected)[raw] += av * bv;
        else
            (*expected)[raw - P::n] -= av * bv;
    };
    const std::array<std::pair<std::uint32_t, double>, 4> lhs_terms = {
        {{0, (*lhs)[0]}, {1, (*lhs)[1]}, {17, (*lhs)[17]},
         {P::n - 2, (*lhs)[P::n - 2]}}};
    const std::array<std::pair<std::uint32_t, double>, 4> rhs_terms = {
        {{0, (*rhs)[0]}, {2, (*rhs)[2]}, {19, (*rhs)[19]},
         {P::n - 3, (*rhs)[P::n - 3]}}};
    for (const auto &[ai, av] : lhs_terms)
        for (const auto &[bi, bv] : rhs_terms) add_product(ai, av, bi, bv);

    auto lhs_ct = std::make_unique<Ct>();
    auto rhs_ct = std::make_unique<Ct>();
    TFHEpp::ckksEncrypt<P, LogQ, log_delta>(*lhs_ct, *lhs, *key, noise);
    TFHEpp::ckksEncrypt<P, LogQ, log_delta>(*rhs_ct, *rhs, *key, noise);

    auto relinkey = TFHEpp::makeCKKSRelinKey<P, ProductCt::log_q>(*key, noise);
    auto seeded_dd_relinkey =
        TFHEpp::makeCKKSSeededDDRelinKey<P, ProductCt::log_q, RelinBgbit,
                                         RelinBbarbit>(*key, noise);

    auto standard = std::make_unique<ProductCt>();
    auto seeded_dd = std::make_unique<ProductCt>();
    TFHEpp::CKKSMult<P>(*standard, *lhs_ct, *rhs_ct, *relinkey);
    TFHEpp::CKKSMult<P>(*seeded_dd, *lhs_ct, *rhs_ct, *seeded_dd_relinkey);

    auto standard_decoded = std::make_unique<std::array<double, P::n>>();
    auto seeded_dd_decoded = std::make_unique<std::array<double, P::n>>();
    TFHEpp::ckksDecrypt<P>(*standard_decoded, *standard, *key);
    TFHEpp::ckksDecrypt<P>(*seeded_dd_decoded, *seeded_dd, *key);

    double standard_error = 0.0;
    double seeded_dd_error = 0.0;
    double max_delta = 0.0;
    for (std::uint32_t i = 0; i < P::n; i++) {
        standard_error =
            std::max(standard_error,
                     std::abs((*standard_decoded)[i] - (*expected)[i]));
        seeded_dd_error =
            std::max(seeded_dd_error,
                     std::abs((*seeded_dd_decoded)[i] - (*expected)[i]));
        max_delta = std::max(
            max_delta,
            std::abs((*standard_decoded)[i] - (*seeded_dd_decoded)[i]));
    }
    std::cout << label << "_log_q=" << LogQ
              << " product_log_q=" << ProductCt::log_q
              << " dd_bgbit=" << RelinBgbit
              << " dd_bbarbit=" << RelinBbarbit
              << " standard_error=" << standard_error
              << " seeded_dd_error=" << seeded_dd_error
              << " max_delta=" << max_delta << '\n';
    require(standard_error < 0.01, "lvl6 CKKS multiply standard relin");
    require(seeded_dd_error < 0.01, "lvl6 CKKS multiply seeded DD relin");
    require(max_delta < 0.01, "lvl6 CKKS multiply seeded DD matches standard");
}

}  // namespace

int main(int argc, char **argv)
{
    if (argc == 2 && std::string(argv[1]) == "--lvl6-dd-relin-low") {
        run_lvl6_seeded_dd_relin_diagnostic_at<168>("lvl6_dd_relin_low");
        std::cout << "Passed" << std::endl;
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--lvl6-dd-relin-mid") {
        run_lvl6_seeded_dd_relin_diagnostic_at<328>("lvl6_dd_relin_mid");
        std::cout << "Passed" << std::endl;
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--lvl6-dd-relin-high") {
        run_lvl6_seeded_dd_relin_diagnostic_at<728>("lvl6_dd_relin_high");
        std::cout << "Passed" << std::endl;
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--lvl6-dd-relin-high-tuned") {
        run_lvl6_seeded_dd_relin_diagnostic_at<
            728, TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule::dd_relin_bgbit,
            TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule::dd_relin_bbarbit>(
            "lvl6_dd_relin_high_tuned");
        std::cout << "Passed" << std::endl;
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--lvl6-dd-mult-high") {
        run_lvl6_seeded_dd_mult_diagnostic_at<768>("lvl6_dd_mult_high",
                                                   {0.0, 0});
        std::cout << "Passed" << std::endl;
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--lvl6-dd-mult-high-noisy") {
        run_lvl6_seeded_dd_mult_diagnostic_at<768>(
            "lvl6_dd_mult_high_noisy", {TFHEpp::lvl6param::α, 0});
        std::cout << "Passed" << std::endl;
        return 0;
    }
    if (argc == 2 &&
        std::string(argv[1]) == "--lvl6-dd-mult-high-noisy-tuned") {
        run_lvl6_seeded_dd_mult_diagnostic_at<
            768, TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule::dd_relin_bgbit,
            TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule::dd_relin_bbarbit>(
            "lvl6_dd_mult_high_noisy_tuned", {TFHEpp::lvl6param::α, 0});
        std::cout << "Passed" << std::endl;
        return 0;
    }

    using Lvl6T = TFHEpp::lvl6param::T;
    static_assert(TFHEpp::lvl5param::nbit == 14);
    static_assert(TFHEpp::lvl6param::nbit == 15);
    static_assert(std::numeric_limits<Lvl6T>::digits == 896);
    static_assert(TFHEpp::lvl6param::l̅ * TFHEpp::lvl6param::B̅gbit == 896);
    static_assert(TFHEpp::lvl6param::Q_mod_t ==
                  pow2_mod(896, TFHEpp::lvl6param::plain_modulus_u64));
    require(std::abs(std::log2(TFHEpp::lvl6param::α) + 886.0) < 1e-9,
            "lvl6 CKKS security noise exponent");

    {
        using P = TFHEpp::lvl6param;
        using T = typename P::T;
        using Ct = TFHEpp::CKKSCiphertext<P, 896, 50>;
        static_assert(Ct::log_q == 896);
        static_assert(Ct::log_budget == 846);
        static_assert(TFHEpp::is_multilimb_digit_fft_compatible_v<P>);

        const T level_neg_one = TFHEpp::ckks_detail::levelMask<P, 880>();
        const T centered_neg_one =
            TFHEpp::ckks_detail::centeredLevelToTorus<P, 880>(level_neg_one);
        require(centered_neg_one == T{-1}, "lvl6 active-level -1 sign extend");

        const T sign = T{1} << 879;
        const T centered_sign =
            TFHEpp::ckks_detail::centeredLevelToTorus<P, 880>(sign);
        require((centered_sign >> 895) == T{1},
                "lvl6 active sign extends to storage top bit");
        require(TFHEpp::ckks_detail::reduceToLevel<P, 880>(centered_sign) ==
                    sign,
                "lvl6 centered value reduces back to active level");

        const T uniform = TFHEpp::ckks_detail::uniformAtLevel<P, 880>();
        require((uniform >> 880) == T{0}, "lvl6 uniform sample is level-bound");

        const TFHEpp::CKKSNoise default_noise{P::α, 0};
        const long double input_noise =
            TFHEpp::ckks_detail::effectiveNoiseStddevAtLevel<P, 60>(
                default_noise);
        require(std::abs(input_noise - P::ckks_min_noise_stddev) < 1e-12L,
                "lvl6 low active-level CKKS noise floor");
        const long double boot_noise =
            TFHEpp::ckks_detail::effectiveNoiseStddevAtLevel<P, 888>(
                default_noise);
        require(std::abs(std::log2(boot_noise) - 2.0L) < 1e-9L,
                "lvl6 boot-level CKKS noise scale");
        const T large_noise =
            TFHEpp::ckks_detail::signedLongDoubleToLevel<P, 888>(
                std::ldexp(1.0L, 2));
        require((large_noise >> 2) == T{1},
                "lvl6 large CKKS noise conversion keeps high bits");

        auto lhs = std::make_unique<TFHEpp::Polynomial<P>>();
        auto rhs = std::make_unique<TFHEpp::Polynomial<P>>();
        auto product = std::make_unique<TFHEpp::Polynomial<P>>();
        (*lhs)[0] = T::from_signed_i64(3);
        (*lhs)[1] = T::from_signed_i64(-1);
        (*lhs)[P::n - 1] = T::from_signed_i64(2);
        (*rhs)[0] = T::from_signed_i64(2);
        (*rhs)[2] = T::from_signed_i64(5);
        (*rhs)[P::n - 1] = T::from_signed_i64(-4);
        TFHEpp::PolyMulDigit<P>(*product, *lhs, *rhs);

        require((*product)[0] == T::from_signed_i64(2),
                "lvl6 digit FFT coefficient 0");
        require((*product)[1] == T::from_signed_i64(-12),
                "lvl6 digit FFT coefficient 1");
        require((*product)[2] == T::from_signed_i64(15),
                "lvl6 digit FFT coefficient 2");
        require((*product)[3] == T::from_signed_i64(-5),
                "lvl6 digit FFT coefficient 3");
        require((*product)[P::n - 2] == T::from_signed_i64(8),
                "lvl6 digit FFT coefficient n-2");
        require((*product)[P::n - 1] == T::from_signed_i64(-8),
                "lvl6 digit FFT coefficient n-1");
        for (std::uint32_t i = 4; i + 2 < P::n; i++)
            require((*product)[i] == T{0}, "lvl6 digit FFT zero coefficient");

        auto dense_lhs = std::make_unique<TFHEpp::Polynomial<P>>();
        auto dense_rhs = std::make_unique<TFHEpp::Polynomial<P>>();
        auto dense_product = std::make_unique<TFHEpp::Polynomial<P>>();
        constexpr int64_t bg_limit =
            (int64_t{1} << (P::Bgbit - 1)) - int64_t{1};
        constexpr int64_t bbar_limit =
            (int64_t{1} << (P::B̅gbit - 1)) - int64_t{1};
        for (std::uint32_t i = 0; i < P::n; i++) {
            const int64_t lhs_value = (i & 1) == 0 ? bg_limit : -bg_limit;
            const int64_t rhs_value =
                (i % 3) == 0 ? bbar_limit : -bbar_limit;
            (*dense_lhs)[i] = T::from_signed_i64(lhs_value);
            (*dense_rhs)[i] = T::from_signed_i64(rhs_value);
        }
        TFHEpp::PolyMulDigit<P>(*dense_product, *dense_lhs, *dense_rhs);

        auto expected_digit_product = [&](std::uint32_t i) {
            __int128 expected = 0;
            for (std::uint32_t j = 0; j <= i; j++) {
                expected +=
                    static_cast<__int128>(TFHEpp::multilimb_to_signed_i64(
                        (*dense_lhs)[j])) *
                    TFHEpp::multilimb_to_signed_i64((*dense_rhs)[i - j]);
            }
            for (std::uint32_t j = i + 1; j < P::n; j++) {
                expected -=
                    static_cast<__int128>(TFHEpp::multilimb_to_signed_i64(
                        (*dense_lhs)[j])) *
                    TFHEpp::multilimb_to_signed_i64(
                        (*dense_rhs)[P::n + i - j]);
            }
            return static_cast<int64_t>(expected);
        };
        for (const std::uint32_t i :
             {0U, 1U, 2U, 17U, P::n / 2U, P::n - 1U}) {
            require(TFHEpp::multilimb_to_signed_i64((*dense_product)[i]) ==
                        expected_digit_product(i),
                    "lvl6 DD relin digit FFT dense coefficient");
        }
    }

    {
        using P = TinyCKKSMLParam;
        constexpr std::uint32_t log_q = 180;
        constexpr std::uint32_t log_delta = 20;
        constexpr std::uint32_t out_log_q = log_q - log_delta;

        auto lhs = std::make_unique<TFHEpp::TRLWE<P>>();
        auto rhs = std::make_unique<TFHEpp::TRLWE<P>>();
        (*lhs)[P::k][0] = encode_integer<P, log_q, log_delta>(3.0);
        (*lhs)[P::k][1] = encode_integer<P, log_q, log_delta>(1.0);
        (*rhs)[P::k][0] = encode_integer<P, log_q, log_delta>(-2.0);
        (*rhs)[P::k][2] = encode_integer<P, log_q, log_delta>(4.0);

        auto product = std::make_unique<TFHEpp::TRLWE3<P>>();
        TFHEpp::CKKSTensorProductRescale<P, log_q, log_q, log_delta>(
            *product, *lhs, *rhs);

        require((*product)[1][0] ==
                    encode_integer<P, out_log_q, log_delta>(-6.0),
                "multi-limb CKKS tensor coefficient 0");
        require((*product)[1][1] ==
                    encode_integer<P, out_log_q, log_delta>(-2.0),
                "multi-limb CKKS tensor coefficient 1");
        require((*product)[1][2] ==
                    encode_integer<P, out_log_q, log_delta>(12.0),
                "multi-limb CKKS tensor coefficient 2");
        require((*product)[1][3] ==
                    encode_integer<P, out_log_q, log_delta>(4.0),
                "multi-limb CKKS tensor coefficient 3");
        for (std::uint32_t i = 0; i < P::n; i++) {
            require((*product)[0][i] == typename P::T{0},
                    "multi-limb CKKS tensor cross term");
            require((*product)[2][i] == typename P::T{0},
                    "multi-limb CKKS tensor square term");
        }

        TFHEpp::Key<P> key{};
        auto relinkey = TFHEpp::makeCKKSRelinKey<P, out_log_q>(key, {0.0, 0});
        TFHEpp::TRLWE<P> switched{};
        TFHEpp::Polynomial<P> zero{};
        TFHEpp::CKKSRelinKeySwitch<P>(switched, zero, *relinkey);
    }

    {
        using P = TinyCKKSMLParam;
        constexpr std::uint32_t log_q = 180;
        constexpr std::uint32_t log_delta = 20;
        using Ct = TFHEpp::CKKSCiphertext<P, log_q, log_delta>;
        using ProductCt =
            TFHEpp::CKKSMultResult<P, log_q, log_delta, log_q, log_delta>;
        static_assert(ProductCt::log_q == 160);
        static_assert(TFHEpp::CKKSDDRelinKey<P, ProductCt::log_q>::
                          primary_rows == 10);
        static_assert(
            TFHEpp::CKKSDDRelinKey<P, ProductCt::log_q>::bbar_rows == 10);

        TFHEpp::Key<P> key{};
        for (std::uint32_t i = 0; i < P::n; i++)
            key[i] =
                P::T::from_signed_i64(static_cast<int>(i % 3) - 1);

        auto relinkey =
            TFHEpp::makeCKKSRelinKey<P, ProductCt::log_q>(key, {0.0, 0});
        auto dd_relinkey =
            TFHEpp::makeCKKSDDRelinKey<P, ProductCt::log_q>(key, {0.0, 0});
        auto dd_chain =
            std::make_unique<TFHEpp::CKKSDDRelinKeyChain<P, log_q, log_delta,
                                                         1>>();
        TFHEpp::CKKSDDRelinKeyChainGen<P, log_q, log_delta, 1>(*dd_chain, key,
                                                               {0.0, 0});

        auto lhs = std::make_unique<std::array<double, P::n>>();
        auto rhs = std::make_unique<std::array<double, P::n>>();
        auto expected = std::make_unique<std::array<double, P::n>>();
        lhs->fill(0.0);
        rhs->fill(0.0);
        expected->fill(0.0);
        (*lhs)[0] = 0.25;
        (*lhs)[1] = -0.125;
        (*lhs)[3] = 0.0625;
        (*rhs)[0] = -0.5;
        (*rhs)[2] = 0.25;

        auto add_product = [&](std::uint32_t ai, double av, std::uint32_t bi,
                               double bv) {
            const std::uint32_t raw = ai + bi;
            if (raw < P::n)
                (*expected)[raw] += av * bv;
            else
                (*expected)[raw - P::n] -= av * bv;
        };
        add_product(0, (*lhs)[0], 0, (*rhs)[0]);
        add_product(0, (*lhs)[0], 2, (*rhs)[2]);
        add_product(1, (*lhs)[1], 0, (*rhs)[0]);
        add_product(1, (*lhs)[1], 2, (*rhs)[2]);
        add_product(3, (*lhs)[3], 0, (*rhs)[0]);
        add_product(3, (*lhs)[3], 2, (*rhs)[2]);

        auto lhs_ct = std::make_unique<Ct>();
        auto rhs_ct = std::make_unique<Ct>();
        TFHEpp::ckksEncrypt<P, log_q, log_delta>(*lhs_ct, *lhs, key, {0.0, 0});
        TFHEpp::ckksEncrypt<P, log_q, log_delta>(*rhs_ct, *rhs, key, {0.0, 0});

        auto product = std::make_unique<ProductCt>();
        auto dd_product = std::make_unique<ProductCt>();
        auto dd_chain_product = std::make_unique<ProductCt>();
        TFHEpp::CKKSMult<P>(*product, *lhs_ct, *rhs_ct, *relinkey);
        TFHEpp::CKKSMult<P>(*dd_product, *lhs_ct, *rhs_ct, *dd_relinkey);
        TFHEpp::CKKSMult<P>(*dd_chain_product, *lhs_ct, *rhs_ct,
                            dd_chain->template get<0>());

        auto got = std::make_unique<std::array<double, P::n>>();
        auto got_dd = std::make_unique<std::array<double, P::n>>();
        auto got_dd_chain = std::make_unique<std::array<double, P::n>>();
        TFHEpp::ckksDecrypt<P>(*got, *product, key);
        TFHEpp::ckksDecrypt<P>(*got_dd, *dd_product, key);
        TFHEpp::ckksDecrypt<P>(*got_dd_chain, *dd_chain_product, key);
        for (std::uint32_t i = 0; i < P::n; i++) {
            require(std::abs((*got)[i] - (*expected)[i]) < 0.05,
                    "tiny CKKS standard relin encrypted multiply");
            require(std::abs((*got_dd)[i] - (*expected)[i]) < 0.05,
                    "tiny CKKS DD relin encrypted multiply");
            require(std::abs((*got_dd)[i] - (*got)[i]) < 0.05,
                    "tiny CKKS DD relin matches standard relin");
            require(std::abs((*got_dd_chain)[i] - (*got)[i]) < 0.05,
                    "tiny CKKS DD relin chain matches standard relin");
        }
    }

    {
        using P = TinyCKKSMLParam;
        constexpr std::uint32_t log_q = 120;
        constexpr std::uint32_t log_delta = 20;
        const auto encoded = TFHEpp::ckksEncodeCoeff<P, log_q, log_delta>(-1.25);
        const double decoded =
            TFHEpp::ckksDecodeCoeff<P, log_q, log_delta>(encoded);
        require(std::abs(decoded + 1.25) < 1e-9,
                "small multi-limb CKKS encode/decode");
    }

    std::cout << "Passed" << std::endl;
}

#endif
