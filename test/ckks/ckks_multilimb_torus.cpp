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

}  // namespace

int main()
{
    using Lvl6T = TFHEpp::lvl6param::T;
    static_assert(TFHEpp::lvl5param::nbit == 14);
    static_assert(TFHEpp::lvl6param::nbit == 15);
    static_assert(std::numeric_limits<Lvl6T>::digits == 1152);
    static_assert(TFHEpp::lvl6param::l̅ * TFHEpp::lvl6param::B̅gbit == 1152);
    static_assert(TFHEpp::lvl6param::Q_mod_t ==
                  pow2_mod(1152, TFHEpp::lvl6param::plain_modulus_u64));
    require(std::abs(std::log2(TFHEpp::lvl6param::α) + 850.0) < 1e-9,
            "lvl6 CKKS security noise exponent");

    {
        using P = TFHEpp::lvl6param;
        using T = typename P::T;
        using Ct = TFHEpp::CKKSCiphertext<P, 1152, 50>;
        static_assert(Ct::log_q == 1152);
        static_assert(Ct::log_budget == 1102);
        static_assert(TFHEpp::is_multilimb_digit_fft_compatible_v<P>);

        const T level_neg_one = TFHEpp::ckks_detail::levelMask<P, 880>();
        const T centered_neg_one =
            TFHEpp::ckks_detail::centeredLevelToTorus<P, 880>(level_neg_one);
        require(centered_neg_one == T{-1}, "lvl6 active-level -1 sign extend");

        const T sign = T{1} << 879;
        const T centered_sign =
            TFHEpp::ckks_detail::centeredLevelToTorus<P, 880>(sign);
        require((centered_sign >> 1151) == T{1},
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
            TFHEpp::ckks_detail::effectiveNoiseStddevAtLevel<P, 1152>(
                default_noise);
        require(std::abs(std::log2(boot_noise) - 302.0L) < 1e-9L,
                "lvl6 boot-level CKKS noise scale");
        const T large_noise =
            TFHEpp::ckks_detail::signedLongDoubleToLevel<P, 1152>(
                std::ldexp(1.0L, 302));
        require((large_noise >> 302) == T{1},
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
        TFHEpp::CKKSRelinKeySwitch<P, out_log_q>(switched, zero, *relinkey);
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
