#if defined(USE_80BIT_SECURITY) || defined(USE_COMPRESS) || \
    defined(USE_CGGI19) || defined(USE_CONCRETE) || defined(USE_TFHE_RS) || \
    defined(USE_TERNARY)

#include <iostream>

int main()
{
    std::cout << "BFV multi-limb torus scaffold test skipped for non-default "
                 "parameters"
              << std::endl;
    return 0;
}

#else

#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>

#include <tfhe++.hpp>

namespace {

using U128x = TFHEpp::MultiLimbUInt<2>;

struct MLTestParam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 3;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 2;
    static constexpr std::uint32_t lₐ = 2;
    static constexpr std::uint32_t Bgbit = 16;
    static constexpr std::uint32_t Bgₐbit = 16;
    using T = U128x;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr T Bgₐ = T{1} << Bgₐbit;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static constexpr double α = 0.0;
    static constexpr T μ = T{1} << 125;
    static constexpr uint64_t plain_modulus_u64 = 257;
    static constexpr uint32_t plain_modulusbit = 9;
    static constexpr T plain_modulus = T{plain_modulus_u64};
    static constexpr double Δ = 0.0;
    static constexpr T delta_int =
        std::numeric_limits<T>::max() / plain_modulus_u64;
    static constexpr uint64_t Q_mod_t =
        (std::numeric_limits<T>::max() % plain_modulus_u64) + 1;
    static constexpr std::uint32_t l̅ = 8;
    static constexpr std::uint32_t l̅ₐ = 8;
    static constexpr std::uint32_t B̅gbit = 16;
    static constexpr std::uint32_t B̅gₐbit = 16;
};

constexpr uint64_t pow2_mod(unsigned bits, uint64_t mod)
{
    uint64_t r = 1 % mod;
    for (unsigned i = 0; i < bits; i++)
        r = static_cast<uint64_t>((static_cast<unsigned __int128>(r) * 2) % mod);
    return r;
}

__uint128_t to_u128(const U128x &x)
{
    return (static_cast<__uint128_t>(x.limb[1]) << 64) | x.limb[0];
}

}  // namespace

int main()
{
    using U = TFHEpp::MultiLimbUInt<7>;
    static_assert(std::numeric_limits<U>::digits == 448);
    static_assert(TFHEpp::is_multilimb_uint_v<U>);
    static_assert(TFHEpp::lvl5param::nbit == 14);
    static_assert(TFHEpp::lvl5param::l̅ * TFHEpp::lvl5param::B̅gbit == 448);
    static_assert(TFHEpp::lvl5param::Q_mod_t ==
                  pow2_mod(448, TFHEpp::lvl5param::plain_modulus_u64));

    using BP = TFHEpp::bfvboot::PrimePower2Param<TFHEpp::lvl5param>;
    static_assert(BP::base_plain_modulus == TFHEpp::lvl5param::plain_modulus_u64);
    static_assert(BP::plain_modulus_u64 ==
                  TFHEpp::lvl5param::plain_modulus_u64 *
                      TFHEpp::lvl5param::plain_modulus_u64);
    static_assert(BP::Q_mod_t == pow2_mod(448, BP::plain_modulus_u64));
    [[maybe_unused]] auto lvl5_make_relin =
        &TFHEpp::makeRelinKeyFFT<TFHEpp::lvl5param>;
    [[maybe_unused]] auto lvl5_relin_fftgen =
        &TFHEpp::relinKeyFFTgen<TFHEpp::lvl5param>;

    U one = 1;
    U top = one << 447;
    auto require = [](bool ok, const char *message) {
        if (!ok) {
            std::cerr << "FAIL: " << message << std::endl;
            std::exit(1);
        }
    };
    require(top.bit(447), "top bit shift");
    require(static_cast<uint64_t>(top >> 447) == 1, "top bit shift back");
    require((U{-1} + U{1}) == U{0}, "wraparound add");

    U a = (U{1} << 64) + U{3};
    U prod = a * U{5};
    require(prod.limb[0] == 15, "multiply low limb");
    require(prod.limb[1] == 5, "multiply carry limb");

    {
        TFHEpp::WideSignedLimbAccumulator<4> acc;
        acc.add_shifted_i64(12345, 90);
        acc.add_shifted_i64(-77, 20);
        const __int128_t expected =
            ((static_cast<__int128_t>(12345) << 90) -
             (static_cast<__int128_t>(77) << 20)) /
            257;
        require(to_u128(acc.div_to_torus<2>(U128x{257})) ==
                    static_cast<__uint128_t>(expected),
                "wide accumulator positive division");
    }

    {
        TFHEpp::WideSignedLimbAccumulator<4> acc;
        acc.add_shifted_i64(-12345, 90);
        acc.add_shifted_i64(77, 20);
        const __int128_t expected =
            (-(static_cast<__int128_t>(12345) << 90) +
             (static_cast<__int128_t>(77) << 20)) /
            257;
        require(to_u128(acc.div_to_torus<2>(U128x{257})) ==
                    static_cast<__uint128_t>(expected),
                "wide accumulator negative division");
    }

    {
        using P = MLTestParam;
        TFHEpp::Polynomial<P> a{}, b{}, expected{}, actual{};
        a[0] = U128x{3};
        a[1] = U128x{-2};
        a[3] = U128x{11};
        b[0] = U128x{-1};
        b[1] = U128x{2};
        b[4] = U128x{1};
        TFHEpp::PolyMulNaive<P>(expected, a, b);
        TFHEpp::PolyMulDigit<P>(actual, a, b);
        for (std::size_t i = 0; i < P::n; i++)
            require(actual[i] == expected[i], "multi-limb digit PolyMul");
    }

    {
        using P = MLTestParam;
        TFHEpp::Polynomial<P> torus{}, digit{}, expected{}, actual{};
        torus[0] = (U128x{0x1234} << 96) + U128x{7};
        torus[1] = -(U128x{0x55} << 80) + U128x{5};
        torus[3] = U128x{0xabcdef};
        digit[0] = U128x{1};
        digit[1] = U128x{-1};
        digit[2] = U128x{2};
        TFHEpp::PolyMulNaive<P>(expected, torus, digit);
        TFHEpp::PolyMulTorusByDigit<P>(actual, torus, digit);
        for (std::size_t i = 0; i < P::n; i++)
            require(actual[i] == expected[i],
                    "multi-limb torus-by-digit PolyMul");
    }

    {
        using P = MLTestParam;
        TFHEpp::Polynomial<P> torus{}, plain_poly{}, expected{}, actual{};
        std::array<uint64_t, P::n> plain{};
        torus[0] = (U128x{0x1234} << 96) + U128x{7};
        torus[2] = -(U128x{0x31} << 88) + U128x{11};
        plain[0] = 257;
        plain[1] = 129;
        plain[5] = 42;
        for (std::size_t i = 0; i < P::n; i++) plain_poly[i] = U128x{plain[i]};
        TFHEpp::PolyMulNaive<P>(expected, torus, plain_poly);
        TFHEpp::PolyMulTorusByUnsigned<P>(actual, torus, plain);
        for (std::size_t i = 0; i < P::n; i++)
            require(actual[i] == expected[i],
                    "multi-limb torus-by-unsigned PolyMul");
    }

    const std::array<uint64_t, 4> decode_cases{
        0, 1, 42, MLTestParam::plain_modulus_u64 - 1};
    for (uint64_t m : decode_cases) {
        const U128x encoded = TFHEpp::bfvEncodeCoeff<MLTestParam>(m);
        require(TFHEpp::bfvDecodeCoeff<MLTestParam>(encoded) == m,
                "multi-limb BFV coefficient decode");
    }

    TFHEpp::TRLWE<MLTestParam> ct{};
    ct[0][0] = U128x{0x1234} << 112;
    ct[1][1] = -(U128x{0x1234} << 112);

    std::array<TFHEpp::TRLWE<MLTestParam>, MLTestParam::l̅> dec{};
    TFHEpp::TRLWEBaseBbarDecompose<MLTestParam>(dec, ct);
    require(static_cast<uint64_t>(dec[0][0][0]) == 0x1234,
            "positive top Bbar digit");
    require(TFHEpp::multilimb_to_signed_i64(dec[0][1][1]) == -0x1234,
            "negative top Bbar digit");
    for (std::size_t j = 1; j < MLTestParam::l̅; j++) {
        require(dec[j][0][0] == U128x{0}, "positive lower Bbar digit");
        require(dec[j][1][1] == U128x{0}, "negative lower Bbar digit");
    }

    {
        using P = TFHEpp::lvl5param;
        auto a = std::make_unique<TFHEpp::Polynomial<P>>();
        auto b = std::make_unique<TFHEpp::Polynomial<P>>();
        auto afd = std::make_unique<TFHEpp::PolynomialInFD<P>>();
        auto bfd = std::make_unique<TFHEpp::PolynomialInFD<P>>();
        auto prod_fd = std::make_unique<TFHEpp::PolynomialInFD<P>>();
        auto prod = std::make_unique<TFHEpp::Polynomial<P>>();

        (*a)[0] = static_cast<typename P::T>(3);
        (*a)[1] = static_cast<typename P::T>(-2);
        (*b)[0] = static_cast<typename P::T>(5);
        (*b)[1] = static_cast<typename P::T>(7);
        TFHEpp::TwistIFFTDigit<P>(*afd, *a);
        TFHEpp::TwistIFFTDigit<P>(*bfd, *b);
        TFHEpp::MulInFD<P::n>(*prod_fd, *afd, *bfd);
        TFHEpp::TwistFFTDigitProduct<P>(*prod, *prod_fd);

        const std::array<int64_t, 4> expected_digit_product{15, 11, -14, 0};
        for (std::size_t i = 0; i < expected_digit_product.size(); i++)
            require(TFHEpp::multilimb_to_signed_i64((*prod)[i]) ==
                        expected_digit_product[i],
                    "lvl5 digit FFT product");
        for (std::size_t i = expected_digit_product.size(); i < P::n; i++)
            require((*prod)[i] == static_cast<typename P::T>(0),
                    "lvl5 digit FFT zero tail");
    }

    {
        using P = TFHEpp::lvl5param;
        constexpr int target_i = 0;
        constexpr int target_j = 0;
        constexpr int shift =
            std::numeric_limits<typename P::T>::digits -
            (target_j + 1) * static_cast<int>(P::B̅gbit);

        auto relinkeyfft = std::make_unique<TFHEpp::relinKeyFFT<P>>();
        auto identity = std::make_unique<TFHEpp::Polynomial<P>>();
        (*identity)[0] = static_cast<typename P::T>(1);
        TFHEpp::TwistIFFTDigit<P>(
            (*relinkeyfft)[target_i * P::l̅ + target_j][P::k], *identity);

        auto poly = std::make_unique<TFHEpp::Polynomial<P>>();
        (*poly)[0] = static_cast<typename P::T>(3)
                     << (std::numeric_limits<typename P::T>::digits -
                         P::Bgbit);
        (*poly)[1] = static_cast<typename P::T>(-2)
                     << (std::numeric_limits<typename P::T>::digits -
                         P::Bgbit);

        auto switched = std::make_unique<TFHEpp::TRLWE<P>>();
        TFHEpp::relinKeySwitch<P>(*switched, *poly, *relinkeyfft);

        auto dec = std::make_unique<TFHEpp::DecomposedPolynomial<P>>();
        TFHEpp::Decomposition<P>(*dec, *poly);
        for (std::size_t i = 0; i < P::n; i++) {
            const typename P::T expected = (*dec)[target_i][i] << shift;
            require((*switched)[P::k][i] == expected,
                    "lvl5 relinKeySwitch body recombination");
            require((*switched)[0][i] == static_cast<typename P::T>(0),
                    "lvl5 relinKeySwitch mask zero");
        }
    }

    {
        using P = TFHEpp::lvl5param;
        auto slots = std::make_unique<std::array<uint64_t, P::n>>();
        auto decoded = std::make_unique<std::array<uint64_t, P::n>>();
        auto poly = std::make_unique<TFHEpp::Polynomial<P>>();
        for (std::size_t i = 0; i < P::n; i++)
            (*slots)[i] = (17 * i + 3) % P::plain_modulus_u64;
        TFHEpp::SlotEncode<P>(*poly, *slots);
        TFHEpp::SlotDecode<P>(*decoded, *poly);
        for (std::size_t i = 0; i < P::n; i++)
            require((*decoded)[i] == (*slots)[i], "lvl5 SlotEncode/SlotDecode");
    }

    if (std::getenv("TFHEPP_BFV_LVL5_KEYGEN_TEST") != nullptr) {
        using P = TFHEpp::lvl5param;
        auto key = std::make_unique<TFHEpp::Key<P>>();
        for (std::size_t i = 0; i < P::n; i++) {
            const int v = static_cast<int>(i % 3) - 1;
            (*key)[i] = static_cast<typename P::T>(v);
        }
        auto relinkey = TFHEpp::makeRelinKeyFFT<P>(*key);
        require(relinkey != nullptr, "lvl5 makeRelinKeyFFT allocation");

        auto slots_a = std::make_unique<std::array<uint64_t, P::n>>();
        auto slots_b = std::make_unique<std::array<uint64_t, P::n>>();
        auto decrypted = std::make_unique<std::array<uint64_t, P::n>>();
        for (std::size_t i = 0; i < P::n; i++) {
            (*slots_a)[i] = (5 * i + 1) % P::plain_modulus_u64;
            (*slots_b)[i] = (7 * i + 3) % P::plain_modulus_u64;
        }

        auto ct_a = std::make_unique<TFHEpp::TRLWE<P>>();
        auto ct_b = std::make_unique<TFHEpp::TRLWE<P>>();
        auto ct_mul = std::make_unique<TFHEpp::TRLWE<P>>();
        TFHEpp::trlweSlotEncrypt<P>(*ct_a, *slots_a, *key);
        TFHEpp::trlweSlotDecrypt<P>(*decrypted, *ct_a, *key);
        for (std::size_t i = 0; i < P::n; i++)
            require((*decrypted)[i] == (*slots_a)[i],
                    "lvl5 slot encrypt/decrypt");

        TFHEpp::trlweSlotEncrypt<P>(*ct_b, *slots_b, *key);
        TFHEpp::TRLWEMultFullDD<P>(*ct_mul, *ct_a, *ct_b, *relinkey);
        TFHEpp::trlweSlotDecrypt<P>(*decrypted, *ct_mul, *key);
        for (std::size_t i = 0; i < P::n; i++) {
            const uint64_t expected =
                ((*slots_a)[i] * (*slots_b)[i]) % P::plain_modulus_u64;
            require((*decrypted)[i] == expected, "lvl5 encrypted slot multiply");
        }

        using BP = TFHEpp::bfvboot::PrimePower2Param<P>;
        auto boot_key = std::make_unique<TFHEpp::Key<BP>>();
        TFHEpp::bfvboot::ConvertSecretKey<P, BP>(*boot_key, *key);
        auto sk_plain = std::make_unique<std::array<uint64_t, BP::n>>();
        TFHEpp::bfvboot::SecretKeyAsPlaintext<BP, P>(*sk_plain, *key);
        auto enc_sk = std::make_unique<TFHEpp::TRLWE<BP>>();
        TFHEpp::bfvboot::BfvPolyEncrypt<BP>(*enc_sk, *sk_plain, *boot_key);
        auto noisy = std::make_unique<TFHEpp::TRLWE<BP>>();
        TFHEpp::bfvboot::NoisyDecrypt<P, BP>(*noisy, *ct_a, *enc_sk);
        auto noisy_plain = std::make_unique<std::array<uint64_t, BP::n>>();
        TFHEpp::bfvboot::BfvPolyDecrypt<BP>(*noisy_plain, *noisy, *boot_key);

        if (std::getenv("TFHEPP_BFV_LVL5_BOOTSTRAP_TEST") != nullptr) {
            relinkey.reset();
            ct_b.reset();
            ct_mul.reset();
            sk_plain.reset();
            enc_sk.reset();
            noisy.reset();
            noisy_plain.reset();

            auto gk = std::make_unique<TFHEpp::GaloisKey<P>>();
            TFHEpp::GaloisKeyGen<P>(*gk, *key);

            constexpr int transform_terms = 9;
            std::vector<std::array<uint64_t, P::n>> diagonals(transform_terms);
            std::vector<int> offsets(transform_terms);
            for (int r = 0; r < transform_terms; r++) {
                diagonals[r].fill(1);
                offsets[r] = r;
            }

            auto transformed = std::make_unique<TFHEpp::TRLWE<P>>();
            TFHEpp::LinearTransformBSGS<P>(*transformed, *ct_a, diagonals,
                                           offsets, transform_terms, *gk);
            TFHEpp::trlweSlotDecrypt<P>(*decrypted, *transformed, *key);

            auto expected = std::make_unique<std::array<uint64_t, P::n>>();
            TFHEpp::c2s::PlainLinearTransform<P>(*expected, *slots_a,
                                                 diagonals);
            for (std::size_t i = 0; i < P::n; i++)
                require((*decrypted)[i] == (*expected)[i],
                        "lvl5 fused BSGS linear transform");

            transformed.reset();
            gk.reset();

            if (std::getenv("TFHEPP_BFV_LVL5_BOOTSTRAP_FULL_TEST") != nullptr) {
                auto bk = std::make_unique<TFHEpp::bfvboot::BootstrapKey<P>>(
                    TFHEpp::bfvboot::MakeBootstrapKey<P>(*key, true));
                auto noisy_slots = std::make_unique<TFHEpp::TRLWE<BP>>();
                TFHEpp::bfvboot::BootstrapNoisySlots<P>(*noisy_slots, *ct_a,
                                                        *bk);
                auto full_noisy_plain =
                    std::make_unique<std::array<uint64_t, BP::n>>();
                TFHEpp::trlweSlotDecrypt<BP>(*full_noisy_plain, *noisy_slots,
                                             *boot_key);
            }
        }
    }

    {
        using P = MLTestParam;
        TFHEpp::TRLWE<P> lhs{}, rhs{};
        lhs[1][0] = TFHEpp::bfvEncodeCoeff<P>(3);
        lhs[1][1] = TFHEpp::bfvEncodeCoeff<P>(2);
        rhs[1][0] = TFHEpp::bfvEncodeCoeff<P>(5);
        rhs[1][1] = TFHEpp::bfvEncodeCoeff<P>(7);

        TFHEpp::TRLWE3<P> mult{};
        TFHEpp::TRLWEMultWithoutRelinearizationFullDD<P>(mult, lhs, rhs);

        const uint64_t expected[P::n] = {15, 31, 14, 0, 0, 0, 0, 0};
        for (std::size_t i = 0; i < P::n; i++) {
            require(TFHEpp::bfvDecodeCoeff<P>(mult[1][i]) == expected[i],
                    "multi-limb FullDD body multiplication");
            require(TFHEpp::bfvDecodeCoeff<P>(mult[0][i]) == 0,
                    "multi-limb FullDD cross term zero");
            require(TFHEpp::bfvDecodeCoeff<P>(mult[2][i]) == 0,
                    "multi-limb FullDD square term zero");
        }
    }

    {
        using P = MLTestParam;
        TFHEpp::Polynomial<P> torus{}, got{};
        std::array<uint64_t, P::n> neg_one{};
        neg_one[0] = P::plain_modulus_u64 - 1;
        torus[0] = P::T{5};
        torus[1] = P::T{7};
        TFHEpp::PolyMulTorusByCenteredPlain<P>(got, torus, neg_one);
        require(got[0] == (P::T{0} - P::T{5}),
                "multi-limb centered constant multiply coefficient 0");
        require(got[1] == (P::T{0} - P::T{7}),
                "multi-limb centered constant multiply coefficient 1");
    }

    std::cout << "BFV multi-limb torus scaffold test" << std::endl;
    std::cout << "  lvl5 n=" << TFHEpp::lvl5param::n
              << " q_bits=" << std::numeric_limits<U>::digits
              << " p=" << TFHEpp::lvl5param::plain_modulus_u64 << std::endl;
    std::cout << "PASS" << std::endl;
    return 0;
}

#endif
