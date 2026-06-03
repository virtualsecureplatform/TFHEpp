#include <cstdint>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <tfhe++.hpp>

namespace {

using Clock = std::chrono::steady_clock;

template <class P, int validbit>
typename P::T decode_digits(const std::array<int, P::n> &digits)
{
    return TFHEpp::decodeHatEncoderInt2T<P, validbit, P::n>(digits);
}

template <class P, int validbit>
unsigned __int128 decode_digits_u128(const std::array<int, P::n> &digits)
{
    return TFHEpp::decodeHatEncoderInt2U128<P, validbit, P::n>(digits);
}

std::string to_string_u128(unsigned __int128 value)
{
    if (value == 0) return "0";
    std::string result;
    while (value != 0) {
        result.insert(result.begin(), static_cast<char>('0' + value % 10));
        value /= 10;
    }
    return result;
}

double elapsed_ms(const Clock::time_point start, const Clock::time_point end)
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

template <class P, int validbit>
bool check_digits(const std::array<int, P::n> &digits,
                  const typename P::T expected, const char *label)
{
    const auto decoded = decode_digits<P, validbit>(digits);
    if (decoded != expected) {
        std::cerr << label << " decoded as "
                  << static_cast<std::uint64_t>(decoded) << ", expected "
                  << static_cast<std::uint64_t>(expected) << "; first digits:";
        for (int i = 0; i < validbit; i++) std::cerr << ' ' << digits[i];
        std::cerr << std::endl;
        return false;
    }
    return true;
}

template <class P, int validbit>
bool check_digits_u128(const std::array<int, P::n> &digits,
                       const unsigned __int128 expected, const char *label)
{
    const auto decoded = decode_digits_u128<P, validbit>(digits);
    if (decoded != expected) {
        std::cerr << label << " decoded as " << to_string_u128(decoded)
                  << ", expected " << to_string_u128(expected)
                  << "; first digits:";
        for (int i = 0; i < validbit; i++) std::cerr << ' ' << digits[i];
        std::cerr << std::endl;
        return false;
    }
    return true;
}

template <class P, int validbit>
bool check_clpx_roundtrip(const typename P::T value, const TFHEpp::Key<P> &key)
{
    const auto poly = TFHEpp::EncodeHatEncoderP<P>(value);
    const auto ct = TFHEpp::clpxSymIntEncrypt<P>(poly, key);
    const auto digits = TFHEpp::clpxSymIntDecrypt<P>(ct, key);
    return check_digits<P, validbit>(digits, value, "CLPX roundtrip");
}

template <class P, int validbit>
bool check_clpx_mult(const typename P::T lhs, const typename P::T rhs,
                     const TFHEpp::Key<P> &key,
                     const TFHEpp::relinKeyFFT<P> &relinkeyfft)
{
    const auto lhs_poly = TFHEpp::EncodeHatEncoderP<P>(lhs);
    const auto rhs_poly = TFHEpp::EncodeHatEncoderP<P>(rhs);
    const auto lhs_ct = TFHEpp::clpxSymIntEncrypt<P>(lhs_poly, key);
    const auto rhs_ct = TFHEpp::clpxSymIntEncrypt<P>(rhs_poly, key);

    TFHEpp::TRLWE<P> product_ct;
    TFHEpp::CLPXMult<P>(product_ct, lhs_ct, rhs_ct, relinkeyfft);

    const auto product_digits = TFHEpp::clpxSymIntDecrypt<P>(product_ct, key);
    return check_digits<P, validbit>(product_digits, lhs * rhs, "CLPXMult");
}

template <class iksP, class bkP, class sskP, int num_multi, int shift, int w,
          int validbit>
bool switch_tlwes_to_clpx(typename bkP::targetP::T value,
                          TFHEpp::TRLWE<typename bkP::targetP> &switched,
                          const TFHEpp::SecretKey &sk, TFHEpp::EvalKey &ek,
                          const TFHEpp::AnnihilateKey<typename bkP::targetP> &ahk)
{
    std::vector<uint8_t> digits(validbit);
    for (int i = 0; i < validbit; i++)
        digits[i] = static_cast<uint8_t>((value >> i) & 1);
    std::vector<TFHEpp::TLWE<typename iksP::domainP>> tlwes;
    TFHEpp::bootsSymEncrypt<typename iksP::domainP>(tlwes, digits, sk);
    if (tlwes.size() != digits.size()) {
        std::cerr << "TLWE input size mismatch" << std::endl;
        return false;
    }

    TFHEpp::TLWES2CLPXIKS<iksP, bkP, sskP, num_multi, shift, w>(
        switched, tlwes, ahk, ek);
    return true;
}

template <int validbit>
bool check_switched_clpx_mult(
    const typename TFHEpp::SS2CLPXlvl2param::T lhs_value,
    const typename TFHEpp::SS2CLPXlvl2param::T rhs_value,
    const TFHEpp::SecretKey &sk, TFHEpp::EvalKey &ek,
    const TFHEpp::AnnihilateKey<TFHEpp::SS2CLPXlvl2param> &ahk,
    const TFHEpp::Key<TFHEpp::SS2CLPXlvl2param> &key,
    const TFHEpp::relinKeyFFT<TFHEpp::SS2CLPXlvl2param> &relinkeyfft)
{
    using iksP = TFHEpp::lvl1hparam;
    using bkP = TFHEpp::SS2CLPXlvlh2param;
    using sskP = TFHEpp::SS2CLPXlvl22param;
    using P = TFHEpp::SS2CLPXlvl2param;
    constexpr int lutnum = 4;
    constexpr int shiftnum = 5;
    constexpr int shift = shiftnum - 1;
    constexpr int w = 20;

    TFHEpp::TRLWE<P> lhs_ct;
    TFHEpp::TRLWE<P> rhs_ct;
    if (!switch_tlwes_to_clpx<iksP, bkP, sskP, lutnum, shift, w, validbit>(
            lhs_value, lhs_ct, sk, ek, ahk))
        return false;
    if (!switch_tlwes_to_clpx<iksP, bkP, sskP, lutnum, shift, w, validbit>(
            rhs_value, rhs_ct, sk, ek, ahk))
        return false;

    const auto lhs_digits = TFHEpp::clpxSymIntDecrypt<P>(lhs_ct, key);
    if (!check_digits<P, validbit>(lhs_digits, lhs_value, "TLWES2CLPXIKS lhs"))
        return false;
    const auto rhs_digits = TFHEpp::clpxSymIntDecrypt<P>(rhs_ct, key);
    if (!check_digits<P, validbit>(rhs_digits, rhs_value, "TLWES2CLPXIKS rhs"))
        return false;

    TFHEpp::TRLWE<P> product_ct;
    TFHEpp::CLPXMult<P>(product_ct, lhs_ct, rhs_ct, relinkeyfft);

    const auto product_digits = TFHEpp::clpxSymIntDecrypt<P>(product_ct, key);
    const auto product_expected = static_cast<unsigned __int128>(lhs_value) *
                                  static_cast<unsigned __int128>(rhs_value);
    return check_digits_u128<P, validbit * 2>(
        product_digits, product_expected, "TLWES2CLPXIKS + CLPXMult");
}

bool check_switched_clpx_mult_with_runtime(
    const typename TFHEpp::SS2CLPXlvl2param::T lhs_value,
    const typename TFHEpp::SS2CLPXlvl2param::T rhs_value,
    const TFHEpp::SecretKey &sk, TFHEpp::EvalKey &ek,
    const TFHEpp::AnnihilateKey<TFHEpp::SS2CLPXlvl2param> &ahk,
    const TFHEpp::Key<TFHEpp::SS2CLPXlvl2param> &key,
    const TFHEpp::relinKeyFFT<TFHEpp::SS2CLPXlvl2param> &relinkeyfft)
{
    using iksP = TFHEpp::lvl1hparam;
    using bkP = TFHEpp::SS2CLPXlvlh2param;
    using sskP = TFHEpp::SS2CLPXlvl22param;
    using P = TFHEpp::SS2CLPXlvl2param;
    using c2tIksP10 = TFHEpp::lvl1hparam;
    using c2tIksP21 = TFHEpp::lvl21param;
    using c2tBkP01 = TFHEpp::lvlh1param;
    using c2tBkP02 = TFHEpp::lvlh2param;
    using c2tIksP20 = TFHEpp::lvl2hparam;
    constexpr int validbit = 64;
    constexpr int productbit = validbit * 2;
    constexpr int lutnum = 4;
    constexpr int shiftnum = 5;
    constexpr int shift = shiftnum - 1;
    constexpr int w = 20;

    TFHEpp::TRLWE<P> lhs_ct;
    TFHEpp::TRLWE<P> rhs_ct;
    const auto tfhe2clpx_start = Clock::now();
    if (!switch_tlwes_to_clpx<iksP, bkP, sskP, lutnum, shift, w, validbit>(
            lhs_value, lhs_ct, sk, ek, ahk))
        return false;
    if (!switch_tlwes_to_clpx<iksP, bkP, sskP, lutnum, shift, w, validbit>(
            rhs_value, rhs_ct, sk, ek, ahk))
        return false;
    const auto tfhe2clpx_end = Clock::now();

    const auto lhs_digits = TFHEpp::clpxSymIntDecrypt<P>(lhs_ct, key);
    if (!check_digits<P, validbit>(lhs_digits, lhs_value, "TLWES2CLPXIKS lhs"))
        return false;
    const auto rhs_digits = TFHEpp::clpxSymIntDecrypt<P>(rhs_ct, key);
    if (!check_digits<P, validbit>(rhs_digits, rhs_value, "TLWES2CLPXIKS rhs"))
        return false;

    TFHEpp::TRLWE<P> product_ct;
    const auto mult_start = Clock::now();
    TFHEpp::CLPXMult<P>(product_ct, lhs_ct, rhs_ct, relinkeyfft);
    const auto mult_end = Clock::now();

    const auto product_digits = TFHEpp::clpxSymIntDecrypt<P>(product_ct, key);
    const auto product_expected = static_cast<unsigned __int128>(lhs_value) *
                                  static_cast<unsigned __int128>(rhs_value);
    if (!check_digits_u128<P, productbit>(
            product_digits, product_expected, "TLWES2CLPXIKS + CLPXMult"))
        return false;

    std::vector<TFHEpp::TLWE<typename c2tBkP01::targetP>> output(productbit);
    const auto clpx2tfhe_start = Clock::now();
    TFHEpp::CLPX2TLWESIKSanybit<c2tIksP10, c2tIksP21, c2tBkP01, c2tBkP02,
                                c2tIksP20, 4, 2>(output, product_ct, ek, sk);
    const auto clpx2tfhe_end = Clock::now();

    bool output_nonzero = false;
    for (const auto &tlwe : output)
        for (const auto value : tlwe)
            output_nonzero = output_nonzero || (value != 0);
    if (!output_nonzero) {
        std::cerr << "CLPX2TFHE produced all-zero TLWEs" << std::endl;
        return false;
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "TFHE2CLPX runtime (2 x 64-bit operands): "
              << elapsed_ms(tfhe2clpx_start, tfhe2clpx_end) << " ms"
              << std::endl;
    std::cout << "CLPX mult runtime (64-bit x 64-bit): "
              << elapsed_ms(mult_start, mult_end) << " ms" << std::endl;
    std::cout << "CLPX2TFHE runtime (128 output TLWEs): "
              << elapsed_ms(clpx2tfhe_start, clpx2tfhe_end) << " ms"
              << std::endl;
    return true;
}

}  // namespace

int main()
{
    using P = TFHEpp::lvl2param;
    constexpr int validbit = 8;

    TFHEpp::SecretKey sk;
    const auto &key = sk.key.get<P>();
    const auto relinkeyfft = TFHEpp::relinKeyFFTgen<P>(key);

    if (!check_clpx_roundtrip<P, validbit>(2, key)) return 1;
    if (!check_clpx_roundtrip<P, validbit>(3, key)) return 1;

    if (!check_clpx_mult<P, validbit>(2, 3, key, relinkeyfft)) return 1;
    if (!check_clpx_mult<P, validbit>(5, 6, key, relinkeyfft)) return 1;
    if (!check_clpx_mult<P, validbit>(7, 7, key, relinkeyfft)) return 1;

    using ssP = TFHEpp::SS2CLPXlvl2param;
    TFHEpp::EvalKey ss_ek;
    ss_ek.emplaceiksk<TFHEpp::lvl1hparam>(sk);
    ss_ek.emplaceiksk<TFHEpp::lvl21param>(sk);
    ss_ek.emplaceiksk<TFHEpp::lvl2hparam>(sk);
    ss_ek.emplacebkfft<TFHEpp::lvlh1param>(sk);
    ss_ek.emplacebkfft<TFHEpp::lvlh2param>(sk);
    const auto ss_key = sk.key.get<TFHEpp::lvl2param>();
    TFHEpp::AnnihilateKey<ssP> ss_ahk;
    TFHEpp::annihilatekeygen<ssP>(ss_ahk, ss_key);
    const auto ss_relinkeyfft = TFHEpp::relinKeyFFTgen<ssP>(ss_key);

    if (!check_switched_clpx_mult<8>(2, 3, sk, ss_ek, ss_ahk, ss_key,
                                     ss_relinkeyfft))
        return 1;
    if (!check_switched_clpx_mult<16>((1ULL << 15) - 1, 3, sk, ss_ek, ss_ahk,
                                      ss_key, ss_relinkeyfft))
        return 1;
    if (!check_switched_clpx_mult<32>((1ULL << 31) - 1, 3, sk, ss_ek, ss_ahk,
                                      ss_key, ss_relinkeyfft))
        return 1;
    if (!check_switched_clpx_mult_with_runtime(
            std::numeric_limits<std::uint64_t>::max(),
            std::numeric_limits<std::uint64_t>::max(), sk, ss_ek, ss_ahk,
            ss_key, ss_relinkeyfft))
        return 1;

    std::cout << "Passed" << std::endl;
}
