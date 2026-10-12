#include <ckks/prepared_plain.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <memory>
#include <random>

namespace {
struct P : TFHEpp::lvl5param {
    static constexpr std::uint32_t B̅gbit = 4;
    static constexpr std::uint32_t l̅ = T::digits / B̅gbit;
};
using Integer = boost::multiprecision::cpp_int;

Integer integer(const P::T &value)
{
    Integer result = 0;
    for (std::size_t i = P::T::limbs; i-- > 0;) {
        result <<= 64;
        result += value.limb[i];
    }
    return result;
}

P::T torus(Integer value)
{
    const Integer mask = (Integer{1} << P::T::digits) - 1;
    value &= mask;
    P::T result{};
    for (auto &word : result.limb) {
        word = (value & ((Integer{1} << 64) - 1)).convert_to<std::uint64_t>();
        value >>= 64;
    }
    return result;
}

// A sparse plaintext permits direct arbitrary-precision negacyclic convolution
// at N=16384, without using production decomposition, FFT, or reduction helpers.
template <unsigned Q>
void check(const TFHEpp::Polynomial<P> &input, bool zero = false)
{
    auto plain = std::make_unique<TFHEpp::Polynomial<P>>();
    auto output = std::make_unique<TFHEpp::Polynomial<P>>();
    const Integer mask = (Integer{1} << Q) - 1;
    const std::array<std::size_t, 4> positions{0, 1, P::n / 2 + 3, P::n - 1};
    const std::array<Integer, 4> values{
        (Integer{1} << 64) - 1, -((Integer{1} << 63) + 3), Integer{1} << 39, -1};
    if (!zero)
        for (std::size_t j = 0; j < positions.size(); ++j)
            (*plain)[positions[j]] = torus(values[j] & mask);
    TFHEpp::ckks_detail::polyMulTorusBySmallSigned<P, Q>(*output, input, *plain);
    auto prepared = TFHEpp::makeCKKSPreparedSmallPlainFFT<P, Q>(*plain);
    auto cached_output = std::make_unique<TFHEpp::Polynomial<P>>();
    TFHEpp::CKKSPreparedSmallPlainMul<P, Q>(*cached_output, input, *prepared);
    for (std::size_t i = 0; i < P::n; ++i) {
        Integer expected = 0;
        if (!zero)
            for (std::size_t j = 0; j < positions.size(); ++j) {
                const bool wrap = i < positions[j];
                const auto index = (i + P::n - positions[j]) % P::n;
                const Integer contribution = integer(input[index]) * values[j];
                if (wrap) expected -= contribution;
                else expected += contribution;
            }
        if ((integer((*output)[i]) & mask) != (expected & mask) || (integer((*cached_output)[i]) & mask) != (expected & mask))
            throw std::runtime_error("signed plaintext product mismatch Q=" +
                                     std::to_string(Q) + " coefficient=" + std::to_string(i));
    }
    std::cout << "signed plaintext exact N=" << P::n << " Q=" << Q
              << " zero=" << zero << " coefficients=" << P::n << '\n';
}


template <unsigned Q>
void checkDense(const TFHEpp::Polynomial<P> &input, bool alternating)
{
    auto plain = std::make_unique<TFHEpp::Polynomial<P>>();
    auto output = std::make_unique<TFHEpp::Polynomial<P>>();
    const Integer magnitude = (Integer{1} << 64) - 1;
    Integer total = 0;
    for (std::size_t i = 0; i < P::n; ++i) {
        const bool negative = alternating && (i & 1);
        (*plain)[i] = torus(negative ? -magnitude : magnitude);
        total += negative ? -integer(input[i]) : integer(input[i]);
    }
    auto cache = TFHEpp::makeCKKSPreparedSmallPlainFFT<P, Q>(*plain);
    TFHEpp::CKKSPreparedSmallPlainMul<P, Q>(*output, input, *cache);
    const Integer mask = (Integer{1} << Q) - 1;
    Integer prefix = 0;
    for (std::size_t i = 0; i < P::n; ++i) {
        const bool negative = alternating && (i & 1);
        prefix += negative ? -integer(input[i]) : integer(input[i]);
        Integer expected = magnitude * (2 * prefix - total);
        if (negative) expected = -expected;
        if ((integer((*output)[i]) & mask) != (expected & mask))
            throw std::runtime_error("prepared dense signed product mismatch");
    }
    std::cout << "prepared dense signed exact N=" << P::n << " Q=" << Q
              << " alternating=" << alternating << " coefficients=" << P::n << '\n';
}
} // namespace

int main()
{
    try {
        auto input = std::make_unique<TFHEpp::Polynomial<P>>();
        std::mt19937_64 random(0x2119);
        for (auto &coefficient : *input)
            for (auto &word : coefficient.limb) word = random();
        // Exercise carries and sign boundaries as well as native high-bit data
        // above the active modulus. Every test checks every ring coefficient.
        (*input)[0] = P::T{};
        (*input)[1] = P::T{} - P::T{1};
        (*input)[2] = P::T{1} << 80;
        (*input)[3] = P::T{1} << 412;
        check<81>(*input);
        check<413>(*input);
        check<448>(*input);
        check<81>(*input, true);
        checkDense<413>(*input, false);
        checkDense<413>(*input, true);
        checkDense<448>(*input, true);
    }
    catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
