#include <ckks/ckks.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>

namespace {
struct P : TFHEpp::lvl5param {
    static constexpr std::uint32_t B̅gbit = 4, l̅ = T::digits / B̅gbit;
};
using Wide = TFHEpp::MultiLimbUInt<16>;
using Spectra = std::vector<std::vector<std::uint64_t>>;
constexpr TFHEpp::modular_ntt::PrimeModulus prime{4611686018427322369ULL, 7};

template <std::uint32_t Drop>
void checkAccumulatorRescale()
{
    using boost::multiprecision::cpp_int;
    using Acc = TFHEpp::WideSignedLimbAccumulator<16>;
    std::mt19937_64 random(20262119 + Drop);
    for (std::uint32_t sample = 0; sample < 272; ++sample) {
        Acc value{};
        if (sample < 256) {
            // Leave a sign limb free so rounding cannot overflow this oracle.
            for (std::size_t i = 0; i + 1 < Acc::limbs; ++i)
                value.limb[i] = random();
        }
        else {
            value.limb[7] = sample % 4;
            if constexpr (Drop > 0) value.add_shifted_i64(1, Drop - 1);
            value.add_shifted_i64(static_cast<int>(sample % 3) - 1, 0);
        }
        if (sample % 2) TFHEpp::detail::twos_complement(value.limb);
        cpp_int original = 0;
        for (std::size_t i = Acc::limbs; i-- > 0;) {
            original <<= 64;
            original += value.limb[i];
        }
        if (value.is_negative()) original -= cpp_int{1} << Acc::digits;
        cpp_int magnitude = original < 0 ? -original : original;
        if constexpr (Drop > 0) magnitude += cpp_int{1} << (Drop - 1);
        magnitude >>= Drop;
        cpp_int expected = original < 0 ? -magnitude : magnitude;
        const cpp_int modulus = cpp_int{1} << P::T::digits;
        expected %= modulus;
        if (expected < 0) expected += modulus;
        const auto actual = TFHEpp::ckks_detail::rescaleLimbAccumulator<P, Drop>(value);
        cpp_int decoded = 0;
        for (std::size_t i = P::T::limbs; i-- > 0;) {
            decoded <<= 64;
            decoded += actual.limb[i];
        }
        if (decoded != expected)
            throw std::runtime_error("limb rescale differs from independent integer rounding");
    }
    std::cout << "exact accumulator rescale passed drop=" << Drop
              << " samples=272\n";
}

// Independent unsigned 16-bit radix, with a signed top chunk, represents the
// same integer as balanced 4-bit rows. Their asymmetric [-8,7] range moves the
// wrap threshold away from Q's sign bit; preserve this at signed rescale ties.
template <std::uint32_t Q>
Spectra prepare(const TFHEpp::Polynomial<P> &input,
                const TFHEpp::modular_ntt::NegacyclicNTTPlan &plan)
{
    constexpr auto rows = (Q + 15) / 16;
    typename P::T offset{};
    for (std::uint32_t shift = 0; shift < Q; shift += 4) offset += typename P::T{8} << shift;
    // A partial top nibble has spare bits, so centered active-Q integers fit
    // without the balanced-digit wrap that occurs at exact nibble boundaries.
    const auto threshold = Q % 4 == 0 ? (typename P::T{1} << Q) - offset
                                    : typename P::T{1} << (Q - 1);
    Spectra output(rows, std::vector<std::uint64_t>(P::n));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t i = 0; i < P::n; ++i) {
            const auto value = TFHEpp::ckks_detail::reduceToLevel<P, Q>(input[i]);
            auto digit = static_cast<std::int64_t>((value >> (16 * row)).limb[0] & 65535);
            if (row + 1 == rows) {
                constexpr auto top_bits = Q - 16 * (rows - 1);
                if (value >= threshold) digit -= std::int64_t{1} << top_bits;
            }
            output[row][i] = digit < 0 ? prime.value - static_cast<std::uint64_t>(-digit) : digit;
        }
        plan.forward(output[row]);
    }
    return output;
}

template <std::uint32_t LQ, std::uint32_t RQ>
auto reference(const TFHEpp::TRLWE<P> &lhs, const TFHEpp::TRLWE<P> &rhs)
{
    using namespace TFHEpp;
    constexpr auto lr = (LQ + 15) / 16, rr = (RQ + 15) / 16;
    // Two cross products, at most min(lr,rr) pairs per radix coefficient.
    // This strict bound makes centered single-prime reconstruction exact.
    static_assert(static_cast<unsigned __int128>(2) * std::min(lr, rr) * P::n * 65535 * 65535 < prime.value / 2);
    static_assert(LQ + RQ + P::nbit + 3 < Wide::digits);
    const modular_ntt::NegacyclicNTTPlan plan(P::n, prime);
    const std::array<Spectra, 2> left{prepare<LQ>(lhs[0], plan), prepare<LQ>(lhs[1], plan)};
    const std::array<Spectra, 2> right{prepare<RQ>(rhs[0], plan), prepare<RQ>(rhs[1], plan)};
    std::array<std::vector<std::vector<modular_ntt::ShoupMultiplier>>, 2> weights;
    for (std::size_t c = 0; c < 2; ++c) {
        weights[c].resize(rr);
        for (std::size_t row = 0; row < rr; ++row) {
            weights[c][row].reserve(P::n);
            for (const auto value : right[c][row]) weights[c][row].emplace_back(value, prime.value);
        }
    }
    auto exact = std::make_unique<std::array<std::array<Wide, P::n>, 3>>();
    std::vector<std::uint64_t> sums(P::n);
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t digit = 0; digit < lr + rr - 1; ++digit) {
            std::fill(sums.begin(), sums.end(), 0);
            auto add = [&](std::size_t a, std::size_t b, std::size_t row, std::size_t other) {
                for (std::size_t i = 0; i < P::n; ++i) {
                    sums[i] += weights[b][other][i].apply(left[a][row][i], prime.value);
                    if (sums[i] >= prime.value) sums[i] -= prime.value;
                }
            };
            for (std::size_t row = 0; row < lr; ++row) {
                if (digit < row || digit - row >= rr) continue;
                const auto other = digit - row;
                if (component == 0) { add(0, 1, row, other); add(1, 0, row, other); }
                else if (component == 1) add(1, 1, row, other);
                else add(0, 0, row, other);
            }
            plan.inverse(sums);
            for (std::size_t i = 0; i < P::n; ++i) {
                const auto signed_value = sums[i] <= prime.value / 2
                    ? static_cast<std::int64_t>(sums[i]) : -static_cast<std::int64_t>(prime.value - sums[i]);
                (*exact)[component][i] += Wide::from_signed_i64(signed_value) << (16 * digit);
            }
        }
    return exact;
}

template <std::uint32_t Q, std::uint32_t Drop>
typename P::T rounded(Wide value)
{
    const bool negative = (value.limb.back() >> 63) != 0;
    if (negative) value = Wide{} - value;
    if constexpr (Drop) value += Wide{1} << (Drop - 1);
    value >>= Drop;
    if (negative) value = Wide{} - value;
    typename P::T result{};
    std::copy_n(value.limb.begin(), P::T::limbs, result.limb.begin());
    return TFHEpp::ckks_detail::reduceToLevel<P, Q - Drop>(result);
}

template <std::uint32_t LQ, std::uint32_t RQ, std::uint32_t Drop>
void compare(const TFHEpp::TRLWE<P> &lhs, const TFHEpp::TRLWE<P> &rhs,
             const std::array<std::array<Wide, P::n>, 3> &exact)
{
    auto actual = std::make_unique<TFHEpp::TRLWE3<P>>();
    TFHEpp::CKKSTensorProductRescale<P, LQ, RQ, Drop>(*actual, lhs, rhs);
    for (std::size_t c = 0; c < 3; ++c)
        for (std::size_t i = 0; i < P::n; ++i)
            if ((*actual)[c][i] != rounded<std::min(LQ, RQ), Drop>(exact[c][i]))
                throw std::runtime_error("full-ring FFT tensor differs from independent exact NTT product");
    std::cout << "exact full-ring tensor passed N=" << P::n << " lhsQ=" << LQ
              << " rhsQ=" << RQ << " drop=" << Drop << " coefficients=" << 3 * P::n << '\n';
}

template <std::uint32_t LQ, std::uint32_t RQ, std::uint32_t Drop = 48>
void run(int boundary)
{
    auto lhs = std::make_unique<TFHEpp::TRLWE<P>>();
    auto rhs = std::make_unique<TFHEpp::TRLWE<P>>();
    if (!boundary) {
        std::mt19937_64 random(20262119 + LQ + RQ);
        for (auto *input : {lhs.get(), rhs.get()})
            for (auto &poly : *input) for (auto &value : poly) for (auto &limb : value.limb) limb = random();
    }
    else if (boundary == 2) {
        // Dense maximum positive balanced 16-bit digits stress FFT sums, not
        // just sparse edges or random coefficients. No encryption is needed.
        typename P::T value{};
        for (std::uint32_t shift = 0; shift < P::T::digits; shift += 16)
            value += typename P::T{32767} << shift;
        for (auto *input : {lhs.get(), rhs.get()})
            for (auto &poly : *input) poly.fill(value);
    }
    else {
        (*lhs)[0][0] = typename P::T{1} << (LQ - 1);
        (*rhs)[0][0] = typename P::T{1} << (RQ - 1);
        (*rhs)[1][0] = typename P::T{1};
        (*rhs)[0][1] = typename P::T{1};
        (*lhs)[1][1] = typename P::T{1} << (Drop - 1);
        (*lhs)[1][2] = typename P::T{} - (typename P::T{1} << (Drop - 1));
        (*lhs)[1][P::n - 1] = typename P::T{} - typename P::T{1};
        typename P::T offset{};
        for (std::uint32_t shift = 0; shift < LQ; shift += 4) offset += typename P::T{8} << shift;
        const auto threshold = LQ % 4 == 0 ? (typename P::T{1} << LQ) - offset
                                          : typename P::T{1} << (LQ - 1);
        (*lhs)[1][4] = threshold - typename P::T{1};
        (*lhs)[1][5] = threshold;
        (*lhs)[1][6] = threshold + typename P::T{1};
        (*rhs)[1][8] = typename P::T{1} << (Drop - 1);
    }
    const auto exact = reference<LQ, RQ>(*lhs, *rhs);
    compare<LQ, RQ, 0>(*lhs, *rhs, *exact);
    compare<LQ, RQ, Drop>(*lhs, *rhs, *exact);
}
} // namespace

int main(int argc, char **argv)
{
    if (argc == 1) { std::cout << "Use --full, --packed or --rescale for independent arithmetic checks\n"; return 0; }
    if (argc == 2 && std::string(argv[1]) == "--rescale") {
        checkAccumulatorRescale<0>();
        checkAccumulatorRescale<1>();
        checkAccumulatorRescale<48>();
        checkAccumulatorRescale<49>();
        checkAccumulatorRescale<63>();
        checkAccumulatorRescale<64>();
        checkAccumulatorRescale<65>();
        checkAccumulatorRescale<127>();
        checkAccumulatorRescale<128>();
        checkAccumulatorRescale<256>();
        checkAccumulatorRescale<447>();
        return 0;
    }
    if (argc == 2 && std::string(argv[1]) == "--packed") {
        run<345, 345, 49>(false);
        run<345, 296, 49>(false);
        run<296, 296, 49>(false);
        run<247, 247, 49>(false);
        run<198, 198, 49>(false);
        run<149, 149, 49>(false);
        run<345, 345, 49>(true);
        run<247, 247, 49>(true);
        run<198, 198, 49>(true);
        run<149, 149, 49>(true);
        return 0;
    }
    if (argc != 2 || std::string(argv[1]) != "--full") return 2;
    run<348, 348>(false);
    run<348, 300>(false);
    run<300, 300>(false);
    run<252, 252>(false);
    run<204, 204>(false);
    run<156, 156>(false);
    run<108, 108>(false);
    run<348, 348>(true);
    run<300, 300>(true);
    run<252, 252>(true);
    run<204, 204>(true);
    run<156, 156>(true);
    run<108, 108>(true);
    run<448, 448>(2);
}
