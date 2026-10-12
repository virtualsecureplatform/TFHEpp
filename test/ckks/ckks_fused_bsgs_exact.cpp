#include <ckks/ckks.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <string>

namespace {
struct P : TFHEpp::lvl5param {
    static constexpr std::uint32_t B̅gbit = 4, l̅ = T::digits / B̅gbit;
};
struct P8192 : P {
    static constexpr std::uint32_t nbit = 13, n = 1U << nbit;
};
struct P32768 : P {
    static constexpr std::uint32_t nbit = 15, n = 1U << nbit;
};
struct Basis576 : P32768 {
    using T = TFHEpp::MultiLimbUInt<9>;
    static constexpr std::uint32_t Bgbit = 18, B̅gbit = 8;
    static constexpr T Bg = T{1} << Bgbit;
    static constexpr std::uint32_t l̅ = T::digits / B̅gbit;
};
struct Basis640 : P32768 {
    using T = TFHEpp::MultiLimbUInt<10>;
    static constexpr std::uint32_t l̅ = T::digits / B̅gbit;
};
using Integer = boost::multiprecision::cpp_int;

template <class Torus>
Integer integer(const Torus &value)
{
    Integer result = 0;
    for (std::size_t i = Torus::limbs; i-- > 0;) {
        result <<= 64;
        result += value.limb[i];
    }
    return result;
}

template <class Torus>
Torus torus(Integer value)
{
    value &= (Integer{1} << Torus::digits) - 1;
    Torus result{};
    for (auto &word : result.limb) {
        word = (value & ((Integer{1} << 64) - 1)).convert_to<std::uint64_t>();
        value >>= 64;
    }
    return result;
}

template <unsigned Q, unsigned Drop, class Ring = P>
void check(unsigned mode)
{
    using P = Ring;
    using namespace TFHEpp;
    auto input = std::make_unique<CKKSCiphertext<P, Q, 48>>();
    std::mt19937_64 random(20262119 + Q);
    for (auto &poly : input->ct)
        for (auto &value : poly)
            for (auto &word : value.limb) word = random();

    auto group = std::make_unique<CKKSLinearTransformGiantStep<P, Q, Drop>>();
    const unsigned terms = mode == 2 ? 1 : 33;
    group->terms.resize(terms);
    const std::array<std::size_t, 4> positions{0, 1, P::n / 2 + 3, P::n - 1};
    const std::array<Integer, 4> values{
        (Integer{1} << 64) - 1, -((Integer{1} << 63) + 3),
        Integer{1} << 39, -1};
    std::array<Integer, 4> combined{};
    Integer dense_sum = 0;
    for (unsigned term = 0; term < terms; ++term) {
        auto &plain = group->terms[term].plain.poly;
        if (mode == 0) {
            for (std::size_t j = 0; j < positions.size(); ++j) {
                Integer value = values[j];
                if (term % 3 == 1) value = -value;
                plain[positions[j]] = torus<typename P::T>(value & ((Integer{1} << Q) - 1));
                combined[j] += value;
            }
        }
        else if (mode == 1) {
            const Integer value = (Integer{1} << 64) - 1;
            dense_sum += value;
            plain.fill(torus<typename P::T>(value));
        }
        else {
            plain[0] = typename P::T{1};
            combined[0] = 1;
        }
    }
    if (mode == 1) {
        typename P::T value{};
        for (unsigned shift = 0; shift < P::T::digits; shift += 16)
            value += typename P::T{32767} << shift;
        for (auto &poly : input->ct) poly.fill(value);
    }
    if (mode == 2) {
        for (auto &poly : input->ct) {
            poly[0] = typename P::T{1} << (Drop - 1);
            poly[1] = typename P::T{} - poly[0];
            poly[2] = typename P::T{1} << (Q - 1);
            poly[3] = (typename P::T{1} << (Q - 1)) - typename P::T{1};
        }
    }

    // A single unrotated group needs no keys. Repeating baby zero stresses
    // fused batches (including their final partial batch) without encryption.
    CKKSSparseGaloisKey<P, Q> input_keys;
    CKKSSparseGaloisKey<P, Q - Drop> output_keys;
    auto provider = [&](std::size_t) -> const auto & { return *group; };
    auto actual = std::make_unique<CKKSPlainMulResult<P, Q, 48, Drop>>();
    const auto begin = std::chrono::steady_clock::now();
    CKKSLinearTransformBSGSWithGroupProvider<P, Q, 48, Drop>(
        *actual, *input, 1, std::vector<bool>{true}, 1,
        provider, input_keys, output_keys);
    const auto finish = std::chrono::steady_clock::now();

    // These direct cpp_int convolutions and rounding do not use production
    // decomposition, FFT, modulus reduction, or rescale helpers.
    const Integer mask = (Integer{1} << Q) - 1;
    const Integer sign = Integer{1} << (Q - 1);
    const Integer modulus = Integer{1} << (Q - Drop);
    for (std::size_t component = 0; component <= P::k; ++component) {
        for (std::size_t i = 0; i < P::n; ++i) {
            Integer expected = 0;
            if (mode == 1) {
                expected = (Integer{2} * i + 2 - P::n) *
                    integer(input->ct[component][0]) * dense_sum;
            }
            else {
                for (std::size_t j = 0; j < positions.size(); ++j) {
                    const auto index = (i + P::n - positions[j]) % P::n;
                    const Integer contribution =
                        integer(input->ct[component][index]) * combined[j];
                    if (i < positions[j]) expected -= contribution;
                    else expected += contribution;
                }
            }
            expected &= mask;
            if (expected >= sign) expected -= Integer{1} << Q;
            const bool negative = expected < 0;
            if (negative) expected = -expected;
            expected += Integer{1} << (Drop - 1);
            expected >>= Drop;
            if (negative) expected = -expected;
            expected %= modulus;
            if (expected < 0) expected += modulus;
            if (integer(actual->ct[component][i]) != expected)
                throw std::runtime_error("fused group differs from independent integer convolution and rounding");
        }
    }
    std::cout << "exact fused BSGS passed N=" << P::n << " Q=" << Q
              << " storage_bits=" << P::T::digits
              << " drop=" << Drop << " mode=" << mode << " terms=" << terms
              << " coefficients=" << (P::k + 1) * P::n
              << " seconds=" << std::chrono::duration<double>(finish - begin).count()
              << '\n';
}
}  // namespace

template <class Ring>
void checkRing()
{
    for (unsigned mode = 0; mode < 3; ++mode) {
        check<81, 20, Ring>(mode);
        check<93, 20, Ring>(mode);
        check<413, 48, Ring>(mode);
        check<448, 48, Ring>(mode);
    }
}

int main(int argc, char **argv)
{
    try {
        if (argc == 1) checkRing<P>();
        else if (argc == 2 && std::string(argv[1]) == "--all-rings") {
            checkRing<P8192>();
            checkRing<P>();
            checkRing<P32768>();
        }
        else if (argc == 2 && std::string(argv[1]) == "--basis-rings") {
            for (unsigned mode = 0; mode < 3; ++mode) {
                check<129, 37, Basis576>(mode);
                check<518, 47, Basis576>(mode);
                check<576, 48, Basis576>(mode);
                check<129, 37, Basis640>(mode);
                check<536, 49, Basis640>(mode);
                check<640, 48, Basis640>(mode);
            }
        }
        else {
            std::cerr << "Usage: ckks_fused_bsgs_exact [--all-rings|--basis-rings]\n";
            return 2;
        }
    }
    catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
