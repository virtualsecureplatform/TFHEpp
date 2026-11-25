#pragma once

#include "narrow_refresh.hpp"
#include "prepared_compact_linear.hpp"
#include <modular_ntt32.hpp>

#include <chrono>
#include <numbers>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Multi-limb torus for TUTTI schedules. Key noise is continuous sigma 3.6
// truncated toward zero, whose integer RMS (3.23) is not below the sigma 3.2
// assumed by the published modulus caps. PrimaryBits is the relinearization
// digit; SecondaryBits is the digit of every other dense-key switch.
template <std::uint32_t Nbit, std::size_t Limbs, std::uint32_t PrimaryBits = 12,
          std::uint32_t SecondaryBits = 8>
struct CKKSTuttiParam {
    static_assert(64 * Limbs % SecondaryBits == 0);
    static constexpr int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = Nbit, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t l = 1, l_a = 1, lₐ = 1;
    static constexpr std::uint32_t Bgbit = PrimaryBits, Bg_abit = PrimaryBits, Bgₐbit = PrimaryBits;
    using T = MultiLimbUInt<Limbs>;
    static constexpr T Bg = T{1} << Bgbit, Bg_a = Bg, Bgₐ = Bg;
    static constexpr ErrorDistribution errordist = ErrorDistribution::ModularGaussian;
    static const inline double α = std::ldexp(1.0, -static_cast<int>(64 * Limbs));
    static constexpr double ckks_min_noise_stddev = 3.6;
    static constexpr T μ = T{1} << (64 * Limbs - 3);
    static constexpr std::uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t B̅gbit = SecondaryBits, B̅gₐbit = SecondaryBits;
    static constexpr std::uint32_t l̅ = 64 * Limbs / SecondaryBits, l̅ₐ = l̅;
};

namespace tutti_detail {
// Same torus and ring with another key-switch digit. Ciphertext and key
// storage are shared, so only the gadget of the sparse switch changes.
template <class P, std::uint32_t DigitBits>
struct SwitchParam : P {
    static_assert(std::numeric_limits<typename P::T>::digits % DigitBits == 0);
    static constexpr std::uint32_t B̅gbit = DigitBits, B̅gₐbit = DigitBits;
    static constexpr std::uint32_t l̅ = std::numeric_limits<typename P::T>::digits / DigitBits;
    static constexpr std::uint32_t l̅ₐ = l̅;
};

// c for which c*g*sin(v/g)-v equioscillates on |v| <= bound, g = gap/(2 pi).
inline double minimaxSineScale(double gap, double bound)
{
    const double g = gap / (2 * std::numbers::pi), x = bound / g;
    if (!(x > 0 && x < std::numbers::pi / 2))
        throw std::invalid_argument("TUTTI message bound exceeds the sine range");
    double low = 1, high = x / std::sin(x);
    for (int i = 0; i < 200; ++i) {
        const double c = (low + high) / 2;
        const double interior = g * (std::sqrt(c * c - 1) - std::acos(1 / c));
        (interior < bound - c * g * std::sin(x) ? low : high) = c;
    }
    return high;
}
}  // namespace tutti_detail

// Full real-slot TUTTI (ePrint 2026/2848, Sections 4 and 5) over TFHEpp's
// power-of-two moduli. Config supplies:
//   Param                 multi-limb ring parameter
//   base_log_q            q0, the only modulus exposing the sparse key
//   input_log_delta       input scale; gap = 2^(base_log_q - input_log_delta)
//   s2c_plain_log_delta   level spent by the opening slot-to-coefficient map
//   log_delta             scale of the product tree, one per tree level
//   indicator_log_delta   scale of the encrypted block indicators
//   plain_log_delta       scale of the public root encodings
//   output_log_q, output_log_delta   refreshed ciphertext
//   blocks, subring_stride           block-binary key s'(X^stride)
//   auxiliary_log_scale   lift of every dense-key automorphism
//   sparse_digit_bits     digit of the switch into the sparse key
//   message_bound, minimax_sine      sine normalization
// AccMul sums at selector_log_q and divides once by
// 2^(plain + indicator - log_delta), the role of the paper's auxiliary prime.
template <class Config>
struct CKKSTuttiSchedule : CKKSRefreshLinearSchedule<
    typename Config::Param, Config::base_log_q, Config::input_log_delta,
    Config::s2c_plain_log_delta, 0, true, Config::auxiliary_log_scale> {
    using Param = typename Config::Param;
    using P = Param;
    using SparseParam = tutti_detail::SwitchParam<P, Config::sparse_digit_bits>;
    static constexpr std::size_t blocks = Config::blocks;
    static constexpr std::size_t stride = Config::subring_stride;
    static_assert(P::k == 1);
    static_assert(blocks > 1 && blocks <= P::n && std::has_single_bit(blocks));
    static constexpr std::size_t block_width = P::n / blocks;
    static_assert(std::has_single_bit(stride) && stride <= block_width);
    static constexpr std::size_t candidates = block_width / stride;
    static constexpr std::size_t indicator_count = (blocks - 1) * candidates;
    static constexpr std::uint32_t tree_depth = std::countr_zero(blocks);
    static constexpr std::uint32_t log_delta = Config::log_delta;
    static constexpr std::uint32_t indicator_log_delta = Config::indicator_log_delta;
    static constexpr std::uint32_t plain_log_delta = Config::plain_log_delta;
    static constexpr std::uint32_t output_log_q = Config::output_log_q;
    static constexpr std::uint32_t output_log_delta = Config::output_log_delta;
    static_assert(output_log_delta <= log_delta && indicator_log_delta >= log_delta);
    static_assert(plain_log_delta > 0 && plain_log_delta <= 50);
    static constexpr std::uint32_t output_drop = log_delta - output_log_delta;
    static constexpr std::uint32_t tree_log_q = output_log_q + output_drop;
    static constexpr std::uint32_t factor_log_q = tree_log_q + tree_depth * log_delta;
    static constexpr std::uint32_t accmul_drop =
        plain_log_delta + indicator_log_delta - log_delta;
    static constexpr std::uint32_t selector_log_q = factor_log_q + accmul_drop;
    static_assert(selector_log_q <= std::numeric_limits<typename P::T>::digits);
    static_assert(output_log_q > Config::base_log_q + Config::s2c_plain_log_delta,
                  "the refreshed modulus must exceed the refresh input modulus");
    static constexpr std::uint32_t output_key_log_q =
        output_log_q + Config::auxiliary_log_scale;
    static constexpr std::uint32_t relin_chain_start_log_q = factor_log_q + log_delta;
    static constexpr bool relinearize_before_rescale = true;
    using Selector = CKKSCiphertext<P, selector_log_q, indicator_log_delta>;
    using Factor = CKKSCiphertext<P, factor_log_q, log_delta>;
    using TreeResult = CKKSCiphertext<P, tree_log_q, log_delta>;
    using Output = CKKSCiphertext<P, output_log_q, output_log_delta>;
    using RelinKeys = CKKSCompactDDRelinKeyChain<P, relin_chain_start_log_q, log_delta, tree_depth>;
    using PreparedRelinKeys = CKKSPreparedDDRelinKeyChainFFT<
        P, relin_chain_start_log_q, log_delta, tree_depth>;
    // Retained S2C spectra: plaintext diagonals and both Galois key sets.
    static constexpr bool prepared_s2c = use_multilimb_digit_fft_v<P> &&
        P::B̅gbit + CKKSPreparedAutoSecondaryBits<P, CKKSTuttiSchedule::evaluation_input_log_q>() +
        P::nbit + std::bit_width(CKKSKeySwitchRowCountForLevel<
            P, CKKSTuttiSchedule::evaluation_input_log_q>()) + 3 <
            std::numeric_limits<double>::digits;
    static constexpr std::uint64_t prepared_s2c_budget_bytes = std::uint64_t{4} << 30;
    struct Unavailable {};
    using PreparedS2C = std::conditional_t<prepared_s2c,
        CKKSPreparedCompactRealLinearTransformFFT<P, CKKSTuttiSchedule::evaluation_input_log_q,
            CKKSTuttiSchedule::evaluation_input_log_delta, Config::s2c_plain_log_delta>,
        Unavailable>;
    using PreparedS2CInputKeys = std::conditional_t<prepared_s2c,
        CKKSRetainedPreparedGaloisKeyFFT<P, CKKSTuttiSchedule::evaluation_input_log_q>,
        Unavailable>;
    using PreparedS2COutputKeys = std::conditional_t<prepared_s2c,
        CKKSRetainedPreparedGaloisKeyFFT<P, CKKSTuttiSchedule::evaluation_coefficient_log_q>,
        Unavailable>;

    // Four-prime CRT range of one accumulator: `candidates` negacyclic
    // products of a root encoding with one unsigned selector digit.
    static constexpr std::size_t ntt_digit_bits = std::min<std::size_t>(60,
        116 - plain_log_delta - P::nbit - std::bit_width(candidates));
    static_assert(ntt_digit_bits >= 16, "root encodings exceed the selector CRT range");
    static constexpr std::size_t ntt_digits =
        (selector_log_q + ntt_digit_bits - 1) / ntt_digit_bits;

    // Amplitude of the product before the conjugation, and of one factor.
    static double sine_scale()
    {
        static const double value = Config::minimax_sine ?
            tutti_detail::minimaxSineScale(CKKSTuttiSchedule::gap(), Config::message_bound) : 1.0;
        return value;
    }
    static double amplitude() { return sine_scale() * CKKSTuttiSchedule::gap() / (4 * std::numbers::pi); }
    static double factor_amplitude() { return std::pow(amplitude(), 1.0 / blocks); }
};

namespace tutti_detail {
template <class Schedule>
struct IndicatorNTT {
    // Montgomery residues, row (component * digits + digit) of each prime.
    std::array<std::vector<std::uint32_t>, 4> spectra;
};

template <class P>
inline const auto &plans()
{
    static const std::array<modular_ntt32::NegacyclicNTTPlan, 4> value{
        modular_ntt32::NegacyclicNTTPlan(P::n, modular_ntt32::primes[0]),
        modular_ntt32::NegacyclicNTTPlan(P::n, modular_ntt32::primes[1]),
        modular_ntt32::NegacyclicNTTPlan(P::n, modular_ntt32::primes[2]),
        modular_ntt32::NegacyclicNTTPlan(P::n, modular_ntt32::primes[3])};
    return value;
}

// Position of each CKKS slot in the backend's negacyclic FFT, found by
// transforming X, so slot encodings need one half-size transform.
template <class P>
struct SlotLayout {
    // Offsets of the real and imaginary part of each slot's evaluation.
    std::array<std::uint32_t, P::n / 2> position{}, imaginary{};
    std::array<bool, P::n / 2> conjugate{};
    double scale = 1;  // backend normalization of a forward/inverse pair
    SlotLayout()
    {
        constexpr std::uint32_t half = P::n / 2;
        constexpr std::uint32_t invalid = 2 * P::n;
        const auto exponent = [](std::complex<double> z) {
            const double turns = std::arg(z) / std::numbers::pi * P::n;
            const auto value = static_cast<std::int64_t>(std::llround(turns));
            if (std::abs(turns - value) > 1e-3 || !(std::abs(std::abs(z) - 1) < 1e-6))
                return invalid;
            return static_cast<std::uint32_t>(((value % (2 * P::n)) + 2 * P::n) % (2 * P::n));
        };
        auto monomial = std::make_unique<std::array<std::int64_t, P::n>>();
        auto spectrum = std::make_unique<PolynomialInFD<P>>();
        monomial->fill(0);
        (*monomial)[1] = 1;
        TwistIFFTDigit<P>(*spectrum, *monomial);
        {
            auto copy = std::make_unique<PolynomialInFD<P>>(*spectrum);
            auto back = std::make_unique<std::array<std::uint64_t, P::n>>();
            TwistFFTDigitProduct64<P>(*back, *copy);
            scale = 1.0 / static_cast<double>(static_cast<std::int64_t>((*back)[1]));
        }
        // Backends store the half-size complex spectrum split or interleaved.
        std::vector<std::uint32_t> where(2 * P::n, P::n);
        std::uint32_t stride = 1, offset = half;
        for (std::uint32_t i = 0; i < half; ++i)
            if (exponent({(*spectrum)[i], (*spectrum)[i + half]}) == invalid) {
                stride = 2;
                offset = 1;
                break;
            }
        for (std::uint32_t i = 0; i < half; ++i) {
            const auto e = exponent({(*spectrum)[stride * i], (*spectrum)[stride * i + offset]});
            if (e == invalid) throw std::runtime_error("unknown negacyclic FFT layout");
            where[e] = stride * i;
        }
        auto polynomial = std::make_unique<Polynomial<P>>();
        polynomial->fill(typename P::T{0});
        (*polynomial)[1] = typename P::T{1} << 20;
        auto slots = std::make_unique<CKKSSlotVector<P>>();
        ckksSlotDecode<P, 40, 20>(*slots, *polynomial);
        for (std::uint32_t j = 0; j < half; ++j) {
            const auto e = exponent((*slots)[j]);
            if (e == invalid) throw std::runtime_error("unexpected CKKS evaluation point");
            conjugate[j] = where[e] == P::n;
            position[j] = conjugate[j] ? where[2 * P::n - e] : where[e];
            if (position[j] == P::n) throw std::runtime_error("unmatched CKKS slot");
            imaginary[j] = position[j] + offset;
        }
    }
};

template <class P>
inline const SlotLayout<P> &slotLayout()
{
    static const SlotLayout<P> value;
    return value;
}

// Integer polynomial whose slot values are `slots`, which already carry
// their scale. Rings without a digit FFT use the reference encoder.
template <class P>
inline void encodeSlots(std::array<std::int64_t, P::n> &encoded, const CKKSSlotVector<P> &slots)
{
    if constexpr (use_multilimb_digit_fft_v<P>) {
        constexpr std::uint32_t half = P::n / 2;
        const auto &layout = slotLayout<P>();
        auto spectrum = std::make_unique<PolynomialInFD<P>>();
        for (std::uint32_t i = 0; i < half; ++i) {
            const auto value = layout.conjugate[i] ? std::conj(slots[i]) : slots[i];
            (*spectrum)[layout.position[i]] = value.real() * layout.scale;
            (*spectrum)[layout.imaginary[i]] = value.imag() * layout.scale;
        }
        auto raw = std::make_unique<std::array<std::uint64_t, P::n>>();
        TwistFFTDigitProduct64<P>(*raw, *spectrum);
        for (std::uint32_t i = 0; i < P::n; ++i)
            encoded[i] = static_cast<std::int64_t>((*raw)[i]);
    }
    else {
        auto polynomial = std::make_unique<Polynomial<P>>();
        ckksSlotEncode<P, 62, 1>(*polynomial, slots);
        for (std::uint32_t i = 0; i < P::n; ++i)
            encoded[i] = static_cast<std::int64_t>(std::llround(
                ckks_detail::levelToLongDouble<P, 62>((*polynomial)[i]) / 2));
    }
}
}  // namespace tutti_detail

template <class Schedule>
struct CKKSTuttiKey {
    using P = typename Schedule::Param;
    typename Schedule::S2CPlan s2c_plan;
    CKKSGaloisKey<P, Schedule::evaluation_input_log_q> s2c_input_galois;
    CKKSGaloisKey<P, Schedule::evaluation_coefficient_log_q> s2c_output_galois;
    // Encrypted under the sparse key, at the base modulus only.
    CKKSSecretKeySwitchKey<typename Schedule::SparseParam, Schedule::coefficient_log_q> to_sparse;
    std::vector<typename Schedule::Selector> indicators;
    typename Schedule::RelinKeys relin;
    CKKSAutoKey<P, Schedule::output_key_log_q> conjugation;
    // Derived from the public members above; never serialized. The S2C
    // caches refer to the raw keys, so a prepared key must not be relocated.
    std::vector<tutti_detail::IndicatorNTT<Schedule>> indicator_ntt;
    std::unique_ptr<typename Schedule::PreparedRelinKeys> relin_fft;
    std::unique_ptr<typename Schedule::PreparedS2C> s2c_fft;
    std::unique_ptr<typename Schedule::PreparedS2CInputKeys> s2c_input_fft;
    std::unique_ptr<typename Schedule::PreparedS2COutputKeys> s2c_output_fft;
    CKKSTuttiKey() = default;
    CKKSTuttiKey(const CKKSTuttiKey &) = delete;
    CKKSTuttiKey &operator=(const CKKSTuttiKey &) = delete;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(s2c_plan.has_direct, s2c_plan.has_conjugate, s2c_plan.direct,
            s2c_plan.conjugate, s2c_input_galois, s2c_output_galois, to_sparse,
            indicators, relin, conjugation);
    }
};

struct CKKSTuttiTimings {
    double linear_seconds = 0, accmul_seconds = 0, tree_seconds = 0, finish_seconds = 0;
};

// Anchored block-binary key in the subring Z[X^stride]: one coefficient per
// block, on the stride grid, with s_0 = 1 public.
template <class Schedule>
inline void CKKSTuttiSparseKeyGen(Key<typename Schedule::Param> &secret)
{
    using T = typename Schedule::Param::T;
    secret.fill(T{0});
    secret[0] = T{1};
    std::uniform_int_distribution<std::size_t> choose(0, Schedule::candidates - 1);
    for (std::size_t block = 1; block < Schedule::blocks; ++block)
        secret[block * Schedule::block_width + choose(generator) * Schedule::stride] = T{1};
}

// Ternary key with `weight` nonzero coefficients and balanced signs.
// Rescaling error grows with the squared norm of this key, so a sparse
// dense key trades modulus budget for precision. Zero selects a uniform key.
template <class P>
inline void CKKSTuttiDenseKeyGen(Key<P> &secret, std::size_t weight = 0)
{
    using T = typename P::T;
    if (weight > P::n) throw std::invalid_argument("invalid TUTTI dense key weight");
    if (weight == 0) {
        std::uniform_int_distribution<int> ternary(-1, 1);
        for (auto &coefficient : secret) {
            const int value = ternary(generator);
            coefficient = value < 0 ? T{0} - T{1} : T{static_cast<std::uint64_t>(value)};
        }
        return;
    }
    secret.fill(T{0});
    for (std::size_t placed = 0; placed < weight;) {
        std::uniform_int_distribution<std::size_t> choose(0, P::n - 1);
        const auto position = choose(generator);
        if (secret[position] != T{0}) continue;
        secret[position] = placed++ % 2 ? T{0} - T{1} : T{1};
    }
}

// Rebuild the derived NTT and FFT caches from the public key material.
template <class Schedule>
inline void CKKSTuttiPrepareKey(CKKSTuttiKey<Schedule> &key)
{
    using P = typename Schedule::Param;
    constexpr std::size_t rows = (P::k + 1) * Schedule::ntt_digits;
    if (key.indicators.size() != Schedule::indicator_count)
        throw std::invalid_argument("incorrect TUTTI indicator count");
    const auto &plans = tutti_detail::plans<P>();
    key.indicator_ntt.resize(key.indicators.size());
#pragma omp parallel for schedule(static)
    for (std::size_t j = 0; j < key.indicators.size(); ++j)
        for (std::size_t prime = 0; prime < 4; ++prime) {
            auto &spectrum = key.indicator_ntt[j].spectra[prime];
            spectrum.resize(rows * P::n);
            for (std::size_t component = 0; component <= P::k; ++component)
                for (std::size_t digit = 0; digit < Schedule::ntt_digits; ++digit) {
                    auto row = std::span(spectrum).subspan(
                        (component * Schedule::ntt_digits + digit) * P::n, P::n);
                    for (std::size_t i = 0; i < P::n; ++i)
                        row[i] = plans[prime].arithmetic().encode(
                            narrow_detail::selectorDigit<P, Schedule::selector_log_q,
                                Schedule::ntt_digit_bits>(key.indicators[j].ct[component][i], digit));
                    plans[prime].forward(row);
                }
        }
    if constexpr (use_multilimb_digit_fft_v<P>)
        key.relin_fft = makeCKKSPreparedDDRelinKeyChainFFT(key.relin);
    if constexpr (Schedule::prepared_s2c) {
        key.s2c_fft = std::make_unique<typename Schedule::PreparedS2C>(
            key.s2c_plan, Schedule::prepared_s2c_budget_bytes);
        key.s2c_input_fft = std::make_unique<typename Schedule::PreparedS2CInputKeys>(
            key.s2c_input_galois);
        key.s2c_output_fft = std::make_unique<typename Schedule::PreparedS2COutputKeys>(
            key.s2c_output_galois);
    }
}

template <class Schedule>
inline void CKKSTuttiKeyGen(CKKSTuttiKey<Schedule> &key,
    const Key<typename Schedule::Param> &dense,
    const Key<typename Schedule::Param> &sparse,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    using T = typename P::T;
    for (std::size_t block = 0; block < Schedule::blocks; ++block) {
        std::size_t ones = 0;
        for (std::size_t offset = 0; offset < Schedule::block_width; ++offset) {
            const auto value = sparse[block * Schedule::block_width + offset];
            if (value == T{0}) continue;
            if (value != T{1} || offset % Schedule::stride || (block == 0 && offset))
                throw std::invalid_argument("invalid TUTTI block-binary key");
            ++ones;
        }
        if (ones != 1) throw std::invalid_argument("TUTTI requires one bit per block");
    }

    CKKSBuildPackedSlotToCoeffPlan<P, Schedule::evaluation_input_log_q,
        Schedule::evaluation_input_log_delta, Schedule::s2c_plain_log_delta>(
        key.s2c_plan, Schedule::bsgs_step);
    CKKSGaloisKeyGen<P, Schedule::evaluation_input_log_q>(key.s2c_input_galois, dense, noise);
    CKKSGaloisKeyGen<P, Schedule::evaluation_coefficient_log_q>(key.s2c_output_galois, dense, noise);
    CKKSSecretKeySwitchKeyGen<typename Schedule::SparseParam, Schedule::coefficient_log_q>(
        key.to_sparse, dense, sparse, noise);

    key.indicators.resize(Schedule::indicator_count);
    auto plain = std::make_unique<Polynomial<P>>();
    plain->fill(T{0});
    for (std::size_t row = 0; row < key.indicators.size(); ++row) {
        const auto position = Schedule::block_width + row * Schedule::stride;
        (*plain)[0] = sparse[position] << Schedule::indicator_log_delta;
        ckksEncryptPolynomial<P>(key.indicators[row], *plain, dense, noise);
    }
    CKKSDDRelinKeyChainGen<P, Schedule::relin_chain_start_log_q,
        Schedule::log_delta, Schedule::tree_depth>(key.relin, dense, noise);
    CKKSAutoKeyGen<P, Schedule::output_key_log_q>(key.conjugation, 2 * P::n - 1, dense, noise);
    CKKSTuttiPrepareKey(key);
}

// Opening slot-to-coefficient map under the dense key, lifted by the
// auxiliary scale, then the switch into the sparse key at the base modulus.
template <class Schedule>
inline void CKKSTuttiPrepareCoefficients(typename Schedule::Coefficients &coefficients,
    const typename Schedule::Input &input, const CKKSTuttiKey<Schedule> &key)
{
    using P = typename Schedule::Param;
    using SP = typename Schedule::SparseParam;
    constexpr auto aux = Schedule::auxiliary_log_scale;
    const auto transform = [&](auto &result, const auto &source) {
        if constexpr (Schedule::prepared_s2c) {
            if (key.s2c_fft) {
                CKKSRealLinearTransform(result, source, *key.s2c_fft,
                    *key.s2c_input_fft, *key.s2c_output_fft);
                return;
            }
        }
        narrow_detail::realTransformWithPreparedGalois<P, Schedule::evaluation_input_log_q,
            Schedule::evaluation_input_log_delta, Schedule::s2c_plain_log_delta>(
            result, source, key.s2c_plan, key.s2c_input_galois, key.s2c_output_galois);
    };
    auto dense = std::make_unique<CKKSCiphertext<SP, Schedule::coefficient_log_q,
        Schedule::input_log_delta>>();
    if constexpr (aux != 0) {
        auto lifted = std::make_unique<typename Schedule::EvaluationInput>();
        for (std::size_t c = 0; c <= P::k; ++c)
            for (std::size_t i = 0; i < P::n; ++i)
                lifted->ct[c][i] = ckks_detail::reduceToLevel<P, Schedule::input_log_q>(
                    input.ct[c][i]) << aux;
        auto raised = std::make_unique<typename Schedule::EvaluationCoefficients>();
        transform(*raised, *lifted);
        typename Schedule::Coefficients lowered;
        CKKSScaleDown<P, Schedule::evaluation_coefficient_log_q,
            Schedule::evaluation_input_log_delta, aux>(lowered, *raised);
        dense->ct = lowered.ct;
    }
    else {
        typename Schedule::Coefficients transformed;
        transform(transformed, input);
        dense->ct = transformed.ct;
    }
    auto sparse = std::make_unique<CKKSCiphertext<SP, Schedule::coefficient_log_q,
        Schedule::input_log_delta>>();
    CKKSSecretKeySwitch<SP>(*sparse, *dense, key.to_sparse);
    coefficients.ct = sparse->ct;
}

namespace tutti_detail {
// AccMul for one block: sum_j indicator_j * Ecd(roots shifted to candidate j),
// accumulated in the NTT domain and divided once by 2^accmul_drop.
template <class Schedule>
inline void accumulate(typename Schedule::Factor &factor, std::size_t block,
    std::span<const std::complex<double>> roots, const CKKSTuttiKey<Schedule> &key)
{
    using P = typename Schedule::Param;
    constexpr std::size_t half = P::n / 2, digits = Schedule::ntt_digits;
    constexpr std::size_t rows = (P::k + 1) * digits;
    const auto &plans = tutti_detail::plans<P>();
    const double scale = std::ldexp(Schedule::factor_amplitude(), Schedule::plain_log_delta);
    const auto limit = static_cast<std::int64_t>(std::ldexp(1.0, Schedule::plain_log_delta + 1));
    std::array<std::vector<std::uint32_t>, 4> sums;
    for (auto &sum : sums) sum.assign(rows * P::n, 0);
    auto encoded = std::make_unique<std::array<std::int64_t, P::n>>();
    auto table = std::make_unique<CKKSSlotVector<P>>();
    std::vector<std::uint32_t> residues(P::n);
    for (std::size_t row = 0; row < Schedule::candidates; ++row) {
        const auto position = block * Schedule::block_width + row * Schedule::stride;
        for (std::size_t i = 0; i < half; ++i) {
            // Coefficient i of X^position * a is a_{i-position}, negated on wrap.
            (*table)[i] = scale * (i >= position ? roots[i - position] :
                std::conj(roots[P::n + i - position]));
        }
        encodeSlots<P>(*encoded, *table);
        const auto &indicator = key.indicator_ntt[(block - 1) * Schedule::candidates + row];
        for (std::size_t prime = 0; prime < 4; ++prime) {
            const auto &arithmetic = plans[prime].arithmetic();
            const std::int64_t modulus = arithmetic.modulus();
            for (std::size_t i = 0; i < P::n; ++i) {
                const auto coefficient = (*encoded)[i];
                if (coefficient >= limit || coefficient <= -limit)
                    throw std::runtime_error("TUTTI root encoding exceeds its CRT bound");
                const auto reduced = coefficient % modulus;
                residues[i] = static_cast<std::uint32_t>(reduced < 0 ? reduced + modulus : reduced);
            }
            plans[prime].forward(residues);
            const std::uint32_t twice = 2 * arithmetic.modulus();
            const auto *selector = indicator.spectra[prime].data();
            auto *sum = sums[prime].data();
            for (std::size_t r = 0; r < rows; ++r, selector += P::n, sum += P::n)
                for (std::size_t i = 0; i < P::n; ++i) {
                    const auto value = sum[i] + arithmetic.multiplyLazy(residues[i], selector[i]);
                    sum[i] = value >= twice ? value - twice : value;
                }
        }
    }
    const modular_ntt32::FourPrimeCRT crt;
    auto accumulated = std::make_unique<Polynomial<P>>();
    for (std::size_t component = 0; component <= P::k; ++component) {
        accumulated->fill(typename P::T{0});
        for (std::size_t digit = 0; digit < digits; ++digit) {
            std::array<std::span<std::uint32_t>, 4> row;
            for (std::size_t prime = 0; prime < 4; ++prime) {
                row[prime] = std::span(sums[prime]).subspan((component * digits + digit) * P::n, P::n);
                for (auto &value : row[prime]) value = plans[prime].arithmetic().canonical(value);
                plans[prime].inverse(row[prime]);
            }
            for (std::size_t i = 0; i < P::n; ++i)
                (*accumulated)[i] += ckks_detail::signedToLevel<P, Schedule::selector_log_q>(
                    crt.reconstructSigned({row[0][i], row[1][i], row[2][i], row[3][i]}))
                    << (digit * Schedule::ntt_digit_bits);
        }
        for (std::size_t i = 0; i < P::n; ++i)
            factor.ct[component][i] = ckks_detail::roundedLevelRightShift<P,
                Schedule::selector_log_q, Schedule::accmul_drop>(
                ckks_detail::reduceToLevel<P, Schedule::selector_log_q>((*accumulated)[i]));
    }
}
}  // namespace tutti_detail

// The input is at the schedule's bottom level with real values in every slot.
// Evaluation takes only public ciphertexts and evaluation keys.
template <class Schedule>
inline void CKKSTuttiBootstrap(typename Schedule::Output &out,
    const typename Schedule::Input &input, const CKKSTuttiKey<Schedule> &key,
    CKKSTuttiTimings *timings = nullptr)
{
    using P = typename Schedule::Param;
    using Clock = std::chrono::steady_clock;
    constexpr auto base_q = Schedule::coefficient_log_q;
    constexpr double tau = 2 * std::numbers::pi;
    if (key.indicator_ntt.size() != Schedule::indicator_count)
        throw std::invalid_argument("TUTTI key caches are not prepared");
    auto record = [&, last = Clock::now()](double CKKSTuttiTimings::*field) mutable {
        const auto now = Clock::now();
        if (timings) timings->*field += std::chrono::duration<double>(now - last).count();
        last = now;
    };

    auto coeffs = std::make_unique<typename Schedule::Coefficients>();
    CKKSTuttiPrepareCoefficients<Schedule>(*coeffs, input, key);
    record(&CKKSTuttiTimings::linear_seconds);

    // One exponentiation per coefficient of the mask (Section 5.1).
    const auto phase = [](const auto &value) {
        return tau * static_cast<double>(std::ldexp(
            ckks_detail::levelToLongDouble<P, base_q>(value), -static_cast<int>(base_q)));
    };
    std::vector<std::complex<double>> roots(P::n);
    for (std::size_t i = 0; i < P::n; ++i)
        roots[i] = std::polar(1.0, -phase(coeffs->ct[0][i]));

    std::vector<typename Schedule::Factor> factors(Schedule::blocks);
    {
        // The anchored block contributes b - a at offset zero and carries
        // the -I that turns the final conjugate sum into a sine.
        auto table = std::make_unique<CKKSSlotVector<P>>();
        for (std::size_t i = 0; i < P::n / 2; ++i)
            (*table)[i] = std::complex<double>(0, -Schedule::factor_amplitude()) *
                std::polar(1.0, phase(coeffs->ct[1][i] - coeffs->ct[0][i]));
        ckksSlotEncode<P, Schedule::factor_log_q, Schedule::log_delta>(
            factors[0].ct[P::k], *table);
    }
    std::vector<std::exception_ptr> errors(Schedule::blocks);
#pragma omp parallel for schedule(dynamic)
    for (std::size_t block = 1; block < Schedule::blocks; ++block) {
        try { tutti_detail::accumulate<Schedule>(factors[block], block, roots, key); }
        catch (...) { errors[block] = std::current_exception(); }
    }
    for (const auto &error : errors)
        if (error) std::rethrow_exception(error);
    record(&CKKSTuttiTimings::accmul_seconds);

    auto tree = std::make_unique<typename Schedule::TreeResult>();
    narrow_detail::productTree<Schedule, 0, typename Schedule::PreparedRelinKeys, true>(
        *tree, std::move(factors), key.relin, key.relin_fft.get());
    record(&CKKSTuttiTimings::tree_seconds);

    if constexpr (Schedule::output_drop != 0)
        CKKSScaleDown<P, Schedule::tree_log_q, Schedule::log_delta, Schedule::output_drop>(out, *tree);
    else out = *tree;
    auto conjugate = std::make_unique<typename Schedule::Output>();
    narrow_detail::conjugateWithKey<Schedule, false>(conjugate->ct, out.ct, key.conjugation);
    CKKSAddInPlace<P>(out, *conjugate);
    record(&CKKSTuttiTimings::finish_seconds);
}

// N = 2^13, 64 blocks, subring stride 4, q0 = 2^24, gap 16: the shape of the
// paper's LL13 set, with scales chosen by Parameter-Selection
// (python/TUTTIparams.py --profile ll13). The dense key must have Hamming
// weight 512 (CKKSTuttiDenseKeyGen<P>(key, dense_weight)); indicators are the
// largest evaluation key, at 2^212. One multiplication at scale 2^20 returns
// the output to the refresh input modulus 2^38.
struct CKKSTuttiLL13Config {
    using Param = CKKSTuttiParam<13, 4>;
    static constexpr std::uint32_t base_log_q = 24, input_log_delta = 20;
    static constexpr std::uint32_t s2c_plain_log_delta = 14;
    static constexpr std::uint32_t log_delta = 22, indicator_log_delta = 24;
    static constexpr std::uint32_t plain_log_delta = 18;
    static constexpr std::uint32_t output_log_q = 58, output_log_delta = 20;
    static constexpr std::size_t blocks = 64, subring_stride = 4;
    static constexpr std::uint32_t auxiliary_log_scale = 20;
    static constexpr std::uint32_t sparse_digit_bits = 1;
    static constexpr std::size_t dense_weight = 512;
    static constexpr double message_bound = 1.0;
    static constexpr bool minimax_sine = true;
};

} // namespace TFHEpp
#endif
