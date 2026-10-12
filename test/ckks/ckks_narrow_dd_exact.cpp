#include <ckks/compact_dd.hpp>
#include <ckks/narrow_bootstrapping.hpp>
#include <iostream>
#include <memory>
#include <random>
#include <cereal/archives/portable_binary.hpp>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <chrono>

namespace {
struct P : TFHEpp::lvl5param {
    static constexpr std::uint32_t Bgbit = 4, B̅gbit = 4;
    static constexpr std::uint32_t l̅ = T::digits / B̅gbit;
    static const inline double α = std::ldexp(1.0, -448);
    static constexpr double ckks_min_noise_stddev = 3.2;
};
constexpr TFHEpp::modular_ntt::PrimeModulus prime{4611686018427322369ULL, 7};

template <std::uint32_t Q>
P::T lowBits(P::T value)
{
    constexpr auto words = Q / 64, bits = Q % 64;
    if constexpr (bits) value.limb[words] &= (std::uint64_t{1} << bits) - 1;
    for (std::size_t i = words + (bits != 0); i < P::T::limbs; ++i) value.limb[i] = 0;
    return value;
}

template <std::uint32_t Q>
P::T magnitude(P::T value)
{
    value = lowBits<Q>(value);
    return ((value >> (Q - 1)).limb[0] & 1) ? lowBits<Q>(P::T{} - value) : value;
}

// Independent exact secret multiplication. Unsigned 16-bit limbs suffice
// modulo 2^Q; signed centering is needed only for each bounded NTT product.
struct SecretProducts {
    TFHEpp::modular_ntt::NegacyclicNTTPlan plan{P::n, prime};
    std::vector<TFHEpp::modular_ntt::ShoupMultiplier> secret, square_weights;
    std::vector<std::int64_t> square;

    explicit SecretProducts(const TFHEpp::Key<P> &key)
    {
        std::vector<std::uint64_t> spectrum(P::n);
        for (std::size_t n = 0; n < P::n; ++n) {
            const auto value = std::bit_cast<std::int64_t>(key[n].limb[0]);
            if (value < -1 || value > 1 || P::T::from_signed_i64(value) != key[n])
                throw std::runtime_error("nonternary exact oracle secret");
            spectrum[n] = value < 0 ? prime.value - 1 : value;
        }
        plan.forward(spectrum);
        secret.reserve(P::n);
        square_weights.reserve(P::n);
        for (std::size_t n = 0; n < P::n; ++n) {
            secret.emplace_back(spectrum[n], prime.value);
            spectrum[n] = secret.back().apply(spectrum[n], prime.value);
            square_weights.emplace_back(spectrum[n], prime.value);
        }
        plan.inverse(spectrum);
        square.resize(P::n);
        for (std::size_t n = 0; n < P::n; ++n) {
            square[n] = centered(spectrum[n]);
            if (std::abs(square[n]) > P::n) throw std::runtime_error("exact secret-square bound failed");
        }
    }

    static std::int64_t centered(std::uint64_t value)
    {
        return value <= prime.value / 2 ? static_cast<std::int64_t>(value)
            : -static_cast<std::int64_t>(prime.value - value);
    }

    template <std::uint32_t Q>
    auto multiply(const TFHEpp::Polynomial<P> &input,
                  const std::vector<TFHEpp::modular_ntt::ShoupMultiplier> &weights) const
    {
        // Both secret and secret-square coefficients have magnitude <= N.
        static_assert(static_cast<std::uint64_t>(P::n) * P::n * 65535 < prime.value / 2);
        auto output = std::make_unique<TFHEpp::Polynomial<P>>();
        std::vector<std::uint64_t> spectrum(P::n);
        for (std::uint32_t bit = 0; bit < Q; bit += 16) {
            const auto mask = (std::uint64_t{1} << std::min<std::uint32_t>(16, Q - bit)) - 1;
            for (std::size_t n = 0; n < P::n; ++n) spectrum[n] = (input[n] >> bit).limb[0] & mask;
            plan.forward(spectrum);
            for (std::size_t n = 0; n < P::n; ++n) spectrum[n] = weights[n].apply(spectrum[n], prime.value);
            plan.inverse(spectrum);
            for (std::size_t n = 0; n < P::n; ++n)
                (*output)[n] += P::T::from_signed_i64(centered(spectrum[n])) << bit;
        }
        for (auto &value : *output) value = lowBits<Q>(value);
        return output;
    }
};

template <std::uint32_t Q, std::uint32_t Bg>
std::uint64_t checkKeyPhase(const TFHEpp::CKKSCompactDDRelinKey<P, Q, Bg, 4> &key,
                           const SecretProducts &products)
{
    constexpr auto rows = (Q + Bg - 1) / Bg;
    constexpr auto secondary_rows = (Q + 3) / 4;
    auto row = std::make_unique<TFHEpp::TRLWE<P>>();
    std::uint64_t maximum = 0;
    for (std::size_t i = 0; i < rows; ++i) {
        for (auto &poly : *row) for (auto &value : poly) value = P::T{};
        for (std::size_t j = 0; j < secondary_rows; ++j)
            for (std::size_t c = 0; c < 2; ++c)
                for (std::size_t n = 0; n < P::n; ++n)
                    (*row)[c][n] += P::T::from_signed_i64(key.data[i][j][c][n]) << (4 * (secondary_rows - j - 1));
        const auto mask_secret = products.multiply<Q>((*row)[0], products.secret);
        for (std::size_t n = 0; n < P::n; ++n) {
            const auto message = P::T::from_signed_i64(products.square[n]) << (Bg * (rows - i - 1));
            const auto error = magnitude<Q>((*row)[1][n] - (*mask_secret)[n] - message);
            // A diagnostic 12-sigma envelope, not a cryptographic noise proof.
            if (error > P::T{39}) throw std::runtime_error("exact DD key row phase outside 12-sigma envelope");
            maximum = std::max(maximum, error.limb[0]);
        }
    }
    std::cout << "exact DD key phase passed Q=" << Q << " coefficients=" << rows * P::n
              << " max_noise=" << maximum << std::endl;
    return maximum;
}

// Independently decompose the low Q bits using balanced integer carries. No production
// decomposition, recomposition, or floating-point transform is used here.
template <std::uint32_t Q, std::uint32_t Bg>
auto reference(const TFHEpp::Polynomial<P> &input,
               const TFHEpp::CKKSCompactDDRelinKey<P, Q, Bg, 4> &key)
{
    using namespace TFHEpp;
    constexpr auto rows = (Q + Bg - 1) / Bg;
    constexpr auto secondary_rows = (Q + 3) / 4;
    constexpr auto half = std::uint64_t{1} << (Bg - 1);
    constexpr auto mask = 2 * half - 1;
    static_assert(static_cast<std::uint64_t>(rows) * P::n * half * 8 < prime.value / 2);
    const modular_ntt::NegacyclicNTTPlan plan(P::n, prime);
    std::vector<std::vector<std::uint64_t>> primary(rows, std::vector<std::uint64_t>(P::n));
    typename P::T offset{};
    for (std::uint32_t bit = 0; bit < Q; bit += Bg) offset += typename P::T{half} << bit;
    for (std::size_t n = 0; n < P::n; ++n) {
        auto centered = lowBits<Q>(input[n]);
        if ((centered >> (Q - 1)).limb[0] & 1) centered -= P::T{1} << Q;
        const auto adjusted = centered + offset;
        for (std::size_t i = 0; i < rows; ++i) {
            const auto digit = static_cast<std::int64_t>((adjusted >> (Bg * (rows - i - 1))).limb[0] & mask) - static_cast<std::int64_t>(half);
            primary[i][n] = digit < 0 ? prime.value - static_cast<std::uint64_t>(-digit) : digit;
        }
    }
    for (auto &row : primary) plan.forward(row);
    auto exact = std::make_unique<TRLWE<P>>();
    std::vector<std::uint64_t> spectrum(P::n), sums(P::n);
    for (std::size_t j = 0; j < secondary_rows; ++j)
        for (std::size_t c = 0; c < 2; ++c) {
            std::fill(sums.begin(), sums.end(), 0);
            for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t n = 0; n < P::n; ++n) {
                    const auto digit = static_cast<std::int64_t>(key.data[i][j][c][n]);
                    if (digit < -8 || digit > 7) throw std::runtime_error("invalid compact key digit");
                    spectrum[n] = digit < 0 ? prime.value - static_cast<std::uint64_t>(-digit) : digit;
                }
                plan.forward(spectrum);
                for (std::size_t n = 0; n < P::n; ++n) {
                    const modular_ntt::ShoupMultiplier weight(spectrum[n], prime.value);
                    sums[n] += weight.apply(primary[i][n], prime.value);
                    if (sums[n] >= prime.value) sums[n] -= prime.value;
                }
            }
            plan.inverse(sums);
            for (std::size_t n = 0; n < P::n; ++n) {
                const auto coefficient = sums[n] <= prime.value / 2
                    ? static_cast<std::int64_t>(sums[n]) : -static_cast<std::int64_t>(prime.value - sums[n]);
                (*exact)[c][n] += P::T::from_signed_i64(coefficient) << (4 * (secondary_rows - j - 1));
            }
        }
    // Independent power-of-two reduction, including garbage above active Q.
    for (auto &poly : *exact) for (auto &value : poly) value = lowBits<Q>(value);
    return exact;
}

template <std::uint32_t Q, std::uint32_t Bg>
void check(const TFHEpp::Key<P> &secret,
           const TFHEpp::CKKSCompactDDRelinKey<P, Q, Bg, 4> &key,
           std::mt19937_64 &random)
{
    using namespace TFHEpp;
    const SecretProducts products(secret);
    const auto maximum_key_noise = checkKeyPhase<Q>(key, products);
    auto input = std::make_unique<Polynomial<P>>();
    for (auto &value : *input) for (auto &word : value.limb) word = random();
    // Carry boundaries, signed wrap, negacyclic wrap, and native-width garbage
    // share the same full-ring sample as the uniformly random coefficients.
    (*input)[0] = P::T{};
    (*input)[1] = P::T{} - P::T{1};
    (*input)[2] = P::T{1} << (Q - 1);
    typename P::T offset{};
    for (std::uint32_t bit = 0; bit < Q; bit += Bg) offset += typename P::T{std::uint64_t{1} << (Bg - 1)} << bit;
    const auto threshold = Q % Bg == 0 ? (P::T{1} << Q) - offset : P::T{1} << (Q - 1);
    (*input)[3] = threshold - P::T{1};
    (*input)[4] = threshold;
    (*input)[5] = threshold + P::T{1};
    (*input)[P::n - 1] = P::T{1};
    std::cout << "exact full-ring DD reference starting Q=" << Q << std::endl;
    const auto exact = reference<Q>(*input, key);
    const auto fft = makeCKKSDDRelinKeyFFT(key);
    auto actual = std::make_unique<TRLWE<P>>();
    CKKSRelinKeySwitch<P, Q>(*actual, *input, *fft);
    for (std::size_t c = 0; c < 2; ++c)
        for (std::size_t n = 0; n < P::n; ++n)
            if ((*actual)[c][n] != (*exact)[c][n])
                throw std::runtime_error("full-ring FFT DD differs from independent exact NTT switch at Q=" +
                    std::to_string(Q) + " component=" + std::to_string(c) + " coefficient=" + std::to_string(n));
    std::cout << "exact full-ring DD passed N=" << P::n << " Q=" << Q
              << " primary_bits=" << Bg << " coefficients=" << 2 * P::n << std::endl;
    const auto mask_secret = products.multiply<Q>((*actual)[0], products.secret);
    const auto expected = products.multiply<Q>(*input, products.square_weights);
    const auto bound = static_cast<std::uint64_t>((Q + Bg - 1) / Bg) * P::n * (std::uint64_t{1} << (Bg - 1)) * maximum_key_noise;
    std::uint64_t maximum = 0;
    for (std::size_t n = 0; n < P::n; ++n) {
        const auto error = magnitude<Q>((*actual)[1][n] - (*mask_secret)[n] - (*expected)[n]);
        if (error > P::T{bound}) throw std::runtime_error("exact switched DD phase exceeds key-noise convolution bound");
        maximum = std::max(maximum, error.limb[0]);
    }
    std::cout << "exact DD switched phase passed Q=" << Q << " coefficients=" << P::n
              << " max_error=" << maximum << " bound=" << bound << std::endl;
}

template <std::uint32_t Q>
void run()
{
    std::mt19937_64 random(20262119 + Q);
    auto secret = std::make_unique<TFHEpp::Key<P>>();
    for (auto &value : *secret) value = static_cast<int>(random() % 3) - 1;
    const auto key = TFHEpp::makeCKKSCompactDDRelinKey<P, Q, 4, 4>(*secret, {P::α, 0});
    check<Q>(*secret, *key, random);
}

using RetainedChain = TFHEpp::CKKSCompactDDRelinKeyChain<P, 396, 48, 5, 4, 4>;

void timing(const TFHEpp::CKKSCompactDDRelinKey<P, 348, 4, 4> &source)
{
    using namespace TFHEpp;
    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();
    auto old_coefficients = makeCKKSCompactDDRelinKeyWithSecondaryBase<16>(source);
    const auto old_fft = makeCKKSDDRelinKeyFFT(*old_coefficients);
    old_coefficients.reset();
    const auto old_prepare = std::chrono::duration<double>(Clock::now()-start).count();
    const auto new_start = Clock::now();
    const auto projected = makeCKKSCompactDDRelinKeyWithPrimaryBase<8>(source);
    auto new_coefficients = makeCKKSCompactDDRelinKeyWithSecondaryBase<16>(*projected);
    const auto new_fft = makeCKKSDDRelinKeyFFT(*new_coefficients);
    new_coefficients.reset();
    const auto new_prepare = std::chrono::duration<double>(Clock::now()-new_start).count();
    std::mt19937_64 random(20262119);
    auto input = std::make_unique<Polynomial<P>>();
    for (auto &value : *input) for (auto &word : value.limb) word = random();
    (*input)[0] = P::T{};
    (*input)[1] = P::T{}-P::T{1};
    (*input)[2] = P::T{1} << 347;
    (*input)[P::n-1] = P::T{1};
    const auto old_exact = reference<348>(*input, source);
    const auto new_exact = reference<348>(*input, *projected);
    auto output = std::make_unique<TRLWE<P>>();
    std::cout << std::setprecision(17)
              << "{\"timing_profile\":\"retained-dd\",\"n\":" << P::n
              << ",\"log_q\":348,\"secondary_bits\":16,\"old_primary_rows\":87,\"new_primary_rows\":44"
              << ",\"old_cache_bytes\":" << sizeof(*old_fft) << ",\"new_cache_bytes\":" << sizeof(*new_fft)
              << ",\"old_cache_prepare_seconds\":" << old_prepare
              << ",\"new_cache_prepare_seconds\":" << new_prepare << ",\"fresh_keys_generated\":false}" << std::endl;
    // Warm both paths, then alternate the paired order to reduce warm-cache bias.
    CKKSRelinKeySwitch<P, 348>(*output, *input, *old_fft);
    if (*output != *old_exact) throw std::runtime_error("secondary16 baseline differs from exact DD oracle");
    CKKSRelinKeySwitch<P, 348>(*output, *input, *new_fft);
    if (*output != *new_exact) throw std::runtime_error("primary8/secondary16 differs from exact DD oracle");
    for (const bool wider : {false, true, true, false, false, true, true, false}) {
        const auto begin = Clock::now();
        constexpr std::size_t repetitions = 8;
        for (std::size_t i = 0; i < repetitions; i++) {
            if (wider) CKKSRelinKeySwitch<P, 348>(*output, *input, *new_fft);
            else CKKSRelinKeySwitch<P, 348>(*output, *input, *old_fft);
        }
        const auto elapsed = std::chrono::duration<double>(Clock::now()-begin).count();
        if (*output != *(wider ? new_exact : old_exact))
            throw std::runtime_error("timed DD result differs from independent exact oracle");
        std::cout << "{\"primary_bits\":" << (wider ? 8 : 4)
                  << ",\"repetitions\":" << repetitions << ",\"seconds\":" << elapsed
                  << ",\"exact_coefficients\":" << 2*P::n << "}" << std::endl;
    }
}

template <std::uint32_t AuxiliaryBits = 28>
void tree(const TFHEpp::Key<P> &secret, const RetainedChain &chain,
          bool primary8 = false, const std::string &public_path = {})
{
    using namespace TFHEpp;
    using Schedule = CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, true, AuxiliaryBits, true>;
    static_assert(std::is_same_v<typename Schedule::RelinKeys, RetainedChain>);
    std::cout << "preparing retained full-ring DD product-tree keys" << std::endl;
    std::ofstream public_file;
    std::unique_ptr<cereal::PortableBinaryOutputArchive> public_archive;
    if (!public_path.empty()) {
        if (std::filesystem::exists(public_path)) throw std::runtime_error("public tree output must be new");
        public_file.open(public_path, std::ios::binary);
        if (!public_file) throw std::runtime_error("cannot create public tree output");
        public_archive = std::make_unique<cereal::PortableBinaryOutputArchive>(public_file);
        const std::array<std::uint64_t, 5> header{0x5452454538303031ULL, P::n, 108, 48, 3};
        (*public_archive)(header);
    }
    auto evaluate = [&](const auto &prepared) {
    std::mt19937_64 random(20262119);
    std::uniform_real_distribution<double> cycles(-0.5, 0.5);
    const std::array<std::complex<double>, 3> weights{{{0.125, -0.25}, {1e-6, 2e-6}, {1, 0.5}}};
    for (std::size_t kind = 0; kind < weights.size(); ++kind) {
        std::vector<typename Schedule::Factor> leaves(Schedule::blocks + 1);
        auto slots = std::make_unique<CKKSSlotVector<P>>();
        auto expected = std::make_unique<CKKSSlotVector<P>>();
        auto active = [&](std::size_t i) {
            return kind == 2 || (kind == 0 ? i < Schedule::head_slots : i >= Schedule::head_slots);
        };
        for (std::size_t i = 0; i < Schedule::slots; ++i)
            (*expected)[i] = active(i) ? weights[kind] : std::complex<double>{};
        for (std::size_t leaf = 0; leaf < Schedule::blocks; ++leaf) {
            for (std::size_t i = 0; i < Schedule::slots; ++i) {
                (*slots)[i] = active(i) ? std::polar(1.0, 2 * std::acos(-1.0) * cycles(random))
                                      : std::complex<double>{};
                if (active(i)) (*expected)[i] *= (*slots)[i];
            }
            ckksSlotEncrypt<P>(leaves[leaf], *slots, secret, {P::α, 0});
        }
        for (std::size_t i = 0; i < Schedule::slots; ++i)
            (*slots)[i] = active(i) ? weights[kind] : std::complex<double>{};
        ckksSlotEncode<P, Schedule::factor_log_q, Schedule::log_delta>(leaves.back().ct[P::k], *slots);
        auto output = std::make_unique<typename Schedule::TreeResult>();
        narrow_detail::productTree<Schedule>(*output, std::move(leaves), chain, prepared.get());
        auto decoded = std::make_unique<CKKSSlotVector<P>>();
        ckksSlotDecrypt<P>(*decoded, *output, secret);
        double maximum = 0, squared = 0;
        for (std::size_t i = 0; i < Schedule::slots; ++i) {
            const auto error = std::abs((*decoded)[i] - (*expected)[i]);
            if (!std::isfinite(error)) throw std::runtime_error("nonfinite retained full-ring DD tree error");
            maximum = std::max(maximum, error);
            squared += error * error;
        }
        const auto rms = std::sqrt(squared / Schedule::slots);
        if (!std::isfinite(rms) || maximum > std::ldexp(1.0, -23))
            throw std::runtime_error("retained full-ring DD tree failed fixed 23-bit accuracy gate");
        if (public_archive) {
            const std::uint64_t case_index = kind;
            (*public_archive)(case_index, *expected, output->ct);
        }
        std::cout << std::setprecision(17) << "{\"tree_kind\":" << kind
                  << ",\"projected_primary8_dd\":" << (primary8 ? "true" : "false")
                  << ",\"fresh_keys_generated\":false"
                  << ",\"encrypted_leaves\":31,\"transparent_beta_leaves\":1,\"depth\":5"
                  << ",\"output_log_q\":108,\"output_log_delta\":48,\"slots_checked\":" << Schedule::slots
                  << ",\"maximum_error\":" << maximum << ",\"rms_error\":" << rms
                  << ",\"accuracy_passed\":true}" << std::endl;
    }
    };
    if (primary8) {
        auto projected = makeCKKSCompactDDRelinKeyChainWithPrimaryBase<8>(chain);
        const auto prepared = makeCKKSPreparedDDRelinKeyChainFFTWithSecondaryBase<
            Schedule::prepared_relin_bbarbit>(*projected);
        projected.reset();
        evaluate(prepared);
    }
    else {
        const auto prepared = CKKSNarrowMakePreparedRelinearization<Schedule>(chain);
        evaluate(prepared);
    }
    public_archive.reset();
    if (public_file.is_open()) {
        public_file.close();
        if (!public_file) throw std::runtime_error("cannot finish public tree output");
    }
}

template <std::uint32_t AuxiliaryBits>
void outputConjugation(const TFHEpp::Key<P> &secret,
                       const TFHEpp::CKKSAutoKey<P, 108 + AuxiliaryBits> &key)
{
    using namespace TFHEpp;
    constexpr auto Q = 108 + AuxiliaryBits;
    constexpr auto rows = CKKSAutoKey<P, Q>::key_switch_rows;
    static_assert((AuxiliaryBits == 8 || AuxiliaryBits == 16 || AuxiliaryBits == 28) && P::k == 1);
    const SecretProducts products(secret);
    std::uint64_t maximum_noise = 0;
    for (std::size_t row = 0; row < rows; ++row) {
        const auto mask_secret = products.multiply<Q>(key[0][row][0], products.secret);
        for (std::size_t n = 0; n < P::n; ++n) {
            // X -> X^-1 maps s[0] to itself and s[n] to -s[N-n].
            // Derive this directly, independently of production Automorphism.
            const auto conjugated = n == 0 ? secret[0] : P::T{} - secret[P::n - n];
            const auto message = conjugated << (Q - 4 - 4 * row);
            const auto error = magnitude<Q>(key[0][row][1][n] - (*mask_secret)[n] - message);
            if (error > P::T{39}) throw std::runtime_error("retained output conjugation key phase outside 12-sigma envelope");
            maximum_noise = std::max(maximum_noise, error.limb[0]);
        }
    }
    std::cout << "retained output conjugation key phase passed coefficients=" << rows * P::n
              << " max_noise=" << maximum_noise << std::endl;
    const auto prepared = makeCKKSPreparedAutoKeyFFT<P, Q>(key);
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto decoded = std::make_unique<CKKSSlotVector<P>>();
    auto input = std::make_unique<CKKSCiphertext<P, 108, 48>>();
    auto output = std::make_unique<CKKSCiphertext<P, 108, 48>>();
    std::mt19937_64 random(20262119);
    std::uniform_real_distribution<double> values(-0.5, 0.5);
    for (std::size_t kind = 0; kind < 3; ++kind) {
        for (std::size_t i = 0; i < slots->size(); ++i)
            (*slots)[i] = kind == 0 ? std::complex<double>{values(random), values(random)}
                       : kind == 1 ? std::polar(1.0, 2 * std::acos(-1.0) * values(random))
                                   : std::complex<double>{};
        if (kind == 2) {
            (*slots)[0] = {0.5, -0.5};
            (*slots)[slots->size() / 2] = {-0.5, 0.5};
            (*slots)[slots->size() - 1] = {0.25, 0.375};
        }
        ckksSlotEncrypt<P>(*input, *slots, secret, {P::α, 0});
        ckks_detail::applyWithAuxiliaryScale<P, 108, AuxiliaryBits>(output->ct, input->ct,
            [&](auto &out, const auto &in) { CKKSEvalAuto<P, Q>(out, in, 2 * P::n - 1, *prepared); });
        ckksSlotDecrypt<P>(*decoded, *output, secret);
        double maximum = 0, squared = 0;
        for (std::size_t i = 0; i < slots->size(); ++i) {
            const auto error = std::abs((*decoded)[i] - std::conj((*slots)[i]));
            if (!std::isfinite(error)) throw std::runtime_error("nonfinite retained output conjugation error");
            maximum = std::max(maximum, error);
            squared += error * error;
        }
        const auto rms = std::sqrt(squared / slots->size());
        if (!std::isfinite(rms) || maximum > std::ldexp(1.0, -23))
            throw std::runtime_error("retained output conjugation failed fixed 23-bit accuracy");
        std::cout << std::setprecision(17) << "{\"output_conjugation_kind\":" << kind
                  << ",\"slots_checked\":" << slots->size() << ",\"maximum_error\":" << maximum
                  << ",\"rms_error\":" << rms << ",\"accuracy_passed\":true}" << std::endl;
    }
}

template <std::uint32_t AuxiliaryBits>
void checkFactorGalois(const TFHEpp::CKKSAutoKey<P, 348 + AuxiliaryBits> &key,
                      const TFHEpp::Key<P> &secret, const SecretProducts &products,
                      std::uint64_t exponent, std::size_t index, bool check_evaluation = false)
{
    using namespace TFHEpp;
    constexpr auto Q = 348 + AuxiliaryBits;
    using Key = CKKSAutoKey<P, Q>;
    static_assert((AuxiliaryBits == 8 || AuxiliaryBits == 16 || AuxiliaryBits == 28) && P::k == 1);
    auto transformed = std::make_unique<Polynomial<P>>();
    // Independently substitute X -> X^d in the negacyclic ring. Odd d
    // permutes all indices; the high bit of n*d gives the negacyclic sign.
    for (std::size_t n = 0; n < P::n; ++n) {
        const auto position = (n * exponent) % (2 * P::n);
        (*transformed)[position % P::n] = position >= P::n ? P::T{} - secret[n] : secret[n];
    }
    std::uint64_t maximum = 0;
    for (std::size_t row = 0; row < Key::key_switch_rows; ++row) {
        const auto mask_secret = products.multiply<Q>(key[0][row][0], products.secret);
        for (std::size_t n = 0; n < P::n; ++n) {
            const auto message = (*transformed)[n] << (Q - 4 - 4 * row);
            const auto error = magnitude<Q>(key[0][row][1][n] - (*mask_secret)[n] - message);
            if (error > P::T{39}) throw std::runtime_error("factor Galois key phase outside 12-sigma envelope");
            maximum = std::max(maximum, error.limb[0]);
        }
    }
    std::cout << "{\"factor_key_index\":" << index << ",\"automorphism_exponent\":" << exponent
              << ",\"key_log_q\":" << Q
              << ",\"phase_coefficients_checked\":" << Key::key_switch_rows * P::n
              << ",\"maximum_noise\":" << maximum << ",\"phase_passed\":true}" << std::endl;
    if (!check_evaluation) return;
    const auto prepared = makeCKKSPreparedAutoKeyFFT<P, Q>(key);
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto expected = std::make_unique<CKKSSlotVector<P>>();
    auto decoded = std::make_unique<CKKSSlotVector<P>>();
    auto plaintext = std::make_unique<Polynomial<P>>();
    auto input = std::make_unique<CKKSCiphertext<P, 348, 48>>();
    auto output = std::make_unique<CKKSCiphertext<P, 348, 48>>();
    std::mt19937_64 random(20262119 + index);
    std::uniform_real_distribution<double> values(-0.5, 0.5);
    for (std::size_t kind = 0; kind < 3; ++kind) {
        for (std::size_t i = 0; i < slots->size(); ++i)
            (*slots)[i] = kind == 0 ? std::complex<double>{values(random), values(random)}
                       : kind == 1 ? std::polar(1.0, 2 * std::acos(-1.0) * values(random))
                                   : std::complex<double>{};
        if (kind == 2) {
            (*slots)[0] = {0.5, -0.5};
            (*slots)[slots->size() / 2] = {-0.5, 0.5};
            (*slots)[slots->size() - 1] = {0.25, 0.375};
        }
        ckksSlotEncode<P, 348, 48>(*plaintext, *slots);
        // Reference substitution acts directly on plaintext coefficients,
        // without production automorphism or key-switch routines.
        for (std::size_t n = 0; n < P::n; ++n) {
            const auto position = (n * exponent) % (2 * P::n);
            (*transformed)[position % P::n] = position >= P::n ? P::T{} - (*plaintext)[n] : (*plaintext)[n];
        }
        ckksSlotDecode<P, 348, 48>(*expected, *transformed);
        ckksSlotEncrypt<P>(*input, *slots, secret, {P::α, 0});
        ckks_detail::applyWithAuxiliaryScale<P, 348, AuxiliaryBits>(output->ct, input->ct,
            [&](auto &out, const auto &in) { CKKSEvalAuto<P, Q>(out, in, exponent, *prepared); });
        ckksSlotDecrypt<P>(*decoded, *output, secret);
        double max_error = 0, squared = 0;
        for (std::size_t i = 0; i < slots->size(); ++i) {
            const auto error = std::abs((*decoded)[i] - (*expected)[i]);
            if (!std::isfinite(error)) throw std::runtime_error("nonfinite factor Galois evaluation error");
            max_error = std::max(max_error, error);
            squared += error * error;
        }
        const auto rms = std::sqrt(squared / slots->size());
        if (!std::isfinite(rms) || max_error > std::ldexp(1.0, -23))
            throw std::runtime_error("factor Galois evaluation failed fixed 23-bit accuracy");
        std::cout << std::setprecision(17) << "{\"factor_key_index\":" << index
                  << ",\"evaluation_kind\":" << kind << ",\"slots_checked\":" << slots->size()
                  << ",\"maximum_error\":" << max_error << ",\"rms_error\":" << rms
                  << ",\"accuracy_passed\":true}" << std::endl;
    }
}

template <std::uint32_t AuxiliaryBits>
void factorGalois(std::ifstream &source, cereal::PortableBinaryInputArchive &archive,
                  const TFHEpp::Key<P> &secret, std::streamoff begin, std::streamoff end)
{
    using namespace TFHEpp;
    const SecretProducts products(secret);
    auto key = std::make_unique<CKKSAutoKey<P, 348 + AuxiliaryBits>>();
    std::uint64_t d = 5;
    source.seekg(begin);
    for (std::size_t index = 0; index <= P::nbit; ++index) {
        archive(*key);
        const auto exponent = index == P::nbit ? 2 * P::n - 1 : d;
        checkFactorGalois<AuxiliaryBits>(*key, secret, products, exponent, index);
        d = d * d % (2 * P::n);
    }
    if (source.tellg() != end) throw std::runtime_error("retained factor Galois payload length mismatch");
}

void freshAux16FactorGalois()
{
    using namespace TFHEpp;
    auto secret = std::make_unique<Key<P>>();
    std::mt19937_64 random(20262119);
    for (auto &value : *secret) value = static_cast<int>(random() % 3) - 1;
    const SecretProducts products(*secret);
    auto key = std::make_unique<CKKSAutoKey<P, 364>>();
    std::uint64_t d = 5;
    for (std::size_t index = 0; index <= P::nbit; ++index) {
        const auto exponent = index == P::nbit ? 2 * P::n - 1 : d;
        CKKSAutoKeyGen<P, 364>(*key, exponent, *secret, {P::α, 0});
        checkFactorGalois<16>(*key, *secret, products, exponent, index, true);
        d = d * d % (2 * P::n);
    }
}

template <std::uint32_t AuxiliaryBits = 28>
void checkpoint(const char *path, bool check_tree = false, bool check_output = false,
                bool check_factor = false, bool check_square = false, bool check_primary8 = false,
                const std::string &public_path = {}, bool check_timing = false)
{
    using namespace TFHEpp;
    using Chain = RetainedChain;
    std::ifstream source(path, std::ios::binary);
    if (!source) throw std::runtime_error("cannot open DD checkpoint");
    cereal::PortableBinaryInputArchive archive(source);
    std::string function;
    if constexpr (AuxiliaryBits == 8) {
        std::array<std::uint64_t, 10> metadata{};
        bool refresh_finished = false;
        archive(metadata, function, refresh_finished);
        const std::array<std::uint64_t, 10> expected{
            0x4155583850523031ULL, P::n, 85, 44, 385, 48, 37, 108, 8, 448};
        if (metadata != expected || !refresh_finished)
            throw std::runtime_error("unsupported or unfinished auxiliary8 probe state");
    }
    else {
        std::array<std::uint64_t, 9> metadata{};
        archive(metadata, function);
        const std::array<std::uint64_t, 9> expected{
            AuxiliaryBits == 16 ? 0x4e41525243503136ULL : 0x4e41525243503031ULL,
            P::n, 85, 44, 385, 48, 37, 108, 448};
        if (metadata != expected)
            throw std::runtime_error("unsupported full plain37 DD checkpoint profile");
    }
    if (function != "sigmoid" && function != "gelu" && function != "tanh")
        throw std::runtime_error("unsupported checkpoint function");
    auto secret = std::make_unique<Key<P>>();
    archive(*secret);
    const auto prefix_bytes = source.tellg();
    source.seekg(0, std::ios::end);
    const auto physical_file_bytes = source.tellg();
    auto file_bytes = physical_file_bytes;
    if constexpr (AuxiliaryBits == 8) {
        // The probe stores its retained-square key after the refresh key.
        // Establish that suffix boundary from the checked compact-key header
        // and payload, preserving this archive's original byte order.
        auto square = std::make_unique<CKKSCompactDDRelinKey<P, 108>>();
        const auto square_bytes = 7 * sizeof(std::uint64_t) + sizeof(square->data);
        if (file_bytes < prefix_bytes + static_cast<std::streamoff>(square_bytes))
            throw std::runtime_error("truncated auxiliary8 retained-square suffix");
        file_bytes -= static_cast<std::streamoff>(square_bytes);
        source.seekg(file_bytes);
        archive(*square);
        if (source.tellg() != physical_file_bytes || source.peek() != std::char_traits<char>::eof())
            throw std::runtime_error("invalid auxiliary8 retained-square termination");
        source.clear();
        std::cout << "validated auxiliary8 retained-square archive boundary at "
                  << file_bytes << std::endl;
        if (check_square) {
            std::mt19937_64 random(20262119);
            check<108>(*secret, *square, random);
            return;
        }
    }
    // Portable arrays have fixed payloads here: each compact key has seven
    // uint64 metadata words and int8 digit data. The narrow cache-mode uint32
    // is the sole field after the chain in the refresh-key archive. Seek using
    // the existing archive to preserve its original endianness conversion.
    auto chain = std::make_unique<Chain>();
    const auto chain_bytes = std::apply([](const auto &...key) {
        return ((7 * sizeof(std::uint64_t) + sizeof(key.data)) + ...);
    }, chain->keys);
    constexpr auto mode_bytes = sizeof(std::uint32_t);
    if (file_bytes < prefix_bytes + static_cast<std::streamoff>(chain_bytes + mode_bytes))
        throw std::runtime_error("truncated DD checkpoint suffix");
    const auto offset = file_bytes - static_cast<std::streamoff>(chain_bytes + mode_bytes);
    source.seekg(offset);
    if (!source) throw std::runtime_error("cannot seek DD checkpoint suffix");
    if (check_output || check_factor) {
        // Validate the adjoining compact-key header and terminal cache mode
        // before seeking back into the fixed-width raw automorphism payload.
        std::array<std::uint64_t, 7> anchor{};
        archive(anchor);
        const std::array<std::uint64_t, 7> expected_anchor{
            0x434b4b5343444431ULL, P::n, P::k, 348, 4, 4, 1};
        if (anchor != expected_anchor) throw std::runtime_error("invalid DD anchor after output conjugation key");
        source.seekg(file_bytes - static_cast<std::streamoff>(mode_bytes));
        std::uint32_t mode{};
        archive(mode);
        if (mode > 2 || source.tellg() != file_bytes ||
            (AuxiliaryBits != 8 && source.peek() != std::char_traits<char>::eof()))
            throw std::runtime_error("invalid output conjugation checkpoint termination");
        source.clear();
        chain.reset();
        if (check_factor) {
            const auto end = offset - static_cast<std::streamoff>(sizeof(CKKSAutoKey<P, 108 + AuxiliaryBits>::data));
            const auto bytes = sizeof(CKKSGaloisKey<P, 348 + AuxiliaryBits>);
            const auto begin = end - static_cast<std::streamoff>(bytes);
            if (begin < prefix_bytes) throw std::runtime_error("truncated factor Galois checkpoint");
            std::cout << "loaded retained factor Galois prefix_bytes=" << prefix_bytes
                      << " key_offset=" << begin << " key_bytes=" << bytes
                      << " cache_mode=" << mode << std::endl;
            factorGalois<AuxiliaryBits>(source, archive, *secret, begin, end);
            return;
        }
        auto key = std::make_unique<CKKSAutoKey<P, 108 + AuxiliaryBits>>();
        const auto key_offset = offset - static_cast<std::streamoff>(sizeof(key->data));
        if (key_offset < prefix_bytes) throw std::runtime_error("truncated output conjugation checkpoint");
        source.seekg(key_offset);
        archive(*key);
        if (source.tellg() != offset) throw std::runtime_error("output conjugation payload length mismatch");
        std::cout << "loaded retained output conjugation prefix_bytes=" << prefix_bytes
                  << " key_offset=" << key_offset << " key_bytes=" << sizeof(key->data)
                  << " cache_mode=" << mode << std::endl;
        outputConjugation<AuxiliaryBits>(*secret, *key);
        return;
    }
    archive(*chain);
    std::uint32_t mode{};
    archive(mode);
    if (mode > 2 || source.tellg() != file_bytes ||
        (AuxiliaryBits != 8 && source.peek() != std::char_traits<char>::eof()))
        throw std::runtime_error("invalid DD checkpoint suffix termination");
    std::cout << "loaded retained DD chain prefix_bytes=" << prefix_bytes << " suffix_offset=" << offset
              << " chain_bytes=" << chain_bytes << " cache_mode=" << mode << std::endl;
    if (check_tree) { tree<AuxiliaryBits>(*secret, *chain, check_primary8, public_path); return; }
    std::mt19937_64 random(20262119);
    if (check_timing) {
        auto original = std::make_unique<CKKSCompactDDRelinKey<P, 348, 4, 4>>(
            std::move(chain->template get<0>()));
        chain.reset();
        timing(*original);
        return;
    }
    if (check_primary8) {
        auto projected = makeCKKSCompactDDRelinKeyChainWithPrimaryBase<8>(*chain);
        [&]<std::size_t... I>(std::index_sequence<I...>) {
            auto compare = [&]<std::size_t J>() {
                const auto &old = chain->template get<J>();
                const auto &wide = projected->template get<J>();
                for (std::size_t row = 0; row < wide.primary_rows; row++) {
                    const auto exponent = (wide.primary_rows-row-1)*8;
                    if (wide.data[row] != old.data[old.primary_rows-1-exponent/4])
                        throw std::runtime_error("projected primary DD row differs from retained key");
                }
            };
            (compare.template operator()<I>(), ...);
        }(std::make_index_sequence<5>{});
        chain.reset(); // Bound memory while the independent full-ring oracle runs.
        [&]<std::size_t... I>(std::index_sequence<I...>) {
            (check<348 - 48 * I>(*secret, projected->template get<I>(), random), ...);
        }(std::make_index_sequence<5>{});
        return;
    }
    [&]<std::size_t... I>(std::index_sequence<I...>) {
        (check<348 - 48 * I>(*secret, chain->template get<I>(), random), ...);
    }(std::make_index_sequence<5>{});
}
} // namespace

int main(int argc, char **argv)
{
    if (argc == 1) {
        std::cout << "Use --full, --packed, --checkpoint CHECKPOINT.bin, --checkpoint-tree CHECKPOINT.bin, --checkpoint-output-conjugation CHECKPOINT.bin, or --checkpoint-factor-galois CHECKPOINT.bin for N16384 checks. Aux16 modes: --aux16-factor-galois, --aux16-output-conjugation and --aux16-checkpoint[-tree|-output-conjugation|-factor-galois]\n";
        return 0;
    }
    try {
        if (argc == 3 && std::string(argv[1]) == "--aux8-checkpoint-dd-timing") {
            checkpoint<8>(argv[2], false, false, false, false, false, "", true);
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--aux16-factor-galois") { freshAux16FactorGalois(); return 0; }
        if (argc == 2 && std::string(argv[1]) == "--aux16-output-conjugation") {
            auto secret = std::make_unique<TFHEpp::Key<P>>();
            std::mt19937_64 random(20262119);
            for (auto &value : *secret) value = static_cast<int>(random() % 3) - 1;
            auto key = std::make_unique<TFHEpp::CKKSAutoKey<P, 124>>();
            TFHEpp::CKKSAutoKeyGen<P, 124>(*key, 2 * P::n - 1, *secret, {P::α, 0});
            outputConjugation<16>(*secret, *key);
            return 0;
        }
        if (argc == 3 && std::string(argv[1]) == "--aux16-checkpoint") { checkpoint<16>(argv[2]); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux16-checkpoint-tree") { checkpoint<16>(argv[2], true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux16-checkpoint-output-conjugation") { checkpoint<16>(argv[2], false, true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux16-checkpoint-factor-galois") { checkpoint<16>(argv[2], false, false, true); return 0; }
        if (argc == 4 && std::string(argv[1]) == "--aux8-checkpoint-primary8-tree") { checkpoint<8>(argv[2], true, false, false, false, true, argv[3]); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux8-checkpoint-primary8") { checkpoint<8>(argv[2], false, false, false, false, true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux8-checkpoint") { checkpoint<8>(argv[2]); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux8-checkpoint-tree") { checkpoint<8>(argv[2], true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux8-checkpoint-output-conjugation") { checkpoint<8>(argv[2], false, true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux8-checkpoint-factor-galois") { checkpoint<8>(argv[2], false, false, true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--aux8-checkpoint-square") { checkpoint<8>(argv[2], false, false, false, true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--checkpoint") { checkpoint(argv[2]); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--checkpoint-tree") { checkpoint(argv[2], true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--checkpoint-output-conjugation") { checkpoint(argv[2], false, true); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--checkpoint-factor-galois") { checkpoint(argv[2], false, false, true); return 0; }
        if (argc == 2 && std::string(argv[1]) == "--packed") {
            run<296>();
            run<247>();
            run<198>();
            run<149>();
            run<100>();
            return 0;
        }
        if (argc != 2 || std::string(argv[1]) != "--full") return 2;
        run<348>();
        run<300>();
        run<252>();
        run<204>();
        run<156>();
        run<108>();
    }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
}
