#include <ckks/ckks.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <filesystem>

namespace {
bool plain37MetadataMatches(std::array<std::uint64_t, 9> metadata,
                            const std::array<std::uint64_t, 9> &expected)
{
    // The prefix and saved output have identical coefficient levels for the
    // aux16 candidate. No evaluation-key suffix is read by these modes.
    if (metadata[0] == 0x4e41525243503136ULL) metadata[0] = 0x4e41525243503031ULL;
    return metadata == expected;
}

template <std::uint32_t Bits> struct Param : TFHEpp::lvl5param {
    static constexpr std::uint32_t nbit = Bits, n = 1U << Bits;
};

// Inspect fresh expanded plain37 selectors without loading any evaluation-key
// suffix or a complete selector group. Secret multiplication is exact NTT,
// independent of the production FFT decryption path.
int inspectSelectorPhase(const char *checkpoint, const char *selector_file,
                         bool auxiliary8_probe = false)
{
    using namespace TFHEpp;
    using P = Param<14>;
    constexpr std::uint32_t Q = 385;
    static_assert(P::T::limbs == 7 && Q == 6 * 64 + 1);
    constexpr std::size_t blocks = 31, width = (P::n + blocks - 1) / blocks;
    constexpr std::size_t tail_width = width - 1, slots = P::n / 2;
    constexpr std::size_t tail_per_pack = slots / (2 * tail_width);
    constexpr modular_ntt::PrimeModulus prime{4611686018427322369ULL, 7};
    static_assert(static_cast<std::uint64_t>(P::n) * 65535 < prime.value / 2);
    const auto low = [](P::T value) {
        value.limb[6] &= 1;
        return value;
    };
    std::ifstream source(checkpoint, std::ios::binary);
    if (!source) throw std::runtime_error("cannot open checkpoint");
    cereal::PortableBinaryInputArchive archive(source);
    std::string function;
    if (auxiliary8_probe) {
        std::array<std::uint64_t, 10> metadata{};
        bool finished = false;
        archive(metadata, function, finished);
        const std::array<std::uint64_t, 10> expected{
            0x4155583850523031ULL, P::n, 85, 44, Q, 48, 37, 108, 8, 448};
        if (metadata != expected || !finished)
            throw std::runtime_error("selector phase inspector requires finished full auxiliary8 probe state");
    }
    else {
        std::array<std::uint64_t, 9> metadata{};
        archive(metadata, function);
        const std::array<std::uint64_t, 9> expected{
            0x4e41525243503031ULL, P::n, 85, 44, Q, 48, 37, 108, 448};
        if (!plain37MetadataMatches(metadata, expected))
            throw std::runtime_error("selector phase inspector requires full plain37 checkpoint");
    }
    if (function != "sigmoid")
        throw std::runtime_error("selector phase inspector requires full plain37 sigmoid checkpoint");
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    archive(*dense, *sparse);
    const auto prefix_bytes = source.tellg();
    source.close();
    std::vector<std::uint64_t> spectrum(P::n);
    for (std::size_t n = 0; n < P::n; ++n) {
        const auto digit = std::bit_cast<std::int64_t>((*dense)[n].limb[0]);
        if (digit < -1 || digit > 1 || P::T::from_signed_i64(digit) != (*dense)[n])
            throw std::runtime_error("nonternary checkpoint dense secret");
        spectrum[n] = digit < 0 ? prime.value - 1 : digit;
    }
    const modular_ntt::NegacyclicNTTPlan plan(P::n, prime);
    plan.forward(spectrum);
    std::vector<modular_ntt::ShoupMultiplier> weights;
    weights.reserve(P::n);
    for (const auto value : spectrum) weights.emplace_back(value, prime.value);
    for (std::size_t x = 0; x < blocks; ++x) {
        std::size_t count = 0;
        for (std::size_t y = 0; y < width && x * width + y < P::n; ++y) {
            const auto value = (*sparse)[x * width + y];
            if (value == P::T{1} || value == P::T{} - P::T{1}) ++count;
            else if (value != P::T{}) throw std::runtime_error("nonternary checkpoint sparse secret");
        }
        if (count != 1) throw std::runtime_error("checkpoint sparse block does not have one secret");
    }
    std::ifstream file(selector_file, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open selector group");
    cereal::PortableBinaryInputArchive selectors(file);
    std::array<std::uint64_t, 9> header{};
    selectors(header);
    const bool tail = header[7] != 0;
    const auto group = header[8];
    const std::array<std::uint64_t, 9> wanted_header{
        0x4e4152524f573031ULL, P::n, Q, 48, 37, width, 448, tail, group};
    if (header != wanted_header || group >= (tail ? (blocks + tail_per_pack - 1) / tail_per_pack : blocks) ||
        std::filesystem::file_size(selector_file) != 1 + sizeof(header) + width * sizeof(CKKSCiphertext<P, Q, 48>))
        throw std::runtime_error("unsupported selector group metadata or length");
    auto indicator = [&](std::size_t x, std::size_t y) -> std::complex<double> {
        if (x >= blocks || x * width + y >= P::n) return {};
        const auto value = (*sparse)[x * width + y];
        return value == P::T{1} ? std::complex<double>{1, 0}
            : value == P::T{} - P::T{1} ? std::complex<double>{0, 1} : std::complex<double>{};
    };
    std::uint64_t maximum = 0;
    std::size_t failed = 0;
    std::size_t failed_coefficients = 0;
    std::uint32_t maximum_error_bits = 0;
    std::vector<std::size_t> failed_records;
    for (std::size_t begin = 0; begin < width; begin += 8) {
        const auto count = std::min<std::size_t>(8, width - begin);
        std::vector<CKKSCiphertext<P, Q, 48>> chunk(count);
        selectors(cereal::binary_data(reinterpret_cast<std::uint64_t *>(chunk.data()), count * sizeof(chunk[0])));
        std::vector<Polynomial<P>> messages(count);
        for (std::size_t j = 0; j < count; ++j) {
            auto mask = std::make_unique<CKKSSlotVector<P>>();
            if (!tail) std::fill_n(mask->begin(), slots - tail_width, indicator(group, begin + j));
            else for (std::size_t ell = 0; ell < tail_per_pack; ++ell)
                std::fill_n(mask->begin() + (2 * ell + 1) * tail_width, tail_width,
                            indicator(group * tail_per_pack + ell, begin + j));
            ckksSlotEncode<P, Q, 48>(messages[j], *mask);
        }
        std::vector<std::uint64_t> maxima(count);
        std::vector<unsigned char> bad(count);
        std::vector<std::size_t> bad_coefficients(count);
        std::vector<std::uint32_t> error_bits(count);
#pragma omp parallel for schedule(static)
        for (std::size_t j = 0; j < count; ++j) {
            auto product = std::make_unique<Polynomial<P>>();
            std::vector<std::uint64_t> row(P::n);
            for (std::uint32_t bit = 0; bit < Q; bit += 16) {
                const auto mask = (std::uint64_t{1} << std::min<std::uint32_t>(16, Q - bit)) - 1;
                for (std::size_t n = 0; n < P::n; ++n) row[n] = (chunk[j].ct[0][n] >> bit).limb[0] & mask;
                plan.forward(row);
                for (std::size_t n = 0; n < P::n; ++n) row[n] = weights[n].apply(row[n], prime.value);
                plan.inverse(row);
                for (std::size_t n = 0; n < P::n; ++n) {
                    const auto signed_value = row[n] <= prime.value / 2 ? static_cast<std::int64_t>(row[n])
                        : -static_cast<std::int64_t>(prime.value - row[n]);
                    (*product)[n] += P::T::from_signed_i64(signed_value) << bit;
                }
            }
            for (std::size_t n = 0; n < P::n; ++n) {
                auto error = low(chunk[j].ct[1][n] - (*product)[n] - messages[j][n]);
                if (error.limb[6] & 1) error = low(P::T{} - error);
                if (error > P::T{39}) {
                    bad[j] = 1;
                    ++bad_coefficients[j];
                    for (std::size_t limb = 0; limb < P::T::limbs; ++limb)
                        if (error.limb[limb]) error_bits[j] = std::max(error_bits[j],
                            static_cast<std::uint32_t>(64 * limb + std::bit_width(error.limb[limb])));
                }
                else maxima[j] = std::max(maxima[j], error.limb[0]);
            }
        }
        for (std::size_t j = 0; j < count; ++j) {
            maximum = std::max(maximum, maxima[j]); failed += bad[j];
            failed_coefficients += bad_coefficients[j];
            maximum_error_bits = std::max(maximum_error_bits, error_bits[j]);
            if (bad[j]) failed_records.push_back(begin + j);
        }
    }
    if (file.peek() != std::char_traits<char>::eof()) throw std::runtime_error("trailing selector group data");
    std::cout << "{\"n\":" << P::n << ",\"group\":" << group << ",\"tail\":" << (tail ? "true" : "false")
              << ",\"selectors_checked\":" << width << ",\"phase_coefficients_checked\":" << width * P::n
              << ",\"checkpoint_prefix_bytes_read\":" << prefix_bytes << ",\"evaluation_keys_loaded\":false"
              << ",\"auxiliary8_probe_state\":" << (auxiliary8_probe ? "true" : "false")
              << ",\"failed_selectors\":" << failed << ",\"maximum_bounded_noise\":" << maximum
              << ",\"failed_phase_coefficients\":" << failed_coefficients
              << ",\"maximum_error_bits\":" << maximum_error_bits << ",\"failed_records\":[";
    for (std::size_t i = 0; i < failed_records.size(); ++i)
        std::cout << (i ? "," : "") << failed_records[i];
    std::cout << "],\"phase_accuracy_passed\":" << (failed ? "false" : "true") << "}\n";
    return failed ? 1 : 0;
}

// Full-ring phase statistics use exact digit NTT multiplication, independently
// of the production FFT secret product. Incorrect keys are explicit controls.
int inspectPhaseSpectrum(const char *checkpoint, const char *ciphertext)
{
    using namespace TFHEpp;
    using P = Param<14>;
    constexpr std::uint32_t Q = 108, Delta = 48;
    constexpr modular_ntt::PrimeModulus prime{4611686018427322369ULL, 7};
    static_assert(static_cast<std::uint64_t>(P::n) * 65535 < prime.value / 2);
    std::ifstream source(checkpoint, std::ios::binary);
    if (!source) throw std::runtime_error("cannot open phase checkpoint");
    cereal::PortableBinaryInputArchive archive(source);
    std::array<std::uint64_t, 9> metadata{};
    std::string function;
    archive(metadata, function);
    const std::array<std::uint64_t, 9> expected{
        0x4e41525243503031ULL, P::n, 85, 44, 385, Delta, 37, Q, 448};
    if (!plain37MetadataMatches(metadata, expected) || function != "sigmoid")
        throw std::runtime_error("phase spectrum requires full plain37 sigmoid checkpoint");
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto input = std::make_unique<CKKSCiphertext<P, 85, 44>>();
    archive(*dense, *sparse, *slots, *input);
    const auto prefix_bytes = source.tellg();
    source.close();
    auto output = std::make_unique<CKKSCiphertext<P, Q, Delta>>();
    std::ifstream result(ciphertext, std::ios::binary);
    if (!result) throw std::runtime_error("cannot open phase ciphertext");
    { cereal::PortableBinaryInputArchive saved(result); saved(*output); }
    if (result.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("trailing phase ciphertext data");
    auto changed = std::make_unique<Key<P>>(*dense);
    (*changed)[0] = (*changed)[0] == P::T{} ? P::T{1} : P::T{};
    const modular_ntt::NegacyclicNTTPlan plan(P::n, prime);
    const long double expected_rms = std::ldexp(std::sqrt(static_cast<long double>(P::n) / 12), Q - Delta);
    bool passed = true;
    std::cout << std::setprecision(17) << "{\"n\":" << P::n
              << ",\"log_q\":" << Q << ",\"log_delta\":" << Delta
              << ",\"checkpoint_prefix_bytes_read\":" << prefix_bytes
              << ",\"evaluation_keys_loaded\":false,\"uniform_phase_predicted_slot_rms\":"
              << expected_rms << ",\"cases\":[";
    const std::array<const Key<P> *, 3> keys{dense.get(), sparse.get(), changed.get()};
    const std::array<const char *, 3> names{"correct_dense", "incorrect_sparse", "one_dense_coefficient_changed"};
    for (std::size_t case_index = 0; case_index < keys.size(); ++case_index) {
        std::vector<std::uint64_t> spectrum(P::n), row(P::n);
        for (std::size_t i = 0; i < P::n; ++i) {
            const auto digit = std::bit_cast<std::int64_t>((*keys[case_index])[i].limb[0]);
            if (digit < -1 || digit > 1 || P::T::from_signed_i64(digit) != (*keys[case_index])[i])
                throw std::runtime_error("nonternary phase diagnostic secret");
            spectrum[i] = digit < 0 ? prime.value - 1 : digit;
        }
        plan.forward(spectrum);
        std::vector<modular_ntt::ShoupMultiplier> weights;
        for (const auto value : spectrum) weights.emplace_back(value, prime.value);
        auto product = std::make_unique<Polynomial<P>>();
        for (std::uint32_t bit = 0; bit < Q; bit += 16) {
            const auto mask = (std::uint64_t{1} << std::min<std::uint32_t>(16, Q - bit)) - 1;
            for (std::size_t i = 0; i < P::n; ++i) row[i] = (output->ct[0][i] >> bit).limb[0] & mask;
            plan.forward(row);
            for (std::size_t i = 0; i < P::n; ++i) row[i] = weights[i].apply(row[i], prime.value);
            plan.inverse(row);
            for (std::size_t i = 0; i < P::n; ++i) {
                const auto digit = row[i] <= prime.value / 2 ? static_cast<std::int64_t>(row[i])
                    : -static_cast<std::int64_t>(prime.value - row[i]);
                (*product)[i] += P::T::from_signed_i64(digit) << bit;
            }
        }
        auto phase = std::make_unique<Polynomial<P>>();
        std::vector<long double> normalized(P::n);
        long double sum = 0, squared_coefficients = 0;
        for (std::size_t i = 0; i < P::n; ++i) {
            (*phase)[i] = output->ct[1][i] - (*product)[i];
            const auto value = ckks_detail::levelToLongDouble<P, Q>((*phase)[i]);
            normalized[i] = std::ldexp(value, -static_cast<int>(Q));
            sum += normalized[i];
            squared_coefficients += normalized[i] * normalized[i];
        }
        const auto mean = sum / P::n;
        const auto variance = squared_coefficients / P::n - mean * mean;
        std::sort(normalized.begin(), normalized.end());
        long double ks = 0;
        for (std::size_t i = 0; i < P::n; ++i) {
            const auto cdf = normalized[i] + 0.5L;
            ks = std::max(ks, std::max(cdf - static_cast<long double>(i) / P::n,
                static_cast<long double>(i + 1) / P::n - cdf));
        }
        auto decoded = std::make_unique<CKKSSlotVector<P>>();
        ckksSlotDecode<P, Q, Delta>(*decoded, *phase);
        long double squared_slots = 0, squared_errors = 0;
        double maximum = 0;
        for (std::size_t i = 0; i < decoded->size(); ++i) {
            const auto target = std::complex<double>{1 / (1 + std::exp(-(*slots)[i].real())),
                1 / (1 + std::exp(-(*slots)[i].imag()))};
            const auto error = std::abs((*decoded)[i] - target);
            if (!std::isfinite(error) || !std::isfinite(std::norm((*decoded)[i])))
                throw std::runtime_error("nonfinite phase spectrum diagnostic");
            maximum = std::max(maximum, error);
            squared_errors += static_cast<long double>(error) * error;
            squared_slots += std::norm((*decoded)[i]);
        }
        const auto rms = std::sqrt(squared_errors / decoded->size());
        const auto slot_norm = std::sqrt(squared_slots / decoded->size());
        const auto coefficient_norm = std::ldexp(std::sqrt(squared_coefficients), Q - Delta);
        const auto parseval_error = std::abs(slot_norm - coefficient_norm) / std::max(slot_norm, coefficient_norm);
        if (!std::isfinite(parseval_error) || parseval_error > 1e-12L)
            throw std::runtime_error("phase Parseval check failed");
        const bool accurate = maximum < std::ldexp(1.0, -23);
        passed &= case_index == 0 ? accurate : (!accurate && rms / expected_rms > 0.9L && rms / expected_rms < 1.1L);
        if (case_index) std::cout << ',';
        std::cout << "{\"key\":\"" << names[case_index] << "\",\"normalized_phase_mean\":" << mean
                  << ",\"normalized_phase_variance\":" << variance << ",\"uniform_cdf_max_deviation\":" << ks
                  << ",\"slot_rms_error\":" << rms << ",\"maximum_error\":" << maximum
                  << ",\"rms_to_uniform_prediction\":" << rms / expected_rms
                  << ",\"parseval_relative_error\":" << parseval_error
                  << ",\"accuracy_passed\":" << (accurate ? "true" : "false") << '}';
    }
    std::cout << "],\"control_checks_passed\":" << (passed ? "true" : "false") << "}\n";
    return passed ? 0 : 1;
}

template <std::uint32_t Bits>
int inspect(const char *checkpoint, const char *ciphertext)
{
    using namespace TFHEpp;
    using P = Param<Bits>;
    std::ifstream source(checkpoint, std::ios::binary);
    if (!source) throw std::runtime_error("cannot open checkpoint");
    cereal::PortableBinaryInputArchive archive(source);
    std::array<std::uint64_t, 9> metadata{};
    std::string function;
    archive(metadata, function);
    const std::array<std::uint64_t, 9> expected{
        0x4e41525243503031ULL, P::n, 85, 44, 385, 48, 37, 108, 448};
    if (!plain37MetadataMatches(metadata, expected) ||
        (function != "sigmoid" && function != "gelu" && function != "tanh"))
        throw std::runtime_error("unsupported checkpoint profile/function");
    auto dense = std::make_unique<Key<P>>();
    auto sparse = std::make_unique<Key<P>>();
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto input = std::make_unique<CKKSCiphertext<P, 85, 44>>();
    archive(*dense, *sparse, *slots, *input);
    // The evaluation-key suffix is intentionally not deserialized. This is a
    // prefix inspector, not a validator of the entire checkpoint archive.
    const auto prefix_bytes = source.tellg();
    source.close();
    auto refreshed = std::make_unique<CKKSCiphertext<P, 108, 48>>();
    std::ifstream result(ciphertext, std::ios::binary);
    if (!result) throw std::runtime_error("cannot open refresh ciphertext");
    { cereal::PortableBinaryInputArchive output(result); output(*refreshed); }
    if (result.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("refresh ciphertext has trailing data");
    auto fresh = std::make_unique<CKKSSlotVector<P>>();
    auto decoded = std::make_unique<CKKSSlotVector<P>>();
    ckksSlotDecrypt<P>(*fresh, *input, *dense);
    ckksSlotDecrypt<P>(*decoded, *refreshed, *dense);
    const auto target = [&](double x) {
        if (function == "sigmoid") return 1 / (1 + std::exp(-x));
        if (function == "gelu") return 0.5 * x * (1 + std::erf(x / std::sqrt(2.0)));
        return std::tanh(x);
    };
    double fresh_error = 0, maximum = 0, squared = 0;
    for (std::size_t i = 0; i < slots->size(); ++i) {
        for (const auto z : {(*slots)[i], (*fresh)[i], (*decoded)[i]})
            if (!std::isfinite(z.real()) || !std::isfinite(z.imag()))
                throw std::runtime_error("nonfinite saved input/decrypted slot");
        const std::complex<double> wanted{target((*slots)[i].real()), target((*slots)[i].imag())};
        const auto error = std::abs((*decoded)[i] - wanted);
        const auto initial = std::abs((*fresh)[i] - (*slots)[i]);
        if (!std::isfinite(error) || !std::isfinite(initial))
            throw std::runtime_error("nonfinite checkpoint error");
        maximum = std::max(maximum, error);
        fresh_error = std::max(fresh_error, initial);
        squared += error * error;
    }
    const auto rms = std::sqrt(squared / slots->size());
    if (!std::isfinite(rms)) throw std::runtime_error("nonfinite checkpoint RMS error");
    const bool passed = maximum <= std::ldexp(1.0, -23);
    std::cout << std::setprecision(17)
              << "{\"n\":" << P::n << ",\"function\":\"" << function
              << "\",\"checkpoint_prefix_bytes_read\":" << prefix_bytes
              << ",\"evaluation_keys_loaded\":false,\"fresh_input_max_error\":" << fresh_error
              << ",\"refresh_max_error\":" << maximum << ",\"refresh_rms_error\":" << rms
              << ",\"target_tolerance\":" << std::ldexp(1.0, -23)
              << ",\"refresh_accuracy_passed\":" << (passed ? "true" : "false") << "}\n";
    return passed ? 0 : 1;
}
} // namespace

int main(int argc, char **argv)
{
    if (argc == 1) {
        std::cout << "Usage: ckks_narrow_checkpoint_inspect checkpoint.bin refresh.ciphertext\n";
        std::cout << "       ckks_narrow_checkpoint_inspect --phase-spectrum checkpoint.bin refresh.ciphertext\n";
        std::cout << "       ckks_narrow_checkpoint_inspect --selector-phase checkpoint.bin expanded-selector-group.bin\n";
        std::cout << "       ckks_narrow_checkpoint_inspect --aux8-selector-phase probe-state.bin expanded-selector-group.bin\n";
        return 0;
    }
    if (argc != 3 && argc != 4) return 2;
    try {
        if (argc == 4) {
            if (std::string(argv[1]) == "--phase-spectrum") return inspectPhaseSpectrum(argv[2], argv[3]);
            if (std::string(argv[1]) == "--aux8-selector-phase") return inspectSelectorPhase(argv[2], argv[3], true);
            if (std::string(argv[1]) != "--selector-phase") return 2;
            return inspectSelectorPhase(argv[2], argv[3]);
        }
        std::ifstream file(argv[1], std::ios::binary);
        if (!file) throw std::runtime_error("cannot open checkpoint");
        cereal::PortableBinaryInputArchive archive(file);
        std::array<std::uint64_t, 9> metadata{};
        archive(metadata);
        if (metadata[1] == 512) return inspect<9>(argv[1], argv[2]);
        if (metadata[1] == 16384) return inspect<14>(argv[1], argv[2]);
        throw std::runtime_error("unsupported checkpoint ring degree");
    }
    catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
