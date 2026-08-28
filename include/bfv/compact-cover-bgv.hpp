#pragma once

#include <algorithm>
#include <array>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "modular_ntt.hpp"
#ifdef USE_BLAKE3
#include "BLAKE3PRNG.hpp"
#endif

namespace TFHEpp::compact_cover_bgv {

inline constexpr std::size_t degree = 65536;
inline constexpr std::size_t max_frontier_width = 368;
inline constexpr std::size_t automorphism_key_count = 362;
inline constexpr std::uint32_t manifest_version = 2;
inline constexpr std::uint64_t manifest_magic = UINT64_C(0x4343424756363535);
using TransitionSeed = std::array<std::uint64_t, 4>;

struct ThinSchedule65536 {
    static constexpr std::uint32_t cyclotomic_index = 2 * degree;
    static constexpr std::uint32_t plaintext_prime = 65537;
    static constexpr std::size_t frobenius_order = 2;
    static constexpr std::size_t slot_count = 32768;
    static constexpr std::array<std::size_t, 2> dimensions{16384, 2};
    static constexpr std::array<std::uint32_t, 2> generators{5, 2 * degree - 1};
    static constexpr std::array<std::size_t, 2> baby_dimensions{2, 91};
    static constexpr std::array<std::size_t, 2> giant_dimensions{1, 181};
    static constexpr std::size_t baby_product = 182;
    static constexpr std::size_t giant_product = 181;
    static constexpr std::size_t peak_live_ciphertexts = max_frontier_width;
};

inline std::uint32_t powModCyclotomic(std::uint32_t base, std::size_t exponent)
{
    std::uint32_t result = 1;
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = static_cast<std::uint32_t>(
                (static_cast<std::uint64_t>(result) * base) %
                ThinSchedule65536::cyclotomic_index);
        base = static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(base) * base) %
            ThinSchedule65536::cyclotomic_index);
        exponent >>= 1;
    }
    return result;
}

// Exact switch-key manifest extracted from the non-cyclic thin-BGV
// baby-step/giant-step schedule.  It is generated here, instead of storing a
// 362-entry literal, so the implementation and the schedule formulas cannot
// silently diverge.
inline std::vector<std::uint32_t> thinScheduleSwitchExponents()
{
    std::vector<std::uint32_t> result;
    result.reserve(automorphism_key_count + 1);
    const auto append_rotation = [&result](const std::size_t minus_one_power,
                                           const std::size_t five_power) {
        std::uint32_t exponent = powModCyclotomic(5, five_power);
        if ((minus_one_power & 1U) != 0)
            exponent = ThinSchedule65536::cyclotomic_index - exponent;
        result.push_back(exponent);
    };

    // One-based mixed-radix indices 2..182 in dimensions (2,91).
    for (std::size_t index = 2; index <= ThinSchedule65536::baby_product;
         ++index) {
        const std::size_t value = index - 1;
        append_rotation(value % 2, value / 2);
    }
    // Giant rotations: baby_dimensions * mixed_radix(index,(1,181)).
    for (std::size_t index = 2; index <= ThinSchedule65536::giant_product;
         ++index)
        append_rotation(0, ThinSchedule65536::baby_dimensions[1] * (index - 1));

    // Bad-dimension backward key and the d=2 trace key.  At this modulus,
    // 5^(-16384) = 5^16384.
    result.push_back(powModCyclotomic(5, ThinSchedule65536::dimensions[0]));
    result.push_back(ThinSchedule65536::plaintext_prime);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    if (result.size() != automorphism_key_count)
        throw std::logic_error("thin-BGV switch manifest size mismatch");
    return result;
}

inline std::size_t generatedSubgroupSize(
    const std::span<const std::uint32_t> generators)
{
    std::vector<bool> seen(ThinSchedule65536::cyclotomic_index, false);
    std::vector<std::uint32_t> frontier{1};
    seen[1] = true;
    for (std::size_t cursor = 0; cursor < frontier.size(); ++cursor) {
        for (const auto generator : generators) {
            const auto next = static_cast<std::uint32_t>(
                (static_cast<std::uint64_t>(frontier[cursor]) * generator) %
                ThinSchedule65536::cyclotomic_index);
            if (!seen[next]) {
                seen[next] = true;
                frontier.push_back(next);
            }
        }
    }
    return frontier.size();
}

inline std::size_t checkedProduct(
    const std::initializer_list<std::size_t> values)
{
    std::size_t product = 1;
    for (const std::size_t value : values) {
        if (value != 0 &&
            product > std::numeric_limits<std::size_t>::max() / value)
            throw std::overflow_error("compact-cover size overflow");
        product *= value;
    }
    return product;
}

struct ResourceEstimate {
    std::size_t element_bytes{};
    std::size_t ciphertext_bytes{};
    std::size_t seeded_key_body_bytes{};
    std::size_t seeded_key_seed_bytes{};
    std::size_t uniform_key_bytes{};
};

inline ResourceEstimate estimateResources(const std::size_t width,
                                          const std::size_t limbs,
                                          const std::size_t rows)
{
    if (width == 0 || width > max_frontier_width)
        throw std::invalid_argument("invalid compact-cover frontier width");
    if (limbs == 0 || limbs > modular_ntt::degree65536_primes.size())
        throw std::invalid_argument("invalid compact-cover RNS limb count");
    const std::size_t residues = checkedProduct({degree, width, limbs});
    const std::size_t element_bytes =
        checkedProduct({residues, sizeof(std::uint64_t)});
    return {
        element_bytes,
        checkedProduct({2, element_bytes}),
        checkedProduct({rows, element_bytes}),
        checkedProduct({rows, width, limbs, sizeof(TransitionSeed)}),
        checkedProduct({2, rows, element_bytes}),
    };
}

struct FrontierManifest {
    std::size_t width{};
    std::size_t limbs{};
    std::vector<std::uint32_t> labels{};  // odd exponents modulo 2N

    void validate() const
    {
        if (width == 0 || width > max_frontier_width || labels.size() != width)
            throw std::invalid_argument("invalid compact-cover manifest width");
        if (limbs == 0 || limbs > modular_ntt::degree65536_primes.size())
            throw std::invalid_argument("invalid compact-cover manifest limbs");
        std::vector<std::uint32_t> sorted = labels;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::invalid_argument(
                "duplicate compact-cover frontier label");
        for (const auto label : labels)
            if ((label & 1U) == 0 || label >= 2 * degree)
                throw std::invalid_argument(
                    "invalid cyclotomic automorphism label");
    }
};

class FrontierElement {
public:
    FrontierElement() = default;
    FrontierElement(const std::size_t width, const std::size_t limbs)
        : width_(width),
          limbs_(limbs),
          data_(checkedProduct({width, limbs, degree}))
    {
        if (width == 0 || width > max_frontier_width)
            throw std::invalid_argument("invalid frontier width");
        if (limbs == 0 || limbs > modular_ntt::degree65536_primes.size())
            throw std::invalid_argument("invalid frontier limb count");
    }

    std::size_t width() const { return width_; }
    std::size_t limbs() const { return limbs_; }
    std::size_t bytes() const { return data_.size() * sizeof(std::uint64_t); }

    std::span<std::uint64_t> slice(const std::size_t limb,
                                   const std::size_t label)
    {
        return std::span<std::uint64_t>(data_.data() + offset(limb, label),
                                        degree);
    }
    std::span<const std::uint64_t> slice(const std::size_t limb,
                                         const std::size_t label) const
    {
        return std::span<const std::uint64_t>(
            data_.data() + offset(limb, label), degree);
    }

    std::vector<std::uint64_t> &storage() { return data_; }
    const std::vector<std::uint64_t> &storage() const { return data_; }

private:
    std::size_t offset(const std::size_t limb, const std::size_t label) const
    {
        if (limb >= limbs_ || label >= width_)
            throw std::out_of_range("frontier slice index");
        return (limb * width_ + label) * degree;
    }
    std::size_t width_{};
    std::size_t limbs_{};
    std::vector<std::uint64_t> data_{};
};

struct FrontierCiphertext {
    FrontierElement mask;
    FrontierElement body;

    FrontierCiphertext(const std::size_t width, const std::size_t limbs)
        : mask(width, limbs), body(width, limbs)
    {
    }
};

inline void requireSameShape(const FrontierElement &left,
                             const FrontierElement &right)
{
    if (left.width() != right.width() || left.limbs() != right.limbs())
        throw std::invalid_argument("compact-cover shape mismatch");
}

inline void add(FrontierElement &result, const FrontierElement &left,
                const FrontierElement &right)
{
    requireSameShape(left, right);
    requireSameShape(result, left);
    for (std::size_t limb = 0; limb < left.limbs(); ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        for (std::size_t label = 0; label < left.width(); ++label) {
            const auto lhs = left.slice(limb, label);
            const auto rhs = right.slice(limb, label);
            auto out = result.slice(limb, label);
#ifdef USE_HEXL
            intel::hexl::EltwiseAddMod(out.data(), lhs.data(), rhs.data(),
                                       degree, modulus);
#else
            for (std::size_t i = 0; i < degree; ++i)
                out[i] = modular_ntt::add(lhs[i], rhs[i], modulus);
#endif
        }
    }
}

inline void multiply(FrontierElement &result, const FrontierElement &left,
                     const FrontierElement &right)
{
    requireSameShape(left, right);
    requireSameShape(result, left);
    for (std::size_t limb = 0; limb < left.limbs(); ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        for (std::size_t label = 0; label < left.width(); ++label) {
            const auto lhs = left.slice(limb, label);
            const auto rhs = right.slice(limb, label);
            auto out = result.slice(limb, label);
#ifdef USE_HEXL
            intel::hexl::EltwiseMultMod(out.data(), lhs.data(), rhs.data(),
                                        degree, modulus, 1);
#else
            for (std::size_t i = 0; i < degree; ++i)
                out[i] = modular_ntt::multiply(lhs[i], rhs[i], modulus);
#endif
        }
    }
}

inline std::size_t automorphismSlot(const std::uint32_t exponent,
                                    const std::size_t natural_slot)
{
    const std::uint64_t odd = 2 * natural_slot + 1;
    const std::uint64_t mapped =
        (static_cast<std::uint64_t>(exponent) * odd) % (2 * degree);
    return static_cast<std::size_t>((mapped - 1) / 2);
}

struct RelabelEntry {
    std::size_t source_index{};
    std::uint32_t automorphism{};
};

inline void relabel(FrontierElement &target, const FrontierElement &source,
                    const std::span<const RelabelEntry> mapping)
{
    if (mapping.size() != target.width() || source.limbs() != target.limbs())
        throw std::invalid_argument("invalid compact-cover relabel mapping");
    for (std::size_t target_label = 0; target_label < target.width();
         ++target_label) {
        const auto entry = mapping[target_label];
        if (entry.source_index >= source.width() ||
            (entry.automorphism & 1U) == 0)
            throw std::invalid_argument("invalid compact-cover relabel entry");
        for (std::size_t limb = 0; limb < target.limbs(); ++limb) {
            const auto input = source.slice(limb, entry.source_index);
            auto output = target.slice(limb, target_label);
            for (std::size_t slot = 0; slot < degree; ++slot)
                output[slot] =
                    input[automorphismSlot(entry.automorphism, slot)];
        }
    }
}

inline void forwardNTT(FrontierElement &value)
{
    for (std::size_t limb = 0; limb < value.limbs(); ++limb) {
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        for (std::size_t label = 0; label < value.width(); ++label)
            plan.forward(value.slice(limb, label));
    }
}

inline void inverseNTT(FrontierElement &value)
{
    for (std::size_t limb = 0; limb < value.limbs(); ++limb) {
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        for (std::size_t label = 0; label < value.width(); ++label)
            plan.inverse(value.slice(limb, label));
    }
}

inline std::uint64_t splitmix64(std::uint64_t value)
{
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

inline void expandSeededUniformSlice(const TransitionSeed &seed,
                                     const std::uint64_t modulus,
                                     const std::span<std::uint64_t> output)
{
    // Reject the short prefix so the accepted 64-bit domain is an exact
    // multiple of modulus.  This removes the visible bias of x % modulus for
    // the 61-bit RNS primes.
    const std::uint64_t rejection_threshold = -modulus % modulus;
#ifdef USE_BLAKE3
    std::array<std::uint8_t, 32> seed_bytes{};
    for (std::size_t word = 0; word < seed.size(); ++word)
        for (std::size_t byte = 0; byte < 8; ++byte)
            seed_bytes[8 * word + byte] =
                static_cast<std::uint8_t>(seed[word] >> (8 * byte));
    BLAKE3PRNG::BLAKE3PRNG<std::uint64_t> engine(seed_bytes);
    for (auto &value : output) {
        std::uint64_t sample;
        do sample = engine();
        while (sample < rejection_threshold);
        value = sample % modulus;
    }
#else
    std::uint64_t state = splitmix64(seed[0]) ^ splitmix64(seed[1] + 1) ^
                          splitmix64(seed[2] + 2) ^ splitmix64(seed[3] + 3);
    for (auto &value : output) {
        std::uint64_t sample;
        do {
            state += UINT64_C(0x9e3779b97f4a7c15);
            sample = splitmix64(state);
        } while (sample < rejection_threshold);
        value = sample % modulus;
    }
#endif
}

struct TransitionKeyHeader {
    std::uint64_t magic{manifest_magic};
    std::uint32_t version{manifest_version};
    std::uint32_t reserved{};
    std::uint64_t width{};
    std::uint64_t limbs{};
    std::uint64_t rows{};
    std::uint64_t degree_value{degree};
};

inline constexpr std::uint32_t seeded_transition_key_format = 1;
inline constexpr std::uint32_t uniform_transition_key_format = 2;

class SeededTransitionKeyWriter {
public:
    SeededTransitionKeyWriter(const std::filesystem::path &path,
                              const std::size_t width, const std::size_t limbs,
                              const std::size_t rows,
                              const std::uintmax_t maximum_bytes)
        : path_(path),
          header_{manifest_magic,
                  manifest_version,
                  seeded_transition_key_format,
                  width,
                  limbs,
                  rows,
                  degree}
    {
        const auto estimate = estimateResources(width, limbs, rows);
        const auto required = sizeof(TransitionKeyHeader) +
                              estimate.seeded_key_seed_bytes +
                              estimate.seeded_key_body_bytes;
        if (required > maximum_bytes)
            throw std::runtime_error(
                "transition key exceeds configured disk budget");
        const auto available =
            std::filesystem::space(path.parent_path()).available;
        if (required > available)
            throw std::runtime_error(
                "insufficient disk space for transition key");
        stream_.open(path, std::ios::binary | std::ios::trunc);
        if (!stream_)
            throw std::runtime_error("cannot create transition key file");
        stream_.write(reinterpret_cast<const char *>(&header_),
                      sizeof(header_));
    }

    void writeSeed(const TransitionSeed &seed)
    {
        stream_.write(reinterpret_cast<const char *>(seed.data()),
                      sizeof(seed));
        if (!stream_)
            throw std::runtime_error("failed writing transition seed");
    }

    void writeBodySlice(const std::span<const std::uint64_t> body)
    {
        if (body.size() != degree)
            throw std::invalid_argument(
                "transition body slice has wrong degree");
        stream_.write(
            reinterpret_cast<const char *>(body.data()),
            static_cast<std::streamsize>(degree * sizeof(std::uint64_t)));
        if (!stream_)
            throw std::runtime_error("failed writing transition body");
    }

private:
    std::filesystem::path path_;
    TransitionKeyHeader header_{};
    std::ofstream stream_{};
};

class SeededTransitionKeyProvider {
public:
    explicit SeededTransitionKeyProvider(const std::filesystem::path &path)
        : path_(path)
    {
        stream_.open(path, std::ios::binary);
        if (!stream_)
            throw std::runtime_error("cannot open transition key file");
        stream_.read(reinterpret_cast<char *>(&header_), sizeof(header_));
        if (!stream_ || header_.magic != manifest_magic ||
            header_.version != manifest_version ||
            header_.reserved != seeded_transition_key_format ||
            header_.degree_value != degree || header_.width == 0 ||
            header_.width > max_frontier_width || header_.limbs == 0 ||
            header_.limbs > modular_ntt::degree65536_primes.size())
            throw std::runtime_error("invalid transition key manifest");
        seeds_offset_ = sizeof(TransitionKeyHeader);
        bodies_offset_ = seeds_offset_ + header_.rows * header_.limbs *
                                             header_.width *
                                             sizeof(TransitionSeed);
        const auto expected_bytes =
            bodies_offset_ + header_.rows * header_.limbs * header_.width *
                                 degree * sizeof(std::uint64_t);
        if (std::filesystem::file_size(path_) != expected_bytes)
            throw std::runtime_error("truncated transition key file");
    }

    const TransitionKeyHeader &header() const { return header_; }

    TransitionSeed seed(const std::size_t row, const std::size_t limb,
                        const std::size_t label) const
    {
        check(row, limb, label);
        const auto index = (row * header_.limbs + limb) * header_.width + label;
        TransitionSeed value{};
        readAt(seeds_offset_ + index * sizeof(value), value.data(),
               sizeof(value));
        return value;
    }

    void expandMaskSlice(const std::size_t row, const std::size_t limb,
                         const std::size_t label,
                         const std::span<std::uint64_t> output) const
    {
        if (output.size() != degree)
            throw std::invalid_argument("mask output slice has wrong degree");
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        expandSeededUniformSlice(seed(row, limb, label), modulus, output);
    }

    void loadBodySlice(const std::size_t row, const std::size_t limb,
                       const std::size_t label,
                       const std::span<std::uint64_t> output) const
    {
        check(row, limb, label);
        if (output.size() != degree)
            throw std::invalid_argument("body output slice has wrong degree");
        const auto index = (row * header_.limbs + limb) * header_.width + label;
        readAt(bodies_offset_ + index * degree * sizeof(std::uint64_t),
               output.data(), degree * sizeof(std::uint64_t));
    }

private:
    void check(const std::size_t row, const std::size_t limb,
               const std::size_t label) const
    {
        if (row >= header_.rows || limb >= header_.limbs ||
            label >= header_.width)
            throw std::out_of_range("transition key slice index");
    }

    void readAt(const std::uint64_t offset, void *output,
                const std::size_t bytes) const
    {
        stream_.clear();
        stream_.seekg(static_cast<std::streamoff>(offset));
        stream_.read(reinterpret_cast<char *>(output),
                     static_cast<std::streamsize>(bytes));
        if (!stream_)
            throw std::runtime_error("failed reading transition key slice");
    }

    std::filesystem::path path_;
    mutable std::ifstream stream_{};
    TransitionKeyHeader header_{};
    std::uint64_t seeds_offset_{};
    std::uint64_t bodies_offset_{};
};

// Information-theoretic key format used by the certified compiler.  Masks
// and bodies are stored interleaved as complete NTT slices in
// [row][limb][label][component] order, so evaluation needs only one row slice
// in memory.  Unlike the seeded format, this adds no PRG assumption.
class UniformTransitionKeyWriter {
public:
    UniformTransitionKeyWriter(const std::filesystem::path &path,
                               const std::size_t width, const std::size_t limbs,
                               const std::size_t rows,
                               const std::uintmax_t maximum_bytes)
        : header_{manifest_magic,
                  manifest_version,
                  uniform_transition_key_format,
                  width,
                  limbs,
                  rows,
                  degree}
    {
        const auto estimate = estimateResources(width, limbs, rows);
        const auto required =
            sizeof(TransitionKeyHeader) + estimate.uniform_key_bytes;
        if (required > maximum_bytes)
            throw std::runtime_error(
                "uniform transition key exceeds configured disk budget");
        const auto parent = path.parent_path().empty()
                                ? std::filesystem::current_path()
                                : path.parent_path();
        if (required > std::filesystem::space(parent).available)
            throw std::runtime_error(
                "insufficient disk space for uniform transition key");
        stream_.open(path, std::ios::binary | std::ios::trunc);
        if (!stream_)
            throw std::runtime_error(
                "cannot create uniform transition key file");
        stream_.write(reinterpret_cast<const char *>(&header_),
                      sizeof(header_));
    }

    void writeSlice(const std::span<const std::uint64_t> mask,
                    const std::span<const std::uint64_t> body)
    {
        if (mask.size() != degree || body.size() != degree)
            throw std::invalid_argument(
                "uniform transition slice has wrong degree");
        stream_.write(reinterpret_cast<const char *>(mask.data()),
                      degree * sizeof(std::uint64_t));
        stream_.write(reinterpret_cast<const char *>(body.data()),
                      degree * sizeof(std::uint64_t));
        if (!stream_)
            throw std::runtime_error("failed writing uniform transition slice");
        ++written_slices_;
    }

    void finish()
    {
        const auto expected = header_.rows * header_.limbs * header_.width;
        if (written_slices_ != expected)
            throw std::runtime_error(
                "uniform transition key has incomplete slice count");
        stream_.flush();
        if (!stream_)
            throw std::runtime_error("failed flushing uniform transition key");
    }

private:
    TransitionKeyHeader header_{};
    std::ofstream stream_{};
    std::size_t written_slices_{};
};

class UniformTransitionKeyProvider {
public:
    explicit UniformTransitionKeyProvider(const std::filesystem::path &path)
        : path_(path)
    {
        stream_.open(path, std::ios::binary);
        if (!stream_)
            throw std::runtime_error("cannot open uniform transition key file");
        stream_.read(reinterpret_cast<char *>(&header_), sizeof(header_));
        if (!stream_ || header_.magic != manifest_magic ||
            header_.version != manifest_version ||
            header_.reserved != uniform_transition_key_format ||
            header_.degree_value != degree || header_.width == 0 ||
            header_.width > max_frontier_width || header_.limbs == 0 ||
            header_.limbs > modular_ntt::degree65536_primes.size() ||
            header_.rows == 0)
            throw std::runtime_error("invalid uniform transition key manifest");
        const auto expected =
            sizeof(TransitionKeyHeader) +
            estimateResources(header_.width, header_.limbs, header_.rows)
                .uniform_key_bytes;
        if (std::filesystem::file_size(path_) != expected)
            throw std::runtime_error("truncated uniform transition key file");
    }

    const TransitionKeyHeader &header() const { return header_; }

    void loadSlice(const std::size_t row, const std::size_t limb,
                   const std::size_t label, const std::span<std::uint64_t> mask,
                   const std::span<std::uint64_t> body) const
    {
        if (row >= header_.rows || limb >= header_.limbs ||
            label >= header_.width || mask.size() != degree ||
            body.size() != degree)
            throw std::out_of_range("uniform transition key slice index");
        const auto index = (row * header_.limbs + limb) * header_.width + label;
        const auto offset = sizeof(TransitionKeyHeader) +
                            index * 2 * degree * sizeof(std::uint64_t);
        readAt(offset, mask.data(), degree * sizeof(std::uint64_t));
        readAt(offset + degree * sizeof(std::uint64_t), body.data(),
               degree * sizeof(std::uint64_t));
    }

private:
    void readAt(const std::uint64_t offset, void *output,
                const std::size_t bytes) const
    {
        stream_.clear();
        stream_.seekg(static_cast<std::streamoff>(offset));
        stream_.read(reinterpret_cast<char *>(output),
                     static_cast<std::streamsize>(bytes));
        if (!stream_)
            throw std::runtime_error(
                "failed reading uniform transition key slice");
    }

    std::filesystem::path path_;
    mutable std::ifstream stream_{};
    TransitionKeyHeader header_{};
};

inline std::size_t gadgetRows(const std::uint64_t modulus,
                              const unsigned gadget_bits)
{
    if (gadget_bits == 0 || gadget_bits > 32)
        throw std::invalid_argument("gadget bits must be in [1,32]");
    unsigned modulus_bits = 0;
    for (std::uint64_t value = modulus - 1; value != 0; value >>= 1)
        ++modulus_bits;
    return (modulus_bits + gadget_bits - 1) / gadget_bits;
}

inline void addProductMod(const std::span<std::uint64_t> accumulator,
                          const std::span<const std::uint64_t> left,
                          const std::span<const std::uint64_t> right,
                          const std::uint64_t modulus,
                          const std::span<std::uint64_t> scratch)
{
    if (accumulator.size() != degree || left.size() != degree ||
        right.size() != degree || scratch.size() != degree)
        throw std::invalid_argument("invalid compact-cover product slice");
#ifdef USE_HEXL
    intel::hexl::EltwiseMultMod(scratch.data(), left.data(), right.data(),
                                degree, modulus, 1);
    intel::hexl::EltwiseAddMod(accumulator.data(), accumulator.data(),
                               scratch.data(), degree, modulus);
#else
    for (std::size_t slot = 0; slot < degree; ++slot)
        accumulator[slot] = modular_ntt::add(
            accumulator[slot],
            modular_ntt::multiply(left[slot], right[slot], modulus), modulus);
#endif
}

// Streamed, limb-local gadget transition in the NTT domain.
//
// Each key row must have phase
//     key_body - key_mask * target_secret = -B^row * source_secret.
// Consequently the returned ciphertext has exactly the input phase, under
// target_secret, apart from evaluation-key error.  Decomposition is applied
// to coefficient-domain residues independently in every RNS limb.  Only two
// key slices and four scratch polynomials are resident beyond input/output.
inline void applyTransition(FrontierCiphertext &output,
                            const FrontierCiphertext &input,
                            const SeededTransitionKeyProvider &key,
                            const unsigned gadget_bits)
{
    requireSameShape(input.mask, input.body);
    requireSameShape(output.mask, output.body);
    requireSameShape(output.mask, input.mask);
    const auto &header = key.header();
    if (header.width != input.mask.width() ||
        header.limbs != input.mask.limbs())
        throw std::invalid_argument("transition key shape mismatch");

    std::vector<std::uint64_t> coefficients(degree);
    std::vector<std::uint64_t> digit(degree);
    std::vector<std::uint64_t> key_mask(degree);
    std::vector<std::uint64_t> key_body(degree);
    std::vector<std::uint64_t> product(degree);
    const std::uint64_t digit_mask = (UINT64_C(1) << gadget_bits) - UINT64_C(1);

    for (std::size_t limb = 0; limb < input.mask.limbs(); ++limb) {
        const auto modulus = modular_ntt::degree65536_primes[limb].value;
        const auto rows = gadgetRows(modulus, gadget_bits);
        if (header.rows != rows)
            throw std::invalid_argument("transition key gadget-row mismatch");
        modular_ntt::NegacyclicNTTPlan plan(
            degree, modular_ntt::degree65536_primes[limb]);
        for (std::size_t label = 0; label < input.mask.width(); ++label) {
            const auto input_mask = input.mask.slice(limb, label);
            const auto input_body = input.body.slice(limb, label);
            auto output_mask = output.mask.slice(limb, label);
            auto output_body = output.body.slice(limb, label);
            std::fill(output_mask.begin(), output_mask.end(), 0);
            std::copy(input_body.begin(), input_body.end(),
                      output_body.begin());
            std::copy(input_mask.begin(), input_mask.end(),
                      coefficients.begin());
            plan.inverse(coefficients);

            for (std::size_t row = 0; row < rows; ++row) {
                const unsigned shift = static_cast<unsigned>(row) * gadget_bits;
                for (std::size_t coefficient = 0; coefficient < degree;
                     ++coefficient)
                    digit[coefficient] =
                        (coefficients[coefficient] >> shift) & digit_mask;
                plan.forward(digit);
                key.expandMaskSlice(row, limb, label, key_mask);
                key.loadBodySlice(row, limb, label, key_body);
                addProductMod(output_mask, digit, key_mask, modulus, product);
                addProductMod(output_body, digit, key_body, modulus, product);
            }
        }
    }
}

class FullModulusGadget {
public:
    FullModulusGadget(const std::size_t limbs, const std::size_t rows)
        : limbs_(limbs), rows_(rows)
    {
        if (limbs == 0 || limbs > modular_ntt::degree65536_primes.size() ||
            rows == 0)
            throw std::invalid_argument("invalid full-modulus gadget shape");
        modulus_ = 1;
        for (std::size_t limb = 0; limb < limbs_; ++limb)
            modulus_ *= modular_ntt::degree65536_primes[limb].value;
        modulus_bits_ =
            static_cast<unsigned>(boost::multiprecision::msb(modulus_) + 1);
        gadget_bits_ =
            static_cast<unsigned>((modulus_bits_ + rows_ - 1) / rows_);
        base_ = boost::multiprecision::cpp_int{1} << gadget_bits_;
    }

    std::size_t limbs() const { return limbs_; }
    std::size_t rows() const { return rows_; }
    unsigned modulusBits() const { return modulus_bits_; }
    unsigned gadgetBits() const { return gadget_bits_; }
    const boost::multiprecision::cpp_int &modulus() const { return modulus_; }
    const boost::multiprecision::cpp_int &base() const { return base_; }

    boost::multiprecision::cpp_int reconstructCentered(
        const std::span<const std::uint64_t> residues) const
    {
        using boost::multiprecision::cpp_int;
        if (residues.size() != limbs_)
            throw std::invalid_argument("CRT residue count mismatch");
        cpp_int value = residues[0];
        cpp_int accumulated = modular_ntt::degree65536_primes[0].value;
        for (std::size_t limb = 1; limb < limbs_; ++limb) {
            const auto prime = modular_ntt::degree65536_primes[limb].value;
            const auto current = static_cast<std::uint64_t>(
                (value % prime).convert_to<std::uint64_t>());
            const auto difference =
                modular_ntt::subtract(residues[limb], current, prime);
            const auto accumulated_mod = static_cast<std::uint64_t>(
                (accumulated % prime).convert_to<std::uint64_t>());
            const auto quotient = modular_ntt::multiply(
                difference, modular_ntt::invert(accumulated_mod, prime), prime);
            value += accumulated * quotient;
            accumulated *= prime;
        }
        if (value > modulus_ / 2) value -= modulus_;
        return value;
    }

    static std::uint64_t residue(const boost::multiprecision::cpp_int &value,
                                 const std::uint64_t modulus)
    {
        boost::multiprecision::cpp_int reduced = value % modulus;
        if (reduced < 0) reduced += modulus;
        return reduced.convert_to<std::uint64_t>();
    }

    boost::multiprecision::cpp_int takeBalancedDigit(
        boost::multiprecision::cpp_int &remaining, const bool final_row) const
    {
        using boost::multiprecision::cpp_int;
        if (final_row) {
            cpp_int digit = remaining;
            remaining = 0;
            return digit;
        }
        cpp_int digit = remaining % base_;
        if (digit < 0) digit += base_;
        if (digit > base_ / 2) digit -= base_;
        remaining = (remaining - digit) / base_;
        return digit;
    }

private:
    std::size_t limbs_{};
    std::size_t rows_{};
    unsigned modulus_bits_{};
    unsigned gadget_bits_{};
    boost::multiprecision::cpp_int modulus_{};
    boost::multiprecision::cpp_int base_{};
};

// Certified transition path.  Unlike `applyTransition`, this decomposes one
// centered coefficient modulo the complete CRT product, so every digit is a
// genuine small polynomial shared by all RNS limbs.
inline void applyFullModulusTransition(FrontierCiphertext &output,
                                       const FrontierCiphertext &input,
                                       const UniformTransitionKeyProvider &key)
{
    requireSameShape(input.mask, input.body);
    requireSameShape(output.mask, output.body);
    requireSameShape(output.mask, input.mask);
    const auto &header = key.header();
    if (header.width != input.mask.width() ||
        header.limbs != input.mask.limbs())
        throw std::invalid_argument("uniform transition key shape mismatch");
    FullModulusGadget gadget(header.limbs, header.rows);

    std::vector<std::vector<std::uint64_t>> coefficients(
        header.limbs, std::vector<std::uint64_t>(degree));
    std::vector<std::vector<std::uint64_t>> digit_ntt(
        header.limbs, std::vector<std::uint64_t>(degree));
    std::vector<boost::multiprecision::cpp_int> remaining(degree);
    std::vector<std::uint64_t> residues(header.limbs);
    std::vector<std::uint64_t> key_mask(degree), key_body(degree),
        product(degree);

    for (std::size_t label = 0; label < input.mask.width(); ++label) {
        for (std::size_t limb = 0; limb < header.limbs; ++limb) {
            std::copy(input.mask.slice(limb, label).begin(),
                      input.mask.slice(limb, label).end(),
                      coefficients[limb].begin());
            modular_ntt::NegacyclicNTTPlan plan(
                degree, modular_ntt::degree65536_primes[limb]);
            plan.inverse(coefficients[limb]);
            std::fill(output.mask.slice(limb, label).begin(),
                      output.mask.slice(limb, label).end(), 0);
            std::copy(input.body.slice(limb, label).begin(),
                      input.body.slice(limb, label).end(),
                      output.body.slice(limb, label).begin());
        }
        for (std::size_t coefficient = 0; coefficient < degree; ++coefficient) {
            for (std::size_t limb = 0; limb < header.limbs; ++limb)
                residues[limb] = coefficients[limb][coefficient];
            remaining[coefficient] = gadget.reconstructCentered(residues);
        }

        for (std::size_t row = 0; row < header.rows; ++row) {
            for (std::size_t coefficient = 0; coefficient < degree;
                 ++coefficient) {
                const auto digit = gadget.takeBalancedDigit(
                    remaining[coefficient], row + 1 == header.rows);
                for (std::size_t limb = 0; limb < header.limbs; ++limb)
                    digit_ntt[limb][coefficient] = FullModulusGadget::residue(
                        digit, modular_ntt::degree65536_primes[limb].value);
            }
            for (std::size_t limb = 0; limb < header.limbs; ++limb) {
                const auto modulus =
                    modular_ntt::degree65536_primes[limb].value;
                modular_ntt::NegacyclicNTTPlan plan(
                    degree, modular_ntt::degree65536_primes[limb]);
                plan.forward(digit_ntt[limb]);
                key.loadSlice(row, limb, label, key_mask, key_body);
                addProductMod(output.mask.slice(limb, label), digit_ntt[limb],
                              key_mask, modulus, product);
                addProductMod(output.body.slice(limb, label), digit_ntt[limb],
                              key_body, modulus, product);
            }
        }
        if (std::any_of(remaining.begin(), remaining.end(),
                        [](const auto &value) { return value != 0; }))
            throw std::runtime_error(
                "full-modulus gadget did not cover the CRT modulus");
    }
}

}  // namespace TFHEpp::compact_cover_bgv
