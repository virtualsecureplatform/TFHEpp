#pragma once

#include "narrow_selector_file.hpp"
#include <bit>

#if defined(TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS) && defined(USE_BLAKE3)
#include <blake3.h>

namespace TFHEpp {

// A derived public cache, not a key archive. Unkeyed BLAKE3 detects accidental
// payload corruption and binds each record to its raw selector; it does not
// authenticate a directory supplied by an adversary. Fixed dimensions avoid
// allocations controlled by serialized vector lengths.
template <class Schedule>
class CKKSNarrowNTTSelectorFile final : public CKKSNarrowSelectorSource<Schedule> {
    using P = typename Schedule::Param;
    using Cache = CKKSNarrowNTTSelector<Schedule>;
    using Raw = CKKSNarrowSelectorFile<Schedule>;
    using Digest = std::array<std::uint8_t, 32>;
    std::filesystem::path directory_;
    std::shared_ptr<const Raw> raw_;
    static constexpr std::size_t words = (P::k + 1) * Cache::digits * P::n;
    static constexpr std::size_t record_bytes = 64 + Cache::prime_count * words * sizeof(std::uint64_t);
    static constexpr std::size_t header_bytes = 1 + 17 * sizeof(std::uint64_t);

    static std::filesystem::path path(const std::filesystem::path &directory, bool tail, std::size_t group)
    {
        if (group >= (tail ? Schedule::tail_packs : Schedule::blocks))
            throw std::out_of_range("prepared selector group");
        return directory / ((tail ? "tail-" : "head-") + std::to_string(group) + ".bin");
    }
    static std::array<std::uint64_t, 17> header(bool tail, std::size_t group)
    {
        return {0x4e41524e54543031ULL, P::n, Schedule::selector_log_q, Schedule::log_delta,
            Schedule::plain_log_delta, Schedule::block_width, std::numeric_limits<typename P::T>::digits,
            Cache::digit_bits, Cache::digits, Cache::prime_count,
            modular_ntt::degree65536_primes[0].value, modular_ntt::degree65536_primes[0].primitive_root,
            modular_ntt::degree65536_primes[1].value, modular_ntt::degree65536_primes[1].primitive_root,
            1, tail, group}; // format1: natural negacyclic NTT index order
    }
    static void hashWords(blake3_hasher &hasher, std::span<const std::uint64_t> input)
    {
        if constexpr (std::endian::native == std::endian::little)
            blake3_hasher_update(&hasher, input.data(), input.size_bytes());
        else for (auto word : input) {
            std::array<std::uint8_t, 8> bytes{};
            for (std::size_t i = 0; i < bytes.size(); ++i) bytes[i] = word >> (8 * i);
            blake3_hasher_update(&hasher, bytes.data(), bytes.size());
        }
    }
    static Digest rawDigest(const typename Schedule::Selector &selector)
    {
        static_assert(std::is_trivially_copyable_v<typename Schedule::Selector>);
        static_assert(sizeof(selector) == (P::k + 1) * P::n * sizeof(typename P::T));
        static_assert(sizeof(typename P::T) == sizeof(std::declval<typename P::T>().limb));
        // Hash typed limbs, avoiding aliasing unrelated aggregate storage.
        blake3_hasher hasher; blake3_hasher_init(&hasher);
        for (const auto &component : selector.ct) {
            if constexpr (std::endian::native == std::endian::little)
                blake3_hasher_update(&hasher, component.data(), sizeof(component));
            else for (const auto &coefficient : component) hashWords(hasher, coefficient.limb);
        }
        Digest result{}; blake3_hasher_finalize(&hasher, result.data(), result.size());
        return result;
    }
    static Digest cacheDigest(const Cache &selector)
    {
        blake3_hasher hasher; blake3_hasher_init(&hasher);
        for (const auto &spectrum : selector.spectra) hashWords(hasher, spectrum);
        Digest result{}; blake3_hasher_finalize(&hasher, result.data(), result.size());
        return result;
    }
    static void validateHeader(const std::filesystem::path &file_path, bool tail, std::size_t group,
                               cereal::PortableBinaryInputArchive &archive)
    {
        std::array<std::uint64_t, 17> actual{}; archive(actual);
        if (actual != header(tail, group) || std::filesystem::file_size(file_path) !=
            header_bytes + Schedule::block_width * record_bytes)
            throw std::runtime_error("prepared selector profile or length mismatch");
    }
    static void writeRecords(cereal::PortableBinaryOutputArchive &archive,
                             std::span<const typename Schedule::Selector> raw,
                             std::span<const Cache> prepared)
    {
        if (raw.size() != prepared.size()) throw std::invalid_argument("prepared record count mismatch");
        for (std::size_t i = 0; i < raw.size(); ++i) {
            auto source_hash = rawDigest(raw[i]);
            auto payload_hash = cacheDigest(prepared[i]);
            archive(cereal::binary_data(source_hash.data(), source_hash.size()));
            for (const auto &spectrum : prepared[i].spectra) {
                if (spectrum.size() != words) throw std::invalid_argument("prepared record length mismatch");
                archive(cereal::binary_data(spectrum.data(), spectrum.size()*sizeof(std::uint64_t)));
            }
            archive(cereal::binary_data(payload_hash.data(), payload_hash.size()));
        }
    }
    static std::vector<Cache> readRecords(cereal::PortableBinaryInputArchive &archive,
                                         std::span<const typename Schedule::Selector> raw)
    {
        std::vector<Cache> loaded(raw.size());
        for (std::size_t i = 0; i < raw.size(); ++i) {
            Digest source_hash{}, payload_hash{};
            archive(cereal::binary_data(source_hash.data(), source_hash.size()));
            for (std::size_t prime = 0; prime < Cache::prime_count; ++prime) {
                auto &spectrum = loaded[i].spectra[prime]; spectrum.resize(words);
                archive(cereal::binary_data(spectrum.data(), spectrum.size()*sizeof(std::uint64_t)));
                for (const auto value : spectrum)
                    if (value >= modular_ntt::degree65536_primes[prime].value)
                        throw std::runtime_error("noncanonical prepared selector residue");
            }
            archive(cereal::binary_data(payload_hash.data(), payload_hash.size()));
            if (source_hash != rawDigest(raw[i]) || payload_hash != cacheDigest(loaded[i]))
                throw std::runtime_error("prepared selector source or payload digest mismatch");
        }
        return loaded;
    }

public:
    struct Existing {};
    CKKSNarrowNTTSelectorFile(std::filesystem::path directory, std::shared_ptr<const Raw> raw, Existing)
        : directory_(std::move(directory)), raw_(std::move(raw))
    {
        if (!raw_) throw std::invalid_argument("prepared selectors require their raw source");
        for (const bool tail : {false, true})
            for (std::size_t group = 0; group < (tail ? Schedule::tail_packs : Schedule::blocks); ++group) {
                std::ifstream file(path(directory_, tail, group), std::ios::binary);
                if (!file) throw std::runtime_error("missing prepared selector group");
                cereal::PortableBinaryInputArchive archive(file);
                validateHeader(path(directory_, tail, group), tail, group, archive);
            }
    }
    static std::shared_ptr<CKKSNarrowNTTSelectorFile> prepare(
        const std::filesystem::path &directory, std::shared_ptr<const Raw> raw, std::size_t chunk_size)
    {
        if (!raw || chunk_size == 0) throw std::invalid_argument("invalid prepared selector source or chunk");
        if (!std::filesystem::create_directory(directory))
            throw std::runtime_error("prepared selector directory must be new");
        std::filesystem::permissions(directory, std::filesystem::perms::owner_all);
        struct Cleanup {
            std::filesystem::path directory;
            bool keep = false;
            ~Cleanup() {
                if (keep) return;
                std::error_code error;
                for (bool tail : {false, true})
                    for (std::size_t group = 0; group < (tail ? Schedule::tail_packs : Schedule::blocks); ++group) {
                        const auto file = path(directory, tail, group);
                        std::filesystem::remove(file, error);
                        std::filesystem::remove(file.string()+".pending", error);
                    }
                std::filesystem::remove(directory, error);
            }
        } cleanup{directory};
        for (bool tail : {false, true})
            for (std::size_t group = 0; group < (tail ? Schedule::tail_packs : Schedule::blocks); ++group) {
                const auto destination = path(directory, tail, group);
                const auto pending = destination.string()+".pending";
                std::ofstream file(pending, std::ios::binary);
                if (!file) throw std::runtime_error("cannot create prepared selector group");
                std::filesystem::permissions(pending, std::filesystem::perms::owner_read | std::filesystem::perms::owner_write);
                {
                    cereal::PortableBinaryOutputArchive archive(file); archive(header(tail, group));
                    for (std::size_t first = 0; first < Schedule::block_width;) {
                        const auto count = std::min(chunk_size, Schedule::block_width-first);
                        std::vector<typename Schedule::Selector> selectors;
                        raw->loadChunk(tail, group, first, count, selectors);
                        std::vector<Cache> prepared;
                        narrow_detail::prepareSelectors<Schedule>(prepared, selectors);
                        writeRecords(archive, selectors, prepared);
                        first += count;
                    }
                }
                file.flush(); if (!file) throw std::runtime_error("cannot flush prepared selector group");
                file.close(); if (!file) throw std::runtime_error("cannot close prepared selector group");
                std::filesystem::rename(pending, destination);
            }
        auto result = std::make_shared<CKKSNarrowNTTSelectorFile>(directory, std::move(raw), Existing{});
        cleanup.keep = true;
        return result;
    }
    std::size_t size(bool tail) const override { return raw_->size(tail); }
    // A separately identified partial archive supports bounded preparation and
    // inspection. It cannot be attached as a complete evaluation source.
    static std::filesystem::path prepareChunk(const std::filesystem::path &directory,
        const Raw &source, bool tail, std::size_t group, std::size_t first, std::size_t count)
    {
        if (count == 0 || first > Schedule::block_width || count > Schedule::block_width-first)
            throw std::out_of_range("prepared selector chunk range");
        std::vector<typename Schedule::Selector> raw;
        source.loadChunk(tail, group, first, count, raw);
        std::vector<Cache> prepared;
        narrow_detail::prepareSelectors<Schedule>(prepared, raw);
        if (!std::filesystem::create_directory(directory))
            throw std::runtime_error("prepared chunk directory must be new");
        std::filesystem::permissions(directory, std::filesystem::perms::owner_all);
        const auto pending = directory/"chunk.pending", destination = directory/"chunk.bin";
        try {
            std::ofstream file(pending, std::ios::binary);
            if (!file) throw std::runtime_error("cannot create prepared chunk");
            std::filesystem::permissions(pending, std::filesystem::perms::owner_read | std::filesystem::perms::owner_write);
            {
                cereal::PortableBinaryOutputArchive archive(file);
                auto metadata = header(tail, group); metadata[0] = 0x4e41524e54433031ULL;
                const std::array<std::uint64_t,2> range{first,count}; archive(metadata,range);
                writeRecords(archive, raw, prepared);
            }
            file.flush(); if (!file) throw std::runtime_error("cannot flush prepared chunk");
            file.close(); if (!file) throw std::runtime_error("cannot close prepared chunk");
            std::filesystem::rename(pending,destination);
        }
        catch (...) {
            std::error_code error;
            std::filesystem::remove(pending,error); std::filesystem::remove(directory,error);
            throw;
        }
        return destination;
    }
    static void loadChunkArchive(const std::filesystem::path &file_path, const Raw &source,
        bool tail, std::size_t group, std::size_t first, std::size_t count, std::vector<Cache> &result)
    {
        if (count == 0 || first > Schedule::block_width || count > Schedule::block_width-first)
            throw std::out_of_range("prepared selector chunk range");
        std::ifstream file(file_path,std::ios::binary);
        if (!file) throw std::runtime_error("cannot read prepared chunk");
        cereal::PortableBinaryInputArchive archive(file);
        std::array<std::uint64_t,17> metadata{};std::array<std::uint64_t,2> range{};
        archive(metadata,range);
        auto expected = header(tail,group); expected[0] = 0x4e41524e54433031ULL;
        if (metadata != expected || range != std::array<std::uint64_t,2>{first,count} ||
            std::filesystem::file_size(file_path) != header_bytes + sizeof(range) + count * record_bytes)
            throw std::runtime_error("prepared chunk profile, range or length mismatch");
        std::vector<typename Schedule::Selector> raw;
        source.loadChunk(tail,group,first,count,raw);
        auto loaded = readRecords(archive,raw);
        result.swap(loaded);
    }
    bool hasPreparedChunks() const override { return true; }
    void load(bool tail, std::size_t group, std::vector<typename Schedule::Selector> &result) const override
    { raw_->load(tail, group, result); }
    void loadPreparedChunk(bool tail, std::size_t group, std::size_t first, std::size_t count,
                           std::vector<Cache> &result) const override
    {
        if (count == 0 || first > Schedule::block_width || count > Schedule::block_width-first)
            throw std::out_of_range("prepared selector chunk range");
        const auto file_path = path(directory_, tail, group);
        std::ifstream file(file_path, std::ios::binary);
        if (!file) throw std::runtime_error("cannot read prepared selector group");
        cereal::PortableBinaryInputArchive archive(file);
        validateHeader(file_path, tail, group, archive);
        file.seekg(header_bytes + first * record_bytes);
        if (!file) throw std::runtime_error("cannot seek prepared selector chunk");
        std::vector<typename Schedule::Selector> raw;
        raw_->loadChunk(tail, group, first, count, raw);
        auto loaded = readRecords(archive, raw);
        result.swap(loaded);
    }
};

} // namespace TFHEpp
#endif
