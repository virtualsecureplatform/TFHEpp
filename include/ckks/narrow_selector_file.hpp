#pragma once

#include "narrow_bootstrapping.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// A caller-owned directory name must be new. The cache owns only its known
// files and removes them on destruction. Portable fixed-size records avoid
// native-layout assumptions and untrusted vector-length allocations.
template <class Schedule>
class CKKSNarrowSelectorFile final : public CKKSNarrowSelectorSource<Schedule> {
    using P = typename Schedule::Param;
    using Selector = typename Schedule::Selector;
    std::filesystem::path directory_;
    std::vector<bool> ready_;
    bool complete_ = false;
    mutable bool remove_on_destroy_ = true;

    static constexpr std::size_t groups = Schedule::blocks + Schedule::tail_packs;
    static std::size_t index(bool tail, std::size_t group)
    {
        if (group >= (tail ? Schedule::tail_packs : Schedule::blocks))
            throw std::out_of_range("narrow selector file group");
        return (tail ? Schedule::blocks : 0) + group;
    }
    std::filesystem::path path(bool tail, std::size_t group) const
    {
        return directory_ / ((tail ? "tail-" : "head-") + std::to_string(group) + ".bin");
    }
    static std::array<std::uint64_t, 9> header(bool tail, std::size_t group)
    {
        return {0x4e4152524f573031ULL, P::n, Schedule::selector_log_q,
                Schedule::log_delta, Schedule::plain_log_delta, Schedule::block_width,
                std::numeric_limits<typename P::T>::digits, tail, group};
    }

public:
    struct Existing {};
    explicit CKKSNarrowSelectorFile(std::filesystem::path directory)
        : directory_(std::move(directory)), ready_(groups, false)
    {
        if (!std::filesystem::create_directory(directory_))
            throw std::runtime_error("narrow selector cache directory must not exist");
    }
    // Reopen a retained cache without taking deletion ownership. Validate all
    // group headers and fixed lengths before publishing a complete source;
    // load() still validates the payload and rejects trailing data.
    CKKSNarrowSelectorFile(std::filesystem::path directory, Existing)
        : directory_(std::move(directory)), ready_(groups, true), remove_on_destroy_(false)
    {
        for (bool tail : {false, true})
            for (std::size_t group = 0; group < (tail ? Schedule::tail_packs : Schedule::blocks); group++) {
                std::ifstream file(path(tail, group), std::ios::binary);
                if (!file) throw std::runtime_error("cannot reopen narrow selector group");
                cereal::PortableBinaryInputArchive archive(file);
                std::array<std::uint64_t, 9> metadata;
                archive(metadata);
                if (metadata != header(tail, group) ||
                    std::filesystem::file_size(path(tail, group)) !=
                        1 + sizeof(metadata) + std::size_t{Schedule::block_width} * sizeof(Selector))
                    throw std::runtime_error("retained narrow selector group metadata or length mismatch");
            }
        complete_ = true;
    }
    CKKSNarrowSelectorFile(const CKKSNarrowSelectorFile &) = delete;
    CKKSNarrowSelectorFile &operator=(const CKKSNarrowSelectorFile &) = delete;
    ~CKKSNarrowSelectorFile() override
    {
        if (!remove_on_destroy_) return;
        std::error_code error;
        for (bool tail : {false, true})
            for (std::size_t group = 0; group < (tail ? Schedule::tail_packs : Schedule::blocks); group++)
                std::filesystem::remove(path(tail, group), error);
        std::filesystem::remove(directory_, error); // Preserve any unrelated files.
    }
    // Diagnostic checkpoints retain completed public selectors for replay.
    // Incomplete caches keep their normal cleanup behavior.
    void preserve() const
    {
        if (!complete_) throw std::logic_error("cannot preserve incomplete narrow selector cache");
        remove_on_destroy_ = false;
    }
    void write(bool tail, std::size_t group, std::span<const Selector> selectors)
    {
        const auto offset = index(tail, group);
        if (complete_ || ready_[offset] || selectors.size() != Schedule::block_width)
            throw std::invalid_argument("invalid narrow selector file write");
        std::ofstream file(path(tail, group), std::ios::binary | std::ios::trunc);
        if (!file) throw std::runtime_error("cannot create narrow selector group file");
        {
            cereal::PortableBinaryOutputArchive archive(file);
            const auto metadata = header(tail, group);
            archive(metadata);
            if constexpr (is_multilimb_uint_v<typename P::T>) {
                // The legacy archive writes each coefficient's uint64 limbs
                // in this order. Transfer the same word sequence in one call;
                // PortableBinary still converts endianness per uint64 word.
                static_assert(std::is_trivially_copyable_v<Selector> &&
                    std::is_standard_layout_v<Selector>);
                static_assert(sizeof(Selector) == (P::k + 1) * P::n * sizeof(typename P::T));
                static_assert(sizeof(typename P::T) == sizeof(selectors.front().ct[0][0].limb));
                archive(cereal::binary_data(reinterpret_cast<const std::uint64_t *>(selectors.data()),
                                            selectors.size_bytes()));
            }
            else for (const auto &selector : selectors) archive(selector);
        }
        file.flush();
        if (!file) throw std::runtime_error("cannot write narrow selector group file");
        ready_[offset] = true;
    }
    void finish()
    {
        if (!std::all_of(ready_.begin(), ready_.end(), [](bool value) { return value; }))
            throw std::runtime_error("incomplete narrow selector file cache");
        complete_ = true;
    }
    std::size_t size(bool tail) const override
    {
        return complete_ ? (tail ? Schedule::tail_packs : Schedule::blocks) * Schedule::block_width : 0;
    }
    void load(bool tail, std::size_t group, std::vector<Selector> &result) const override
    { loadChunk(tail, group, 0, Schedule::block_width, result); }
    bool hasRawChunks() const override { return true; }
    bool supportsRawChunkPrefetch() const override { return true; }
    void loadChunk(bool tail, std::size_t group, std::size_t first, std::size_t count,
                   std::vector<Selector> &result) const override
    {
        index(tail, group);
        if (first > Schedule::block_width || count > Schedule::block_width - first)
            throw std::out_of_range("raw selector chunk range");
        if (!complete_) throw std::runtime_error("incomplete narrow selector file cache");
        std::ifstream file(path(tail, group), std::ios::binary);
        if (!file) throw std::runtime_error("cannot read narrow selector group file");
        cereal::PortableBinaryInputArchive archive(file);
        std::array<std::uint64_t, 9> metadata;
        archive(metadata);
        if (metadata != header(tail, group) || std::filesystem::file_size(path(tail, group)) !=
            1 + sizeof(metadata) + std::size_t{Schedule::block_width} * sizeof(Selector))
            throw std::runtime_error("narrow selector group metadata mismatch");
        result.clear();
        result.resize(count);
        file.seekg(1 + sizeof(metadata) + first * sizeof(Selector));
        if (!file) throw std::runtime_error("cannot seek raw selector chunk");
        if constexpr (is_multilimb_uint_v<typename P::T>) {
            archive(cereal::binary_data(reinterpret_cast<std::uint64_t *>(result.data()),
                                        result.size() * sizeof(Selector)));
        }
        else for (auto &selector : result) archive(selector);
        if (first + count == Schedule::block_width && file.peek() != std::char_traits<char>::eof())
            throw std::runtime_error("unexpected trailing narrow selector file data");
    }
};

// Export existing public ciphertexts without re-encryption. The source is
// committed only after all groups have been written; raw vector storage and
// its NTT caches are then released. External files are absent from key archives.
template <class Schedule>
inline void CKKSNarrowStoreSelectors(CKKSNarrowBootstrapKey<Schedule> &key,
                                    const std::filesystem::path &directory)
{
    if (key.selector_source || key.head.size() != Schedule::blocks * Schedule::block_width ||
        key.tail.size() != Schedule::tail_packs * Schedule::block_width)
        throw std::invalid_argument("invalid expanded narrow selector key");
    auto storage = std::make_shared<CKKSNarrowSelectorFile<Schedule>>(directory);
    for (bool tail : {false, true}) {
        const auto &selectors = tail ? key.tail : key.head;
        for (std::size_t group = 0; group < (tail ? Schedule::tail_packs : Schedule::blocks); group++)
            storage->write(tail, group, std::span<const typename Schedule::Selector>(selectors).subspan(
                group * Schedule::block_width, Schedule::block_width));
    }
    storage->finish();
    key.selector_source = std::move(storage);
    decltype(key.head){}.swap(key.head);
    decltype(key.tail){}.swap(key.tail);
    decltype(key.head_ntt){}.swap(key.head_ntt);
    decltype(key.tail_ntt){}.swap(key.tail_ntt);
}

// Accumulate public packed-extraction terms in fixed-size disk records. Each
// record has its own portable-archive header, so arbitrary coefficient order
// requires only one selector in RAM. Finalization uses one complete group.
// The new source is committed only after extraction and all writes succeed.
template <class Schedule, class PackedSchedule = Schedule>
inline std::shared_ptr<CKKSNarrowSelectorFile<Schedule>> CKKSNarrowReconstructSelectorFile(
    const CKKSNarrowPackedKey<PackedSchedule> &packed,
    const std::filesystem::path &directory,
    const std::function<void(std::size_t)> &term_observer = {},
    const std::filesystem::path &scratch_directory = {})
{
    using P = typename Schedule::Param;
    using Selector = typename Schedule::Selector;
    auto storage = std::make_shared<CKKSNarrowSelectorFile<Schedule>>(directory);
    // A caller may put the temporary random-access accumulator on local disk.
    // Claim a new scratch directory and remove only its known file and directory.
    if (!scratch_directory.empty() && !std::filesystem::create_directory(scratch_directory))
        throw std::runtime_error("packed selector scratch directory must not exist");
    const auto accumulator_path = (scratch_directory.empty() ? directory : scratch_directory) /
        "packed-accumulator.bin";
    struct Cleanup {
        std::filesystem::path path, directory;
        ~Cleanup() {
            std::error_code error;
            std::filesystem::remove(path, error);
            if (!directory.empty()) std::filesystem::remove(directory, error);
        }
    } cleanup{accumulator_path, scratch_directory};
    std::fstream accumulator(accumulator_path,
        std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
    if (!accumulator) throw std::runtime_error("cannot create packed selector accumulator");
    auto transfer = [](auto &archive, auto &selector) {
        if constexpr (is_multilimb_uint_v<typename P::T>) {
            static_assert(std::is_trivially_copyable_v<Selector> &&
                          std::is_standard_layout_v<Selector>);
            static_assert(sizeof(Selector) == (P::k + 1) * P::n * sizeof(typename P::T));
            archive(cereal::binary_data(reinterpret_cast<std::uint64_t *>(&selector), sizeof(Selector)));
        }
        else archive(selector);
    };
    auto value = std::make_unique<Selector>();
    std::stringstream sample;
    { cereal::PortableBinaryOutputArchive archive(sample); transfer(archive, *value); }
    const auto record_bytes = static_cast<std::streamoff>(sample.tellp());
    constexpr std::size_t head_count = Schedule::blocks * Schedule::block_width;
    constexpr std::size_t count = head_count + Schedule::tail_packs * Schedule::block_width;
    std::vector<bool> initialized(count, false);
    auto read = [&](std::size_t index, Selector &selector) {
        if (!initialized.at(index)) {
            for (auto &component : selector.ct) component.fill(typename P::T{});
            return;
        }
        accumulator.seekg(static_cast<std::streamoff>(index) * record_bytes);
        if (!accumulator) throw std::runtime_error("cannot seek packed selector accumulator");
        cereal::PortableBinaryInputArchive archive(accumulator);
        transfer(archive, selector);
    };
    std::size_t terms_written = 0;
    CKKSNarrowExtractSelectorTerms<Schedule>(packed,
        [&](bool tail, std::size_t index, const Selector &term) {
            const auto flat = (tail ? head_count : 0) + index;
            if (index >= (tail ? count - head_count : head_count))
                throw std::out_of_range("packed selector accumulator index");
            read(flat, *value);
            CKKSAddInPlace<P>(*value, term);
            const auto offset = static_cast<std::streamoff>(flat) * record_bytes;
            accumulator.seekp(offset);
            { cereal::PortableBinaryOutputArchive archive(accumulator); transfer(archive, *value); }
            if (!accumulator || accumulator.tellp() != offset + record_bytes)
                throw std::runtime_error("cannot write fixed-size packed selector record");
            initialized[flat] = true;
            if (term_observer) term_observer(++terms_written);
        });
    std::vector<Selector> group(Schedule::block_width);
    for (bool tail : {false, true})
        for (std::size_t index = 0; index < (tail ? Schedule::tail_packs : Schedule::blocks); index++) {
            for (std::size_t y = 0; y < Schedule::block_width; y++)
                read((tail ? head_count : 0) + index * Schedule::block_width + y, group[y]);
            storage->write(tail, index, group);
        }
    storage->finish();
    accumulator.close();
    if (!accumulator) throw std::runtime_error("cannot close packed selector accumulator");
    // Remove the temporary before publishing the completed selector source.
    std::filesystem::remove(accumulator_path);
    return storage;
}

// Attach a completed reconstruction to an existing compatible refresh key.
template <class Schedule>
inline void CKKSNarrowReconstructSelectorsToFile(
    CKKSNarrowBootstrapKey<Schedule> &key,
    const CKKSNarrowPackedKey<Schedule> &packed,
    const std::filesystem::path &directory)
{
    key.selector_source = CKKSNarrowReconstructSelectorFile<Schedule>(packed, directory);
    decltype(key.head){}.swap(key.head);
    decltype(key.tail){}.swap(key.tail);
    decltype(key.head_ntt){}.swap(key.head_ntt);
    decltype(key.tail_ntt){}.swap(key.tail_ntt);
}

// Setup-side generation needs only one selector group in RAM. Reusable
// automorphism/relinearization setup must already have been generated.
template <class Schedule>
inline void CKKSNarrowEncryptSelectorsToFile(CKKSNarrowBootstrapKey<Schedule> &key,
    const Key<typename Schedule::Param> &source,
    const Key<typename Schedule::Param> &evaluation,
    const std::filesystem::path &directory,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    const auto positions = narrow_detail::structuredPositions<Schedule>(source);
    auto storage = std::make_shared<CKKSNarrowSelectorFile<Schedule>>(directory);
    std::vector<typename Schedule::Selector> group(Schedule::block_width);
    auto indicator = [&](std::size_t x, std::size_t y) -> std::complex<double> {
        if (y >= Schedule::block_size(x) || positions[x] != static_cast<int>(y)) return {};
        return source[x * Schedule::block_width + y] == typename P::T{1}
            ? std::complex<double>{1, 0} : std::complex<double>{0, 1};
    };
    for (bool tail : {false, true})
        for (std::size_t index = 0; index < (tail ? Schedule::tail_packs : Schedule::blocks); index++) {
#pragma omp parallel for schedule(static)
            for (std::size_t y = 0; y < Schedule::block_width; y++) {
                CKKSSlotVector<P> mask{};
                if (!tail) std::fill_n(mask.begin(), Schedule::head_slots, indicator(index, y));
                else for (std::size_t ell = 0; ell < Schedule::tail_blocks_per_pack; ell++) {
                    const auto x = index * Schedule::tail_blocks_per_pack + ell;
                    if (x == Schedule::blocks) break;
                    std::fill_n(mask.begin() + (2 * ell + 1) * Schedule::tail_width,
                                Schedule::tail_width, indicator(x, y));
                }
                ckksSlotEncrypt<P>(group[y], mask, evaluation, noise);
            }
            storage->write(tail, index, group);
        }
    storage->finish();
    key.selector_source = std::move(storage);
    decltype(key.head){}.swap(key.head);
    decltype(key.tail){}.swap(key.tail);
    decltype(key.head_ntt){}.swap(key.head_ntt);
    decltype(key.tail_ntt){}.swap(key.tail_ntt);
}

} // namespace TFHEpp
#endif
