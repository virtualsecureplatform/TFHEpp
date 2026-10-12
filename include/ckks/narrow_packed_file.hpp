#pragma once

#include "narrow_bootstrapping.hpp"
#include <cereal/archives/portable_binary.hpp>
#include <filesystem>
#include <fstream>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

template <class Schedule>
constexpr std::array<std::uint64_t, 9> CKKSNarrowPackedFileMetadata()
{
    using P = typename Schedule::Param;
    return {0x5041434b43463031ULL, P::n, Schedule::selector_log_q,
            Schedule::log_delta, Schedule::plain_log_delta, Schedule::blocks,
            std::numeric_limits<typename P::T>::digits, P::Bgbit, P::B̅gbit};
}

// Public encrypted packed indicators and substitution keys only. Derived FFT
// caches are rebuilt after loading; no secret key is stored in this archive.
template <class Schedule>
inline void CKKSNarrowSavePackedKey(const CKKSNarrowPackedKey<Schedule> &key,
                                   const std::filesystem::path &path)
{
    const auto pending = std::filesystem::path(path.string() + ".pending");
    if (std::filesystem::exists(path)) throw std::runtime_error("packed archive already exists");
    if (!std::filesystem::create_directory(pending))
        throw std::runtime_error("packed archive staging directory already exists");
    struct Cleanup {
        std::filesystem::path path;
        ~Cleanup() { std::error_code error; std::filesystem::remove_all(path, error); }
    } cleanup{pending};
    std::filesystem::permissions(pending, std::filesystem::perms::owner_all,
        std::filesystem::perm_options::replace);
    const auto payload = pending / "payload.bin";
    std::ofstream file(payload, std::ios::binary);
    if (!file) throw std::runtime_error("cannot create packed archive");
    std::filesystem::permissions(payload,
        std::filesystem::perms::owner_read | std::filesystem::perms::owner_write,
        std::filesystem::perm_options::replace);
    {
        cereal::PortableBinaryOutputArchive archive(file);
        archive(CKKSNarrowPackedFileMetadata<Schedule>(), key);
    }
    file.flush();
    if (!file) throw std::runtime_error("cannot flush packed archive");
    file.close();
    if (!file) throw std::runtime_error("cannot close packed archive");
    // Publish without replacing an archive created by another writer.
    std::filesystem::create_hard_link(payload, path);
}

template <class Schedule>
inline std::unique_ptr<CKKSNarrowPackedKey<Schedule>> CKKSNarrowLoadPackedKey(
    const std::filesystem::path &path)
{
    using Packed = CKKSNarrowPackedKey<Schedule>;
    // All archived fields are fixed-size arrays, with one endian marker and
    // nine uint64 profile words. Reject truncated/extended input before allocation.
    constexpr std::uintmax_t bytes = 1 + 9 * sizeof(std::uint64_t) +
        sizeof(decltype(Packed::indicators)) + sizeof(decltype(Packed::substitutions));
    if (std::filesystem::file_size(path) != bytes)
        throw std::runtime_error("packed archive length mismatch");
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open packed archive");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t, 9> actual{};
    archive(actual);
    if (actual != CKKSNarrowPackedFileMetadata<Schedule>())
        throw std::runtime_error("packed archive profile mismatch");
    auto result = std::make_unique<Packed>();
    archive(*result);
    if (file.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("packed archive has trailing data");
    return result;
}

} // namespace TFHEpp
#endif
