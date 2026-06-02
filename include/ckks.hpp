#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <tuple>
#include <utility>
#include <vector>

#include "bfv++.hpp"
#include "bfv-slots.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
struct CKKSCiphertext {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);
    static_assert(LogDelta < LogQ);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t log_budget = LogQ - LogDelta;

    TRLWE<P> ct{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(ct);
    }
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
struct CKKSPlaintext {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);
    static_assert(LogDelta < LogQ);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;

    Polynomial<P> poly{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(poly);
    }
};

template <class P, std::uint32_t LogQ>
constexpr std::uint32_t CKKSKeySwitchRowCountForLevel()
{
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);
    constexpr std::uint32_t rows =
        (LogQ + P::B̅gbit - 1) / P::B̅gbit;
    static_assert(rows > 0);
    static_assert(rows <= P::l̅);
    return rows;
}

template <class P, std::uint32_t LogQ>
constexpr std::uint32_t CKKSKeySwitchFirstRowForLevel()
{
    return P::l̅ - CKKSKeySwitchRowCountForLevel<P, LogQ>();
}

template <class P, std::uint32_t LogQ>
struct CKKSRelinKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t key_switch_rows =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    using Rows = std::array<TRLWE<P>, key_switch_rows>;

    Rows data{};

    TRLWE<P> &operator[](std::size_t i) { return data[i]; }
    const TRLWE<P> &operator[](std::size_t i) const { return data[i]; }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(data);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSAutoKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t key_switch_rows =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    using Rows = std::array<TRLWE<P>, key_switch_rows>;

    std::array<Rows, P::k> data{};

    Rows &operator[](std::size_t i) { return data[i]; }
    const Rows &operator[](std::size_t i) const
    {
        return data[i];
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(data);
    }
};

template <class P, std::uint32_t LogQ>
using CKKSGaloisKey = std::array<CKKSAutoKey<P, LogQ>, P::nbit + 1>;

using CKKSSeed = std::array<std::uint64_t, 4>;

template <class P, std::uint32_t LogQ>
struct CKKSSeededKeySwitchRow {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    std::array<CKKSSeed, P::k> mask_seeds{};
    Polynomial<P> body{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(mask_seeds, body);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSSeededRelinKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t key_switch_rows =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    using Rows = std::array<CKKSSeededKeySwitchRow<P, LogQ>, key_switch_rows>;

    Rows data{};

    CKKSSeededKeySwitchRow<P, LogQ> &operator[](std::size_t i)
    {
        return data[i];
    }
    const CKKSSeededKeySwitchRow<P, LogQ> &operator[](std::size_t i) const
    {
        return data[i];
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(data);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSSeededAutoKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t key_switch_rows =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    using Rows = std::array<CKKSSeededKeySwitchRow<P, LogQ>, key_switch_rows>;

    std::array<Rows, P::k> data{};

    Rows &operator[](std::size_t i) { return data[i]; }
    const Rows &operator[](std::size_t i) const { return data[i]; }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(data);
    }
};

template <class P, std::uint32_t LogQ>
using CKKSSeededSecretKeySwitchKey = CKKSSeededAutoKey<P, LogQ>;

template <class P, std::uint32_t LogQ>
using CKKSSeededGaloisKey =
    std::array<CKKSSeededAutoKey<P, LogQ>, P::nbit + 1>;

namespace ckks_detail {

template <class T>
struct is_seeded_key_switch_row_array : std::false_type {};

template <class P, std::uint32_t LogQ, std::size_t N>
struct is_seeded_key_switch_row_array<
    std::array<CKKSSeededKeySwitchRow<P, LogQ>, N>> : std::true_type {};

template <class T>
inline constexpr bool is_seeded_key_switch_row_array_v =
    is_seeded_key_switch_row_array<std::remove_cvref_t<T>>::value;

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ>
using CKKSSecretKeySwitchKey = CKKSAutoKey<P, LogQ>;

template <class P>
constexpr std::size_t CKKSRelinKeySwitchRowCount()
{
    return P::l̅;
}

template <class P, std::uint32_t LogQ>
constexpr std::size_t CKKSRelinKeySwitchRowCount()
{
    return CKKSKeySwitchRowCountForLevel<P, LogQ>();
}

template <class P>
constexpr std::size_t CKKSAutoKeySwitchRowCount()
{
    return P::k * P::l̅;
}

template <class P, std::uint32_t LogQ>
constexpr std::size_t CKKSAutoKeySwitchRowCount()
{
    return P::k * CKKSKeySwitchRowCountForLevel<P, LogQ>();
}

template <class P>
constexpr std::size_t CKKSKeySwitchRowByteSize()
{
    return sizeof(TRLWE<P>);
}

template <class P>
constexpr std::size_t CKKSSeededKeySwitchRowByteSize()
{
    return sizeof(Polynomial<P>) + P::k * sizeof(CKKSSeed);
}

template <class P, std::uint32_t LogQ>
constexpr std::size_t CKKSRelinKeyByteEstimate()
{
    return CKKSRelinKeySwitchRowCount<P, LogQ>() *
           CKKSKeySwitchRowByteSize<P>();
}

template <class P, std::uint32_t LogQ>
constexpr std::size_t CKKSSeededRelinKeyByteEstimate()
{
    return CKKSRelinKeySwitchRowCount<P, LogQ>() *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class P, std::uint32_t LogQ>
constexpr std::size_t CKKSAutoKeyByteEstimate()
{
    return CKKSAutoKeySwitchRowCount<P, LogQ>() *
           CKKSKeySwitchRowByteSize<P>();
}

template <class P, std::uint32_t LogQ>
constexpr std::size_t CKKSSeededAutoKeyByteEstimate()
{
    return CKKSAutoKeySwitchRowCount<P, LogQ>() *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class P>
using CKKSRotationKeyIndexSet = std::array<bool, P::nbit + 1>;

template <class P>
using CKKSDirectRotationKeyIndexSet = std::array<bool, P::n / 2>;

template <class P>
inline void CKKSClearRotationKeyIndexSet(CKKSRotationKeyIndexSet<P> &keys)
{
    keys.fill(false);
}

template <class P>
inline void CKKSClearDirectRotationKeyIndexSet(
    CKKSDirectRotationKeyIndexSet<P> &keys)
{
    keys.fill(false);
}

template <class P>
inline void CKKSMarkRotationPowerKeyIndices(CKKSRotationKeyIndexSet<P> &keys,
                                            int steps)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    for (int i = 0; i < static_cast<int>(P::nbit) - 1; i++)
        if ((steps >> i) & 1) keys[static_cast<std::size_t>(i)] = true;
}

template <class P>
inline void CKKSMarkDirectRotationKeyIndex(
    CKKSDirectRotationKeyIndexSet<P> &keys, int steps)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    if (steps != 0) keys[static_cast<std::size_t>(steps)] = true;
}

template <class P>
inline void CKKSMarkConjugationKeyIndex(CKKSRotationKeyIndexSet<P> &keys)
{
    keys[P::nbit] = true;
}

template <class P>
inline std::size_t
CKKSRotationKeyIndexSetCount(const CKKSRotationKeyIndexSet<P> &keys)
{
    return static_cast<std::size_t>(
        std::count(keys.begin(), keys.end(), true));
}

template <class P>
inline std::size_t CKKSDirectRotationKeyIndexSetCount(
    const CKKSDirectRotationKeyIndexSet<P> &keys)
{
    return static_cast<std::size_t>(
        std::count(keys.begin(), keys.end(), true));
}

template <class P, std::uint32_t LogQ>
struct CKKSSparseGaloisKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    CKKSRotationKeyIndexSet<P> available{};
    std::array<std::unique_ptr<CKKSAutoKey<P, LogQ>>, P::nbit + 1> keys{};

    bool has(std::size_t i) const
    {
        return i < keys.size() && available[i] && keys[i] != nullptr;
    }

    CKKSAutoKey<P, LogQ> &get(std::size_t i)
    {
        assert(has(i));
        return *keys[i];
    }

    const CKKSAutoKey<P, LogQ> &get(std::size_t i) const
    {
        assert(has(i));
        return *keys[i];
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(available, keys);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSDirectSparseGaloisKeyEntry {
    int steps = 0;
    CKKSAutoKey<P, LogQ> key{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(steps, key);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSDirectSparseGaloisKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    std::vector<CKKSDirectSparseGaloisKeyEntry<P, LogQ>> keys{};

    static constexpr int normalize_steps(int steps)
    {
        constexpr int half = static_cast<int>(P::n) / 2;
        return ((steps % half) + half) % half;
    }

    bool has(int steps) const
    {
        steps = normalize_steps(steps);
        if (steps == 0) return true;
        return std::any_of(keys.begin(), keys.end(), [&](const auto &entry) {
            return entry.steps == steps;
        });
    }

    CKKSAutoKey<P, LogQ> &get(int steps)
    {
        steps = normalize_steps(steps);
        assert(steps != 0);
        auto it = std::find_if(keys.begin(), keys.end(), [&](auto &entry) {
            return entry.steps == steps;
        });
        assert(it != keys.end());
        return it->key;
    }

    const CKKSAutoKey<P, LogQ> &get(int steps) const
    {
        steps = normalize_steps(steps);
        assert(steps != 0);
        auto it =
            std::find_if(keys.begin(), keys.end(), [&](const auto &entry) {
                return entry.steps == steps;
            });
        assert(it != keys.end());
        return it->key;
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSHybridSparseGaloisKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    CKKSSparseGaloisKey<P, LogQ> binary{};
    CKKSDirectSparseGaloisKey<P, LogQ> direct{};

    bool has(int steps) const
    {
        return direct.has(steps) || binary.has(steps);
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(binary, direct);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSSeededSparseGaloisKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    CKKSRotationKeyIndexSet<P> available{};
    std::array<std::unique_ptr<CKKSSeededAutoKey<P, LogQ>>, P::nbit + 1>
        keys{};

    bool has(std::size_t i) const
    {
        return i < keys.size() && available[i] && keys[i] != nullptr;
    }

    CKKSSeededAutoKey<P, LogQ> &get(std::size_t i)
    {
        assert(has(i));
        return *keys[i];
    }

    const CKKSSeededAutoKey<P, LogQ> &get(std::size_t i) const
    {
        assert(has(i));
        return *keys[i];
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(available, keys);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSSeededDirectSparseGaloisKeyEntry {
    int steps = 0;
    CKKSSeededAutoKey<P, LogQ> key{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(steps, key);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSSeededDirectSparseGaloisKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    std::vector<CKKSSeededDirectSparseGaloisKeyEntry<P, LogQ>> keys{};

    static constexpr int normalize_steps(int steps)
    {
        constexpr int half = static_cast<int>(P::n) / 2;
        return ((steps % half) + half) % half;
    }

    bool has(int steps) const
    {
        steps = normalize_steps(steps);
        if (steps == 0) return true;
        return std::any_of(keys.begin(), keys.end(), [&](const auto &entry) {
            return entry.steps == steps;
        });
    }

    CKKSSeededAutoKey<P, LogQ> &get(int steps)
    {
        steps = normalize_steps(steps);
        assert(steps != 0);
        auto it = std::find_if(keys.begin(), keys.end(), [&](auto &entry) {
            return entry.steps == steps;
        });
        assert(it != keys.end());
        return it->key;
    }

    const CKKSSeededAutoKey<P, LogQ> &get(int steps) const
    {
        steps = normalize_steps(steps);
        assert(steps != 0);
        auto it =
            std::find_if(keys.begin(), keys.end(), [&](const auto &entry) {
                return entry.steps == steps;
            });
        assert(it != keys.end());
        return it->key;
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSSeededHybridSparseGaloisKey {
    static_assert(std::is_same_v<typename P::T, __uint128_t> ||
                  is_multilimb_uint_v<typename P::T>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);

    static constexpr std::uint32_t log_q = LogQ;

    CKKSSeededSparseGaloisKey<P, LogQ> binary{};
    CKKSSeededDirectSparseGaloisKey<P, LogQ> direct{};

    bool has(int steps) const { return direct.has(steps) || binary.has(steps); }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(binary, direct);
    }
};

struct CKKSNoise {
    double modular_stdev = 0.0;
    std::uint32_t uniform_bits = 0;
};

template <class T>
inline void CKKSSavePortableBinary(const std::filesystem::path &path,
                                   const T &value)
{
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) throw std::runtime_error("failed to open CKKS key file for write");
    cereal::PortableBinaryOutputArchive archive(ofs);
    archive(value);
}

template <class T>
inline void CKKSSavePortableBinaryAtomic(const std::filesystem::path &path,
                                         const T &value)
{
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    std::filesystem::path tmp = path;
    tmp += ".tmp";

    try {
        {
            std::ofstream ofs(tmp, std::ios::binary | std::ios::trunc);
            if (!ofs) {
                throw std::runtime_error(
                    "failed to open temporary CKKS key file for write");
            }
            cereal::PortableBinaryOutputArchive archive(ofs);
            archive(value);
            ofs.flush();
            if (!ofs) {
                throw std::runtime_error(
                    "failed to flush temporary CKKS key file");
            }
        }
        std::filesystem::rename(tmp, path);
    }
    catch (...) {
        std::error_code ec;
        std::filesystem::remove(tmp, ec);
        throw;
    }
}

template <class T>
inline void CKKSLoadPortableBinary(T &value, const std::filesystem::path &path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("failed to open CKKS key file for read: " +
                                 path.string());
    }
    cereal::PortableBinaryInputArchive archive(ifs);
    archive(value);
}

inline void CKKSRequireKeyFileReadable(const std::filesystem::path &path,
                                       const char *label)
{
    std::error_code ec;
    const bool exists = std::filesystem::exists(path, ec);
    if (ec || !exists) {
        throw std::runtime_error(std::string(label) +
                                 " CKKS key file missing: " + path.string());
    }
    const bool regular = std::filesystem::is_regular_file(path, ec);
    if (ec || !regular) {
        throw std::runtime_error(std::string(label) +
                                 " CKKS key file is not a regular file: " +
                                 path.string());
    }
}

struct CKKSDenseBootstrapKeyDirectoryOptions {
    bool overwrite_existing = true;
};

struct CKKSDenseBootstrapKeyFileByteEstimate {
    std::filesystem::path path{};
    std::size_t estimated_bytes = 0;
};

struct CKKSDenseBootstrapTimings {
    double normalize_ms = 0.0;
    double input_secret_switch_ms = 0.0;
    double modraise_ms = 0.0;
    double coeff_to_slot_ms = 0.0;
    double split_ms = 0.0;
    double real_evalmod_ms = 0.0;
    double imag_evalmod_ms = 0.0;
    double slot_to_coeff_ms = 0.0;
    double output_secret_switch_ms = 0.0;

    double total_ms() const
    {
        return normalize_ms + input_secret_switch_ms + modraise_ms +
               coeff_to_slot_ms + split_ms + real_evalmod_ms +
               imag_evalmod_ms + slot_to_coeff_ms + output_secret_switch_ms;
    }
};

struct CKKSDenseBootstrapProductTimings {
    double multiply_ms = 0.0;
    CKKSDenseBootstrapTimings bootstrap{};

    double total_ms() const { return multiply_ms + bootstrap.total_ms(); }
};

struct CKKSDenseBootstrapProgress {
    void (*stage_begin)(const char *stage, const void *context) = nullptr;
    void (*stage_end)(const char *stage, double elapsed_ms,
                      const void *context) = nullptr;
    const void *context = nullptr;
};

struct CKKSDenseBootstrapKeyDirectoryManifest {
    static constexpr std::uint32_t current_version = 3;
    static constexpr std::uint32_t sparse_format = 0;
    static constexpr std::uint32_t hybrid_giant_format = 1;
    static constexpr std::uint32_t seeded_hybrid_giant_format = 2;
    static constexpr std::uint32_t hybrid_giant_tuned_format_base = 1000;
    static constexpr std::uint32_t seeded_hybrid_giant_tuned_format_base = 2000;
    static constexpr std::uint32_t seeded_hybrid_giant_streamed_format = 3;
    static constexpr std::uint32_t seeded_hybrid_giant_streamed_tuned_format_base =
        3000;

    std::uint32_t version = current_version;
    std::uint32_t key_format = sparse_format;
    std::uint32_t n = 0;
    std::uint32_t nbit = 0;
    std::uint32_t torus_bits = 0;
    std::uint32_t log_delta = 0;
    std::uint32_t log_message_ratio = 0;
    std::uint32_t input_log_q = 0;
    std::uint32_t boot_log_q = 0;
    std::uint32_t output_log_q = 0;
    std::uint32_t linear_plain_log_delta = 0;
    std::uint32_t linear_fuse_radix = 0;
    std::uint32_t coeff_to_slot_plain_log_delta = 0;
    std::uint32_t component_split_plain_log_delta = 0;
    std::uint32_t slot_to_coeff_plain_log_delta = 0;
    std::uint32_t coeff_to_slot_fuse_radix = 0;
    std::uint32_t slot_to_coeff_fuse_radix = 0;
    std::int32_t linear_bsgs_step = 0;
    std::uint32_t coeff_to_slot_level_count = 0;
    std::uint32_t slot_to_coeff_level_count = 0;
    std::uint32_t evalmod_degree = 0;
    std::uint32_t evalmod_k = 0;
    std::uint32_t evalmod_double_angle = 0;
    std::uint32_t evalmod_inv_degree = 0;
    std::uint32_t evalmod_log_scale = 0;
    std::uint32_t evalmod_depth = 0;
    std::uint32_t modraise_mask_bound = 0;
    std::uint64_t expected_file_count = 0;
    std::uint64_t sparse_key_rows = 0;
    std::uint64_t streamed_peak_key_rows = 0;
    std::uint64_t hybrid_key_rows = 0;
    std::uint64_t hybrid_streamed_peak_key_rows = 0;
    std::uint64_t full_key_rows = 0;

    bool operator==(const CKKSDenseBootstrapKeyDirectoryManifest &) const =
        default;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(version, key_format, n, nbit, torus_bits, log_delta,
                log_message_ratio, input_log_q, boot_log_q, output_log_q,
                linear_plain_log_delta, linear_fuse_radix,
                coeff_to_slot_plain_log_delta,
                component_split_plain_log_delta,
                slot_to_coeff_plain_log_delta, coeff_to_slot_fuse_radix,
                slot_to_coeff_fuse_radix, linear_bsgs_step,
                coeff_to_slot_level_count, slot_to_coeff_level_count,
                evalmod_degree, evalmod_k, evalmod_double_angle,
                evalmod_inv_degree, evalmod_log_scale, evalmod_depth,
                modraise_mask_bound, expected_file_count, sparse_key_rows,
                streamed_peak_key_rows, hybrid_key_rows,
                hybrid_streamed_peak_key_rows, full_key_rows);
    }
};

namespace ckks_detail {

inline bool ckks_linear_trace_enabled()
{
    static const bool enabled =
        std::getenv("TFHEPP_CKKS_TRACE_LINEAR") != nullptr;
    return enabled;
}

inline bool ckks_evalmod_trace_enabled()
{
    static const bool enabled =
        std::getenv("TFHEPP_CKKS_TRACE_EVALMOD") != nullptr;
    return enabled;
}

inline double ckks_trace_elapsed_ms(
    const std::chrono::steady_clock::time_point &start)
{
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - start)
        .count();
}

inline void ckks_trace_linear_event(const char *name, double elapsed_ms)
{
    if (!ckks_linear_trace_enabled()) return;
    std::cerr << "ckks_linear_trace name=" << name << " ms=" << elapsed_ms
              << '\n';
}

inline void ckks_trace_linear_group_event(const char *name, std::size_t group,
                                          double elapsed_ms)
{
    if (!ckks_linear_trace_enabled()) return;
    std::cerr << "ckks_linear_trace name=" << name << " group=" << group
              << " ms=" << elapsed_ms << '\n';
}

inline void ckks_trace_evalmod_event(const char *name, double elapsed_ms)
{
    if (!ckks_evalmod_trace_enabled()) return;
    std::cerr << "ckks_evalmod_trace name=" << name << " ms=" << elapsed_ms
              << '\n';
}

inline void ckks_trace_evalmod_power_event(const char *name, std::size_t power,
                                           std::uint32_t depth,
                                           double elapsed_ms)
{
    if (!ckks_evalmod_trace_enabled()) return;
    std::cerr << "ckks_evalmod_trace name=" << name << " power=" << power
              << " depth=" << depth << " ms=" << elapsed_ms << '\n';
}

inline void ckks_trace_evalmod_mult_event(const char *name,
                                          std::uint32_t lhs_log_q,
                                          std::uint32_t rhs_log_q,
                                          std::uint32_t out_log_q,
                                          double elapsed_ms)
{
    if (!ckks_evalmod_trace_enabled()) return;
    std::cerr << "ckks_evalmod_trace name=" << name
              << " lhs_logQ=" << lhs_log_q << " rhs_logQ=" << rhs_log_q
              << " out_logQ=" << out_log_q << " ms=" << elapsed_ms << '\n';
}

template <std::uint32_t A, std::uint32_t B>
inline constexpr std::uint32_t static_min_v = A < B ? A : B;

template <std::uint32_t A, std::uint32_t B>
inline constexpr std::uint32_t static_max_v = A < B ? B : A;

constexpr std::uint32_t ceil_div(std::uint32_t num, std::uint32_t den)
{
    return (num + den - 1) / den;
}

constexpr std::uint32_t bit_width_u64(std::uint64_t value)
{
    std::uint32_t width = 0;
    while (value != 0) {
        width++;
        value >>= 1;
    }
    return width;
}

constexpr double exp2_double(std::uint32_t exponent)
{
    double value = 1.0;
    for (std::uint32_t i = 0; i < exponent; i++) value *= 2.0;
    return value;
}

template <class P>
inline constexpr std::uint32_t torus_width_v =
    std::numeric_limits<typename P::T>::digits;

template <class P>
inline constexpr bool supported_torus_v =
    std::is_same_v<typename P::T, __uint128_t> ||
    is_multilimb_uint_v<typename P::T>;

inline std::uint64_t splitmix64Scramble(std::uint64_t x)
{
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

struct CKKSSeedStream {
    std::uint64_t state = UINT64_C(0x6a09e667f3bcc909);

    explicit CKKSSeedStream(const CKKSSeed &seed)
    {
        for (std::uint64_t limb : seed)
            state = splitmix64Scramble(state ^ limb);
    }

    std::uint64_t next()
    {
        state += UINT64_C(0x9e3779b97f4a7c15);
        return splitmix64Scramble(state);
    }
};

inline CKKSSeed randomSeed()
{
    std::uniform_int_distribution<std::uint64_t> dist64(
        0, std::numeric_limits<std::uint64_t>::max());
    CKKSSeed seed{};
    for (auto &word : seed) word = dist64(generator);
    return seed;
}

template <class P>
inline __int128_t unsignedToI128(typename P::T value)
{
    using T = typename P::T;
    if constexpr (std::is_same_v<T, __uint128_t>) {
        return static_cast<__int128_t>(value);
    }
    else {
        static_assert(is_multilimb_uint_v<T>);
        __uint128_t low = value.limb[0];
        if constexpr (T::limbs > 1)
            low |= static_cast<__uint128_t>(value.limb[1]) << 64;
        return static_cast<__int128_t>(low);
    }
}

template <class P>
inline long double unsignedToLongDouble(typename P::T value)
{
    using T = typename P::T;
    if constexpr (std::is_same_v<T, __uint128_t>) {
        constexpr long double limb_scale = 18446744073709551616.0L;
        const auto lo = static_cast<std::uint64_t>(value);
        const auto hi = static_cast<std::uint64_t>(value >> 64);
        return static_cast<long double>(hi) * limb_scale +
               static_cast<long double>(lo);
    }
    else {
        static_assert(is_multilimb_uint_v<T>);
        long double out = 0.0L;
        for (std::size_t i = T::limbs; i-- > 0;)
            out = std::ldexp(out, 64) +
                  static_cast<long double>(value.limb[i]);
        return out;
    }
}

template <class P>
inline __int128_t torusToSigned(typename P::T value)
{
    static_assert(supported_torus_v<P>);
    if constexpr (std::is_same_v<typename P::T, __uint128_t>) {
        return static_cast<__int128_t>(value);
    }
    else {
        static_assert(torus_width_v<P> <= 128,
                      "Use levelToSigned only for <=128-bit multi-limb levels");
        return unsignedToI128<P>(value);
    }
}

template <class P>
inline typename P::T signedToTorus(__int128_t value)
{
    static_assert(supported_torus_v<P>);
    return static_cast<typename P::T>(value);
}

template <class P, std::uint32_t LogQ>
inline constexpr typename P::T levelMask()
{
    static_assert(supported_torus_v<P>);
    static_assert(LogQ > 0);
    static_assert(LogQ <= torus_width_v<P>);
    if constexpr (LogQ == torus_width_v<P>)
        return std::numeric_limits<typename P::T>::max();
    else
        return (typename P::T{1} << LogQ) - typename P::T{1};
}

template <class P, std::uint32_t LogQ>
inline typename P::T reduceToLevel(typename P::T value)
{
    return value & levelMask<P, LogQ>();
}

template <class P, std::uint32_t LogQ>
inline typename P::T signedToLevel(__int128_t value)
{
    return reduceToLevel<P, LogQ>(signedToTorus<P>(value));
}

template <class P, std::uint32_t LogQ>
inline __int128_t levelToSigned(typename P::T value)
{
    value = reduceToLevel<P, LogQ>(value);
    if constexpr (std::is_same_v<typename P::T, __uint128_t> &&
                  LogQ == torus_width_v<P>) {
        return torusToSigned<P>(value);
    }
    else {
        static_assert(LogQ < 128,
                      "levelToSigned returns __int128_t and cannot represent "
                      "larger active CKKS levels");
        const typename P::T sign = typename P::T{1} << (LogQ - 1);
        if ((value & sign) == typename P::T{0}) return unsignedToI128<P>(value);

        const typename P::T modulus = typename P::T{1} << LogQ;
        const typename P::T magnitude = modulus - value;
        return -unsignedToI128<P>(magnitude);
    }
}

template <class P, std::uint32_t LogQ>
inline long double levelToLongDouble(typename P::T value)
{
    using T = typename P::T;
    static_assert(LogQ > 0);
    static_assert(LogQ <= torus_width_v<P>);

    value = reduceToLevel<P, LogQ>(value);
    const T sign = T{1} << (LogQ - 1);
    if ((value & sign) == T{0}) return unsignedToLongDouble<P>(value);

    T magnitude;
    if constexpr (LogQ == torus_width_v<P>)
        magnitude = T{0} - value;
    else
        magnitude = (T{1} << LogQ) - value;
    return -unsignedToLongDouble<P>(magnitude);
}

template <class P, std::uint32_t LogQ>
inline typename P::T centeredLevelToTorus(typename P::T value)
{
    value = reduceToLevel<P, LogQ>(value);
    if constexpr (LogQ == torus_width_v<P>) {
        return value;
    }
    else {
        const typename P::T sign = typename P::T{1} << (LogQ - 1);
        if ((value & sign) == typename P::T{0}) return value;
        return value | ~levelMask<P, LogQ>();
    }
}

template <class P, std::uint32_t LogQ>
inline typename P::T uniformAtLevel()
{
    if constexpr (LogQ == torus_width_v<P>) {
        return UniformTorusRandom<P>();
    }
    else {
        std::uniform_int_distribution<std::uint64_t> dist64(
            0, std::numeric_limits<std::uint64_t>::max());
        typename P::T value = 0;
        std::uint32_t generated = 0;
        while (generated < LogQ) {
            value |= static_cast<typename P::T>(dist64(generator)) << generated;
            generated += 64;
        }
        return reduceToLevel<P, LogQ>(value);
    }
}

template <class P, std::uint32_t LogQ>
inline typename P::T seededUniformAtLevel(CKKSSeedStream &stream)
{
    static_assert(LogQ > 0);
    static_assert(LogQ <= torus_width_v<P>);
    typename P::T value = 0;
    std::uint32_t generated = 0;
    while (generated < LogQ) {
        value |= static_cast<typename P::T>(stream.next()) << generated;
        generated += 64;
    }
    return reduceToLevel<P, LogQ>(value);
}

template <class P, std::uint32_t LogQ>
inline void fillSeededUniformPolynomialAtLevel(Polynomial<P> &poly,
                                               const CKKSSeed &seed)
{
    CKKSSeedStream stream(seed);
    for (std::uint32_t i = 0; i < P::n; i++)
        poly[i] = seededUniformAtLevel<P, LogQ>(stream);
}

inline __int128_t randomSignedBits(std::uint32_t bits)
{
    assert(bits > 0);
    assert(bits < 127);

    std::uniform_int_distribution<std::uint64_t> dist64(
        0, std::numeric_limits<std::uint64_t>::max());
    __uint128_t value = 0;
    std::uint32_t generated = 0;
    while (generated < bits) {
        value |= static_cast<__uint128_t>(dist64(generator)) << generated;
        generated += 64;
    }
    value &= (static_cast<__uint128_t>(1) << bits) - 1;

    const __uint128_t sign = static_cast<__uint128_t>(1) << (bits - 1);
    if ((value & sign) == 0) return static_cast<__int128_t>(value);
    return static_cast<__int128_t>(value) -
           static_cast<__int128_t>(static_cast<__uint128_t>(1) << bits);
}

inline void fftInPlace(std::vector<std::complex<long double>> &a, bool inverse)
{
    const std::size_t n = a.size();
    assert(n != 0 && (n & (n - 1)) == 0);

    for (std::size_t i = 1, j = 0; i < n; i++) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const long double angle =
            (inverse ? 2.0L : -2.0L) * pi / static_cast<long double>(len);
        const std::complex<long double> wlen(std::cos(angle),
                                             std::sin(angle));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<long double> w = 1.0L;
            for (std::size_t j = 0; j < len / 2; j++) {
                const std::complex<long double> u = a[i + j];
                const std::complex<long double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (inverse)
        for (auto &x : a) x /= static_cast<long double>(n);
}

template <class P, std::uint32_t LogQ>
inline typename P::T sampleNoise(const CKKSNoise &noise)
{
    typename P::T value = 0;
    if (noise.modular_stdev != 0.0) {
        std::normal_distribution<long double> distribution(0.0L, 1.0L);
        const long double scaled =
            std::ldexp(static_cast<long double>(noise.modular_stdev) *
                           distribution(generator),
                       LogQ);
        value += signedToLevel<P, LogQ>(longDoubleToI128(scaled));
    }
    if (noise.uniform_bits > 0) {
        assert(noise.uniform_bits < LogQ);
        value += signedToLevel<P, LogQ>(randomSignedBits(noise.uniform_bits));
    }
    return reduceToLevel<P, LogQ>(value);
}

template <class P, std::uint32_t LogQ>
inline void reducePolynomialToLevel(Polynomial<P> &poly)
{
    for (std::uint32_t i = 0; i < P::n; i++)
        poly[i] = reduceToLevel<P, LogQ>(poly[i]);
}

template <class P, std::uint32_t LogQ>
inline void reduceTRLWEToLevel(TRLWE<P> &ct)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        reducePolynomialToLevel<P, LogQ>(ct[c]);
}

template <class P, std::uint32_t LogQ>
inline void reduceTRLWE3ToLevel(TRLWE3<P> &ct)
{
    for (int c = 0; c < 3; c++) reducePolynomialToLevel<P, LogQ>(ct[c]);
}

template <class P>
inline bool fitsU64(typename P::T value)
{
    if constexpr (std::is_same_v<typename P::T, __uint128_t>) {
        return value <= static_cast<__uint128_t>(
                            std::numeric_limits<std::uint64_t>::max());
    }
    else {
        using T = typename P::T;
        static_assert(is_multilimb_uint_v<T>);
        for (std::size_t i = 1; i < T::limbs; i++)
            if (value.limb[i] != 0) return false;
        return true;
    }
}

template <class P>
inline std::uint64_t lowU64(typename P::T value)
{
    if constexpr (std::is_same_v<typename P::T, __uint128_t>)
        return static_cast<std::uint64_t>(value);
    else {
        static_assert(is_multilimb_uint_v<typename P::T>);
        return value.limb[0];
    }
}

template <class P, std::uint32_t LogQ>
inline std::pair<bool, std::uint64_t> smallSignedMagnitude(typename P::T value)
{
    using T = typename P::T;
    static_assert(LogQ > 0);
    static_assert(LogQ <= torus_width_v<P>);

    value = reduceToLevel<P, LogQ>(value);
    const T sign = T{1} << (LogQ - 1);
    const bool negative = (value & sign) != T{0};
    T magnitude = value;
    if (negative) {
        if constexpr (LogQ == torus_width_v<P>)
            magnitude = T{0} - value;
        else
            magnitude = (T{1} << LogQ) - value;
    }

    assert(fitsU64<P>(magnitude));
    return {negative, lowU64<P>(magnitude)};
}

template <class P, std::uint32_t LogQ, std::uint32_t DropBits>
inline typename P::T roundedLevelRightShift(typename P::T value)
{
    using T = typename P::T;
    static_assert(DropBits > 0);
    static_assert(DropBits < LogQ);
    static_assert(LogQ <= torus_width_v<P>);
    constexpr std::uint32_t out_log_q = LogQ - DropBits;

    value = reduceToLevel<P, LogQ>(value);
    const T sign = T{1} << (LogQ - 1);
    const bool negative = (value & sign) != T{0};
    const T half = T{1} << (DropBits - 1);

    if (!negative)
        return reduceToLevel<P, out_log_q>((value + half) >> DropBits);

    T magnitude;
    if constexpr (LogQ == torus_width_v<P>)
        magnitude = T{0} - value;
    else
        magnitude = (T{1} << LogQ) - value;
    const T rounded = (magnitude + half) >> DropBits;
    return reduceToLevel<P, out_log_q>(T{0} - rounded);
}

template <class P, std::uint32_t PlainLogQ>
inline void polyMulTorusBySmallSigned(Polynomial<P> &res,
                                      const Polynomial<P> &torus,
                                      const Polynomial<P> &plain)
{
    using T = typename P::T;
    static_assert(PlainLogQ > 0);
    static_assert(PlainLogQ <= torus_width_v<P>);

    if constexpr (is_multilimb_uint_v<T>) {
        auto pos = std::make_unique<std::array<std::uint64_t, P::n>>();
        auto neg = std::make_unique<std::array<std::uint64_t, P::n>>();
        bool has_pos = false;
        bool has_neg = false;
        for (std::uint32_t i = 0; i < P::n; i++) {
            const auto [negative, magnitude] =
                smallSignedMagnitude<P, PlainLogQ>(plain[i]);
            (*pos)[i] = negative ? 0 : magnitude;
            (*neg)[i] = negative ? magnitude : 0;
            has_pos = has_pos || (!negative && magnitude != 0);
            has_neg = has_neg || (negative && magnitude != 0);
        }

        if (has_pos)
            PolyMulTorusByUnsignedBits<P, 64>(res, torus, *pos);
        else
            for (std::uint32_t i = 0; i < P::n; i++) res[i] = 0;

        if (has_neg) {
            auto tmp = std::make_unique<Polynomial<P>>();
            PolyMulTorusByUnsignedBits<P, 64>(*tmp, torus, *neg);
            for (std::uint32_t i = 0; i < P::n; i++) res[i] -= (*tmp)[i];
        }
    }
    else {
        static_assert(std::is_same_v<T, __uint128_t>);
        for (int i = 0; i < static_cast<int>(P::n); i++) {
            T ri = 0;
            for (int j = 0; j <= i; j++)
                ri += torus[j] * plain[i - j];
            for (int j = i + 1; j < static_cast<int>(P::n); j++)
                ri -= torus[j] * plain[P::n + i - j];
            res[i] = reduceToLevel<P, PlainLogQ>(ri);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void centeredTRLWEAtLevel(TRLWE<P> &out, const TRLWE<P> &in)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++)
            out[c][n] = centeredLevelToTorus<P, LogQ>(in[c][n]);
}

template <class P, std::uint32_t LogQ>
inline void centeredPolynomialAtLevel(Polynomial<P> &out,
                                      const Polynomial<P> &in)
{
    for (std::uint32_t n = 0; n < P::n; n++)
        out[n] = centeredLevelToTorus<P, LogQ>(in[n]);
}

template <class P, std::uint32_t LogScale>
inline typename P::T rescaleAccumulator(Wide384 acc)
{
    static_assert(std::is_same_v<typename P::T, __uint128_t>,
                  "CKKS v1 expects a 128-bit torus parameter set");
    static_assert(LogScale < torus_width_v<P>);
    if constexpr (LogScale > 0) {
        const bool negative = (acc.w[5] >> 63) != 0;
        acc.add_shifted(negative ? -static_cast<__int128_t>(1)
                                 : static_cast<__int128_t>(1),
                        static_cast<int>(LogScale) - 1);
    }
    const typename P::T divisor = LogScale == 0
                                      ? static_cast<typename P::T>(1)
                                      : (static_cast<typename P::T>(1)
                                         << LogScale);
    return acc.div128(divisor);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void encodeRealPolynomial(Polynomial<P> &poly,
                                 const std::array<double, P::n> &coeffs)
{
    static_assert(LogDelta < LogQ);
    for (std::uint32_t i = 0; i < P::n; i++) {
        const long double scaled =
            std::ldexp(static_cast<long double>(coeffs[i]), LogDelta);
        const auto rounded = static_cast<__int128_t>(std::round(scaled));
        poly[i] = signedToLevel<P, LogQ>(rounded);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void decodeRealPolynomial(std::array<double, P::n> &coeffs,
                                 const Polynomial<P> &poly)
{
    static_assert(LogDelta < LogQ);
    for (std::uint32_t i = 0; i < P::n; i++) {
        const long double value =
            std::ldexp(levelToLongDouble<P, LogQ>(poly[i]),
                       -static_cast<int>(LogDelta));
        coeffs[i] = static_cast<double>(value);
    }
}

template <class P>
inline const std::array<std::uint32_t, P::n / 2> &ckksSlotToEvalIndex()
{
    static const auto table = [] {
        std::array<std::uint32_t, P::n / 2> tbl{};
        std::uint64_t g = 1;
        constexpr std::uint64_t twon = 2 * P::n;
        for (std::uint32_t i = 0; i < P::n / 2; i++) {
            tbl[i] = static_cast<std::uint32_t>((g - 1) / 2);
            g = (g * 5) % twon;
        }
        return tbl;
    }();
    return table;
}

inline std::uint32_t bitReverse(std::uint32_t value, std::uint32_t bits)
{
    std::uint32_t reversed = 0;
    for (std::uint32_t i = 0; i < bits; i++) {
        reversed = (reversed << 1) | (value & 1U);
        value >>= 1;
    }
    return reversed;
}

template <std::uint32_t LogDelta, class P>
inline void fillConjugateSymmetricEvaluations(
    std::array<std::complex<double>, P::n> &evals,
    const std::array<std::complex<double>, P::n / 2> &slots)
{
    const double scale = std::ldexp(1.0, LogDelta);
    const auto &slot_to_eval = ckksSlotToEvalIndex<P>();
    for (std::uint32_t i = 0; i < P::n / 2; i++) {
        const std::uint32_t eval_index = slot_to_eval[i];
        evals[eval_index] = slots[i] * scale;
        evals[P::n - 1 - eval_index] = std::conj(evals[eval_index]);
    }
}

template <class P, std::uint32_t LogQ>
inline void encryptPolynomialAtLevel(TRLWE<P> &ct, const Polynomial<P> &poly,
                                     const Key<P> &key,
                                     CKKSNoise noise = {P::α, 0})
{
    static_assert(supported_torus_v<P>);

    for (std::uint32_t i = 0; i < P::n; i++)
        ct[P::k][i] =
            reduceToLevel<P, LogQ>(poly[i] + sampleNoise<P, LogQ>(noise));

    auto partkey = std::make_unique<Polynomial<P>>();
    auto mask_phase = std::make_unique<Polynomial<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        for (std::uint32_t i = 0; i < P::n; i++)
            ct[k][i] = uniformAtLevel<P, LogQ>();

        for (std::uint32_t i = 0; i < P::n; i++)
            (*partkey)[i] = key[k * P::n + i];

        if constexpr (is_multilimb_uint_v<typename P::T>)
            PolyMulTorusByDigit<P>(*mask_phase, ct[k], *partkey);
        else
            PolyMul<P>(*mask_phase, ct[k], *partkey);
        for (std::uint32_t i = 0; i < P::n; i++)
            ct[P::k][i] =
                reduceToLevel<P, LogQ>(ct[P::k][i] + (*mask_phase)[i]);
    }
}

template <class P, std::uint32_t LogQ>
inline void expandSeededKeySwitchRow(
    TRLWE<P> &ct, const CKKSSeededKeySwitchRow<P, LogQ> &row)
{
    static_assert(supported_torus_v<P>);
    for (int k = 0; k < static_cast<int>(P::k); k++)
        fillSeededUniformPolynomialAtLevel<P, LogQ>(
            ct[static_cast<std::size_t>(k)],
            row.mask_seeds[static_cast<std::size_t>(k)]);
    ct[P::k] = row.body;
}

template <class P, std::uint32_t LogQ>
inline void encryptPolynomialAtLevel(CKKSSeededKeySwitchRow<P, LogQ> &row,
                                     const Polynomial<P> &poly,
                                     const Key<P> &key,
                                     CKKSNoise noise = {P::α, 0})
{
    static_assert(supported_torus_v<P>);

    for (std::uint32_t i = 0; i < P::n; i++)
        row.body[i] =
            reduceToLevel<P, LogQ>(poly[i] + sampleNoise<P, LogQ>(noise));

    auto partkey = std::make_unique<Polynomial<P>>();
    auto mask = std::make_unique<Polynomial<P>>();
    auto mask_phase = std::make_unique<Polynomial<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        row.mask_seeds[static_cast<std::size_t>(k)] = randomSeed();
        fillSeededUniformPolynomialAtLevel<P, LogQ>(
            *mask, row.mask_seeds[static_cast<std::size_t>(k)]);

        for (std::uint32_t i = 0; i < P::n; i++)
            (*partkey)[i] = key[k * P::n + i];

        if constexpr (is_multilimb_uint_v<typename P::T>)
            PolyMulTorusByDigit<P>(*mask_phase, *mask, *partkey);
        else
            PolyMul<P>(*mask_phase, *mask, *partkey);
        for (std::uint32_t i = 0; i < P::n; i++)
            row.body[i] =
                reduceToLevel<P, LogQ>(row.body[i] + (*mask_phase)[i]);
    }
}

template <class P>
inline void polyMulTorusByBbarDigit(Polynomial<P> &res,
                                    const Polynomial<P> &torus,
                                    const Polynomial<P> &digit)
{
    if constexpr (is_multilimb_uint_v<typename P::T>) {
        PolyMulTorusByDigit<P>(res, torus, digit);
    }
    else {
        static_assert(std::is_same_v<typename P::T, __uint128_t>);
        static_assert(P::l̅ * P::B̅gbit ==
                          std::numeric_limits<typename P::T>::digits,
                      "CKKS Bbar digits must cover the torus");

        using T = typename P::T;
        constexpr int width = std::numeric_limits<T>::digits;
        constexpr T half = static_cast<T>(1) << (P::B̅gbit - 1);
        constexpr T mask = (static_cast<T>(1) << P::B̅gbit) - 1;
        constexpr T offset = [] {
            constexpr int local_width = std::numeric_limits<T>::digits;
            constexpr T local_half = static_cast<T>(1) << (P::B̅gbit - 1);
            T value = 0;
            for (int j = 0; j < static_cast<int>(P::l̅); j++)
                value += local_half << (local_width - (j + 1) * P::B̅gbit);
            return value;
        }();

        for (std::uint32_t n = 0; n < P::n; n++) res[n] = 0;

        auto torus_digit = std::make_unique<Polynomial<P>>();
        auto product = std::make_unique<Polynomial<P>>();
        auto digit_fft = std::make_unique<PolynomialInFD<P>>();
        auto torus_digit_fft = std::make_unique<PolynomialInFD<P>>();
        auto product_fft = std::make_unique<PolynomialInFD<P>>();
        TwistIFFT<P>(*digit_fft, digit);

        for (int j = 0; j < static_cast<int>(P::l̅); j++) {
            const int shift = width - (j + 1) * static_cast<int>(P::B̅gbit);
            for (std::uint32_t n = 0; n < P::n; n++)
                (*torus_digit)[n] =
                    (((torus[n] + offset) >> shift) & mask) - half;

            TwistIFFT<P>(*torus_digit_fft, *torus_digit);
            MulInFD<P::n>(*product_fft, *torus_digit_fft, *digit_fft);
            TwistFFT<P>(*product, *product_fft);

            for (std::uint32_t n = 0; n < P::n; n++)
                res[n] += (*product)[n] << shift;
        }
    }
}

template <class P, std::uint32_t FirstRow, std::size_t RowCount>
inline void baseBbarDecomposeRows(std::array<TRLWE<P>, RowCount> &result,
                                  const TRLWE<P> &input)
{
    static_assert(P::l̅ * P::B̅gbit ==
                      std::numeric_limits<typename P::T>::digits,
                  "CKKS Bbar digits must cover the torus");
    static_assert(RowCount > 0);
    static_assert(FirstRow + RowCount <= P::l̅);

    using T = typename P::T;
    constexpr int width = std::numeric_limits<T>::digits;
    constexpr T half = T{1} << (P::B̅gbit - 1);
    constexpr T mask = (T{1} << P::B̅gbit) - T{1};
    constexpr T offset = [] {
        constexpr int local_width = std::numeric_limits<T>::digits;
        constexpr T local_half = T{1} << (P::B̅gbit - 1);
        T value = 0;
        for (int j = 0; j < static_cast<int>(P::l̅); j++)
            value += local_half
                     << (local_width - (j + 1) * P::B̅gbit);
        return value;
    }();

    for (int c = 0; c <= static_cast<int>(P::k); c++) {
        for (std::uint32_t n = 0; n < P::n; n++) {
            const T a = input[c][n] + offset;
            for (std::size_t row = 0; row < RowCount; row++) {
                const std::uint32_t full_row =
                    FirstRow + static_cast<std::uint32_t>(row);
                const int shift =
                    width - (full_row + 1) * static_cast<int>(P::B̅gbit);
                result[row][c][n] = ((a >> shift) & mask) - half;
            }
        }
    }
}

template <class P, std::uint32_t FirstRow, std::size_t RowCount>
inline void baseBbarDecomposePolynomialRows(
    std::array<Polynomial<P>, RowCount> &result, const Polynomial<P> &input)
{
    static_assert(P::l̅ * P::B̅gbit ==
                      std::numeric_limits<typename P::T>::digits,
                  "CKKS Bbar digits must cover the torus");
    static_assert(RowCount > 0);
    static_assert(FirstRow + RowCount <= P::l̅);

    using T = typename P::T;
    constexpr int width = std::numeric_limits<T>::digits;
    constexpr T half = T{1} << (P::B̅gbit - 1);
    constexpr T mask = (T{1} << P::B̅gbit) - T{1};
    constexpr T offset = [] {
        constexpr int local_width = std::numeric_limits<T>::digits;
        constexpr T local_half = T{1} << (P::B̅gbit - 1);
        T value = 0;
        for (int j = 0; j < static_cast<int>(P::l̅); j++)
            value += local_half
                     << (local_width - (j + 1) * P::B̅gbit);
        return value;
    }();

    for (std::uint32_t n = 0; n < P::n; n++) {
        const T a = input[n] + offset;
        for (std::size_t row = 0; row < RowCount; row++) {
            const std::uint32_t full_row =
                FirstRow + static_cast<std::uint32_t>(row);
            const int shift =
                width - (full_row + 1) * static_cast<int>(P::B̅gbit);
            result[row][n] = ((a >> shift) & mask) - half;
        }
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline typename P::T ckksEncodeCoeff(double value)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta < LogQ);
    const long double scaled =
        std::ldexp(static_cast<long double>(value), LogDelta);
    return ckks_detail::signedToLevel<P, LogQ>(
        static_cast<__int128_t>(std::round(scaled)));
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline double ckksDecodeCoeff(typename P::T value)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta < LogQ);
    const long double decoded =
        std::ldexp(ckks_detail::levelToLongDouble<P, LogQ>(value),
                   -static_cast<int>(LogDelta));
    return static_cast<double>(decoded);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksEncodePolynomial(Polynomial<P> &poly,
                                 const std::array<double, P::n> &coeffs)
{
    ckks_detail::encodeRealPolynomial<P, LogQ, LogDelta>(poly, coeffs);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksDecodePolynomial(std::array<double, P::n> &coeffs,
                                 const Polynomial<P> &poly)
{
    ckks_detail::decodeRealPolynomial<P, LogQ, LogDelta>(coeffs, poly);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotEncode(
    Polynomial<P> &poly,
    const std::array<std::complex<double>, P::n / 2> &slots)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta < LogQ);
    std::array<std::complex<double>, P::n> evals{};
    ckks_detail::fillConjugateSymmetricEvaluations<LogDelta, P>(evals, slots);

    constexpr double pi = 3.141592653589793238462643383279502884;
    std::vector<std::complex<long double>> spectrum(P::n);
    for (std::uint32_t k = 0; k < P::n; k++)
        spectrum[k] = static_cast<std::complex<long double>>(evals[k]);
    ckks_detail::fftInPlace(spectrum, false);

    for (std::uint32_t j = 0; j < P::n; j++) {
        const long double angle =
            -pi * static_cast<long double>(j) / static_cast<long double>(P::n);
        const auto coeff =
            spectrum[j] *
            std::complex<long double>(std::cos(angle), std::sin(angle)) /
            static_cast<long double>(P::n);
        poly[j] = ckks_detail::signedToLevel<P, LogQ>(
            static_cast<__int128_t>(std::round(coeff.real())));
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotDecode(
    std::array<std::complex<double>, P::n / 2> &slots,
    const Polynomial<P> &poly)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta < LogQ);
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double inv_scale = std::ldexp(1.0, -static_cast<int>(LogDelta));
    const auto &slot_to_eval = ckks_detail::ckksSlotToEvalIndex<P>();
    std::vector<std::complex<long double>> values(P::n);
    for (std::uint32_t j = 0; j < P::n; j++) {
        const long double angle =
            pi * static_cast<long double>(j) / static_cast<long double>(P::n);
        const long double coeff =
            ckks_detail::levelToLongDouble<P, LogQ>(poly[j]);
        values[j] =
            coeff * std::complex<long double>(std::cos(angle), std::sin(angle));
    }
    ckks_detail::fftInPlace(values, true);

    for (std::uint32_t k = 0; k < P::n / 2; k++) {
        const std::uint32_t eval_index = slot_to_eval[k];
        const auto value = values[eval_index] * static_cast<long double>(P::n);
        slots[k] = static_cast<std::complex<double>>(value) * inv_scale;
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksEncrypt(CKKSCiphertext<P, LogQ, LogDelta> &out,
                        const std::array<double, P::n> &coeffs,
                        const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    Polynomial<P> poly;
    ckksEncodePolynomial<P, LogQ, LogDelta>(poly, coeffs);
    ckks_detail::encryptPolynomialAtLevel<P, LogQ>(out.ct, poly, key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksEncryptPolynomial(CKKSCiphertext<P, LogQ, LogDelta> &out,
                                  const Polynomial<P> &poly,
                                  const Key<P> &key,
                                  CKKSNoise noise = {P::α, 0})
{
    ckks_detail::encryptPolynomialAtLevel<P, LogQ>(out.ct, poly, key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksDecrypt(std::array<double, P::n> &out,
                        const CKKSCiphertext<P, LogQ, LogDelta> &ct,
                        const Key<P> &key)
{
    Polynomial<P> phase = trlwePhase<P>(ct.ct, key);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(phase);
    ckksDecodePolynomial<P, LogQ, LogDelta>(out, phase);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotEncrypt(
    CKKSCiphertext<P, LogQ, LogDelta> &out,
    const std::array<std::complex<double>, P::n / 2> &slots,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    Polynomial<P> poly;
    ckksSlotEncode<P, LogQ, LogDelta>(poly, slots);
    ckksEncryptPolynomial<P, LogQ, LogDelta>(out, poly, key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void ckksSlotDecrypt(
    std::array<std::complex<double>, P::n / 2> &slots,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct, const Key<P> &key)
{
    Polynomial<P> phase = trlwePhase<P>(ct.ct, key);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(phase);
    ckksSlotDecode<P, LogQ, LogDelta>(slots, phase);
}

template <class P>
using CKKSSlotVector = std::array<std::complex<double>, P::n / 2>;

template <class P>
inline void rotateCKKSSlotVector(CKKSSlotVector<P> &out,
                                 const CKKSSlotVector<P> &in, int steps)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    for (int i = 0; i < half; i++) out[i] = in[(i + steps) % half];
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
struct CKKSPlainMulTraits {
    static_assert(PlainLogDelta < LogQ);
    static constexpr std::uint32_t log_q = LogQ - PlainLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
using CKKSPlainMulResult =
    typename CKKSPlainMulTraits<P, LogQ, LogDelta,
                                PlainLogDelta>::Ciphertext;

template <class P, std::uint32_t InLogQ, std::uint32_t OutLogQ,
          std::uint32_t LogDelta>
inline void CKKSModRaise(CKKSCiphertext<P, OutLogQ, LogDelta> &res,
                         const CKKSCiphertext<P, InLogQ, LogDelta> &ct)
{
    static_assert(InLogQ < OutLogQ);
    static_assert(OutLogQ <= ckks_detail::torus_width_v<P>);
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            res.ct[c][i] = ckks_detail::reduceToLevel<P, OutLogQ>(
                ckks_detail::centeredLevelToTorus<P, InLogQ>(ct.ct[c][i]));
}

template <class P, std::uint32_t InLogQ, std::uint32_t OutLogQ,
          std::uint32_t LogDelta>
inline void CKKSModRaiseRandomized(CKKSCiphertext<P, OutLogQ, LogDelta> &res,
                                   const CKKSCiphertext<P, InLogQ, LogDelta> &ct)
{
    static_assert(InLogQ < OutLogQ);
    static_assert(OutLogQ <= ckks_detail::torus_width_v<P>);
    constexpr std::uint32_t high_log_q = OutLogQ - InLogQ;
    for (int c = 0; c <= static_cast<int>(P::k); c++) {
        for (std::uint32_t i = 0; i < P::n; i++) {
            const typename P::T low =
                ckks_detail::centeredLevelToTorus<P, InLogQ>(ct.ct[c][i]);
            const typename P::T high =
                ckks_detail::uniformAtLevel<P, high_log_q>() << InLogQ;
            res.ct[c][i] = ckks_detail::reduceToLevel<P, OutLogQ>(low + high);
        }
    }
}

template <class P, std::uint32_t InLogQ, std::uint32_t OutLogQ,
          std::uint32_t LogDelta, std::uint32_t MaskBound>
inline void CKKSModRaiseBoundedPhaseRandomized(
    CKKSCiphertext<P, OutLogQ, LogDelta> &res,
    const CKKSCiphertext<P, InLogQ, LogDelta> &ct)
{
    static_assert(InLogQ < OutLogQ);
    static_assert(OutLogQ <= ckks_detail::torus_width_v<P>);
    CKKSModRaise<P, InLogQ, OutLogQ, LogDelta>(res, ct);

    if constexpr (MaskBound > 0) {
        std::uniform_int_distribution<std::int64_t> dist(
            -static_cast<std::int64_t>(MaskBound),
            static_cast<std::int64_t>(MaskBound));
        for (std::uint32_t i = 0; i < P::n; i++) {
            const std::int64_t mask = dist(generator);
            const typename P::T high_mask =
                ckks_detail::reduceToLevel<P, OutLogQ>(
                    ckks_detail::signedToLevel<P, OutLogQ>(
                        static_cast<__int128_t>(mask))
                    << InLogQ);
            res.ct[P::k][i] = ckks_detail::reduceToLevel<P, OutLogQ>(
                res.ct[P::k][i] + high_mask);
        }
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t PlainLogDelta>
inline void CKKSPlainMulRescaleTRLWE(TRLWE<P> &res, const TRLWE<P> &ct,
                                     const Polynomial<P> &plain)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    auto product = std::make_unique<Polynomial<P>>();
    for (int c = 0; c <= static_cast<int>(P::k); c++) {
        ckks_detail::polyMulTorusBySmallSigned<P, LogQ>(*product, ct[c], plain);
        for (std::uint32_t i = 0; i < P::n; i++)
            res[c][i] =
                ckks_detail::roundedLevelRightShift<P, LogQ, PlainLogDelta>(
                    (*product)[i]);
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSPlainMulRescale(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSPlaintext<P, LogQ, PlainLogDelta> &plain)
{
    CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(res.ct, ct.ct,
                                                     plain.poly);
}

template <class P, std::uint32_t LogQ>
inline void CKKSAddTRLWEInPlace(TRLWE<P> &acc, const TRLWE<P> &term)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            acc[c][i] = ckks_detail::reduceToLevel<P, LogQ>(acc[c][i] +
                                                            term[c][i]);
}

template <class P, std::uint32_t LogQ>
inline void CKKSSubTRLWEInPlace(TRLWE<P> &acc, const TRLWE<P> &term)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            acc[c][i] = ckks_detail::reduceToLevel<P, LogQ>(acc[c][i] -
                                                            term[c][i]);
}

template <class P, std::uint32_t InLogQ, std::uint32_t OutLogQ,
          std::uint32_t LogDelta>
inline void CKKSLevelReduce(CKKSCiphertext<P, OutLogQ, LogDelta> &res,
                            const CKKSCiphertext<P, InLogQ, LogDelta> &ct)
{
    static_assert(OutLogQ <= InLogQ);
    static_assert(LogDelta < OutLogQ);
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            res.ct[c][i] = ckks_detail::reduceToLevel<P, OutLogQ>(ct.ct[c][i]);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t DropBits>
struct CKKSScaleDownTraits {
    static_assert(DropBits > 0);
    static_assert(DropBits < LogQ);
    static_assert(DropBits <= LogDelta);
    static constexpr std::uint32_t log_q = LogQ - DropBits;
    static constexpr std::uint32_t log_delta = LogDelta - DropBits;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t DropBits>
using CKKSScaleDownResult =
    typename CKKSScaleDownTraits<P, LogQ, LogDelta, DropBits>::Ciphertext;

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t DropBits>
inline void CKKSScaleDown(
    CKKSScaleDownResult<P, LogQ, LogDelta, DropBits> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct)
{
    static_assert(DropBits > 0);
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            res.ct[c][i] =
                ckks_detail::roundedLevelRightShift<P, LogQ, DropBits>(
                    ct.ct[c][i]);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSAddInPlace(CKKSCiphertext<P, LogQ, LogDelta> &acc,
                           const CKKSCiphertext<P, LogQ, LogDelta> &term)
{
    CKKSAddTRLWEInPlace<P, LogQ>(acc.ct, term.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSAdd(CKKSCiphertext<P, LogQ, LogDelta> &res,
                    const CKKSCiphertext<P, LogQ, LogDelta> &lhs,
                    const CKKSCiphertext<P, LogQ, LogDelta> &rhs)
{
    res = lhs;
    CKKSAddInPlace<P, LogQ, LogDelta>(res, rhs);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSubInPlace(CKKSCiphertext<P, LogQ, LogDelta> &acc,
                           const CKKSCiphertext<P, LogQ, LogDelta> &term)
{
    CKKSSubTRLWEInPlace<P, LogQ>(acc.ct, term.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSub(CKKSCiphertext<P, LogQ, LogDelta> &res,
                    const CKKSCiphertext<P, LogQ, LogDelta> &lhs,
                    const CKKSCiphertext<P, LogQ, LogDelta> &rhs)
{
    res = lhs;
    CKKSSubInPlace<P, LogQ, LogDelta>(res, rhs);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSMulIntegerInPlace(CKKSCiphertext<P, LogQ, LogDelta> &ct,
                                  std::int64_t scalar)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            ct.ct[c][i] =
                ckks_detail::reduceToLevel<P, LogQ>(ct.ct[c][i] * scalar);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSMulInteger(CKKSCiphertext<P, LogQ, LogDelta> &res,
                           const CKKSCiphertext<P, LogQ, LogDelta> &ct,
                           std::int64_t scalar)
{
    res = ct;
    CKKSMulIntegerInPlace<P, LogQ, LogDelta>(res, scalar);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSPlainMulByConstant(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    std::complex<double> scalar)
{
    CKKSSlotVector<P> slots{};
    slots.fill(scalar);
    CKKSPlaintext<P, LogQ, PlainLogDelta> plain;
    ckksSlotEncode<P, LogQ, PlainLogDelta>(plain.poly, slots);
    CKKSPlainMulRescale<P, LogQ, LogDelta, PlainLogDelta>(res, ct, plain);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSPlainMulByReal(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct, double scalar)
{
    CKKSPlainMulByConstant<P, LogQ, LogDelta, PlainLogDelta>(
        res, ct, {scalar, 0.0});
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSetTransparentReal(CKKSCiphertext<P, LogQ, LogDelta> &res,
                                   double value)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t i = 0; i < P::n; i++) res.ct[c][i] = 0;

    CKKSSlotVector<P> slots{};
    slots.fill({value, 0.0});
    ckksSlotEncode<P, LogQ, LogDelta>(res.ct[P::k], slots);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSAddPlainRealInPlace(CKKSCiphertext<P, LogQ, LogDelta> &acc,
                                    double value)
{
    CKKSSlotVector<P> slots{};
    slots.fill({value, 0.0});
    Polynomial<P> plain;
    ckksSlotEncode<P, LogQ, LogDelta>(plain, slots);
    for (std::uint32_t i = 0; i < P::n; i++)
        acc.ct[P::k][i] =
            ckks_detail::reduceToLevel<P, LogQ>(acc.ct[P::k][i] + plain[i]);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSubPlainRealInPlace(CKKSCiphertext<P, LogQ, LogDelta> &acc,
                                    double value)
{
    CKKSAddPlainRealInPlace<P, LogQ, LogDelta>(acc, -value);
}

template <class P, std::uint32_t LogQ, class Rows>
    requires(!ckks_detail::is_seeded_key_switch_row_array_v<Rows>)
inline void CKKSKeySwitchRows(TRLWE<P> &res, const Polynomial<P> &poly,
                              const Rows &rows)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(P::l̅ * P::B̅gbit ==
                      std::numeric_limits<typename P::T>::digits,
                  "CKKS key switching requires full Bbar decomposition");
    constexpr std::uint32_t row_count =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row =
        CKKSKeySwitchFirstRowForLevel<P, LogQ>();
    static_assert(std::tuple_size_v<std::remove_reference_t<Rows>> ==
                  row_count);

    auto centered = std::make_unique<Polynomial<P>>();
    ckks_detail::centeredPolynomialAtLevel<P, LogQ>(*centered, poly);

    auto decomposed = std::make_unique<std::array<Polynomial<P>, row_count>>();
    ckks_detail::baseBbarDecomposePolynomialRows<P, first_row>(*decomposed,
                                                               *centered);

    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++) res[c][n] = 0;

#pragma omp parallel
    {
        auto local = std::make_unique<TRLWE<P>>();
        auto product = std::make_unique<Polynomial<P>>();

        for (int c = 0; c <= static_cast<int>(P::k); c++)
            for (std::uint32_t n = 0; n < P::n; n++) (*local)[c][n] = 0;

#pragma omp for schedule(dynamic)
        for (int j = 0; j < static_cast<int>(row_count); j++) {
            const Polynomial<P> &digit = (*decomposed)[j];
            for (int c = 0; c <= static_cast<int>(P::k); c++) {
                ckks_detail::polyMulTorusByBbarDigit<P>(*product, rows[j][c],
                                                        digit);
                for (std::uint32_t n = 0; n < P::n; n++)
                    (*local)[c][n] += (*product)[n];
            }
        }

#pragma omp critical
        {
            for (int c = 0; c <= static_cast<int>(P::k); c++)
                for (std::uint32_t n = 0; n < P::n; n++)
                    res[c][n] =
                        ckks_detail::reduceToLevel<P, LogQ>(res[c][n] +
                                                            (*local)[c][n]);
        }
    }
}

template <class P, std::uint32_t LogQ, std::size_t RowCount>
inline void CKKSKeySwitchRows(
    TRLWE<P> &res, const Polynomial<P> &poly,
    const std::array<CKKSSeededKeySwitchRow<P, LogQ>, RowCount> &rows)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(P::l̅ * P::B̅gbit ==
                      std::numeric_limits<typename P::T>::digits,
                  "CKKS key switching requires full Bbar decomposition");
    constexpr std::uint32_t row_count =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row =
        CKKSKeySwitchFirstRowForLevel<P, LogQ>();
    static_assert(RowCount == row_count);

    auto centered = std::make_unique<Polynomial<P>>();
    ckks_detail::centeredPolynomialAtLevel<P, LogQ>(*centered, poly);

    auto decomposed = std::make_unique<std::array<Polynomial<P>, row_count>>();
    ckks_detail::baseBbarDecomposePolynomialRows<P, first_row>(*decomposed,
                                                               *centered);

    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++) res[c][n] = 0;

#pragma omp parallel
    {
        auto local = std::make_unique<TRLWE<P>>();
        auto mask = std::make_unique<Polynomial<P>>();
        auto product = std::make_unique<Polynomial<P>>();

        for (int c = 0; c <= static_cast<int>(P::k); c++)
            for (std::uint32_t n = 0; n < P::n; n++) (*local)[c][n] = 0;

#pragma omp for schedule(dynamic)
        for (int j = 0; j < static_cast<int>(row_count); j++) {
            const Polynomial<P> &digit = (*decomposed)[j];
            const auto &row = rows[static_cast<std::size_t>(j)];
            for (int c = 0; c < static_cast<int>(P::k); c++) {
                ckks_detail::fillSeededUniformPolynomialAtLevel<P, LogQ>(
                    *mask, row.mask_seeds[static_cast<std::size_t>(c)]);
                ckks_detail::polyMulTorusByBbarDigit<P>(*product, *mask,
                                                        digit);
                for (std::uint32_t n = 0; n < P::n; n++)
                    (*local)[c][n] += (*product)[n];
            }
            ckks_detail::polyMulTorusByBbarDigit<P>(*product, row.body, digit);
            for (std::uint32_t n = 0; n < P::n; n++)
                (*local)[P::k][n] += (*product)[n];
        }

#pragma omp critical
        {
            for (int c = 0; c <= static_cast<int>(P::k); c++)
                for (std::uint32_t n = 0; n < P::n; n++)
                    res[c][n] =
                        ckks_detail::reduceToLevel<P, LogQ>(res[c][n] +
                                                            (*local)[c][n]);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSSecretKeySwitchKeyGen(CKKSSecretKeySwitchKey<P, LogQ> &switch_key,
                                      const Key<P> &source_key,
                                      const Key<P> &target_key,
                                      CKKSNoise noise = {P::α, 0})
{
    static_assert(ckks_detail::supported_torus_v<P>);
    constexpr int width = std::numeric_limits<typename P::T>::digits;
    constexpr std::uint32_t row_count =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row =
        CKKSKeySwitchFirstRowForLevel<P, LogQ>();

    auto partkey = std::make_unique<Polynomial<P>>();
    auto gadget = std::make_unique<Polynomial<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        for (std::uint32_t i = 0; i < P::n; i++)
            (*partkey)[i] = source_key[k * P::n + i];

        for (int j = 0; j < static_cast<int>(row_count); j++) {
            const std::uint32_t full_row = first_row + j;
            const int shift =
                width - (full_row + 1) * static_cast<int>(P::B̅gbit);
            for (std::uint32_t n = 0; n < P::n; n++)
                (*gadget)[n] =
                    ckks_detail::reduceToLevel<P, LogQ>((*partkey)[n] << shift);
            ckks_detail::encryptPolynomialAtLevel<P, LogQ>(
                switch_key[static_cast<std::size_t>(k)]
                          [static_cast<std::size_t>(j)],
                *gadget, target_key, noise);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSSecretKeySwitchKeyGen(
    CKKSSeededSecretKeySwitchKey<P, LogQ> &switch_key,
    const Key<P> &source_key, const Key<P> &target_key,
    CKKSNoise noise = {P::α, 0})
{
    static_assert(ckks_detail::supported_torus_v<P>);
    constexpr int width = std::numeric_limits<typename P::T>::digits;
    constexpr std::uint32_t row_count =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row =
        CKKSKeySwitchFirstRowForLevel<P, LogQ>();

    auto partkey = std::make_unique<Polynomial<P>>();
    auto gadget = std::make_unique<Polynomial<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        for (std::uint32_t i = 0; i < P::n; i++)
            (*partkey)[i] = source_key[k * P::n + i];

        for (int j = 0; j < static_cast<int>(row_count); j++) {
            const std::uint32_t full_row = first_row + j;
            const int shift =
                width - (full_row + 1) * static_cast<int>(P::B̅gbit);
            for (std::uint32_t n = 0; n < P::n; n++)
                (*gadget)[n] =
                    ckks_detail::reduceToLevel<P, LogQ>((*partkey)[n] << shift);
            ckks_detail::encryptPolynomialAtLevel<P, LogQ>(
                switch_key[static_cast<std::size_t>(k)]
                          [static_cast<std::size_t>(j)],
                *gadget, target_key, noise);
        }
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSecretKeySwitch(
    CKKSCiphertext<P, LogQ, LogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSSecretKeySwitchKey<P, LogQ> &switch_key)
{
    auto switched = std::make_unique<TRLWE<P>>();
    auto out = std::make_unique<TRLWE<P>>();
    for (std::uint32_t n = 0; n < P::n; n++)
        (*out)[P::k][n] = ckks_detail::reduceToLevel<P, LogQ>(ct.ct[P::k][n]);

    for (int k = 0; k < static_cast<int>(P::k); k++) {
        CKKSKeySwitchRows<P, LogQ>(
            *switched, ct.ct[static_cast<std::size_t>(k)],
            switch_key[static_cast<std::size_t>(k)]);
        CKKSSubTRLWEInPlace<P, LogQ>(*out, *switched);
    }
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(*out);
    res.ct = *out;
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void CKKSSecretKeySwitch(
    CKKSCiphertext<P, LogQ, LogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSSeededSecretKeySwitchKey<P, LogQ> &switch_key)
{
    auto switched = std::make_unique<TRLWE<P>>();
    auto out = std::make_unique<TRLWE<P>>();
    for (std::uint32_t n = 0; n < P::n; n++)
        (*out)[P::k][n] = ckks_detail::reduceToLevel<P, LogQ>(ct.ct[P::k][n]);

    for (int k = 0; k < static_cast<int>(P::k); k++) {
        CKKSKeySwitchRows<P, LogQ>(
            *switched, ct.ct[static_cast<std::size_t>(k)],
            switch_key[static_cast<std::size_t>(k)]);
        CKKSSubTRLWEInPlace<P, LogQ>(*out, *switched);
    }
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(*out);
    res.ct = *out;
}

template <class P, std::uint32_t LogQ>
inline void CKKSAutoKeyGen(CKKSAutoKey<P, LogQ> &autokey, const uint d,
                           const Key<P> &key,
                           CKKSNoise noise = {P::α, 0})
{
    static_assert(ckks_detail::supported_torus_v<P>);
    constexpr int width = std::numeric_limits<typename P::T>::digits;
    constexpr std::uint32_t row_count =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row =
        CKKSKeySwitchFirstRowForLevel<P, LogQ>();

    auto partkey = std::make_unique<Polynomial<P>>();
    auto autokey_poly = std::make_unique<Polynomial<P>>();
    auto gadget = std::make_unique<Polynomial<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        for (std::uint32_t i = 0; i < P::n; i++)
            (*partkey)[i] = key[k * P::n + i];

        Automorphism<P>(*autokey_poly, *partkey, d);

        for (int j = 0; j < static_cast<int>(row_count); j++) {
            const std::uint32_t full_row = first_row + j;
            const int shift =
                width - (full_row + 1) * static_cast<int>(P::B̅gbit);
            for (std::uint32_t n = 0; n < P::n; n++)
                (*gadget)[n] = ckks_detail::reduceToLevel<P, LogQ>(
                    (*autokey_poly)[n] << shift);
            ckks_detail::encryptPolynomialAtLevel<P, LogQ>(
                autokey[k][j], *gadget, key, noise);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSAutoKeyGen(CKKSSeededAutoKey<P, LogQ> &autokey, const uint d,
                           const Key<P> &key,
                           CKKSNoise noise = {P::α, 0})
{
    static_assert(ckks_detail::supported_torus_v<P>);
    constexpr int width = std::numeric_limits<typename P::T>::digits;
    constexpr std::uint32_t row_count =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row =
        CKKSKeySwitchFirstRowForLevel<P, LogQ>();

    auto partkey = std::make_unique<Polynomial<P>>();
    auto autokey_poly = std::make_unique<Polynomial<P>>();
    auto gadget = std::make_unique<Polynomial<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        for (std::uint32_t i = 0; i < P::n; i++)
            (*partkey)[i] = key[k * P::n + i];

        Automorphism<P>(*autokey_poly, *partkey, d);

        for (int j = 0; j < static_cast<int>(row_count); j++) {
            const std::uint32_t full_row = first_row + j;
            const int shift =
                width - (full_row + 1) * static_cast<int>(P::B̅gbit);
            for (std::uint32_t n = 0; n < P::n; n++)
                (*gadget)[n] = ckks_detail::reduceToLevel<P, LogQ>(
                    (*autokey_poly)[n] << shift);
            ckks_detail::encryptPolynomialAtLevel<P, LogQ>(
                autokey[k][j], *gadget, key, noise);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSGaloisKeyGen(CKKSGaloisKey<P, LogQ> &gk, const Key<P> &key,
                             CKKSNoise noise = {P::α, 0})
{
    std::uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit); i++) {
        CKKSAutoKeyGen<P, LogQ>(gk[i], static_cast<uint>(d), key, noise);
        d = d * d % (2 * P::n);
    }
    CKKSAutoKeyGen<P, LogQ>(gk[P::nbit], 2 * P::n - 1, key, noise);
}

template <class P, std::uint32_t LogQ>
inline void CKKSSeededGaloisKeyGen(CKKSSeededGaloisKey<P, LogQ> &gk,
                                   const Key<P> &key,
                                   CKKSNoise noise = {P::α, 0})
{
    std::uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit); i++) {
        CKKSAutoKeyGen<P, LogQ>(gk[i], static_cast<uint>(d), key, noise);
        d = d * d % (2 * P::n);
    }
    CKKSAutoKeyGen<P, LogQ>(gk[P::nbit], 2 * P::n - 1, key, noise);
}

namespace ckks_detail {

template <class P>
inline int rotationAutomorphismForSteps(int steps)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    constexpr std::uint64_t mod = 2 * static_cast<std::uint64_t>(P::n);
    std::uint64_t d = 1;
    std::uint64_t base = 5;
    std::uint64_t exp = static_cast<std::uint64_t>(steps);
    while (exp != 0) {
        if ((exp & 1U) != 0) d = (d * base) % mod;
        base = (base * base) % mod;
        exp >>= 1U;
    }
    return static_cast<int>(d);
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ>
inline void CKKSSparseGaloisKeyGen(CKKSSparseGaloisKey<P, LogQ> &gk,
                                   const Key<P> &key,
                                   const CKKSRotationKeyIndexSet<P> &indices,
                                   CKKSNoise noise = {P::α, 0})
{
    gk.available = indices;
    std::uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit); i++) {
        if (indices[static_cast<std::size_t>(i)]) {
            gk.keys[static_cast<std::size_t>(i)] =
                std::make_unique<CKKSAutoKey<P, LogQ>>();
            CKKSAutoKeyGen<P, LogQ>(
                *gk.keys[static_cast<std::size_t>(i)], static_cast<uint>(d),
                key, noise);
        }
        else {
            gk.keys[static_cast<std::size_t>(i)].reset();
        }
        d = d * d % (2 * P::n);
    }

    if (indices[P::nbit]) {
        gk.keys[P::nbit] = std::make_unique<CKKSAutoKey<P, LogQ>>();
        CKKSAutoKeyGen<P, LogQ>(*gk.keys[P::nbit], 2 * P::n - 1, key,
                                noise);
    }
    else {
        gk.keys[P::nbit].reset();
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSDirectSparseGaloisKeyGen(
    CKKSDirectSparseGaloisKey<P, LogQ> &gk, const Key<P> &key,
    const CKKSDirectRotationKeyIndexSet<P> &indices,
    CKKSNoise noise = {P::α, 0})
{
    gk.keys.clear();
    gk.keys.reserve(CKKSDirectRotationKeyIndexSetCount<P>(indices));
    for (int steps = 1; steps < static_cast<int>(P::n / 2); steps++) {
        if (!indices[static_cast<std::size_t>(steps)]) continue;
        auto &entry = gk.keys.emplace_back();
        entry.steps = steps;
        CKKSAutoKeyGen<P, LogQ>(
            entry.key, ckks_detail::rotationAutomorphismForSteps<P>(steps),
            key, noise);
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSHybridSparseGaloisKeyGen(
    CKKSHybridSparseGaloisKey<P, LogQ> &gk, const Key<P> &key,
    const CKKSRotationKeyIndexSet<P> &binary_indices,
    const CKKSDirectRotationKeyIndexSet<P> &direct_indices,
    CKKSNoise noise = {P::α, 0})
{
    CKKSSparseGaloisKeyGen<P, LogQ>(gk.binary, key, binary_indices, noise);
    CKKSDirectSparseGaloisKeyGen<P, LogQ>(gk.direct, key, direct_indices,
                                          noise);
}

template <class P, std::uint32_t LogQ>
inline void CKKSSparseGaloisKeyGen(
    CKKSSeededSparseGaloisKey<P, LogQ> &gk, const Key<P> &key,
    const CKKSRotationKeyIndexSet<P> &indices, CKKSNoise noise = {P::α, 0})
{
    gk.available = indices;
    std::uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit); i++) {
        if (indices[static_cast<std::size_t>(i)]) {
            gk.keys[static_cast<std::size_t>(i)] =
                std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
            CKKSAutoKeyGen<P, LogQ>(
                *gk.keys[static_cast<std::size_t>(i)], static_cast<uint>(d),
                key, noise);
        }
        else {
            gk.keys[static_cast<std::size_t>(i)].reset();
        }
        d = d * d % (2 * P::n);
    }

    if (indices[P::nbit]) {
        gk.keys[P::nbit] = std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
        CKKSAutoKeyGen<P, LogQ>(*gk.keys[P::nbit], 2 * P::n - 1, key,
                                noise);
    }
    else {
        gk.keys[P::nbit].reset();
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSDirectSparseGaloisKeyGen(
    CKKSSeededDirectSparseGaloisKey<P, LogQ> &gk, const Key<P> &key,
    const CKKSDirectRotationKeyIndexSet<P> &indices,
    CKKSNoise noise = {P::α, 0})
{
    gk.keys.clear();
    gk.keys.reserve(CKKSDirectRotationKeyIndexSetCount<P>(indices));
    for (int steps = 1; steps < static_cast<int>(P::n / 2); steps++) {
        if (!indices[static_cast<std::size_t>(steps)]) continue;
        auto &entry = gk.keys.emplace_back();
        entry.steps = steps;
        CKKSAutoKeyGen<P, LogQ>(
            entry.key, ckks_detail::rotationAutomorphismForSteps<P>(steps),
            key, noise);
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSHybridSparseGaloisKeyGen(
    CKKSSeededHybridSparseGaloisKey<P, LogQ> &gk, const Key<P> &key,
    const CKKSRotationKeyIndexSet<P> &binary_indices,
    const CKKSDirectRotationKeyIndexSet<P> &direct_indices,
    CKKSNoise noise = {P::α, 0})
{
    CKKSSparseGaloisKeyGen<P, LogQ>(gk.binary, key, binary_indices, noise);
    CKKSDirectSparseGaloisKeyGen<P, LogQ>(gk.direct, key, direct_indices,
                                          noise);
}

template <class P, std::uint32_t LogQ>
inline void CKKSEvalAuto(TRLWE<P> &res, const TRLWE<P> &trlwe, const int d,
                         const CKKSAutoKey<P, LogQ> &autokey)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++) res[c][n] = 0;
    Automorphism<P>(res[P::k], trlwe[P::k], d);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(res[P::k]);

    auto temppoly = std::make_unique<Polynomial<P>>();
    auto switched = std::make_unique<TRLWE<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        Automorphism<P>(*temppoly, trlwe[k], d);
        CKKSKeySwitchRows<P, LogQ>(*switched, *temppoly, autokey[k]);
        CKKSSubTRLWEInPlace<P, LogQ>(res, *switched);
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSEvalAuto(TRLWE<P> &res, const TRLWE<P> &trlwe, const int d,
                         const CKKSSeededAutoKey<P, LogQ> &autokey)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (std::uint32_t n = 0; n < P::n; n++) res[c][n] = 0;
    Automorphism<P>(res[P::k], trlwe[P::k], d);
    ckks_detail::reducePolynomialToLevel<P, LogQ>(res[P::k]);

    auto temppoly = std::make_unique<Polynomial<P>>();
    auto switched = std::make_unique<TRLWE<P>>();
    for (int k = 0; k < static_cast<int>(P::k); k++) {
        Automorphism<P>(*temppoly, trlwe[k], d);
        CKKSKeySwitchRows<P, LogQ>(*switched, *temppoly, autokey[k]);
        CKKSSubTRLWEInPlace<P, LogQ>(res, *switched);
    }
}

namespace ckks_detail {

template <class P, std::uint32_t LogQ, class KeyForBit>
inline void CKKSRotateSlotsBinaryChain(TRLWE<P> &res, const TRLWE<P> &ct,
                                       int steps, KeyForBit &&key_for_bit)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    if (steps == 0) {
        res = ct;
        reduceTRLWEToLevel<P, LogQ>(res);
        return;
    }

    auto first = std::make_unique<TRLWE<P>>();
    auto second = std::make_unique<TRLWE<P>>();
    const TRLWE<P> *cur = &ct;
    TRLWE<P> *next = first.get();

    std::uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit) - 1; i++) {
        if ((steps >> i) & 1) {
            CKKSEvalAuto<P, LogQ>(
                *next, *cur, static_cast<int>(d),
                key_for_bit(static_cast<std::size_t>(i)));
            cur = next;
            next = next == first.get() ? second.get() : first.get();
        }
        d = d * d % (2 * P::n);
    }
    res = *cur;
    reduceTRLWEToLevel<P, LogQ>(res);
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSGaloisKey<P, LogQ> &gk)
{
    ckks_detail::CKKSRotateSlotsBinaryChain<P, LogQ>(
        res, ct, steps,
        [&](std::size_t i) -> const CKKSAutoKey<P, LogQ> & { return gk[i]; });
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSSparseGaloisKey<P, LogQ> &gk)
{
    ckks_detail::CKKSRotateSlotsBinaryChain<P, LogQ>(
        res, ct, steps, [&](std::size_t i) -> const CKKSAutoKey<P, LogQ> & {
            return gk.get(i);
        });
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSSeededSparseGaloisKey<P, LogQ> &gk)
{
    ckks_detail::CKKSRotateSlotsBinaryChain<P, LogQ>(
        res, ct, steps,
        [&](std::size_t i) -> const CKKSSeededAutoKey<P, LogQ> & {
            return gk.get(i);
        });
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSDirectSparseGaloisKey<P, LogQ> &gk)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    if (steps == 0) {
        res = ct;
        ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
        return;
    }

    CKKSEvalAuto<P, LogQ>(
        res, ct, ckks_detail::rotationAutomorphismForSteps<P>(steps),
        gk.get(steps));
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSSeededDirectSparseGaloisKey<P, LogQ> &gk)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    if (steps == 0) {
        res = ct;
        ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
        return;
    }

    CKKSEvalAuto<P, LogQ>(
        res, ct, ckks_detail::rotationAutomorphismForSteps<P>(steps),
        gk.get(steps));
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSHybridSparseGaloisKey<P, LogQ> &gk)
{
    if (gk.direct.has(steps) &&
        CKKSDirectSparseGaloisKey<P, LogQ>::normalize_steps(steps) != 0) {
        CKKSRotateSlots<P, LogQ>(res, ct, steps, gk.direct);
    }
    else {
        CKKSRotateSlots<P, LogQ>(res, ct, steps, gk.binary);
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSSeededHybridSparseGaloisKey<P, LogQ> &gk)
{
    if (gk.direct.has(steps) &&
        CKKSSeededDirectSparseGaloisKey<P, LogQ>::normalize_steps(steps) != 0) {
        CKKSRotateSlots<P, LogQ>(res, ct, steps, gk.direct);
    }
    else {
        CKKSRotateSlots<P, LogQ>(res, ct, steps, gk.binary);
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSConjugateSlots(TRLWE<P> &res, const TRLWE<P> &ct,
                               const CKKSGaloisKey<P, LogQ> &gk)
{
    CKKSEvalAuto<P, LogQ>(res, ct, 2 * P::n - 1, gk[P::nbit]);
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ>
inline void CKKSConjugateSlots(TRLWE<P> &res, const TRLWE<P> &ct,
                               const CKKSSeededGaloisKey<P, LogQ> &gk)
{
    CKKSEvalAuto<P, LogQ>(res, ct, 2 * P::n - 1, gk[P::nbit]);
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ>
inline void CKKSConjugateSlots(TRLWE<P> &res, const TRLWE<P> &ct,
                               const CKKSSparseGaloisKey<P, LogQ> &gk)
{
    CKKSEvalAuto<P, LogQ>(res, ct, 2 * P::n - 1, gk.get(P::nbit));
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ>
inline void CKKSConjugateSlots(TRLWE<P> &res, const TRLWE<P> &ct,
                               const CKKSSeededSparseGaloisKey<P, LogQ> &gk)
{
    CKKSEvalAuto<P, LogQ>(res, ct, 2 * P::n - 1, gk.get(P::nbit));
    ckks_detail::reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ>
inline void CKKSConjugateSlots(TRLWE<P> &res, const TRLWE<P> &ct,
                               const CKKSHybridSparseGaloisKey<P, LogQ> &gk)
{
    CKKSConjugateSlots<P, LogQ>(res, ct, gk.binary);
}

template <class P, std::uint32_t LogQ>
inline void CKKSConjugateSlots(TRLWE<P> &res, const TRLWE<P> &ct,
                               const CKKSSeededHybridSparseGaloisKey<P, LogQ> &gk)
{
    CKKSConjugateSlots<P, LogQ>(res, ct, gk.binary);
}

namespace ckks_detail {

inline int highestPowerOfTwoLE(int value)
{
    assert(value > 0);
    int bit = 1;
    while ((bit << 1) <= value) bit <<= 1;
    return bit;
}

inline std::size_t popcountNonnegativeInt(int value)
{
    assert(value >= 0);
    std::size_t count = 0;
    auto bits = static_cast<unsigned int>(value);
    while (bits != 0) {
        count += bits & 1U;
        bits >>= 1;
    }
    return count;
}

template <class P>
using CKKSSparseBabyStepTable = std::vector<std::unique_ptr<TRLWE<P>>>;

template <class P>
inline TRLWE<P> &CKKSEnsureSparseBabyStep(
    CKKSSparseBabyStepTable<P> &baby, std::size_t index)
{
    assert(index < baby.size());
    if (!baby[index]) baby[index] = std::make_unique<TRLWE<P>>();
    return *baby[index];
}

template <class P, std::uint32_t LogQ, class GaloisKey>
inline void CKKSBuildSparseBabyStepRotationTable(
    CKKSSparseBabyStepTable<P> &baby, const TRLWE<P> &ct,
    const std::vector<bool> &baby_used, const GaloisKey &gk)
{
    assert(!baby.empty());
    assert(baby.size() == baby_used.size());
    CKKSEnsureSparseBabyStep<P>(baby, 0) = ct;

    std::vector<bool> needed(baby.size(), false);
    needed[0] = true;
    for (int target = 1; target < static_cast<int>(baby_used.size());
         target++) {
        if (!baby_used[static_cast<std::size_t>(target)]) continue;

        int step = target;
        while (step != 0) {
            needed[static_cast<std::size_t>(step)] = true;
            step -= highestPowerOfTwoLE(step);
        }
    }

    for (std::size_t i = 1; i < needed.size(); i++)
        if (needed[i]) CKKSEnsureSparseBabyStep<P>(baby, i);

    for (int depth = 1; depth <= static_cast<int>(P::nbit); depth++) {
        for (int current = 1; current < static_cast<int>(needed.size());
             current++) {
            if (!needed[static_cast<std::size_t>(current)] ||
                static_cast<int>(popcountNonnegativeInt(current)) != depth)
                continue;
            const int bit = highestPowerOfTwoLE(current);
            const int parent = current - bit;
            assert(baby[static_cast<std::size_t>(parent)]);
            assert(baby[static_cast<std::size_t>(current)]);
            CKKSRotateSlots<P, LogQ>(
                *baby[static_cast<std::size_t>(current)],
                *baby[static_cast<std::size_t>(parent)], bit, gk);
        }
    }
}

template <class P, std::uint32_t LogQ, class GaloisKey>
inline void CKKSBuildSparseDirectBabyStepRotationTable(
    CKKSSparseBabyStepTable<P> &baby, const TRLWE<P> &ct,
    const std::vector<bool> &baby_used, const GaloisKey &gk)
{
    assert(!baby.empty());
    assert(baby.size() == baby_used.size());
    CKKSEnsureSparseBabyStep<P>(baby, 0) = ct;
    for (int j1 = 1; j1 < static_cast<int>(baby_used.size()); j1++) {
        if (baby_used[static_cast<std::size_t>(j1)])
            CKKSRotateSlots<P, LogQ>(
                CKKSEnsureSparseBabyStep<P>(
                    baby, static_cast<std::size_t>(j1)),
                ct, j1, gk);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class GaloisKey>
inline void CKKSExtractRealSlots(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct, const GaloisKey &gk)
{
    auto conjugated = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
    CKKSConjugateSlots<P, LogQ>(conjugated->ct, ct.ct, gk);

    auto sum = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
    CKKSAdd<P, LogQ, LogDelta>(*sum, ct, *conjugated);
    conjugated.reset();
    CKKSPlainMulByReal<P, LogQ, LogDelta, PlainLogDelta>(res, *sum, 0.5);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class GaloisKey>
inline void CKKSExtractImagSlots(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct, const GaloisKey &gk)
{
    auto conjugated = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
    CKKSConjugateSlots<P, LogQ>(conjugated->ct, ct.ct, gk);

    auto diff = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
    CKKSSub<P, LogQ, LogDelta>(*diff, ct, *conjugated);
    conjugated.reset();
    CKKSPlainMulByConstant<P, LogQ, LogDelta, PlainLogDelta>(
        res, *diff, {0.0, -0.5});
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class GaloisKey>
inline void CKKSExtractRealImagSlots(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &real,
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &imag,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct, const GaloisKey &gk)
{
    auto conjugated = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
    CKKSConjugateSlots<P, LogQ>(conjugated->ct, ct.ct, gk);

    auto sum = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
    CKKSAdd<P, LogQ, LogDelta>(*sum, ct, *conjugated);
    CKKSPlainMulByReal<P, LogQ, LogDelta, PlainLogDelta>(real, *sum, 0.5);
    sum.reset();

    auto diff = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
    CKKSSub<P, LogQ, LogDelta>(*diff, ct, *conjugated);
    conjugated.reset();
    CKKSPlainMulByConstant<P, LogQ, LogDelta, PlainLogDelta>(
        imag, *diff, {0.0, -0.5});
}

template <class P, std::uint32_t LogQ, std::uint32_t PlainLogDelta>
struct CKKSLinearTransformTerm {
    int baby_step = 0;
    CKKSPlaintext<P, LogQ, PlainLogDelta> plain{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(baby_step, plain);
    }
};

template <class P, std::uint32_t LogQ, std::uint32_t PlainLogDelta>
struct CKKSLinearTransformGiantStep {
    int giant_step = 0;
    std::vector<CKKSLinearTransformTerm<P, LogQ, PlainLogDelta>> terms{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(giant_step, terms);
    }
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
struct CKKSLinearTransformPlan {
    static_assert(PlainLogDelta < LogQ);
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t plain_log_delta = PlainLogDelta;
    static constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;

    int k_step = 0;
    std::vector<CKKSLinearTransformGiantStep<P, LogQ, PlainLogDelta>> groups{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(k_step, groups);
    }
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSLinearTransformPlanUsedBabySteps(
    std::vector<bool> &used,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan)
{
    assert(plan.k_step > 0);
    used.assign(static_cast<std::size_t>(plan.k_step), false);
    for (const auto &group : plan.groups) {
        for (const auto &entry : group.terms) {
            assert(entry.baby_step >= 0 && entry.baby_step < plan.k_step);
            used[static_cast<std::size_t>(entry.baby_step)] = true;
        }
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline std::size_t CKKSLinearTransformPlanBabyRotationCount(
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan)
{
    std::vector<bool> used;
    CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta, PlainLogDelta>(
        used, plan);
    std::size_t count = 0;
    for (std::size_t i = 1; i < used.size(); i++)
        if (used[i]) count++;
    return count;
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline std::size_t CKKSLinearTransformPlanGiantRotationCount(
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan)
{
    std::size_t count = 0;
    for (const auto &group : plan.groups)
        if (group.giant_step != 0) count++;
    return count;
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSCollectLinearTransformPlanRotationKeyIndices(
    CKKSRotationKeyIndexSet<P> &input_keys,
    CKKSRotationKeyIndexSet<P> &output_keys,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan)
{
    std::vector<bool> baby_used;
    CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta, PlainLogDelta>(
        baby_used, plan);
    for (int j1 = 1; j1 < plan.k_step; j1++) {
        if (baby_used[static_cast<std::size_t>(j1)])
            CKKSMarkRotationPowerKeyIndices<P>(input_keys, j1);
    }
    for (const auto &group : plan.groups) {
        if (group.giant_step != 0)
            CKKSMarkRotationPowerKeyIndices<P>(
                output_keys, plan.k_step * group.giant_step);
    }
}

template <class P>
struct CKKSLinearTransformStage {
    std::vector<CKKSSlotVector<P>> diagonals{};
    std::vector<int> rotation_offsets{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(diagonals, rotation_offsets);
    }
};

template <class P>
using CKKSLinearTransformStages = std::vector<CKKSLinearTransformStage<P>>;

template <class P>
inline void CKKSCollectLinearTransformStageRotationKeyIndices(
    CKKSRotationKeyIndexSet<P> &input_keys,
    CKKSRotationKeyIndexSet<P> &output_keys,
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    int max_j2 = 0;
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        max_j2 = std::max(max_j2, r / k_step);
    }

    std::vector<bool> baby_used(static_cast<std::size_t>(k_step), false);
    std::vector<bool> giant_used(static_cast<std::size_t>(max_j2 + 1), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        baby_used[static_cast<std::size_t>(r % k_step)] = true;
        giant_used[static_cast<std::size_t>(r / k_step)] = true;
    }

    for (int j1 = 1; j1 < k_step; j1++) {
        if (baby_used[static_cast<std::size_t>(j1)])
            CKKSMarkRotationPowerKeyIndices<P>(input_keys, j1);
    }
    for (int j2 = 1; j2 <= max_j2; j2++) {
        if (giant_used[static_cast<std::size_t>(j2)])
            CKKSMarkRotationPowerKeyIndices<P>(output_keys, k_step * j2);
    }
}

template <class P>
inline void CKKSCollectLinearTransformStageDirectRotationKeyIndices(
    CKKSDirectRotationKeyIndexSet<P> &input_keys,
    CKKSDirectRotationKeyIndexSet<P> &output_keys,
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    int max_j2 = 0;
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        max_j2 = std::max(max_j2, r / k_step);
    }

    std::vector<bool> baby_used(static_cast<std::size_t>(k_step), false);
    std::vector<bool> giant_used(static_cast<std::size_t>(max_j2 + 1), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        baby_used[static_cast<std::size_t>(r % k_step)] = true;
        giant_used[static_cast<std::size_t>(r / k_step)] = true;
    }

    for (int j1 = 1; j1 < k_step; j1++) {
        if (baby_used[static_cast<std::size_t>(j1)])
            CKKSMarkDirectRotationKeyIndex<P>(input_keys, j1);
    }
    for (int j2 = 1; j2 <= max_j2; j2++) {
        if (giant_used[static_cast<std::size_t>(j2)])
            CKKSMarkDirectRotationKeyIndex<P>(output_keys, k_step * j2);
    }
}

template <class P>
inline void CKKSCollectLinearTransformStageHybridGiantRotationKeyIndices(
    CKKSRotationKeyIndexSet<P> &input_binary_keys,
    CKKSRotationKeyIndexSet<P> &output_binary_keys,
    CKKSDirectRotationKeyIndexSet<P> &output_direct_keys,
    const CKKSLinearTransformStage<P> &stage, int k_step,
    int direct_giant_popcount_threshold = 4)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    int max_j2 = 0;
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        max_j2 = std::max(max_j2, r / k_step);
    }

    std::vector<bool> baby_used(static_cast<std::size_t>(k_step), false);
    std::vector<bool> giant_used(static_cast<std::size_t>(max_j2 + 1), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        baby_used[static_cast<std::size_t>(r % k_step)] = true;
        giant_used[static_cast<std::size_t>(r / k_step)] = true;
    }

    for (int j1 = 1; j1 < k_step; j1++) {
        if (baby_used[static_cast<std::size_t>(j1)])
            CKKSMarkRotationPowerKeyIndices<P>(input_binary_keys, j1);
    }
    for (int j2 = 1; j2 <= max_j2; j2++) {
        if (!giant_used[static_cast<std::size_t>(j2)]) continue;
        const int steps = k_step * j2;
        if (ckks_detail::popcountNonnegativeInt(steps) >=
            static_cast<std::size_t>(direct_giant_popcount_threshold))
            CKKSMarkDirectRotationKeyIndex<P>(output_direct_keys, steps);
        else
            CKKSMarkRotationPowerKeyIndices<P>(output_binary_keys, steps);
    }
}

template <class P>
inline std::size_t CKKSLinearTransformStageBabyRotationCount(
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    std::vector<bool> used(static_cast<std::size_t>(k_step), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        used[static_cast<std::size_t>(r % k_step)] = true;
    }

    std::size_t count = 0;
    for (int j1 = 1; j1 < k_step; j1++)
        if (used[static_cast<std::size_t>(j1)]) count++;
    return count;
}

template <class P>
inline std::size_t CKKSLinearTransformStageBabyRotationTableEvalAutoCount(
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    std::vector<bool> used(static_cast<std::size_t>(k_step), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        used[static_cast<std::size_t>(r % k_step)] = true;
    }

    std::vector<bool> ready(static_cast<std::size_t>(k_step), false);
    ready[0] = true;
    std::vector<int> stack;
    std::size_t count = 0;
    for (int target = 1; target < k_step; target++) {
        if (!used[static_cast<std::size_t>(target)]) continue;

        stack.clear();
        int step = target;
        while (!ready[static_cast<std::size_t>(step)]) {
            stack.push_back(step);
            step -= ckks_detail::highestPowerOfTwoLE(step);
        }

        while (!stack.empty()) {
            const int current = stack.back();
            stack.pop_back();
            if (ready[static_cast<std::size_t>(current)]) continue;
            ready[static_cast<std::size_t>(current)] = true;
            count++;
        }
    }
    return count;
}

template <class P>
inline std::size_t CKKSLinearTransformStageBabyRotationBinaryEvalAutoCount(
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    std::vector<bool> used(static_cast<std::size_t>(k_step), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        used[static_cast<std::size_t>(r % k_step)] = true;
    }

    std::size_t count = 0;
    for (int j1 = 1; j1 < k_step; j1++)
        if (used[static_cast<std::size_t>(j1)])
            count += ckks_detail::popcountNonnegativeInt(j1);
    return count;
}

template <class P>
inline std::size_t CKKSLinearTransformStageGiantRotationCount(
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    int max_j2 = 0;
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        max_j2 = std::max(max_j2, r / k_step);
    }

    std::vector<bool> used(static_cast<std::size_t>(max_j2 + 1), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        used[static_cast<std::size_t>(r / k_step)] = true;
    }

    std::size_t count = 0;
    for (int j2 = 1; j2 <= max_j2; j2++)
        if (used[static_cast<std::size_t>(j2)]) count++;
    return count;
}

template <class P>
inline std::size_t CKKSLinearTransformStageGiantRotationBinaryEvalAutoCount(
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    int max_j2 = 0;
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        max_j2 = std::max(max_j2, r / k_step);
    }

    std::vector<bool> used(static_cast<std::size_t>(max_j2 + 1), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        used[static_cast<std::size_t>(r / k_step)] = true;
    }

    std::size_t count = 0;
    for (int j2 = 1; j2 <= max_j2; j2++) {
        if (used[static_cast<std::size_t>(j2)])
            count += ckks_detail::popcountNonnegativeInt(k_step * j2);
    }
    return count;
}

template <class P>
inline std::size_t CKKSLinearTransformStageRotationEvalAutoCount(
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    return CKKSLinearTransformStageBabyRotationTableEvalAutoCount<P>(
               stage, k_step) +
           CKKSLinearTransformStageGiantRotationBinaryEvalAutoCount<P>(
               stage, k_step);
}

template <class P>
inline std::size_t CKKSLinearTransformStageDirectRotationEvalAutoCount(
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    return CKKSLinearTransformStageBabyRotationCount<P>(stage, k_step) +
           CKKSLinearTransformStageGiantRotationCount<P>(stage, k_step);
}

template <class P>
inline std::size_t CKKSLinearTransformStageHybridGiantRotationEvalAutoCount(
    const CKKSLinearTransformStage<P> &stage, int k_step,
    int direct_giant_popcount_threshold = 4)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(stage.diagonals.size() == stage.rotation_offsets.size());
    assert(!stage.rotation_offsets.empty());
    assert(k_step > 0 && k_step <= half);

    int max_j2 = 0;
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        max_j2 = std::max(max_j2, r / k_step);
    }

    std::vector<bool> giant_used(static_cast<std::size_t>(max_j2 + 1), false);
    for (const int offset : stage.rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        giant_used[static_cast<std::size_t>(r / k_step)] = true;
    }

    std::size_t giant_evalautos = 0;
    for (int j2 = 1; j2 <= max_j2; j2++) {
        if (!giant_used[static_cast<std::size_t>(j2)]) continue;
        const int steps = k_step * j2;
        const std::size_t binary_evalautos =
            ckks_detail::popcountNonnegativeInt(steps);
        giant_evalautos +=
            binary_evalautos >=
                    static_cast<std::size_t>(direct_giant_popcount_threshold)
                ? 1
                : binary_evalautos;
    }

    return CKKSLinearTransformStageBabyRotationTableEvalAutoCount<P>(
               stage, k_step) +
           giant_evalautos;
}

template <class P>
inline std::size_t CKKSLinearTransformStagesRotationEvalAutoCount(
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    std::size_t stage_count, int k_step)
{
    assert(first_stage + stage_count <= stages.size());
    std::size_t count = 0;
    for (std::size_t i = 0; i < stage_count; i++)
        count += CKKSLinearTransformStageRotationEvalAutoCount<P>(
            stages[first_stage + i], k_step);
    return count;
}

template <class P>
inline std::size_t CKKSLinearTransformStagesDirectRotationEvalAutoCount(
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    std::size_t stage_count, int k_step)
{
    assert(first_stage + stage_count <= stages.size());
    std::size_t count = 0;
    for (std::size_t i = 0; i < stage_count; i++)
        count += CKKSLinearTransformStageDirectRotationEvalAutoCount<P>(
            stages[first_stage + i], k_step);
    return count;
}

template <class P>
inline std::size_t CKKSLinearTransformStagesHybridGiantRotationEvalAutoCount(
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    std::size_t stage_count, int k_step,
    int direct_giant_popcount_threshold = 4)
{
    assert(first_stage + stage_count <= stages.size());
    std::size_t count = 0;
    for (std::size_t i = 0; i < stage_count; i++)
        count += CKKSLinearTransformStageHybridGiantRotationEvalAutoCount<P>(
            stages[first_stage + i], k_step,
            direct_giant_popcount_threshold);
    return count;
}

template <class P>
inline std::size_t
CKKSLinearTransformStagesDualInputSharedTailRotationEvalAutoCount(
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    std::size_t stage_count, int k_step)
{
    assert(stage_count > 0);
    assert(first_stage + stage_count <= stages.size());
    std::size_t count =
        2 * CKKSLinearTransformStageBabyRotationTableEvalAutoCount<P>(
                stages[first_stage], k_step) +
        CKKSLinearTransformStageGiantRotationBinaryEvalAutoCount<P>(
            stages[first_stage], k_step);
    for (std::size_t i = 1; i < stage_count; i++)
        count += CKKSLinearTransformStageRotationEvalAutoCount<P>(
            stages[first_stage + i], k_step);
    return count;
}

template <class P>
inline std::size_t
CKKSLinearTransformStagesDualInputSharedTailDirectRotationEvalAutoCount(
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    std::size_t stage_count, int k_step)
{
    assert(stage_count > 0);
    assert(first_stage + stage_count <= stages.size());
    std::size_t count =
        2 * CKKSLinearTransformStageBabyRotationCount<P>(
                stages[first_stage], k_step) +
        CKKSLinearTransformStageGiantRotationCount<P>(stages[first_stage],
                                                     k_step);
    for (std::size_t i = 1; i < stage_count; i++)
        count += CKKSLinearTransformStageDirectRotationEvalAutoCount<P>(
            stages[first_stage + i], k_step);
    return count;
}

template <class P>
inline std::size_t
CKKSLinearTransformStagesDualInputSharedTailHybridGiantRotationEvalAutoCount(
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    std::size_t stage_count, int k_step,
    int direct_giant_popcount_threshold = 4)
{
    assert(stage_count > 0);
    assert(first_stage + stage_count <= stages.size());
    std::size_t count =
        2 * CKKSLinearTransformStageBabyRotationTableEvalAutoCount<P>(
                stages[first_stage], k_step) +
        (CKKSLinearTransformStageHybridGiantRotationEvalAutoCount<P>(
             stages[first_stage], k_step, direct_giant_popcount_threshold) -
         CKKSLinearTransformStageBabyRotationTableEvalAutoCount<P>(
             stages[first_stage], k_step));
    for (std::size_t i = 1; i < stage_count; i++)
        count += CKKSLinearTransformStageHybridGiantRotationEvalAutoCount<P>(
            stages[first_stage + i], k_step,
            direct_giant_popcount_threshold);
    return count;
}

template <class P>
inline void CKKSComposeLinearTransformStages(
    CKKSLinearTransformStage<P> &out, const CKKSLinearTransformStage<P> &first,
    const CKKSLinearTransformStage<P> &second)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(first.diagonals.size() == first.rotation_offsets.size());
    assert(second.diagonals.size() == second.rotation_offsets.size());

    CKKSLinearTransformStage<P> composed;
    for (std::size_t b = 0; b < second.diagonals.size(); b++) {
        const int second_offset =
            ((second.rotation_offsets[b] % half) + half) % half;
        for (std::size_t a = 0; a < first.diagonals.size(); a++) {
            const int first_offset =
                ((first.rotation_offsets[a] % half) + half) % half;
            const int offset = (second_offset + first_offset) % half;
            auto it = std::find(composed.rotation_offsets.begin(),
                                composed.rotation_offsets.end(), offset);
            if (it == composed.rotation_offsets.end()) {
                composed.rotation_offsets.push_back(offset);
                composed.diagonals.emplace_back();
                composed.diagonals.back().fill({0.0, 0.0});
                it = composed.rotation_offsets.end() - 1;
            }
            auto &diagonal = composed.diagonals[static_cast<std::size_t>(
                it - composed.rotation_offsets.begin())];
            for (int i = 0; i < half; i++) {
                const int shifted = (i + second_offset) % half;
                diagonal[i] += second.diagonals[b][i] *
                               first.diagonals[a][shifted];
            }
        }
    }
    out = std::move(composed);
}

template <class P>
inline void CKKSFuseLinearTransformStages(CKKSLinearTransformStages<P> &out,
                                          const CKKSLinearTransformStages<P> &in,
                                          std::size_t stages_per_fused_level)
{
    assert(stages_per_fused_level > 0);
    out.clear();
    out.reserve((in.size() + stages_per_fused_level - 1) /
                stages_per_fused_level);
    for (std::size_t start = 0; start < in.size();
         start += stages_per_fused_level) {
        CKKSLinearTransformStage<P> fused = in[start];
        const std::size_t end =
            std::min(in.size(), start + stages_per_fused_level);
        for (std::size_t i = start + 1; i < end; i++) {
            CKKSLinearTransformStage<P> next;
            CKKSComposeLinearTransformStages<P>(next, fused, in[i]);
            fused = std::move(next);
        }
        out.push_back(std::move(fused));
    }
}

template <class P>
inline void CKKSScaleLinearTransformStages(CKKSLinearTransformStages<P> &stages,
                                           std::complex<double> scalar)
{
    for (auto &stage : stages) {
        if (stage.diagonals.empty()) continue;
        for (auto &diagonal : stage.diagonals)
            for (auto &entry : diagonal) entry *= scalar;
        return;
    }
}

namespace ckks_detail {

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          class IndexSequence>
struct CKKSGaloisKeyChainTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t... Is>
struct CKKSGaloisKeyChainTuple<P, StartLogQ, PlainLogDelta,
                               std::index_sequence<Is...>> {
    using type = std::tuple<CKKSGaloisKey<P, StartLogQ - Is * PlainLogDelta>...>;
};

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
struct CKKSGaloisKeyChain {
    static_assert(StartLogQ > StageCount * PlainLogDelta);

    using Tuple = typename ckks_detail::CKKSGaloisKeyChainTuple<
        P, StartLogQ, PlainLogDelta,
        std::make_index_sequence<StageCount + 1>>::type;

    Tuple keys{};

    template <std::size_t I>
    auto &get()
    {
        return std::get<I>(keys);
    }

    template <std::size_t I>
    const auto &get() const
    {
        return std::get<I>(keys);
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
inline void CKKSGaloisKeyChainGenImpl(
    CKKSGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount> &chain,
    const Key<P> &key, CKKSNoise noise)
{
    if constexpr (I <= StageCount) {
        constexpr std::uint32_t log_q = StartLogQ - I * PlainLogDelta;
        CKKSGaloisKeyGen<P, log_q>(chain.template get<I>(), key, noise);
        CKKSGaloisKeyChainGenImpl<I + 1, P, StartLogQ, PlainLogDelta,
                                  StageCount>(chain, key, noise);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
inline void CKKSGaloisKeyChainGen(
    CKKSGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount> &chain,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    ckks_detail::CKKSGaloisKeyChainGenImpl<0, P, StartLogQ, PlainLogDelta,
                                           StageCount>(chain, key, noise);
}

namespace ckks_detail {

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          class IndexSequence>
struct CKKSSparseGaloisKeyChainTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t... Is>
struct CKKSSparseGaloisKeyChainTuple<P, StartLogQ, PlainLogDelta,
                                     std::index_sequence<Is...>> {
    using type = std::tuple<
        CKKSSparseGaloisKey<P, StartLogQ - Is * PlainLogDelta>...>;
};

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
struct CKKSSparseGaloisKeyChain {
    static_assert(StartLogQ > StageCount * PlainLogDelta);

    using Tuple = typename ckks_detail::CKKSSparseGaloisKeyChainTuple<
        P, StartLogQ, PlainLogDelta,
        std::make_index_sequence<StageCount + 1>>::type;

    Tuple keys{};

    template <std::size_t I>
    auto &get()
    {
        return std::get<I>(keys);
    }

    template <std::size_t I>
    const auto &get() const
    {
        return std::get<I>(keys);
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
inline void CKKSSparseGaloisKeyChainGenImpl(
    CKKSSparseGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount> &chain,
    const std::array<CKKSRotationKeyIndexSet<P>, StageCount + 1> &usage,
    const Key<P> &key, CKKSNoise noise)
{
    if constexpr (I <= StageCount) {
        constexpr std::uint32_t log_q = StartLogQ - I * PlainLogDelta;
        CKKSSparseGaloisKeyGen<P, log_q>(chain.template get<I>(), key,
                                         usage[I], noise);
        CKKSSparseGaloisKeyChainGenImpl<I + 1, P, StartLogQ, PlainLogDelta,
                                        StageCount>(chain, usage, key, noise);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
inline void CKKSSparseGaloisKeyChainGen(
    CKKSSparseGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount> &chain,
    const std::array<CKKSRotationKeyIndexSet<P>, StageCount + 1> &usage,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    ckks_detail::CKKSSparseGaloisKeyChainGenImpl<0, P, StartLogQ,
                                                 PlainLogDelta, StageCount>(
        chain, usage, key, noise);
}

namespace ckks_detail {

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          class IndexSequence>
struct CKKSDirectSparseGaloisKeyChainTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t... Is>
struct CKKSDirectSparseGaloisKeyChainTuple<P, StartLogQ, PlainLogDelta,
                                           std::index_sequence<Is...>> {
    using type = std::tuple<
        CKKSDirectSparseGaloisKey<P, StartLogQ - Is * PlainLogDelta>...>;
};

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
struct CKKSDirectSparseGaloisKeyChain {
    static_assert(StartLogQ > StageCount * PlainLogDelta);

    using Tuple = typename ckks_detail::CKKSDirectSparseGaloisKeyChainTuple<
        P, StartLogQ, PlainLogDelta,
        std::make_index_sequence<StageCount + 1>>::type;

    Tuple keys{};

    template <std::size_t I>
    auto &get()
    {
        return std::get<I>(keys);
    }

    template <std::size_t I>
    const auto &get() const
    {
        return std::get<I>(keys);
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
inline void CKKSDirectSparseGaloisKeyChainGenImpl(
    CKKSDirectSparseGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount>
        &chain,
    const std::array<CKKSDirectRotationKeyIndexSet<P>, StageCount + 1> &usage,
    const Key<P> &key, CKKSNoise noise)
{
    if constexpr (I <= StageCount) {
        constexpr std::uint32_t log_q = StartLogQ - I * PlainLogDelta;
        CKKSDirectSparseGaloisKeyGen<P, log_q>(chain.template get<I>(), key,
                                               usage[I], noise);
        CKKSDirectSparseGaloisKeyChainGenImpl<I + 1, P, StartLogQ,
                                              PlainLogDelta, StageCount>(
            chain, usage, key, noise);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
inline void CKKSDirectSparseGaloisKeyChainGen(
    CKKSDirectSparseGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount>
        &chain,
    const std::array<CKKSDirectRotationKeyIndexSet<P>, StageCount + 1> &usage,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    ckks_detail::CKKSDirectSparseGaloisKeyChainGenImpl<
        0, P, StartLogQ, PlainLogDelta, StageCount>(chain, usage, key, noise);
}

namespace ckks_detail {

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          class IndexSequence>
struct CKKSHybridSparseGaloisKeyChainTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t... Is>
struct CKKSHybridSparseGaloisKeyChainTuple<P, StartLogQ, PlainLogDelta,
                                           std::index_sequence<Is...>> {
    using type = std::tuple<
        CKKSHybridSparseGaloisKey<P, StartLogQ - Is * PlainLogDelta>...>;
};

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
struct CKKSHybridSparseGaloisKeyChain {
    static_assert(StartLogQ > StageCount * PlainLogDelta);

    using Tuple = typename ckks_detail::CKKSHybridSparseGaloisKeyChainTuple<
        P, StartLogQ, PlainLogDelta,
        std::make_index_sequence<StageCount + 1>>::type;

    Tuple keys{};

    template <std::size_t I>
    auto &get()
    {
        return std::get<I>(keys);
    }

    template <std::size_t I>
    const auto &get() const
    {
        return std::get<I>(keys);
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
inline void CKKSHybridSparseGaloisKeyChainGenImpl(
    CKKSHybridSparseGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount>
        &chain,
    const std::array<CKKSRotationKeyIndexSet<P>, StageCount + 1>
        &binary_usage,
    const std::array<CKKSDirectRotationKeyIndexSet<P>, StageCount + 1>
        &direct_usage,
    const Key<P> &key, CKKSNoise noise)
{
    if constexpr (I <= StageCount) {
        constexpr std::uint32_t log_q = StartLogQ - I * PlainLogDelta;
        CKKSHybridSparseGaloisKeyGen<P, log_q>(
            chain.template get<I>(), key, binary_usage[I], direct_usage[I],
            noise);
        CKKSHybridSparseGaloisKeyChainGenImpl<I + 1, P, StartLogQ,
                                              PlainLogDelta, StageCount>(
            chain, binary_usage, direct_usage, key, noise);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t PlainLogDelta,
          std::size_t StageCount>
inline void CKKSHybridSparseGaloisKeyChainGen(
    CKKSHybridSparseGaloisKeyChain<P, StartLogQ, PlainLogDelta, StageCount>
        &chain,
    const std::array<CKKSRotationKeyIndexSet<P>, StageCount + 1>
        &binary_usage,
    const std::array<CKKSDirectRotationKeyIndexSet<P>, StageCount + 1>
        &direct_usage,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    ckks_detail::CKKSHybridSparseGaloisKeyChainGenImpl<
        0, P, StartLogQ, PlainLogDelta, StageCount>(
        chain, binary_usage, direct_usage, key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildLinearTransformBSGSPlan(
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const std::vector<CKKSSlotVector<P>> &diagonals,
    const std::vector<int> &rotation_offsets, int k_step)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(diagonals.size() == rotation_offsets.size());
    assert(!diagonals.empty());
    assert(k_step > 0 && k_step <= half);

    int max_j2 = 0;
    for (const int offset : rotation_offsets) {
        const int r = ((offset % half) + half) % half;
        max_j2 = std::max(max_j2, r / k_step);
    }

    std::vector<std::vector<CKKSSlotVector<P>>> accumulated(
        max_j2 + 1, std::vector<CKKSSlotVector<P>>(k_step));
    std::vector<std::vector<bool>> used(max_j2 + 1,
                                        std::vector<bool>(k_step, false));
    for (auto &group : accumulated)
        for (auto &diag : group) diag.fill({0.0, 0.0});

    auto rotated = std::make_unique<CKKSSlotVector<P>>();
    for (std::size_t i = 0; i < rotation_offsets.size(); i++) {
        const int r = ((rotation_offsets[i] % half) + half) % half;
        const int j2 = r / k_step;
        const int j1 = r % k_step;
        rotateCKKSSlotVector<P>(*rotated, diagonals[i], -k_step * j2);
        for (int s = 0; s < half; s++) accumulated[j2][j1][s] += (*rotated)[s];
        used[j2][j1] = true;
    }

    plan.k_step = k_step;
    plan.groups.clear();
    plan.groups.reserve(max_j2 + 1);
    for (int j2 = 0; j2 <= max_j2; j2++) {
        CKKSLinearTransformGiantStep<P, LogQ, PlainLogDelta> group;
        group.giant_step = j2;
        for (int j1 = 0; j1 < k_step; j1++) {
            if (!used[j2][j1]) continue;
            CKKSLinearTransformTerm<P, LogQ, PlainLogDelta> term;
            term.baby_step = j1;
            ckksSlotEncode<P, LogQ, PlainLogDelta>(term.plain.poly,
                                                   accumulated[j2][j1]);
            group.terms.push_back(std::move(term));
        }
        if (!group.terms.empty()) plan.groups.push_back(std::move(group));
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildLinearTransformBSGSPlan(
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const CKKSLinearTransformStage<P> &stage, int k_step)
{
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, stage.diagonals, stage.rotation_offsets, k_step);
}

namespace ckks_detail {

template <class GaloisKey>
inline void maybe_clear_cached_keys(const GaloisKey &gk)
{
    if constexpr (requires { gk.clear_cached_keys(); }) {
        gk.clear_cached_keys();
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSLinearTransformBSGS(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const InputGaloisKey &input_gk, const OutputGaloisKey &output_gk)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    assert(plan.k_step > 0 && plan.k_step <= static_cast<int>(P::n / 2));
    assert(!plan.groups.empty());
    const bool trace_linear = ckks_detail::ckks_linear_trace_enabled();
    const auto trace_total_start = std::chrono::steady_clock::now();
    auto trace_section_start = trace_total_start;

    std::vector<bool> baby_used;
    CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta, PlainLogDelta>(
        baby_used, plan);

    auto baby =
        std::make_unique<ckks_detail::CKKSSparseBabyStepTable<P>>(plan.k_step);
    ckks_detail::CKKSBuildSparseBabyStepRotationTable<P, LogQ>(
        *baby, ct.ct, baby_used, input_gk);
    ckks_detail::maybe_clear_cached_keys(input_gk);
    if (trace_linear) {
        ckks_detail::ckks_trace_linear_event(
            "baby_rotation_table",
            ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
        trace_section_start = std::chrono::steady_clock::now();
    }

    if constexpr (is_multilimb_uint_v<typename P::T> &&
                  use_multilimb_digit_fft_v<P>) {
        using T = typename P::T;
        constexpr int width = std::numeric_limits<T>::digits;
        constexpr int plain_digits =
            (64 + static_cast<int>(P::B̅gbit) - 1) /
            static_cast<int>(P::B̅gbit);
        constexpr uint64_t plain_mask =
            P::B̅gbit == 64 ? std::numeric_limits<uint64_t>::max()
                            : ((uint64_t{1} << P::B̅gbit) - 1);
        constexpr T half_digit = T{1} << (P::B̅gbit - 1);
        constexpr T torus_mask = (T{1} << P::B̅gbit) - T{1};
        constexpr T offset = [] {
            constexpr int local_width = std::numeric_limits<T>::digits;
            constexpr T local_half_digit = T{1} << (P::B̅gbit - 1);
            T value = 0;
            for (int j = 0; j < static_cast<int>(P::l̅); j++)
                value += local_half_digit
                         << (local_width - (j + 1) * P::B̅gbit);
            return value;
        }();
        constexpr int fd_margin =
            std::numeric_limits<double>::digits -
            (2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3);
        static_assert(fd_margin > 0,
                      "fused CKKS BSGS digit products must fit in double");
        constexpr int max_fused_terms = 1 << (fd_margin - 1);

        using BabyComponentFFT = std::array<PolynomialInFD<P>, P::l̅>;
        using BabyFFT = std::array<BabyComponentFFT, P::k + 1>;
        auto baby_fft =
            std::make_unique<std::vector<std::unique_ptr<BabyFFT>>>(
                plan.k_step);
        for (int j1 = 0; j1 < static_cast<int>(baby_used.size()); j1++) {
            if (!baby_used[static_cast<std::size_t>(j1)]) continue;
            (*baby_fft)[static_cast<std::size_t>(j1)] =
                std::make_unique<BabyFFT>();
        }
#pragma omp parallel
        {
            auto torus_digit = std::make_unique<Polynomial<P>>();
#pragma omp for schedule(dynamic)
            for (int j1 = 0; j1 < static_cast<int>(baby_used.size()); j1++) {
                if (!baby_used[static_cast<std::size_t>(j1)]) continue;
                const auto &baby_ct =
                    *(*baby)[static_cast<std::size_t>(j1)];
                auto &baby_out =
                    *(*baby_fft)[static_cast<std::size_t>(j1)];
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                        const int torus_shift =
                            width - (j + 1) * static_cast<int>(P::B̅gbit);
                        for (uint32_t n = 0; n < P::n; n++) {
                            const T adjusted = baby_ct[c][n] + offset;
                            (*torus_digit)[n] =
                                ((adjusted >> torus_shift) & torus_mask) -
                                half_digit;
                        }
                        TwistIFFTDigit<P>(baby_out[c][j], *torus_digit);
                    }
                }
            }
        }
        baby.reset();
        if (trace_linear) {
            ckks_detail::ckks_trace_linear_event(
                "baby_digit_fft",
                ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
            trace_section_start = std::chrono::steady_clock::now();
        }

        bool res_initialized = false;
        auto inner = std::make_unique<TRLWE<P>>();
        auto giant_rotated = std::make_unique<TRLWE<P>>();

        std::size_t group_index = 0;
        for (const auto &group : plan.groups) {
            using PlainTermFFT = std::array<PolynomialInFD<P>, plain_digits>;
            auto plain_fft =
                std::make_unique<std::vector<PlainTermFFT>>(group.terms.size());
            const int term_count = static_cast<int>(group.terms.size());
#pragma omp parallel
            {
                auto plain_digit = std::make_unique<Polynomial<P>>();
#pragma omp for collapse(2) schedule(dynamic)
                for (int t = 0; t < term_count; t++) {
                    for (int d = 0; d < plain_digits; d++) {
                        const auto &plain =
                            group.terms[static_cast<std::size_t>(t)]
                                .plain.poly;
                        const int plain_shift =
                            d * static_cast<int>(P::B̅gbit);
                        for (uint32_t n = 0; n < P::n; n++) {
                            const auto [negative, magnitude] =
                                ckks_detail::smallSignedMagnitude<P, LogQ>(
                                    plain[n]);
                            const uint64_t digit =
                                (magnitude >> plain_shift) & plain_mask;
                            const T digit_value{digit};
                            (*plain_digit)[n] =
                                negative ? T{0} - digit_value : digit_value;
                        }
                        TwistIFFTDigit<P>(
                            (*plain_fft)[static_cast<std::size_t>(t)][d],
                            *plain_digit);
                    }
                }
            }
            if (trace_linear) {
                ckks_detail::ckks_trace_linear_group_event(
                    "plain_digit_fft", group_index,
                    ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
                trace_section_start = std::chrono::steady_clock::now();
            }

            for (int c = 0; c <= static_cast<int>(P::k); c++)
                for (uint32_t n = 0; n < P::n; n++) (*inner)[c][n] = 0;

            const int component_count = static_cast<int>(P::k) + 1;
            const int aux_row_count = static_cast<int>(P::l̅);
#pragma omp parallel
            {
                auto local_inner = std::make_unique<TRLWE<P>>();
                auto product_acc = std::make_unique<PolynomialInFD<P>>();
                auto digit_prod = std::make_unique<Polynomial<P>>();

#pragma omp for collapse(3) schedule(dynamic)
                for (int c = 0; c < component_count; c++) {
                    for (int j = 0; j < aux_row_count; j++) {
                        for (int d = 0; d < plain_digits; d++) {
                            const int torus_shift =
                                width -
                                (j + 1) * static_cast<int>(P::B̅gbit);
                            const int plain_shift =
                                d * static_cast<int>(P::B̅gbit);
                            const int shift = torus_shift + plain_shift;
                            if (shift >= width) continue;

                            for (std::size_t batch_begin = 0;
                                 batch_begin < group.terms.size();
                                 batch_begin += max_fused_terms) {
                                const std::size_t batch_end =
                                    std::min<std::size_t>(
                                        group.terms.size(),
                                        batch_begin + max_fused_terms);
                                const int first_j1 =
                                    group.terms[batch_begin].baby_step;
                                assert(first_j1 >= 0 &&
                                       first_j1 < plan.k_step);
                                const auto &first_fft =
                                    *(*baby_fft)[static_cast<std::size_t>(
                                        first_j1)];
                                MulInFD<P::n>(
                                    *product_acc, first_fft[c][j],
                                    (*plain_fft)[batch_begin][d]);
                                for (std::size_t t = batch_begin + 1;
                                     t < batch_end; t++) {
                                    const int j1 = group.terms[t].baby_step;
                                    assert(j1 >= 0 && j1 < plan.k_step);
                                    const auto &fft =
                                        *(*baby_fft)[static_cast<std::size_t>(
                                            j1)];
                                    FMAInFD<P::n>(*product_acc, fft[c][j],
                                                  (*plain_fft)[t][d]);
                                }
                                TwistFFTDigitProduct<P>(*digit_prod,
                                                        *product_acc);
                                for (uint32_t n = 0; n < P::n; n++)
                                    (*local_inner)[c][n] +=
                                        (*digit_prod)[n] << shift;
                            }
                        }
                    }
                }
#pragma omp critical
                {
                    for (int c = 0; c < component_count; c++)
                        for (uint32_t n = 0; n < P::n; n++)
                            (*inner)[c][n] += (*local_inner)[c][n];
                }
            }
            if (trace_linear) {
                ckks_detail::ckks_trace_linear_group_event(
                    "fused_digit_products", group_index,
                    ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
                trace_section_start = std::chrono::steady_clock::now();
            }

            for (int c = 0; c <= static_cast<int>(P::k); c++)
                for (uint32_t n = 0; n < P::n; n++)
                    (*inner)[c][n] =
                        ckks_detail::roundedLevelRightShift<
                            P, LogQ, PlainLogDelta>((*inner)[c][n]);
            ckks_detail::reduceTRLWEToLevel<P, out_log_q>(*inner);

            const TRLWE<P> *to_add = inner.get();
            if (group.giant_step != 0) {
                CKKSRotateSlots<P, out_log_q>(
                    *giant_rotated, *inner, plan.k_step * group.giant_step,
                    output_gk);
                to_add = giant_rotated.get();
            }

            if (!res_initialized) {
                res.ct = *to_add;
                res_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, *to_add);
            }
            if (trace_linear) {
                ckks_detail::ckks_trace_linear_group_event(
                    "group_finalize", group_index,
                    ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
                trace_section_start = std::chrono::steady_clock::now();
            }
            group_index++;
        }
        ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
        if (trace_linear)
            ckks_detail::ckks_trace_linear_event(
                "linear_transform_total",
                ckks_detail::ckks_trace_elapsed_ms(trace_total_start));
        return;
    }

    bool res_initialized = false;
    auto term = std::make_unique<TRLWE<P>>();
    auto inner = std::make_unique<TRLWE<P>>();
    auto giant_rotated = std::make_unique<TRLWE<P>>();

    for (const auto &group : plan.groups) {
        bool inner_initialized = false;
        for (const auto &entry : group.terms) {
            assert((*baby)[static_cast<std::size_t>(entry.baby_step)]);
            CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(
                *term, *(*baby)[static_cast<std::size_t>(entry.baby_step)],
                entry.plain.poly);

            if (!inner_initialized) {
                *inner = *term;
                inner_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(*inner, *term);
            }
        }

        const TRLWE<P> *to_add = inner.get();
        if (group.giant_step != 0) {
            CKKSRotateSlots<P, out_log_q>(
                *giant_rotated, *inner, plan.k_step * group.giant_step,
                output_gk);
            to_add = giant_rotated.get();
        }

        if (!res_initialized) {
            res.ct = *to_add;
            res_initialized = true;
        }
        else {
            CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, *to_add);
        }
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSLinearTransformBSGSDirect(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const InputGaloisKey &input_gk, const OutputGaloisKey &output_gk)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    assert(plan.k_step > 0 && plan.k_step <= static_cast<int>(P::n / 2));
    assert(!plan.groups.empty());

    std::vector<bool> baby_used;
    CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta, PlainLogDelta>(
        baby_used, plan);

    auto baby =
        std::make_unique<ckks_detail::CKKSSparseBabyStepTable<P>>(plan.k_step);
    ckks_detail::CKKSBuildSparseDirectBabyStepRotationTable<P, LogQ>(
        *baby, ct.ct, baby_used, input_gk);
    ckks_detail::maybe_clear_cached_keys(input_gk);

    bool res_initialized = false;
    auto term = std::make_unique<TRLWE<P>>();
    auto inner = std::make_unique<TRLWE<P>>();
    auto giant_rotated = std::make_unique<TRLWE<P>>();

    for (const auto &group : plan.groups) {
        bool inner_initialized = false;
        for (const auto &entry : group.terms) {
            assert((*baby)[static_cast<std::size_t>(entry.baby_step)]);
            CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(
                *term, *(*baby)[static_cast<std::size_t>(entry.baby_step)],
                entry.plain.poly);

            if (!inner_initialized) {
                *inner = *term;
                inner_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(*inner, *term);
            }
        }

        const TRLWE<P> *to_add = inner.get();
        if (group.giant_step != 0) {
            CKKSRotateSlots<P, out_log_q>(
                *giant_rotated, *inner, plan.k_step * group.giant_step,
                output_gk);
            to_add = giant_rotated.get();
        }

        if (!res_initialized) {
            res.ct = *to_add;
            res_initialized = true;
        }
        else {
            CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, *to_add);
        }
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSLinearTransformBSGSDualInput(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &lhs,
    const CKKSCiphertext<P, LogQ, LogDelta> &rhs,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &lhs_plan,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &rhs_plan,
    const InputGaloisKey &input_gk, const OutputGaloisKey &output_gk)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    assert(lhs_plan.k_step > 0 &&
           lhs_plan.k_step <= static_cast<int>(P::n / 2));
    assert(rhs_plan.k_step == lhs_plan.k_step);
    assert(!lhs_plan.groups.empty());
    assert(lhs_plan.groups.size() == rhs_plan.groups.size());

    std::vector<bool> lhs_baby_used;
    std::vector<bool> rhs_baby_used;
    CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta, PlainLogDelta>(
        lhs_baby_used, lhs_plan);
    CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta, PlainLogDelta>(
        rhs_baby_used, rhs_plan);
    for (std::size_t i = 0; i < lhs_baby_used.size(); i++)
        lhs_baby_used[i] = lhs_baby_used[i] || rhs_baby_used[i];

    auto lhs_baby =
        std::make_unique<ckks_detail::CKKSSparseBabyStepTable<P>>(
            lhs_plan.k_step);
    auto rhs_baby =
        std::make_unique<ckks_detail::CKKSSparseBabyStepTable<P>>(
            lhs_plan.k_step);
    ckks_detail::CKKSBuildSparseBabyStepRotationTable<P, LogQ>(
        *lhs_baby, lhs.ct, lhs_baby_used, input_gk);
    ckks_detail::CKKSBuildSparseBabyStepRotationTable<P, LogQ>(
        *rhs_baby, rhs.ct, lhs_baby_used, input_gk);
    ckks_detail::maybe_clear_cached_keys(input_gk);

    const bool trace_linear = ckks_detail::ckks_linear_trace_enabled();
    const auto trace_total_start = std::chrono::steady_clock::now();
    auto trace_section_start = trace_total_start;

    if (trace_linear) {
        ckks_detail::ckks_trace_linear_event(
            "dual_baby_rotation_table",
            ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
        trace_section_start = std::chrono::steady_clock::now();
    }

    if constexpr (is_multilimb_uint_v<typename P::T> &&
                  use_multilimb_digit_fft_v<P>) {
        using T = typename P::T;
        constexpr int width = std::numeric_limits<T>::digits;
        constexpr int plain_digits =
            (64 + static_cast<int>(P::B̅gbit) - 1) /
            static_cast<int>(P::B̅gbit);
        constexpr uint64_t plain_mask =
            P::B̅gbit == 64 ? std::numeric_limits<uint64_t>::max()
                            : ((uint64_t{1} << P::B̅gbit) - 1);
        constexpr T half_digit = T{1} << (P::B̅gbit - 1);
        constexpr T torus_mask = (T{1} << P::B̅gbit) - T{1};
        constexpr T offset = [] {
            constexpr int local_width = std::numeric_limits<T>::digits;
            constexpr T local_half_digit = T{1} << (P::B̅gbit - 1);
            T value = 0;
            for (int j = 0; j < static_cast<int>(P::l̅); j++)
                value += local_half_digit
                         << (local_width - (j + 1) * P::B̅gbit);
            return value;
        }();
        constexpr int fd_margin =
            std::numeric_limits<double>::digits -
            (2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3);
        static_assert(fd_margin > 0,
                      "fused CKKS dual-input BSGS digit products must fit in double");
        constexpr int max_fused_terms = 1 << (fd_margin - 1);

        using BabyComponentFFT = std::array<PolynomialInFD<P>, P::l̅>;
        using BabyFFT = std::array<BabyComponentFFT, P::k + 1>;
        auto lhs_baby_fft =
            std::make_unique<std::vector<std::unique_ptr<BabyFFT>>>(
                lhs_plan.k_step);
        auto rhs_baby_fft =
            std::make_unique<std::vector<std::unique_ptr<BabyFFT>>>(
                lhs_plan.k_step);
        for (int j1 = 0; j1 < static_cast<int>(lhs_baby_used.size()); j1++) {
            if (!lhs_baby_used[static_cast<std::size_t>(j1)]) continue;
            (*lhs_baby_fft)[static_cast<std::size_t>(j1)] =
                std::make_unique<BabyFFT>();
            (*rhs_baby_fft)[static_cast<std::size_t>(j1)] =
                std::make_unique<BabyFFT>();
        }
#pragma omp parallel
        {
            auto torus_digit = std::make_unique<Polynomial<P>>();
#pragma omp for schedule(dynamic)
            for (int j1 = 0; j1 < static_cast<int>(lhs_baby_used.size());
                 j1++) {
                if (!lhs_baby_used[static_cast<std::size_t>(j1)]) continue;
                const auto &lhs_baby_ct =
                    *(*lhs_baby)[static_cast<std::size_t>(j1)];
                const auto &rhs_baby_ct =
                    *(*rhs_baby)[static_cast<std::size_t>(j1)];
                auto &lhs_baby_out =
                    *(*lhs_baby_fft)[static_cast<std::size_t>(j1)];
                auto &rhs_baby_out =
                    *(*rhs_baby_fft)[static_cast<std::size_t>(j1)];
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    for (int j = 0; j < static_cast<int>(P::l̅); j++) {
                        const int torus_shift =
                            width - (j + 1) * static_cast<int>(P::B̅gbit);
                        for (uint32_t n = 0; n < P::n; n++) {
                            const T adjusted = lhs_baby_ct[c][n] + offset;
                            (*torus_digit)[n] =
                                ((adjusted >> torus_shift) & torus_mask) -
                                half_digit;
                        }
                        TwistIFFTDigit<P>(lhs_baby_out[c][j], *torus_digit);
                        for (uint32_t n = 0; n < P::n; n++) {
                            const T adjusted = rhs_baby_ct[c][n] + offset;
                            (*torus_digit)[n] =
                                ((adjusted >> torus_shift) & torus_mask) -
                                half_digit;
                        }
                        TwistIFFTDigit<P>(rhs_baby_out[c][j], *torus_digit);
                    }
                }
            }
        }
        lhs_baby.reset();
        rhs_baby.reset();
        if (trace_linear) {
            ckks_detail::ckks_trace_linear_event(
                "dual_baby_digit_fft",
                ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
            trace_section_start = std::chrono::steady_clock::now();
        }

        bool res_initialized = false;
        auto inner = std::make_unique<TRLWE<P>>();
        auto giant_rotated = std::make_unique<TRLWE<P>>();

        for (std::size_t group_index = 0;
             group_index < lhs_plan.groups.size(); group_index++) {
            const auto &lhs_group = lhs_plan.groups[group_index];
            const auto &rhs_group = rhs_plan.groups[group_index];
            assert(lhs_group.giant_step == rhs_group.giant_step);

            using PlainTermFFT = std::array<PolynomialInFD<P>, plain_digits>;
            auto lhs_plain_fft =
                std::make_unique<std::vector<PlainTermFFT>>(
                    lhs_group.terms.size());
            auto rhs_plain_fft =
                std::make_unique<std::vector<PlainTermFFT>>(
                    rhs_group.terms.size());
            const int lhs_term_count =
                static_cast<int>(lhs_group.terms.size());
            const int rhs_term_count =
                static_cast<int>(rhs_group.terms.size());
            const int total_term_count = lhs_term_count + rhs_term_count;
            assert(total_term_count > 0);
#pragma omp parallel
            {
                auto plain_digit = std::make_unique<Polynomial<P>>();
#pragma omp for collapse(2) schedule(dynamic)
                for (int t = 0; t < total_term_count; t++) {
                    for (int d = 0; d < plain_digits; d++) {
                        const auto &plain =
                            t < lhs_term_count
                                ? lhs_group.terms[static_cast<std::size_t>(t)]
                                      .plain.poly
                                : rhs_group
                                      .terms[static_cast<std::size_t>(
                                          t - lhs_term_count)]
                                      .plain.poly;
                        auto &plain_fft =
                            t < lhs_term_count
                                ? (*lhs_plain_fft)[static_cast<std::size_t>(t)]
                                : (*rhs_plain_fft)[static_cast<std::size_t>(
                                      t - lhs_term_count)];
                        const int plain_shift =
                            d * static_cast<int>(P::B̅gbit);
                        for (uint32_t n = 0; n < P::n; n++) {
                            const auto [negative, magnitude] =
                                ckks_detail::smallSignedMagnitude<P, LogQ>(
                                    plain[n]);
                            const uint64_t digit =
                                (magnitude >> plain_shift) & plain_mask;
                            const T digit_value{digit};
                            (*plain_digit)[n] =
                                negative ? T{0} - digit_value : digit_value;
                        }
                        TwistIFFTDigit<P>(plain_fft[d], *plain_digit);
                    }
                }
            }
            if (trace_linear) {
                ckks_detail::ckks_trace_linear_group_event(
                    "dual_plain_digit_fft", group_index,
                    ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
                trace_section_start = std::chrono::steady_clock::now();
            }

            for (int c = 0; c <= static_cast<int>(P::k); c++)
                for (uint32_t n = 0; n < P::n; n++) (*inner)[c][n] = 0;

            const int component_count = static_cast<int>(P::k) + 1;
            const int aux_row_count = static_cast<int>(P::l̅);
#pragma omp parallel
            {
                auto local_inner = std::make_unique<TRLWE<P>>();
                auto product_acc = std::make_unique<PolynomialInFD<P>>();
                auto digit_prod = std::make_unique<Polynomial<P>>();

#pragma omp for collapse(3) schedule(dynamic)
                for (int c = 0; c < component_count; c++) {
                    for (int j = 0; j < aux_row_count; j++) {
                        for (int d = 0; d < plain_digits; d++) {
                            const int torus_shift =
                                width -
                                (j + 1) * static_cast<int>(P::B̅gbit);
                            const int plain_shift =
                                d * static_cast<int>(P::B̅gbit);
                            const int shift = torus_shift + plain_shift;
                            if (shift >= width) continue;

                            for (std::size_t batch_begin = 0;
                                 batch_begin <
                                 static_cast<std::size_t>(total_term_count);
                                 batch_begin += max_fused_terms) {
                                const std::size_t batch_end =
                                    std::min<std::size_t>(
                                        static_cast<std::size_t>(
                                            total_term_count),
                                        batch_begin + max_fused_terms);
                                for (std::size_t t = batch_begin;
                                     t < batch_end; t++) {
                                    const bool is_lhs =
                                        t <
                                        static_cast<std::size_t>(
                                            lhs_term_count);
                                    const std::size_t local_t =
                                        is_lhs
                                            ? t
                                            : t - static_cast<std::size_t>(
                                                      lhs_term_count);
                                    const auto &entry =
                                        is_lhs
                                            ? lhs_group.terms[local_t]
                                            : rhs_group.terms[local_t];
                                    assert(entry.baby_step >= 0 &&
                                           entry.baby_step < lhs_plan.k_step);
                                    const auto &baby_fft =
                                        is_lhs
                                            ? *(*lhs_baby_fft)
                                                   [static_cast<std::size_t>(
                                                       entry.baby_step)]
                                            : *(*rhs_baby_fft)
                                                   [static_cast<std::size_t>(
                                                       entry.baby_step)];
                                    const auto &plain_fft =
                                        is_lhs ? (*lhs_plain_fft)[local_t]
                                               : (*rhs_plain_fft)[local_t];
                                    if (t == batch_begin)
                                        MulInFD<P::n>(*product_acc,
                                                      baby_fft[c][j],
                                                      plain_fft[d]);
                                    else
                                        FMAInFD<P::n>(*product_acc,
                                                     baby_fft[c][j],
                                                     plain_fft[d]);
                                }
                                TwistFFTDigitProduct<P>(*digit_prod,
                                                        *product_acc);
                                for (uint32_t n = 0; n < P::n; n++)
                                    (*local_inner)[c][n] +=
                                        (*digit_prod)[n] << shift;
                            }
                        }
                    }
                }
#pragma omp critical
                {
                    for (int c = 0; c < component_count; c++)
                        for (uint32_t n = 0; n < P::n; n++)
                            (*inner)[c][n] += (*local_inner)[c][n];
                }
            }
            if (trace_linear) {
                ckks_detail::ckks_trace_linear_group_event(
                    "dual_fused_digit_products", group_index,
                    ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
                trace_section_start = std::chrono::steady_clock::now();
            }

            for (int c = 0; c <= static_cast<int>(P::k); c++)
                for (uint32_t n = 0; n < P::n; n++)
                    (*inner)[c][n] =
                        ckks_detail::roundedLevelRightShift<
                            P, LogQ, PlainLogDelta>((*inner)[c][n]);
            ckks_detail::reduceTRLWEToLevel<P, out_log_q>(*inner);

            const TRLWE<P> *to_add = inner.get();
            if (lhs_group.giant_step != 0) {
                CKKSRotateSlots<P, out_log_q>(
                    *giant_rotated, *inner,
                    lhs_plan.k_step * lhs_group.giant_step, output_gk);
                to_add = giant_rotated.get();
            }

            if (!res_initialized) {
                res.ct = *to_add;
                res_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, *to_add);
            }
            if (trace_linear) {
                ckks_detail::ckks_trace_linear_group_event(
                    "dual_group_finalize", group_index,
                    ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
                trace_section_start = std::chrono::steady_clock::now();
            }
        }
        ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
        if (trace_linear)
            ckks_detail::ckks_trace_linear_event(
                "dual_linear_transform_total",
                ckks_detail::ckks_trace_elapsed_ms(trace_total_start));
        return;
    }

    bool res_initialized = false;
    auto term = std::make_unique<TRLWE<P>>();
    auto inner = std::make_unique<TRLWE<P>>();
    auto giant_rotated = std::make_unique<TRLWE<P>>();

    for (std::size_t group_index = 0; group_index < lhs_plan.groups.size();
         group_index++) {
        const auto &lhs_group = lhs_plan.groups[group_index];
        const auto &rhs_group = rhs_plan.groups[group_index];
        assert(lhs_group.giant_step == rhs_group.giant_step);

        bool inner_initialized = false;
        for (const auto &entry : lhs_group.terms) {
            assert((*lhs_baby)[static_cast<std::size_t>(entry.baby_step)]);
            CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(
                *term,
                *(*lhs_baby)[static_cast<std::size_t>(entry.baby_step)],
                entry.plain.poly);

            if (!inner_initialized) {
                *inner = *term;
                inner_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(*inner, *term);
            }
        }
        for (const auto &entry : rhs_group.terms) {
            assert((*rhs_baby)[static_cast<std::size_t>(entry.baby_step)]);
            CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(
                *term,
                *(*rhs_baby)[static_cast<std::size_t>(entry.baby_step)],
                entry.plain.poly);

            if (!inner_initialized) {
                *inner = *term;
                inner_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(*inner, *term);
            }
        }

        const TRLWE<P> *to_add = inner.get();
        if (lhs_group.giant_step != 0) {
            CKKSRotateSlots<P, out_log_q>(
                *giant_rotated, *inner, lhs_plan.k_step * lhs_group.giant_step,
                output_gk);
            to_add = giant_rotated.get();
        }

        if (!res_initialized) {
            res.ct = *to_add;
            res_initialized = true;
        }
        else {
            CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, *to_add);
        }
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSLinearTransformBSGSDualInputDirect(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &lhs,
    const CKKSCiphertext<P, LogQ, LogDelta> &rhs,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &lhs_plan,
    const CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &rhs_plan,
    const InputGaloisKey &input_gk, const OutputGaloisKey &output_gk)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    assert(lhs_plan.k_step > 0 &&
           lhs_plan.k_step <= static_cast<int>(P::n / 2));
    assert(rhs_plan.k_step == lhs_plan.k_step);
    assert(!lhs_plan.groups.empty());
    assert(lhs_plan.groups.size() == rhs_plan.groups.size());

    std::vector<bool> lhs_baby_used;
    std::vector<bool> rhs_baby_used;
    CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta, PlainLogDelta>(
        lhs_baby_used, lhs_plan);
    CKKSLinearTransformPlanUsedBabySteps<P, LogQ, LogDelta, PlainLogDelta>(
        rhs_baby_used, rhs_plan);
    for (std::size_t i = 0; i < lhs_baby_used.size(); i++)
        lhs_baby_used[i] = lhs_baby_used[i] || rhs_baby_used[i];

    auto lhs_baby =
        std::make_unique<ckks_detail::CKKSSparseBabyStepTable<P>>(
            lhs_plan.k_step);
    auto rhs_baby =
        std::make_unique<ckks_detail::CKKSSparseBabyStepTable<P>>(
            lhs_plan.k_step);
    ckks_detail::CKKSBuildSparseDirectBabyStepRotationTable<P, LogQ>(
        *lhs_baby, lhs.ct, lhs_baby_used, input_gk);
    ckks_detail::CKKSBuildSparseDirectBabyStepRotationTable<P, LogQ>(
        *rhs_baby, rhs.ct, lhs_baby_used, input_gk);
    ckks_detail::maybe_clear_cached_keys(input_gk);

    bool res_initialized = false;
    auto term = std::make_unique<TRLWE<P>>();
    auto inner = std::make_unique<TRLWE<P>>();
    auto giant_rotated = std::make_unique<TRLWE<P>>();

    for (std::size_t group_index = 0; group_index < lhs_plan.groups.size();
         group_index++) {
        const auto &lhs_group = lhs_plan.groups[group_index];
        const auto &rhs_group = rhs_plan.groups[group_index];
        assert(lhs_group.giant_step == rhs_group.giant_step);

        bool inner_initialized = false;
        for (const auto &entry : lhs_group.terms) {
            assert((*lhs_baby)[static_cast<std::size_t>(entry.baby_step)]);
            CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(
                *term,
                *(*lhs_baby)[static_cast<std::size_t>(entry.baby_step)],
                entry.plain.poly);

            if (!inner_initialized) {
                *inner = *term;
                inner_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(*inner, *term);
            }
        }
        for (const auto &entry : rhs_group.terms) {
            assert((*rhs_baby)[static_cast<std::size_t>(entry.baby_step)]);
            CKKSPlainMulRescaleTRLWE<P, LogQ, PlainLogDelta>(
                *term,
                *(*rhs_baby)[static_cast<std::size_t>(entry.baby_step)],
                entry.plain.poly);

            if (!inner_initialized) {
                *inner = *term;
                inner_initialized = true;
            }
            else {
                CKKSAddTRLWEInPlace<P, out_log_q>(*inner, *term);
            }
        }

        const TRLWE<P> *to_add = inner.get();
        if (lhs_group.giant_step != 0) {
            CKKSRotateSlots<P, out_log_q>(
                *giant_rotated, *inner, lhs_plan.k_step * lhs_group.giant_step,
                output_gk);
            to_add = giant_rotated.get();
        }

        if (!res_initialized) {
            res.ct = *to_add;
            res_initialized = true;
        }
        else {
            CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, *to_add);
        }
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSLinearTransformBSGS(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const std::vector<CKKSSlotVector<P>> &diagonals,
    const std::vector<int> &rotation_offsets, int k_step,
    const InputGaloisKey &input_gk, const OutputGaloisKey &output_gk)
{
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> plan;
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, diagonals, rotation_offsets, k_step);
    CKKSLinearTransformBSGS<P, LogQ, LogDelta, PlainLogDelta>(
        res, ct, plan, input_gk, output_gk);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSLinearTransformBSGSDirect(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const std::vector<CKKSSlotVector<P>> &diagonals,
    const std::vector<int> &rotation_offsets, int k_step,
    const InputGaloisKey &input_gk, const OutputGaloisKey &output_gk)
{
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> plan;
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, diagonals, rotation_offsets, k_step);
    CKKSLinearTransformBSGSDirect<P, LogQ, LogDelta, PlainLogDelta>(
        res, ct, plan, input_gk, output_gk);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
struct CKKSStagedPlainMulTraits {
    static_assert(LogQ > StageCount * PlainLogDelta);
    static constexpr std::uint32_t log_q = LogQ - StageCount * PlainLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount>
using CKKSStagedPlainMulResult =
    typename CKKSStagedPlainMulTraits<P, LogQ, LogDelta, PlainLogDelta,
                                      StageCount>::Ciphertext;

namespace ckks_detail {

template <std::size_t I, class KeyProvider>
inline void maybe_release_key(const KeyProvider &keys)
{
    if constexpr (requires { keys.template release<I>(); }) {
        keys.template release<I>();
    }
}

template <std::size_t KeyOffset, std::size_t I, class P,
          std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::size_t StageCount, class GaloisKeyChain>
inline void CKKSLinearTransformStagesBSGSImpl(
    CKKSStagedPlainMulResult<P, StartLogQ, LogDelta, PlainLogDelta,
                             StageCount> &res,
    const CKKSCiphertext<P, StartLogQ - I * PlainLogDelta, LogDelta> &ct,
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    int k_step, const GaloisKeyChain &gk_chain)
{
    if constexpr (I == StageCount) {
        res = ct;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * PlainLogDelta;
        CKKSLinearTransformPlan<P, log_q, LogDelta, PlainLogDelta> plan;
        CKKSBuildLinearTransformBSGSPlan<P, log_q, LogDelta, PlainLogDelta>(
            plan, stages[first_stage + I], k_step);

        auto next =
            std::make_unique<CKKSPlainMulResult<P, log_q, LogDelta,
                                                PlainLogDelta>>();
        CKKSLinearTransformBSGS<P, log_q, LogDelta, PlainLogDelta>(
            *next, ct, plan, gk_chain.template get<KeyOffset + I>(),
            gk_chain.template get<KeyOffset + I + 1>());
        maybe_release_key<KeyOffset + I>(gk_chain);
        if constexpr (I + 1 == StageCount) {
            maybe_release_key<KeyOffset + I + 1>(gk_chain);
        }
        CKKSLinearTransformStagesBSGSImpl<KeyOffset, I + 1, P, StartLogQ,
                                          LogDelta, PlainLogDelta, StageCount>(
            res, *next, stages, first_stage, k_step, gk_chain);
    }
}

template <std::size_t KeyOffset, std::size_t I, class P,
          std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::size_t StageCount, class GaloisKeyChain>
inline void CKKSLinearTransformStagesBSGSDirectImpl(
    CKKSStagedPlainMulResult<P, StartLogQ, LogDelta, PlainLogDelta,
                             StageCount> &res,
    const CKKSCiphertext<P, StartLogQ - I * PlainLogDelta, LogDelta> &ct,
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    int k_step, const GaloisKeyChain &gk_chain)
{
    if constexpr (I == StageCount) {
        res = ct;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * PlainLogDelta;
        CKKSLinearTransformPlan<P, log_q, LogDelta, PlainLogDelta> plan;
        CKKSBuildLinearTransformBSGSPlan<P, log_q, LogDelta, PlainLogDelta>(
            plan, stages[first_stage + I], k_step);

        auto next =
            std::make_unique<CKKSPlainMulResult<P, log_q, LogDelta,
                                                PlainLogDelta>>();
        CKKSLinearTransformBSGSDirect<P, log_q, LogDelta, PlainLogDelta>(
            *next, ct, plan, gk_chain.template get<KeyOffset + I>(),
            gk_chain.template get<KeyOffset + I + 1>());
        maybe_release_key<KeyOffset + I>(gk_chain);
        if constexpr (I + 1 == StageCount) {
            maybe_release_key<KeyOffset + I + 1>(gk_chain);
        }
        CKKSLinearTransformStagesBSGSDirectImpl<
            KeyOffset, I + 1, P, StartLogQ, LogDelta, PlainLogDelta,
            StageCount>(res, *next, stages, first_stage, k_step, gk_chain);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount,
          class GaloisKeyChain>
inline void CKKSLinearTransformStagesBSGS(
    CKKSStagedPlainMulResult<P, LogQ, LogDelta, PlainLogDelta, StageCount>
        &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    int k_step, const GaloisKeyChain &gk_chain)
{
    assert(first_stage + StageCount <= stages.size());
    ckks_detail::CKKSLinearTransformStagesBSGSImpl<0, 0, P, LogQ, LogDelta,
                                                   PlainLogDelta, StageCount>(
        res, ct, stages, first_stage, k_step, gk_chain);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount,
          class GaloisKeyChain>
inline void CKKSLinearTransformStagesBSGSDirect(
    CKKSStagedPlainMulResult<P, LogQ, LogDelta, PlainLogDelta, StageCount>
        &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSLinearTransformStages<P> &stages, std::size_t first_stage,
    int k_step, const GaloisKeyChain &gk_chain)
{
    assert(first_stage + StageCount <= stages.size());
    ckks_detail::CKKSLinearTransformStagesBSGSDirectImpl<
        0, 0, P, LogQ, LogDelta, PlainLogDelta, StageCount>(
        res, ct, stages, first_stage, k_step, gk_chain);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount,
          class GaloisKeyChain>
inline void CKKSLinearTransformStagesBSGSDualInputSharedTail(
    CKKSStagedPlainMulResult<P, LogQ, LogDelta, PlainLogDelta, StageCount>
        &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &lhs,
    const CKKSCiphertext<P, LogQ, LogDelta> &rhs,
    const CKKSLinearTransformStages<P> &lhs_stages,
    const CKKSLinearTransformStages<P> &rhs_stages, std::size_t first_stage,
    int k_step, const GaloisKeyChain &gk_chain)
{
    static_assert(StageCount > 0);
    assert(first_stage + StageCount <= lhs_stages.size());
    assert(first_stage + StageCount <= rhs_stages.size());
#ifndef NDEBUG
    for (std::size_t i = first_stage + 1; i < first_stage + StageCount; i++) {
        assert(lhs_stages[i].rotation_offsets == rhs_stages[i].rotation_offsets);
        assert(lhs_stages[i].diagonals == rhs_stages[i].diagonals);
    }
#endif

    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> lhs_plan;
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> rhs_plan;
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        lhs_plan, lhs_stages[first_stage], k_step);
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        rhs_plan, rhs_stages[first_stage], k_step);

    auto next =
        std::make_unique<CKKSPlainMulResult<P, LogQ, LogDelta,
                                            PlainLogDelta>>();
    CKKSLinearTransformBSGSDualInput<P, LogQ, LogDelta, PlainLogDelta>(
        *next, lhs, rhs, lhs_plan, rhs_plan, gk_chain.template get<0>(),
        gk_chain.template get<1>());
    ckks_detail::maybe_release_key<0>(gk_chain);

    if constexpr (StageCount == 1) {
        res = *next;
        ckks_detail::maybe_release_key<1>(gk_chain);
    }
    else {
        constexpr std::uint32_t tail_log_q = LogQ - PlainLogDelta;
        ckks_detail::CKKSLinearTransformStagesBSGSImpl<
            1, 0, P, tail_log_q, LogDelta, PlainLogDelta, StageCount - 1>(
            res, *next, lhs_stages, first_stage + 1, k_step, gk_chain);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, std::size_t StageCount,
          class GaloisKeyChain>
inline void CKKSLinearTransformStagesBSGSDualInputSharedTailDirect(
    CKKSStagedPlainMulResult<P, LogQ, LogDelta, PlainLogDelta, StageCount>
        &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &lhs,
    const CKKSCiphertext<P, LogQ, LogDelta> &rhs,
    const CKKSLinearTransformStages<P> &lhs_stages,
    const CKKSLinearTransformStages<P> &rhs_stages, std::size_t first_stage,
    int k_step, const GaloisKeyChain &gk_chain)
{
    static_assert(StageCount > 0);
    assert(first_stage + StageCount <= lhs_stages.size());
    assert(first_stage + StageCount <= rhs_stages.size());
#ifndef NDEBUG
    for (std::size_t i = first_stage + 1; i < first_stage + StageCount; i++) {
        assert(lhs_stages[i].rotation_offsets == rhs_stages[i].rotation_offsets);
        assert(lhs_stages[i].diagonals == rhs_stages[i].diagonals);
    }
#endif

    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> lhs_plan;
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> rhs_plan;
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        lhs_plan, lhs_stages[first_stage], k_step);
    CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
        rhs_plan, rhs_stages[first_stage], k_step);

    auto next =
        std::make_unique<CKKSPlainMulResult<P, LogQ, LogDelta,
                                            PlainLogDelta>>();
    CKKSLinearTransformBSGSDualInputDirect<P, LogQ, LogDelta, PlainLogDelta>(
        *next, lhs, rhs, lhs_plan, rhs_plan, gk_chain.template get<0>(),
        gk_chain.template get<1>());
    ckks_detail::maybe_release_key<0>(gk_chain);

    if constexpr (StageCount == 1) {
        res = *next;
        ckks_detail::maybe_release_key<1>(gk_chain);
    }
    else {
        constexpr std::uint32_t tail_log_q = LogQ - PlainLogDelta;
        ckks_detail::CKKSLinearTransformStagesBSGSDirectImpl<
            1, 0, P, tail_log_q, LogDelta, PlainLogDelta, StageCount - 1>(
            res, *next, lhs_stages, first_stage + 1, k_step, gk_chain);
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
struct CKKSRealLinearTransformPlan {
    static_assert(PlainLogDelta < LogQ);
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t plain_log_delta = PlainLogDelta;
    static constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;

    bool has_direct = false;
    bool has_conjugate = false;
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> direct{};
    CKKSLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> conjugate{};
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildRealLinearTransformPlan(
    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const std::vector<CKKSSlotVector<P>> &direct_diagonals,
    const std::vector<int> &direct_offsets,
    const std::vector<CKKSSlotVector<P>> &conjugate_diagonals,
    const std::vector<int> &conjugate_offsets, int k_step)
{
    assert(!direct_diagonals.empty() || !conjugate_diagonals.empty());
    plan.has_direct = !direct_diagonals.empty();
    plan.has_conjugate = !conjugate_diagonals.empty();
    if (plan.has_direct)
        CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
            plan.direct, direct_diagonals, direct_offsets, k_step);
    if (plan.has_conjugate)
        CKKSBuildLinearTransformBSGSPlan<P, LogQ, LogDelta, PlainLogDelta>(
            plan.conjugate, conjugate_diagonals, conjugate_offsets, k_step);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta, class InputGaloisKey,
          class OutputGaloisKey>
inline void CKKSRealLinearTransform(
    CKKSPlainMulResult<P, LogQ, LogDelta, PlainLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    const InputGaloisKey &input_gk, const OutputGaloisKey &output_gk)
{
    static_assert(PlainLogDelta < LogQ);
    constexpr std::uint32_t out_log_q = LogQ - PlainLogDelta;
    assert(plan.has_direct || plan.has_conjugate);

    bool initialized = false;
    auto term = std::make_unique<CKKSPlainMulResult<P, LogQ, LogDelta,
                                                   PlainLogDelta>>();
    if (plan.has_direct) {
        CKKSLinearTransformBSGS<P, LogQ, LogDelta, PlainLogDelta>(
            *term, ct, plan.direct, input_gk, output_gk);
        res = *term;
        initialized = true;
    }

    if (plan.has_conjugate) {
        auto conjugated = std::make_unique<CKKSCiphertext<P, LogQ, LogDelta>>();
        CKKSConjugateSlots<P, LogQ>(conjugated->ct, ct.ct, input_gk);
        CKKSLinearTransformBSGS<P, LogQ, LogDelta, PlainLogDelta>(
            *term, *conjugated, plan.conjugate, input_gk, output_gk);
        if (!initialized) {
            res = *term;
            initialized = true;
        }
        else {
            CKKSAddTRLWEInPlace<P, out_log_q>(res.ct, term->ct);
        }
    }
    ckks_detail::reduceTRLWEToLevel<P, out_log_q>(res.ct);
}

namespace ckks_detail {

template <class P>
inline void addCKKSStageDiagonal(CKKSLinearTransformStage<P> &stage, int offset,
                                 const CKKSSlotVector<P> &diag)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    offset = ((offset % half) + half) % half;
    for (std::size_t i = 0; i < stage.rotation_offsets.size(); i++) {
        if (stage.rotation_offsets[i] != offset) continue;
        for (int j = 0; j < half; j++) stage.diagonals[i][j] += diag[j];
        return;
    }
    stage.rotation_offsets.push_back(offset);
    stage.diagonals.push_back(diag);
}

template <class P>
inline std::vector<std::complex<long double>> ckksPowerOfFiveTable()
{
    constexpr int half = static_cast<int>(P::n) / 2;
    std::vector<std::complex<long double>> roots(4 * half + 1);
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    for (int i = 0; i <= 4 * half; i++) {
        const long double angle =
            2.0L * pi * static_cast<long double>(i) /
            static_cast<long double>(4 * half);
        roots[i] = {std::cos(angle), std::sin(angle)};
    }
    return roots;
}

template <class P>
inline std::vector<int> ckksPowerOfFiveExponents()
{
    constexpr int half = static_cast<int>(P::n) / 2;
    std::vector<int> pow5(2 * half + 1);
    pow5[0] = 1;
    for (int i = 1; i <= 2 * half; i++)
        pow5[i] = (pow5[i - 1] * 5) & (4 * half - 1);
    return pow5;
}

}  // namespace ckks_detail

template <class P>
inline std::uint32_t CKKSBitReverseSlotIndex(std::uint32_t index)
{
    static_assert(P::nbit > 0);
    assert(index < P::n / 2);
    return ckks_detail::bitReverse(index, P::nbit - 1);
}

// Maps CKKS slots to packed coefficients in bit-reversed order:
// out[bitreverse(i)] = coeff[i] + i*coeff[i + N/2].
template <class P>
inline void CKKSBuildCoeffToPackedSlotStages(
    CKKSLinearTransformStages<P> &stages)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    constexpr int log_slots = static_cast<int>(P::nbit) - 1;
    static_assert(log_slots > 0);

    stages.clear();
    stages.reserve(log_slots);
    const auto roots = ckks_detail::ckksPowerOfFiveTable<P>();
    const auto pow5 = ckks_detail::ckksPowerOfFiveExponents<P>();
    const long double stage_scale = 0.5L;

    for (int m = half; m >= 2; m >>= 1) {
        CKKSLinearTransformStage<P> stage;
        CKKSSlotVector<P> a{};
        CKKSSlotVector<P> b{};
        CKKSSlotVector<P> c{};
        a.fill({0.0, 0.0});
        b.fill({0.0, 0.0});
        c.fill({0.0, 0.0});

        const int tt = m >> 1;
        const int gap = half / m;
        const int mask = (m << 2) - 1;
        for (int i = 0; i < half; i += m) {
            for (int j = 0; j < tt; j++) {
                const int k = ((m << 2) - (pow5[j] & mask)) * gap;
                const int idx1 = i + j;
                const int idx2 = idx1 + tt;
                a[idx1] = {static_cast<double>(stage_scale), 0.0};
                a[idx2] = static_cast<std::complex<double>>(
                    -stage_scale * roots[k]);
                b[idx1] = {static_cast<double>(stage_scale), 0.0};
                c[idx2] = static_cast<std::complex<double>>(
                    stage_scale * roots[k]);
            }
        }

        ckks_detail::addCKKSStageDiagonal<P>(stage, 0, a);
        ckks_detail::addCKKSStageDiagonal<P>(stage, tt, b);
        ckks_detail::addCKKSStageDiagonal<P>(stage, half - tt, c);
        stages.push_back(std::move(stage));
    }
}

// Inverse of CKKSBuildCoeffToPackedSlotStages.  The input is the same
// bit-reversed packed coefficient layout emitted by CoeffToPackedSlot.
template <class P>
inline void CKKSBuildPackedSlotToCoeffStages(
    CKKSLinearTransformStages<P> &stages)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    constexpr int log_slots = static_cast<int>(P::nbit) - 1;
    static_assert(log_slots > 0);

    stages.clear();
    stages.reserve(log_slots);
    const auto roots = ckks_detail::ckksPowerOfFiveTable<P>();
    const auto pow5 = ckks_detail::ckksPowerOfFiveExponents<P>();

    for (int m = 2; m <= half; m <<= 1) {
        CKKSLinearTransformStage<P> stage;
        CKKSSlotVector<P> a{};
        CKKSSlotVector<P> b{};
        CKKSSlotVector<P> c{};
        a.fill({0.0, 0.0});
        b.fill({0.0, 0.0});
        c.fill({0.0, 0.0});

        const int tt = m >> 1;
        const int gap = half / m;
        const int mask = (m << 2) - 1;
        for (int i = 0; i < half; i += m) {
            for (int j = 0; j < tt; j++) {
                const int k = (pow5[j] & mask) * gap;
                const int idx1 = i + j;
                const int idx2 = idx1 + tt;
                a[idx1] = {1.0, 0.0};
                a[idx2] =
                    static_cast<std::complex<double>>(-roots[k]);
                b[idx1] = static_cast<std::complex<double>>(roots[k]);
                c[idx2] = {1.0, 0.0};
            }
        }

        ckks_detail::addCKKSStageDiagonal<P>(stage, 0, a);
        ckks_detail::addCKKSStageDiagonal<P>(stage, tt, b);
        ckks_detail::addCKKSStageDiagonal<P>(stage, half - tt, c);
        stages.push_back(std::move(stage));
    }
}

template <class P>
inline void CKKSBuildCoeffToPackedSlotDiagonals(
    std::vector<CKKSSlotVector<P>> &direct_diagonals,
    std::vector<int> &direct_offsets,
    std::vector<CKKSSlotVector<P>> &conjugate_diagonals,
    std::vector<int> &conjugate_offsets)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const auto &slot_to_eval = ckks_detail::ckksSlotToEvalIndex<P>();

    direct_diagonals.assign(half, CKKSSlotVector<P>{});
    conjugate_diagonals.assign(half, CKKSSlotVector<P>{});
    direct_offsets.resize(half);
    conjugate_offsets.resize(half);
    for (int d = 0; d < half; d++) {
        direct_offsets[d] = d;
        conjugate_offsets[d] = d;
        direct_diagonals[d].fill({0.0, 0.0});
        conjugate_diagonals[d].fill({0.0, 0.0});
    }

    for (int r = 0; r < half; r++) {
        for (int d = 0; d < half; d++) {
            const int c = (r + d) % half;
            const long double h =
                static_cast<long double>(2 * slot_to_eval[c] + 1);
            const long double theta = pi * h / static_cast<long double>(P::n);
            const auto e0 = std::complex<long double>(
                std::cos(theta * r), -std::sin(theta * r));
            const auto e1 = std::complex<long double>(
                std::cos(theta * (r + half)),
                -std::sin(theta * (r + half)));
            const auto f0 = std::conj(e0);
            const auto f1 = std::conj(e1);
            direct_diagonals[d][r] = static_cast<std::complex<double>>(
                (e0 + std::complex<long double>(0.0L, 1.0L) * e1) /
                static_cast<long double>(P::n));
            conjugate_diagonals[d][r] = static_cast<std::complex<double>>(
                (f0 + std::complex<long double>(0.0L, 1.0L) * f1) /
                static_cast<long double>(P::n));
        }
    }
}

template <class P>
inline void CKKSBuildPackedSlotToCoeffDiagonals(
    std::vector<CKKSSlotVector<P>> &direct_diagonals,
    std::vector<int> &direct_offsets,
    std::vector<CKKSSlotVector<P>> &conjugate_diagonals,
    std::vector<int> &conjugate_offsets)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    const auto &slot_to_eval = ckks_detail::ckksSlotToEvalIndex<P>();

    direct_diagonals.assign(half, CKKSSlotVector<P>{});
    conjugate_diagonals.assign(half, CKKSSlotVector<P>{});
    direct_offsets.resize(half);
    conjugate_offsets.resize(half);
    for (int d = 0; d < half; d++) {
        direct_offsets[d] = d;
        conjugate_offsets[d] = d;
        direct_diagonals[d].fill({0.0, 0.0});
        conjugate_diagonals[d].fill({0.0, 0.0});
    }

    for (int i = 0; i < half; i++) {
        const long double h =
            static_cast<long double>(2 * slot_to_eval[i] + 1);
        const long double theta = pi * h / static_cast<long double>(P::n);
        for (int d = 0; d < half; d++) {
            const int c = (i + d) % half;
            const auto e0 = std::complex<long double>(
                std::cos(theta * c), std::sin(theta * c));
            const auto e1 = std::complex<long double>(
                std::cos(theta * (c + half)),
                std::sin(theta * (c + half)));
            direct_diagonals[d][i] = static_cast<std::complex<double>>(
                0.5L *
                (e0 - std::complex<long double>(0.0L, 1.0L) * e1));
            conjugate_diagonals[d][i] = static_cast<std::complex<double>>(
                0.5L *
                (e0 + std::complex<long double>(0.0L, 1.0L) * e1));
        }
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildCoeffToPackedSlotPlan(
    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    int k_step)
{
    std::vector<CKKSSlotVector<P>> direct_diagonals;
    std::vector<CKKSSlotVector<P>> conjugate_diagonals;
    std::vector<int> direct_offsets;
    std::vector<int> conjugate_offsets;
    CKKSBuildCoeffToPackedSlotDiagonals<P>(
        direct_diagonals, direct_offsets, conjugate_diagonals,
        conjugate_offsets);
    CKKSBuildRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, direct_diagonals, direct_offsets, conjugate_diagonals,
        conjugate_offsets, k_step);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PlainLogDelta>
inline void CKKSBuildPackedSlotToCoeffPlan(
    CKKSRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta> &plan,
    int k_step)
{
    std::vector<CKKSSlotVector<P>> direct_diagonals;
    std::vector<CKKSSlotVector<P>> conjugate_diagonals;
    std::vector<int> direct_offsets;
    std::vector<int> conjugate_offsets;
    CKKSBuildPackedSlotToCoeffDiagonals<P>(
        direct_diagonals, direct_offsets, conjugate_diagonals,
        conjugate_offsets);
    CKKSBuildRealLinearTransformPlan<P, LogQ, LogDelta, PlainLogDelta>(
        plan, direct_diagonals, direct_offsets, conjugate_diagonals,
        conjugate_offsets, k_step);
}

template <
    class P, std::uint32_t LogDelta = 50,
    std::uint32_t LogMessageRatio = 8, std::uint32_t BootLogQ = 880,
    std::uint32_t LinearPlainLogDelta = 50,
    std::uint32_t LinearFuseRadix = 5, std::uint32_t EvalModDegree = 30,
    std::uint32_t EvalModK = 16, std::uint32_t EvalModDoubleAngle = 3,
    std::uint32_t EvalModInvDegree = 0,
    std::uint32_t EvalModLogScale = LogDelta, int LinearBSGSStep = 128,
    std::uint32_t ModRaiseMaskBound = 0,
    std::uint32_t CoeffToSlotPlainLogDelta = LinearPlainLogDelta,
    std::uint32_t ComponentSplitPlainLogDelta = CoeffToSlotPlainLogDelta,
    std::uint32_t SlotToCoeffPlainLogDelta =
        (LinearPlainLogDelta == 50 ? 35 : LinearPlainLogDelta),
    std::uint32_t CoeffToSlotFuseRadix = LinearFuseRadix,
    std::uint32_t SlotToCoeffFuseRadix = LinearFuseRadix,
    int HybridGiantDirectPopcountThreshold = 4>
struct CKKSDenseBootstrapSchedule {
    static_assert(ckks_detail::supported_torus_v<P>);
    static_assert(LogDelta > 0);
    static_assert(LogMessageRatio > 0);
    static_assert(LinearPlainLogDelta > 0);
    static_assert(LinearFuseRadix > 0);
    static_assert(CoeffToSlotPlainLogDelta > 0);
    static_assert(ComponentSplitPlainLogDelta > 0);
    static_assert(SlotToCoeffPlainLogDelta > 0);
    static_assert(CoeffToSlotFuseRadix > 0);
    static_assert(SlotToCoeffFuseRadix > 0);
    static_assert(LinearBSGSStep > 0);
    static_assert(HybridGiantDirectPopcountThreshold > 0);
    static_assert(EvalModK > 0);
    static_assert(ModRaiseMaskBound < EvalModK);
    static_assert(EvalModInvDegree == 0 || EvalModInvDegree > 1);
    static_assert(EvalModDegree >= 2 * (EvalModK - 1),
                  "cos-discrete EvalMod requires degree at least 2*(K-1)");

    using Param = P;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t log_message_ratio = LogMessageRatio;
    static constexpr std::uint32_t input_log_q = LogDelta + LogMessageRatio;
    static constexpr std::uint32_t boot_log_q = BootLogQ;
    static constexpr std::uint32_t linear_plain_log_delta = LinearPlainLogDelta;
    static constexpr std::uint32_t linear_fuse_radix = LinearFuseRadix;
    static constexpr std::uint32_t coeff_to_slot_plain_log_delta =
        CoeffToSlotPlainLogDelta;
    static constexpr std::uint32_t component_split_plain_log_delta =
        ComponentSplitPlainLogDelta;
    static constexpr std::uint32_t slot_to_coeff_plain_log_delta =
        SlotToCoeffPlainLogDelta;
    static constexpr std::uint32_t coeff_to_slot_fuse_radix =
        CoeffToSlotFuseRadix;
    static constexpr std::uint32_t slot_to_coeff_fuse_radix =
        SlotToCoeffFuseRadix;
    static constexpr int linear_bsgs_step = LinearBSGSStep;
    static constexpr int hybrid_giant_direct_popcount_threshold =
        HybridGiantDirectPopcountThreshold;
    template <std::size_t I>
    static consteval int coeff_to_slot_bsgs_step()
    {
        return LinearBSGSStep;
    }
    template <std::size_t I>
    static consteval int slot_to_coeff_bsgs_step()
    {
        return LinearBSGSStep;
    }
    static constexpr std::uint32_t raw_linear_stage_count = P::nbit - 1;
    static constexpr std::uint32_t coeff_to_slot_level_count =
        ckks_detail::ceil_div(raw_linear_stage_count, CoeffToSlotFuseRadix);
    static constexpr std::uint32_t slot_to_coeff_level_count =
        ckks_detail::ceil_div(raw_linear_stage_count, SlotToCoeffFuseRadix);
    static constexpr std::uint32_t evalmod_degree = EvalModDegree;
    static constexpr std::uint32_t evalmod_k = EvalModK;
    static constexpr std::uint32_t evalmod_double_angle = EvalModDoubleAngle;
    static constexpr std::uint32_t evalmod_inv_degree = EvalModInvDegree;
    static constexpr std::uint32_t evalmod_log_scale = EvalModLogScale;
    static constexpr std::uint32_t modraise_mask_bound = ModRaiseMaskBound;
    static constexpr std::uint32_t evalmod_polynomial_depth =
        ckks_detail::bit_width_u64(
            ckks_detail::static_max_v<EvalModDegree, 2 * (EvalModK - 1)>) +
        1;
    static constexpr std::uint32_t evalmod_polynomial_power_depth =
        evalmod_polynomial_depth - 1;
    static constexpr std::uint32_t evalmod_inv_depth =
        EvalModInvDegree == 0
            ? 0
            : ckks_detail::bit_width_u64(EvalModInvDegree - 1) + 1;
    static constexpr std::uint32_t evalmod_inv_power_depth =
        evalmod_inv_depth == 0 ? 0 : evalmod_inv_depth - 1;
    static constexpr std::uint32_t evalmod_depth =
        evalmod_polynomial_depth + EvalModDoubleAngle + evalmod_inv_depth;
    static constexpr std::uint32_t evalmod_coeff_rescale_count =
        1 + (EvalModInvDegree == 0 ? 0 : 1);
    static constexpr std::uint32_t evalmod_log_q_consumption =
        (evalmod_polynomial_power_depth + EvalModDoubleAngle +
         evalmod_inv_power_depth) *
            LogDelta +
        evalmod_coeff_rescale_count * EvalModLogScale;
    static constexpr double message_ratio =
        ckks_detail::exp2_double(LogMessageRatio);
    // Modulus raising already introduces key-dependent multiples of the input
    // level.  Any extra explicit phase mask must fit in the remaining EvalMod
    // interval, so it is opt-in rather than tied to EvalModK by default.
    static constexpr std::uint32_t evalmod_mask_bound =
        EvalModK == 0 ? 0 : EvalModK - 1;
    static constexpr std::uint32_t bounded_sparse_secret_key_weight =
        EvalModK > ModRaiseMaskBound + 2 ? EvalModK - ModRaiseMaskBound - 2
                                         : 0;
    static constexpr double coeff_to_slot_scaling_factor =
        1.0 / (static_cast<double>(EvalModK) * message_ratio);
    static constexpr double slot_to_coeff_scaling_factor = message_ratio;

    static constexpr std::uint32_t after_coeff_to_slot_log_q =
        BootLogQ - coeff_to_slot_level_count * CoeffToSlotPlainLogDelta;
    static constexpr std::uint32_t after_component_split_log_q =
        after_coeff_to_slot_log_q - ComponentSplitPlainLogDelta;
    static constexpr std::uint32_t after_evalmod_log_q =
        after_component_split_log_q - evalmod_log_q_consumption;
    static constexpr std::uint32_t output_log_q =
        after_evalmod_log_q - slot_to_coeff_level_count * SlotToCoeffPlainLogDelta;
    static constexpr std::uint32_t post_bootstrap_product_log_q =
        output_log_q - LogDelta;
    static constexpr bool supports_post_bootstrap_product =
        post_bootstrap_product_log_q >= input_log_q;
    static constexpr std::uint32_t post_bootstrap_product_slack =
        supports_post_bootstrap_product
            ? post_bootstrap_product_log_q - input_log_q
            : 0;

    static_assert(input_log_q < BootLogQ);
    static_assert(BootLogQ <= ckks_detail::torus_width_v<P>);
    static_assert(after_coeff_to_slot_log_q > ComponentSplitPlainLogDelta);
    static_assert(after_component_split_log_q > evalmod_log_q_consumption);
    static_assert(output_log_q > LogDelta);

    using InputCiphertext = CKKSCiphertext<P, input_log_q, log_delta>;
    using BootstrapCiphertext = CKKSCiphertext<P, boot_log_q, log_delta>;
    using CoeffToSlotCiphertext =
        CKKSStagedPlainMulResult<P, boot_log_q, log_delta,
                                 coeff_to_slot_plain_log_delta,
                                 coeff_to_slot_level_count>;
    using ComponentCiphertext =
        CKKSPlainMulResult<P, after_coeff_to_slot_log_q, log_delta,
                           component_split_plain_log_delta>;
    using EvalModCiphertext = CKKSCiphertext<P, after_evalmod_log_q, log_delta>;
    using OutputCiphertext =
        CKKSStagedPlainMulResult<P, after_evalmod_log_q, log_delta,
                                 slot_to_coeff_plain_log_delta,
                                 slot_to_coeff_level_count>;
};

template <int HybridGiantDirectPopcountThreshold = 3>
using lvl6CKKSDenseBootstrapRobustHybridSchedule =
    CKKSDenseBootstrapSchedule<lvl6param, 50, 8, 888, 50, 5, 34, 18, 3, 0, 50,
                               128, 0, 50, 50, 25, 5, 5,
                               HybridGiantDirectPopcountThreshold>;

using lvl6CKKSDenseBootstrapFastSchedule =
    lvl6CKKSDenseBootstrapRobustHybridSchedule<3>;
using lvl6CKKSDenseBootstrapCompactSchedule =
    lvl6CKKSDenseBootstrapRobustHybridSchedule<4>;
using lvl6CKKSDenseBootstrapSchedule =
    lvl6CKKSDenseBootstrapFastSchedule;
using lvl6CKKSDenseBootstrapInverseSchedule =
    CKKSDenseBootstrapSchedule<lvl6param, 40, 8, 880, 40, 5, 52, 18, 4, 3, 40,
                               128, 0, 50, 50, 20, 5, 5, 3>;
struct lvl6CKKSDenseBootstrapTunedSchedule
    : CKKSDenseBootstrapSchedule<lvl6param, 52, 8, 1152, 52, 7, 34, 18, 4, 5,
                                 52, 128, 0, 52, 52, 30, 7, 7, 5> {
    template <std::size_t I>
    static consteval int coeff_to_slot_bsgs_step()
    {
        return I == 0 ? 2048 : 32;
    }

    template <std::size_t I>
    static consteval int slot_to_coeff_bsgs_step()
    {
        return I + 1 == slot_to_coeff_level_count ? 1024 : 64;
    }
};

struct CKKSBoundedCosEvalModPolynomial {
    std::uint32_t k = 0;
    std::uint32_t degree = 0;
    std::uint32_t log_message_ratio = 0;
    std::uint32_t double_angle = 0;
    double message_ratio = 1.0;
    double q_diff = 1.0;
    double sqrt_coeff = 1.0;
    double domain_offset = 0.0;
    std::vector<double> chebyshev_coeffs{};
    std::vector<double> power_coeffs{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(k, degree, log_message_ratio, double_angle, message_ratio,
                q_diff, sqrt_coeff, domain_offset, chebyshev_coeffs,
                power_coeffs);
    }
};

inline double CKKSEvaluateChebyshevUnit(const std::vector<double> &coeffs,
                                        double x)
{
    if (coeffs.empty()) return 0.0;
    double b_next = 0.0;
    double b_next_next = 0.0;
    for (std::size_t i = coeffs.size() - 1; i > 0; i--) {
        const double b = 2.0 * x * b_next - b_next_next + coeffs[i];
        b_next_next = b_next;
        b_next = b;
    }
    return coeffs[0] + x * b_next - b_next_next;
}

inline double CKKSEvaluatePowerPolynomial(const std::vector<double> &coeffs,
                                          double x)
{
    double y = 0.0;
    for (std::size_t i = coeffs.size(); i-- > 0;)
        y = y * x + coeffs[i];
    return y;
}

inline std::vector<double> CKKSChebyshevToPowerCoefficients(
    const std::vector<double> &chebyshev_coeffs)
{
    const std::size_t degree = chebyshev_coeffs.empty()
                                   ? 0
                                   : chebyshev_coeffs.size() - 1;
    std::vector<double> power_coeffs(degree + 1, 0.0);
    if (chebyshev_coeffs.empty()) return {};

    std::vector<double> prev{1.0};
    power_coeffs[0] += chebyshev_coeffs[0];
    if (chebyshev_coeffs.size() == 1) return power_coeffs;

    std::vector<double> cur{0.0, 1.0};
    power_coeffs[1] += chebyshev_coeffs[1];
    for (std::size_t i = 2; i < chebyshev_coeffs.size(); i++) {
        std::vector<double> next(i + 1, 0.0);
        for (std::size_t j = 0; j < cur.size(); j++) next[j + 1] += 2.0 * cur[j];
        for (std::size_t j = 0; j < prev.size(); j++) next[j] -= prev[j];
        for (std::size_t j = 0; j < next.size(); j++)
            power_coeffs[j] += chebyshev_coeffs[i] * next[j];
        prev = std::move(cur);
        cur = std::move(next);
    }
    return power_coeffs;
}

inline CKKSBoundedCosEvalModPolynomial CKKSBuildBoundedCosEvalModPolynomial(
    std::uint32_t k, std::uint32_t degree, std::uint32_t log_message_ratio,
    std::uint32_t double_angle, double q_diff = 1.0,
    bool use_inverse_correction = false)
{
    assert(k > 0);
    constexpr long double pi =
        3.141592653589793238462643383279502884L;
    CKKSBoundedCosEvalModPolynomial poly;
    poly.k = k;
    poly.degree = degree;
    poly.log_message_ratio = log_message_ratio;
    poly.double_angle = double_angle;
    poly.message_ratio =
        ckks_detail::exp2_double(log_message_ratio);
    poly.q_diff = q_diff;

    const double shrink =
        ckks_detail::exp2_double(double_angle);
    poly.sqrt_coeff =
        use_inverse_correction
            ? 1.0
            : std::pow(q_diff / (2.0 * static_cast<double>(pi)), 1.0 / shrink);
    poly.domain_offset =
        1.0 / (4.0 * static_cast<double>(k) * shrink);

    const std::size_t nodes = static_cast<std::size_t>(degree) + 1;
    std::vector<long double> values(nodes);
    for (std::size_t j = 0; j < nodes; j++) {
        const long double theta =
            pi * (static_cast<long double>(j) + 0.5L) /
            static_cast<long double>(nodes);
        const long double x = std::cos(theta);
        const long double unnormalized =
            static_cast<long double>(k) *
            (x + static_cast<long double>(poly.domain_offset));
        const long double angle =
            2.0L * pi * (unnormalized - 0.25L) /
            static_cast<long double>(shrink);
        values[j] = static_cast<long double>(poly.sqrt_coeff) * std::cos(angle);
    }

    poly.chebyshev_coeffs.resize(nodes);
    for (std::size_t i = 0; i < nodes; i++) {
        long double acc = 0.0L;
        for (std::size_t j = 0; j < nodes; j++) {
            const long double theta =
                pi * static_cast<long double>(i) *
                (static_cast<long double>(j) + 0.5L) /
                static_cast<long double>(nodes);
            acc += values[j] * std::cos(theta);
        }
        acc *= 2.0L / static_cast<long double>(nodes);
        if (i == 0) acc *= 0.5L;
        poly.chebyshev_coeffs[i] = static_cast<double>(acc);
    }
    poly.power_coeffs =
        CKKSChebyshevToPowerCoefficients(poly.chebyshev_coeffs);
    return poly;
}

template <class Schedule>
inline CKKSBoundedCosEvalModPolynomial
CKKSBuildBoundedCosEvalModPolynomial()
{
    return CKKSBuildBoundedCosEvalModPolynomial(
        Schedule::evalmod_k, Schedule::evalmod_degree,
        Schedule::log_message_ratio, Schedule::evalmod_double_angle, 1.0,
        Schedule::evalmod_inv_degree != 0);
}

inline std::vector<double> CKKSBuildEvalModInversePowerCoefficients(
    const CKKSBoundedCosEvalModPolynomial &poly, std::uint32_t degree)
{
    if (degree == 0) return {};
    constexpr double pi = 3.141592653589793238462643383279502884;
    std::vector<double> coeffs(static_cast<std::size_t>(degree) + 1, 0.0);
    coeffs[1] = poly.q_diff / (2.0 * pi);
    for (std::uint32_t i = 3; i <= degree; i += 2) {
        const double numerator =
            static_cast<double>(i * i - 4 * i + 4);
        const double denominator = static_cast<double>(i * i - i);
        coeffs[i] = coeffs[i - 2] * numerator / denominator;
    }
    return coeffs;
}

inline double CKKSApplyEvalModInverseCorrection(
    const CKKSBoundedCosEvalModPolynomial &poly, double value,
    std::uint32_t degree)
{
    if (degree == 0) return value;
    return CKKSEvaluatePowerPolynomial(
        CKKSBuildEvalModInversePowerCoefficients(poly, degree), value);
}

inline double CKKSPlainEvalModBoundedCosNormalized(
    const CKKSBoundedCosEvalModPolynomial &poly, double normalized,
    std::uint32_t inverse_degree = 0)
{
    const double x = normalized - poly.domain_offset;
    double y = CKKSEvaluateChebyshevUnit(poly.chebyshev_coeffs, x);
    double sqrt_coeff = poly.sqrt_coeff;
    for (std::uint32_t i = 0; i < poly.double_angle; i++) {
        sqrt_coeff *= sqrt_coeff;
        y = 2.0 * y * y - sqrt_coeff;
    }
    y = CKKSApplyEvalModInverseCorrection(poly, y, inverse_degree);
    return poly.message_ratio * y;
}

inline double CKKSPlainEvalModBoundedCosNormalizedPower(
    const CKKSBoundedCosEvalModPolynomial &poly, double normalized,
    std::uint32_t inverse_degree = 0)
{
    const double x = normalized - poly.domain_offset;
    double y = CKKSEvaluatePowerPolynomial(poly.power_coeffs, x);
    double sqrt_coeff = poly.sqrt_coeff;
    for (std::uint32_t i = 0; i < poly.double_angle; i++) {
        sqrt_coeff *= sqrt_coeff;
        y = 2.0 * y * y - sqrt_coeff;
    }
    y = CKKSApplyEvalModInverseCorrection(poly, y, inverse_degree);
    return poly.message_ratio * y;
}

inline double CKKSPlainEvalModBoundedCos(
    const CKKSBoundedCosEvalModPolynomial &poly, double masked_value,
    std::uint32_t inverse_degree = 0)
{
    const double normalized =
        masked_value /
        (static_cast<double>(poly.k) * poly.message_ratio * poly.q_diff);
    return CKKSPlainEvalModBoundedCosNormalized(poly, normalized,
                                                inverse_degree);
}

inline double CKKSPlainEvalModBoundedCosPower(
    const CKKSBoundedCosEvalModPolynomial &poly, double masked_value,
    std::uint32_t inverse_degree = 0)
{
    const double normalized =
        masked_value /
        (static_cast<double>(poly.k) * poly.message_ratio * poly.q_diff);
    return CKKSPlainEvalModBoundedCosNormalizedPower(poly, normalized,
                                                     inverse_degree);
}

template <class Schedule>
struct CKKSDenseBootstrapLinearPlan {
    using P = typename Schedule::Param;
    CKKSLinearTransformStages<P> coeff_to_slot_stages{};
    CKKSLinearTransformStages<P> slot_to_coeff_stages{};
    CKKSLinearTransformStages<P> slot_to_coeff_imag_stages{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(coeff_to_slot_stages, slot_to_coeff_stages,
                slot_to_coeff_imag_stages);
    }
};

namespace ckks_detail {

template <std::size_t I, class Schedule>
consteval int CKKSDenseBootstrapCoeffToSlotBSGSStep()
{
    if constexpr (requires { Schedule::template coeff_to_slot_bsgs_step<I>(); })
        return Schedule::template coeff_to_slot_bsgs_step<I>();
    else
        return Schedule::linear_bsgs_step;
}

template <std::size_t I, class Schedule>
consteval int CKKSDenseBootstrapSlotToCoeffBSGSStep()
{
    if constexpr (requires { Schedule::template slot_to_coeff_bsgs_step<I>(); })
        return Schedule::template slot_to_coeff_bsgs_step<I>();
    else
        return Schedule::linear_bsgs_step;
}

template <class Schedule, std::size_t I = 0>
constexpr int CKKSDenseBootstrapCoeffToSlotBSGSStepAt(std::size_t index)
{
    if constexpr (I >= Schedule::coeff_to_slot_level_count) {
        return Schedule::linear_bsgs_step;
    }
    else {
        return index == I
                   ? CKKSDenseBootstrapCoeffToSlotBSGSStep<I, Schedule>()
                   : CKKSDenseBootstrapCoeffToSlotBSGSStepAt<Schedule, I + 1>(
                         index);
    }
}

template <class Schedule, std::size_t I = 0>
constexpr int CKKSDenseBootstrapSlotToCoeffBSGSStepAt(std::size_t index)
{
    if constexpr (I >= Schedule::slot_to_coeff_level_count) {
        return Schedule::linear_bsgs_step;
    }
    else {
        return index == I
                   ? CKKSDenseBootstrapSlotToCoeffBSGSStep<I, Schedule>()
                   : CKKSDenseBootstrapSlotToCoeffBSGSStepAt<Schedule, I + 1>(
                         index);
    }
}

}  // namespace ckks_detail

template <class Schedule>
struct CKKSDenseBootstrapRotationKeyUsage {
    using P = typename Schedule::Param;
    std::array<CKKSRotationKeyIndexSet<P>,
               Schedule::coeff_to_slot_level_count + 1>
        coeff_to_slot{};
    CKKSRotationKeyIndexSet<P> packed_conjugate{};
    std::array<CKKSRotationKeyIndexSet<P>,
               Schedule::slot_to_coeff_level_count + 1>
        slot_to_coeff{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(coeff_to_slot, packed_conjugate, slot_to_coeff);
    }
};

template <class Schedule>
struct CKKSDenseBootstrapDirectRotationKeyUsage {
    using P = typename Schedule::Param;
    std::array<CKKSDirectRotationKeyIndexSet<P>,
               Schedule::coeff_to_slot_level_count + 1>
        coeff_to_slot{};
    CKKSRotationKeyIndexSet<P> packed_conjugate{};
    std::array<CKKSDirectRotationKeyIndexSet<P>,
               Schedule::slot_to_coeff_level_count + 1>
        slot_to_coeff{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(coeff_to_slot, packed_conjugate, slot_to_coeff);
    }
};

template <class Schedule>
struct CKKSDenseBootstrapHybridGiantRotationKeyUsage {
    using P = typename Schedule::Param;
    std::array<CKKSRotationKeyIndexSet<P>,
               Schedule::coeff_to_slot_level_count + 1>
        coeff_to_slot_binary{};
    std::array<CKKSDirectRotationKeyIndexSet<P>,
               Schedule::coeff_to_slot_level_count + 1>
        coeff_to_slot_direct{};
    CKKSRotationKeyIndexSet<P> packed_conjugate{};
    std::array<CKKSRotationKeyIndexSet<P>,
               Schedule::slot_to_coeff_level_count + 1>
        slot_to_coeff_binary{};
    std::array<CKKSDirectRotationKeyIndexSet<P>,
               Schedule::slot_to_coeff_level_count + 1>
        slot_to_coeff_direct{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(coeff_to_slot_binary, coeff_to_slot_direct,
                packed_conjugate, slot_to_coeff_binary,
                slot_to_coeff_direct);
    }
};

template <class P, std::size_t N>
inline void CKKSClearRotationKeyIndexSets(
    std::array<CKKSRotationKeyIndexSet<P>, N> &levels)
{
    for (auto &level : levels) CKKSClearRotationKeyIndexSet<P>(level);
}

template <class P, std::size_t N>
inline void CKKSClearDirectRotationKeyIndexSets(
    std::array<CKKSDirectRotationKeyIndexSet<P>, N> &levels)
{
    for (auto &level : levels) CKKSClearDirectRotationKeyIndexSet<P>(level);
}

template <class P, std::size_t N>
inline std::size_t CKKSRotationKeyIndexSetsCount(
    const std::array<CKKSRotationKeyIndexSet<P>, N> &levels)
{
    std::size_t count = 0;
    for (const auto &level : levels)
        count += CKKSRotationKeyIndexSetCount<P>(level);
    return count;
}

template <class P, std::size_t N>
inline std::size_t CKKSDirectRotationKeyIndexSetsCount(
    const std::array<CKKSDirectRotationKeyIndexSet<P>, N> &levels)
{
    std::size_t count = 0;
    for (const auto &level : levels)
        count += CKKSDirectRotationKeyIndexSetCount<P>(level);
    return count;
}

template <class Schedule>
inline void CKKSBuildDenseBootstrapRotationKeyUsage(
    CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan)
{
    using P = typename Schedule::Param;
    CKKSClearRotationKeyIndexSets<P>(usage.coeff_to_slot);
    CKKSClearRotationKeyIndexSet<P>(usage.packed_conjugate);
    CKKSClearRotationKeyIndexSets<P>(usage.slot_to_coeff);

    for (std::size_t i = 0; i < linear_plan.coeff_to_slot_stages.size(); i++) {
        CKKSCollectLinearTransformStageRotationKeyIndices<P>(
            usage.coeff_to_slot[i], usage.coeff_to_slot[i + 1],
            linear_plan.coeff_to_slot_stages[i],
            ckks_detail::CKKSDenseBootstrapCoeffToSlotBSGSStepAt<Schedule>(
                i));
    }
    CKKSMarkConjugationKeyIndex<P>(usage.packed_conjugate);
    for (std::size_t i = 0; i < linear_plan.slot_to_coeff_stages.size(); i++) {
        CKKSCollectLinearTransformStageRotationKeyIndices<P>(
            usage.slot_to_coeff[i], usage.slot_to_coeff[i + 1],
            linear_plan.slot_to_coeff_stages[i],
            ckks_detail::CKKSDenseBootstrapSlotToCoeffBSGSStepAt<Schedule>(
                i));
        CKKSCollectLinearTransformStageRotationKeyIndices<P>(
            usage.slot_to_coeff[i], usage.slot_to_coeff[i + 1],
            linear_plan.slot_to_coeff_imag_stages[i],
            ckks_detail::CKKSDenseBootstrapSlotToCoeffBSGSStepAt<Schedule>(
                i));
    }
}

template <class Schedule>
inline void CKKSBuildDenseBootstrapDirectRotationKeyUsage(
    CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage,
    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan)
{
    using P = typename Schedule::Param;
    CKKSClearDirectRotationKeyIndexSets<P>(usage.coeff_to_slot);
    CKKSClearRotationKeyIndexSet<P>(usage.packed_conjugate);
    CKKSClearDirectRotationKeyIndexSets<P>(usage.slot_to_coeff);

    for (std::size_t i = 0; i < linear_plan.coeff_to_slot_stages.size(); i++) {
        CKKSCollectLinearTransformStageDirectRotationKeyIndices<P>(
            usage.coeff_to_slot[i], usage.coeff_to_slot[i + 1],
            linear_plan.coeff_to_slot_stages[i],
            ckks_detail::CKKSDenseBootstrapCoeffToSlotBSGSStepAt<Schedule>(
                i));
    }
    CKKSMarkConjugationKeyIndex<P>(usage.packed_conjugate);
    for (std::size_t i = 0; i < linear_plan.slot_to_coeff_stages.size(); i++) {
        CKKSCollectLinearTransformStageDirectRotationKeyIndices<P>(
            usage.slot_to_coeff[i], usage.slot_to_coeff[i + 1],
            linear_plan.slot_to_coeff_stages[i],
            ckks_detail::CKKSDenseBootstrapSlotToCoeffBSGSStepAt<Schedule>(
                i));
        CKKSCollectLinearTransformStageDirectRotationKeyIndices<P>(
            usage.slot_to_coeff[i], usage.slot_to_coeff[i + 1],
            linear_plan.slot_to_coeff_imag_stages[i],
            ckks_detail::CKKSDenseBootstrapSlotToCoeffBSGSStepAt<Schedule>(
                i));
    }
}

template <class Schedule>
inline void CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage(
    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan)
{
    using P = typename Schedule::Param;
    CKKSClearRotationKeyIndexSets<P>(usage.coeff_to_slot_binary);
    CKKSClearDirectRotationKeyIndexSets<P>(usage.coeff_to_slot_direct);
    CKKSClearRotationKeyIndexSet<P>(usage.packed_conjugate);
    CKKSClearRotationKeyIndexSets<P>(usage.slot_to_coeff_binary);
    CKKSClearDirectRotationKeyIndexSets<P>(usage.slot_to_coeff_direct);

    for (std::size_t i = 0; i < linear_plan.coeff_to_slot_stages.size(); i++) {
        CKKSCollectLinearTransformStageHybridGiantRotationKeyIndices<P>(
            usage.coeff_to_slot_binary[i], usage.coeff_to_slot_binary[i + 1],
            usage.coeff_to_slot_direct[i + 1],
            linear_plan.coeff_to_slot_stages[i],
            ckks_detail::CKKSDenseBootstrapCoeffToSlotBSGSStepAt<Schedule>(i),
            Schedule::hybrid_giant_direct_popcount_threshold);
    }
    CKKSMarkConjugationKeyIndex<P>(usage.packed_conjugate);
    for (std::size_t i = 0; i < linear_plan.slot_to_coeff_stages.size(); i++) {
        CKKSCollectLinearTransformStageHybridGiantRotationKeyIndices<P>(
            usage.slot_to_coeff_binary[i], usage.slot_to_coeff_binary[i + 1],
            usage.slot_to_coeff_direct[i + 1],
            linear_plan.slot_to_coeff_stages[i],
            ckks_detail::CKKSDenseBootstrapSlotToCoeffBSGSStepAt<Schedule>(i),
            Schedule::hybrid_giant_direct_popcount_threshold);
        CKKSCollectLinearTransformStageHybridGiantRotationKeyIndices<P>(
            usage.slot_to_coeff_binary[i], usage.slot_to_coeff_binary[i + 1],
            usage.slot_to_coeff_direct[i + 1],
            linear_plan.slot_to_coeff_imag_stages[i],
            ckks_detail::CKKSDenseBootstrapSlotToCoeffBSGSStepAt<Schedule>(i),
            Schedule::hybrid_giant_direct_popcount_threshold);
    }
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapRotationKeyUsageCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSRotationKeyIndexSetsCount<P>(usage.coeff_to_slot) +
           CKKSRotationKeyIndexSetCount<P>(usage.packed_conjugate) +
           CKKSRotationKeyIndexSetsCount<P>(usage.slot_to_coeff);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectRotationKeyUsageCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDirectRotationKeyIndexSetsCount<P>(usage.coeff_to_slot) +
           CKKSRotationKeyIndexSetCount<P>(usage.packed_conjugate) +
           CKKSDirectRotationKeyIndexSetsCount<P>(usage.slot_to_coeff);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantRotationKeyUsageCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSRotationKeyIndexSetsCount<P>(usage.coeff_to_slot_binary) +
           CKKSDirectRotationKeyIndexSetsCount<P>(usage.coeff_to_slot_direct) +
           CKKSRotationKeyIndexSetCount<P>(usage.packed_conjugate) +
           CKKSRotationKeyIndexSetsCount<P>(usage.slot_to_coeff_binary) +
           CKKSDirectRotationKeyIndexSetsCount<P>(usage.slot_to_coeff_direct);
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapFullGaloisKeyIndexCount()
{
    using P = typename Schedule::Param;
    return (Schedule::coeff_to_slot_level_count + 1 + 1 +
            Schedule::slot_to_coeff_level_count + 1) *
           (P::nbit + 1);
}

template <class Schedule>
inline void CKKSBuildDenseBootstrapLinearPlan(
    CKKSDenseBootstrapLinearPlan<Schedule> &plan)
{
    using P = typename Schedule::Param;
    CKKSLinearTransformStages<P> raw_coeff_to_slot;
    CKKSLinearTransformStages<P> raw_slot_to_coeff;
    CKKSBuildCoeffToPackedSlotStages<P>(raw_coeff_to_slot);
    CKKSBuildPackedSlotToCoeffStages<P>(raw_slot_to_coeff);

    CKKSFuseLinearTransformStages<P>(plan.coeff_to_slot_stages,
                                     raw_coeff_to_slot,
                                     Schedule::coeff_to_slot_fuse_radix);
    CKKSFuseLinearTransformStages<P>(plan.slot_to_coeff_stages,
                                     raw_slot_to_coeff,
                                     Schedule::slot_to_coeff_fuse_radix);
    CKKSFuseLinearTransformStages<P>(plan.slot_to_coeff_imag_stages,
                                     raw_slot_to_coeff,
                                     Schedule::slot_to_coeff_fuse_radix);
    assert(plan.coeff_to_slot_stages.size() ==
           Schedule::coeff_to_slot_level_count);
    assert(plan.slot_to_coeff_stages.size() ==
           Schedule::slot_to_coeff_level_count);
    assert(plan.slot_to_coeff_imag_stages.size() ==
           Schedule::slot_to_coeff_level_count);

    CKKSScaleLinearTransformStages<P>(
        plan.coeff_to_slot_stages,
        {Schedule::coeff_to_slot_scaling_factor, 0.0});
    CKKSScaleLinearTransformStages<P>(
        plan.slot_to_coeff_stages,
        {Schedule::slot_to_coeff_scaling_factor, 0.0});
    CKKSScaleLinearTransformStages<P>(
        plan.slot_to_coeff_imag_stages,
        {0.0, Schedule::slot_to_coeff_scaling_factor});
}

template <class P, std::uint32_t LogQ>
inline std::unique_ptr<CKKSRelinKey<P, LogQ>> makeCKKSRelinKey(
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    static_assert(ckks_detail::supported_torus_v<P>);

    auto relinkey = std::make_unique<CKKSRelinKey<P, LogQ>>();
    auto keysquare = std::make_unique<Polynomial<P>>();
    auto partkey = std::make_unique<Polynomial<P>>();
    for (std::uint32_t i = 0; i < P::n; i++) (*partkey)[i] = key[i];
    if constexpr (is_multilimb_uint_v<typename P::T>)
        PolyMulDigit<P>(*keysquare, *partkey, *partkey);
    else
        PolyMul<P>(*keysquare, *partkey, *partkey);

    constexpr int width = std::numeric_limits<typename P::T>::digits;
    constexpr std::uint32_t row_count =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row =
        CKKSKeySwitchFirstRowForLevel<P, LogQ>();
    auto gadget = std::make_unique<Polynomial<P>>();
    for (int j = 0; j < static_cast<int>(row_count); j++) {
        const std::uint32_t full_row = first_row + j;
        const int shift =
            width - (full_row + 1) * static_cast<int>(P::B̅gbit);
        for (std::uint32_t n = 0; n < P::n; n++)
            (*gadget)[n] = ckks_detail::reduceToLevel<P, LogQ>((*keysquare)[n]
                                                               << shift);
        ckks_detail::encryptPolynomialAtLevel<P, LogQ>((*relinkey)[j], *gadget,
                                                       key, noise);
    }

    return relinkey;
}

template <class P, std::uint32_t LogQ>
inline std::unique_ptr<CKKSSeededRelinKey<P, LogQ>> makeCKKSSeededRelinKey(
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    static_assert(ckks_detail::supported_torus_v<P>);

    auto relinkey = std::make_unique<CKKSSeededRelinKey<P, LogQ>>();
    auto keysquare = std::make_unique<Polynomial<P>>();
    auto partkey = std::make_unique<Polynomial<P>>();
    for (std::uint32_t i = 0; i < P::n; i++) (*partkey)[i] = key[i];
    if constexpr (is_multilimb_uint_v<typename P::T>)
        PolyMulDigit<P>(*keysquare, *partkey, *partkey);
    else
        PolyMul<P>(*keysquare, *partkey, *partkey);

    constexpr int width = std::numeric_limits<typename P::T>::digits;
    constexpr std::uint32_t row_count =
        CKKSKeySwitchRowCountForLevel<P, LogQ>();
    constexpr std::uint32_t first_row =
        CKKSKeySwitchFirstRowForLevel<P, LogQ>();
    auto gadget = std::make_unique<Polynomial<P>>();
    for (int j = 0; j < static_cast<int>(row_count); j++) {
        const std::uint32_t full_row = first_row + j;
        const int shift =
            width - (full_row + 1) * static_cast<int>(P::B̅gbit);
        for (std::uint32_t n = 0; n < P::n; n++)
            (*gadget)[n] = ckks_detail::reduceToLevel<P, LogQ>((*keysquare)[n]
                                                               << shift);
        ckks_detail::encryptPolynomialAtLevel<P, LogQ>((*relinkey)[j], *gadget,
                                                       key, noise);
    }

    return relinkey;
}

namespace ckks_detail {

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          class Seq>
struct CKKSRelinKeyChainTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t... Is>
struct CKKSRelinKeyChainTuple<P, StartLogQ, LogDelta,
                              std::index_sequence<Is...>> {
    using type =
        std::tuple<CKKSRelinKey<P, StartLogQ - (Is + 1) * LogDelta>...>;
};

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t Depth>
struct CKKSRelinKeyChain {
    static_assert(Depth > 0);
    static_assert(StartLogQ > Depth * LogDelta);

    using Tuple = typename ckks_detail::CKKSRelinKeyChainTuple<
        P, StartLogQ, LogDelta, std::make_index_sequence<Depth>>::type;

    Tuple keys{};

    template <std::size_t I>
    auto &get()
    {
        static_assert(I < Depth);
        return std::get<I>(keys);
    }

    template <std::size_t I>
    const auto &get() const
    {
        static_assert(I < Depth);
        return std::get<I>(keys);
    }

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t I>
using CKKSRelinKeyChainElement =
    CKKSRelinKey<P, StartLogQ - (I + 1) * LogDelta>;

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::size_t Depth>
inline void CKKSRelinKeyChainElementGen(
    CKKSRelinKeyChainElement<P, StartLogQ, LogDelta, I> &key_out,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    static_assert(I < Depth);
    key_out = *makeCKKSRelinKey<
        P, StartLogQ - (I + 1) * LogDelta>(key, noise);
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::size_t Depth>
inline void CKKSSeededRelinKeyChainElementGen(
    CKKSSeededRelinKey<P, StartLogQ - (I + 1) * LogDelta> &key_out,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    static_assert(I < Depth);
    key_out = *makeCKKSSeededRelinKey<
        P, StartLogQ - (I + 1) * LogDelta>(key, noise);
}

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::size_t Depth>
inline void CKKSRelinKeyChainGenImpl(
    CKKSRelinKeyChain<P, StartLogQ, LogDelta, Depth> &chain,
    const Key<P> &key, CKKSNoise noise)
{
    if constexpr (I < Depth) {
        constexpr std::uint32_t log_q = StartLogQ - (I + 1) * LogDelta;
        chain.template get<I>() = *makeCKKSRelinKey<P, log_q>(key, noise);
        CKKSRelinKeyChainGenImpl<I + 1, P, StartLogQ, LogDelta, Depth>(
            chain, key, noise);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t Depth>
inline void CKKSRelinKeyChainGen(
    CKKSRelinKeyChain<P, StartLogQ, LogDelta, Depth> &chain,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    ckks_detail::CKKSRelinKeyChainGenImpl<0, P, StartLogQ, LogDelta, Depth>(
        chain, key, noise);
}

template <class P, std::uint32_t LogQ, class RelinKey>
    requires(RelinKey::log_q == LogQ)
inline void CKKSRelinKeySwitch(TRLWE<P> &res, const Polynomial<P> &poly,
                               const RelinKey &relinkey)
{
    CKKSKeySwitchRows<P, LogQ>(res, poly, relinkey.data);
}

template <class P, std::uint32_t LogQ, class RelinKey>
    requires(RelinKey::log_q == LogQ)
inline void CKKSRelinearization(TRLWE<P> &res, const TRLWE3<P> &mult,
                                const RelinKey &relinkey)
{
    auto squareterm = std::make_unique<TRLWE<P>>();
    CKKSRelinKeySwitch<P, LogQ>(*squareterm, mult[2], relinkey);
    for (std::uint32_t i = 0; i < P::n; i++)
        res[0][i] = ckks_detail::reduceToLevel<P, LogQ>(mult[0][i] +
                                                        (*squareterm)[0][i]);
    for (std::uint32_t i = 0; i < P::n; i++)
        res[1][i] = ckks_detail::reduceToLevel<P, LogQ>(mult[1][i] +
                                                        (*squareterm)[1][i]);
}

template <class P, std::uint32_t LhsLogQ, std::uint32_t RhsLogQ,
          std::uint32_t LogScale>
inline void CKKSTensorProductRescale(TRLWE3<P> &res, const TRLWE<P> &a,
                                     const TRLWE<P> &b)
{
    static_assert(ckks_detail::supported_torus_v<P>);
    constexpr int width = std::numeric_limits<typename P::T>::digits;
    static_assert(P::l̅ * P::B̅gbit == width,
                  "CKKS FullDD multiply requires Bbar digits to cover T");
    static_assert(LogScale < static_cast<std::uint32_t>(width));
    constexpr std::uint32_t base_log_q =
        ckks_detail::static_min_v<LhsLogQ, RhsLogQ>;
    static_assert(base_log_q > LogScale);
    constexpr std::uint32_t out_log_q = base_log_q - LogScale;
    constexpr std::uint32_t lhs_row_count =
        CKKSKeySwitchRowCountForLevel<P, LhsLogQ>();
    constexpr std::uint32_t rhs_row_count =
        CKKSKeySwitchRowCountForLevel<P, RhsLogQ>();
    constexpr std::uint32_t lhs_first_row =
        CKKSKeySwitchFirstRowForLevel<P, LhsLogQ>();
    constexpr std::uint32_t rhs_first_row =
        CKKSKeySwitchFirstRowForLevel<P, RhsLogQ>();

    auto a_centered = std::make_unique<TRLWE<P>>();
    auto b_centered = std::make_unique<TRLWE<P>>();
    ckks_detail::centeredTRLWEAtLevel<P, LhsLogQ>(*a_centered, a);
    ckks_detail::centeredTRLWEAtLevel<P, RhsLogQ>(*b_centered, b);

    auto a_dec = std::make_unique<std::array<TRLWE<P>, lhs_row_count>>();
    auto b_dec = std::make_unique<std::array<TRLWE<P>, rhs_row_count>>();
    ckks_detail::baseBbarDecomposeRows<P, lhs_first_row>(*a_dec, *a_centered);
    ckks_detail::baseBbarDecomposeRows<P, rhs_first_row>(*b_dec, *b_centered);

    if constexpr (std::is_same_v<typename P::T, __uint128_t>) {
        constexpr bool use_fft_digits =
            2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3 <
            std::numeric_limits<double>::digits;

        auto acc = std::make_unique<std::array<std::array<Wide384, P::n>, 3>>();

        if constexpr (use_fft_digits) {
            auto a_fd =
                std::make_unique<std::array<TRLWEInFD<P>, lhs_row_count>>();
            auto b_fd =
                std::make_unique<std::array<TRLWEInFD<P>, rhs_row_count>>();
            for (int i = 0; i < static_cast<int>(lhs_row_count); i++) {
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    TwistIFFT<P>((*a_fd)[i][c], (*a_dec)[i][c]);
                }
            }
            for (int i = 0; i < static_cast<int>(rhs_row_count); i++) {
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    TwistIFFT<P>((*b_fd)[i][c], (*b_dec)[i][c]);
                }
            }

            alignas(64) PolynomialInFD<P> prod_fd;
            Polynomial<P> prod;
            for (int i = 0; i < static_cast<int>(lhs_row_count); i++) {
                for (int j = 0; j < static_cast<int>(rhs_row_count); j++) {
                    const int digit_sum =
                        static_cast<int>(lhs_first_row + rhs_first_row) + i +
                        j;
                    const int shift =
                        2 * width -
                        (digit_sum + 2) * static_cast<int>(P::B̅gbit);
                    if (shift <= -width) continue;

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][0], (*b_fd)[j][0]);
                    TwistFFT<P>(prod, prod_fd);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[2][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][1], (*b_fd)[j][1]);
                    TwistFFT<P>(prod, prod_fd);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[1][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);

                    MulInFD<P::n>(prod_fd, (*a_fd)[i][0], (*b_fd)[j][1]);
                    FMAInFD<P::n>(prod_fd, (*a_fd)[i][1], (*b_fd)[j][0]);
                    TwistFFT<P>(prod, prod_fd);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                }
            }
        }
        else {
            Polynomial<P> prod;
            for (int i = 0; i < static_cast<int>(lhs_row_count); i++) {
                for (int j = 0; j < static_cast<int>(rhs_row_count); j++) {
                    const int digit_sum =
                        static_cast<int>(lhs_first_row + rhs_first_row) + i +
                        j;
                    const int shift =
                        2 * width -
                        (digit_sum + 2) * static_cast<int>(P::B̅gbit);
                    if (shift <= -width) continue;

                    PolyMulNaive<P>(prod, (*a_dec)[i][0], (*b_dec)[j][0]);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[2][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][1], (*b_dec)[j][1]);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[1][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][0], (*b_dec)[j][1]);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                    PolyMulNaive<P>(prod, (*a_dec)[i][1], (*b_dec)[j][0]);
                    for (std::uint32_t n = 0; n < P::n; n++)
                        (*acc)[0][n].add_shifted(
                            static_cast<__int128_t>(prod[n]), shift);
                }
            }
        }

        for (int c = 0; c < 3; c++)
            for (std::uint32_t n = 0; n < P::n; n++)
                res[c][n] = ckks_detail::reduceToLevel<P, out_log_q>(
                    ckks_detail::rescaleAccumulator<P, LogScale>((*acc)[c][n]));
    }
    else {
        using Torus = typename P::T;
        using Acc = WideSignedLimbAccumulator<2 * Torus::limbs + 2>;
        auto acc = std::make_unique<std::array<std::array<Acc, P::n>, 3>>();

        auto accumulate = [](std::array<std::array<Acc, P::n>, 3> &accum,
                             int component, const Polynomial<P> &prod,
                             int shift) {
            for (std::uint32_t n = 0; n < P::n; n++)
                accum[component][n].add_shifted_i64(
                    multilimb_to_signed_i64(prod[n]), shift);
        };

        auto merge_accumulator =
            [](std::array<std::array<Acc, P::n>, 3> &dst,
               const std::array<std::array<Acc, P::n>, 3> &src) {
                for (int c = 0; c < 3; c++)
                    for (std::uint32_t n = 0; n < P::n; n++)
                        dst[c][n].add(src[c][n]);
            };

        constexpr bool use_fft_digits =
            use_multilimb_digit_fft_v<P> &&
            2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 4 <
                std::numeric_limits<double>::digits;

        if constexpr (use_fft_digits) {
            constexpr int digit_product_bits =
                2 * static_cast<int>(P::B̅gbit) + static_cast<int>(P::nbit) + 3;
            constexpr int cross_term_product_bits = digit_product_bits + 1;
            constexpr int fd_batch_slack =
                std::numeric_limits<double>::digits - cross_term_product_bits;
            static_assert(fd_batch_slack > 0);
            constexpr int fd_batch_size =
                std::min(static_cast<int>(ckks_detail::static_min_v<
                             lhs_row_count, rhs_row_count>),
                         1 << fd_batch_slack);

            auto a_fd =
                std::make_unique<std::array<TRLWEInFD<P>, lhs_row_count>>();
            auto b_fd =
                std::make_unique<std::array<TRLWEInFD<P>, rhs_row_count>>();
            for (int i = 0; i < static_cast<int>(lhs_row_count); i++) {
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    TwistIFFTDigit<P>((*a_fd)[i][c], (*a_dec)[i][c]);
                }
            }
            for (int i = 0; i < static_cast<int>(rhs_row_count); i++) {
                for (int c = 0; c <= static_cast<int>(P::k); c++) {
                    TwistIFFTDigit<P>((*b_fd)[i][c], (*b_dec)[i][c]);
                }
            }

#pragma omp parallel
            {
                auto local_acc =
                    std::make_unique<std::array<std::array<Acc, P::n>, 3>>();
                auto sum_fd =
                    std::make_unique<std::array<PolynomialInFD<P>, 3>>();
                auto prod = std::make_unique<Polynomial<P>>();
#pragma omp for schedule(dynamic)
                for (int digit_sum = 0;
                     digit_sum <=
                         static_cast<int>(lhs_row_count + rhs_row_count) - 2;
                     digit_sum++) {
                    const int i_begin = std::max(
                        0,
                        digit_sum - static_cast<int>(rhs_row_count) + 1);
                    const int i_end =
                        std::min(static_cast<int>(lhs_row_count) - 1,
                                 digit_sum);
                    const int full_digit_sum =
                        static_cast<int>(lhs_first_row + rhs_first_row) +
                        digit_sum;
                    const int shift =
                        2 * width -
                        (full_digit_sum + 2) * static_cast<int>(P::B̅gbit);

                    for (int batch_begin = i_begin; batch_begin <= i_end;
                         batch_begin += fd_batch_size) {
                        for (int c = 0; c < 3; c++) (*sum_fd)[c].fill(0.0);
                        const int batch_end =
                            std::min(i_end, batch_begin + fd_batch_size - 1);
                        for (int i = batch_begin; i <= batch_end; i++) {
                            const int j = digit_sum - i;
                            FMAInFD<P::n>((*sum_fd)[2], (*a_fd)[i][0],
                                          (*b_fd)[j][0]);
                            FMAInFD<P::n>((*sum_fd)[1], (*a_fd)[i][1],
                                          (*b_fd)[j][1]);
                            FMAInFD<P::n>((*sum_fd)[0], (*a_fd)[i][0],
                                          (*b_fd)[j][1]);
                            FMAInFD<P::n>((*sum_fd)[0], (*a_fd)[i][1],
                                          (*b_fd)[j][0]);
                        }
                        for (int c = 0; c < 3; c++) {
                            TwistFFTDigitProduct<P>(*prod, (*sum_fd)[c]);
                            accumulate(*local_acc, c, *prod, shift);
                        }
                    }
                }
#pragma omp critical
                { merge_accumulator(*acc, *local_acc); }
            }
        }
        else {
#pragma omp parallel
            {
                auto local_acc =
                    std::make_unique<std::array<std::array<Acc, P::n>, 3>>();
                auto prod = std::make_unique<Polynomial<P>>();
#pragma omp for collapse(2) schedule(dynamic)
                for (int i = 0; i < static_cast<int>(lhs_row_count); i++) {
                    for (int j = 0; j < static_cast<int>(rhs_row_count); j++) {
                        const int digit_sum =
                            static_cast<int>(lhs_first_row + rhs_first_row) +
                            i + j;
                        const int shift =
                            2 * width -
                            (digit_sum + 2) * static_cast<int>(P::B̅gbit);

                        PolyMulDigit<P>(*prod, (*a_dec)[i][0],
                                        (*b_dec)[j][0]);
                        accumulate(*local_acc, 2, *prod, shift);
                        PolyMulDigit<P>(*prod, (*a_dec)[i][1],
                                        (*b_dec)[j][1]);
                        accumulate(*local_acc, 1, *prod, shift);
                        PolyMulDigit<P>(*prod, (*a_dec)[i][0],
                                        (*b_dec)[j][1]);
                        accumulate(*local_acc, 0, *prod, shift);
                        PolyMulDigit<P>(*prod, (*a_dec)[i][1],
                                        (*b_dec)[j][0]);
                        accumulate(*local_acc, 0, *prod, shift);
                    }
                }
#pragma omp critical
                { merge_accumulator(*acc, *local_acc); }
            }
        }

        const Torus divisor =
            LogScale == 0 ? Torus{1} : (Torus{1} << LogScale);
        for (int c = 0; c < 3; c++) {
            for (std::uint32_t n = 0; n < P::n; n++) {
                if constexpr (LogScale > 0) {
                    (*acc)[c][n].add_shifted_i64(
                        (*acc)[c][n].is_negative() ? -1 : 1,
                        static_cast<int>(LogScale) - 1);
                }
                res[c][n] = ckks_detail::reduceToLevel<P, out_log_q>(
                    (*acc)[c][n].template div_to_torus<Torus::limbs>(divisor));
            }
        }
    }
}

template <class P, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
struct CKKSMultTraits {
    static constexpr std::uint32_t log_scale =
        ckks_detail::static_max_v<LhsLogDelta, RhsLogDelta>;
    static constexpr std::uint32_t base_log_q =
        ckks_detail::static_min_v<LhsLogQ, RhsLogQ>;
    static_assert(base_log_q > log_scale);

    static constexpr std::uint32_t log_q = base_log_q - log_scale;
    static constexpr std::uint32_t log_delta =
        ckks_detail::static_min_v<LhsLogDelta, RhsLogDelta>;

    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
using CKKSMultResult =
    typename CKKSMultTraits<P, LhsLogQ, LhsLogDelta, RhsLogQ,
                            RhsLogDelta>::Ciphertext;

template <class P, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta, class RelinKey>
    requires(RelinKey::log_q ==
             CKKSMultTraits<P, LhsLogQ, LhsLogDelta, RhsLogQ,
                            RhsLogDelta>::log_q)
inline void CKKSMult(
    CKKSMultResult<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta> &res,
    const CKKSCiphertext<P, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<P, RhsLogQ, RhsLogDelta> &rhs,
    const RelinKey &relinkey)
{
    using Traits =
        CKKSMultTraits<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta>;
    const bool trace_evalmod = ckks_detail::ckks_evalmod_trace_enabled();
    std::chrono::steady_clock::time_point trace_total_start{};
    std::chrono::steady_clock::time_point trace_section_start{};
    if (trace_evalmod) {
        trace_total_start = std::chrono::steady_clock::now();
        trace_section_start = trace_total_start;
    }

    auto mult = std::make_unique<TRLWE3<P>>();
    CKKSTensorProductRescale<P, LhsLogQ, RhsLogQ, Traits::log_scale>(
        *mult, lhs.ct, rhs.ct);
    if (trace_evalmod) {
        ckks_detail::ckks_trace_evalmod_mult_event(
            "mult_tensor_product", LhsLogQ, RhsLogQ, Traits::log_q,
            ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
        trace_section_start = std::chrono::steady_clock::now();
    }
    CKKSRelinearization<P, Traits::log_q>(res.ct, *mult, relinkey);
    if (trace_evalmod) {
        ckks_detail::ckks_trace_evalmod_mult_event(
            "mult_relinearization", LhsLogQ, RhsLogQ, Traits::log_q,
            ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
        ckks_detail::ckks_trace_evalmod_mult_event(
            "mult_total", LhsLogQ, RhsLogQ, Traits::log_q,
            ckks_detail::ckks_trace_elapsed_ms(trace_total_start));
    }
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta, class RelinKey>
    requires(RelinKey::log_q ==
             CKKSMultTraits<P, LogQ, LogDelta, LogQ, LogDelta>::log_q)
inline void CKKSSquare(
    CKKSMultResult<P, LogQ, LogDelta, LogQ, LogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct, const RelinKey &relinkey)
{
    CKKSMult<P>(res, ct, ct, relinkey);
}

namespace ckks_detail {

constexpr std::uint32_t power_tree_depth(std::size_t power)
{
    return power <= 1 ? 0 : bit_width_u64(power - 1);
}

constexpr std::size_t choose_power_bsgs_baby_step(std::size_t degree)
{
    if (degree < 16) return 0;

    const std::uint32_t total_depth = power_tree_depth(degree);
    std::size_t best_step = 0;
    std::size_t best_cost = degree - 1;
    for (std::size_t step = 2; step < degree; step <<= 1) {
        const std::size_t giant_steps = degree / step;
        if (giant_steps == 0 || giant_steps > total_depth) continue;

        const std::uint32_t block_input_depth =
            total_depth - static_cast<std::uint32_t>(giant_steps);
        if (power_tree_depth(step - 1) > block_input_depth) continue;

        const std::size_t cost = (step - 1) + giant_steps;
        if (cost < best_cost) {
            best_step = step;
            best_cost = cost;
        }
    }
    return best_step;
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t Degree, class Seq>
struct CKKSPowerBasisTuple;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::size_t Degree, std::size_t... Is>
struct CKKSPowerBasisTuple<P, StartLogQ, LogDelta, Degree,
                           std::index_sequence<Is...>> {
    using type = std::tuple<
        std::unique_ptr<CKKSCiphertext<
            P, StartLogQ - power_tree_depth(Is + 1) * LogDelta, LogDelta>>...>;
};

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
struct CKKSPowerPolynomialEvaluatorTraits {
    static_assert(Degree > 1);
    static constexpr std::uint32_t power_depth =
        ckks_detail::power_tree_depth(Degree);
    static constexpr std::size_t bsgs_baby_step =
        CoeffLogDelta == LogDelta
            ? ckks_detail::choose_power_bsgs_baby_step(Degree)
            : 0;
    static constexpr bool bsgs_enabled = bsgs_baby_step != 0;
    static constexpr std::uint32_t relin_depth =
        power_depth + (bsgs_enabled ? 1 : 0);
    static_assert(StartLogQ > power_depth * LogDelta + CoeffLogDelta);

    static constexpr std::uint32_t term_input_log_q =
        StartLogQ - power_depth * LogDelta;
    static constexpr std::uint32_t log_q = term_input_log_q - CoeffLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;

    using PowerBasis = typename ckks_detail::CKKSPowerBasisTuple<
        P, StartLogQ, LogDelta, Degree,
        std::make_index_sequence<Degree>>::type;
    using RelinKeyChain = CKKSRelinKeyChain<P, StartLogQ, LogDelta, relin_depth>;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
using CKKSPowerPolynomialResult =
    typename CKKSPowerPolynomialEvaluatorTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>::Ciphertext;

namespace ckks_detail {

template <std::size_t I, class Traits, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t CoeffLogDelta,
          std::size_t Degree, class RelinKeyProvider>
inline void CKKSBuildPowerBasisImpl(
    typename Traits::PowerBasis &powers,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const RelinKeyProvider &keys)
{
    if constexpr (I <= Degree) {
        if constexpr (I == 1) {
            using Ct = CKKSCiphertext<P, StartLogQ, LogDelta>;
            std::get<0>(powers) = std::make_unique<Ct>(ct);
        }
        else {
            constexpr std::size_t lhs_power = I / 2;
            constexpr std::size_t rhs_power = I - lhs_power;
            constexpr std::uint32_t depth = power_tree_depth(I);
            using Ct = CKKSCiphertext<P, StartLogQ - depth * LogDelta,
                                      LogDelta>;
            auto &dst = std::get<I - 1>(powers);
            dst = std::make_unique<Ct>();
            const bool trace_evalmod = ckks_evalmod_trace_enabled();
            std::chrono::steady_clock::time_point trace_start{};
            if (trace_evalmod) trace_start = std::chrono::steady_clock::now();
            CKKSMult<P>(*dst, *std::get<lhs_power - 1>(powers),
                        *std::get<rhs_power - 1>(powers),
                        keys.template get<depth - 1>());
            if (trace_evalmod)
                ckks_trace_evalmod_power_event(
                    "power_basis_mult", I, depth,
                    ckks_trace_elapsed_ms(trace_start));
            if constexpr (I == Degree || power_tree_depth(I + 1) > depth) {
                maybe_release_key<depth - 1>(keys);
            }
        }
        CKKSBuildPowerBasisImpl<I + 1, Traits, P, StartLogQ, LogDelta,
                                CoeffLogDelta, Degree>(powers, ct, keys);
    }
}

template <std::size_t I, class Traits, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t CoeffLogDelta,
          std::size_t Degree>
inline void CKKSAddPowerPolynomialTermsImpl(
    typename Traits::Ciphertext &res, const typename Traits::PowerBasis &powers,
    const std::vector<double> &coeffs)
{
    if constexpr (I <= Degree) {
        if (I < coeffs.size() && coeffs[I] != 0.0) {
            using PowerPtr =
                std::tuple_element_t<I - 1, typename Traits::PowerBasis>;
            using PowerCt = typename PowerPtr::element_type;
            const auto &power = std::get<I - 1>(powers);
            assert(power != nullptr);
            auto reduced = std::make_unique<
                CKKSCiphertext<P, Traits::term_input_log_q, LogDelta>>();
            CKKSLevelReduce<P, PowerCt::log_q, Traits::term_input_log_q,
                            LogDelta>(*reduced, *power);

            auto term = std::make_unique<
                CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta,
                                   CoeffLogDelta>>();
            CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta,
                               CoeffLogDelta>(*term, *reduced, coeffs[I]);
            CKKSAddInPlace<P, Traits::log_q, LogDelta>(res, *term);
        }
        CKKSAddPowerPolynomialTermsImpl<I + 1, Traits, P, StartLogQ, LogDelta,
                                        CoeffLogDelta, Degree>(res, powers,
                                                               coeffs);
    }
}

template <class EvalTraits, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t CoeffLogDelta,
          std::size_t Degree>
struct CKKSPowerPolynomialBSGSTraits {
    static_assert(EvalTraits::bsgs_enabled);
    static_assert(CoeffLogDelta == LogDelta);
    static constexpr std::size_t baby_step = EvalTraits::bsgs_baby_step;
    static constexpr std::size_t giant_steps = Degree / baby_step;
    static_assert(giant_steps > 0);
    static constexpr std::uint32_t block_input_depth =
        EvalTraits::power_depth - static_cast<std::uint32_t>(giant_steps);
    static_assert(power_tree_depth(baby_step - 1) <= block_input_depth);

    static constexpr std::uint32_t block_input_log_q =
        StartLogQ - block_input_depth * LogDelta;
    static constexpr std::uint32_t block_log_q =
        block_input_log_q - CoeffLogDelta;

    using PowerBasis = typename CKKSPowerBasisTuple<
        P, StartLogQ, LogDelta, baby_step,
        std::make_index_sequence<baby_step>>::type;
    using BlockCiphertext = CKKSCiphertext<P, block_log_q, LogDelta>;
    using GiantPtr = std::tuple_element_t<baby_step - 1, PowerBasis>;
    using GiantCiphertext = typename GiantPtr::element_type;
    using FinalCiphertext = typename EvalTraits::Ciphertext;
};

template <std::size_t Block, std::size_t I, class BSGS, class P,
          std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
inline void CKKSAddPowerPolynomialBSGSBlockTermsImpl(
    typename BSGS::BlockCiphertext &res,
    const typename BSGS::PowerBasis &powers, const std::vector<double> &coeffs)
{
    if constexpr (I < BSGS::baby_step) {
        constexpr std::size_t coeff_index = Block * BSGS::baby_step + I;
        if constexpr (coeff_index <= Degree) {
            if (coeff_index < coeffs.size() && coeffs[coeff_index] != 0.0) {
                using PowerPtr =
                    std::tuple_element_t<I - 1, typename BSGS::PowerBasis>;
                using PowerCt = typename PowerPtr::element_type;
                static_assert(BSGS::block_input_log_q <= PowerCt::log_q);
                const auto &power = std::get<I - 1>(powers);
                assert(power != nullptr);
                auto reduced = std::make_unique<
                    CKKSCiphertext<P, BSGS::block_input_log_q, LogDelta>>();
                CKKSLevelReduce<P, PowerCt::log_q, BSGS::block_input_log_q,
                                LogDelta>(*reduced, *power);

                auto term = std::make_unique<CKKSPlainMulResult<
                    P, BSGS::block_input_log_q, LogDelta, CoeffLogDelta>>();
                CKKSPlainMulByReal<P, BSGS::block_input_log_q, LogDelta,
                                   CoeffLogDelta>(*term, *reduced,
                                                  coeffs[coeff_index]);
                CKKSAddInPlace<P, BSGS::block_log_q, LogDelta>(res, *term);
            }
        }
        CKKSAddPowerPolynomialBSGSBlockTermsImpl<
            Block, I + 1, BSGS, P, StartLogQ, LogDelta, CoeffLogDelta,
            Degree>(res, powers, coeffs);
    }
}

template <std::size_t Block, class BSGS, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t CoeffLogDelta,
          std::size_t Degree>
inline void CKKSEvalPowerPolynomialBSGSBlock(
    typename BSGS::BlockCiphertext &res,
    const typename BSGS::PowerBasis &powers, const std::vector<double> &coeffs)
{
    constexpr std::size_t constant_index = Block * BSGS::baby_step;
    CKKSSetTransparentReal<P, BSGS::block_log_q, LogDelta>(
        res, constant_index < coeffs.size() ? coeffs[constant_index] : 0.0);
    CKKSAddPowerPolynomialBSGSBlockTermsImpl<
        Block, 1, BSGS, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
        res, powers, coeffs);
}

template <std::size_t Remaining, class EvalTraits, class BSGS, class P,
          std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree, class CurrentCt,
          class RelinKeyProvider>
inline void CKKSPowerPolynomialBSGSHornerImpl(
    typename BSGS::FinalCiphertext &res, const CurrentCt &current,
    const typename BSGS::GiantCiphertext &giant,
    const typename BSGS::PowerBasis &powers, const std::vector<double> &coeffs,
    const RelinKeyProvider &keys)
{
    if constexpr (Remaining == 0) {
        static_assert(CurrentCt::log_q == BSGS::FinalCiphertext::log_q);
        res = current;
    }
    else {
        using ProductCt =
            CKKSMultResult<P, CurrentCt::log_q, LogDelta,
                           BSGS::GiantCiphertext::log_q, LogDelta>;
        constexpr std::uint32_t consumed = StartLogQ - ProductCt::log_q;
        static_assert(consumed % LogDelta == 0);
        constexpr std::size_t key_index = consumed / LogDelta - 1;
        static_assert(key_index < EvalTraits::relin_depth);

        auto product = std::make_unique<ProductCt>();
        const bool trace_evalmod = ckks_evalmod_trace_enabled();
        std::chrono::steady_clock::time_point trace_start{};
        if (trace_evalmod) trace_start = std::chrono::steady_clock::now();
        CKKSMult<P>(*product, current, giant, keys.template get<key_index>());
        if (trace_evalmod)
            ckks_trace_evalmod_power_event(
                "power_bsgs_giant_mult",
                (BSGS::giant_steps - Remaining + 1) * BSGS::baby_step,
                consumed / LogDelta, ckks_trace_elapsed_ms(trace_start));
        maybe_release_key<key_index>(keys);

        constexpr std::size_t block_index = Remaining - 1;
        auto block = std::make_unique<typename BSGS::BlockCiphertext>();
        CKKSEvalPowerPolynomialBSGSBlock<
            block_index, BSGS, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
            *block, powers, coeffs);
        auto reduced = std::make_unique<ProductCt>();
        CKKSLevelReduce<P, BSGS::block_log_q, ProductCt::log_q, LogDelta>(
            *reduced, *block);
        CKKSAddInPlace<P, ProductCt::log_q, LogDelta>(*product, *reduced);

        CKKSPowerPolynomialBSGSHornerImpl<
            Remaining - 1, EvalTraits, BSGS, P, StartLogQ, LogDelta,
            CoeffLogDelta, Degree>(res, *product, giant, powers, coeffs, keys);
    }
}

template <class EvalTraits, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t CoeffLogDelta,
          std::size_t Degree, class RelinKeyProvider>
inline void CKKSEvalPowerPolynomialBSGSWithKeyProvider(
    typename EvalTraits::Ciphertext &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const std::vector<double> &coeffs, const RelinKeyProvider &keys)
{
    using BSGS = CKKSPowerPolynomialBSGSTraits<
        EvalTraits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>;
    typename BSGS::PowerBasis powers;
    CKKSBuildPowerBasisImpl<1, BSGS, P, StartLogQ, LogDelta, CoeffLogDelta,
                            BSGS::baby_step>(powers, ct, keys);

    auto current = std::make_unique<typename BSGS::BlockCiphertext>();
    CKKSEvalPowerPolynomialBSGSBlock<
        BSGS::giant_steps, BSGS, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
        *current, powers, coeffs);
    const auto &giant = std::get<BSGS::baby_step - 1>(powers);
    assert(giant != nullptr);
    CKKSPowerPolynomialBSGSHornerImpl<
        BSGS::giant_steps, EvalTraits, BSGS, P, StartLogQ, LogDelta,
        CoeffLogDelta, Degree>(res, *current, *giant, powers, coeffs, keys);
}

template <std::size_t I, class Traits, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t CoeffLogDelta,
          std::size_t Degree, class RelinKeyProvider>
inline void CKKSBuildChebyshevBasisImpl(
    typename Traits::PowerBasis &powers,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const RelinKeyProvider &keys)
{
    if constexpr (I <= Degree) {
        if constexpr (I == 1) {
            using Ct = CKKSCiphertext<P, StartLogQ, LogDelta>;
            std::get<0>(powers) = std::make_unique<Ct>(ct);
        }
        else {
            constexpr std::size_t lhs_power = I / 2;
            constexpr std::size_t rhs_power = I - lhs_power;
            constexpr std::size_t correction_power =
                rhs_power - lhs_power;
            constexpr std::uint32_t depth = power_tree_depth(I);
            using Ct = CKKSCiphertext<P, StartLogQ - depth * LogDelta,
                                      LogDelta>;
            auto &dst = std::get<I - 1>(powers);
            dst = std::make_unique<Ct>();
            const bool trace_evalmod = ckks_evalmod_trace_enabled();
            std::chrono::steady_clock::time_point trace_start{};
            if (trace_evalmod) trace_start = std::chrono::steady_clock::now();
            CKKSMult<P>(*dst, *std::get<lhs_power - 1>(powers),
                        *std::get<rhs_power - 1>(powers),
                        keys.template get<depth - 1>());
            if (trace_evalmod)
                ckks_trace_evalmod_power_event(
                    "chebyshev_basis_mult", I, depth,
                    ckks_trace_elapsed_ms(trace_start));
            CKKSMulIntegerInPlace<P, Ct::log_q, LogDelta>(*dst, 2);
            if constexpr (correction_power == 0) {
                CKKSSubPlainRealInPlace<P, Ct::log_q, LogDelta>(*dst, 1.0);
            }
            else {
                using CorrectionPtr = std::tuple_element_t<
                    correction_power - 1, typename Traits::PowerBasis>;
                using CorrectionCt = typename CorrectionPtr::element_type;
                auto correction =
                    std::make_unique<CKKSCiphertext<P, Ct::log_q, LogDelta>>();
                CKKSLevelReduce<P, CorrectionCt::log_q, Ct::log_q, LogDelta>(
                    *correction, *std::get<correction_power - 1>(powers));
                CKKSSubInPlace<P, Ct::log_q, LogDelta>(*dst, *correction);
            }
            if constexpr (I == Degree || power_tree_depth(I + 1) > depth) {
                maybe_release_key<depth - 1>(keys);
            }
        }
        CKKSBuildChebyshevBasisImpl<I + 1, Traits, P, StartLogQ, LogDelta,
                                    CoeffLogDelta, Degree>(powers, ct, keys);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          class RelinKeyProvider>
inline void CKKSEvalPowerPolynomialWithKeyProvider(
    CKKSPowerPolynomialResult<P, StartLogQ, LogDelta, CoeffLogDelta, Degree>
        &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const std::vector<double> &coeffs, const RelinKeyProvider &keys)
{
    using Traits = CKKSPowerPolynomialEvaluatorTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>;
    if constexpr (Traits::bsgs_enabled) {
        ckks_detail::CKKSEvalPowerPolynomialBSGSWithKeyProvider<
            Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
            res, ct, coeffs, keys);
    }
    else {
        typename Traits::PowerBasis powers;
        ckks_detail::CKKSBuildPowerBasisImpl<1, Traits, P, StartLogQ, LogDelta,
                                             CoeffLogDelta, Degree>(powers, ct,
                                                                    keys);
        CKKSSetTransparentReal<P, Traits::log_q, LogDelta>(
            res, coeffs.empty() ? 0.0 : coeffs[0]);
        ckks_detail::CKKSAddPowerPolynomialTermsImpl<
            1, Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(
            res, powers, coeffs);
    }
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree>
inline void CKKSEvalPowerPolynomial(
    CKKSPowerPolynomialResult<P, StartLogQ, LogDelta, CoeffLogDelta, Degree>
        &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const std::vector<double> &coeffs,
    const typename CKKSPowerPolynomialEvaluatorTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>::RelinKeyChain &keys)
{
    CKKSEvalPowerPolynomialWithKeyProvider<P, StartLogQ, LogDelta,
                                           CoeffLogDelta, Degree>(
        res, ct, coeffs, keys);
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          class RelinKeyProvider>
inline void CKKSEvalChebyshevPolynomialWithKeyProvider(
    CKKSPowerPolynomialResult<P, StartLogQ, LogDelta, CoeffLogDelta, Degree>
        &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const std::vector<double> &coeffs, const RelinKeyProvider &keys)
{
    using Traits = CKKSPowerPolynomialEvaluatorTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree>;
    typename Traits::PowerBasis powers;
    ckks_detail::CKKSBuildChebyshevBasisImpl<
        1, Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(powers, ct,
                                                                  keys);
    CKKSSetTransparentReal<P, Traits::log_q, LogDelta>(
        res, coeffs.empty() ? 0.0 : coeffs[0]);
    ckks_detail::CKKSAddPowerPolynomialTermsImpl<
        1, Traits, P, StartLogQ, LogDelta, CoeffLogDelta, Degree>(res, powers,
                                                                  coeffs);
}

namespace ckks_detail {

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t InvDegree,
          bool Enabled>
struct CKKSEvalModInverseCorrectionTraitsImpl;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t InvDegree>
struct CKKSEvalModInverseCorrectionTraitsImpl<
    P, StartLogQ, LogDelta, CoeffLogDelta, InvDegree, false> {
    static constexpr bool enabled = false;
    static constexpr std::uint32_t power_depth = 0;
    static constexpr std::uint32_t log_q = StartLogQ;
    static constexpr std::uint32_t log_delta = LogDelta;

    using RelinKeyChain = std::tuple<>;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t InvDegree>
struct CKKSEvalModInverseCorrectionTraitsImpl<
    P, StartLogQ, LogDelta, CoeffLogDelta, InvDegree, true> {
    static_assert(InvDegree > 1);
    using PolynomialTraits =
        CKKSPowerPolynomialEvaluatorTraits<P, StartLogQ, LogDelta,
                                           CoeffLogDelta, InvDegree>;
    static constexpr bool enabled = true;
    static constexpr std::uint32_t power_depth =
        PolynomialTraits::power_depth;
    static constexpr std::uint32_t log_q = PolynomialTraits::log_q;
    static constexpr std::uint32_t log_delta = LogDelta;

    using RelinKeyChain = typename PolynomialTraits::RelinKeyChain;
    using Ciphertext = typename PolynomialTraits::Ciphertext;
};

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t InvDegree>
using CKKSEvalModInverseCorrectionTraits =
    CKKSEvalModInverseCorrectionTraitsImpl<
        P, StartLogQ, LogDelta, CoeffLogDelta, InvDegree, InvDegree != 0>;

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
struct CKKSEvalModBoundedCosTraits {
    static_assert(DoubleAngle > 0);
    static_assert(InvDegree == 0 || InvDegree > 1);
    using PolynomialTraits =
        CKKSPowerPolynomialEvaluatorTraits<P, StartLogQ, LogDelta,
                                           CoeffLogDelta, Degree>;
    static constexpr std::uint32_t polynomial_log_q =
        PolynomialTraits::log_q;
    static_assert(polynomial_log_q > DoubleAngle * LogDelta);
    static constexpr std::uint32_t after_double_angle_log_q =
        polynomial_log_q - DoubleAngle * LogDelta;
    using InverseTraits = ckks_detail::CKKSEvalModInverseCorrectionTraits<
        P, after_double_angle_log_q, LogDelta, CoeffLogDelta, InvDegree>;
    static constexpr bool inverse_enabled = InverseTraits::enabled;
    static constexpr std::uint32_t log_q = InverseTraits::log_q;
    static constexpr std::uint32_t log_delta = LogDelta;

    using PolynomialCiphertext = typename PolynomialTraits::Ciphertext;
    using DoubleAngleCiphertext =
        CKKSCiphertext<P, after_double_angle_log_q, log_delta>;
    using PolynomialRelinKeyChain = typename PolynomialTraits::RelinKeyChain;
    using DoubleAngleRelinKeyChain =
        CKKSRelinKeyChain<P, polynomial_log_q, LogDelta, DoubleAngle>;
    using InverseRelinKeyChain = typename InverseTraits::RelinKeyChain;
    using Ciphertext = typename InverseTraits::Ciphertext;
};

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
using CKKSEvalModBoundedCosResult =
    typename CKKSEvalModBoundedCosTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree,
        DoubleAngle, InvDegree>::Ciphertext;

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
struct CKKSEvalModBoundedCosRelinKeys {
    using Traits = CKKSEvalModBoundedCosTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree, DoubleAngle,
        InvDegree>;

    typename Traits::PolynomialRelinKeyChain polynomial{};
    typename Traits::DoubleAngleRelinKeyChain double_angle{};
    typename Traits::InverseRelinKeyChain inverse{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        if constexpr (Traits::inverse_enabled)
            archive(polynomial, double_angle, inverse);
        else
            archive(polynomial, double_angle);
    }
};

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
inline void CKKSEvalModBoundedCosKeyGen(
    CKKSEvalModBoundedCosRelinKeys<P, StartLogQ, LogDelta, CoeffLogDelta,
                                   Degree, DoubleAngle, InvDegree> &keys,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    using Traits = CKKSEvalModBoundedCosTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree, DoubleAngle,
        InvDegree>;
    CKKSRelinKeyChainGen<P, StartLogQ, LogDelta,
                         Traits::PolynomialTraits::relin_depth>(
        keys.polynomial, key, noise);
    CKKSRelinKeyChainGen<P, Traits::polynomial_log_q, LogDelta, DoubleAngle>(
        keys.double_angle, key, noise);
    if constexpr (Traits::inverse_enabled) {
        CKKSRelinKeyChainGen<P, Traits::after_double_angle_log_q, LogDelta,
                             Traits::InverseTraits::power_depth>(
            keys.inverse, key, noise);
    }
}

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::uint32_t DoubleAngle,
          class RelinKeyProvider>
inline void CKKSBoundedCosDoubleAngleImpl(
    CKKSCiphertext<P, StartLogQ - DoubleAngle * LogDelta, LogDelta> &res,
    const CKKSCiphertext<P, StartLogQ - I * LogDelta, LogDelta> &ct,
    const RelinKeyProvider &keys, double sqrt_coeff)
{
    if constexpr (I == DoubleAngle) {
        res = ct;
    }
    else {
        constexpr std::uint32_t cur_log_q = StartLogQ - I * LogDelta;
        using NextCt = CKKSMultResult<P, cur_log_q, LogDelta, cur_log_q,
                                      LogDelta>;
        auto next = std::make_unique<NextCt>();
        CKKSSquare<P>(*next, ct, keys.template get<I>());
        maybe_release_key<I>(keys);
        CKKSMulIntegerInPlace<P, NextCt::log_q, LogDelta>(*next, 2);

        const double squared_coeff = sqrt_coeff * sqrt_coeff;
        CKKSSubPlainRealInPlace<P, NextCt::log_q, LogDelta>(*next,
                                                            squared_coeff);
        CKKSBoundedCosDoubleAngleImpl<I + 1, P, StartLogQ, LogDelta,
                                      DoubleAngle>(res, *next, keys,
                                                   squared_coeff);
    }
}

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
struct CKKSEvalModBoundedCosInMemoryKeyProvider {
    using RelinKeys =
        CKKSEvalModBoundedCosRelinKeys<P, StartLogQ, LogDelta, CoeffLogDelta,
                                       Degree, DoubleAngle, InvDegree>;
    const RelinKeys *keys = nullptr;

    explicit CKKSEvalModBoundedCosInMemoryKeyProvider(const RelinKeys &keys_)
        : keys(&keys_)
    {}

    template <std::size_t I>
    const auto &polynomial_relin() const
    {
        assert(keys != nullptr);
        return keys->polynomial.template get<I>();
    }

    template <std::size_t I>
    const auto &double_angle_relin() const
    {
        assert(keys != nullptr);
        return keys->double_angle.template get<I>();
    }

    template <std::size_t I>
    const auto &inverse_relin() const
    {
        assert(keys != nullptr);
        return keys->inverse.template get<I>();
    }
};

namespace ckks_detail {

template <class KeyProvider, bool Polynomial>
struct CKKSEvalModBoundedCosRelinKeyProviderChain {
    const KeyProvider &provider;

    template <std::size_t I>
    decltype(auto) get() const
    {
        if constexpr (Polynomial)
            return provider.template polynomial_relin<I>();
        else
            return provider.template double_angle_relin<I>();
    }

    template <std::size_t I>
    void release() const
    {
        if constexpr (Polynomial) {
            if constexpr (requires {
                              provider.template release_polynomial_relin<I>();
                          }) {
                provider.template release_polynomial_relin<I>();
            }
        }
        else {
            if constexpr (requires {
                              provider.template release_double_angle_relin<I>();
                          }) {
                provider.template release_double_angle_relin<I>();
            }
        }
    }
};

template <class KeyProvider>
struct CKKSEvalModBoundedCosInverseRelinKeyProviderChain {
    const KeyProvider &provider;

    template <std::size_t I>
    decltype(auto) get() const
    {
        return provider.template inverse_relin<I>();
    }

    template <std::size_t I>
    void release() const
    {
        if constexpr (requires { provider.template release_inverse_relin<I>(); }) {
            provider.template release_inverse_relin<I>();
        }
    }
};

}  // namespace ckks_detail

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, class KeyProvider,
          std::size_t InvDegree = 0>
inline void CKKSEvalModBoundedCosNormalizedWithKeyProvider(
    CKKSEvalModBoundedCosResult<P, StartLogQ, LogDelta, CoeffLogDelta, Degree,
                                DoubleAngle, InvDegree> &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const CKKSBoundedCosEvalModPolynomial &poly,
    const KeyProvider &key_provider)
{
    assert(poly.chebyshev_coeffs.size() <= Degree + 1);
    assert(poly.double_angle == DoubleAngle);

    using Traits = CKKSEvalModBoundedCosTraits<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree, DoubleAngle,
        InvDegree>;
    const bool trace_evalmod = ckks_detail::ckks_evalmod_trace_enabled();
    std::chrono::steady_clock::time_point trace_total_start{};
    std::chrono::steady_clock::time_point trace_section_start{};
    if (trace_evalmod) {
        trace_total_start = std::chrono::steady_clock::now();
        trace_section_start = trace_total_start;
    }
    auto shifted =
        std::make_unique<CKKSCiphertext<P, StartLogQ, LogDelta>>(ct);
    CKKSSubPlainRealInPlace<P, StartLogQ, LogDelta>(*shifted,
                                                    poly.domain_offset);

    auto polynomial = std::make_unique<typename Traits::PolynomialCiphertext>();
    const ckks_detail::CKKSEvalModBoundedCosRelinKeyProviderChain<KeyProvider,
                                                                  true>
        polynomial_keys{key_provider};
    if constexpr (Traits::PolynomialTraits::bsgs_enabled) {
        assert(poly.power_coeffs.size() <= Degree + 1);
        CKKSEvalPowerPolynomialWithKeyProvider<P, StartLogQ, LogDelta,
                                               CoeffLogDelta, Degree>(
            *polynomial, *shifted, poly.power_coeffs, polynomial_keys);
    }
    else {
        CKKSEvalChebyshevPolynomialWithKeyProvider<P, StartLogQ, LogDelta,
                                                   CoeffLogDelta, Degree>(
            *polynomial, *shifted, poly.chebyshev_coeffs, polynomial_keys);
    }
    shifted.reset();
    if (trace_evalmod) {
        ckks_detail::ckks_trace_evalmod_event(
            Traits::PolynomialTraits::bsgs_enabled
                ? "evalmod_power_bsgs_polynomial"
                : "evalmod_chebyshev_polynomial",
            ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
        trace_section_start = std::chrono::steady_clock::now();
    }

    const ckks_detail::CKKSEvalModBoundedCosRelinKeyProviderChain<KeyProvider,
                                                                  false>
        double_angle_keys{key_provider};
    if constexpr (Traits::inverse_enabled) {
        auto double_angle =
            std::make_unique<typename Traits::DoubleAngleCiphertext>();
        ckks_detail::CKKSBoundedCosDoubleAngleImpl<
            0, P, Traits::polynomial_log_q, LogDelta, DoubleAngle>(
            *double_angle, *polynomial, double_angle_keys, poly.sqrt_coeff);
        polynomial.reset();
        if (trace_evalmod) {
            ckks_detail::ckks_trace_evalmod_event(
                "evalmod_double_angle",
                ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
            trace_section_start = std::chrono::steady_clock::now();
        }

        const ckks_detail::CKKSEvalModBoundedCosInverseRelinKeyProviderChain<
            KeyProvider>
            inverse_keys{key_provider};
        CKKSEvalPowerPolynomialWithKeyProvider<
            P, Traits::after_double_angle_log_q, LogDelta, CoeffLogDelta,
            InvDegree>(res, *double_angle,
                       CKKSBuildEvalModInversePowerCoefficients(poly, InvDegree),
                       inverse_keys);
        if (trace_evalmod) {
            ckks_detail::ckks_trace_evalmod_event(
                "evalmod_inverse_correction",
                ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
            ckks_detail::ckks_trace_evalmod_event(
                "evalmod_total",
                ckks_detail::ckks_trace_elapsed_ms(trace_total_start));
        }
    }
    else {
        ckks_detail::CKKSBoundedCosDoubleAngleImpl<
            0, P, Traits::polynomial_log_q, LogDelta, DoubleAngle>(
            res, *polynomial, double_angle_keys, poly.sqrt_coeff);
        if (trace_evalmod) {
            ckks_detail::ckks_trace_evalmod_event(
                "evalmod_double_angle",
                ckks_detail::ckks_trace_elapsed_ms(trace_section_start));
            ckks_detail::ckks_trace_evalmod_event(
                "evalmod_total",
                ckks_detail::ckks_trace_elapsed_ms(trace_total_start));
        }
    }
}

template <class P, std::uint32_t StartLogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta, std::size_t Degree,
          std::uint32_t DoubleAngle, std::size_t InvDegree = 0>
inline void CKKSEvalModBoundedCosNormalized(
    CKKSEvalModBoundedCosResult<P, StartLogQ, LogDelta, CoeffLogDelta, Degree,
                                DoubleAngle, InvDegree> &res,
    const CKKSCiphertext<P, StartLogQ, LogDelta> &ct,
    const CKKSBoundedCosEvalModPolynomial &poly,
    const CKKSEvalModBoundedCosRelinKeys<P, StartLogQ, LogDelta,
                                         CoeffLogDelta, Degree, DoubleAngle,
                                         InvDegree>
        &keys)
{
    const CKKSEvalModBoundedCosInMemoryKeyProvider<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree, DoubleAngle, InvDegree>
        key_provider(keys);
    CKKSEvalModBoundedCosNormalizedWithKeyProvider<
        P, StartLogQ, LogDelta, CoeffLogDelta, Degree, DoubleAngle,
        decltype(key_provider), InvDegree>(
        res, ct, poly, key_provider);
}

template <class Schedule>
using CKKSDenseEvalModBoundedCosTraits = CKKSEvalModBoundedCosTraits<
    typename Schedule::Param, Schedule::after_component_split_log_q,
    Schedule::log_delta, Schedule::evalmod_log_scale, Schedule::evalmod_degree,
    Schedule::evalmod_double_angle, Schedule::evalmod_inv_degree>;

template <class Schedule>
using CKKSDenseEvalModBoundedCosResult =
    typename CKKSDenseEvalModBoundedCosTraits<Schedule>::Ciphertext;

template <class Schedule>
using CKKSDenseEvalModBoundedCosRelinKeys =
    CKKSEvalModBoundedCosRelinKeys<
        typename Schedule::Param, Schedule::after_component_split_log_q,
        Schedule::log_delta, Schedule::evalmod_log_scale,
        Schedule::evalmod_degree, Schedule::evalmod_double_angle,
        Schedule::evalmod_inv_degree>;

template <class Schedule>
using CKKSDenseEvalModBoundedCosInMemoryKeyProvider =
    CKKSEvalModBoundedCosInMemoryKeyProvider<
        typename Schedule::Param, Schedule::after_component_split_log_q,
        Schedule::log_delta, Schedule::evalmod_log_scale,
        Schedule::evalmod_degree, Schedule::evalmod_double_angle,
        Schedule::evalmod_inv_degree>;

template <class Schedule, std::size_t I>
using CKKSDenseEvalModPolynomialRelinKey =
    CKKSRelinKeyChainElement<typename Schedule::Param,
                             Schedule::after_component_split_log_q,
                             Schedule::log_delta, I>;

template <class Schedule, std::size_t I>
using CKKSDenseEvalModDoubleAngleRelinKey =
    CKKSRelinKeyChainElement<
        typename Schedule::Param,
        CKKSDenseEvalModBoundedCosTraits<Schedule>::polynomial_log_q,
        Schedule::log_delta, I>;

template <class Schedule, std::size_t I>
using CKKSDenseEvalModInverseRelinKey =
    CKKSRelinKeyChainElement<
        typename Schedule::Param,
        CKKSDenseEvalModBoundedCosTraits<Schedule>::after_double_angle_log_q,
        Schedule::log_delta, I>;

template <class Schedule, std::size_t I>
using CKKSDenseSeededEvalModPolynomialRelinKey =
    CKKSSeededRelinKey<typename Schedule::Param,
                       Schedule::after_component_split_log_q -
                           (I + 1) * Schedule::log_delta>;

template <class Schedule, std::size_t I>
using CKKSDenseSeededEvalModDoubleAngleRelinKey =
    CKKSSeededRelinKey<
        typename Schedule::Param,
        CKKSDenseEvalModBoundedCosTraits<Schedule>::polynomial_log_q -
            (I + 1) * Schedule::log_delta>;

template <class Schedule, std::size_t I>
using CKKSDenseSeededEvalModInverseRelinKey =
    CKKSSeededRelinKey<
        typename Schedule::Param,
        CKKSDenseEvalModBoundedCosTraits<Schedule>::after_double_angle_log_q -
            (I + 1) * Schedule::log_delta>;

namespace ckks_detail {

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
inline std::size_t CKKSRotationUsageKeySwitchRows(
    const std::array<CKKSRotationKeyIndexSet<P>, LevelCount> &usage)
{
    if constexpr (I == LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * LevelStep;
        return CKKSRotationKeyIndexSetCount<P>(usage[I]) *
                   CKKSAutoKeySwitchRowCount<P, log_q>() +
               CKKSRotationUsageKeySwitchRows<I + 1, P, StartLogQ, LevelStep>(
                   usage);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
inline std::size_t CKKSDirectRotationUsageKeySwitchRows(
    const std::array<CKKSDirectRotationKeyIndexSet<P>, LevelCount> &usage)
{
    if constexpr (I == LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * LevelStep;
        return CKKSDirectRotationKeyIndexSetCount<P>(usage[I]) *
                   CKKSAutoKeySwitchRowCount<P, log_q>() +
               CKKSDirectRotationUsageKeySwitchRows<I + 1, P, StartLogQ,
                                                    LevelStep>(usage);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
inline std::size_t CKKSRotationUsagePeakKeySwitchRows(
    const std::array<CKKSRotationKeyIndexSet<P>, LevelCount> &usage)
{
    if constexpr (I == LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * LevelStep;
        const std::size_t here =
            CKKSRotationKeyIndexSetCount<P>(usage[I]) *
            CKKSAutoKeySwitchRowCount<P, log_q>();
        const std::size_t tail =
            CKKSRotationUsagePeakKeySwitchRows<I + 1, P, StartLogQ,
                                               LevelStep>(usage);
        return std::max(here, tail);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
inline std::size_t CKKSDirectRotationUsagePeakKeySwitchRows(
    const std::array<CKKSDirectRotationKeyIndexSet<P>, LevelCount> &usage)
{
    if constexpr (I == LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * LevelStep;
        const std::size_t here =
            CKKSDirectRotationKeyIndexSetCount<P>(usage[I]) *
            CKKSAutoKeySwitchRowCount<P, log_q>();
        const std::size_t tail =
            CKKSDirectRotationUsagePeakKeySwitchRows<I + 1, P, StartLogQ,
                                                     LevelStep>(usage);
        return std::max(here, tail);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
inline std::size_t CKKSRotationUsageAdjacentPeakKeySwitchRows(
    const std::array<CKKSRotationKeyIndexSet<P>, LevelCount> &usage)
{
    if constexpr (I + 1 >= LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t input_log_q = StartLogQ - I * LevelStep;
        constexpr std::uint32_t output_log_q =
            StartLogQ - (I + 1) * LevelStep;
        const std::size_t here =
            CKKSRotationKeyIndexSetCount<P>(usage[I]) *
                CKKSAutoKeySwitchRowCount<P, input_log_q>() +
            CKKSRotationKeyIndexSetCount<P>(usage[I + 1]) *
                CKKSAutoKeySwitchRowCount<P, output_log_q>();
        const std::size_t tail =
            CKKSRotationUsageAdjacentPeakKeySwitchRows<I + 1, P, StartLogQ,
                                                       LevelStep>(usage);
        return std::max(here, tail);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
inline std::size_t CKKSDirectRotationUsageAdjacentPeakKeySwitchRows(
    const std::array<CKKSDirectRotationKeyIndexSet<P>, LevelCount> &usage)
{
    if constexpr (I + 1 >= LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t input_log_q = StartLogQ - I * LevelStep;
        constexpr std::uint32_t output_log_q =
            StartLogQ - (I + 1) * LevelStep;
        const std::size_t here =
            CKKSDirectRotationKeyIndexSetCount<P>(usage[I]) *
                CKKSAutoKeySwitchRowCount<P, input_log_q>() +
            CKKSDirectRotationKeyIndexSetCount<P>(usage[I + 1]) *
                CKKSAutoKeySwitchRowCount<P, output_log_q>();
        const std::size_t tail =
            CKKSDirectRotationUsageAdjacentPeakKeySwitchRows<
                I + 1, P, StartLogQ, LevelStep>(usage);
        return std::max(here, tail);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
inline std::size_t CKKSHybridRotationUsageAdjacentPeakKeySwitchRows(
    const std::array<CKKSRotationKeyIndexSet<P>, LevelCount> &binary_usage,
    const std::array<CKKSDirectRotationKeyIndexSet<P>, LevelCount>
        &direct_usage)
{
    if constexpr (I + 1 >= LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t input_log_q = StartLogQ - I * LevelStep;
        constexpr std::uint32_t output_log_q =
            StartLogQ - (I + 1) * LevelStep;
        const std::size_t input_count =
            CKKSRotationKeyIndexSetCount<P>(binary_usage[I]) +
            CKKSDirectRotationKeyIndexSetCount<P>(direct_usage[I]);
        const std::size_t output_count =
            CKKSRotationKeyIndexSetCount<P>(binary_usage[I + 1]) +
            CKKSDirectRotationKeyIndexSetCount<P>(direct_usage[I + 1]);
        const std::size_t here =
            input_count * CKKSAutoKeySwitchRowCount<P, input_log_q>() +
            output_count * CKKSAutoKeySwitchRowCount<P, output_log_q>();
        const std::size_t tail =
            CKKSHybridRotationUsageAdjacentPeakKeySwitchRows<
                I + 1, P, StartLogQ, LevelStep>(binary_usage, direct_usage);
        return std::max(here, tail);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
inline std::size_t CKKSHybridRotationUsagePeakKeySwitchRows(
    const std::array<CKKSRotationKeyIndexSet<P>, LevelCount> &binary_usage,
    const std::array<CKKSDirectRotationKeyIndexSet<P>, LevelCount>
        &direct_usage)
{
    if constexpr (I == LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * LevelStep;
        const std::size_t count =
            CKKSRotationKeyIndexSetCount<P>(binary_usage[I]) +
            CKKSDirectRotationKeyIndexSetCount<P>(direct_usage[I]);
        const std::size_t here = count * CKKSAutoKeySwitchRowCount<P, log_q>();
        const std::size_t tail =
            CKKSHybridRotationUsagePeakKeySwitchRows<
                I + 1, P, StartLogQ, LevelStep>(binary_usage, direct_usage);
        return std::max(here, tail);
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LevelStep, std::size_t LevelCount>
constexpr std::size_t CKKSFullGaloisKeySwitchRows()
{
    if constexpr (I == LevelCount) {
        return 0;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - I * LevelStep;
        return (P::nbit + 1) * CKKSAutoKeySwitchRowCount<P, log_q>() +
               CKKSFullGaloisKeySwitchRows<I + 1, P, StartLogQ, LevelStep,
                                           LevelCount>();
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::size_t Depth>
constexpr std::size_t CKKSRelinChainKeySwitchRows()
{
    if constexpr (I == Depth) {
        return 0;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - (I + 1) * LogDelta;
        return CKKSRelinKeySwitchRowCount<P, log_q>() +
               CKKSRelinChainKeySwitchRows<I + 1, P, StartLogQ, LogDelta,
                                           Depth>();
    }
}

template <std::size_t I, class P, std::uint32_t StartLogQ,
          std::uint32_t LogDelta, std::size_t Depth>
constexpr std::size_t CKKSRelinChainPeakKeySwitchRows()
{
    if constexpr (I == Depth) {
        return 0;
    }
    else {
        constexpr std::uint32_t log_q = StartLogQ - (I + 1) * LogDelta;
        return std::max(CKKSRelinKeySwitchRowCount<P, log_q>(),
                        CKKSRelinChainPeakKeySwitchRows<
                            I + 1, P, StartLogQ, LogDelta, Depth>());
    }
}

}  // namespace ckks_detail

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapEvalModRelinKeyCount()
{
    return CKKSDenseEvalModBoundedCosTraits<
               Schedule>::PolynomialTraits::relin_depth +
           Schedule::evalmod_double_angle +
           CKKSDenseEvalModBoundedCosTraits<
               Schedule>::InverseTraits::power_depth;
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapEvalModKeySwitchRowCount()
{
    using P = typename Schedule::Param;
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    return ckks_detail::CKKSRelinChainKeySwitchRows<
               0, P, Schedule::after_component_split_log_q,
               Schedule::log_delta,
               Traits::PolynomialTraits::relin_depth>() +
           ckks_detail::CKKSRelinChainKeySwitchRows<
               0, P, Traits::polynomial_log_q, Schedule::log_delta,
               Schedule::evalmod_double_angle>() +
           ckks_detail::CKKSRelinChainKeySwitchRows<
               0, P, Traits::after_double_angle_log_q, Schedule::log_delta,
               Traits::InverseTraits::power_depth>();
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapEvalModPeakKeySwitchRowCount()
{
    using P = typename Schedule::Param;
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    const std::size_t polynomial_peak =
        ckks_detail::CKKSRelinChainPeakKeySwitchRows<
            0, P, Schedule::after_component_split_log_q,
            Schedule::log_delta, Traits::PolynomialTraits::relin_depth>();
    const std::size_t double_angle_peak =
        ckks_detail::CKKSRelinChainPeakKeySwitchRows<
            0, P, Traits::polynomial_log_q, Schedule::log_delta,
            Schedule::evalmod_double_angle>();
    const std::size_t inverse_peak =
        ckks_detail::CKKSRelinChainPeakKeySwitchRows<
            0, P, Traits::after_double_angle_log_q, Schedule::log_delta,
            Traits::InverseTraits::power_depth>();
    return std::max(std::max(polynomial_peak, double_angle_peak), inverse_peak);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseCoeffToSlotKeySwitchRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSRotationUsageKeySwitchRows<
        0, P, Schedule::boot_log_q,
        Schedule::coeff_to_slot_plain_log_delta>(
        usage.coeff_to_slot);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectCoeffToSlotKeySwitchRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSDirectRotationUsageKeySwitchRows<
        0, P, Schedule::boot_log_q,
        Schedule::coeff_to_slot_plain_log_delta>(
        usage.coeff_to_slot);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantCoeffToSlotKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSRotationUsageKeySwitchRows<
               0, P, Schedule::boot_log_q,
               Schedule::coeff_to_slot_plain_log_delta>(
               usage.coeff_to_slot_binary) +
           ckks_detail::CKKSDirectRotationUsageKeySwitchRows<
               0, P, Schedule::boot_log_q,
               Schedule::coeff_to_slot_plain_log_delta>(
               usage.coeff_to_slot_direct);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseCoeffToSlotPeakKeySwitchRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSRotationUsageAdjacentPeakKeySwitchRows<
        0, P, Schedule::boot_log_q,
        Schedule::coeff_to_slot_plain_log_delta>(
        usage.coeff_to_slot);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantCoeffToSlotPeakKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSHybridRotationUsageAdjacentPeakKeySwitchRows<
        0, P, Schedule::boot_log_q,
        Schedule::coeff_to_slot_plain_log_delta>(
        usage.coeff_to_slot_binary, usage.coeff_to_slot_direct);
}

template <class Schedule>
inline std::size_t
CKKSDenseBootstrapHybridGiantCoeffToSlotLowCachePeakKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSHybridRotationUsagePeakKeySwitchRows<
        0, P, Schedule::boot_log_q,
        Schedule::coeff_to_slot_plain_log_delta>(
        usage.coeff_to_slot_binary, usage.coeff_to_slot_direct);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectCoeffToSlotPeakKeySwitchRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSDirectRotationUsageAdjacentPeakKeySwitchRows<
        0, P, Schedule::boot_log_q,
        Schedule::coeff_to_slot_plain_log_delta>(
        usage.coeff_to_slot);
}

template <class Schedule>
constexpr std::size_t
CKKSDenseBootstrapPackedConjugateKeySwitchRowCount()
{
    using P = typename Schedule::Param;
    return CKKSAutoKeySwitchRowCount<P,
                                     Schedule::after_coeff_to_slot_log_q>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparsePackedConjugateKeySwitchRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSRotationKeyIndexSetCount<P>(usage.packed_conjugate) *
           CKKSDenseBootstrapPackedConjugateKeySwitchRowCount<Schedule>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectPackedConjugateKeySwitchRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSRotationKeyIndexSetCount<P>(usage.packed_conjugate) *
           CKKSDenseBootstrapPackedConjugateKeySwitchRowCount<Schedule>();
}

template <class Schedule>
inline std::size_t
CKKSDenseBootstrapHybridGiantPackedConjugateKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSRotationKeyIndexSetCount<P>(usage.packed_conjugate) *
           CKKSDenseBootstrapPackedConjugateKeySwitchRowCount<Schedule>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseSlotToCoeffKeySwitchRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSRotationUsageKeySwitchRows<
        0, P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta>(usage.slot_to_coeff);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantSlotToCoeffKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSRotationUsageKeySwitchRows<
               0, P, Schedule::after_evalmod_log_q,
               Schedule::slot_to_coeff_plain_log_delta>(
               usage.slot_to_coeff_binary) +
           ckks_detail::CKKSDirectRotationUsageKeySwitchRows<
               0, P, Schedule::after_evalmod_log_q,
               Schedule::slot_to_coeff_plain_log_delta>(
               usage.slot_to_coeff_direct);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectSlotToCoeffKeySwitchRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSDirectRotationUsageKeySwitchRows<
        0, P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta>(usage.slot_to_coeff);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseSlotToCoeffPeakKeySwitchRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSRotationUsageAdjacentPeakKeySwitchRows<
        0, P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta>(usage.slot_to_coeff);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantSlotToCoeffPeakKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSHybridRotationUsageAdjacentPeakKeySwitchRows<
        0, P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta>(
        usage.slot_to_coeff_binary, usage.slot_to_coeff_direct);
}

template <class Schedule>
inline std::size_t
CKKSDenseBootstrapHybridGiantSlotToCoeffLowCachePeakKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSHybridRotationUsagePeakKeySwitchRows<
        0, P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta>(
        usage.slot_to_coeff_binary, usage.slot_to_coeff_direct);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectSlotToCoeffPeakKeySwitchRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSDirectRotationUsageAdjacentPeakKeySwitchRows<
        0, P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta>(usage.slot_to_coeff);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseGaloisKeySwitchRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    return CKKSDenseBootstrapSparseCoeffToSlotKeySwitchRowCount<Schedule>(
               usage) +
           CKKSDenseBootstrapSparsePackedConjugateKeySwitchRowCount<Schedule>(
               usage) +
           CKKSDenseBootstrapSparseSlotToCoeffKeySwitchRowCount<Schedule>(
           usage);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectGaloisKeySwitchRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    return CKKSDenseBootstrapDirectCoeffToSlotKeySwitchRowCount<Schedule>(
               usage) +
           CKKSDenseBootstrapDirectPackedConjugateKeySwitchRowCount<Schedule>(
               usage) +
           CKKSDenseBootstrapDirectSlotToCoeffKeySwitchRowCount<Schedule>(
           usage);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantGaloisKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    return CKKSDenseBootstrapHybridGiantCoeffToSlotKeySwitchRowCount<
               Schedule>(usage) +
           CKKSDenseBootstrapHybridGiantPackedConjugateKeySwitchRowCount<
               Schedule>(usage) +
           CKKSDenseBootstrapHybridGiantSlotToCoeffKeySwitchRowCount<
               Schedule>(usage);
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseGaloisPeakKeySwitchRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    return std::max(
        {CKKSDenseBootstrapSparseCoeffToSlotPeakKeySwitchRowCount<Schedule>(
             usage),
         CKKSDenseBootstrapSparsePackedConjugateKeySwitchRowCount<Schedule>(
             usage),
         CKKSDenseBootstrapSparseSlotToCoeffPeakKeySwitchRowCount<Schedule>(
             usage)});
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantGaloisPeakKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    return std::max(
        {CKKSDenseBootstrapHybridGiantCoeffToSlotPeakKeySwitchRowCount<
             Schedule>(usage),
         CKKSDenseBootstrapHybridGiantPackedConjugateKeySwitchRowCount<
             Schedule>(usage),
         CKKSDenseBootstrapHybridGiantSlotToCoeffPeakKeySwitchRowCount<
             Schedule>(usage)});
}

template <class Schedule>
inline std::size_t
CKKSDenseBootstrapHybridGiantGaloisLowCachePeakKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    return std::max(
        {CKKSDenseBootstrapHybridGiantCoeffToSlotLowCachePeakKeySwitchRowCount<
             Schedule>(usage),
         CKKSDenseBootstrapHybridGiantPackedConjugateKeySwitchRowCount<
             Schedule>(usage),
         CKKSDenseBootstrapHybridGiantSlotToCoeffLowCachePeakKeySwitchRowCount<
             Schedule>(usage)});
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectGaloisPeakKeySwitchRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    return std::max(
        {CKKSDenseBootstrapDirectCoeffToSlotPeakKeySwitchRowCount<Schedule>(
             usage),
         CKKSDenseBootstrapDirectPackedConjugateKeySwitchRowCount<Schedule>(
             usage),
         CKKSDenseBootstrapDirectSlotToCoeffPeakKeySwitchRowCount<Schedule>(
             usage)});
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapFullGaloisKeySwitchRowCount()
{
    using P = typename Schedule::Param;
    return ckks_detail::CKKSFullGaloisKeySwitchRows<
               0, P, Schedule::boot_log_q,
               Schedule::coeff_to_slot_plain_log_delta,
               Schedule::coeff_to_slot_level_count + 1>() +
           (P::nbit + 1) *
               CKKSAutoKeySwitchRowCount<
                   P, Schedule::after_coeff_to_slot_log_q>() +
           ckks_detail::CKKSFullGaloisKeySwitchRows<
               0, P, Schedule::after_evalmod_log_q,
               Schedule::slot_to_coeff_plain_log_delta,
               Schedule::slot_to_coeff_level_count + 1>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseKeySwitchRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    return CKKSDenseBootstrapSparseGaloisKeySwitchRowCount<Schedule>(usage) +
           CKKSDenseBootstrapEvalModKeySwitchRowCount<Schedule>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectKeySwitchRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    return CKKSDenseBootstrapDirectGaloisKeySwitchRowCount<Schedule>(usage) +
           CKKSDenseBootstrapEvalModKeySwitchRowCount<Schedule>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantKeySwitchRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    return CKKSDenseBootstrapHybridGiantGaloisKeySwitchRowCount<Schedule>(
               usage) +
           CKKSDenseBootstrapEvalModKeySwitchRowCount<Schedule>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapStreamedKeySwitchPeakRowCount(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    return std::max(
        CKKSDenseBootstrapSparseGaloisPeakKeySwitchRowCount<Schedule>(usage),
        CKKSDenseBootstrapEvalModPeakKeySwitchRowCount<Schedule>());
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantStreamedKeySwitchPeakRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    return std::max(
        CKKSDenseBootstrapHybridGiantGaloisPeakKeySwitchRowCount<Schedule>(
            usage),
        CKKSDenseBootstrapEvalModPeakKeySwitchRowCount<Schedule>());
}

template <class Schedule>
inline std::size_t
CKKSDenseBootstrapHybridGiantStreamedLowCacheKeySwitchPeakRowCount(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    return std::max(
        CKKSDenseBootstrapHybridGiantGaloisLowCachePeakKeySwitchRowCount<
            Schedule>(usage),
        CKKSDenseBootstrapEvalModPeakKeySwitchRowCount<Schedule>());
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectStreamedKeySwitchPeakRowCount(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    return std::max(
        CKKSDenseBootstrapDirectGaloisPeakKeySwitchRowCount<Schedule>(usage),
        CKKSDenseBootstrapEvalModPeakKeySwitchRowCount<Schedule>());
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapFullKeySwitchRowCount()
{
    return CKKSDenseBootstrapFullGaloisKeySwitchRowCount<Schedule>() +
           CKKSDenseBootstrapEvalModKeySwitchRowCount<Schedule>();
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapEncapsulationKeySwitchRowCount()
{
    using P = typename Schedule::Param;
    return CKKSAutoKeySwitchRowCount<P, Schedule::input_log_q>() +
           CKKSAutoKeySwitchRowCount<P, Schedule::output_log_q>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseKeyByteEstimate(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapSparseKeySwitchRowCount<Schedule>(usage) *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapSparseSeededKeyByteEstimate(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapSparseKeySwitchRowCount<Schedule>(usage) *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectKeyByteEstimate(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapDirectKeySwitchRowCount<Schedule>(usage) *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectSeededKeyByteEstimate(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapDirectKeySwitchRowCount<Schedule>(usage) *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantKeyByteEstimate(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapHybridGiantKeySwitchRowCount<Schedule>(usage) *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantSeededKeyByteEstimate(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapHybridGiantKeySwitchRowCount<Schedule>(usage) *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapStreamedPeakKeyByteEstimate(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapStreamedKeySwitchPeakRowCount<Schedule>(usage) *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapStreamedPeakSeededKeyByteEstimate(
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapStreamedKeySwitchPeakRowCount<Schedule>(usage) *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapHybridGiantStreamedPeakKeyByteEstimate(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapHybridGiantStreamedKeySwitchPeakRowCount<
               Schedule>(usage) *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t
CKKSDenseBootstrapHybridGiantStreamedPeakSeededKeyByteEstimate(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapHybridGiantStreamedKeySwitchPeakRowCount<
               Schedule>(usage) *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t
CKKSDenseBootstrapHybridGiantStreamedLowCachePeakKeyByteEstimate(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapHybridGiantStreamedLowCacheKeySwitchPeakRowCount<
               Schedule>(usage) *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t
CKKSDenseBootstrapHybridGiantStreamedLowCachePeakSeededKeyByteEstimate(
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapHybridGiantStreamedLowCacheKeySwitchPeakRowCount<
               Schedule>(usage) *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectStreamedPeakKeyByteEstimate(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapDirectStreamedKeySwitchPeakRowCount<Schedule>(
               usage) *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline std::size_t CKKSDenseBootstrapDirectStreamedPeakSeededKeyByteEstimate(
    const CKKSDenseBootstrapDirectRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapDirectStreamedKeySwitchPeakRowCount<Schedule>(
               usage) *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapFullKeyByteEstimate()
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapFullKeySwitchRowCount<Schedule>() *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapFullSeededKeyByteEstimate()
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapFullKeySwitchRowCount<Schedule>() *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapEncapsulationKeyByteEstimate()
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapEncapsulationKeySwitchRowCount<Schedule>() *
           CKKSKeySwitchRowByteSize<P>();
}

template <class Schedule>
constexpr std::size_t CKKSDenseBootstrapEncapsulationSeededKeyByteEstimate()
{
    using P = typename Schedule::Param;
    return CKKSDenseBootstrapEncapsulationKeySwitchRowCount<Schedule>() *
           CKKSSeededKeySwitchRowByteSize<P>();
}

template <class Schedule>
inline void CKKSDenseEvalModBoundedCosKeyGen(
    CKKSDenseEvalModBoundedCosRelinKeys<Schedule> &keys,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    static_assert(CKKSDenseEvalModBoundedCosTraits<Schedule>::log_q ==
                  Schedule::after_evalmod_log_q);
    CKKSEvalModBoundedCosKeyGen<
        typename Schedule::Param, Schedule::after_component_split_log_q,
        Schedule::log_delta, Schedule::evalmod_log_scale,
        Schedule::evalmod_degree, Schedule::evalmod_double_angle,
        Schedule::evalmod_inv_degree>(keys, key, noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseEvalModPolynomialRelinKeyGen(
    CKKSDenseEvalModPolynomialRelinKey<Schedule, I> &relinkey,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    CKKSRelinKeyChainElementGen<
        I, typename Schedule::Param, Schedule::after_component_split_log_q,
        Schedule::log_delta, Traits::PolynomialTraits::relin_depth>(relinkey, key,
                                                                    noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseEvalModDoubleAngleRelinKeyGen(
    CKKSDenseEvalModDoubleAngleRelinKey<Schedule, I> &relinkey,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    CKKSRelinKeyChainElementGen<I, typename Schedule::Param,
                                Traits::polynomial_log_q, Schedule::log_delta,
                                Schedule::evalmod_double_angle>(relinkey, key,
                                                                noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseEvalModInverseRelinKeyGen(
    CKKSDenseEvalModInverseRelinKey<Schedule, I> &relinkey,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    CKKSRelinKeyChainElementGen<I, typename Schedule::Param,
                                Traits::after_double_angle_log_q,
                                Schedule::log_delta,
                                Traits::InverseTraits::power_depth>(relinkey,
                                                                        key, noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseSeededEvalModPolynomialRelinKeyGen(
    CKKSDenseSeededEvalModPolynomialRelinKey<Schedule, I> &relinkey,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    CKKSSeededRelinKeyChainElementGen<
        I, typename Schedule::Param, Schedule::after_component_split_log_q,
        Schedule::log_delta, Traits::PolynomialTraits::relin_depth>(relinkey, key,
                                                                    noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseSeededEvalModDoubleAngleRelinKeyGen(
    CKKSDenseSeededEvalModDoubleAngleRelinKey<Schedule, I> &relinkey,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    CKKSSeededRelinKeyChainElementGen<
        I, typename Schedule::Param, Traits::polynomial_log_q,
        Schedule::log_delta, Schedule::evalmod_double_angle>(relinkey, key,
                                                             noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseSeededEvalModInverseRelinKeyGen(
    CKKSDenseSeededEvalModInverseRelinKey<Schedule, I> &relinkey,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    CKKSSeededRelinKeyChainElementGen<
        I, typename Schedule::Param, Traits::after_double_angle_log_q,
        Schedule::log_delta, Traits::InverseTraits::power_depth>(relinkey, key,
                                                                 noise);
}

template <class Schedule, class KeyProvider>
inline void CKKSDenseEvalModBoundedCosNormalizedWithKeyProvider(
    CKKSDenseEvalModBoundedCosResult<Schedule> &res,
    const typename Schedule::ComponentCiphertext &ct,
    const CKKSBoundedCosEvalModPolynomial &poly,
    const KeyProvider &key_provider);

template <class Schedule>
inline void CKKSDenseEvalModBoundedCosNormalized(
    CKKSDenseEvalModBoundedCosResult<Schedule> &res,
    const typename Schedule::ComponentCiphertext &ct,
    const CKKSBoundedCosEvalModPolynomial &poly,
    const CKKSDenseEvalModBoundedCosRelinKeys<Schedule> &keys)
{
    const CKKSDenseEvalModBoundedCosInMemoryKeyProvider<Schedule> key_provider(
        keys);
    CKKSDenseEvalModBoundedCosNormalizedWithKeyProvider<Schedule>(
        res, ct, poly, key_provider);
}

template <class Schedule, class KeyProvider>
inline void CKKSDenseEvalModBoundedCosNormalizedWithKeyProvider(
    CKKSDenseEvalModBoundedCosResult<Schedule> &res,
    const typename Schedule::ComponentCiphertext &ct,
    const CKKSBoundedCosEvalModPolynomial &poly,
    const KeyProvider &key_provider)
{
    static_assert(CKKSDenseEvalModBoundedCosTraits<Schedule>::log_q ==
                  Schedule::after_evalmod_log_q);
    CKKSEvalModBoundedCosNormalizedWithKeyProvider<
        typename Schedule::Param, Schedule::after_component_split_log_q,
        Schedule::log_delta, Schedule::evalmod_log_scale,
        Schedule::evalmod_degree, Schedule::evalmod_double_angle, KeyProvider,
        Schedule::evalmod_inv_degree>(res, ct, poly, key_provider);
}

template <class Schedule>
struct CKKSDenseBootstrapKey {
    using P = typename Schedule::Param;
    using CoeffToSlotGaloisKeyChain =
        CKKSSparseGaloisKeyChain<P, Schedule::boot_log_q,
                                 Schedule::coeff_to_slot_plain_log_delta,
                                 Schedule::coeff_to_slot_level_count>;
    using SlotToCoeffGaloisKeyChain =
        CKKSSparseGaloisKeyChain<P, Schedule::after_evalmod_log_q,
                                 Schedule::slot_to_coeff_plain_log_delta,
                                 Schedule::slot_to_coeff_level_count>;

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan{};
    CKKSBoundedCosEvalModPolynomial evalmod_polynomial{};
    CoeffToSlotGaloisKeyChain coeff_to_slot_galois{};
    CKKSSparseGaloisKey<P, Schedule::after_coeff_to_slot_log_q>
        packed_conjugate_galois{};
    CKKSDenseEvalModBoundedCosRelinKeys<Schedule> evalmod_relin{};
    SlotToCoeffGaloisKeyChain slot_to_coeff_galois{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(linear_plan, evalmod_polynomial, coeff_to_slot_galois,
                packed_conjugate_galois, evalmod_relin,
                slot_to_coeff_galois);
    }
};

template <class Schedule>
struct CKKSDenseBootstrapHybridGiantKey {
    using P = typename Schedule::Param;
    using CoeffToSlotGaloisKeyChain =
        CKKSHybridSparseGaloisKeyChain<P, Schedule::boot_log_q,
                                       Schedule::coeff_to_slot_plain_log_delta,
                                       Schedule::coeff_to_slot_level_count>;
    using SlotToCoeffGaloisKeyChain =
        CKKSHybridSparseGaloisKeyChain<P, Schedule::after_evalmod_log_q,
                                       Schedule::slot_to_coeff_plain_log_delta,
                                       Schedule::slot_to_coeff_level_count>;

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan{};
    CKKSBoundedCosEvalModPolynomial evalmod_polynomial{};
    CoeffToSlotGaloisKeyChain coeff_to_slot_galois{};
    CKKSSparseGaloisKey<P, Schedule::after_coeff_to_slot_log_q>
        packed_conjugate_galois{};
    CKKSDenseEvalModBoundedCosRelinKeys<Schedule> evalmod_relin{};
    SlotToCoeffGaloisKeyChain slot_to_coeff_galois{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(linear_plan, evalmod_polynomial, coeff_to_slot_galois,
                packed_conjugate_galois, evalmod_relin,
                slot_to_coeff_galois);
    }
};

template <class Schedule>
struct CKKSDenseBootstrapEncapsulationKey {
    using P = typename Schedule::Param;

    CKKSSecretKeySwitchKey<P, Schedule::input_log_q> input_to_bootstrap{};
    CKKSSecretKeySwitchKey<P, Schedule::output_log_q> bootstrap_to_output{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(input_to_bootstrap, bootstrap_to_output);
    }
};

template <class Schedule>
struct CKKSDenseBootstrapSeededEncapsulationKey {
    using P = typename Schedule::Param;

    CKKSSeededSecretKeySwitchKey<P, Schedule::input_log_q>
        input_to_bootstrap{};
    CKKSSeededSecretKeySwitchKey<P, Schedule::output_log_q>
        bootstrap_to_output{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(input_to_bootstrap, bootstrap_to_output);
    }
};

template <class Schedule>
inline void CKKSDenseBootstrapEncapsulationKeyGen(
    CKKSDenseBootstrapEncapsulationKey<Schedule> &encapsulation_key,
    const Key<typename Schedule::Param> &input_output_key,
    const Key<typename Schedule::Param> &bootstrap_key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    CKKSSecretKeySwitchKeyGen<P, Schedule::input_log_q>(
        encapsulation_key.input_to_bootstrap, input_output_key, bootstrap_key,
        noise);
    CKKSSecretKeySwitchKeyGen<P, Schedule::output_log_q>(
        encapsulation_key.bootstrap_to_output, bootstrap_key, input_output_key,
        noise);
}

template <class Schedule>
inline void CKKSDenseBootstrapSeededEncapsulationKeyGen(
    CKKSDenseBootstrapSeededEncapsulationKey<Schedule> &encapsulation_key,
    const Key<typename Schedule::Param> &input_output_key,
    const Key<typename Schedule::Param> &bootstrap_key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    CKKSSecretKeySwitchKeyGen<P, Schedule::input_log_q>(
        encapsulation_key.input_to_bootstrap, input_output_key, bootstrap_key,
        noise);
    CKKSSecretKeySwitchKeyGen<P, Schedule::output_log_q>(
        encapsulation_key.bootstrap_to_output, bootstrap_key, input_output_key,
        noise);
}

template <class Schedule, std::size_t I>
using CKKSDenseBootstrapCoeffToSlotGaloisKey = CKKSSparseGaloisKey<
    typename Schedule::Param,
    Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta>;

template <class Schedule, std::size_t I>
using CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKey =
    CKKSHybridSparseGaloisKey<
        typename Schedule::Param,
        Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta>;

template <class Schedule, std::size_t I>
using CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKey =
    CKKSSeededHybridSparseGaloisKey<
        typename Schedule::Param,
        Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta>;

template <class Schedule>
using CKKSDenseBootstrapPackedConjugateGaloisKey =
    CKKSSparseGaloisKey<typename Schedule::Param,
                        Schedule::after_coeff_to_slot_log_q>;

template <class Schedule>
using CKKSDenseBootstrapSeededPackedConjugateGaloisKey =
    CKKSSeededSparseGaloisKey<typename Schedule::Param,
                              Schedule::after_coeff_to_slot_log_q>;

template <class Schedule, std::size_t I>
using CKKSDenseBootstrapSlotToCoeffGaloisKey = CKKSSparseGaloisKey<
    typename Schedule::Param,
    Schedule::after_evalmod_log_q -
        I * Schedule::slot_to_coeff_plain_log_delta>;

template <class Schedule, std::size_t I>
using CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKey =
    CKKSHybridSparseGaloisKey<
        typename Schedule::Param,
        Schedule::after_evalmod_log_q -
            I * Schedule::slot_to_coeff_plain_log_delta>;

template <class Schedule, std::size_t I>
using CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKey =
    CKKSSeededHybridSparseGaloisKey<
        typename Schedule::Param,
        Schedule::after_evalmod_log_q -
            I * Schedule::slot_to_coeff_plain_log_delta>;

template <class Schedule, std::size_t I>
inline void CKKSDenseBootstrapCoeffToSlotGaloisKeyGen(
    CKKSDenseBootstrapCoeffToSlotGaloisKey<Schedule, I> &gk,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    static_assert(I <= Schedule::coeff_to_slot_level_count);
    constexpr std::uint32_t log_q =
        Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta;
    CKKSSparseGaloisKeyGen<P, log_q>(gk, key, usage.coeff_to_slot[I],
                                     noise);
}

template <class Schedule>
inline void CKKSDenseBootstrapPackedConjugateGaloisKeyGen(
    CKKSDenseBootstrapPackedConjugateGaloisKey<Schedule> &gk,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
        gk, key, usage.packed_conjugate, noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseBootstrapSlotToCoeffGaloisKeyGen(
    CKKSDenseBootstrapSlotToCoeffGaloisKey<Schedule, I> &gk,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    static_assert(I <= Schedule::slot_to_coeff_level_count);
    constexpr std::uint32_t log_q =
        Schedule::after_evalmod_log_q -
        I * Schedule::slot_to_coeff_plain_log_delta;
    CKKSSparseGaloisKeyGen<P, log_q>(gk, key, usage.slot_to_coeff[I],
                                     noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKeyGen(
    CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKey<Schedule, I> &gk,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    static_assert(I <= Schedule::coeff_to_slot_level_count);
    constexpr std::uint32_t log_q =
        Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta;
    CKKSHybridSparseGaloisKeyGen<P, log_q>(
        gk, key, usage.coeff_to_slot_binary[I],
        usage.coeff_to_slot_direct[I], noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKeyGen(
    CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKey<Schedule, I> &gk,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    static_assert(I <= Schedule::coeff_to_slot_level_count);
    constexpr std::uint32_t log_q =
        Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta;
    CKKSHybridSparseGaloisKeyGen<P, log_q>(
        gk, key, usage.coeff_to_slot_binary[I],
        usage.coeff_to_slot_direct[I], noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKeyGen(
    CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKey<Schedule, I> &gk,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    static_assert(I <= Schedule::slot_to_coeff_level_count);
    constexpr std::uint32_t log_q =
        Schedule::after_evalmod_log_q -
        I * Schedule::slot_to_coeff_plain_log_delta;
    CKKSHybridSparseGaloisKeyGen<P, log_q>(
        gk, key, usage.slot_to_coeff_binary[I],
        usage.slot_to_coeff_direct[I], noise);
}

template <class Schedule, std::size_t I>
inline void CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKeyGen(
    CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKey<Schedule, I> &gk,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    static_assert(I <= Schedule::slot_to_coeff_level_count);
    constexpr std::uint32_t log_q =
        Schedule::after_evalmod_log_q -
        I * Schedule::slot_to_coeff_plain_log_delta;
    CKKSHybridSparseGaloisKeyGen<P, log_q>(
        gk, key, usage.slot_to_coeff_binary[I],
        usage.slot_to_coeff_direct[I], noise);
}

namespace ckks_detail {

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapCoeffToSlotGaloisKeyChainGenImpl(
    typename CKKSDenseBootstrapKey<Schedule>::CoeffToSlotGaloisKeyChain &chain,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        CKKSDenseBootstrapCoeffToSlotGaloisKeyGen<Schedule, I>(
            chain.template get<I>(), usage, key, noise);
        CKKSDenseBootstrapCoeffToSlotGaloisKeyChainGenImpl<I + 1, Schedule>(
            chain, usage, key, noise);
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapSlotToCoeffGaloisKeyChainGenImpl(
    typename CKKSDenseBootstrapKey<Schedule>::SlotToCoeffGaloisKeyChain &chain,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        CKKSDenseBootstrapSlotToCoeffGaloisKeyGen<Schedule, I>(
            chain.template get<I>(), usage, key, noise);
        CKKSDenseBootstrapSlotToCoeffGaloisKeyChainGenImpl<I + 1, Schedule>(
            chain, usage, key, noise);
    }
}

}  // namespace ckks_detail

template <class Schedule>
inline void CKKSDenseBootstrapCoeffToSlotGaloisKeyChainGen(
    typename CKKSDenseBootstrapKey<Schedule>::CoeffToSlotGaloisKeyChain &chain,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    ckks_detail::CKKSDenseBootstrapCoeffToSlotGaloisKeyChainGenImpl<
        0, Schedule>(chain, usage, key, noise);
}

template <class Schedule>
inline void CKKSDenseBootstrapSlotToCoeffGaloisKeyChainGen(
    typename CKKSDenseBootstrapKey<Schedule>::SlotToCoeffGaloisKeyChain &chain,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    ckks_detail::CKKSDenseBootstrapSlotToCoeffGaloisKeyChainGenImpl<
        0, Schedule>(chain, usage, key, noise);
}

template <class Schedule, class BootstrapKey = CKKSDenseBootstrapKey<Schedule>>
struct CKKSDenseBootstrapInMemoryKeyProvider {
    const BootstrapKey *key = nullptr;

    explicit CKKSDenseBootstrapInMemoryKeyProvider(
        const BootstrapKey &bootstrap_key)
        : key(&bootstrap_key)
    {}

    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan() const
    {
        assert(key != nullptr);
        return key->linear_plan;
    }

    const CKKSBoundedCosEvalModPolynomial &evalmod_polynomial() const
    {
        assert(key != nullptr);
        return key->evalmod_polynomial;
    }

    template <std::size_t I>
    const auto &coeff_to_slot_galois() const
    {
        assert(key != nullptr);
        return key->coeff_to_slot_galois.template get<I>();
    }

    const auto &packed_conjugate_galois() const
    {
        assert(key != nullptr);
        return key->packed_conjugate_galois;
    }

    const CKKSDenseEvalModBoundedCosRelinKeys<Schedule> &evalmod_relin() const
    {
        assert(key != nullptr);
        return key->evalmod_relin;
    }

    template <std::size_t I>
    const auto &polynomial_relin() const
    {
        assert(key != nullptr);
        return key->evalmod_relin.polynomial.template get<I>();
    }

    template <std::size_t I>
    const auto &double_angle_relin() const
    {
        assert(key != nullptr);
        return key->evalmod_relin.double_angle.template get<I>();
    }

    template <std::size_t I>
    const auto &inverse_relin() const
    {
        assert(key != nullptr);
        return key->evalmod_relin.inverse.template get<I>();
    }

    template <std::size_t I>
    const auto &slot_to_coeff_galois() const
    {
        assert(key != nullptr);
        return key->slot_to_coeff_galois.template get<I>();
    }
};

template <class Schedule>
inline void CKKSDenseBootstrapKeyGen(
    CKKSDenseBootstrapKey<Schedule> &bootstrap_key,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(bootstrap_key.linear_plan);
    CKKSDenseBootstrapRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapRotationKeyUsage<Schedule>(
        rotation_usage, bootstrap_key.linear_plan);
    bootstrap_key.evalmod_polynomial =
        CKKSBuildBoundedCosEvalModPolynomial<Schedule>();
    CKKSDenseBootstrapCoeffToSlotGaloisKeyChainGen<Schedule>(
        bootstrap_key.coeff_to_slot_galois, rotation_usage, key, noise);
    CKKSDenseBootstrapPackedConjugateGaloisKeyGen<Schedule>(
        bootstrap_key.packed_conjugate_galois, rotation_usage, key, noise);
    CKKSDenseEvalModBoundedCosKeyGen<Schedule>(bootstrap_key.evalmod_relin, key,
                                               noise);
    CKKSDenseBootstrapSlotToCoeffGaloisKeyChainGen<Schedule>(
        bootstrap_key.slot_to_coeff_galois, rotation_usage, key, noise);
}

template <class Schedule>
inline void CKKSDenseBootstrapHybridGiantKeyGen(
    CKKSDenseBootstrapHybridGiantKey<Schedule> &bootstrap_key,
    const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(bootstrap_key.linear_plan);
    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        rotation_usage, bootstrap_key.linear_plan);
    bootstrap_key.evalmod_polynomial =
        CKKSBuildBoundedCosEvalModPolynomial<Schedule>();
    CKKSHybridSparseGaloisKeyChainGen<
        P, Schedule::boot_log_q, Schedule::coeff_to_slot_plain_log_delta,
        Schedule::coeff_to_slot_level_count>(
        bootstrap_key.coeff_to_slot_galois,
        rotation_usage.coeff_to_slot_binary,
        rotation_usage.coeff_to_slot_direct, key, noise);
    CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
        bootstrap_key.packed_conjugate_galois, key,
        rotation_usage.packed_conjugate, noise);
    CKKSDenseEvalModBoundedCosKeyGen<Schedule>(bootstrap_key.evalmod_relin, key,
                                               noise);
    CKKSHybridSparseGaloisKeyChainGen<
        P, Schedule::after_evalmod_log_q,
        Schedule::slot_to_coeff_plain_log_delta,
        Schedule::slot_to_coeff_level_count>(
        bootstrap_key.slot_to_coeff_galois,
        rotation_usage.slot_to_coeff_binary,
        rotation_usage.slot_to_coeff_direct, key, noise);
}

namespace ckks_detail {

inline std::filesystem::path CKKSDenseBootstrapNamedPath(
    const std::filesystem::path &root, const std::string &name)
{
    return root / (name + ".bin");
}

inline std::filesystem::path CKKSDenseBootstrapIndexedPath(
    const std::filesystem::path &root, const std::string &name,
    std::size_t index)
{
    return root / (name + "_" + std::to_string(index) + ".bin");
}

inline bool CKKSDenseBootstrapShouldWriteKeyFile(
    const std::filesystem::path &path,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    return options.overwrite_existing || !std::filesystem::exists(path);
}

inline std::string CKKSDenseBootstrapLeveledKeyPrefix(
    const std::string &name, std::size_t level, const std::string &kind)
{
    return name + "_" + std::to_string(level) + "_" + kind;
}

template <class P, std::uint32_t LogQ>
struct CKKSSeededSparseGaloisKeyFileView {
    std::filesystem::path root{};
    std::string prefix{};
    CKKSRotationKeyIndexSet<P> available{};
    mutable std::array<std::unique_ptr<CKKSSeededAutoKey<P, LogQ>>,
                       P::nbit + 1>
        keys{};

    bool has(std::size_t i) const
    {
        return i < keys.size() && available[i];
    }

    const CKKSSeededAutoKey<P, LogQ> &get(std::size_t i) const
    {
        assert(has(i));
        auto &entry = keys[i];
        if (!entry) {
            entry = std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
            CKKSLoadPortableBinary(*entry,
                                   CKKSDenseBootstrapIndexedPath(root, prefix,
                                                                 i));
        }
        return *entry;
    }

    void clear_cached_keys() const
    {
        for (auto &entry : keys) entry.reset();
    }
};

template <class P, std::uint32_t LogQ>
struct CKKSSeededDirectSparseGaloisKeyFileView {
    std::filesystem::path root{};
    std::string prefix{};
    CKKSDirectRotationKeyIndexSet<P> available{};
    mutable std::vector<CKKSSeededDirectSparseGaloisKeyEntry<P, LogQ>> keys{};

    static constexpr int normalize_steps(int steps)
    {
        constexpr int half = static_cast<int>(P::n) / 2;
        return ((steps % half) + half) % half;
    }

    bool has(int steps) const
    {
        steps = normalize_steps(steps);
        return steps == 0 || available[static_cast<std::size_t>(steps)];
    }

    const CKKSSeededAutoKey<P, LogQ> &get(int steps) const
    {
        steps = normalize_steps(steps);
        assert(steps != 0);
        assert(available[static_cast<std::size_t>(steps)]);
        auto it = std::find_if(keys.begin(), keys.end(),
                               [&](const auto &entry) {
                                   return entry.steps == steps;
                               });
        if (it == keys.end()) {
            auto &entry = keys.emplace_back();
            entry.steps = steps;
            CKKSLoadPortableBinary(entry.key,
                                   CKKSDenseBootstrapIndexedPath(
                                       root, prefix,
                                       static_cast<std::size_t>(steps)));
            return entry.key;
        }
        return it->key;
    }

    void clear_cached_keys() const { keys.clear(); }
};

template <class P, std::uint32_t LogQ>
struct CKKSSeededHybridSparseGaloisKeyFileView {
    CKKSSeededSparseGaloisKeyFileView<P, LogQ> binary{};
    CKKSSeededDirectSparseGaloisKeyFileView<P, LogQ> direct{};

    bool has(int steps) const { return direct.has(steps) || binary.has(steps); }

    void clear_cached_keys() const
    {
        binary.clear_cached_keys();
        direct.clear_cached_keys();
    }
};

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                            const CKKSSeededSparseGaloisKeyFileView<P, LogQ>
                                &gk)
{
    CKKSRotateSlotsBinaryChain<P, LogQ>(
        res, ct, steps,
        [&](std::size_t i) -> const CKKSSeededAutoKey<P, LogQ> & {
            return gk.get(i);
        });
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(
    TRLWE<P> &res, const TRLWE<P> &ct, int steps,
    const CKKSSeededDirectSparseGaloisKeyFileView<P, LogQ> &gk)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;
    if (steps == 0) {
        res = ct;
        reduceTRLWEToLevel<P, LogQ>(res);
        return;
    }

    CKKSEvalAuto<P, LogQ>(res, ct, rotationAutomorphismForSteps<P>(steps),
                          gk.get(steps));
    reduceTRLWEToLevel<P, LogQ>(res);
}

template <class P, std::uint32_t LogQ>
inline void CKKSRotateSlots(
    TRLWE<P> &res, const TRLWE<P> &ct, int steps,
    const CKKSSeededHybridSparseGaloisKeyFileView<P, LogQ> &gk)
{
    if (gk.direct.has(steps) &&
        CKKSSeededDirectSparseGaloisKeyFileView<P, LogQ>::normalize_steps(
            steps) != 0) {
        CKKSRotateSlots<P, LogQ>(res, ct, steps, gk.direct);
    }
    else {
        CKKSRotateSlots<P, LogQ>(res, ct, steps, gk.binary);
    }
}

template <class P, std::uint32_t LogQ>
inline void CKKSSeededSparseGaloisKeyFilePaths(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root,
    const std::string &prefix, const CKKSRotationKeyIndexSet<P> &indices)
{
    for (std::size_t i = 0; i < indices.size(); i++)
        if (indices[i]) paths.push_back(CKKSDenseBootstrapIndexedPath(root,
                                                                      prefix,
                                                                      i));
}

template <class P, std::uint32_t LogQ>
inline void CKKSSeededDirectSparseGaloisKeyFilePaths(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root,
    const std::string &prefix, const CKKSDirectRotationKeyIndexSet<P> &indices)
{
    for (std::size_t steps = 1; steps < indices.size(); steps++)
        if (indices[steps])
            paths.push_back(CKKSDenseBootstrapIndexedPath(root, prefix, steps));
}

template <class P, std::uint32_t LogQ>
inline void CKKSSeededSparseGaloisKeyGenToDirectory(
    const std::filesystem::path &root, const std::string &prefix,
    const CKKSRotationKeyIndexSet<P> &indices, const Key<P> &key,
    CKKSNoise noise, const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    std::uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit); i++) {
        if (indices[static_cast<std::size_t>(i)]) {
            const std::filesystem::path path =
                CKKSDenseBootstrapIndexedPath(root, prefix,
                                              static_cast<std::size_t>(i));
            if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
                auto autokey =
                    std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
                CKKSAutoKeyGen<P, LogQ>(*autokey, static_cast<uint>(d), key,
                                         noise);
                CKKSSavePortableBinaryAtomic(path, *autokey);
            }
        }
        d = d * d % (2 * P::n);
    }

    if (indices[P::nbit]) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, prefix, P::nbit);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto autokey = std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
            CKKSAutoKeyGen<P, LogQ>(*autokey, 2 * P::n - 1, key, noise);
            CKKSSavePortableBinaryAtomic(path, *autokey);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline bool CKKSSeededSparseGaloisKeyGenNextMissingToDirectory(
    const std::filesystem::path &root, const std::string &prefix,
    const CKKSRotationKeyIndexSet<P> &indices, const Key<P> &key,
    CKKSNoise noise)
{
    std::uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit); i++) {
        if (indices[static_cast<std::size_t>(i)]) {
            const std::filesystem::path path =
                CKKSDenseBootstrapIndexedPath(root, prefix,
                                              static_cast<std::size_t>(i));
            if (!std::filesystem::exists(path)) {
                auto autokey =
                    std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
                CKKSAutoKeyGen<P, LogQ>(*autokey, static_cast<uint>(d), key,
                                         noise);
                CKKSSavePortableBinaryAtomic(path, *autokey);
                return true;
            }
        }
        d = d * d % (2 * P::n);
    }

    if (indices[P::nbit]) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, prefix, P::nbit);
        if (!std::filesystem::exists(path)) {
            auto autokey = std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
            CKKSAutoKeyGen<P, LogQ>(*autokey, 2 * P::n - 1, key, noise);
            CKKSSavePortableBinaryAtomic(path, *autokey);
            return true;
        }
    }
    return false;
}

template <class P, std::uint32_t LogQ>
inline void CKKSSeededDirectSparseGaloisKeyGenToDirectory(
    const std::filesystem::path &root, const std::string &prefix,
    const CKKSDirectRotationKeyIndexSet<P> &indices, const Key<P> &key,
    CKKSNoise noise, const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    for (int steps = 1; steps < static_cast<int>(P::n / 2); steps++) {
        if (!indices[static_cast<std::size_t>(steps)]) continue;
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, prefix,
                                          static_cast<std::size_t>(steps));
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto autokey = std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
            CKKSAutoKeyGen<P, LogQ>(
                *autokey, rotationAutomorphismForSteps<P>(steps), key, noise);
            CKKSSavePortableBinaryAtomic(path, *autokey);
        }
    }
}

template <class P, std::uint32_t LogQ>
inline bool CKKSSeededDirectSparseGaloisKeyGenNextMissingToDirectory(
    const std::filesystem::path &root, const std::string &prefix,
    const CKKSDirectRotationKeyIndexSet<P> &indices, const Key<P> &key,
    CKKSNoise noise)
{
    for (int steps = 1; steps < static_cast<int>(P::n / 2); steps++) {
        if (!indices[static_cast<std::size_t>(steps)]) continue;
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, prefix,
                                          static_cast<std::size_t>(steps));
        if (!std::filesystem::exists(path)) {
            auto autokey = std::make_unique<CKKSSeededAutoKey<P, LogQ>>();
            CKKSAutoKeyGen<P, LogQ>(
                *autokey, rotationAutomorphismForSteps<P>(steps), key, noise);
            CKKSSavePortableBinaryAtomic(path, *autokey);
            return true;
        }
    }
    return false;
}

template <std::size_t I, class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyFilePathsImpl(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        constexpr std::uint32_t log_q =
            Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta;
        CKKSSeededSparseGaloisKeyFilePaths<P, log_q>(
            paths, root,
            CKKSDenseBootstrapLeveledKeyPrefix("coeff_to_slot_galois", I,
                                               "binary"),
            usage.coeff_to_slot_binary[I]);
        CKKSSeededDirectSparseGaloisKeyFilePaths<P, log_q>(
            paths, root,
            CKKSDenseBootstrapLeveledKeyPrefix("coeff_to_slot_galois", I,
                                               "direct"),
            usage.coeff_to_slot_direct[I]);
        CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyFilePathsImpl<
            I + 1, Schedule>(paths, root, usage);
    }
}

template <std::size_t I, class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyFilePathsImpl(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        constexpr std::uint32_t log_q =
            Schedule::after_evalmod_log_q -
            I * Schedule::slot_to_coeff_plain_log_delta;
        CKKSSeededSparseGaloisKeyFilePaths<P, log_q>(
            paths, root,
            CKKSDenseBootstrapLeveledKeyPrefix("slot_to_coeff_galois", I,
                                               "binary"),
            usage.slot_to_coeff_binary[I]);
        CKKSSeededDirectSparseGaloisKeyFilePaths<P, log_q>(
            paths, root,
            CKKSDenseBootstrapLeveledKeyPrefix("slot_to_coeff_galois", I,
                                               "direct"),
            usage.slot_to_coeff_direct[I]);
        CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyFilePathsImpl<
            I + 1, Schedule>(paths, root, usage);
    }
}

template <std::size_t I, class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyGenToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    using P = typename Schedule::Param;
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        constexpr std::uint32_t log_q =
            Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta;
        CKKSSeededSparseGaloisKeyGenToDirectory<P, log_q>(
            root,
            CKKSDenseBootstrapLeveledKeyPrefix("coeff_to_slot_galois", I,
                                               "binary"),
            usage.coeff_to_slot_binary[I], key, noise, options);
        CKKSSeededDirectSparseGaloisKeyGenToDirectory<P, log_q>(
            root,
            CKKSDenseBootstrapLeveledKeyPrefix("coeff_to_slot_galois", I,
                                               "direct"),
            usage.coeff_to_slot_direct[I], key, noise, options);
        CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool
CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    using P = typename Schedule::Param;
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        constexpr std::uint32_t log_q =
            Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta;
        if (CKKSSeededSparseGaloisKeyGenNextMissingToDirectory<P, log_q>(
                root,
                CKKSDenseBootstrapLeveledKeyPrefix("coeff_to_slot_galois", I,
                                                   "binary"),
                usage.coeff_to_slot_binary[I], key, noise))
            return true;
        if (CKKSSeededDirectSparseGaloisKeyGenNextMissingToDirectory<P, log_q>(
                root,
                CKKSDenseBootstrapLeveledKeyPrefix("coeff_to_slot_galois", I,
                                                   "direct"),
                usage.coeff_to_slot_direct[I], key, noise))
            return true;
        return CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyGenToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    using P = typename Schedule::Param;
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        constexpr std::uint32_t log_q =
            Schedule::after_evalmod_log_q -
            I * Schedule::slot_to_coeff_plain_log_delta;
        CKKSSeededSparseGaloisKeyGenToDirectory<P, log_q>(
            root,
            CKKSDenseBootstrapLeveledKeyPrefix("slot_to_coeff_galois", I,
                                               "binary"),
            usage.slot_to_coeff_binary[I], key, noise, options);
        CKKSSeededDirectSparseGaloisKeyGenToDirectory<P, log_q>(
            root,
            CKKSDenseBootstrapLeveledKeyPrefix("slot_to_coeff_galois", I,
                                               "direct"),
            usage.slot_to_coeff_direct[I], key, noise, options);
        CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool
CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    using P = typename Schedule::Param;
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        constexpr std::uint32_t log_q =
            Schedule::after_evalmod_log_q -
            I * Schedule::slot_to_coeff_plain_log_delta;
        if (CKKSSeededSparseGaloisKeyGenNextMissingToDirectory<P, log_q>(
                root,
                CKKSDenseBootstrapLeveledKeyPrefix("slot_to_coeff_galois", I,
                                                   "binary"),
                usage.slot_to_coeff_binary[I], key, noise))
            return true;
        if (CKKSSeededDirectSparseGaloisKeyGenNextMissingToDirectory<P, log_q>(
                root,
                CKKSDenseBootstrapLeveledKeyPrefix("slot_to_coeff_galois", I,
                                                   "direct"),
                usage.slot_to_coeff_direct[I], key, noise))
            return true;
        return CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapCoeffToSlotKeyGenToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "coeff_to_slot_galois", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto gk =
                std::make_unique<CKKSDenseBootstrapCoeffToSlotGaloisKey<
                    Schedule, I>>();
            CKKSDenseBootstrapCoeffToSlotGaloisKeyGen<Schedule, I>(
                *gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
        }
        CKKSDenseBootstrapCoeffToSlotKeyGenToDirectoryImpl<I + 1, Schedule>(
            root, usage, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool CKKSDenseBootstrapCoeffToSlotKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "coeff_to_slot_galois", I);
        if (!std::filesystem::exists(path)) {
            auto gk =
                std::make_unique<CKKSDenseBootstrapCoeffToSlotGaloisKey<
                    Schedule, I>>();
            CKKSDenseBootstrapCoeffToSlotGaloisKeyGen<Schedule, I>(
                *gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
            return true;
        }
        return CKKSDenseBootstrapCoeffToSlotKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapSlotToCoeffKeyGenToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "slot_to_coeff_galois", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto gk =
                std::make_unique<CKKSDenseBootstrapSlotToCoeffGaloisKey<
                    Schedule, I>>();
            CKKSDenseBootstrapSlotToCoeffGaloisKeyGen<Schedule, I>(
                *gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
        }
        CKKSDenseBootstrapSlotToCoeffKeyGenToDirectoryImpl<I + 1, Schedule>(
            root, usage, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool CKKSDenseBootstrapSlotToCoeffKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "slot_to_coeff_galois", I);
        if (!std::filesystem::exists(path)) {
            auto gk =
                std::make_unique<CKKSDenseBootstrapSlotToCoeffGaloisKey<
                    Schedule, I>>();
            CKKSDenseBootstrapSlotToCoeffGaloisKeyGen<Schedule, I>(
                *gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
            return true;
        }
        return CKKSDenseBootstrapSlotToCoeffKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapHybridGiantCoeffToSlotKeyGenToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "coeff_to_slot_galois", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto gk = std::make_unique<
                CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKey<Schedule,
                                                                  I>>();
            CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKeyGen<Schedule,
                                                                 I>(
                *gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
        }
        CKKSDenseBootstrapHybridGiantCoeffToSlotKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool
CKKSDenseBootstrapHybridGiantCoeffToSlotKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "coeff_to_slot_galois", I);
        if (!std::filesystem::exists(path)) {
            auto gk = std::make_unique<
                CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKey<Schedule,
                                                                  I>>();
            CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKeyGen<Schedule,
                                                                 I>(
                *gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
            return true;
        }
        return CKKSDenseBootstrapHybridGiantCoeffToSlotKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapHybridGiantSlotToCoeffKeyGenToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "slot_to_coeff_galois", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto gk = std::make_unique<
                CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKey<Schedule,
                                                                  I>>();
            CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKeyGen<Schedule,
                                                                 I>(
                *gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
        }
        CKKSDenseBootstrapHybridGiantSlotToCoeffKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool
CKKSDenseBootstrapHybridGiantSlotToCoeffKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "slot_to_coeff_galois", I);
        if (!std::filesystem::exists(path)) {
            auto gk = std::make_unique<
                CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKey<Schedule,
                                                                  I>>();
            CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKeyGen<Schedule,
                                                                 I>(
                *gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
            return true;
        }
        return CKKSDenseBootstrapHybridGiantSlotToCoeffKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantCoeffToSlotKeyGenToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "coeff_to_slot_galois", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto gk = std::make_unique<
                CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKey<
                    Schedule, I>>();
            CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKeyGen<
                Schedule, I>(*gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
        }
        CKKSDenseBootstrapSeededHybridGiantCoeffToSlotKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool
CKKSDenseBootstrapSeededHybridGiantCoeffToSlotKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "coeff_to_slot_galois", I);
        if (!std::filesystem::exists(path)) {
            auto gk = std::make_unique<
                CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKey<
                    Schedule, I>>();
            CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKeyGen<
                Schedule, I>(*gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
            return true;
        }
        return CKKSDenseBootstrapSeededHybridGiantCoeffToSlotKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantSlotToCoeffKeyGenToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise,
    const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "slot_to_coeff_galois", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto gk = std::make_unique<
                CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKey<
                    Schedule, I>>();
            CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKeyGen<
                Schedule, I>(*gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
        }
        CKKSDenseBootstrapSeededHybridGiantSlotToCoeffKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool
CKKSDenseBootstrapSeededHybridGiantSlotToCoeffKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage,
    const Key<typename Schedule::Param> &key, CKKSNoise noise)
{
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "slot_to_coeff_galois", I);
        if (!std::filesystem::exists(path)) {
            auto gk = std::make_unique<
                CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKey<
                    Schedule, I>>();
            CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKeyGen<
                Schedule, I>(*gk, usage, key, noise);
            CKKSSavePortableBinaryAtomic(path, *gk);
            return true;
        }
        return CKKSDenseBootstrapSeededHybridGiantSlotToCoeffKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, usage, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapPolynomialRelinKeyGenToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise, const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::PolynomialTraits::relin_depth) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "polynomial_relin", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto relinkey =
                std::make_unique<CKKSDenseEvalModPolynomialRelinKey<Schedule,
                                                                    I>>();
            CKKSDenseEvalModPolynomialRelinKeyGen<Schedule, I>(*relinkey, key,
                                                               noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
        }
        CKKSDenseBootstrapPolynomialRelinKeyGenToDirectoryImpl<I + 1, Schedule>(
            root, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool CKKSDenseBootstrapPolynomialRelinKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::PolynomialTraits::relin_depth) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "polynomial_relin", I);
        if (!std::filesystem::exists(path)) {
            auto relinkey =
                std::make_unique<CKKSDenseEvalModPolynomialRelinKey<Schedule,
                                                                    I>>();
            CKKSDenseEvalModPolynomialRelinKeyGen<Schedule, I>(*relinkey, key,
                                                               noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
            return true;
        }
        return CKKSDenseBootstrapPolynomialRelinKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapSeededPolynomialRelinKeyGenToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise, const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::PolynomialTraits::relin_depth) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "polynomial_relin", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto relinkey = std::make_unique<
                CKKSDenseSeededEvalModPolynomialRelinKey<Schedule, I>>();
            CKKSDenseSeededEvalModPolynomialRelinKeyGen<Schedule, I>(
                *relinkey, key, noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
        }
        CKKSDenseBootstrapSeededPolynomialRelinKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool
CKKSDenseBootstrapSeededPolynomialRelinKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::PolynomialTraits::relin_depth) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "polynomial_relin", I);
        if (!std::filesystem::exists(path)) {
            auto relinkey = std::make_unique<
                CKKSDenseSeededEvalModPolynomialRelinKey<Schedule, I>>();
            CKKSDenseSeededEvalModPolynomialRelinKeyGen<Schedule, I>(
                *relinkey, key, noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
            return true;
        }
        return CKKSDenseBootstrapSeededPolynomialRelinKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapDoubleAngleRelinKeyGenToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise, const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    if constexpr (I < Schedule::evalmod_double_angle) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "double_angle_relin", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto relinkey =
                std::make_unique<CKKSDenseEvalModDoubleAngleRelinKey<
                    Schedule, I>>();
            CKKSDenseEvalModDoubleAngleRelinKeyGen<Schedule, I>(
                *relinkey, key, noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
        }
        CKKSDenseBootstrapDoubleAngleRelinKeyGenToDirectoryImpl<I + 1,
                                                                Schedule>(
            root, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool CKKSDenseBootstrapDoubleAngleRelinKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise)
{
    if constexpr (I < Schedule::evalmod_double_angle) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "double_angle_relin", I);
        if (!std::filesystem::exists(path)) {
            auto relinkey =
                std::make_unique<CKKSDenseEvalModDoubleAngleRelinKey<
                    Schedule, I>>();
            CKKSDenseEvalModDoubleAngleRelinKeyGen<Schedule, I>(
                *relinkey, key, noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
            return true;
        }
        return CKKSDenseBootstrapDoubleAngleRelinKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapSeededDoubleAngleRelinKeyGenToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise, const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    if constexpr (I < Schedule::evalmod_double_angle) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "double_angle_relin", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto relinkey = std::make_unique<
                CKKSDenseSeededEvalModDoubleAngleRelinKey<Schedule, I>>();
            CKKSDenseSeededEvalModDoubleAngleRelinKeyGen<Schedule, I>(
                *relinkey, key, noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
        }
        CKKSDenseBootstrapSeededDoubleAngleRelinKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool
CKKSDenseBootstrapSeededDoubleAngleRelinKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise)
{
    if constexpr (I < Schedule::evalmod_double_angle) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "double_angle_relin", I);
        if (!std::filesystem::exists(path)) {
            auto relinkey = std::make_unique<
                CKKSDenseSeededEvalModDoubleAngleRelinKey<Schedule, I>>();
            CKKSDenseSeededEvalModDoubleAngleRelinKeyGen<Schedule, I>(
                *relinkey, key, noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
            return true;
        }
        return CKKSDenseBootstrapSeededDoubleAngleRelinKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapInverseRelinKeyGenToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise, const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::InverseTraits::power_depth) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "inverse_relin", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto relinkey =
                std::make_unique<CKKSDenseEvalModInverseRelinKey<Schedule,
                                                                 I>>();
            CKKSDenseEvalModInverseRelinKeyGen<Schedule, I>(*relinkey, key,
                                                            noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
        }
        CKKSDenseBootstrapInverseRelinKeyGenToDirectoryImpl<I + 1, Schedule>(
            root, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool CKKSDenseBootstrapInverseRelinKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::InverseTraits::power_depth) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "inverse_relin", I);
        if (!std::filesystem::exists(path)) {
            auto relinkey =
                std::make_unique<CKKSDenseEvalModInverseRelinKey<Schedule,
                                                                 I>>();
            CKKSDenseEvalModInverseRelinKeyGen<Schedule, I>(*relinkey, key,
                                                            noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
            return true;
        }
        return CKKSDenseBootstrapInverseRelinKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapSeededInverseRelinKeyGenToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise, const CKKSDenseBootstrapKeyDirectoryOptions &options)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::InverseTraits::power_depth) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "inverse_relin", I);
        if (CKKSDenseBootstrapShouldWriteKeyFile(path, options)) {
            auto relinkey = std::make_unique<
                CKKSDenseSeededEvalModInverseRelinKey<Schedule, I>>();
            CKKSDenseSeededEvalModInverseRelinKeyGen<Schedule, I>(
                *relinkey, key, noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
        }
        CKKSDenseBootstrapSeededInverseRelinKeyGenToDirectoryImpl<
            I + 1, Schedule>(root, key, noise, options);
    }
}

template <std::size_t I, class Schedule>
inline bool CKKSDenseBootstrapSeededInverseRelinKeyGenNextMissingToDirectoryImpl(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::InverseTraits::power_depth) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "inverse_relin", I);
        if (!std::filesystem::exists(path)) {
            auto relinkey = std::make_unique<
                CKKSDenseSeededEvalModInverseRelinKey<Schedule, I>>();
            CKKSDenseSeededEvalModInverseRelinKeyGen<Schedule, I>(
                *relinkey, key, noise);
            CKKSSavePortableBinaryAtomic(path, *relinkey);
            return true;
        }
        return CKKSDenseBootstrapSeededInverseRelinKeyGenNextMissingToDirectoryImpl<
            I + 1, Schedule>(root, key, noise);
    }
    else {
        return false;
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapCoeffToSlotKeyFilePathsImpl(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root)
{
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        paths.push_back(
            CKKSDenseBootstrapIndexedPath(root, "coeff_to_slot_galois", I));
        CKKSDenseBootstrapCoeffToSlotKeyFilePathsImpl<I + 1, Schedule>(paths,
                                                                       root);
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapSlotToCoeffKeyFilePathsImpl(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root)
{
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        paths.push_back(
            CKKSDenseBootstrapIndexedPath(root, "slot_to_coeff_galois", I));
        CKKSDenseBootstrapSlotToCoeffKeyFilePathsImpl<I + 1, Schedule>(paths,
                                                                       root);
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapPolynomialRelinKeyFilePathsImpl(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::PolynomialTraits::relin_depth) {
        paths.push_back(CKKSDenseBootstrapIndexedPath(root, "polynomial_relin",
                                                      I));
        CKKSDenseBootstrapPolynomialRelinKeyFilePathsImpl<I + 1, Schedule>(
            paths, root);
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapDoubleAngleRelinKeyFilePathsImpl(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root)
{
    if constexpr (I < Schedule::evalmod_double_angle) {
        paths.push_back(
            CKKSDenseBootstrapIndexedPath(root, "double_angle_relin", I));
        CKKSDenseBootstrapDoubleAngleRelinKeyFilePathsImpl<I + 1, Schedule>(
            paths, root);
    }
}

template <std::size_t I, class Schedule>
inline void CKKSDenseBootstrapInverseRelinKeyFilePathsImpl(
    std::vector<std::filesystem::path> &paths, const std::filesystem::path &root)
{
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::InverseTraits::power_depth) {
        paths.push_back(
            CKKSDenseBootstrapIndexedPath(root, "inverse_relin", I));
        CKKSDenseBootstrapInverseRelinKeyFilePathsImpl<I + 1, Schedule>(
            paths, root);
    }
}

template <class P, std::uint32_t LogQ>
inline std::optional<CKKSDenseBootstrapKeyFileByteEstimate>
CKKSSeededSparseGaloisKeyNextMissingFileByteEstimate(
    const std::filesystem::path &root, const std::string &prefix,
    const CKKSRotationKeyIndexSet<P> &indices)
{
    for (std::size_t i = 0; i < P::nbit; i++) {
        if (!indices[i]) continue;
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, prefix, i);
        if (!std::filesystem::exists(path))
            return CKKSDenseBootstrapKeyFileByteEstimate{
                path, CKKSSeededAutoKeyByteEstimate<P, LogQ>()};
    }

    if (indices[P::nbit]) {
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, prefix, P::nbit);
        if (!std::filesystem::exists(path))
            return CKKSDenseBootstrapKeyFileByteEstimate{
                path, CKKSSeededAutoKeyByteEstimate<P, LogQ>()};
    }
    return std::nullopt;
}

template <class P, std::uint32_t LogQ>
inline std::optional<CKKSDenseBootstrapKeyFileByteEstimate>
CKKSSeededDirectSparseGaloisKeyNextMissingFileByteEstimate(
    const std::filesystem::path &root, const std::string &prefix,
    const CKKSDirectRotationKeyIndexSet<P> &indices)
{
    for (std::size_t steps = 1; steps < P::n / 2; steps++) {
        if (!indices[steps]) continue;
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, prefix, steps);
        if (!std::filesystem::exists(path))
            return CKKSDenseBootstrapKeyFileByteEstimate{
                path, CKKSSeededAutoKeyByteEstimate<P, LogQ>()};
    }
    return std::nullopt;
}

template <std::size_t I, class Schedule>
inline std::optional<CKKSDenseBootstrapKeyFileByteEstimate>
CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotNextMissingFileByteEstimateImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    if constexpr (I <= Schedule::coeff_to_slot_level_count) {
        constexpr std::uint32_t log_q =
            Schedule::boot_log_q - I * Schedule::coeff_to_slot_plain_log_delta;
        if (auto next =
                CKKSSeededSparseGaloisKeyNextMissingFileByteEstimate<P, log_q>(
                    root,
                    CKKSDenseBootstrapLeveledKeyPrefix("coeff_to_slot_galois",
                                                       I, "binary"),
                    usage.coeff_to_slot_binary[I]))
            return next;
        if (auto next =
                CKKSSeededDirectSparseGaloisKeyNextMissingFileByteEstimate<
                    P, log_q>(
                    root,
                    CKKSDenseBootstrapLeveledKeyPrefix("coeff_to_slot_galois",
                                                       I, "direct"),
                    usage.coeff_to_slot_direct[I]))
            return next;
        return CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotNextMissingFileByteEstimateImpl<
            I + 1, Schedule>(root, usage);
    }
    else {
        return std::nullopt;
    }
}

template <std::size_t I, class Schedule>
inline std::optional<CKKSDenseBootstrapKeyFileByteEstimate>
CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffNextMissingFileByteEstimateImpl(
    const std::filesystem::path &root,
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> &usage)
{
    using P = typename Schedule::Param;
    if constexpr (I <= Schedule::slot_to_coeff_level_count) {
        constexpr std::uint32_t log_q =
            Schedule::after_evalmod_log_q -
            I * Schedule::slot_to_coeff_plain_log_delta;
        if (auto next =
                CKKSSeededSparseGaloisKeyNextMissingFileByteEstimate<P, log_q>(
                    root,
                    CKKSDenseBootstrapLeveledKeyPrefix("slot_to_coeff_galois",
                                                       I, "binary"),
                    usage.slot_to_coeff_binary[I]))
            return next;
        if (auto next =
                CKKSSeededDirectSparseGaloisKeyNextMissingFileByteEstimate<
                    P, log_q>(
                    root,
                    CKKSDenseBootstrapLeveledKeyPrefix("slot_to_coeff_galois",
                                                       I, "direct"),
                    usage.slot_to_coeff_direct[I]))
            return next;
        return CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffNextMissingFileByteEstimateImpl<
            I + 1, Schedule>(root, usage);
    }
    else {
        return std::nullopt;
    }
}

template <std::size_t I, class Schedule>
inline std::optional<CKKSDenseBootstrapKeyFileByteEstimate>
CKKSDenseBootstrapSeededPolynomialRelinNextMissingFileByteEstimateImpl(
    const std::filesystem::path &root)
{
    using P = typename Schedule::Param;
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::PolynomialTraits::relin_depth) {
        constexpr std::uint32_t log_q =
            Schedule::after_component_split_log_q -
            (I + 1) * Schedule::log_delta;
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "polynomial_relin", I);
        if (!std::filesystem::exists(path))
            return CKKSDenseBootstrapKeyFileByteEstimate{
                path, CKKSSeededRelinKeyByteEstimate<P, log_q>()};
        return CKKSDenseBootstrapSeededPolynomialRelinNextMissingFileByteEstimateImpl<
            I + 1, Schedule>(root);
    }
    else {
        return std::nullopt;
    }
}

template <std::size_t I, class Schedule>
inline std::optional<CKKSDenseBootstrapKeyFileByteEstimate>
CKKSDenseBootstrapSeededDoubleAngleRelinNextMissingFileByteEstimateImpl(
    const std::filesystem::path &root)
{
    using P = typename Schedule::Param;
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Schedule::evalmod_double_angle) {
        constexpr std::uint32_t log_q = Traits::polynomial_log_q -
                                        (I + 1) * Schedule::log_delta;
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "double_angle_relin", I);
        if (!std::filesystem::exists(path))
            return CKKSDenseBootstrapKeyFileByteEstimate{
                path, CKKSSeededRelinKeyByteEstimate<P, log_q>()};
        return CKKSDenseBootstrapSeededDoubleAngleRelinNextMissingFileByteEstimateImpl<
            I + 1, Schedule>(root);
    }
    else {
        return std::nullopt;
    }
}

template <std::size_t I, class Schedule>
inline std::optional<CKKSDenseBootstrapKeyFileByteEstimate>
CKKSDenseBootstrapSeededInverseRelinNextMissingFileByteEstimateImpl(
    const std::filesystem::path &root)
{
    using P = typename Schedule::Param;
    using Traits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    if constexpr (I < Traits::InverseTraits::power_depth) {
        constexpr std::uint32_t log_q = Traits::after_double_angle_log_q -
                                        (I + 1) * Schedule::log_delta;
        const std::filesystem::path path =
            CKKSDenseBootstrapIndexedPath(root, "inverse_relin", I);
        if (!std::filesystem::exists(path))
            return CKKSDenseBootstrapKeyFileByteEstimate{
                path, CKKSSeededRelinKeyByteEstimate<P, log_q>()};
        return CKKSDenseBootstrapSeededInverseRelinNextMissingFileByteEstimateImpl<
            I + 1, Schedule>(root);
    }
    else {
        return std::nullopt;
    }
}

template <class Schedule, class Seq>
struct CKKSDenseBootstrapCoeffToSlotCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapCoeffToSlotCacheTuple<Schedule,
                                               std::index_sequence<Is...>> {
    using type = std::tuple<
        std::unique_ptr<CKKSDenseBootstrapCoeffToSlotGaloisKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapHybridGiantCoeffToSlotCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapHybridGiantCoeffToSlotCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using type = std::tuple<std::unique_ptr<
        CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapSeededHybridGiantCoeffToSlotCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapSeededHybridGiantCoeffToSlotCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using type = std::tuple<std::unique_ptr<
        CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKey<Schedule,
                                                               Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using P = typename Schedule::Param;
    using type = std::tuple<std::unique_ptr<
        CKKSSeededHybridSparseGaloisKeyFileView<
            P, Schedule::boot_log_q -
                   Is * Schedule::coeff_to_slot_plain_log_delta>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapSlotToCoeffCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapSlotToCoeffCacheTuple<Schedule,
                                               std::index_sequence<Is...>> {
    using type = std::tuple<
        std::unique_ptr<CKKSDenseBootstrapSlotToCoeffGaloisKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapHybridGiantSlotToCoeffCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapHybridGiantSlotToCoeffCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using type = std::tuple<std::unique_ptr<
        CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapSeededHybridGiantSlotToCoeffCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapSeededHybridGiantSlotToCoeffCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using type = std::tuple<std::unique_ptr<
        CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKey<Schedule,
                                                               Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using P = typename Schedule::Param;
    using type = std::tuple<std::unique_ptr<
        CKKSSeededHybridSparseGaloisKeyFileView<
            P, Schedule::after_evalmod_log_q -
                   Is * Schedule::slot_to_coeff_plain_log_delta>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapPolynomialRelinCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapPolynomialRelinCacheTuple<Schedule,
                                                   std::index_sequence<Is...>> {
    using type =
        std::tuple<std::unique_ptr<
            CKKSDenseEvalModPolynomialRelinKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapSeededPolynomialRelinCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapSeededPolynomialRelinCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using type =
        std::tuple<std::unique_ptr<
            CKKSDenseSeededEvalModPolynomialRelinKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapDoubleAngleRelinCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapDoubleAngleRelinCacheTuple<Schedule,
                                                    std::index_sequence<Is...>> {
    using type =
        std::tuple<std::unique_ptr<
            CKKSDenseEvalModDoubleAngleRelinKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapSeededDoubleAngleRelinCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapSeededDoubleAngleRelinCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using type =
        std::tuple<std::unique_ptr<
            CKKSDenseSeededEvalModDoubleAngleRelinKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapInverseRelinCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapInverseRelinCacheTuple<Schedule,
                                                std::index_sequence<Is...>> {
    using type =
        std::tuple<std::unique_ptr<
            CKKSDenseEvalModInverseRelinKey<Schedule, Is>>...>;
};

template <class Schedule, class Seq>
struct CKKSDenseBootstrapSeededInverseRelinCacheTuple;

template <class Schedule, std::size_t... Is>
struct CKKSDenseBootstrapSeededInverseRelinCacheTuple<
    Schedule, std::index_sequence<Is...>> {
    using type =
        std::tuple<std::unique_ptr<
            CKKSDenseSeededEvalModInverseRelinKey<Schedule, Is>>...>;
};

}  // namespace ckks_detail

inline std::filesystem::path CKKSDenseBootstrapKeyDirectoryManifestFile(
    const std::filesystem::path &root)
{
    return ckks_detail::CKKSDenseBootstrapNamedPath(root, "manifest");
}

inline std::filesystem::path CKKSDenseBootstrapEncapsulationKeyFile(
    const std::filesystem::path &root)
{
    return ckks_detail::CKKSDenseBootstrapNamedPath(root, "encapsulation_key");
}

inline std::filesystem::path CKKSDenseBootstrapSeededEncapsulationKeyFile(
    const std::filesystem::path &root)
{
    return ckks_detail::CKKSDenseBootstrapNamedPath(
        root, "seeded_encapsulation_key");
}

inline std::filesystem::path CKKSRelinKeyFile(const std::filesystem::path &root,
                                              const std::string &name =
                                                  "relin_key")
{
    return ckks_detail::CKKSDenseBootstrapNamedPath(root, name);
}

template <class P, std::uint32_t LogQ>
inline void CKKSRelinKeyGenToFile(
    const std::filesystem::path &path, const Key<P> &key,
    CKKSNoise noise = {P::α, 0},
    CKKSDenseBootstrapKeyDirectoryOptions options = {})
{
    if (!options.overwrite_existing && std::filesystem::exists(path)) return;

    auto relinkey = makeCKKSRelinKey<P, LogQ>(key, noise);
    CKKSSavePortableBinaryAtomic(path, *relinkey);
}

template <class P, std::uint32_t LogQ>
inline void CKKSSeededRelinKeyGenToFile(
    const std::filesystem::path &path, const Key<P> &key,
    CKKSNoise noise = {P::α, 0},
    CKKSDenseBootstrapKeyDirectoryOptions options = {})
{
    if (!options.overwrite_existing && std::filesystem::exists(path)) return;

    auto relinkey = makeCKKSSeededRelinKey<P, LogQ>(key, noise);
    CKKSSavePortableBinaryAtomic(path, *relinkey);
}

template <class P, std::uint32_t LogQ>
inline bool CKKSRelinKeyGenNextMissingToFile(
    const std::filesystem::path &path, const Key<P> &key,
    CKKSNoise noise = {P::α, 0})
{
    if (std::filesystem::exists(path)) return false;
    CKKSRelinKeyGenToFile<P, LogQ>(path, key, noise);
    return true;
}

template <class P, std::uint32_t LogQ>
inline bool CKKSSeededRelinKeyGenNextMissingToFile(
    const std::filesystem::path &path, const Key<P> &key,
    CKKSNoise noise = {P::α, 0})
{
    if (std::filesystem::exists(path)) return false;
    CKKSSeededRelinKeyGenToFile<P, LogQ>(path, key, noise);
    return true;
}

template <class P, std::uint32_t LogQ>
inline void CKKSLoadRelinKeyFromFile(CKKSRelinKey<P, LogQ> &relinkey,
                                     const std::filesystem::path &path)
{
    CKKSRequireKeyFileReadable(path, "relinearization");
    CKKSLoadPortableBinary(relinkey, path);
}

template <class P, std::uint32_t LogQ>
inline void CKKSLoadSeededRelinKeyFromFile(
    CKKSSeededRelinKey<P, LogQ> &relinkey,
    const std::filesystem::path &path)
{
    CKKSRequireKeyFileReadable(path, "seeded relinearization");
    CKKSLoadPortableBinary(relinkey, path);
}

template <class P, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
inline void CKKSMultWithRelinKeyFile(
    CKKSMultResult<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta> &res,
    const CKKSCiphertext<P, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<P, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &relin_key_file)
{
    using Traits =
        CKKSMultTraits<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta>;
    auto relinkey = std::make_unique<CKKSRelinKey<P, Traits::log_q>>();
    CKKSLoadRelinKeyFromFile<P, Traits::log_q>(*relinkey, relin_key_file);
    CKKSMult<P>(res, lhs, rhs, *relinkey);
}

template <class P, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta>
inline void CKKSMultWithSeededRelinKeyFile(
    CKKSMultResult<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta> &res,
    const CKKSCiphertext<P, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<P, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &relin_key_file)
{
    using Traits =
        CKKSMultTraits<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta>;
    auto relinkey = std::make_unique<CKKSSeededRelinKey<P, Traits::log_q>>();
    CKKSLoadSeededRelinKeyFromFile<P, Traits::log_q>(*relinkey,
                                                     relin_key_file);
    CKKSMult<P>(res, lhs, rhs, *relinkey);
}

template <class Schedule>
inline void CKKSDenseBootstrapEncapsulationKeyGenToFile(
    const std::filesystem::path &path,
    const Key<typename Schedule::Param> &input_output_key,
    const Key<typename Schedule::Param> &bootstrap_key,
    CKKSNoise noise = {Schedule::Param::α, 0},
    CKKSDenseBootstrapKeyDirectoryOptions options = {})
{
    if (!options.overwrite_existing && std::filesystem::exists(path)) return;

    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapEncapsulationKeyGen<Schedule>(
        *encapsulation_key, input_output_key, bootstrap_key, noise);
    CKKSSavePortableBinaryAtomic(path, *encapsulation_key);
}

template <class Schedule>
inline void CKKSDenseBootstrapSeededEncapsulationKeyGenToFile(
    const std::filesystem::path &path,
    const Key<typename Schedule::Param> &input_output_key,
    const Key<typename Schedule::Param> &bootstrap_key,
    CKKSNoise noise = {Schedule::Param::α, 0},
    CKKSDenseBootstrapKeyDirectoryOptions options = {})
{
    if (!options.overwrite_existing && std::filesystem::exists(path)) return;

    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapSeededEncapsulationKeyGen<Schedule>(
        *encapsulation_key, input_output_key, bootstrap_key, noise);
    CKKSSavePortableBinaryAtomic(path, *encapsulation_key);
}

template <class Schedule>
inline bool CKKSDenseBootstrapEncapsulationKeyGenNextMissingToFile(
    const std::filesystem::path &path,
    const Key<typename Schedule::Param> &input_output_key,
    const Key<typename Schedule::Param> &bootstrap_key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    if (std::filesystem::exists(path)) return false;
    CKKSDenseBootstrapEncapsulationKeyGenToFile<Schedule>(
        path, input_output_key, bootstrap_key, noise);
    return true;
}

template <class Schedule>
inline bool CKKSDenseBootstrapSeededEncapsulationKeyGenNextMissingToFile(
    const std::filesystem::path &path,
    const Key<typename Schedule::Param> &input_output_key,
    const Key<typename Schedule::Param> &bootstrap_key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    if (std::filesystem::exists(path)) return false;
    CKKSDenseBootstrapSeededEncapsulationKeyGenToFile<Schedule>(
        path, input_output_key, bootstrap_key, noise);
    return true;
}

template <class Schedule>
inline void CKKSDenseBootstrapLoadEncapsulationKeyFromFile(
    CKKSDenseBootstrapEncapsulationKey<Schedule> &encapsulation_key,
    const std::filesystem::path &path)
{
    CKKSRequireKeyFileReadable(path, "encapsulation");
    CKKSLoadPortableBinary(encapsulation_key, path);
}

template <class Schedule>
inline void CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile(
    CKKSDenseBootstrapSeededEncapsulationKey<Schedule> &encapsulation_key,
    const std::filesystem::path &path)
{
    CKKSRequireKeyFileReadable(path, "seeded encapsulation");
    CKKSLoadPortableBinary(encapsulation_key, path);
}

template <class Schedule>
inline std::vector<std::filesystem::path> CKKSDenseBootstrapKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    using EvalModTraits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    std::vector<std::filesystem::path> paths;
    paths.reserve(6 + Schedule::coeff_to_slot_level_count +
                  Schedule::slot_to_coeff_level_count +
                  EvalModTraits::PolynomialTraits::relin_depth +
                  Schedule::evalmod_double_angle +
                  EvalModTraits::InverseTraits::power_depth);
    paths.push_back(CKKSDenseBootstrapKeyDirectoryManifestFile(root));
    paths.push_back(ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                             "linear_plan"));
    paths.push_back(ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                             "rotation_usage"));
    paths.push_back(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "evalmod_polynomial"));
    ckks_detail::CKKSDenseBootstrapCoeffToSlotKeyFilePathsImpl<0, Schedule>(
        paths, root);
    paths.push_back(ckks_detail::CKKSDenseBootstrapNamedPath(
        root, "packed_conjugate_galois"));
    ckks_detail::CKKSDenseBootstrapPolynomialRelinKeyFilePathsImpl<0,
                                                                    Schedule>(
        paths, root);
    ckks_detail::CKKSDenseBootstrapDoubleAngleRelinKeyFilePathsImpl<0,
                                                                    Schedule>(
        paths, root);
    ckks_detail::CKKSDenseBootstrapInverseRelinKeyFilePathsImpl<0, Schedule>(
        paths, root);
    ckks_detail::CKKSDenseBootstrapSlotToCoeffKeyFilePathsImpl<0, Schedule>(
        paths, root);
    return paths;
}

template <class Schedule>
inline std::vector<std::filesystem::path>
CKKSDenseBootstrapMissingKeyDirectoryFiles(const std::filesystem::path &root)
{
    std::vector<std::filesystem::path> missing;
    for (const std::filesystem::path &path :
         CKKSDenseBootstrapKeyDirectoryFiles<Schedule>(root)) {
        if (!std::filesystem::exists(path)) missing.push_back(path);
    }
    return missing;
}

template <class Schedule>
inline bool CKKSDenseBootstrapKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    return CKKSDenseBootstrapMissingKeyDirectoryFiles<Schedule>(root).empty();
}

namespace ckks_detail {

inline void CKKSDenseBootstrapRequireCompleteKeyDirectoryImpl(
    const std::filesystem::path &root,
    const std::vector<std::filesystem::path> &missing, const char *label)
{
    if (missing.empty()) return;
    std::string message = std::string("incomplete ") + label +
                          " CKKS bootstrap key directory: " + root.string() +
                          " missing=" + std::to_string(missing.size());
    if (!missing.front().empty())
        message += " first_missing=" + missing.front().string();
    throw std::runtime_error(message);
}

}  // namespace ckks_detail

template <class Schedule>
inline void CKKSDenseBootstrapRequireKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    ckks_detail::CKKSDenseBootstrapRequireCompleteKeyDirectoryImpl(
        root, CKKSDenseBootstrapMissingKeyDirectoryFiles<Schedule>(root),
        "sparse");
}

template <class Schedule>
inline std::vector<std::filesystem::path>
CKKSDenseBootstrapHybridGiantKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    return CKKSDenseBootstrapKeyDirectoryFiles<Schedule>(root);
}

template <class Schedule>
inline std::vector<std::filesystem::path>
CKKSDenseBootstrapHybridGiantMissingKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    return CKKSDenseBootstrapMissingKeyDirectoryFiles<Schedule>(root);
}

template <class Schedule>
inline bool CKKSDenseBootstrapHybridGiantKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    return CKKSDenseBootstrapHybridGiantMissingKeyDirectoryFiles<Schedule>(
               root)
        .empty();
}

template <class Schedule>
inline void CKKSDenseBootstrapHybridGiantRequireKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    ckks_detail::CKKSDenseBootstrapRequireCompleteKeyDirectoryImpl(
        root,
        CKKSDenseBootstrapHybridGiantMissingKeyDirectoryFiles<Schedule>(root),
        "hybrid");
}

template <class Schedule>
inline std::vector<std::filesystem::path>
CKKSDenseBootstrapSeededHybridGiantKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    return CKKSDenseBootstrapHybridGiantKeyDirectoryFiles<Schedule>(root);
}

template <class Schedule>
inline std::vector<std::filesystem::path>
CKKSDenseBootstrapSeededHybridGiantMissingKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    return CKKSDenseBootstrapHybridGiantMissingKeyDirectoryFiles<Schedule>(
        root);
}

template <class Schedule>
inline bool CKKSDenseBootstrapSeededHybridGiantKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    return CKKSDenseBootstrapSeededHybridGiantMissingKeyDirectoryFiles<
               Schedule>(root)
        .empty();
}

template <class Schedule>
inline void CKKSDenseBootstrapSeededHybridGiantRequireKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    ckks_detail::CKKSDenseBootstrapRequireCompleteKeyDirectoryImpl(
        root,
        CKKSDenseBootstrapSeededHybridGiantMissingKeyDirectoryFiles<Schedule>(
            root),
        "seeded hybrid");
}

template <class Schedule>
inline std::vector<std::filesystem::path>
CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    using P = typename Schedule::Param;
    using EvalModTraits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    std::vector<std::filesystem::path> paths;
    paths.reserve(8 + Schedule::coeff_to_slot_level_count *
                          (P::nbit + 1) +
                  Schedule::slot_to_coeff_level_count *
                      (P::nbit + 1) +
                  EvalModTraits::PolynomialTraits::relin_depth +
                  Schedule::evalmod_double_angle +
                  EvalModTraits::InverseTraits::power_depth);
    paths.push_back(CKKSDenseBootstrapKeyDirectoryManifestFile(root));
    paths.push_back(ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                             "linear_plan"));
    paths.push_back(ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                             "rotation_usage"));
    paths.push_back(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "evalmod_polynomial"));

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        rotation_usage, linear_plan);

    ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyFilePathsImpl<
            0, Schedule>(paths, root, rotation_usage);
    paths.push_back(ckks_detail::CKKSDenseBootstrapNamedPath(
        root, "packed_conjugate_galois"));
    ckks_detail::CKKSDenseBootstrapPolynomialRelinKeyFilePathsImpl<0,
                                                                    Schedule>(
        paths, root);
    ckks_detail::CKKSDenseBootstrapDoubleAngleRelinKeyFilePathsImpl<0,
                                                                    Schedule>(
        paths, root);
    ckks_detail::CKKSDenseBootstrapInverseRelinKeyFilePathsImpl<0, Schedule>(
        paths, root);
    ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyFilePathsImpl<
            0, Schedule>(paths, root, rotation_usage);
    return paths;
}

template <class Schedule>
inline std::vector<std::filesystem::path>
CKKSDenseBootstrapSeededHybridGiantStreamedMissingKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    std::vector<std::filesystem::path> missing;
    for (const std::filesystem::path &path :
         CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryFiles<Schedule>(
             root)) {
        if (!std::filesystem::exists(path)) missing.push_back(path);
    }
    return missing;
}

template <class Schedule>
inline bool CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    return CKKSDenseBootstrapSeededHybridGiantStreamedMissingKeyDirectoryFiles<
               Schedule>(root)
        .empty();
}

template <class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantStreamedRequireKeyDirectoryComplete(
    const std::filesystem::path &root)
{
    ckks_detail::CKKSDenseBootstrapRequireCompleteKeyDirectoryImpl(
        root,
        CKKSDenseBootstrapSeededHybridGiantStreamedMissingKeyDirectoryFiles<
            Schedule>(root),
        "streamed seeded hybrid");
}

template <class Schedule>
inline std::optional<CKKSDenseBootstrapKeyFileByteEstimate>
CKKSDenseBootstrapSeededHybridGiantStreamedNextMissingGeneratedKeyFileEstimate(
    const std::filesystem::path &root)
{
    using P = typename Schedule::Param;
    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(rotation_usage,
                                                                 linear_plan);

    if (auto next = ckks_detail::
            CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotNextMissingFileByteEstimateImpl<
                0, Schedule>(root, rotation_usage))
        return next;

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (!std::filesystem::exists(packed_path)) {
        return CKKSDenseBootstrapKeyFileByteEstimate{
            packed_path,
            CKKSDenseBootstrapHybridGiantPackedConjugateKeySwitchRowCount<
                Schedule>(rotation_usage) *
                CKKSSeededKeySwitchRowByteSize<P>()};
    }

    if (auto next = ckks_detail::
            CKKSDenseBootstrapSeededPolynomialRelinNextMissingFileByteEstimateImpl<
                0, Schedule>(root))
        return next;
    if (auto next = ckks_detail::
            CKKSDenseBootstrapSeededDoubleAngleRelinNextMissingFileByteEstimateImpl<
                0, Schedule>(root))
        return next;
    if (auto next = ckks_detail::
            CKKSDenseBootstrapSeededInverseRelinNextMissingFileByteEstimateImpl<
                0, Schedule>(root))
        return next;
    return ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffNextMissingFileByteEstimateImpl<
            0, Schedule>(root, rotation_usage);
}

template <class Schedule>
inline CKKSDenseBootstrapKeyDirectoryManifest
CKKSDenseBootstrapBuildKeyDirectoryManifest(const std::filesystem::path &root)
{
    using P = typename Schedule::Param;
    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    CKKSDenseBootstrapRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapRotationKeyUsage<Schedule>(rotation_usage,
                                                      linear_plan);

    CKKSDenseBootstrapKeyDirectoryManifest manifest;
    manifest.key_format =
        CKKSDenseBootstrapKeyDirectoryManifest::sparse_format;
    manifest.n = P::n;
    manifest.nbit = P::nbit;
    manifest.torus_bits = std::numeric_limits<typename P::T>::digits;
    manifest.log_delta = Schedule::log_delta;
    manifest.log_message_ratio = Schedule::log_message_ratio;
    manifest.input_log_q = Schedule::input_log_q;
    manifest.boot_log_q = Schedule::boot_log_q;
    manifest.output_log_q = Schedule::output_log_q;
    manifest.linear_plain_log_delta = Schedule::linear_plain_log_delta;
    manifest.linear_fuse_radix = Schedule::linear_fuse_radix;
    manifest.coeff_to_slot_plain_log_delta =
        Schedule::coeff_to_slot_plain_log_delta;
    manifest.component_split_plain_log_delta =
        Schedule::component_split_plain_log_delta;
    manifest.slot_to_coeff_plain_log_delta =
        Schedule::slot_to_coeff_plain_log_delta;
    manifest.coeff_to_slot_fuse_radix = Schedule::coeff_to_slot_fuse_radix;
    manifest.slot_to_coeff_fuse_radix = Schedule::slot_to_coeff_fuse_radix;
    manifest.linear_bsgs_step = Schedule::linear_bsgs_step;
    manifest.coeff_to_slot_level_count = Schedule::coeff_to_slot_level_count;
    manifest.slot_to_coeff_level_count = Schedule::slot_to_coeff_level_count;
    manifest.evalmod_degree = Schedule::evalmod_degree;
    manifest.evalmod_k = Schedule::evalmod_k;
    manifest.evalmod_double_angle = Schedule::evalmod_double_angle;
    manifest.evalmod_inv_degree = Schedule::evalmod_inv_degree;
    manifest.evalmod_log_scale = Schedule::evalmod_log_scale;
    manifest.evalmod_depth = Schedule::evalmod_depth;
    manifest.modraise_mask_bound = Schedule::modraise_mask_bound;
    manifest.expected_file_count =
        CKKSDenseBootstrapKeyDirectoryFiles<Schedule>(root).size();
    manifest.sparse_key_rows =
        CKKSDenseBootstrapSparseKeySwitchRowCount<Schedule>(rotation_usage);
    manifest.streamed_peak_key_rows =
        CKKSDenseBootstrapStreamedKeySwitchPeakRowCount<Schedule>(
            rotation_usage);
    manifest.hybrid_key_rows = 0;
    manifest.hybrid_streamed_peak_key_rows = 0;
    manifest.full_key_rows = CKKSDenseBootstrapFullKeySwitchRowCount<Schedule>();
    return manifest;
}

template <class Schedule>
inline CKKSDenseBootstrapKeyDirectoryManifest
CKKSDenseBootstrapBuildKeyDirectoryManifest()
{
    return CKKSDenseBootstrapBuildKeyDirectoryManifest<Schedule>(
        std::filesystem::path{});
}

template <class Schedule>
inline CKKSDenseBootstrapKeyDirectoryManifest
CKKSDenseBootstrapBuildHybridGiantKeyDirectoryManifest(
    const std::filesystem::path &root)
{
    CKKSDenseBootstrapKeyDirectoryManifest manifest =
        CKKSDenseBootstrapBuildKeyDirectoryManifest<Schedule>(root);
    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        rotation_usage, linear_plan);
    if constexpr (Schedule::hybrid_giant_direct_popcount_threshold == 4) {
        manifest.key_format =
            CKKSDenseBootstrapKeyDirectoryManifest::hybrid_giant_format;
    }
    else {
        manifest.key_format =
            CKKSDenseBootstrapKeyDirectoryManifest::
                hybrid_giant_tuned_format_base +
            static_cast<std::uint32_t>(
                Schedule::hybrid_giant_direct_popcount_threshold);
    }
    manifest.sparse_key_rows = 0;
    manifest.streamed_peak_key_rows = 0;
    manifest.hybrid_key_rows =
        CKKSDenseBootstrapHybridGiantKeySwitchRowCount<Schedule>(
            rotation_usage);
    manifest.hybrid_streamed_peak_key_rows =
        CKKSDenseBootstrapHybridGiantStreamedKeySwitchPeakRowCount<Schedule>(
            rotation_usage);
    return manifest;
}

template <class Schedule>
inline CKKSDenseBootstrapKeyDirectoryManifest
CKKSDenseBootstrapBuildHybridGiantKeyDirectoryManifest()
{
    return CKKSDenseBootstrapBuildHybridGiantKeyDirectoryManifest<Schedule>(
        std::filesystem::path{});
}

template <class Schedule>
inline CKKSDenseBootstrapKeyDirectoryManifest
CKKSDenseBootstrapBuildSeededHybridGiantKeyDirectoryManifest(
    const std::filesystem::path &root)
{
    CKKSDenseBootstrapKeyDirectoryManifest manifest =
        CKKSDenseBootstrapBuildHybridGiantKeyDirectoryManifest<Schedule>(root);
    if constexpr (Schedule::hybrid_giant_direct_popcount_threshold == 4) {
        manifest.key_format =
            CKKSDenseBootstrapKeyDirectoryManifest::seeded_hybrid_giant_format;
    }
    else {
        manifest.key_format =
            CKKSDenseBootstrapKeyDirectoryManifest::
                seeded_hybrid_giant_tuned_format_base +
            static_cast<std::uint32_t>(
                Schedule::hybrid_giant_direct_popcount_threshold);
    }
    return manifest;
}

template <class Schedule>
inline CKKSDenseBootstrapKeyDirectoryManifest
CKKSDenseBootstrapBuildSeededHybridGiantKeyDirectoryManifest()
{
    return CKKSDenseBootstrapBuildSeededHybridGiantKeyDirectoryManifest<
        Schedule>(std::filesystem::path{});
}

template <class Schedule>
inline CKKSDenseBootstrapKeyDirectoryManifest
CKKSDenseBootstrapBuildSeededHybridGiantStreamedKeyDirectoryManifest(
    const std::filesystem::path &root)
{
    CKKSDenseBootstrapKeyDirectoryManifest manifest =
        CKKSDenseBootstrapBuildHybridGiantKeyDirectoryManifest<Schedule>(root);
    if constexpr (Schedule::hybrid_giant_direct_popcount_threshold == 4) {
        manifest.key_format =
            CKKSDenseBootstrapKeyDirectoryManifest::
                seeded_hybrid_giant_streamed_format;
    }
    else {
        manifest.key_format =
            CKKSDenseBootstrapKeyDirectoryManifest::
                seeded_hybrid_giant_streamed_tuned_format_base +
            static_cast<std::uint32_t>(
                Schedule::hybrid_giant_direct_popcount_threshold);
    }
    manifest.expected_file_count =
        CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryFiles<Schedule>(
            root)
            .size();
    return manifest;
}

template <class Schedule>
inline CKKSDenseBootstrapKeyDirectoryManifest
CKKSDenseBootstrapBuildSeededHybridGiantStreamedKeyDirectoryManifest()
{
    return CKKSDenseBootstrapBuildSeededHybridGiantStreamedKeyDirectoryManifest<
        Schedule>(std::filesystem::path{});
}

template <class Schedule>
inline bool CKKSDenseBootstrapKeyDirectoryManifestMatches(
    const CKKSDenseBootstrapKeyDirectoryManifest &manifest)
{
    return manifest == CKKSDenseBootstrapBuildKeyDirectoryManifest<Schedule>();
}

template <class Schedule>
inline bool CKKSDenseBootstrapHybridGiantKeyDirectoryManifestMatches(
    const CKKSDenseBootstrapKeyDirectoryManifest &manifest)
{
    return manifest ==
           CKKSDenseBootstrapBuildHybridGiantKeyDirectoryManifest<Schedule>();
}

template <class Schedule>
inline bool CKKSDenseBootstrapSeededHybridGiantKeyDirectoryManifestMatches(
    const CKKSDenseBootstrapKeyDirectoryManifest &manifest)
{
    return manifest ==
           CKKSDenseBootstrapBuildSeededHybridGiantKeyDirectoryManifest<
               Schedule>();
}

template <class Schedule>
inline bool
CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryManifestMatches(
    const CKKSDenseBootstrapKeyDirectoryManifest &manifest)
{
    return manifest ==
           CKKSDenseBootstrapBuildSeededHybridGiantStreamedKeyDirectoryManifest<
               Schedule>();
}

template <class Schedule>
inline bool CKKSDenseBootstrapKeyDirectoryManifestMatches(
    const std::filesystem::path &root)
{
    CKKSDenseBootstrapKeyDirectoryManifest manifest;
    CKKSLoadPortableBinary(manifest,
                           CKKSDenseBootstrapKeyDirectoryManifestFile(root));
    return CKKSDenseBootstrapKeyDirectoryManifestMatches<Schedule>(manifest);
}

template <class Schedule>
inline bool CKKSDenseBootstrapHybridGiantKeyDirectoryManifestMatches(
    const std::filesystem::path &root)
{
    CKKSDenseBootstrapKeyDirectoryManifest manifest;
    CKKSLoadPortableBinary(manifest,
                           CKKSDenseBootstrapKeyDirectoryManifestFile(root));
    return CKKSDenseBootstrapHybridGiantKeyDirectoryManifestMatches<Schedule>(
        manifest);
}

template <class Schedule>
inline bool CKKSDenseBootstrapSeededHybridGiantKeyDirectoryManifestMatches(
    const std::filesystem::path &root)
{
    CKKSDenseBootstrapKeyDirectoryManifest manifest;
    CKKSLoadPortableBinary(manifest,
                           CKKSDenseBootstrapKeyDirectoryManifestFile(root));
    return CKKSDenseBootstrapSeededHybridGiantKeyDirectoryManifestMatches<
        Schedule>(manifest);
}

template <class Schedule>
inline bool
CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryManifestMatches(
    const std::filesystem::path &root)
{
    CKKSDenseBootstrapKeyDirectoryManifest manifest;
    CKKSLoadPortableBinary(manifest,
                           CKKSDenseBootstrapKeyDirectoryManifestFile(root));
    return CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryManifestMatches<
        Schedule>(manifest);
}

template <class Schedule>
inline void CKKSDenseBootstrapCheckKeyDirectoryManifestForWrite(
    const std::filesystem::path &root, bool overwrite_existing)
{
    const std::filesystem::path manifest_path =
        CKKSDenseBootstrapKeyDirectoryManifestFile(root);
    if (!std::filesystem::exists(manifest_path)) {
        if (!overwrite_existing) {
            for (const std::filesystem::path &path :
                 CKKSDenseBootstrapKeyDirectoryFiles<Schedule>(root)) {
                if (path != manifest_path && std::filesystem::exists(path)) {
                    throw std::runtime_error(
                        "CKKS bootstrap key directory has key material but no "
                        "manifest; regenerate with overwrite enabled");
                }
            }
        }
        return;
    }

    if (!CKKSDenseBootstrapKeyDirectoryManifestMatches<Schedule>(root) &&
        !overwrite_existing) {
        throw std::runtime_error(
            "CKKS bootstrap key directory manifest does not match schedule; "
            "regenerate with overwrite enabled");
    }
}

template <class Schedule>
inline void CKKSDenseBootstrapHybridGiantCheckKeyDirectoryManifestForWrite(
    const std::filesystem::path &root, bool overwrite_existing)
{
    const std::filesystem::path manifest_path =
        CKKSDenseBootstrapKeyDirectoryManifestFile(root);
    if (!std::filesystem::exists(manifest_path)) {
        if (!overwrite_existing) {
            for (const std::filesystem::path &path :
                 CKKSDenseBootstrapHybridGiantKeyDirectoryFiles<Schedule>(
                     root)) {
                if (path != manifest_path && std::filesystem::exists(path)) {
                    throw std::runtime_error(
                        "CKKS hybrid bootstrap key directory has key material "
                        "but no manifest; regenerate with overwrite enabled");
                }
            }
        }
        return;
    }

    if (!CKKSDenseBootstrapHybridGiantKeyDirectoryManifestMatches<Schedule>(
            root) &&
        !overwrite_existing) {
        throw std::runtime_error(
            "CKKS hybrid bootstrap key directory manifest does not match "
            "schedule; regenerate with overwrite enabled");
    }
}

template <class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantCheckKeyDirectoryManifestForWrite(
    const std::filesystem::path &root, bool overwrite_existing)
{
    const std::filesystem::path manifest_path =
        CKKSDenseBootstrapKeyDirectoryManifestFile(root);
    if (!std::filesystem::exists(manifest_path)) {
        if (!overwrite_existing) {
            for (const std::filesystem::path &path :
                 CKKSDenseBootstrapHybridGiantKeyDirectoryFiles<Schedule>(
                     root)) {
                if (path != manifest_path && std::filesystem::exists(path)) {
                    throw std::runtime_error(
                        "CKKS seeded hybrid bootstrap key directory has key "
                        "material but no manifest; regenerate with overwrite "
                        "enabled");
                }
            }
        }
        return;
    }

    if (!CKKSDenseBootstrapSeededHybridGiantKeyDirectoryManifestMatches<
            Schedule>(root) &&
        !overwrite_existing) {
        throw std::runtime_error(
            "CKKS seeded hybrid bootstrap key directory manifest does not "
            "match schedule; regenerate with overwrite enabled");
    }
}

template <class Schedule>
inline void
CKKSDenseBootstrapSeededHybridGiantStreamedCheckKeyDirectoryManifestForWrite(
    const std::filesystem::path &root, bool overwrite_existing)
{
    const std::filesystem::path manifest_path =
        CKKSDenseBootstrapKeyDirectoryManifestFile(root);
    if (!std::filesystem::exists(manifest_path)) {
        if (!overwrite_existing) {
            for (const std::filesystem::path &path :
                 CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryFiles<
                     Schedule>(root)) {
                if (path != manifest_path && std::filesystem::exists(path)) {
                    throw std::runtime_error(
                        "CKKS streamed seeded hybrid bootstrap key directory "
                        "has key material but no manifest; regenerate with "
                        "overwrite enabled");
                }
            }
        }
        return;
    }

    if (!CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryManifestMatches<
            Schedule>(root) &&
        !overwrite_existing) {
        throw std::runtime_error(
            "CKKS streamed seeded hybrid bootstrap key directory manifest "
            "does not match schedule; regenerate with overwrite enabled");
    }
}

template <class Schedule>
inline void CKKSDenseBootstrapRemoveKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    for (const std::filesystem::path &path :
         CKKSDenseBootstrapKeyDirectoryFiles<Schedule>(root)) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::filesystem::path tmp = path;
        tmp += ".tmp";
        std::filesystem::remove(tmp, ec);
    }
}

template <class Schedule>
inline void CKKSDenseBootstrapSeededHybridGiantStreamedRemoveKeyDirectoryFiles(
    const std::filesystem::path &root)
{
    for (const std::filesystem::path &path :
         CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryFiles<Schedule>(
             root)) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::filesystem::path tmp = path;
        tmp += ".tmp";
        std::filesystem::remove(tmp, ec);
    }
}

template <class Schedule>
inline CKKSDenseBootstrapRotationKeyUsage<Schedule>
CKKSDenseBootstrapWriteKeyDirectoryMetadata(const std::filesystem::path &root,
                                            bool overwrite_existing)
{
    std::filesystem::create_directories(root);
    if (overwrite_existing)
        CKKSDenseBootstrapRemoveKeyDirectoryFiles<Schedule>(root);
    else
        CKKSDenseBootstrapCheckKeyDirectoryManifestForWrite<Schedule>(
            root, overwrite_existing);

    const CKKSDenseBootstrapKeyDirectoryManifest manifest =
        CKKSDenseBootstrapBuildKeyDirectoryManifest<Schedule>(root);
    CKKSSavePortableBinaryAtomic(
        CKKSDenseBootstrapKeyDirectoryManifestFile(root), manifest);

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "linear_plan"),
        linear_plan);

    CKKSDenseBootstrapRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapRotationKeyUsage<Schedule>(rotation_usage,
                                                      linear_plan);
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "rotation_usage"),
        rotation_usage);

    const CKKSBoundedCosEvalModPolynomial evalmod_polynomial =
        CKKSBuildBoundedCosEvalModPolynomial<Schedule>();
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "evalmod_polynomial"),
        evalmod_polynomial);

    return rotation_usage;
}

template <class Schedule>
inline CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
CKKSDenseBootstrapHybridGiantWriteKeyDirectoryMetadata(
    const std::filesystem::path &root, bool overwrite_existing)
{
    std::filesystem::create_directories(root);
    if (overwrite_existing)
        CKKSDenseBootstrapRemoveKeyDirectoryFiles<Schedule>(root);
    else
        CKKSDenseBootstrapHybridGiantCheckKeyDirectoryManifestForWrite<
            Schedule>(root, overwrite_existing);

    const CKKSDenseBootstrapKeyDirectoryManifest manifest =
        CKKSDenseBootstrapBuildHybridGiantKeyDirectoryManifest<Schedule>(root);
    CKKSSavePortableBinaryAtomic(
        CKKSDenseBootstrapKeyDirectoryManifestFile(root), manifest);

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "linear_plan"),
        linear_plan);

    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        rotation_usage, linear_plan);
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "rotation_usage"),
        rotation_usage);

    const CKKSBoundedCosEvalModPolynomial evalmod_polynomial =
        CKKSBuildBoundedCosEvalModPolynomial<Schedule>();
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "evalmod_polynomial"),
        evalmod_polynomial);

    return rotation_usage;
}

template <class Schedule>
inline CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
CKKSDenseBootstrapSeededHybridGiantWriteKeyDirectoryMetadata(
    const std::filesystem::path &root, bool overwrite_existing)
{
    std::filesystem::create_directories(root);
    if (overwrite_existing)
        CKKSDenseBootstrapRemoveKeyDirectoryFiles<Schedule>(root);
    else
        CKKSDenseBootstrapSeededHybridGiantCheckKeyDirectoryManifestForWrite<
            Schedule>(root, overwrite_existing);

    const CKKSDenseBootstrapKeyDirectoryManifest manifest =
        CKKSDenseBootstrapBuildSeededHybridGiantKeyDirectoryManifest<Schedule>(
            root);
    CKKSSavePortableBinaryAtomic(
        CKKSDenseBootstrapKeyDirectoryManifestFile(root), manifest);

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "linear_plan"),
        linear_plan);

    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        rotation_usage, linear_plan);
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "rotation_usage"),
        rotation_usage);

    const CKKSBoundedCosEvalModPolynomial evalmod_polynomial =
        CKKSBuildBoundedCosEvalModPolynomial<Schedule>();
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "evalmod_polynomial"),
        evalmod_polynomial);

    return rotation_usage;
}

template <class Schedule>
inline CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
CKKSDenseBootstrapSeededHybridGiantStreamedWriteKeyDirectoryMetadata(
    const std::filesystem::path &root, bool overwrite_existing)
{
    std::filesystem::create_directories(root);
    if (overwrite_existing)
        CKKSDenseBootstrapSeededHybridGiantStreamedRemoveKeyDirectoryFiles<
            Schedule>(root);
    else
        CKKSDenseBootstrapSeededHybridGiantStreamedCheckKeyDirectoryManifestForWrite<
            Schedule>(root, overwrite_existing);

    const CKKSDenseBootstrapKeyDirectoryManifest manifest =
        CKKSDenseBootstrapBuildSeededHybridGiantStreamedKeyDirectoryManifest<
            Schedule>(root);
    CKKSSavePortableBinaryAtomic(
        CKKSDenseBootstrapKeyDirectoryManifestFile(root), manifest);

    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan;
    CKKSBuildDenseBootstrapLinearPlan<Schedule>(linear_plan);
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "linear_plan"),
        linear_plan);

    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule> rotation_usage;
    CKKSBuildDenseBootstrapHybridGiantRotationKeyUsage<Schedule>(
        rotation_usage, linear_plan);
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "rotation_usage"),
        rotation_usage);

    const CKKSBoundedCosEvalModPolynomial evalmod_polynomial =
        CKKSBuildBoundedCosEvalModPolynomial<Schedule>();
    CKKSSavePortableBinaryAtomic(
        ckks_detail::CKKSDenseBootstrapNamedPath(root, "evalmod_polynomial"),
        evalmod_polynomial);

    return rotation_usage;
}

template <class Schedule>
inline bool CKKSDenseBootstrapKeyGenNextMissingToDirectory(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> rotation_usage =
        CKKSDenseBootstrapWriteKeyDirectoryMetadata<Schedule>(
            root, false);

    if (ckks_detail::CKKSDenseBootstrapCoeffToSlotKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, rotation_usage, key, noise))
        return true;

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (!std::filesystem::exists(packed_path)) {
        auto packed_conjugate =
            std::make_unique<CKKSDenseBootstrapPackedConjugateGaloisKey<
                Schedule>>();
        CKKSDenseBootstrapPackedConjugateGaloisKeyGen<Schedule>(
            *packed_conjugate, rotation_usage, key, noise);
        CKKSSavePortableBinaryAtomic(packed_path, *packed_conjugate);
        return true;
    }

    if (ckks_detail::CKKSDenseBootstrapPolynomialRelinKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, key, noise))
        return true;
    if (ckks_detail::CKKSDenseBootstrapDoubleAngleRelinKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, key, noise))
        return true;
    if (ckks_detail::CKKSDenseBootstrapInverseRelinKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, key, noise))
        return true;
    return ckks_detail::CKKSDenseBootstrapSlotToCoeffKeyGenNextMissingToDirectoryImpl<
        0, Schedule>(root, rotation_usage, key, noise);
}

template <class Schedule>
inline void CKKSDenseBootstrapKeyGenToDirectory(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0},
    CKKSDenseBootstrapKeyDirectoryOptions options = {})
{
    const CKKSDenseBootstrapRotationKeyUsage<Schedule> rotation_usage =
        CKKSDenseBootstrapWriteKeyDirectoryMetadata<Schedule>(
            root, options.overwrite_existing);

    ckks_detail::CKKSDenseBootstrapCoeffToSlotKeyGenToDirectoryImpl<0,
                                                                    Schedule>(
        root, rotation_usage, key, noise, options);

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (ckks_detail::CKKSDenseBootstrapShouldWriteKeyFile(packed_path,
                                                          options)) {
        auto packed_conjugate =
            std::make_unique<CKKSDenseBootstrapPackedConjugateGaloisKey<
                Schedule>>();
        CKKSDenseBootstrapPackedConjugateGaloisKeyGen<Schedule>(
            *packed_conjugate, rotation_usage, key, noise);
        CKKSSavePortableBinaryAtomic(packed_path, *packed_conjugate);
    }

    ckks_detail::CKKSDenseBootstrapPolynomialRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);
    ckks_detail::CKKSDenseBootstrapDoubleAngleRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);
    ckks_detail::CKKSDenseBootstrapInverseRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);

    ckks_detail::CKKSDenseBootstrapSlotToCoeffKeyGenToDirectoryImpl<0,
                                                                    Schedule>(
        root, rotation_usage, key, noise, options);
}

template <class Schedule>
inline bool CKKSDenseBootstrapHybridGiantKeyGenNextMissingToDirectory(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        rotation_usage =
            CKKSDenseBootstrapHybridGiantWriteKeyDirectoryMetadata<Schedule>(
                root, false);

    if (ckks_detail::
            CKKSDenseBootstrapHybridGiantCoeffToSlotKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, rotation_usage, key, noise))
        return true;

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (!std::filesystem::exists(packed_path)) {
        auto packed_conjugate =
            std::make_unique<CKKSDenseBootstrapPackedConjugateGaloisKey<
                Schedule>>();
        CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
            *packed_conjugate, key, rotation_usage.packed_conjugate, noise);
        CKKSSavePortableBinaryAtomic(packed_path, *packed_conjugate);
        return true;
    }

    if (ckks_detail::CKKSDenseBootstrapPolynomialRelinKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, key, noise))
        return true;
    if (ckks_detail::CKKSDenseBootstrapDoubleAngleRelinKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, key, noise))
        return true;
    if (ckks_detail::CKKSDenseBootstrapInverseRelinKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, key, noise))
        return true;
    return ckks_detail::
        CKKSDenseBootstrapHybridGiantSlotToCoeffKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, rotation_usage, key, noise);
}

template <class Schedule>
inline void CKKSDenseBootstrapHybridGiantKeyGenToDirectory(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0},
    CKKSDenseBootstrapKeyDirectoryOptions options = {})
{
    using P = typename Schedule::Param;
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        rotation_usage =
            CKKSDenseBootstrapHybridGiantWriteKeyDirectoryMetadata<Schedule>(
                root, options.overwrite_existing);

    ckks_detail::CKKSDenseBootstrapHybridGiantCoeffToSlotKeyGenToDirectoryImpl<
        0, Schedule>(root, rotation_usage, key, noise, options);

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (ckks_detail::CKKSDenseBootstrapShouldWriteKeyFile(packed_path,
                                                          options)) {
        auto packed_conjugate =
            std::make_unique<CKKSDenseBootstrapPackedConjugateGaloisKey<
                Schedule>>();
        CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
            *packed_conjugate, key, rotation_usage.packed_conjugate, noise);
        CKKSSavePortableBinaryAtomic(packed_path, *packed_conjugate);
    }

    ckks_detail::CKKSDenseBootstrapPolynomialRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);
    ckks_detail::CKKSDenseBootstrapDoubleAngleRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);
    ckks_detail::CKKSDenseBootstrapInverseRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);

    ckks_detail::CKKSDenseBootstrapHybridGiantSlotToCoeffKeyGenToDirectoryImpl<
        0, Schedule>(root, rotation_usage, key, noise, options);
}

template <class Schedule>
inline bool CKKSDenseBootstrapSeededHybridGiantKeyGenNextMissingToDirectory(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        rotation_usage =
            CKKSDenseBootstrapSeededHybridGiantWriteKeyDirectoryMetadata<
                Schedule>(root, false);

    if (ckks_detail::
            CKKSDenseBootstrapSeededHybridGiantCoeffToSlotKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, rotation_usage, key, noise))
        return true;

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (!std::filesystem::exists(packed_path)) {
        auto packed_conjugate =
            std::make_unique<CKKSDenseBootstrapSeededPackedConjugateGaloisKey<
                Schedule>>();
        CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
            *packed_conjugate, key, rotation_usage.packed_conjugate, noise);
        CKKSSavePortableBinaryAtomic(packed_path, *packed_conjugate);
        return true;
    }

    if (ckks_detail::
            CKKSDenseBootstrapSeededPolynomialRelinKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, key, noise))
        return true;
    if (ckks_detail::
            CKKSDenseBootstrapSeededDoubleAngleRelinKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, key, noise))
        return true;
    if (ckks_detail::
            CKKSDenseBootstrapSeededInverseRelinKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, key, noise))
        return true;
    return ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantSlotToCoeffKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, rotation_usage, key, noise);
}

template <class Schedule>
inline void CKKSDenseBootstrapSeededHybridGiantKeyGenToDirectory(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0},
    CKKSDenseBootstrapKeyDirectoryOptions options = {})
{
    using P = typename Schedule::Param;
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        rotation_usage =
            CKKSDenseBootstrapSeededHybridGiantWriteKeyDirectoryMetadata<
                Schedule>(root, options.overwrite_existing);

    ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantCoeffToSlotKeyGenToDirectoryImpl<
            0, Schedule>(root, rotation_usage, key, noise, options);

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (ckks_detail::CKKSDenseBootstrapShouldWriteKeyFile(packed_path,
                                                          options)) {
        auto packed_conjugate =
            std::make_unique<CKKSDenseBootstrapSeededPackedConjugateGaloisKey<
                Schedule>>();
        CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
            *packed_conjugate, key, rotation_usage.packed_conjugate, noise);
        CKKSSavePortableBinaryAtomic(packed_path, *packed_conjugate);
    }

    ckks_detail::CKKSDenseBootstrapSeededPolynomialRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);
    ckks_detail::CKKSDenseBootstrapSeededDoubleAngleRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);
    ckks_detail::CKKSDenseBootstrapSeededInverseRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);

    ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantSlotToCoeffKeyGenToDirectoryImpl<
            0, Schedule>(root, rotation_usage, key, noise, options);
}

template <class Schedule>
inline bool
CKKSDenseBootstrapSeededHybridGiantStreamedKeyGenNextMissingToDirectory(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        rotation_usage =
            CKKSDenseBootstrapSeededHybridGiantStreamedWriteKeyDirectoryMetadata<
                Schedule>(root, false);

    if (ckks_detail::
            CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, rotation_usage, key, noise))
        return true;

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (!std::filesystem::exists(packed_path)) {
        auto packed_conjugate =
            std::make_unique<CKKSDenseBootstrapSeededPackedConjugateGaloisKey<
                Schedule>>();
        CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
            *packed_conjugate, key, rotation_usage.packed_conjugate, noise);
        CKKSSavePortableBinaryAtomic(packed_path, *packed_conjugate);
        return true;
    }

    if (ckks_detail::
            CKKSDenseBootstrapSeededPolynomialRelinKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, key, noise))
        return true;
    if (ckks_detail::
            CKKSDenseBootstrapSeededDoubleAngleRelinKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, key, noise))
        return true;
    if (ckks_detail::
            CKKSDenseBootstrapSeededInverseRelinKeyGenNextMissingToDirectoryImpl<
                0, Schedule>(root, key, noise))
        return true;
    return ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyGenNextMissingToDirectoryImpl<
            0, Schedule>(root, rotation_usage, key, noise);
}

template <class Schedule>
inline void CKKSDenseBootstrapSeededHybridGiantStreamedKeyGenToDirectory(
    const std::filesystem::path &root, const Key<typename Schedule::Param> &key,
    CKKSNoise noise = {Schedule::Param::α, 0},
    CKKSDenseBootstrapKeyDirectoryOptions options = {})
{
    using P = typename Schedule::Param;
    const CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        rotation_usage =
            CKKSDenseBootstrapSeededHybridGiantStreamedWriteKeyDirectoryMetadata<
                Schedule>(root, options.overwrite_existing);

    ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotKeyGenToDirectoryImpl<
            0, Schedule>(root, rotation_usage, key, noise, options);

    const std::filesystem::path packed_path =
        ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                 "packed_conjugate_galois");
    if (ckks_detail::CKKSDenseBootstrapShouldWriteKeyFile(packed_path,
                                                          options)) {
        auto packed_conjugate =
            std::make_unique<CKKSDenseBootstrapSeededPackedConjugateGaloisKey<
                Schedule>>();
        CKKSSparseGaloisKeyGen<P, Schedule::after_coeff_to_slot_log_q>(
            *packed_conjugate, key, rotation_usage.packed_conjugate, noise);
        CKKSSavePortableBinaryAtomic(packed_path, *packed_conjugate);
    }

    ckks_detail::CKKSDenseBootstrapSeededPolynomialRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);
    ckks_detail::CKKSDenseBootstrapSeededDoubleAngleRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);
    ckks_detail::CKKSDenseBootstrapSeededInverseRelinKeyGenToDirectoryImpl<
        0, Schedule>(root, key, noise, options);

    ckks_detail::
        CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffKeyGenToDirectoryImpl<
            0, Schedule>(root, rotation_usage, key, noise, options);
}

template <class Schedule>
struct CKKSDenseBootstrapFilesystemKeyProvider {
    using EvalModTraits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    using CoeffToSlotCache =
        typename ckks_detail::CKKSDenseBootstrapCoeffToSlotCacheTuple<
            Schedule,
            std::make_index_sequence<Schedule::coeff_to_slot_level_count + 1>>::
            type;
    using SlotToCoeffCache =
        typename ckks_detail::CKKSDenseBootstrapSlotToCoeffCacheTuple<
            Schedule,
            std::make_index_sequence<Schedule::slot_to_coeff_level_count + 1>>::
            type;
    using PolynomialRelinCache =
        typename ckks_detail::CKKSDenseBootstrapPolynomialRelinCacheTuple<
            Schedule,
            std::make_index_sequence<EvalModTraits::PolynomialTraits::
                                         relin_depth>>::type;
    using DoubleAngleRelinCache =
        typename ckks_detail::CKKSDenseBootstrapDoubleAngleRelinCacheTuple<
            Schedule,
            std::make_index_sequence<Schedule::evalmod_double_angle>>::type;
    using InverseRelinCache =
        typename ckks_detail::CKKSDenseBootstrapInverseRelinCacheTuple<
            Schedule,
            std::make_index_sequence<EvalModTraits::InverseTraits::
                                         power_depth>>::type;

    std::filesystem::path root;
    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan_cache{};
    CKKSBoundedCosEvalModPolynomial evalmod_polynomial_cache{};
    mutable CoeffToSlotCache coeff_to_slot_cache{};
    mutable std::unique_ptr<CKKSDenseBootstrapPackedConjugateGaloisKey<Schedule>>
        packed_conjugate_cache{};
    mutable PolynomialRelinCache polynomial_relin_cache{};
    mutable DoubleAngleRelinCache double_angle_relin_cache{};
    mutable InverseRelinCache inverse_relin_cache{};
    mutable SlotToCoeffCache slot_to_coeff_cache{};

    explicit CKKSDenseBootstrapFilesystemKeyProvider(
        std::filesystem::path root_)
        : root(std::move(root_))
    {
        try {
            if (!CKKSDenseBootstrapKeyDirectoryManifestMatches<Schedule>(
                    root)) {
                throw std::runtime_error(
                    "CKKS bootstrap key directory manifest does not match "
                    "schedule");
            }
        }
        catch (const std::exception &e) {
            throw std::runtime_error(
                std::string("invalid CKKS bootstrap key directory manifest: ") +
                e.what());
        }
        CKKSDenseBootstrapRequireKeyDirectoryComplete<Schedule>(root);
        CKKSLoadPortableBinary(
            linear_plan_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root, "linear_plan"));
        CKKSLoadPortableBinary(
            evalmod_polynomial_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                     "evalmod_polynomial"));
    }

    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan() const
    {
        return linear_plan_cache;
    }

    const CKKSBoundedCosEvalModPolynomial &evalmod_polynomial() const
    {
        return evalmod_polynomial_cache;
    }

    template <std::size_t I>
    const auto &coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        auto &entry = std::get<I>(coeff_to_slot_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseBootstrapCoeffToSlotGaloisKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "coeff_to_slot_galois", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        std::get<I>(coeff_to_slot_cache).reset();
    }

    const auto &packed_conjugate_galois() const
    {
        if (!packed_conjugate_cache) {
            packed_conjugate_cache =
                std::make_unique<CKKSDenseBootstrapPackedConjugateGaloisKey<
                    Schedule>>();
            CKKSLoadPortableBinary(
                *packed_conjugate_cache,
                ckks_detail::CKKSDenseBootstrapNamedPath(
                    root, "packed_conjugate_galois"));
        }
        return *packed_conjugate_cache;
    }

    void release_packed_conjugate_galois() const
    {
        packed_conjugate_cache.reset();
    }

    template <std::size_t I>
    const auto &polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::relin_depth);
        auto &entry = std::get<I>(polynomial_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseEvalModPolynomialRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "polynomial_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::relin_depth);
        std::get<I>(polynomial_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        auto &entry = std::get<I>(double_angle_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseEvalModDoubleAngleRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "double_angle_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        std::get<I>(double_angle_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        auto &entry = std::get<I>(inverse_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseEvalModInverseRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(root,
                                                           "inverse_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        std::get<I>(inverse_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        auto &entry = std::get<I>(slot_to_coeff_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseBootstrapSlotToCoeffGaloisKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "slot_to_coeff_galois", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        std::get<I>(slot_to_coeff_cache).reset();
    }
};

template <class Schedule>
struct CKKSDenseBootstrapHybridGiantFilesystemKeyProvider {
    using EvalModTraits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    using CoeffToSlotCache =
        typename ckks_detail::
            CKKSDenseBootstrapHybridGiantCoeffToSlotCacheTuple<
                Schedule,
                std::make_index_sequence<
                    Schedule::coeff_to_slot_level_count + 1>>::type;
    using SlotToCoeffCache =
        typename ckks_detail::
            CKKSDenseBootstrapHybridGiantSlotToCoeffCacheTuple<
                Schedule,
                std::make_index_sequence<
                    Schedule::slot_to_coeff_level_count + 1>>::type;
    using PolynomialRelinCache =
        typename ckks_detail::CKKSDenseBootstrapPolynomialRelinCacheTuple<
            Schedule,
            std::make_index_sequence<EvalModTraits::PolynomialTraits::
                                         relin_depth>>::type;
    using DoubleAngleRelinCache =
        typename ckks_detail::CKKSDenseBootstrapDoubleAngleRelinCacheTuple<
            Schedule,
            std::make_index_sequence<Schedule::evalmod_double_angle>>::type;
    using InverseRelinCache =
        typename ckks_detail::CKKSDenseBootstrapInverseRelinCacheTuple<
            Schedule,
            std::make_index_sequence<EvalModTraits::InverseTraits::
                                         power_depth>>::type;

    std::filesystem::path root;
    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan_cache{};
    CKKSBoundedCosEvalModPolynomial evalmod_polynomial_cache{};
    mutable CoeffToSlotCache coeff_to_slot_cache{};
    mutable std::unique_ptr<CKKSDenseBootstrapPackedConjugateGaloisKey<Schedule>>
        packed_conjugate_cache{};
    mutable PolynomialRelinCache polynomial_relin_cache{};
    mutable DoubleAngleRelinCache double_angle_relin_cache{};
    mutable InverseRelinCache inverse_relin_cache{};
    mutable SlotToCoeffCache slot_to_coeff_cache{};

    explicit CKKSDenseBootstrapHybridGiantFilesystemKeyProvider(
        std::filesystem::path root_)
        : root(std::move(root_))
    {
        try {
            if (!CKKSDenseBootstrapHybridGiantKeyDirectoryManifestMatches<
                    Schedule>(root)) {
                throw std::runtime_error(
                    "CKKS hybrid bootstrap key directory manifest does not "
                    "match schedule");
            }
        }
        catch (const std::exception &e) {
            throw std::runtime_error(
                std::string(
                    "invalid CKKS hybrid bootstrap key directory manifest: ") +
                e.what());
        }
        CKKSDenseBootstrapHybridGiantRequireKeyDirectoryComplete<Schedule>(
            root);
        CKKSLoadPortableBinary(
            linear_plan_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root, "linear_plan"));
        CKKSLoadPortableBinary(
            evalmod_polynomial_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                     "evalmod_polynomial"));
    }

    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan() const
    {
        return linear_plan_cache;
    }

    const CKKSBoundedCosEvalModPolynomial &evalmod_polynomial() const
    {
        return evalmod_polynomial_cache;
    }

    template <std::size_t I>
    const auto &coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        auto &entry = std::get<I>(coeff_to_slot_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseBootstrapHybridGiantCoeffToSlotGaloisKey<Schedule,
                                                                  I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "coeff_to_slot_galois", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        std::get<I>(coeff_to_slot_cache).reset();
    }

    const auto &packed_conjugate_galois() const
    {
        if (!packed_conjugate_cache) {
            packed_conjugate_cache =
                std::make_unique<CKKSDenseBootstrapPackedConjugateGaloisKey<
                    Schedule>>();
            CKKSLoadPortableBinary(
                *packed_conjugate_cache,
                ckks_detail::CKKSDenseBootstrapNamedPath(
                    root, "packed_conjugate_galois"));
        }
        return *packed_conjugate_cache;
    }

    void release_packed_conjugate_galois() const
    {
        packed_conjugate_cache.reset();
    }

    template <std::size_t I>
    const auto &polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::relin_depth);
        auto &entry = std::get<I>(polynomial_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseEvalModPolynomialRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "polynomial_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::relin_depth);
        std::get<I>(polynomial_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        auto &entry = std::get<I>(double_angle_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseEvalModDoubleAngleRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "double_angle_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        std::get<I>(double_angle_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        auto &entry = std::get<I>(inverse_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseEvalModInverseRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(root,
                                                           "inverse_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        std::get<I>(inverse_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        auto &entry = std::get<I>(slot_to_coeff_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseBootstrapHybridGiantSlotToCoeffGaloisKey<Schedule,
                                                                  I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "slot_to_coeff_galois", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        std::get<I>(slot_to_coeff_cache).reset();
    }
};

template <class Schedule>
struct CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider {
    using EvalModTraits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    using CoeffToSlotCache =
        typename ckks_detail::
            CKKSDenseBootstrapSeededHybridGiantCoeffToSlotCacheTuple<
                Schedule,
                std::make_index_sequence<
                    Schedule::coeff_to_slot_level_count + 1>>::type;
    using SlotToCoeffCache =
        typename ckks_detail::
            CKKSDenseBootstrapSeededHybridGiantSlotToCoeffCacheTuple<
                Schedule,
                std::make_index_sequence<
                    Schedule::slot_to_coeff_level_count + 1>>::type;
    using PolynomialRelinCache =
        typename ckks_detail::
            CKKSDenseBootstrapSeededPolynomialRelinCacheTuple<
                Schedule,
                std::make_index_sequence<EvalModTraits::PolynomialTraits::
                                             relin_depth>>::type;
    using DoubleAngleRelinCache =
        typename ckks_detail::
            CKKSDenseBootstrapSeededDoubleAngleRelinCacheTuple<
                Schedule,
                std::make_index_sequence<Schedule::evalmod_double_angle>>::type;
    using InverseRelinCache =
        typename ckks_detail::CKKSDenseBootstrapSeededInverseRelinCacheTuple<
            Schedule,
            std::make_index_sequence<EvalModTraits::InverseTraits::
                                         power_depth>>::type;

    std::filesystem::path root;
    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan_cache{};
    CKKSBoundedCosEvalModPolynomial evalmod_polynomial_cache{};
    mutable CoeffToSlotCache coeff_to_slot_cache{};
    mutable std::unique_ptr<
        CKKSDenseBootstrapSeededPackedConjugateGaloisKey<Schedule>>
        packed_conjugate_cache{};
    mutable PolynomialRelinCache polynomial_relin_cache{};
    mutable DoubleAngleRelinCache double_angle_relin_cache{};
    mutable InverseRelinCache inverse_relin_cache{};
    mutable SlotToCoeffCache slot_to_coeff_cache{};

    explicit CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider(
        std::filesystem::path root_)
        : root(std::move(root_))
    {
        try {
            if (!CKKSDenseBootstrapSeededHybridGiantKeyDirectoryManifestMatches<
                    Schedule>(root)) {
                throw std::runtime_error(
                    "CKKS seeded hybrid bootstrap key directory manifest does "
                    "not match schedule");
            }
        }
        catch (const std::exception &e) {
            throw std::runtime_error(
                std::string("invalid CKKS seeded hybrid bootstrap key "
                            "directory manifest: ") +
                e.what());
        }
        CKKSDenseBootstrapSeededHybridGiantRequireKeyDirectoryComplete<
            Schedule>(root);
        CKKSLoadPortableBinary(
            linear_plan_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root, "linear_plan"));
        CKKSLoadPortableBinary(
            evalmod_polynomial_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                     "evalmod_polynomial"));
    }

    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan() const
    {
        return linear_plan_cache;
    }

    const CKKSBoundedCosEvalModPolynomial &evalmod_polynomial() const
    {
        return evalmod_polynomial_cache;
    }

    template <std::size_t I>
    const auto &coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        auto &entry = std::get<I>(coeff_to_slot_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseBootstrapSeededHybridGiantCoeffToSlotGaloisKey<
                    Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "coeff_to_slot_galois", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        std::get<I>(coeff_to_slot_cache).reset();
    }

    const auto &packed_conjugate_galois() const
    {
        if (!packed_conjugate_cache) {
            packed_conjugate_cache = std::make_unique<
                CKKSDenseBootstrapSeededPackedConjugateGaloisKey<Schedule>>();
            CKKSLoadPortableBinary(
                *packed_conjugate_cache,
                ckks_detail::CKKSDenseBootstrapNamedPath(
                    root, "packed_conjugate_galois"));
        }
        return *packed_conjugate_cache;
    }

    void release_packed_conjugate_galois() const
    {
        packed_conjugate_cache.reset();
    }

    template <std::size_t I>
    const auto &polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::relin_depth);
        auto &entry = std::get<I>(polynomial_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseSeededEvalModPolynomialRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "polynomial_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::relin_depth);
        std::get<I>(polynomial_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        auto &entry = std::get<I>(double_angle_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseSeededEvalModDoubleAngleRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "double_angle_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        std::get<I>(double_angle_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        auto &entry = std::get<I>(inverse_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseSeededEvalModInverseRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(root,
                                                           "inverse_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        std::get<I>(inverse_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        auto &entry = std::get<I>(slot_to_coeff_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseBootstrapSeededHybridGiantSlotToCoeffGaloisKey<
                    Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "slot_to_coeff_galois", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        std::get<I>(slot_to_coeff_cache).reset();
    }
};

template <class Schedule>
struct CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider {
    using P = typename Schedule::Param;
    using EvalModTraits = CKKSDenseEvalModBoundedCosTraits<Schedule>;
    using CoeffToSlotCache =
        typename ckks_detail::
            CKKSDenseBootstrapSeededHybridGiantStreamedCoeffToSlotCacheTuple<
                Schedule,
                std::make_index_sequence<
                    Schedule::coeff_to_slot_level_count + 1>>::type;
    using SlotToCoeffCache =
        typename ckks_detail::
            CKKSDenseBootstrapSeededHybridGiantStreamedSlotToCoeffCacheTuple<
                Schedule,
                std::make_index_sequence<
                    Schedule::slot_to_coeff_level_count + 1>>::type;
    using PolynomialRelinCache =
        typename ckks_detail::
            CKKSDenseBootstrapSeededPolynomialRelinCacheTuple<
                Schedule,
                std::make_index_sequence<EvalModTraits::PolynomialTraits::
                                             relin_depth>>::type;
    using DoubleAngleRelinCache =
        typename ckks_detail::
            CKKSDenseBootstrapSeededDoubleAngleRelinCacheTuple<
                Schedule,
                std::make_index_sequence<Schedule::evalmod_double_angle>>::type;
    using InverseRelinCache =
        typename ckks_detail::CKKSDenseBootstrapSeededInverseRelinCacheTuple<
            Schedule,
            std::make_index_sequence<EvalModTraits::InverseTraits::
                                         power_depth>>::type;

    std::filesystem::path root;
    CKKSDenseBootstrapLinearPlan<Schedule> linear_plan_cache{};
    CKKSDenseBootstrapHybridGiantRotationKeyUsage<Schedule>
        rotation_usage_cache{};
    CKKSBoundedCosEvalModPolynomial evalmod_polynomial_cache{};
    mutable CoeffToSlotCache coeff_to_slot_cache{};
    mutable std::unique_ptr<
        CKKSDenseBootstrapSeededPackedConjugateGaloisKey<Schedule>>
        packed_conjugate_cache{};
    mutable PolynomialRelinCache polynomial_relin_cache{};
    mutable DoubleAngleRelinCache double_angle_relin_cache{};
    mutable InverseRelinCache inverse_relin_cache{};
    mutable SlotToCoeffCache slot_to_coeff_cache{};

    explicit CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider(
        std::filesystem::path root_)
        : root(std::move(root_))
    {
        try {
            if (!CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryManifestMatches<
                    Schedule>(root)) {
                throw std::runtime_error(
                    "CKKS streamed seeded hybrid bootstrap key directory "
                    "manifest does not match schedule");
            }
        }
        catch (const std::exception &e) {
            throw std::runtime_error(
                std::string("invalid CKKS streamed seeded hybrid bootstrap key "
                            "directory manifest: ") +
                e.what());
        }
        CKKSDenseBootstrapSeededHybridGiantStreamedRequireKeyDirectoryComplete<
            Schedule>(root);
        CKKSLoadPortableBinary(
            linear_plan_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root, "linear_plan"));
        CKKSLoadPortableBinary(
            rotation_usage_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root, "rotation_usage"));
        CKKSLoadPortableBinary(
            evalmod_polynomial_cache,
            ckks_detail::CKKSDenseBootstrapNamedPath(root,
                                                     "evalmod_polynomial"));
    }

    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan() const
    {
        return linear_plan_cache;
    }

    const CKKSBoundedCosEvalModPolynomial &evalmod_polynomial() const
    {
        return evalmod_polynomial_cache;
    }

    template <std::size_t I>
    const auto &coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        auto &entry = std::get<I>(coeff_to_slot_cache);
        if (!entry) {
            constexpr std::uint32_t log_q =
                Schedule::boot_log_q -
                I * Schedule::coeff_to_slot_plain_log_delta;
            entry = std::make_unique<
                ckks_detail::CKKSSeededHybridSparseGaloisKeyFileView<P,
                                                                     log_q>>();
            entry->binary.root = root;
            entry->binary.prefix =
                ckks_detail::CKKSDenseBootstrapLeveledKeyPrefix(
                    "coeff_to_slot_galois", I, "binary");
            entry->binary.available =
                rotation_usage_cache.coeff_to_slot_binary[I];
            entry->direct.root = root;
            entry->direct.prefix =
                ckks_detail::CKKSDenseBootstrapLeveledKeyPrefix(
                    "coeff_to_slot_galois", I, "direct");
            entry->direct.available =
                rotation_usage_cache.coeff_to_slot_direct[I];
        }
        return *entry;
    }

    template <std::size_t I>
    void release_coeff_to_slot_galois() const
    {
        static_assert(I <= Schedule::coeff_to_slot_level_count);
        std::get<I>(coeff_to_slot_cache).reset();
    }

    const auto &packed_conjugate_galois() const
    {
        if (!packed_conjugate_cache) {
            packed_conjugate_cache = std::make_unique<
                CKKSDenseBootstrapSeededPackedConjugateGaloisKey<Schedule>>();
            CKKSLoadPortableBinary(
                *packed_conjugate_cache,
                ckks_detail::CKKSDenseBootstrapNamedPath(
                    root, "packed_conjugate_galois"));
        }
        return *packed_conjugate_cache;
    }

    void release_packed_conjugate_galois() const
    {
        packed_conjugate_cache.reset();
    }

    template <std::size_t I>
    const auto &polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::relin_depth);
        auto &entry = std::get<I>(polynomial_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseSeededEvalModPolynomialRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "polynomial_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_polynomial_relin() const
    {
        static_assert(I < EvalModTraits::PolynomialTraits::relin_depth);
        std::get<I>(polynomial_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        auto &entry = std::get<I>(double_angle_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseSeededEvalModDoubleAngleRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(
                    root, "double_angle_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_double_angle_relin() const
    {
        static_assert(I < Schedule::evalmod_double_angle);
        std::get<I>(double_angle_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        auto &entry = std::get<I>(inverse_relin_cache);
        if (!entry) {
            entry = std::make_unique<
                CKKSDenseSeededEvalModInverseRelinKey<Schedule, I>>();
            CKKSLoadPortableBinary(
                *entry,
                ckks_detail::CKKSDenseBootstrapIndexedPath(root,
                                                           "inverse_relin", I));
        }
        return *entry;
    }

    template <std::size_t I>
    void release_inverse_relin() const
    {
        static_assert(I < EvalModTraits::InverseTraits::power_depth);
        std::get<I>(inverse_relin_cache).reset();
    }

    template <std::size_t I>
    const auto &slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        auto &entry = std::get<I>(slot_to_coeff_cache);
        if (!entry) {
            constexpr std::uint32_t log_q =
                Schedule::after_evalmod_log_q -
                I * Schedule::slot_to_coeff_plain_log_delta;
            entry = std::make_unique<
                ckks_detail::CKKSSeededHybridSparseGaloisKeyFileView<P,
                                                                     log_q>>();
            entry->binary.root = root;
            entry->binary.prefix =
                ckks_detail::CKKSDenseBootstrapLeveledKeyPrefix(
                    "slot_to_coeff_galois", I, "binary");
            entry->binary.available =
                rotation_usage_cache.slot_to_coeff_binary[I];
            entry->direct.root = root;
            entry->direct.prefix =
                ckks_detail::CKKSDenseBootstrapLeveledKeyPrefix(
                    "slot_to_coeff_galois", I, "direct");
            entry->direct.available =
                rotation_usage_cache.slot_to_coeff_direct[I];
        }
        return *entry;
    }

    template <std::size_t I>
    void release_slot_to_coeff_galois() const
    {
        static_assert(I <= Schedule::slot_to_coeff_level_count);
        std::get<I>(slot_to_coeff_cache).reset();
    }
};

namespace ckks_detail {

template <class KeyProvider, bool CoeffToSlot>
struct CKKSDenseBootstrapLinearKeyProviderChain {
    const KeyProvider &provider;

    template <std::size_t I>
    decltype(auto) get() const
    {
        if constexpr (CoeffToSlot)
            return provider.template coeff_to_slot_galois<I>();
        else
            return provider.template slot_to_coeff_galois<I>();
    }

    template <std::size_t I>
    void release() const
    {
        if constexpr (CoeffToSlot) {
            if constexpr (requires {
                              provider.template release_coeff_to_slot_galois<I>();
                          }) {
                provider.template release_coeff_to_slot_galois<I>();
            }
        }
        else {
            if constexpr (requires {
                              provider.template release_slot_to_coeff_galois<I>();
                          }) {
                provider.template release_slot_to_coeff_galois<I>();
            }
        }
    }
};

template <class KeyProvider>
struct CKKSDenseBootstrapRetainEvalModRelinKeyProvider {
    const KeyProvider &provider;

    template <std::size_t I>
    decltype(auto) polynomial_relin() const
    {
        return provider.template polynomial_relin<I>();
    }

    template <std::size_t I>
    void release_polynomial_relin() const
    {}

    template <std::size_t I>
    decltype(auto) double_angle_relin() const
    {
        return provider.template double_angle_relin<I>();
    }

    template <std::size_t I>
    void release_double_angle_relin() const
    {}

    template <std::size_t I>
    decltype(auto) inverse_relin() const
    {
        return provider.template inverse_relin<I>();
    }

    template <std::size_t I>
    void release_inverse_relin() const
    {}
};

inline void CKKSDenseBootstrapProgressBegin(
    const CKKSDenseBootstrapProgress *progress, const char *stage)
{
    if (progress != nullptr && progress->stage_begin != nullptr)
        progress->stage_begin(stage, progress->context);
}

inline void CKKSDenseBootstrapProgressEnd(
    const CKKSDenseBootstrapProgress *progress, const char *stage,
    double elapsed_ms)
{
    if (progress != nullptr && progress->stage_end != nullptr)
        progress->stage_end(stage, elapsed_ms, progress->context);
}

template <class F>
inline void CKKSTimeBootstrapStage(
    double *duration_ms, F &&stage,
    const CKKSDenseBootstrapProgress *progress, const char *stage_name)
{
    CKKSDenseBootstrapProgressBegin(progress, stage_name);
    const auto start = std::chrono::steady_clock::now();
    std::forward<F>(stage)();
    const auto stop = std::chrono::steady_clock::now();
    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(stop - start).count();
    if (duration_ms != nullptr) *duration_ms = elapsed_ms;
    CKKSDenseBootstrapProgressEnd(progress, stage_name, elapsed_ms);
}

template <class F>
inline void CKKSTimeBootstrapStage(double *duration_ms, F &&stage)
{
    CKKSTimeBootstrapStage(duration_ms, std::forward<F>(stage), nullptr,
                           nullptr);
}

template <std::size_t I, class Schedule, class GaloisKeyChain>
inline void CKKSDenseBootstrapCoeffToSlotStagesBSGSImpl(
    typename Schedule::CoeffToSlotCiphertext &res,
    const CKKSCiphertext<
        typename Schedule::Param,
        Schedule::boot_log_q -
            I * Schedule::coeff_to_slot_plain_log_delta,
        Schedule::log_delta> &ct,
    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan,
    const GaloisKeyChain &gk_chain)
{
    using P = typename Schedule::Param;
    if constexpr (I == Schedule::coeff_to_slot_level_count) {
        res = ct;
    }
    else {
        constexpr std::uint32_t log_q =
            Schedule::boot_log_q -
            I * Schedule::coeff_to_slot_plain_log_delta;
        constexpr int k_step =
            CKKSDenseBootstrapCoeffToSlotBSGSStep<I, Schedule>();
        CKKSLinearTransformPlan<P, log_q, Schedule::log_delta,
                                Schedule::coeff_to_slot_plain_log_delta>
            plan;
        CKKSBuildLinearTransformBSGSPlan<
            P, log_q, Schedule::log_delta,
            Schedule::coeff_to_slot_plain_log_delta>(
            plan, linear_plan.coeff_to_slot_stages[I], k_step);

        auto next = std::make_unique<CKKSPlainMulResult<
            P, log_q, Schedule::log_delta,
            Schedule::coeff_to_slot_plain_log_delta>>();
        CKKSLinearTransformBSGS<
            P, log_q, Schedule::log_delta,
            Schedule::coeff_to_slot_plain_log_delta>(
            *next, ct, plan, gk_chain.template get<I>(),
            gk_chain.template get<I + 1>());
        maybe_release_key<I>(gk_chain);
        if constexpr (I + 1 == Schedule::coeff_to_slot_level_count) {
            maybe_release_key<I + 1>(gk_chain);
        }
        CKKSDenseBootstrapCoeffToSlotStagesBSGSImpl<I + 1, Schedule>(
            res, *next, linear_plan, gk_chain);
    }
}

template <class Schedule, class GaloisKeyChain>
inline void CKKSDenseBootstrapCoeffToSlotStagesBSGS(
    typename Schedule::CoeffToSlotCiphertext &res,
    const typename Schedule::BootstrapCiphertext &ct,
    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan,
    const GaloisKeyChain &gk_chain)
{
    CKKSDenseBootstrapCoeffToSlotStagesBSGSImpl<0, Schedule>(
        res, ct, linear_plan, gk_chain);
}

template <std::size_t I, class Schedule, class GaloisKeyChain>
inline void CKKSDenseBootstrapSlotToCoeffTailStagesBSGSImpl(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<
        typename Schedule::Param,
        Schedule::after_evalmod_log_q -
            I * Schedule::slot_to_coeff_plain_log_delta,
        Schedule::log_delta> &ct,
    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan,
    const GaloisKeyChain &gk_chain)
{
    using P = typename Schedule::Param;
    if constexpr (I == Schedule::slot_to_coeff_level_count) {
        res = ct;
    }
    else {
        constexpr std::uint32_t log_q =
            Schedule::after_evalmod_log_q -
            I * Schedule::slot_to_coeff_plain_log_delta;
        constexpr int k_step =
            CKKSDenseBootstrapSlotToCoeffBSGSStep<I, Schedule>();
        CKKSLinearTransformPlan<P, log_q, Schedule::log_delta,
                                Schedule::slot_to_coeff_plain_log_delta>
            plan;
        CKKSBuildLinearTransformBSGSPlan<
            P, log_q, Schedule::log_delta,
            Schedule::slot_to_coeff_plain_log_delta>(
            plan, linear_plan.slot_to_coeff_stages[I], k_step);

        auto next = std::make_unique<CKKSPlainMulResult<
            P, log_q, Schedule::log_delta,
            Schedule::slot_to_coeff_plain_log_delta>>();
        CKKSLinearTransformBSGS<
            P, log_q, Schedule::log_delta,
            Schedule::slot_to_coeff_plain_log_delta>(
            *next, ct, plan, gk_chain.template get<I>(),
            gk_chain.template get<I + 1>());
        maybe_release_key<I>(gk_chain);
        if constexpr (I + 1 == Schedule::slot_to_coeff_level_count) {
            maybe_release_key<I + 1>(gk_chain);
        }
        CKKSDenseBootstrapSlotToCoeffTailStagesBSGSImpl<I + 1, Schedule>(
            res, *next, linear_plan, gk_chain);
    }
}

template <class Schedule, class GaloisKeyChain>
inline void CKKSDenseBootstrapSlotToCoeffStagesBSGSDualInputSharedTail(
    typename Schedule::OutputCiphertext &res,
    const CKKSDenseEvalModBoundedCosResult<Schedule> &real_evalmod,
    const CKKSDenseEvalModBoundedCosResult<Schedule> &imag_evalmod,
    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan,
    const GaloisKeyChain &gk_chain)
{
    using P = typename Schedule::Param;
    static_assert(Schedule::slot_to_coeff_level_count > 0);
    constexpr std::uint32_t log_q = Schedule::after_evalmod_log_q;
    constexpr int k_step =
        CKKSDenseBootstrapSlotToCoeffBSGSStep<0, Schedule>();

    CKKSLinearTransformPlan<P, log_q, Schedule::log_delta,
                            Schedule::slot_to_coeff_plain_log_delta>
        real_plan;
    CKKSLinearTransformPlan<P, log_q, Schedule::log_delta,
                            Schedule::slot_to_coeff_plain_log_delta>
        imag_plan;
    CKKSBuildLinearTransformBSGSPlan<
        P, log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta>(
        real_plan, linear_plan.slot_to_coeff_stages[0], k_step);
    CKKSBuildLinearTransformBSGSPlan<
        P, log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta>(
        imag_plan, linear_plan.slot_to_coeff_imag_stages[0], k_step);

    auto next = std::make_unique<CKKSPlainMulResult<
        P, log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta>>();
    CKKSLinearTransformBSGSDualInput<
        P, log_q, Schedule::log_delta,
        Schedule::slot_to_coeff_plain_log_delta>(
        *next, real_evalmod, imag_evalmod, real_plan, imag_plan,
        gk_chain.template get<0>(), gk_chain.template get<1>());
    maybe_release_key<0>(gk_chain);

    if constexpr (Schedule::slot_to_coeff_level_count == 1) {
        res = *next;
        maybe_release_key<1>(gk_chain);
    }
    else {
        CKKSDenseBootstrapSlotToCoeffTailStagesBSGSImpl<1, Schedule>(
            res, *next, linear_plan, gk_chain);
    }
}

template <class Schedule, class KeyProvider>
inline void CKKSDenseBootstrapWithKeyProviderImpl(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const KeyProvider &key_provider, CKKSDenseBootstrapTimings *timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    using P = typename Schedule::Param;
    const CKKSDenseBootstrapLinearPlan<Schedule> &linear_plan =
        key_provider.linear_plan();

    auto raised = std::make_unique<typename Schedule::BootstrapCiphertext>();
    CKKSTimeBootstrapStage(timings == nullptr ? nullptr : &timings->modraise_ms,
                           [&] {
                               CKKSModRaiseBoundedPhaseRandomized<
                                   P, Schedule::input_log_q,
                                   Schedule::boot_log_q, Schedule::log_delta,
                                   Schedule::modraise_mask_bound>(*raised, ct);
                           },
                           progress, "modraise");

    auto coeff_to_slot =
        std::make_unique<typename Schedule::CoeffToSlotCiphertext>();
    const CKKSDenseBootstrapLinearKeyProviderChain<KeyProvider, true>
        coeff_to_slot_galois{key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->coeff_to_slot_ms, [&] {
            CKKSDenseBootstrapCoeffToSlotStagesBSGS<Schedule>(
                *coeff_to_slot, *raised, linear_plan, coeff_to_slot_galois);
        },
        progress, "coeff_to_slot");
    raised.reset();

    auto real_component =
        std::make_unique<typename Schedule::ComponentCiphertext>();
    auto imag_component =
        std::make_unique<typename Schedule::ComponentCiphertext>();
    CKKSTimeBootstrapStage(timings == nullptr ? nullptr : &timings->split_ms,
                           [&] {
                               CKKSExtractRealImagSlots<
                                   P, Schedule::after_coeff_to_slot_log_q,
                                   Schedule::log_delta,
                                   Schedule::component_split_plain_log_delta>(
                                   *real_component, *imag_component,
                                   *coeff_to_slot,
                                   key_provider.packed_conjugate_galois());
                           },
                           progress, "split");
    coeff_to_slot.reset();
    if constexpr (requires { key_provider.release_packed_conjugate_galois(); }) {
        key_provider.release_packed_conjugate_galois();
    }

    auto real_evalmod =
        std::make_unique<CKKSDenseEvalModBoundedCosResult<Schedule>>();
    const CKKSDenseBootstrapRetainEvalModRelinKeyProvider<KeyProvider>
        retained_evalmod_relin{key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->real_evalmod_ms, [&] {
            CKKSDenseEvalModBoundedCosNormalizedWithKeyProvider<Schedule>(
                *real_evalmod, *real_component, key_provider.evalmod_polynomial(),
                retained_evalmod_relin);
        },
        progress, "real_evalmod");
    real_component.reset();
    auto imag_evalmod =
        std::make_unique<CKKSDenseEvalModBoundedCosResult<Schedule>>();
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->imag_evalmod_ms, [&] {
            CKKSDenseEvalModBoundedCosNormalizedWithKeyProvider<Schedule>(
                *imag_evalmod, *imag_component, key_provider.evalmod_polynomial(),
                key_provider);
        },
        progress, "imag_evalmod");
    imag_component.reset();

    const CKKSDenseBootstrapLinearKeyProviderChain<KeyProvider, false>
        slot_to_coeff_galois{key_provider};
    CKKSTimeBootstrapStage(
        timings == nullptr ? nullptr : &timings->slot_to_coeff_ms, [&] {
            CKKSDenseBootstrapSlotToCoeffStagesBSGSDualInputSharedTail<
                Schedule>(res, *real_evalmod, *imag_evalmod, linear_plan,
                          slot_to_coeff_galois);
        },
        progress, "slot_to_coeff");
    real_evalmod.reset();
    imag_evalmod.reset();
}

}  // namespace ckks_detail

template <class Schedule, class KeyProvider>
inline void CKKSDenseBootstrapWithKeyProvider(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const KeyProvider &key_provider)
{
    ckks_detail::CKKSDenseBootstrapWithKeyProviderImpl<Schedule>(
        res, ct, key_provider, nullptr);
}

template <class Schedule, class KeyProvider>
inline void CKKSDenseBootstrapWithKeyProviderTimed(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const KeyProvider &key_provider, CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    timings = {};
    ckks_detail::CKKSDenseBootstrapWithKeyProviderImpl<Schedule>(
        res, ct, key_provider, &timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapNormalizeInput(
    typename Schedule::InputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct)
{
    using P = typename Schedule::Param;
    static_assert(InLogDelta >= Schedule::log_delta,
                  "CKKS bootstrap input normalization cannot raise scale");
    constexpr std::uint32_t scale_drop = InLogDelta - Schedule::log_delta;
    static_assert(InLogQ >= Schedule::input_log_q + scale_drop,
                  "CKKS bootstrap input normalization cannot raise a residual "
                  "ciphertext level after scale-down");
    if constexpr (scale_drop == 0) {
        CKKSLevelReduce<P, InLogQ, Schedule::input_log_q,
                        Schedule::log_delta>(res, ct);
    }
    else {
        using ScaledCt =
            CKKSScaleDownResult<P, InLogQ, InLogDelta, scale_drop>;
        auto scaled = std::make_unique<ScaledCt>();
        CKKSScaleDown<P, InLogQ, InLogDelta, scale_drop>(*scaled, ct);
        CKKSLevelReduce<P, ScaledCt::log_q, Schedule::input_log_q,
                        Schedule::log_delta>(res, *scaled);
    }
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta,
          class KeyProvider>
inline void CKKSDenseBootstrapFromLevelWithKeyProvider(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const KeyProvider &key_provider)
{
    auto normalized = std::make_unique<typename Schedule::InputCiphertext>();
    CKKSDenseBootstrapNormalizeInput<Schedule>(*normalized, ct);
    CKKSDenseBootstrapWithKeyProvider<Schedule>(res, *normalized, key_provider);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta,
          class KeyProvider>
inline void CKKSDenseBootstrapFromLevelWithKeyProviderTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const KeyProvider &key_provider, CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    timings = {};
    auto normalized = std::make_unique<typename Schedule::InputCiphertext>();
    ckks_detail::CKKSTimeBootstrapStage(
        &timings.normalize_ms,
        [&] { CKKSDenseBootstrapNormalizeInput<Schedule>(*normalized, ct); },
        progress, "normalize");
    ckks_detail::CKKSDenseBootstrapWithKeyProviderImpl<Schedule>(
        res, *normalized, key_provider, &timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta,
          class KeyProvider, class EncapsulationKey>
inline void CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const KeyProvider &bootstrap_key_provider,
    const EncapsulationKey &encapsulation_key)
{
    auto normalized = std::make_unique<typename Schedule::InputCiphertext>();
    CKKSDenseBootstrapNormalizeInput<Schedule>(*normalized, ct);

    auto bootstrap_input = std::make_unique<typename Schedule::InputCiphertext>();
    CKKSSecretKeySwitch<typename Schedule::Param, Schedule::input_log_q,
                        Schedule::log_delta>(
        *bootstrap_input, *normalized, encapsulation_key.input_to_bootstrap);

    auto bootstrap_output =
        std::make_unique<typename Schedule::OutputCiphertext>();
    CKKSDenseBootstrapWithKeyProvider<Schedule>(
        *bootstrap_output, *bootstrap_input, bootstrap_key_provider);

    CKKSSecretKeySwitch<typename Schedule::Param, Schedule::output_log_q,
                        Schedule::log_delta>(
        res, *bootstrap_output, encapsulation_key.bootstrap_to_output);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta,
          class KeyProvider, class EncapsulationKey>
inline void CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const KeyProvider &bootstrap_key_provider,
    const EncapsulationKey &encapsulation_key,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    timings = {};
    auto normalized = std::make_unique<typename Schedule::InputCiphertext>();
    ckks_detail::CKKSTimeBootstrapStage(
        &timings.normalize_ms,
        [&] { CKKSDenseBootstrapNormalizeInput<Schedule>(*normalized, ct); },
        progress, "normalize");

    auto bootstrap_input = std::make_unique<typename Schedule::InputCiphertext>();
    ckks_detail::CKKSTimeBootstrapStage(&timings.input_secret_switch_ms, [&] {
        CKKSSecretKeySwitch<typename Schedule::Param, Schedule::input_log_q,
                            Schedule::log_delta>(
            *bootstrap_input, *normalized, encapsulation_key.input_to_bootstrap);
    },
                                        progress, "input_secret_switch");

    auto bootstrap_output =
        std::make_unique<typename Schedule::OutputCiphertext>();
    ckks_detail::CKKSDenseBootstrapWithKeyProviderImpl<Schedule>(
        *bootstrap_output, *bootstrap_input, bootstrap_key_provider, &timings,
        progress);

    ckks_detail::CKKSTimeBootstrapStage(&timings.output_secret_switch_ms, [&] {
        CKKSSecretKeySwitch<typename Schedule::Param, Schedule::output_log_q,
                            Schedule::log_delta>(
            res, *bootstrap_output, encapsulation_key.bootstrap_to_output);
    },
                                        progress, "output_secret_switch");
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapEncapsulatedFromLevelWithFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file)
{
    CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider<Schedule>(
        res, ct, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapEncapsulatedFromLevelWithFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, provider, *encapsulation_key, timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapEncapsulatedFromLevelWithHybridGiantFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file)
{
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider<Schedule>(
        res, ct, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithHybridGiantFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, provider, *encapsulation_key, timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapFromLevelWithSeededHybridGiantFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    CKKSDenseBootstrapFromLevelWithKeyProvider<Schedule>(res, ct, provider);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapFromLevelWithSeededHybridGiantFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir, CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    CKKSDenseBootstrapFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, provider, timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithSeededHybridGiantFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider<Schedule>(
        res, ct, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithSeededHybridGiantFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, provider, *encapsulation_key, timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithSeededHybridGiantFilesystemSeededKeys(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider<Schedule>(
        res, ct, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithSeededHybridGiantFilesystemSeededKeysTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, provider, *encapsulation_key, timings, progress);
}

template <class Schedule>
inline void CKKSDenseBootstrapWithSeededHybridGiantStreamedFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const std::filesystem::path &key_dir)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    CKKSDenseBootstrapWithKeyProvider<Schedule>(res, ct, provider);
}

template <class Schedule>
inline void CKKSDenseBootstrapWithSeededHybridGiantStreamedFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const std::filesystem::path &key_dir, CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    CKKSDenseBootstrapWithKeyProviderTimed<Schedule>(
        res, ct, provider, timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapFromLevelWithSeededHybridGiantStreamedFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    CKKSDenseBootstrapFromLevelWithKeyProvider<Schedule>(res, ct, provider);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapFromLevelWithSeededHybridGiantStreamedFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir, CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    CKKSDenseBootstrapFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, provider, timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithSeededHybridGiantStreamedFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider<Schedule>(
        res, ct, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithSeededHybridGiantStreamedFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, provider, *encapsulation_key, timings, progress);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithSeededHybridGiantStreamedFilesystemSeededKeys(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider<Schedule>(
        res, ct, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void
CKKSDenseBootstrapEncapsulatedFromLevelWithSeededHybridGiantStreamedFilesystemSeededKeysTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    CKKSDenseBootstrapTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, provider, *encapsulation_key, timings, progress);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
struct CKKSDenseBootstrapProductTraits {
    using P = typename Schedule::Param;
    using ProductCiphertext =
        CKKSMultResult<P, LhsLogQ, LhsLogDelta, RhsLogQ, RhsLogDelta>;

    static constexpr std::uint32_t scale_drop =
        ProductCiphertext::log_delta >= Schedule::log_delta
            ? ProductCiphertext::log_delta - Schedule::log_delta
            : 0;
    static constexpr bool product_scale_compatible =
        ProductCiphertext::log_delta >= Schedule::log_delta;
    static constexpr bool product_level_compatible =
        ProductCiphertext::log_q >= Schedule::input_log_q + scale_drop;

    static_assert(product_scale_compatible,
                  "CKKS product bootstrap cannot raise product scale");
    static_assert(product_level_compatible,
                  "CKKS product bootstrap needs enough product level for "
                  "bootstrap input normalization");
};

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta, class KeyProvider>
inline void CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProvider(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &relin_key_file,
    const KeyProvider &bootstrap_key_provider,
    const CKKSDenseBootstrapEncapsulationKey<Schedule> &encapsulation_key)
{
    using Traits =
        CKKSDenseBootstrapProductTraits<Schedule, LhsLogQ, LhsLogDelta,
                                        RhsLogQ, RhsLogDelta>;
    auto product = std::make_unique<typename Traits::ProductCiphertext>();
    CKKSMultWithRelinKeyFile<typename Schedule::Param>(
        *product, lhs, rhs, relin_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider<Schedule>(
        res, *product, bootstrap_key_provider, encapsulation_key);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta, class KeyProvider>
inline void CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProviderTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &relin_key_file,
    const KeyProvider &bootstrap_key_provider,
    const CKKSDenseBootstrapEncapsulationKey<Schedule> &encapsulation_key,
    CKKSDenseBootstrapProductTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    using Traits =
        CKKSDenseBootstrapProductTraits<Schedule, LhsLogQ, LhsLogDelta,
                                        RhsLogQ, RhsLogDelta>;
    timings = {};
    auto product = std::make_unique<typename Traits::ProductCiphertext>();
    ckks_detail::CKKSTimeBootstrapStage(&timings.multiply_ms, [&] {
        CKKSMultWithRelinKeyFile<typename Schedule::Param>(
            *product, lhs, rhs, relin_key_file);
    },
                                        progress, "multiply");
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed<Schedule>(
        res, *product, bootstrap_key_provider, encapsulation_key,
        timings.bootstrap, progress);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta, class KeyProvider, class EncapsulationKey>
inline void CKKSDenseBootstrapProductWithSeededRelinKeyFileAndKeyProvider(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &relin_key_file,
    const KeyProvider &bootstrap_key_provider,
    const EncapsulationKey &encapsulation_key)
{
    using Traits =
        CKKSDenseBootstrapProductTraits<Schedule, LhsLogQ, LhsLogDelta,
                                        RhsLogQ, RhsLogDelta>;
    auto product = std::make_unique<typename Traits::ProductCiphertext>();
    CKKSMultWithSeededRelinKeyFile<typename Schedule::Param>(
        *product, lhs, rhs, relin_key_file);
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProvider<Schedule>(
        res, *product, bootstrap_key_provider, encapsulation_key);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta, class KeyProvider, class EncapsulationKey>
inline void CKKSDenseBootstrapProductWithSeededRelinKeyFileAndKeyProviderTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &relin_key_file,
    const KeyProvider &bootstrap_key_provider,
    const EncapsulationKey &encapsulation_key,
    CKKSDenseBootstrapProductTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    using Traits =
        CKKSDenseBootstrapProductTraits<Schedule, LhsLogQ, LhsLogDelta,
                                        RhsLogQ, RhsLogDelta>;
    timings = {};
    auto product = std::make_unique<typename Traits::ProductCiphertext>();
    ckks_detail::CKKSTimeBootstrapStage(&timings.multiply_ms, [&] {
        CKKSMultWithSeededRelinKeyFile<typename Schedule::Param>(
            *product, lhs, rhs, relin_key_file);
    },
                                        progress, "multiply");
    CKKSDenseBootstrapEncapsulatedFromLevelWithKeyProviderTimed<Schedule>(
        res, *product, bootstrap_key_provider, encapsulation_key,
        timings.bootstrap, progress);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void CKKSDenseBootstrapProductWithFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file)
{
    CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProvider<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void CKKSDenseBootstrapProductWithFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file,
    CKKSDenseBootstrapProductTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapFilesystemKeyProvider<Schedule> provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProviderTimed<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key, timings,
        progress);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void CKKSDenseBootstrapProductWithHybridGiantFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file)
{
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProvider<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void CKKSDenseBootstrapProductWithHybridGiantFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file,
    CKKSDenseBootstrapProductTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProviderTimed<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key, timings,
        progress);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void CKKSDenseBootstrapProductWithSeededHybridGiantFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProvider<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void CKKSDenseBootstrapProductWithSeededHybridGiantFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file,
    CKKSDenseBootstrapProductTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProviderTimed<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key, timings,
        progress);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void
CKKSDenseBootstrapProductWithSeededHybridGiantFilesystemSeededKeys(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithSeededRelinKeyFileAndKeyProvider<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void
CKKSDenseBootstrapProductWithSeededHybridGiantFilesystemSeededKeysTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file,
    CKKSDenseBootstrapProductTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<Schedule> provider(
        key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithSeededRelinKeyFileAndKeyProviderTimed<
        Schedule>(res, lhs, rhs, relin_key_file, provider, *encapsulation_key,
                  timings, progress);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void CKKSDenseBootstrapProductWithSeededHybridGiantStreamedFilesystemKey(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProvider<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void
CKKSDenseBootstrapProductWithSeededHybridGiantStreamedFilesystemKeyTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file,
    CKKSDenseBootstrapProductTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithRelinKeyFileAndKeyProviderTimed<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key, timings,
        progress);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void
CKKSDenseBootstrapProductWithSeededHybridGiantStreamedFilesystemSeededKeys(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithSeededRelinKeyFileAndKeyProvider<Schedule>(
        res, lhs, rhs, relin_key_file, provider, *encapsulation_key);
}

template <class Schedule, std::uint32_t LhsLogQ,
          std::uint32_t LhsLogDelta, std::uint32_t RhsLogQ,
          std::uint32_t RhsLogDelta>
inline void
CKKSDenseBootstrapProductWithSeededHybridGiantStreamedFilesystemSeededKeysTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, LhsLogQ, LhsLogDelta> &lhs,
    const CKKSCiphertext<typename Schedule::Param, RhsLogQ, RhsLogDelta> &rhs,
    const std::filesystem::path &key_dir,
    const std::filesystem::path &encapsulation_key_file,
    const std::filesystem::path &relin_key_file,
    CKKSDenseBootstrapProductTimings &timings,
    const CKKSDenseBootstrapProgress *progress = nullptr)
{
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<Schedule>
        provider(key_dir);
    auto encapsulation_key =
        std::make_unique<CKKSDenseBootstrapSeededEncapsulationKey<Schedule>>();
    CKKSDenseBootstrapLoadSeededEncapsulationKeyFromFile<Schedule>(
        *encapsulation_key, encapsulation_key_file);
    CKKSDenseBootstrapProductWithSeededRelinKeyFileAndKeyProviderTimed<
        Schedule>(res, lhs, rhs, relin_key_file, provider, *encapsulation_key,
                  timings, progress);
}

template <class Schedule>
inline void CKKSDenseBootstrap(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSDenseBootstrapKey<Schedule> &bootstrap_key)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<Schedule> key_provider(
        bootstrap_key);
    CKKSDenseBootstrapWithKeyProvider<Schedule>(res, ct, key_provider);
}

template <class Schedule>
inline void CKKSDenseBootstrapTimed(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSDenseBootstrapKey<Schedule> &bootstrap_key,
    CKKSDenseBootstrapTimings &timings)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<Schedule> key_provider(
        bootstrap_key);
    CKKSDenseBootstrapWithKeyProviderTimed<Schedule>(res, ct, key_provider,
                                                     timings);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapFromLevel(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const CKKSDenseBootstrapKey<Schedule> &bootstrap_key)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<Schedule> key_provider(
        bootstrap_key);
    CKKSDenseBootstrapFromLevelWithKeyProvider<Schedule>(res, ct, key_provider);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapFromLevelTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const CKKSDenseBootstrapKey<Schedule> &bootstrap_key,
    CKKSDenseBootstrapTimings &timings)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<Schedule> key_provider(
        bootstrap_key);
    CKKSDenseBootstrapFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, key_provider, timings);
}

template <class Schedule>
inline void CKKSDenseBootstrapHybridGiant(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSDenseBootstrapHybridGiantKey<Schedule> &bootstrap_key)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<
        Schedule, CKKSDenseBootstrapHybridGiantKey<Schedule>>
        key_provider(bootstrap_key);
    CKKSDenseBootstrapWithKeyProvider<Schedule>(res, ct, key_provider);
}

template <class Schedule>
inline void CKKSDenseBootstrapHybridGiantTimed(
    typename Schedule::OutputCiphertext &res,
    const typename Schedule::InputCiphertext &ct,
    const CKKSDenseBootstrapHybridGiantKey<Schedule> &bootstrap_key,
    CKKSDenseBootstrapTimings &timings)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<
        Schedule, CKKSDenseBootstrapHybridGiantKey<Schedule>>
        key_provider(bootstrap_key);
    CKKSDenseBootstrapWithKeyProviderTimed<Schedule>(res, ct, key_provider,
                                                     timings);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapHybridGiantFromLevel(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const CKKSDenseBootstrapHybridGiantKey<Schedule> &bootstrap_key)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<
        Schedule, CKKSDenseBootstrapHybridGiantKey<Schedule>>
        key_provider(bootstrap_key);
    CKKSDenseBootstrapFromLevelWithKeyProvider<Schedule>(res, ct, key_provider);
}

template <class Schedule, std::uint32_t InLogQ, std::uint32_t InLogDelta>
inline void CKKSDenseBootstrapHybridGiantFromLevelTimed(
    typename Schedule::OutputCiphertext &res,
    const CKKSCiphertext<typename Schedule::Param, InLogQ, InLogDelta> &ct,
    const CKKSDenseBootstrapHybridGiantKey<Schedule> &bootstrap_key,
    CKKSDenseBootstrapTimings &timings)
{
    const CKKSDenseBootstrapInMemoryKeyProvider<
        Schedule, CKKSDenseBootstrapHybridGiantKey<Schedule>>
        key_provider(bootstrap_key);
    CKKSDenseBootstrapFromLevelWithKeyProviderTimed<Schedule>(
        res, ct, key_provider, timings);
}

using lvl6CKKSDenseBootstrapInput =
    typename lvl6CKKSDenseBootstrapSchedule::InputCiphertext;
using lvl6CKKSDenseBootstrapOutput =
    typename lvl6CKKSDenseBootstrapSchedule::OutputCiphertext;
using lvl6CKKSDenseBootstrapSparseKey =
    CKKSDenseBootstrapKey<lvl6CKKSDenseBootstrapSchedule>;
using lvl6CKKSDenseBootstrapHybridGiantKey =
    CKKSDenseBootstrapHybridGiantKey<lvl6CKKSDenseBootstrapSchedule>;
using lvl6CKKSDenseBootstrapSparseFilesystemKeyProvider =
    CKKSDenseBootstrapFilesystemKeyProvider<lvl6CKKSDenseBootstrapSchedule>;
using lvl6CKKSDenseBootstrapHybridGiantFilesystemKeyProvider =
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapSchedule>;
using lvl6CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider =
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapSchedule>;
using lvl6CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider =
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapSchedule>;
using lvl6CKKSDenseBootstrapFastHybridGiantKey =
    CKKSDenseBootstrapHybridGiantKey<lvl6CKKSDenseBootstrapFastSchedule>;
using lvl6CKKSDenseBootstrapCompactHybridGiantKey =
    CKKSDenseBootstrapHybridGiantKey<lvl6CKKSDenseBootstrapCompactSchedule>;
using lvl6CKKSDenseBootstrapFastHybridGiantFilesystemKeyProvider =
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapFastSchedule>;
using lvl6CKKSDenseBootstrapCompactHybridGiantFilesystemKeyProvider =
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapCompactSchedule>;
using lvl6CKKSDenseBootstrapInverseInput =
    typename lvl6CKKSDenseBootstrapInverseSchedule::InputCiphertext;
using lvl6CKKSDenseBootstrapInverseOutput =
    typename lvl6CKKSDenseBootstrapInverseSchedule::OutputCiphertext;
using lvl6CKKSDenseBootstrapInverseHybridGiantKey =
    CKKSDenseBootstrapHybridGiantKey<lvl6CKKSDenseBootstrapInverseSchedule>;
using lvl6CKKSDenseBootstrapInverseHybridGiantFilesystemKeyProvider =
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapInverseSchedule>;
using lvl6CKKSDenseBootstrapInverseSeededHybridGiantFilesystemKeyProvider =
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapInverseSchedule>;
using lvl6CKKSDenseBootstrapInverseSeededHybridGiantStreamedFilesystemKeyProvider =
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapInverseSchedule>;
using lvl6CKKSDenseBootstrapTunedHybridGiantKey =
    CKKSDenseBootstrapHybridGiantKey<lvl6CKKSDenseBootstrapTunedSchedule>;
using lvl6CKKSDenseBootstrapTunedHybridGiantFilesystemKeyProvider =
    CKKSDenseBootstrapHybridGiantFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapTunedSchedule>;
using lvl6CKKSDenseBootstrapTunedSeededHybridGiantFilesystemKeyProvider =
    CKKSDenseBootstrapSeededHybridGiantFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapTunedSchedule>;
using lvl6CKKSDenseBootstrapTunedSeededHybridGiantStreamedFilesystemKeyProvider =
    CKKSDenseBootstrapSeededHybridGiantStreamedFilesystemKeyProvider<
        lvl6CKKSDenseBootstrapTunedSchedule>;

inline double CKKSPlainEvalModSineDegree5(double x)
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    const double x2 = x * x;
    return x * (1.0 - (two_pi * two_pi / 6.0) * x2 +
                (two_pi * two_pi * two_pi * two_pi / 120.0) * x2 * x2);
}

inline double CKKSPlainEvalModSineDegree3(double x)
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    return x * (1.0 - (two_pi * two_pi / 6.0) * x * x);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
struct CKKSEvalModSineDegree3Traits {
    static_assert(LogQ > 3 * LogDelta + CoeffLogDelta);
    static constexpr std::uint32_t x2_log_q = LogQ - LogDelta;
    static constexpr std::uint32_t x3_log_q = LogQ - 2 * LogDelta;
    static constexpr std::uint32_t term_input_log_q = x3_log_q;
    static constexpr std::uint32_t log_q = x3_log_q - CoeffLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
using CKKSEvalModSineDegree3Result =
    typename CKKSEvalModSineDegree3Traits<P, LogQ, LogDelta,
                                          CoeffLogDelta>::Ciphertext;

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
struct CKKSEvalModSineDegree3RelinKeys {
    using Traits =
        CKKSEvalModSineDegree3Traits<P, LogQ, LogDelta, CoeffLogDelta>;

    CKKSRelinKey<P, Traits::x2_log_q> x2{};
    CKKSRelinKey<P, Traits::x3_log_q> x3{};
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
inline void CKKSEvalModSineDegree3KeyGen(
    CKKSEvalModSineDegree3RelinKeys<P, LogQ, LogDelta, CoeffLogDelta> &keys,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    using Traits =
        CKKSEvalModSineDegree3Traits<P, LogQ, LogDelta, CoeffLogDelta>;
    keys.x2 = *makeCKKSRelinKey<P, Traits::x2_log_q>(key, noise);
    keys.x3 = *makeCKKSRelinKey<P, Traits::x3_log_q>(key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
inline void CKKSEvalModSineDegree3(
    CKKSEvalModSineDegree3Result<P, LogQ, LogDelta, CoeffLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSEvalModSineDegree3RelinKeys<P, LogQ, LogDelta, CoeffLogDelta>
        &keys)
{
    using Traits =
        CKKSEvalModSineDegree3Traits<P, LogQ, LogDelta, CoeffLogDelta>;
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    constexpr double c1 = -(two_pi * two_pi / 6.0);

    CKKSMultResult<P, LogQ, LogDelta, LogQ, LogDelta> x2;
    CKKSMult<P>(x2, ct, ct, keys.x2);

    CKKSMultResult<P, Traits::x2_log_q, LogDelta, LogQ, LogDelta> x3;
    CKKSMult<P>(x3, x2, ct, keys.x3);

    CKKSCiphertext<P, Traits::term_input_log_q, LogDelta> x_term_input;
    CKKSLevelReduce<P, LogQ, Traits::term_input_log_q, LogDelta>(
        x_term_input, ct);

    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x_term;
    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x3_term;
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x_term, x_term_input, 1.0);
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x3_term, x3, c1);

    CKKSAdd<P, Traits::log_q, LogDelta>(res, x_term, x3_term);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
struct CKKSEvalModSineDegree5Traits {
    static_assert(LogQ > 4 * LogDelta + CoeffLogDelta);
    static constexpr std::uint32_t x2_log_q = LogQ - LogDelta;
    static constexpr std::uint32_t x3_log_q = LogQ - 2 * LogDelta;
    static constexpr std::uint32_t x5_log_q = LogQ - 3 * LogDelta;
    static constexpr std::uint32_t term_input_log_q = x5_log_q;
    static constexpr std::uint32_t log_q = x5_log_q - CoeffLogDelta;
    static constexpr std::uint32_t log_delta = LogDelta;
    using Ciphertext = CKKSCiphertext<P, log_q, log_delta>;
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
using CKKSEvalModSineDegree5Result =
    typename CKKSEvalModSineDegree5Traits<P, LogQ, LogDelta,
                                          CoeffLogDelta>::Ciphertext;

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
struct CKKSEvalModSineDegree5RelinKeys {
    using Traits =
        CKKSEvalModSineDegree5Traits<P, LogQ, LogDelta, CoeffLogDelta>;

    CKKSRelinKey<P, Traits::x2_log_q> x2{};
    CKKSRelinKey<P, Traits::x3_log_q> x3{};
    CKKSRelinKey<P, Traits::x5_log_q> x5{};
};

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
inline void CKKSEvalModSineDegree5KeyGen(
    CKKSEvalModSineDegree5RelinKeys<P, LogQ, LogDelta, CoeffLogDelta> &keys,
    const Key<P> &key, CKKSNoise noise = {P::α, 0})
{
    using Traits =
        CKKSEvalModSineDegree5Traits<P, LogQ, LogDelta, CoeffLogDelta>;
    keys.x2 = *makeCKKSRelinKey<P, Traits::x2_log_q>(key, noise);
    keys.x3 = *makeCKKSRelinKey<P, Traits::x3_log_q>(key, noise);
    keys.x5 = *makeCKKSRelinKey<P, Traits::x5_log_q>(key, noise);
}

template <class P, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t CoeffLogDelta>
inline void CKKSEvalModSineDegree5(
    CKKSEvalModSineDegree5Result<P, LogQ, LogDelta, CoeffLogDelta> &res,
    const CKKSCiphertext<P, LogQ, LogDelta> &ct,
    const CKKSEvalModSineDegree5RelinKeys<P, LogQ, LogDelta, CoeffLogDelta>
        &keys)
{
    using Traits =
        CKKSEvalModSineDegree5Traits<P, LogQ, LogDelta, CoeffLogDelta>;
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    constexpr double c1 = -(two_pi * two_pi / 6.0);
    constexpr double c2 = two_pi * two_pi * two_pi * two_pi / 120.0;

    CKKSMultResult<P, LogQ, LogDelta, LogQ, LogDelta> x2;
    CKKSMult<P>(x2, ct, ct, keys.x2);

    CKKSMultResult<P, Traits::x2_log_q, LogDelta, LogQ, LogDelta> x3;
    CKKSMult<P>(x3, x2, ct, keys.x3);

    CKKSMultResult<P, Traits::x3_log_q, LogDelta, Traits::x2_log_q, LogDelta>
        x5;
    CKKSMult<P>(x5, x3, x2, keys.x5);

    CKKSCiphertext<P, Traits::term_input_log_q, LogDelta> x_term_input;
    CKKSCiphertext<P, Traits::term_input_log_q, LogDelta> x3_term_input;
    CKKSLevelReduce<P, LogQ, Traits::term_input_log_q, LogDelta>(
        x_term_input, ct);
    CKKSLevelReduce<P, Traits::x3_log_q, Traits::term_input_log_q, LogDelta>(
        x3_term_input, x3);

    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x_term;
    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x3_term;
    CKKSPlainMulResult<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>
        x5_term;
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x_term, x_term_input, 1.0);
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x3_term, x3_term_input, c1);
    CKKSPlainMulByReal<P, Traits::term_input_log_q, LogDelta, CoeffLogDelta>(
        x5_term, x5, c2);

    CKKSAdd<P, Traits::log_q, LogDelta>(res, x_term, x3_term);
    CKKSAddInPlace<P, Traits::log_q, LogDelta>(res, x5_term);
}

}  // namespace TFHEpp
