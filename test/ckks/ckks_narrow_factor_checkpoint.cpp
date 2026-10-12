#include <ckks/narrow_selector_file.hpp>
#include <ckks/narrow_refresh.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

namespace {
template <std::uint32_t Nbit>
struct FactorParam : TFHEpp::lvl5param {
    static constexpr std::uint32_t nbit = Nbit, n = 1U << Nbit;
    static constexpr std::uint32_t Bgbit = 4, B̅gbit = 4, l̅ = T::digits / 4;
    static const inline double α = std::ldexp(1.0, -448);
    static constexpr double ckks_min_noise_stddev = 3.2;
};
using P = FactorParam<14>;
using Schedule = TFHEpp::CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, false, 28, true>;

template <class S, std::size_t DigitBits = 0, class FactorKey = TFHEpp::CKKSNarrowBootstrapKey<S>>
void checkGroup(bool tail, std::size_t group, const TFHEpp::Key<typename S::Param> &dense,
                const TFHEpp::Key<typename S::Param> &sparse,
                const TFHEpp::CKKSNarrowSelectorFile<S> &source,
                const FactorKey &key,
                const TFHEpp::CKKSPreparedAutoKeyFFT<typename S::Param, S::factor_key_log_q> *factor_fft,
                std::size_t columns = 12, std::size_t chunk = 40, std::size_t alternate_chunk = 17)
{
    using namespace TFHEpp;
    using Schedule = S;
    using P = typename Schedule::Param;
    using T = typename P::T;
    const auto positions = narrow_detail::structuredPositions<Schedule>(sparse);
    std::vector<typename Schedule::Selector> selectors;
    const auto load_start = std::chrono::steady_clock::now();
    source.load(tail, group, selectors);
    const double load_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - load_start).count();
    std::mt19937_64 random(20262119 + group + 100 * tail);
    std::uniform_real_distribution<double> value(-0.5, 0.5);
    std::vector<CKKSSlotVector<P>> tables(columns);
    auto normalized = std::make_unique<std::array<double, P::n>>();
    for (std::size_t column = 0; column < tables.size(); ++column) {
        for (auto &coefficient : *normalized) coefficient = value(random);
        tables[column] = tail ? narrow_detail::tailTable<Schedule>(*normalized, group, column + 1)
                              : narrow_detail::headTable<Schedule>(*normalized, group, column + 1);
    }
    CKKSNarrowFactorTimings timings;
    const auto factors = narrow_detail::packedFactorsBatch<Schedule, DigitBits>(tables, selectors, {}, 0, key, chunk,
        factor_fft, &timings);
    CKKSNarrowFactorTimings baseline_timings;
    const auto baseline = narrow_detail::packedFactorsBatch<Schedule,16>(tables, selectors, {}, 0, key, chunk,
        factor_fft, &baseline_timings);
    for (std::size_t column = 0; column < tables.size(); ++column)
        if (factors[column].ct != baseline[column].ct)
            throw std::runtime_error("wide NTT limbs changed full-ring factor ciphertext");
    const auto alternate = narrow_detail::packedFactorsBatch<Schedule>(tables, selectors, {}, 0, key, alternate_chunk);
    for (std::size_t column = 0; column < tables.size(); ++column)
        if (factors[column].ct != alternate[column].ct)
            throw std::runtime_error("full-ring selector chunk sizes changed factor ciphertext");
    for (const std::size_t column : {std::size_t{0}, (columns - 1) / 2, columns - 1}) {
        const auto single = narrow_detail::packedFactorsBatch<Schedule>(
            {&tables[column], 1}, selectors, {}, 0, key, chunk);
        if (single.front().ct != factors[column].ct)
            throw std::runtime_error("full-ring batched factor differs from single-column factor");
    }
    auto decoded = std::make_unique<CKKSSlotVector<P>>();
    auto expected = std::make_unique<CKKSSlotVector<P>>();
    auto rotated = std::make_unique<CKKSSlotVector<P>>();
    double maximum = 0;
    for (std::size_t column = 0; column < tables.size(); ++column) {
        expected->fill({});
        const auto count = tail ? std::min(Schedule::tail_blocks_per_pack,
            Schedule::blocks - group * Schedule::tail_blocks_per_pack) : std::size_t{1};
        for (std::size_t ell = 0; ell < count; ++ell) {
            const auto block = tail ? group * Schedule::tail_blocks_per_pack + ell : group;
            const auto position = positions[block];
            rotateCKKSSlotVector<P>(*rotated, tables[column], -position);
            const bool negative = sparse[block * Schedule::block_width + position] == T{} - T{1};
            const auto begin = tail ? (2 * ell + 1) * Schedule::tail_width : 0;
            const auto end = begin + (tail ? Schedule::tail_width : Schedule::head_slots);
            for (std::size_t i = begin; i < end; ++i)
                (*expected)[i] = negative ? std::conj((*rotated)[i]) : (*rotated)[i];
        }
        ckksSlotDecrypt<P>(*decoded, factors[column], dense);
        for (std::size_t i = 0; i < Schedule::slots; ++i) {
            const auto error = std::abs((*decoded)[i] - (*expected)[i]);
            if (!std::isfinite(error)) throw std::runtime_error("nonfinite full-ring factor error");
            maximum = std::max(maximum, error);
        }
    }
    if (maximum > std::ldexp(1.0, -23)) throw std::runtime_error("full-ring packed factor failed fixed 23-bit gate");
    std::cout << std::setprecision(17) << "{\"tail\":" << (tail ? "true" : "false")
              << ",\"group\":" << group << ",\"columns\":" << columns << ",\"slots_checked\":" << columns * Schedule::slots
              << ",\"max_error\":" << maximum << ",\"chunk_size\":" << chunk
              << ",\"alternate_chunk_size\":" << alternate_chunk << ",\"chunks_exact\":true"
              << ",\"single_columns_first_middle_last_exact\":true,\"accuracy_passed\":true"
              << ",\"ntt_digit_bits\":" << CKKSNarrowNTTSelector<Schedule, DigitBits>::digit_bits
              << ",\"ntt_digits\":" << CKKSNarrowNTTSelector<Schedule, DigitBits>::digits
              << ",\"baseline16_ciphertexts_exact\":true"
              << ",\"baseline16_timings\":{\"selector_ntt_seconds\":" << baseline_timings.selector_ntt_seconds
              << ",\"matrix_seconds\":" << baseline_timings.matrix_seconds
              << ",\"inverse_crt_seconds\":" << baseline_timings.inverse_crt_seconds
              << ",\"conjugation_seconds\":" << baseline_timings.conjugation_seconds << "}"
              << ",\"timings\":{\"selector_load_seconds\":" << load_seconds
              << ",\"table_ntt_seconds\":" << timings.table_ntt_seconds
              << ",\"selector_ntt_seconds\":" << timings.selector_ntt_seconds
              << ",\"matrix_seconds\":" << timings.matrix_seconds
              << ",\"inverse_crt_seconds\":" << timings.inverse_crt_seconds
              << ",\"conjugation_seconds\":" << timings.conjugation_seconds << "}}" << std::endl;
}

int run(const char *checkpoint, const char *selector_directory)
{
    using namespace TFHEpp;
    std::ifstream file(checkpoint, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open checkpoint");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t, 9> metadata{};
    std::string function;
    archive(metadata, function);
    const std::array<std::uint64_t, 9> expected{
        0x4e41525243503031ULL, P::n, 85, 44, 385, 48, 37, 108, 448};
    if (metadata != expected || function != "sigmoid") throw std::runtime_error("unsupported factor checkpoint profile");
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    archive(*dense, *sparse);
    file.close();
    for (const auto value : *dense)
        if (value != typename P::T{} && value != typename P::T{1} && value != typename P::T{} - typename P::T{1})
            throw std::runtime_error("nonternary factor checkpoint dense secret");
    narrow_detail::structuredPositions<Schedule>(*sparse);
    CKKSNarrowSelectorFile<Schedule> source(selector_directory,
        typename CKKSNarrowSelectorFile<Schedule>::Existing{});
    // No checkpoint evaluation keys are loaded. The normal key container's
    // unused tree/rotation fields remain zero; only conjugation is exercised.
    auto key = std::make_unique<CKKSNarrowBootstrapKey<Schedule>>();
    std::cout << "generating full-ring factor conjugation key" << std::endl;
    CKKSAutoKeyGen<P, Schedule::factor_key_log_q>(key->factor_galois[P::nbit],
        2 * P::n - 1, *dense, {P::α, 0});
    key->factor_conjugation_fft = makeCKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q>(
        key->factor_galois[P::nbit]);
    checkGroup<Schedule>(false, Schedule::blocks - 1, *dense, *sparse, source, *key, key->factor_conjugation_fft.get());
    checkGroup<Schedule>(true, Schedule::tail_packs - 1, *dense, *sparse, source, *key, key->factor_conjugation_fft.get());
    return 0;
}
// Retained auxiliary8 diagnostic: load the original compatible evaluation key.
// The legacy plain37 mode above remains available for older checkpoints.
template <std::uint32_t Nbit = 14, std::size_t DigitBits = 0>
int runAuxiliary8(const std::string &state, const std::string &selector_directory)
{
    using namespace TFHEpp;
    using P = FactorParam<Nbit>;
    using S = CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, true, 8, true>;
    using R = CKKSNarrowRefreshSchedule<S, 44, 37, 0, true>;
    std::ifstream file(state, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open retained auxiliary8 state");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t,10> metadata{};
    std::string function;
    bool finished = false;
    archive(metadata, function, finished);
    const std::array<std::uint64_t,10> expected{
        0x4155583850523031ULL, P::n, 85, 44, 385, 48, 37, 108, 8, 448};
    if (metadata != expected || function != "sigmoid" || !finished)
        throw std::runtime_error("incompatible retained auxiliary8 state");
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto input = std::make_unique<typename R::Input>();
    auto saved = std::make_unique<typename R::Output>();
    auto key = std::make_unique<CKKSNarrowRefreshKey<R>>();
    auto square = std::make_unique<CKKSCompactDDRelinKey<P,108>>();
    std::cout << "loading retained auxiliary8 evaluation keys; no new keys" << std::endl;
    archive(*dense, *sparse, *slots, *input, *saved, *key, *square);
    if (file.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("unexpected retained state suffix");
    file.close();
    for (const auto &value : *dense)
        if (value != typename P::T{} && value != typename P::T{1} && value != typename P::T{} - typename P::T{1})
            throw std::runtime_error("nonternary retained dense key");
    narrow_detail::structuredPositions<S>(*sparse);
    CKKSNarrowSelectorFile<S> source(selector_directory, typename CKKSNarrowSelectorFile<S>::Existing{});
    key->narrow.factor_conjugation_fft = makeCKKSPreparedAutoKeyFFT<P,S::factor_key_log_q>(
        key->narrow.factor_galois[P::nbit]);
    checkGroup<S, DigitBits>(false, S::blocks - 1, *dense, *sparse, source, key->narrow, key->narrow.factor_conjugation_fft.get());
    checkGroup<S, DigitBits>(true, S::tail_packs - 1, *dense, *sparse, source, key->narrow, key->narrow.factor_conjugation_fft.get());
    return 0;
}

template <class Chain, std::size_t... I>
constexpr std::size_t chainArchiveBytes(std::index_sequence<I...>)
{
    using Keys = decltype(std::declval<Chain>().keys);
    return ((7 * sizeof(std::uint64_t) + sizeof(std::tuple_element_t<I, Keys>::data)) + ...);
}

// The checked square suffix, DD header and cache mode anchor the fixed-size
// conjugation payload. Seek with the original archive to retain its byte order.
template <std::uint32_t Nbit>
int runAuxiliary8Bounded(const std::string &state, const std::string &selector_directory)
{
    using namespace TFHEpp;
    using P = FactorParam<Nbit>;
    using S = CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, true, 8, true>;
    using Chain = typename S::RelinKeys;
    using Keys = decltype(std::declval<Chain>().keys);
    using Square = CKKSCompactDDRelinKey<P, 108>;
    using Auto = CKKSAutoKey<P, S::factor_key_log_q>;
    std::ifstream file(state, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open retained factor state");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t, 10> metadata{};
    std::string function;
    bool finished = false;
    archive(metadata, function, finished);
    const std::array<std::uint64_t, 10> expected{
        0x4155583850523031ULL, P::n, 85, 44, 385, 48, 37, 108, 8, 448};
    if (metadata != expected || function != "sigmoid" || !finished)
        throw std::runtime_error("incompatible retained factor state");
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    archive(*dense, *sparse);
    for (const auto &value : *dense)
        if (value != typename P::T{} && value != typename P::T{1} && value != typename P::T{} - typename P::T{1})
            throw std::runtime_error("nonternary retained factor secret");
    narrow_detail::structuredPositions<S>(*sparse);
    const auto prefix_end = file.tellg();
    file.seekg(0, std::ios::end);
    const auto end = file.tellg();
    constexpr auto square_bytes = 7 * sizeof(std::uint64_t) + sizeof(Square::data);
    constexpr auto chain_bytes = chainArchiveBytes<Chain>(std::make_index_sequence<std::tuple_size_v<Keys>>{});
    const auto square_begin = end - static_cast<std::streamoff>(square_bytes);
    const auto chain_begin = square_begin - static_cast<std::streamoff>(sizeof(std::uint32_t) + chain_bytes);
    const auto auto_end = chain_begin - static_cast<std::streamoff>(sizeof(CKKSAutoKey<P, S::output_key_log_q>::data));
    const auto auto_begin = auto_end - static_cast<std::streamoff>(sizeof(Auto::data));
    if (auto_begin < prefix_end) throw std::runtime_error("truncated retained factor archive");
    std::array<std::uint64_t, 7> anchor{};
    file.seekg(square_begin);
    archive(anchor);
    const std::array<std::uint64_t, 7> square_anchor{0x434b4b5343444431ULL, P::n, 1, 108, 4, 4, 1};
    if (anchor != square_anchor) throw std::runtime_error("invalid retained square anchor");
    file.seekg(chain_begin);
    archive(anchor);
    const std::array<std::uint64_t, 7> chain_anchor{0x434b4b5343444431ULL, P::n, 1, 348, 4, 4, 1};
    if (anchor != chain_anchor) throw std::runtime_error("invalid retained factor DD anchor");
    file.seekg(square_begin - static_cast<std::streamoff>(sizeof(std::uint32_t)));
    std::uint32_t mode{};
    archive(mode);
    if (mode > 2 || file.tellg() != square_begin)
        throw std::runtime_error("invalid retained factor cache mode");
    auto raw = std::make_unique<Auto>();
    file.seekg(auto_begin);
    archive(*raw);
    if (file.tellg() != auto_end) throw std::runtime_error("invalid retained factor key length");
    file.close();
    auto prepared = makeCKKSPreparedAutoKeyFFT<P, S::factor_key_log_q>(*raw);
    raw.reset();
    CKKSNarrowSelectorFile<S> source(selector_directory, typename CKKSNarrowSelectorFile<S>::Existing{});
    std::cout << "loaded only retained factor conjugation key; no new keys" << std::endl;
    checkGroup<S, 60>(false, S::blocks - 1, *dense, *sparse, source, *prepared, prepared.get(), 5, 8, 5);
    checkGroup<S, 60>(true, S::tail_packs - 1, *dense, *sparse, source, *prepared, prepared.get(), 5, 8, 5);
    return 0;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc == 1) {
        std::cout << "Usage: ckks_narrow_factor_checkpoint [--retained-aux8|--retained-aux8-wide60|--retained-small-aux8-wide60] checkpoint.bin selector-directory\n";
        return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "--retained-aux8") {
        try { return runAuxiliary8<>(argv[2], argv[3]); }
        catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
    }
    if (argc == 4 && std::string(argv[1]) == "--retained-small-aux8-wide60") {
        try { return runAuxiliary8<9, 60>(argv[2], argv[3]); }
        catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
    }
    if (argc == 4 && std::string(argv[1]) == "--retained-aux8-wide60") {
        try { return runAuxiliary8<14, 60>(argv[2], argv[3]); }
        catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
    }
    if (argc == 4 && std::string(argv[1]) == "--retained-aux8-wide60-bounded") {
        try { return runAuxiliary8Bounded<14>(argv[2], argv[3]); }
        catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
    }
    if (argc == 4 && std::string(argv[1]) == "--retained-small-aux8-wide60-bounded") {
        try { return runAuxiliary8Bounded<9>(argv[2], argv[3]); }
        catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
    }
    if (argc != 3) return 2;
    try { return run(argv[1], argv[2]); }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
}
