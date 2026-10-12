#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <fstream>
#include <iostream>

namespace {
template <std::uint32_t Nbit>
struct BasisParam {
    static constexpr std::int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = Nbit, n = 1U << Nbit, k = 1;
    static constexpr std::uint32_t l = 1, l_a = 1, lₐ = 1;
    static constexpr std::uint32_t Bgbit = 18, Bg_abit = 4, Bgₐbit = 4;
    using T = TFHEpp::MultiLimbUInt<9>;
    static constexpr T Bg = T{1} << Bgbit, Bg_a = T{1} << Bg_abit, Bgₐ = Bg_a;
    static constexpr TFHEpp::ErrorDistribution errordist = TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::ldexp(1.0, -576);
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << 253;
    static constexpr std::uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0;
    static constexpr std::uint32_t l̅ = 72, l̅ₐ = l̅, B̅gbit = 8, B̅gₐbit = 8;
};
using P = BasisParam<9>;
using Narrow = TFHEpp::CKKSNarrowBootstrapSchedule<P, 48, 518, 47, 39, 31, true, 44, true>;
using Refresh = TFHEpp::CKKSNarrowRefreshSchedule<Narrow, 44, 37, 0, true>;
using Basis = TFHEpp::CKKSNarrowBasisRefreshSchedule<Refresh, 31>;

struct GroupSource : TFHEpp::CKKSNarrowSelectorSource<Narrow> {
    std::shared_ptr<const TFHEpp::CKKSNarrowSelectorFile<Narrow>> raw;
    explicit GroupSource(decltype(raw) source) : raw(std::move(source)) {}
    std::size_t size(bool tail) const override { return raw->size(tail); }
    void load(bool tail, std::size_t group, std::vector<Selector> &result) const override
    { raw->load(tail, group, result); }
};

struct ChunkSource : GroupSource {
    using GroupSource::GroupSource;
    mutable std::size_t calls = 0, maximum_count = 0;
    bool short_chunk = false;
    bool prefetch = false;
    void load(bool, std::size_t, std::vector<Selector> &) const override
    { throw std::runtime_error("streamed evaluation loaded a full raw group"); }
    bool hasRawChunks() const override { return true; }
    bool supportsRawChunkPrefetch() const override { return prefetch; }
    void loadChunk(bool tail, std::size_t group, std::size_t first, std::size_t count,
                   std::vector<Selector> &result) const override
    {
        if (count > 8) throw std::runtime_error("raw chunk exceeded requested bound");
        ++calls; maximum_count = std::max(maximum_count, count);
        raw->loadChunk(tail, group, first, count, result);
        if (short_chunk && !result.empty()) result.pop_back();
    }
};

void run(const std::filesystem::path &directory)
{
    using namespace TFHEpp;
    std::ifstream file(directory / "checkpoint.bin", std::ios::binary);
    if (!file) throw std::runtime_error("cannot open retained small basis checkpoint");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t, 12> metadata{};
    archive(metadata);
    const std::array<std::uint64_t, 12> expected{
        0x4e41525242433031ULL, 512, 576, 518, 47, 31, 85, 44, 18, 8, 44, 37};
    if (metadata != expected) throw std::runtime_error("unsupported retained basis profile");
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto input = std::make_unique<typename Refresh::Input>();
    auto key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    auto basis = std::make_unique<CKKSNarrowFrequencyBasisKey<Basis>>();
    auto square = std::make_unique<CKKSCompactDDRelinKey<P, Basis::output_log_q>>();
    archive(*dense, *sparse, *slots, *input, *key, *basis, *square);
    if (file.peek() != std::char_traits<char>::eof()) throw std::runtime_error("unexpected retained basis suffix");
    file.close();
    narrow_detail::structuredPositions<Narrow>(*sparse);
    auto raw_source = std::make_shared<CKKSNarrowSelectorFile<Narrow>>(
        directory / "selectors", CKKSNarrowSelectorFile<Narrow>::Existing{});
    key->narrow.selector_source = std::make_shared<GroupSource>(raw_source);
    CKKSNarrowUseStagedFFT(key->narrow, 8);
    std::vector<typename Refresh::Coefficients> inputs(3);
    CKKSNarrowPrepareCoefficients<Refresh>(inputs[0], *input, *key);
    narrow_detail::imaginaryCoefficientHalf<P>(inputs[1], inputs[0]);
    inputs[2] = inputs[0];
    inputs[2].ct[1][0] = ckks_detail::reduceToLevel<P, 48>(inputs[2].ct[1][0] + (P::T{1} << 44));
    const CKKSNarrowFrequencyBasis plan(6, {1, 3, 5, 6});
    const std::vector<std::complex<double>> coefficients{
        {0.125, 0}, {0.0625, -0.03125}, {-0.03125, 0.015625},
        {0.015625, 0.0078125}, {-0.0078125, -0.00390625}, {0.00390625, 0}, {0, 0.001953125}};
    std::size_t groups = 0;
    key->narrow.factor_observer = [&](bool, std::size_t, const auto &) { ++groups; };
    std::vector<typename Basis::Output> expected_outputs(inputs.size());
    for (std::size_t i = 0; i < inputs.size(); ++i)
        CKKSNarrowEvalSeriesWithBasis<Basis>(expected_outputs[i], inputs[i], key->narrow,
            *basis, plan, coefficients, 4);
    if (groups != 3 * (Narrow::blocks + Narrow::tail_packs))
        throw std::runtime_error("invalid per-input reference selector count");
    auto chunk_source = std::make_shared<ChunkSource>(raw_source);
    key->narrow.selector_source = chunk_source;
    for (const std::size_t tile : {5, 7, 12}) {
        groups = 0;
        std::vector<typename Basis::Output> actual;
        CKKSNarrowEvalSeriesWithBasisBatch<Basis, 44>(actual, inputs, key->narrow,
            *basis, plan, coefficients, tile);
        const auto expected_groups = ((12 + tile - 1) / tile) * (Narrow::blocks + Narrow::tail_packs);
        if (groups != expected_groups || actual.size() != expected_outputs.size())
            throw std::runtime_error("cross-input basis tile count mismatch");
        for (std::size_t i = 0; i < actual.size(); ++i)
            if (actual[i].ct != expected_outputs[i].ct)
                throw std::runtime_error("cross-input basis tile changed ciphertext or owner");
        std::cout << "{\"tile_size\":" << tile << ",\"inputs\":3,\"frequencies_per_input\":4"
                  << ",\"factor_groups\":" << groups << ",\"expected_factor_groups\":" << expected_groups
                  << ",\"exact_coefficients\":" << 3 * 2 * P::n << ",\"passed\":true}" << std::endl;
    }
    groups = 0;
    std::vector<typename Basis::Output> empty;
    CKKSNarrowEvalSeriesWithBasisBatch<Basis, 44>(empty, {}, key->narrow, *basis, plan, coefficients, 5);
    if (!empty.empty() || groups) throw std::runtime_error("empty basis batch evaluated selectors");

    const std::array<CKKSNarrowFrequencyBasis, 3> plans{
        plan, CKKSNarrowFrequencyBasis(40, {1, 3, 4, 9, 11, 16, 17, 19, 20}), plan};
    const auto frequencies = CKKSNarrowFrequencyBasisUnion(plans);
    if (frequencies != std::vector<std::uint32_t>{1, 3, 4, 5, 6, 9, 11, 16, 17, 19, 20})
        throw std::runtime_error("incorrect shared frequency union");
    std::vector<CKKSNarrowBlindRotationRequest<Narrow, 44>> requests;
    for (const auto &coefficient_input : inputs)
        for (const auto frequency : frequencies)
            requests.push_back({&coefficient_input, frequency, {1, 0}});
    groups = 0;
    std::vector<typename Narrow::Output> reference_roots;
    CKKSNarrowBlindRotateBatch<Narrow, 44>(reference_roots, requests, key->narrow, frequencies.size());
    if (groups != inputs.size() * (Narrow::blocks + Narrow::tail_packs))
        throw std::runtime_error("incorrect per-input union reference pass count");
    for (const std::size_t tile : {22, 33}) {
        groups = 0;
        std::vector<std::map<std::uint32_t, typename Narrow::Output>> actual;
        CKKSNarrowBlindRotateBasisPlansBatch<Narrow, 44>(actual, inputs, plans, key->narrow, tile);
        const auto expected_groups = ((requests.size() + tile - 1) / tile) * (Narrow::blocks + Narrow::tail_packs);
        if (groups != expected_groups || actual.size() != inputs.size())
            throw std::runtime_error("incorrect union input ownership or pass count");
        for (std::size_t input_index = 0; input_index < inputs.size(); ++input_index) {
            if (actual[input_index].size() != frequencies.size())
                throw std::runtime_error("duplicate plan introduced duplicate roots");
            for (std::size_t f = 0; f < frequencies.size(); ++f)
                if (actual[input_index].at(frequencies[f]).ct != reference_roots[input_index * frequencies.size() + f].ct)
                    throw std::runtime_error("union batching changed a root or its input owner");
        }
        std::cout << "{\"union_tile_size\":" << tile << ",\"inputs\":3,\"unique_frequencies\":11"
                  << ",\"factor_groups\":" << groups << ",\"expected_factor_groups\":" << expected_groups
                  << ",\"exact_root_coefficients\":" << requests.size() * 2 * P::n << ",\"passed\":true}" << std::endl;
    }
    groups = 0;
    std::vector<std::map<std::uint32_t, typename Narrow::Output>> union_empty;
    CKKSNarrowBlindRotateBasisPlansBatch<Narrow, 44>(union_empty, {}, plans, key->narrow, 22);
    if (!union_empty.empty() || groups) throw std::runtime_error("empty union inputs evaluated selectors");
    CKKSNarrowBlindRotateBasisPlansBatch<Narrow, 44>(union_empty, inputs, {}, key->narrow, 22);
    if (union_empty.size() != inputs.size() || groups ||
        std::any_of(union_empty.begin(), union_empty.end(), [](const auto &roots) { return !roots.empty(); }))
        throw std::runtime_error("empty union plans evaluated selectors");
    try {
        CKKSNarrowBlindRotateBasisPlansBatch<Narrow, 44>(union_empty, {}, {}, key->narrow, 0);
        throw std::runtime_error("zero union tile accepted");
    }
    catch (const std::invalid_argument &) {}
    if (!chunk_source->calls || chunk_source->maximum_count != 8)
        throw std::runtime_error("raw chunk path was not exercised");
    std::cout << "{\"raw_chunk_calls\":" << chunk_source->calls
              << ",\"maximum_raw_chunk_count\":" << chunk_source->maximum_count
              << ",\"full_group_load_forbidden\":true,\"passed\":true}" << std::endl;
    chunk_source->prefetch = true;
    chunk_source->calls = chunk_source->maximum_count = 0;
    {
        std::vector<typename Basis::Output> actual;
        CKKSNarrowEvalSeriesWithBasisBatch<Basis, 44>(actual, inputs, key->narrow,
            *basis, plan, coefficients, 12);
        if (actual.size() != expected_outputs.size())
            throw std::runtime_error("prefetched basis output count mismatch");
        for (std::size_t i = 0; i < actual.size(); ++i)
            if (actual[i].ct != expected_outputs[i].ct)
                throw std::runtime_error("prefetched basis changed ciphertext or owner");
        if (!chunk_source->calls || chunk_source->maximum_count != 8)
            throw std::runtime_error("prefetched raw chunk bound was not exercised");
        std::cout << "{\"prefetched_exact_coefficients\":" << 3 * 2 * P::n
                  << ",\"raw_chunk_calls\":" << chunk_source->calls
                  << ",\"maximum_raw_chunk_count\":8,\"full_group_load_forbidden\":true,\"passed\":true}"
                  << std::endl;
    }
    chunk_source->short_chunk = true;
    for (const bool prefetch : {false, true}) {
        chunk_source->prefetch = prefetch;
        try {
            std::vector<typename Narrow::Output> rejected;
            CKKSNarrowBlindRotateBatch<Narrow, 44>(rejected,
                std::span<const CKKSNarrowBlindRotationRequest<Narrow, 44>>(requests).first(1), key->narrow, 1);
            throw std::runtime_error("short raw selector chunk accepted");
        }
        catch (const std::invalid_argument &error) {
            if (std::string(error.what()) != "raw selector chunk count mismatch") throw;
        }
        std::cout << "{\"short_raw_chunk_rejected\":true,\"prefetch\":"
                  << (prefetch ? "true" : "false") << ",\"passed\":true}" << std::endl;
    }
}

template <class P>
void fixed16Product(TFHEpp::Polynomial<P> &out, const TFHEpp::Polynomial<P> &mask,
                    const TFHEpp::Polynomial<P> &secret)
{
    using namespace TFHEpp;
    static const modular_ntt::NegacyclicNTTPlan plan(P::n, {4611686018427322369ULL, 7});
    const auto modulus = plan.modulus();
    std::vector<std::uint64_t> lhs(P::n), rhs(P::n);
    for (std::size_t i = 0; i < P::n; ++i) {
        const auto value = multilimb_to_signed_i64(secret[i]);
        rhs[i] = value < 0 ? modulus - (-value) : value;
    }
    plan.forward(rhs); out.fill({});
    for (std::size_t shift = 0; shift < 518; shift += 16) {
        for (std::size_t i = 0; i < P::n; ++i)
            lhs[i] = static_cast<std::uint64_t>(ckks_detail::reduceToLevel<P,518>(mask[i]) >> shift) & 65535;
        plan.forward(lhs);
        for (std::size_t i = 0; i < P::n; ++i) lhs[i] = modular_ntt::multiply(lhs[i], rhs[i], modulus);
        plan.inverse(lhs);
        for (std::size_t i = 0; i < P::n; ++i) {
            const auto value = lhs[i] <= modulus / 2 ? static_cast<std::int64_t>(lhs[i])
                                                    : -static_cast<std::int64_t>(modulus - lhs[i]);
            out[i] += P::T::from_signed_i64(value) << shift;
        }
    }
    for (auto &value : out) value = ckks_detail::reduceToLevel<P,518>(value);
}

template <std::uint32_t Nbit>
void phaseSample(const std::filesystem::path &directory)
{
    using namespace TFHEpp;
    using P = BasisParam<Nbit>;
    using S = CKKSNarrowBootstrapSchedule<P,48,518,47,39,31,true,44,true>;
    std::ifstream file(directory / "checkpoint.bin",std::ios::binary);
    if (!file) throw std::runtime_error("cannot open selector sample checkpoint");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t,12> metadata{};archive(metadata);
    const std::array<std::uint64_t,12> expected{0x4e41525242433031ULL,P::n,576,518,47,31,85,44,18,8,44,37};
    if (metadata != expected) throw std::runtime_error("invalid selector sample profile");
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    archive(*dense,*sparse);file.close();
    for (const auto &value : *dense)
        if (value != typename P::T{} && value != typename P::T{1} && value != typename P::T{} - typename P::T{1})
            throw std::runtime_error("nonternary selector sample secret");
    const auto positions = narrow_detail::structuredPositions<S>(*sparse);
    auto selector = std::make_unique<typename S::Selector>();
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto encoded = std::make_unique<Polynomial<P>>(), product = std::make_unique<Polynomial<P>>();
    auto reference = std::make_unique<Polynomial<P>>();
    for (const bool tail : {false,true})
        for (std::size_t group = 0; group < (tail ? S::tail_packs : S::blocks); ++group) {
            std::vector<std::size_t> records;
            const auto begin = tail ? group * S::tail_blocks_per_pack : group;
            const auto count = tail ? std::min(S::tail_blocks_per_pack,S::blocks - begin) : 1;
            for (std::size_t ell = 0; ell < count; ++ell) records.push_back(positions[begin + ell]);
            std::sort(records.begin(),records.end());records.erase(std::unique(records.begin(),records.end()),records.end());
            std::size_t zero = 0;
            while (std::binary_search(records.begin(),records.end(),zero)) ++zero;
            records.push_back(zero);
            const auto path = directory / "selectors" / ((tail ? "tail-" : "head-") + std::to_string(group) + ".bin");
            std::ifstream source(path,std::ios::binary);
            if (!source) throw std::runtime_error("missing sampled selector group");
            cereal::PortableBinaryInputArchive reader(source);
            std::array<std::uint64_t,9> header{};reader(header);
            const std::array<std::uint64_t,9> expected_header{0x4e4152524f573031ULL,P::n,518,47,39,S::block_width,576,tail,group};
            if (header != expected_header || std::filesystem::file_size(path) != 73 + S::block_width * sizeof(*selector))
                throw std::runtime_error("invalid sampled selector header or length");
            std::uint64_t maximum_noise = 0;
            for (const auto record : records) {
                source.seekg(73 + static_cast<std::streamoff>(record * sizeof(*selector)));
                reader(*selector);slots->fill({});
                for (std::size_t ell = 0; ell < count; ++ell) {
                    const auto block = begin + ell;
                    if (record != static_cast<std::size_t>(positions[block])) continue;
                    const auto positive = (*sparse)[block * S::block_width + record] == typename P::T{1};
                    const auto first = tail ? (2 * ell + 1) * S::tail_width : 0;
                    const auto length = tail ? S::tail_width : S::head_slots;
                    std::fill_n(slots->begin()+first,length,positive ? std::complex<double>{1,0} : std::complex<double>{0,1});
                }
                ckksSlotEncode<P,518,47>(*encoded,*slots);
                fixed16Product<P>(*reference,selector->ct[0],*dense);
                ckks_detail::exactActiveTorusMulBySmall<P,518>(*product,selector->ct[0],*dense);
                if (*product != *reference) throw std::runtime_error("sampled selector adaptive product mismatch");
                for (std::size_t i = 0; i < P::n; ++i) {
                    const auto error = ckks_detail::centeredLevelToTorus<P,518>(selector->ct[1][i]-(*reference)[i]-(*encoded)[i]);
                    const auto signed_error = multilimb_to_signed_i64(error);
                    if (P::T::from_signed_i64(signed_error) != error || signed_error < -39 || signed_error > 39)
                        throw std::runtime_error("sampled selector phase failed 12-sigma noise bound");
                    maximum_noise = std::max(maximum_noise,static_cast<std::uint64_t>(std::abs(signed_error)));
                }
            }
            std::cout << "{\"n\":" << P::n << ",\"tail\":" << (tail ? "true" : "false")
                      << ",\"group\":" << group << ",\"selectors_checked\":" << records.size()
                      << ",\"phase_coefficients_checked\":" << records.size()*P::n
                      << ",\"maximum_noise\":" << maximum_noise << ",\"passed\":true}" << std::endl;
        }
}
}

int main(int argc, char **argv)
{
    if (argc != 2 && argc != 3) { std::cerr << "Usage: ckks_narrow_basis_batch_checkpoint [--phase-sample] CHECKPOINT_DIRECTORY\n"; return 2; }
    try {
        if (argc == 2) run(argv[1]);
        else if (std::string(argv[1]) == "--phase-sample") {
            std::ifstream file(std::filesystem::path(argv[2])/"checkpoint.bin",std::ios::binary);
            cereal::PortableBinaryInputArchive archive(file);
            std::array<std::uint64_t,12> metadata{};archive(metadata);
            if (metadata[1] == 512) phaseSample<9>(argv[2]);
            else if (metadata[1] == 32768) phaseSample<15>(argv[2]);
            else throw std::runtime_error("unsupported selector sample degree");
        }
        else throw std::runtime_error("unknown basis checkpoint mode");
    }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
}
