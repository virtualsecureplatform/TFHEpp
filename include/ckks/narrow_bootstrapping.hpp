#pragma once

#include "ckks.hpp"
#include "prepared_auto.hpp"
#include "prepared_plain.hpp"
#include "prepared_dd.hpp"
#include "compact_dd.hpp"
#include "narrow_matrix.hpp"
#include <functional>
#include <chrono>
#include <future>
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Ouyang et al., ePrint 2026/2119, Sections 3 and 4. The source key has
// exactly one signed nonzero coefficient per public block; the evaluation key may
// be independent. All product-tree relinearizations use Double Decomposition.
template <class P, std::uint32_t InputLogQ, std::uint32_t SelectorLogQ,
          std::uint32_t LogDelta, std::uint32_t PlainLogDelta,
          std::size_t HammingWeight, bool RelinBeforeRescale = false, std::uint32_t AuxiliaryLogScale = 0,
          bool CompactDD = false>
struct CKKSNarrowBootstrapSchedule {
    using Param = P;
    static_assert(P::k == 1);
    static_assert(P::nbit >= 3 && P::nbit <= 16 && P::n == (1U << P::nbit));
    static_assert(InputLogQ > 0 && SelectorLogQ <= std::numeric_limits<typename P::T>::digits);
    static_assert(PlainLogDelta > 0 && PlainLogDelta <= 50);
    static_assert(HammingWeight > 1 && HammingWeight < P::n);
    static constexpr std::size_t blocks = HammingWeight;
    static constexpr std::size_t slots = P::n / 2;
    // Section 8 uses h=31, which does not divide power-of-two ring degrees.
    // Use ceil(N/h)-wide blocks and a shorter final block; its nonexistent
    // candidate positions have zero encrypted selectors. The same head/tail
    // windows apply unchanged to every valid position.
    static constexpr std::size_t block_width = (P::n + blocks - 1) / blocks;
    static_assert((blocks - 1) * block_width < P::n);
    static constexpr std::size_t block_size(std::size_t x)
    {
        return std::min(block_width, P::n - x * block_width);
    }
    static_assert(block_width > 1 && 2 * (block_width - 1) <= slots);
    static constexpr std::size_t tail_width = block_width - 1;
    static constexpr std::size_t head_slots = slots - tail_width;
    static constexpr std::size_t tail_blocks_per_pack = slots / (2 * tail_width);
    static constexpr std::size_t tail_packs =
        (blocks + tail_blocks_per_pack - 1) / tail_blocks_per_pack;
    static constexpr std::size_t selector_count = block_width * (blocks + tail_packs);
    static constexpr std::uint32_t input_log_q = InputLogQ;
    static constexpr std::uint32_t selector_log_q = SelectorLogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    static constexpr std::uint32_t plain_log_delta = PlainLogDelta;
    static constexpr std::uint32_t factor_log_q = SelectorLogQ - PlainLogDelta;
    // Include the public beta as a transparent leaf, as in SHIP Algorithm 1.
    // For the paper's h=31, its 32-leaf tree still consumes five levels and
    // avoids an extra rescale after the tree for multiplying beta.
    static constexpr std::uint32_t tree_depth = ckks_detail::power_tree_depth(blocks + 1);
    static_assert(SelectorLogQ > PlainLogDelta + tree_depth * LogDelta + LogDelta);
    static constexpr std::uint32_t tree_log_q = factor_log_q - tree_depth * LogDelta;
    static constexpr std::uint32_t output_log_q = tree_log_q;
    // Optional public scale lift for automorphisms. Data/tree levels stay
    // fixed; security tuning must count these larger evaluation-key levels.
    static constexpr std::uint32_t auxiliary_log_scale = AuxiliaryLogScale;
    static constexpr std::uint32_t factor_key_log_q = factor_log_q + AuxiliaryLogScale;
    static constexpr std::uint32_t output_key_log_q = output_log_q + AuxiliaryLogScale;
    static_assert(factor_key_log_q <= std::numeric_limits<typename P::T>::digits);
    using Selector = CKKSCiphertext<P, SelectorLogQ, LogDelta>;
    using Factor = CKKSCiphertext<P, factor_log_q, LogDelta>;
    using TreeResult = CKKSCiphertext<P, tree_log_q, LogDelta>;
    using Output = CKKSCiphertext<P, output_log_q, LogDelta>;
    static constexpr bool relinearize_before_rescale = RelinBeforeRescale;
    // Existing DD chains index keys at StartQ-(level+1)*delta. Shift their
    // starting metadata by one scale to switch at each product's input Q.
    static constexpr std::uint32_t relin_chain_start_log_q = factor_log_q +
        (RelinBeforeRescale ? LogDelta : 0);
    static constexpr bool compact_dd = CompactDD;
    using RelinKeys = std::conditional_t<CompactDD,
        CKKSCompactDDRelinKeyChain<P, relin_chain_start_log_q, LogDelta, tree_depth>,
        CKKSDDRelinKeyChain<P, relin_chain_start_log_q, LogDelta, tree_depth>>;
    // Secondary digits represent existing encrypted rows; widening the
    // derived FFT cache does not change the primary gadget or saved keys.
    static constexpr std::uint32_t prepared_relin_bbarbit = CompactDD &&
        P::B̅gbit < 16 && std::numeric_limits<typename P::T>::digits % 16 == 0 &&
        P::Bgbit + 16 + P::nbit + std::bit_width(
            CKKSDDRelinPrimaryRowCountForLevel<P, relin_chain_start_log_q - LogDelta, P::Bgbit>()) +
            3 < std::numeric_limits<double>::digits ? 16 : P::B̅gbit;
    using PreparedRelinKeys = CKKSPreparedDDRelinKeyChainFFT<P, relin_chain_start_log_q,
        LogDelta, tree_depth, P::Bgbit, prepared_relin_bbarbit>;
};

template <class Schedule>
inline auto CKKSNarrowMakePreparedRelinearization(const typename Schedule::RelinKeys &keys)
{
    if constexpr (Schedule::compact_dd &&
        Schedule::prepared_relin_bbarbit != Schedule::Param::B̅gbit)
        return makeCKKSPreparedDDRelinKeyChainFFTWithSecondaryBase<
            Schedule::prepared_relin_bbarbit>(keys);
    else return makeCKKSPreparedDDRelinKeyChainFFT(keys);
}

template <class Schedule, std::size_t DigitBits = 0>
struct CKKSNarrowNTTSelector {
    // Large rings benefit from fewer integer-convolution limbs. Select wide
    // limbs only when the full worst-case sum fits the existing two-prime CRT.
    // Division tests avoid overflowing the bound while screening candidates.
    static constexpr unsigned __int128 crt_limit =
        static_cast<unsigned __int128>(modular_ntt::degree65536_primes[0].value) *
        modular_ntt::degree65536_primes[1].value / 2;
    static constexpr unsigned __int128 table_bound =
        (static_cast<unsigned __int128>(1) << Schedule::plain_log_delta) + 1;
    static constexpr bool fits(std::size_t bits) {
        if (bits == 0 || bits > 60) return false;
        return table_bound <= (crt_limit - 1) / Schedule::Param::n /
            Schedule::block_width / ((std::uint64_t{1} << bits) - 1);
    }
    // A derived schedule may select scratch limbs independently of the
    // serialized ciphertext/key parameters. Explicit per-call widths win.
    static constexpr std::size_t configured_digit_bits = [] {
        if constexpr (requires { Schedule::selector_matrix_digit_bits; })
            return std::size_t{Schedule::selector_matrix_digit_bits};
        else return std::size_t{0};
    }();
    static constexpr std::size_t digit_bits = DigitBits != 0 ? DigitBits :
        configured_digit_bits != 0 ? configured_digit_bits :
        Schedule::Param::n >= 4096 ? (fits(48) ? 48 : fits(32) ? 32 : 16) : 16;
    // Explicit wider widths are useful for profiling; the default remains
    // the qualified 48/32/16-bit policy until full evaluation is validated.
    static_assert(digit_bits > 0 && digit_bits <= 60);
    static_assert(fits(digit_bits), "narrow selector limbs exceed exact CRT range");
    static constexpr std::uint64_t digit_mask = (std::uint64_t{1} << digit_bits) - 1;
    static constexpr std::size_t digits = (Schedule::selector_log_q + digit_bits - 1) / digit_bits;
    // Each prime: [component][unsigned decomposition digit][natural NTT index].
    // Unit-norm transformed root tables encode coefficients of magnitude
    // <= 2^PlainLogDelta + 1. Use one prime when this worst-case convolution
    // bound fits its centered range; otherwise retain the two-prime CRT.
    static constexpr unsigned __int128 convolution_bound =
        ((static_cast<unsigned __int128>(1) << Schedule::plain_log_delta) + 1) *
        digit_mask * Schedule::Param::n * Schedule::block_width;
    static constexpr std::size_t prime_count = convolution_bound <
        modular_ntt::degree65536_primes[0].value / 2 ? 1 : 2;
    std::array<std::vector<std::uint64_t>, prime_count> spectra;
};

// Optional diagnostic counters for a completed factor batch. They are owned
// by the caller and never serialized or used in ciphertext evaluation.
struct CKKSNarrowFactorTimings {
    // Exclusive consumer time: synchronous reads or uncovered prefetch waits.
    // Background I/O overlapping matrix work is not counted a second time.
    double selector_load_seconds = 0;
    double table_ntt_seconds = 0;
    double selector_ntt_seconds = 0;
    double matrix_seconds = 0;
    double inverse_crt_seconds = 0;
    double conjugation_seconds = 0;
};

// Exclusive wall-clock scopes for one root tile. The leaf scope includes
// inline tail alignment when bounded batch rotations are disabled.
struct CKKSNarrowRootTimings {
    double input_and_factor_cache_seconds = 0;
    double leaf_factors_and_inline_alignment_seconds = 0;
    double bounded_tail_alignment_seconds = 0;
    double dd_cache_preparation_seconds = 0;
    double dd_trees_and_beta_seconds = 0;
};

enum class CKKSNarrowCacheMode { Prepared, Streamed, Reference };

// Optional public selector storage. Chunk-capable sources bound raw selector
// memory without changing the ciphertexts or exact NTT arithmetic.
template <class Schedule>
struct CKKSNarrowSelectorSource {
    using Selector = typename Schedule::Selector;
    virtual ~CKKSNarrowSelectorSource() = default;
    virtual std::size_t size(bool tail) const = 0;
    virtual void load(bool tail, std::size_t group, std::vector<Selector> &result) const = 0;
    virtual bool hasRawChunks() const { return false; }
    // Opt in only when loadChunk may run on a worker thread. At most one
    // read is outstanding, and evaluation joins it before leaving the group.
    virtual bool supportsRawChunkPrefetch() const { return false; }
    virtual void loadChunk(bool, std::size_t, std::size_t, std::size_t,
                           std::vector<Selector> &) const
    { throw std::logic_error("selector source has no raw chunks"); }
    virtual bool hasPreparedChunks() const { return false; }
    virtual void loadPreparedChunk(bool, std::size_t, std::size_t, std::size_t,
                                   std::vector<CKKSNarrowNTTSelector<Schedule>> &) const
    { throw std::logic_error("selector source has no prepared chunks"); }
};

template <class Schedule> struct CKKSNarrowBootstrapKey;

// External factor construction can preserve a retained selector scale until
// after the weighted matrix sum. Like selector storage, this is runtime state.
template <class Schedule>
struct CKKSNarrowFactorSource {
    using P = typename Schedule::Param;
    virtual ~CKKSNarrowFactorSource() = default;
    virtual std::size_t size(bool tail) const = 0;
    virtual void adjustBeta(CKKSSlotVector<P> &, CKKSSlotVector<P> &) const {}
    virtual std::vector<typename Schedule::Factor> evaluate(
        std::span<const CKKSSlotVector<P>> tables, bool tail, std::size_t group,
        const CKKSNarrowBootstrapKey<Schedule> &key,
        const CKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q> *cache,
        CKKSNarrowFactorTimings *timings) const = 0;
};

template <class Schedule>
struct CKKSNarrowBootstrapKey {
    using P = typename Schedule::Param;
    std::vector<typename Schedule::Selector> head;
    std::vector<typename Schedule::Selector> tail;
    // External storage is not serialized. Reattach/regenerate it after load.
    std::shared_ptr<const CKKSNarrowSelectorSource<Schedule>> selector_source;
    std::shared_ptr<const CKKSNarrowFactorSource<Schedule>> factor_source;
    CKKSGaloisKey<P, Schedule::factor_key_log_q> factor_galois;
    CKKSAutoKey<P, Schedule::output_key_log_q> output_conjugation;
    typename Schedule::RelinKeys relin;
    // Publicly derived runtime backend; never serialized with raw keys.
    std::function<void(TRLWE<P>&, const TRLWE<P>&, uint32_t, uint32_t)>
        automorphism_evaluator;
    // Optional shared-key batches. Runtime state is excluded from archives.
    std::function<void(std::span<TRLWE<P>>,std::span<const TRLWE<P>>,uint32_t,uint32_t)>
        automorphism_batch_evaluator;
    // Optional prepared backend; derived state is never archived.
    std::function<void(typename Schedule::TreeResult &,
        std::vector<typename Schedule::Factor>, const typename Schedule::RelinKeys &)>
        product_tree_evaluator;
    std::function<void(std::vector<typename Schedule::TreeResult>&,
        std::vector<std::vector<typename Schedule::Factor>>,const typename Schedule::RelinKeys&)>
        product_tree_batch_evaluator;
    bool batch_product_trees = false;
    std::unique_ptr<typename Schedule::PreparedRelinKeys> relin_fft;
    static constexpr std::uint32_t projected_relin_primary_bits = 8 % P::Bgbit == 0 ? 8 : P::Bgbit;
    using ProjectedRelinKeys = CKKSPreparedDDRelinKeyChainFFT<P,
        Schedule::relin_chain_start_log_q, Schedule::log_delta, Schedule::tree_depth,
        projected_relin_primary_bits, Schedule::prepared_relin_bbarbit>;
    std::unique_ptr<ProjectedRelinKeys> relin_projected_fft;
    // Derived Section-5 setup cache; rebuild after loading a serialized key.
    std::vector<CKKSNarrowNTTSelector<Schedule>> head_ntt, tail_ntt;
    std::unique_ptr<CKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q>> factor_conjugation_fft;
    std::unique_ptr<CKKSPreparedAutoKeyFFT<P, Schedule::output_key_log_q>> output_conjugation_fft;
    static constexpr std::uint32_t projected_primary_bits = 8 % P::B̅gbit == 0 ? 8 : P::B̅gbit;
    using ProjectedFactorAuto = CKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q, 16, projected_primary_bits>;
    using ProjectedOutputAuto = CKKSPreparedAutoKeyFFT<P, Schedule::output_key_log_q, 16, projected_primary_bits>;
    std::unique_ptr<ProjectedFactorAuto> factor_conjugation_projected_fft;
    std::unique_ptr<ProjectedOutputAuto> output_conjugation_projected_fft;
    std::unique_ptr<CKKSPreparedRotationKeyFFT<P, Schedule::factor_key_log_q>> tail_rotations_fft;
    // Runtime cache policy, omitted from the legacy archive like derived keys.
    bool bounded_tail_rotation_fft = false;
    // Opt-in parallel calls require a thread-safe public evaluator.
    bool parallel_public_tail_alignment = false;
    // Opt in only with thread-safe beta adjustment and product-tree backends.
    bool parallel_root_columns = false;
    std::size_t selector_ntt_chunk_size = 0;
    bool staged_fft_caches = false;
    // Optional caller-owned diagnostic observer, excluded from key archives.
    std::function<void(bool, std::size_t, const CKKSNarrowFactorTimings &)> factor_observer;
    std::function<void(const CKKSNarrowRootTimings &)> root_stage_observer;
    CKKSNarrowCacheMode cache_mode = CKKSNarrowCacheMode::Prepared;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(head, tail, factor_galois, output_conjugation, relin, cache_mode);
        if constexpr (Archive::is_loading::value) {
            selector_source.reset();
            factor_source.reset();
            product_tree_evaluator = {};
            product_tree_batch_evaluator = {};
            batch_product_trees = false;
            automorphism_evaluator = {};
            automorphism_batch_evaluator = {};
            relin_fft.reset();
            relin_projected_fft.reset();
            factor_conjugation_fft.reset();
            output_conjugation_fft.reset();
            factor_conjugation_projected_fft.reset();
            output_conjugation_projected_fft.reset();
            tail_rotations_fft.reset();
            bounded_tail_rotation_fft = false;
            parallel_public_tail_alignment = false;
            parallel_root_columns = false;
            selector_ntt_chunk_size = 0;
            staged_fft_caches = false;
            factor_observer = {};
            root_stage_observer = {};
            head_ntt.clear();
            tail_ntt.clear();
        }
    }
};

// Tail pack extraction uses public shifts 2*(ell+1)*tail_width modulo N/2.
template <class Schedule>
inline void CKKSNarrowPrepareTailRotations(CKKSNarrowBootstrapKey<Schedule> &key)
{
    using P = typename Schedule::Param;
    if constexpr (use_multilimb_digit_fft_v<P>) {
        std::uint32_t bits = 0;
        for (std::size_t ell = 0;
             ell < std::min(Schedule::blocks, Schedule::tail_blocks_per_pack); ell++)
            bits |= (2 * (ell + 1) * Schedule::tail_width) % Schedule::slots;
        key.tail_rotations_fft = makeCKKSPreparedRotationKeyFFT<P,
            Schedule::factor_key_log_q>(key.factor_galois, bits);
        key.bounded_tail_rotation_fft = false;
        key.staged_fft_caches = false;
    }
    else throw std::invalid_argument("prepared tail rotations require a digit FFT backend");
}

template <class Schedule>
inline void CKKSNarrowUseBoundedTailRotations(CKKSNarrowBootstrapKey<Schedule> &key)
{
    if constexpr (use_multilimb_digit_fft_v<typename Schedule::Param>) {
        key.tail_rotations_fft.reset();
        key.bounded_tail_rotation_fft = true;
    }
    else throw std::invalid_argument("bounded tail rotations require a digit FFT backend");
}

template <class Schedule>
inline void CKKSNarrowUseSelectorChunks(CKKSNarrowBootstrapKey<Schedule> &key,
    std::size_t chunk_size)
{
    if (chunk_size == 0) throw std::invalid_argument("selector NTT chunk size must be positive");
    key.selector_ntt_chunk_size = std::min(chunk_size, Schedule::block_width);
    key.head_ntt.clear();
    key.tail_ntt.clear();
    key.cache_mode = CKKSNarrowCacheMode::Streamed;
}

template <class Schedule>
inline void CKKSNarrowUseStagedFFT(CKKSNarrowBootstrapKey<Schedule> &key,
    std::size_t selector_chunk_size = 40)
{
    if constexpr (use_multilimb_digit_fft_v<typename Schedule::Param>) {
        CKKSNarrowUseSelectorChunks<Schedule>(key, selector_chunk_size);
        CKKSNarrowUseBoundedTailRotations<Schedule>(key);
        key.factor_conjugation_fft.reset();
        key.output_conjugation_fft.reset();
        key.factor_conjugation_projected_fft.reset();
        key.output_conjugation_projected_fft.reset();
        key.relin_fft.reset();
        key.relin_projected_fft.reset();
        key.product_tree_evaluator = {};
        key.product_tree_batch_evaluator = {};
        key.batch_product_trees = false;
        key.automorphism_evaluator = {};
        key.automorphism_batch_evaluator = {};
        key.staged_fft_caches = true;
    }
    else throw std::invalid_argument("staged FFT caches require a digit FFT backend");
}

template <class Schedule>
inline void CKKSNarrowPrepareRelinearization(CKKSNarrowBootstrapKey<Schedule> &key)
{
    if constexpr (use_multilimb_digit_fft_v<typename Schedule::Param>) {
        auto prepared = CKKSNarrowMakePreparedRelinearization<Schedule>(key.relin);
        key.relin_fft = std::move(prepared);
        key.product_tree_evaluator = {};
        key.product_tree_batch_evaluator = {};
        key.batch_product_trees = false;
        key.relin_projected_fft.reset();
        key.staged_fft_caches = false;
    }
    else throw std::invalid_argument("prepared DD relinearization requires a digit FFT backend");
}

// Optional retained-key primary8 DD cache. Wider primary digits change noise;
// this leaves the serialized keys and all other cache policies unchanged.
template <class Schedule>
inline void CKKSNarrowPrepareProjectedRelinearization(CKKSNarrowBootstrapKey<Schedule> &key)
{
    using P = typename Schedule::Param;
    static_assert(Schedule::compact_dd && use_multilimb_digit_fft_v<P>);
    static_assert(P::Bgbit <= 8 && 8 % P::Bgbit == 0);
    auto projected = makeCKKSCompactDDRelinKeyChainWithPrimaryBase<8>(key.relin);
    auto prepared = makeCKKSPreparedDDRelinKeyChainFFTWithSecondaryBase<
        Schedule::prepared_relin_bbarbit>(*projected);
    key.product_tree_evaluator = {};
    key.product_tree_batch_evaluator = {};
    key.batch_product_trees = false;
    key.relin_projected_fft = std::move(prepared);
    key.relin_fft.reset();
}

// Optional derived FFT caches. Rebuild after changing raw keys or loading.
template <class Schedule>
inline void CKKSNarrowPrepareAutomorphisms(CKKSNarrowBootstrapKey<Schedule> &key)
{
    using P = typename Schedule::Param;
    if constexpr (use_multilimb_digit_fft_v<P>) {
        auto factor = makeCKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q>(key.factor_galois[P::nbit]);
        auto output = makeCKKSPreparedAutoKeyFFT<P, Schedule::output_key_log_q>(key.output_conjugation);
        key.automorphism_evaluator = {};
        key.automorphism_batch_evaluator = {};
        key.factor_conjugation_fft = std::move(factor);
        key.output_conjugation_fft = std::move(output);
        key.factor_conjugation_projected_fft.reset();
        key.output_conjugation_projected_fft.reset();
        key.staged_fft_caches = false;
    }
    else throw std::invalid_argument("prepared automorphisms require a digit FFT backend");
}

// Optional primary8 conjugation caches for schedules with an auxiliary lift.
// Raw keys, ciphertext parameters and DD relinearization remain unchanged.
template <class Schedule>
inline void CKKSNarrowPrepareProjectedAutomorphisms(CKKSNarrowBootstrapKey<Schedule> &key)
{
    using P = typename Schedule::Param;
    if constexpr (use_multilimb_digit_fft_v<P> && P::B̅gbit <= 8 && 8 % P::B̅gbit == 0 && Schedule::auxiliary_log_scale >= 8) {
        auto factor = makeCKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q, 16, 8>(key.factor_galois[P::nbit]);
        auto output = makeCKKSPreparedAutoKeyFFT<P, Schedule::output_key_log_q, 16, 8>(key.output_conjugation);
        key.automorphism_evaluator = {};
        key.automorphism_batch_evaluator = {};
        key.factor_conjugation_projected_fft = std::move(factor);
        key.output_conjugation_projected_fft = std::move(output);
        key.factor_conjugation_fft.reset();
        key.output_conjugation_fft.reset();
    }
    else throw std::invalid_argument("projected conjugation requires a digit FFT backend, primary width at most8, and at least8 auxiliary bits");
}

namespace narrow_detail {
template <class Schedule, bool Factor, class AutoKey>
inline void conjugateWithKey(TRLWE<typename Schedule::Param> &result,
    const TRLWE<typename Schedule::Param> &input, const AutoKey &key)
{
    using P = typename Schedule::Param;
    constexpr auto level = Factor ? Schedule::factor_log_q : Schedule::output_log_q;
    if constexpr (Schedule::auxiliary_log_scale != 0)
        ckks_detail::applyWithAuxiliaryScale<P, level, Schedule::auxiliary_log_scale>(result, input,
            [&](auto &out, const auto &in) {
                CKKSEvalAuto<P, level + Schedule::auxiliary_log_scale>(out, in, 2 * P::n - 1, key);
            });
    else CKKSEvalAuto<P, level>(result, input, 2 * P::n - 1, key);
}

// Factor evaluation can reuse a standalone prepared conjugation key without
// allocating the unrelated rotation and product-tree keys.
template <class Schedule, bool Factor>
inline void conjugate(TRLWE<typename Schedule::Param> &result,
    const TRLWE<typename Schedule::Param> &input,
    const CKKSPreparedAutoKeyFFT<typename Schedule::Param,
        Factor ? Schedule::factor_key_log_q : Schedule::output_key_log_q> &key,
    const CKKSPreparedAutoKeyFFT<typename Schedule::Param,
        Factor ? Schedule::factor_key_log_q : Schedule::output_key_log_q> *override_cache = nullptr)
{
    conjugateWithKey<Schedule, Factor>(result, input, override_cache ? *override_cache : key);
}

template <class Schedule, bool Factor, class FactorKey = CKKSNarrowBootstrapKey<Schedule>>
inline void conjugate(TRLWE<typename Schedule::Param> &result,
    const TRLWE<typename Schedule::Param> &input, const FactorKey &key,
    const CKKSPreparedAutoKeyFFT<typename Schedule::Param,
        Factor ? Schedule::factor_key_log_q : Schedule::output_key_log_q> *override_cache = nullptr)
{
    using P = typename Schedule::Param;
    constexpr auto level = Factor ? Schedule::factor_log_q : Schedule::output_log_q;
    auto apply = [&](const auto &auto_key) {
        conjugateWithKey<Schedule, Factor>(result, input, auto_key);
    };
    if constexpr (requires { key.automorphism_evaluator; })
        if (key.automorphism_evaluator) {
            key.automorphism_evaluator(result, input, level, 2*P::n-1);
            return;
        }
    if constexpr (use_multilimb_digit_fft_v<P>) {
        if (override_cache) { apply(*override_cache); return; }
        const auto &projected = [&]() -> const auto & {
            if constexpr (Factor) return key.factor_conjugation_projected_fft;
            else return key.output_conjugation_projected_fft;
        }();
        if (projected) { apply(*projected); return; }
        const auto &cache = [&]() -> const auto & {
            if constexpr (Factor) return key.factor_conjugation_fft;
            else return key.output_conjugation_fft;
        }();
        if (cache) { apply(*cache); return; }
        if (key.staged_fft_caches) {
            const auto &raw = [&]() -> const auto & {
                if constexpr (Factor) return key.factor_galois[P::nbit];
                else return key.output_conjugation;
            }();
            auto temporary = makeCKKSPreparedAutoKeyFFT<P, level + Schedule::auxiliary_log_scale>(raw);
            apply(*temporary);
            return;
        }
    }
    if constexpr (Factor) apply(key.factor_galois[P::nbit]);
    else apply(key.output_conjugation);
}
template <class Schedule, class GaloisKey>
inline void rotate(TRLWE<typename Schedule::Param> &result,
    const TRLWE<typename Schedule::Param> &input, int steps, const GaloisKey &key)
{
    using P = typename Schedule::Param;
    if constexpr (Schedule::auxiliary_log_scale != 0)
        ckks_detail::applyWithAuxiliaryScale<P, Schedule::factor_log_q, Schedule::auxiliary_log_scale>(result, input,
            [&](auto &out, const auto &in) { CKKSRotateSlots<P, Schedule::factor_key_log_q>(out, in, steps, key); });
    else CKKSRotateSlots<P, Schedule::factor_log_q>(result, input, steps, key);
}

}  // namespace narrow_detail

// Section 6's four coefficient-packed indicator ciphertexts. Substitution
// keys are shared by every decompression tree, and are setup data rather
// than part of the online selector payload. Division by N and slot masking
// consume nbit + PlainLogDelta modulus bits during reconstruction.
template <class Schedule>
struct CKKSNarrowPackedKey {
    using P = typename Schedule::Param;
    static constexpr std::uint32_t log_q = Schedule::selector_log_q +
        Schedule::plain_log_delta + P::nbit;
    static_assert(log_q <= std::numeric_limits<typename P::T>::digits);
    using Packed = CKKSCiphertext<P, log_q, Schedule::log_delta>;
    std::array<Packed, 4> indicators;
    std::array<CKKSAutoKey<P, log_q>, P::nbit> substitutions;
    std::array<std::unique_ptr<CKKSPreparedAutoKeyFFT<P, log_q>>, P::nbit> substitutions_fft;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(indicators, substitutions);
        if constexpr (Archive::is_loading::value)
            for (auto &cache : substitutions_fft) cache.reset();
    }
};

// Reconstruction can use more accurate public masks without generating a new
// packed key. This schedule describes selectors only: the evaluator may rescale
// them before its product tree, so its tree budget is specified separately.
template <class PackedSchedule, std::uint32_t PlainLogDelta>
struct CKKSNarrowPackedReconstructionSchedule {
    using Param = typename PackedSchedule::Param;
    static_assert(PlainLogDelta > 0 && PlainLogDelta <= 50);
    static constexpr auto packed_log_q = CKKSNarrowPackedKey<PackedSchedule>::log_q;
    static_assert(packed_log_q > Param::nbit + PlainLogDelta);
    static constexpr auto input_log_q = PackedSchedule::input_log_q;
    static constexpr auto selector_log_q = packed_log_q - Param::nbit - PlainLogDelta;
    static constexpr auto log_delta = PackedSchedule::log_delta;
    static constexpr auto plain_log_delta = PlainLogDelta;
    static constexpr auto blocks = PackedSchedule::blocks;
    static constexpr auto slots = PackedSchedule::slots;
    static constexpr auto block_width = PackedSchedule::block_width;
    static constexpr auto head_slots = PackedSchedule::head_slots;
    static constexpr auto tail_width = PackedSchedule::tail_width;
    static constexpr auto tail_blocks_per_pack = PackedSchedule::tail_blocks_per_pack;
    static constexpr auto tail_packs = PackedSchedule::tail_packs;
    static constexpr auto selector_count = PackedSchedule::selector_count;
    using Selector = CKKSCiphertext<Param, selector_log_q, log_delta>;
};

// Cache the most frequently used depths of the coefficient extraction tree.
// A depth near the leaves is reused by exponentially more nodes. Raw keys
// remain available at uncached depths; zero levels clears the derived cache.
template <class Schedule>
inline void CKKSNarrowPreparePackedSubstitutions(CKKSNarrowPackedKey<Schedule> &key,
                                                std::size_t levels = 2)
{
    using P = typename Schedule::Param;
    if (levels > P::nbit) throw std::invalid_argument("invalid packed substitution cache depth");
    decltype(key.substitutions_fft) caches;
    if (levels != 0) {
        if constexpr (use_multilimb_digit_fft_v<P>) {
            for (std::size_t r = P::nbit - levels; r < P::nbit; r++)
                caches[r] = makeCKKSPreparedAutoKeyFFT<P,
                    CKKSNarrowPackedKey<Schedule>::log_q>(key.substitutions[r]);
        }
        else throw std::invalid_argument("prepared packed substitutions require a digit FFT backend");
    }
    key.substitutions_fft = std::move(caches);
}

// The budget covers FFT spectrum data only. Callers must separately reserve
// memory for raw keys, ciphertexts, allocator metadata and reconstruction.
template <class Schedule>
constexpr std::size_t CKKSNarrowPackedSubstitutionCacheBytesPerDepth()
{
    using P = typename Schedule::Param;
    if constexpr (use_multilimb_digit_fft_v<P>) {
        using Cache = CKKSPreparedAutoKeyFFT<P, CKKSNarrowPackedKey<Schedule>::log_q>;
        static_assert(Cache::spectrum_count <=
            std::numeric_limits<std::size_t>::max() / sizeof(PolynomialInFD<P>));
        return Cache::spectrum_count * sizeof(PolynomialInFD<P>);
    }
    else throw std::invalid_argument("prepared packed substitutions require a digit FFT backend");
}

template <class Schedule>
constexpr std::size_t CKKSNarrowPackedSubstitutionCacheDepthForBudget(std::size_t bytes)
{
    if (bytes == 0) return 0;
    return std::min<std::size_t>(Schedule::Param::nbit,
        bytes / CKKSNarrowPackedSubstitutionCacheBytesPerDepth<Schedule>());
}

template <class Schedule>
inline std::size_t CKKSNarrowPreparePackedSubstitutionsWithinBudget(
    CKKSNarrowPackedKey<Schedule> &key, std::size_t bytes)
{
    const auto levels = CKKSNarrowPackedSubstitutionCacheDepthForBudget<Schedule>(bytes);
    CKKSNarrowPreparePackedSubstitutions<Schedule>(key, levels);
    return levels;
}

namespace narrow_detail {

template <class Schedule>
inline std::array<int, Schedule::blocks> structuredPositions(
    const Key<typename Schedule::Param> &source)
{
    using P = typename Schedule::Param;
    std::array<int, Schedule::blocks> positions{};
    for (std::size_t x = 0; x < Schedule::blocks; x++) {
        std::size_t count = 0;
        for (std::size_t y = 0; y < Schedule::block_size(x); y++) {
            const auto value = source[x * Schedule::block_width + y];
            if (value == typename P::T{0}) continue;
            if (value != typename P::T{1} && value != -typename P::T{1})
                throw std::invalid_argument("narrow bootstrap requires a ternary source key");
            positions[x] = static_cast<int>(y);
            count++;
        }
        if (count != 1)
            throw std::invalid_argument("narrow bootstrap requires one nonzero per block");
    }
    return positions;
}

// a contains a_i/q, reduced to [-1/2,1/2). Negative indices wrap modulo N;
// negacyclic signs are supplied separately by Equation (2).
template <class P>
inline std::complex<double> root(const std::array<double, P::n> &a,
                                 std::int64_t index, int sign,
                                 std::uint32_t frequency)
{
    index = (index % P::n + P::n) % P::n;
    const double cycles = std::remainder(a[static_cast<std::size_t>(index)] *
                                            sign * frequency, 1.0);
    return std::polar(1.0, 2.0 * std::acos(-1.0) * cycles);
}

template <class Schedule>
inline CKKSSlotVector<typename Schedule::Param> headTable(
    const std::array<double, Schedule::Param::n> &a, std::size_t block,
    std::uint32_t frequency)
{
    using P = typename Schedule::Param;
    CKKSSlotVector<P> table;
    const std::int64_t j = block * Schedule::block_width;
    for (std::size_t i = 0; i < Schedule::slots; i++) {
        if (i < Schedule::head_slots)
            table[i] = root<P>(a, static_cast<std::int64_t>(i) - j,
                               j <= static_cast<std::int64_t>(i) ? 1 : -1, frequency);
        else
            table[i] = root<P>(a, static_cast<std::int64_t>(P::n / 2 + i) - j,
                               -1, frequency);
    }
    return table;
}

template <class Schedule>
inline CKKSSlotVector<typename Schedule::Param> tailTable(
    const std::array<double, Schedule::Param::n> &a, std::size_t pack,
    std::uint32_t frequency)
{
    using P = typename Schedule::Param;
    CKKSSlotVector<P> table{};
    const std::size_t begin = pack * Schedule::tail_blocks_per_pack;
    const std::size_t count = std::min(Schedule::tail_blocks_per_pack,
                                     Schedule::blocks - begin);
    for (std::size_t ell = 0; ell < count; ell++) {
        const std::int64_t j = (begin + ell) * Schedule::block_width;
        for (std::size_t t = 0; t < 2 * Schedule::tail_width; t++) {
            const std::int64_t i = Schedule::slots - 2 * Schedule::tail_width + t;
            table[2 * ell * Schedule::tail_width + t] =
                root<P>(a, i - j, j <= i ? 1 : -1, frequency);
        }
    }
    return table;
}

// Equations (19)--(22). Accumulate before rescaling and conjugate only once
// per completed factor. TFHEpp rotates left; the paper's Rot_y rotates right.
template <class P, std::uint32_t LogQ, std::size_t DigitBits>
inline std::uint64_t selectorDigit(const typename P::T &coefficient, std::size_t digit)
{
    static_assert(DigitBits > 0 && DigitBits < 64);
    static_assert(LogQ > 0 && LogQ <= ckks_detail::torus_width_v<P>);
    const auto offset = digit * DigitBits;
    if (offset >= LogQ) return 0;
    const auto bits = std::min<std::size_t>(DigitBits, LogQ - offset);
    const auto mask = (std::uint64_t{1} << bits) - 1;
    if constexpr (requires { P::T::limbs; coefficient.limb; }) {
        // A scratch digit spans at most two native limbs. Avoid copying,
        // masking and shifting the entire torus for every digit and prime.
        const auto limb = offset / 64;
        const auto shift = offset % 64;
        auto value = coefficient.limb[limb] >> shift;
        if (shift != 0 && limb + 1 < P::T::limbs)
            value |= coefficient.limb[limb + 1] << (64 - shift);
        return value & mask;
    }
    else return static_cast<std::uint64_t>(coefficient >> offset) & mask;
}

template <class Schedule, std::size_t DigitBits = 0>
inline void prepareSelectors(
    std::vector<CKKSNarrowNTTSelector<Schedule, DigitBits>> &prepared,
    std::span<const typename Schedule::Selector> selectors)
{
    using P = typename Schedule::Param;
    using Cache = CKKSNarrowNTTSelector<Schedule, DigitBits>;
    // The transform tables depend only on the public degree and primes.
    // Share their immutable tables across selectors, chunks and later runs
    // in this process instead of rebuilding them for every ciphertext.
    static const auto plans = [] {
        std::vector<modular_ntt::NegacyclicNTTPlan> result;
        result.reserve(Cache::prime_count);
        for (std::size_t prime = 0; prime < Cache::prime_count; prime++)
            result.emplace_back(P::n, modular_ntt::degree65536_primes[prime]);
        return result;
    }();
    prepared.resize(selectors.size());
#pragma omp parallel for schedule(static)
    for (std::size_t j = 0; j < selectors.size(); j++)
        for (std::size_t prime = 0; prime < Cache::prime_count; prime++) {
            auto &spectrum = prepared[j].spectra[prime];
            spectrum.resize((P::k + 1) * Cache::digits * P::n);
            const auto &plan = plans[prime];
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t digit = 0; digit < Cache::digits; digit++) {
                    auto row = std::span<std::uint64_t>(spectrum).subspan(
                        (c * Cache::digits + digit) * P::n, P::n);
                    for (std::size_t i = 0; i < P::n; i++) {
                        row[i] = selectorDigit<P, Schedule::selector_log_q,
                            Cache::digit_bits>(selectors[j].ct[c][i], digit);
                    }
                    plan.forward(row);
                }
        }
}

template <class Schedule, class FactorKey = CKKSNarrowBootstrapKey<Schedule>>
inline typename Schedule::Factor packedFactorReference(
    const CKKSSlotVector<typename Schedule::Param> &table,
    std::span<const typename Schedule::Selector> selectors,
    std::size_t begin, const FactorKey &key,
    const CKKSPreparedAutoKeyFFT<typename Schedule::Param, Schedule::factor_key_log_q> *factor_fft = nullptr)
{
    using P = typename Schedule::Param;
    auto accum_storage = std::make_unique<std::array<TRLWE<P>, 2>>();
    auto &accum1 = (*accum_storage)[0];
    auto &accum2 = (*accum_storage)[1];
    for (std::size_t y = 0; y < Schedule::block_width; y++) {
        CKKSSlotVector<P> rotated;
        rotateCKKSSlotVector<P>(rotated, table, -static_cast<int>(y));
        CKKSSlotVector<P> first, second;
        for (std::size_t i = 0; i < Schedule::slots; i++) {
            const auto z = rotated[i];
            first[i] = (z - std::complex<double>(0, 1) * std::conj(z)) / 2.0;
            second[i] = (std::conj(z) - std::complex<double>(0, 1) * z) / 2.0;
        }
        auto plain_storage = std::make_unique<std::array<Polynomial<P>, 3>>();
        auto &plain1 = (*plain_storage)[0];
        auto &plain2 = (*plain_storage)[1];
        auto &product = (*plain_storage)[2];
        ckksSlotEncode<P, Schedule::selector_log_q, Schedule::plain_log_delta>(plain1, first);
        ckksSlotEncode<P, Schedule::selector_log_q, Schedule::plain_log_delta>(plain2, second);
        for (std::size_t c = 0; c <= P::k; c++) {
            ckks_detail::polyMulTorusBySmallSigned<P, Schedule::selector_log_q>(
                product, selectors[begin + y].ct[c], plain1);
            for (std::size_t i = 0; i < P::n; i++) accum1[c][i] += product[i];
            ckks_detail::polyMulTorusBySmallSigned<P, Schedule::selector_log_q>(
                product, selectors[begin + y].ct[c], plain2);
            for (std::size_t i = 0; i < P::n; i++) accum2[c][i] += product[i];
        }
    }
    auto factor_storage = std::make_unique<std::array<typename Schedule::Factor, 3>>();
    auto &first = (*factor_storage)[0];
    auto &second = (*factor_storage)[1];
    auto &conjugated = (*factor_storage)[2];
    for (std::size_t c = 0; c <= P::k; c++)
        for (std::size_t i = 0; i < P::n; i++) {
            first.ct[c][i] = ckks_detail::roundedLevelRightShift<
                P, Schedule::selector_log_q, Schedule::plain_log_delta>(accum1[c][i]);
            second.ct[c][i] = ckks_detail::roundedLevelRightShift<
                P, Schedule::selector_log_q, Schedule::plain_log_delta>(accum2[c][i]);
        }
    conjugate<Schedule, true>(conjugated.ct, second.ct, key, factor_fft);
    CKKSAddInPlace<P>(first, conjugated);
    return first;
}

// Section 5, Equation (25), applied independently to the small limbs of the
// power-of-two torus. One or two auxiliary NTT primes recover exact signed integer
// convolutions before reduction modulo 2^Q. The prime product is an arithmetic
// scratch modulus; it does not change the ciphertext's cryptographic modulus.
template <class Schedule, std::size_t DigitBits = 0,
          class FactorKey = CKKSNarrowBootstrapKey<Schedule>>
inline std::vector<typename Schedule::Factor> packedFactorsBatch(
    std::span<const CKKSSlotVector<typename Schedule::Param>> tables,
    std::span<const typename Schedule::Selector> selectors,
    std::type_identity_t<std::span<const CKKSNarrowNTTSelector<Schedule, DigitBits>>> prepared,
    std::size_t begin, const FactorKey &key,
    std::size_t selector_chunk_size = 0,
    const CKKSPreparedAutoKeyFFT<typename Schedule::Param, Schedule::factor_key_log_q> *factor_fft = nullptr,
    CKKSNarrowFactorTimings *timings = nullptr,
    const std::function<void(std::size_t, std::size_t,
                            std::vector<CKKSNarrowNTTSelector<Schedule, DigitBits>> &)> &load_prepared = {},
    const std::function<void(std::size_t, std::size_t,
                            std::vector<typename Schedule::Selector> &)> &load_raw = {},
    bool prefetch_raw = false)
{
    using P = typename Schedule::Param;
    using Cache = CKKSNarrowNTTSelector<Schedule, DigitBits>;
    const auto columns = tables.size();
    if (columns == 0) return {};
    if (load_prepared && load_raw)
        throw std::invalid_argument("conflicting raw and prepared selector loaders");
    if (prefetch_raw && !load_raw)
        throw std::invalid_argument("raw selector prefetch requires a raw loader");
    if (!load_prepared && !load_raw && (begin > selectors.size() || selectors.size() - begin < Schedule::block_width))
        throw std::invalid_argument("invalid narrow selector group");
    if ((load_prepared || load_raw) && (selector_chunk_size == 0 || !prepared.empty() || begin != 0))
        throw std::invalid_argument("invalid prepared selector chunk policy");
    if (prepared.empty() && selector_chunk_size == 0) {
        std::vector<typename Schedule::Factor> result;
        result.reserve(columns);
        for (const auto &table : tables)
            result.push_back(packedFactorReference<Schedule>(table, selectors, begin, key, factor_fft));
        return result;
    }
    if (selector_chunk_size != 0 && !prepared.empty())
        throw std::invalid_argument("chunked selectors must not retain a full NTT cache");
    if (selector_chunk_size == 0 && prepared.size() != selectors.size())
        throw std::invalid_argument("invalid narrow NTT selector cache");
    for (std::size_t y = 0; selector_chunk_size == 0 && y < Schedule::block_width; y++)
        for (const auto &spectrum : prepared[begin + y].spectra)
            if (spectrum.size() != (P::k + 1) * Cache::digits * P::n)
                throw std::invalid_argument("invalid narrow NTT selector spectrum");
    using Clock = std::chrono::steady_clock;
    auto stage_start = timings ? Clock::now() : Clock::time_point{};
    const modular_ntt::TwoPrimeCRT crt(modular_ntt::degree65536_primes[0], modular_ntt::degree65536_primes[1]);
    // These immutable tables depend only on the public ring degree and
    // primes. Reuse them across factor groups and frequency tiles, as the
    // selector transforms already do, instead of rebuilding both plans.
    static const std::array<modular_ntt::NegacyclicNTTPlan, 2> plans = {
        modular_ntt::NegacyclicNTTPlan(P::n, modular_ntt::degree65536_primes[0]),
        modular_ntt::NegacyclicNTTPlan(P::n, modular_ntt::degree65536_primes[1])};
    // NTT points are rows and batch columns are contiguous. Reusing one
    // selector across columns now streams cache lines rather than striding
    // across separate degree-sized arrays. Selector cache format is unchanged.
    std::array<std::array<std::vector<std::uint64_t>, Cache::prime_count>, 2> base;
    for (auto &signed_tables : base)
        for (auto &spectrum : signed_tables) spectrum.resize(columns * P::n);
    for (std::size_t t = 0; t < 2; t++)
      for (std::size_t column = 0; column < columns; column++) {
        CKKSSlotVector<P> transformed;
        for (std::size_t i = 0; i < Schedule::slots; i++) {
            const auto z = tables[column][i];
            transformed[i] = t == 0
                ? (z - std::complex<double>(0, 1) * std::conj(z)) / 2.0
                : (std::conj(z) - std::complex<double>(0, 1) * z) / 2.0;
        }
        auto polynomial_storage = std::make_unique<Polynomial<P>>();
        auto &polynomial = *polynomial_storage;
        ckksSlotEncode<P, Schedule::selector_log_q, Schedule::plain_log_delta>(polynomial, transformed);
        std::vector<std::uint64_t> spectrum(P::n);
        for (std::size_t prime = 0; prime < Cache::prime_count; prime++) {
            const auto modulus = modular_ntt::degree65536_primes[prime].value;
            for (std::size_t i = 0; i < P::n; i++) {
                const auto [negative, magnitude] =
                    ckks_detail::smallSignedMagnitude<P, Schedule::selector_log_q>(polynomial[i]);
                const unsigned __int128 bound = static_cast<unsigned __int128>(magnitude) *
                    Cache::digit_mask * P::n * Schedule::block_width;
                const unsigned __int128 limit = Cache::prime_count == 1
                    ? modular_ntt::degree65536_primes[0].value / 2 : crt.modulusProduct() / 2;
                if (bound >= limit)
                    throw std::invalid_argument("narrow NTT convolution exceeds CRT bound");
                const auto reduced = magnitude % modulus;
                spectrum[i] = negative && reduced != 0 ? modulus - reduced : reduced;
            }
            plans[prime].forward(spectrum);
            for (std::size_t i = 0; i < P::n; i++)
                base[t][prime][i * columns + column] = spectrum[i];
        }
    }
    std::vector<std::size_t> permutation_exponents(Schedule::block_width);
    for (std::size_t y = 0; y < Schedule::block_width; y++)
        permutation_exponents[y] = modular_ntt::power(5, (Schedule::slots - y) % Schedule::slots, 2 * P::n);
    using Sums = std::array<std::vector<std::uint64_t>, Cache::prime_count>;
    const std::array<SelectorMultiplierFactory, 2> multiplier_factories{
        SelectorMultiplierFactory(modular_ntt::degree65536_primes[0].value),
        SelectorMultiplierFactory(modular_ntt::degree65536_primes[1].value)};
    auto add_rows = [&](Sums &first, Sums &second, std::size_t c, std::size_t digit,
                        std::span<const Cache> cache, std::size_t cache_begin,
                        std::size_t first_y, std::size_t count) {
        for (std::size_t prime = 0; prime < Cache::prime_count; prime++) {
            const auto modulus = modular_ntt::degree65536_primes[prime].value;
            const auto &multiplier_factory = multiplier_factories[prime];
            const auto offset = (c * Cache::digits + digit) * P::n;
            addSelectorMatrixRowsPair<P::n>(first[prime].data(), second[prime].data(),
                base[0][prime].data(), base[1][prime].data(), columns, count,
                [&](std::size_t local) {
                    return SelectorMatrixRow{cache[cache_begin + local].spectra[prime].data() + offset,
                        permutation_exponents[first_y + local]};
                }, multiplier_factory, modulus);
        }
    };
    // Keep exact modular sums across selector chunks. Neither inverse NTT,
    // CRT, torus recomposition nor rounding occurs at a chunk boundary.
    if (timings) {
        timings->table_ntt_seconds += std::chrono::duration<double>(Clock::now() - stage_start).count();
        stage_start = Clock::now();
    }
    double selector_seconds = 0, provider_seconds = 0;
    std::vector<Sums> chunk_sums;
    auto row_index = [](std::size_t t, std::size_t c, std::size_t digit) {
        return (t * (P::k + 1) + c) * Cache::digits + digit;
    };
    chunk_sums.resize(2 * (P::k + 1) * Cache::digits);
    for (auto &row : chunk_sums)
        for (auto &sum : row) sum.assign(columns * P::n, 0);
    if (selector_chunk_size != 0) {
        std::future<std::vector<typename Schedule::Selector>> pending_raw;
        auto start_raw_read = [&](std::size_t first) {
            const auto count = std::min(selector_chunk_size, Schedule::block_width - first);
            pending_raw = std::async(std::launch::async, [&, first, count] {
                std::vector<typename Schedule::Selector> raw;
                load_raw(first, count, raw);
                return raw;
            });
        };
        if (prefetch_raw) start_raw_read(0);
        for (std::size_t first_y = 0; first_y < Schedule::block_width;) {
            const auto count = std::min(selector_chunk_size, Schedule::block_width - first_y);
            std::vector<Cache> chunk;
            const auto prepare_start = timings ? Clock::now() : Clock::time_point{};
            if (load_prepared) {
                load_prepared(first_y, count, chunk);
                if (chunk.size() != count) throw std::invalid_argument("prepared selector chunk count mismatch");
                for (const auto &selector : chunk)
                    for (const auto &spectrum : selector.spectra)
                        if (spectrum.size() != (P::k + 1) * Cache::digits * P::n)
                            throw std::invalid_argument("prepared selector chunk length mismatch");
                if (timings) provider_seconds += std::chrono::duration<double>(Clock::now() - prepare_start).count();
            }
            else if (load_raw) {
                std::vector<typename Schedule::Selector> raw;
                if (prefetch_raw) raw = pending_raw.get();
                else load_raw(first_y, count, raw);
                if (raw.size() != count)
                    throw std::invalid_argument("raw selector chunk count mismatch");
                if (timings) provider_seconds += std::chrono::duration<double>(Clock::now() - prepare_start).count();
                const auto ntt_start = timings ? Clock::now() : Clock::time_point{};
                prepareSelectors<Schedule, DigitBits>(chunk, raw);
                if (timings) selector_seconds += std::chrono::duration<double>(Clock::now() - ntt_start).count();
                // Release this raw chunk before reading the next one. One
                // prepared chunk and at most one raw chunk remain resident.
                if (prefetch_raw) {
                    decltype(raw){}.swap(raw);
                    if (first_y + count < Schedule::block_width)
                        start_raw_read(first_y + count);
                }
            }
            else {
                prepareSelectors<Schedule, DigitBits>(chunk, selectors.subspan(begin + first_y, count));
                if (timings) selector_seconds += std::chrono::duration<double>(Clock::now() - prepare_start).count();
            }
#pragma omp parallel for collapse(2) schedule(static)
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t digit = 0; digit < Cache::digits; digit++)
                    add_rows(chunk_sums[row_index(0, c, digit)], chunk_sums[row_index(1, c, digit)],
                        c, digit, chunk, 0, first_y, count);
            first_y += count;
        }
    }
    else {
#pragma omp parallel for collapse(2) schedule(static)
        for (std::size_t c = 0; c <= P::k; c++)
            for (std::size_t digit = 0; digit < Cache::digits; digit++)
                add_rows(chunk_sums[row_index(0, c, digit)], chunk_sums[row_index(1, c, digit)], c, digit,
                         prepared, begin, 0, Schedule::block_width);
    }
    if (timings) {
        timings->selector_ntt_seconds += selector_seconds;
        timings->selector_load_seconds += provider_seconds;
        timings->matrix_seconds += std::chrono::duration<double>(Clock::now() - stage_start).count() - selector_seconds - provider_seconds;
        stage_start = Clock::now();
    }
    // Base spectra are no longer used; release them before inverse CRT and
    // conjugation allocate their own workspaces.
    for (auto &signed_tables : base)
        for (auto &spectrum : signed_tables) std::vector<std::uint64_t>{}.swap(spectrum);
    std::vector<std::array<TRLWE<P>, 2>> accum(columns);
    // Each (table, component, column) owns its torus accumulator. Reconstruct
    // one digit at a time directly into that accumulator instead of retaining
    // a full torus ciphertext for every digit and reading it again afterwards.
    // Modular sums remain separate until the complete selector group is done.
#pragma omp parallel for collapse(3) schedule(static)
    for (std::size_t t = 0; t < 2; t++)
        for (std::size_t c = 0; c <= P::k; c++)
            for (std::size_t column = 0; column < columns; column++) {
                std::array<std::vector<std::uint64_t>, Cache::prime_count> spectra;
                for (auto &spectrum : spectra) spectrum.resize(P::n);
                for (std::size_t digit = 0; digit < Cache::digits; digit++) {
                    for (std::size_t prime = 0; prime < Cache::prime_count; prime++) {
                        auto &spectrum = spectra[prime];
                        const auto &sum = chunk_sums[row_index(t, c, digit)][prime];
                        for (std::size_t i = 0; i < P::n; i++)
                            spectrum[i] = sum[i * columns + column];
                        plans[prime].inverse(spectrum);
                    }
                    auto &destination = accum[column][t][c];
                    for (std::size_t i = 0; i < P::n; i++) {
                        __int128 exact;
                        if constexpr (Cache::prime_count == 1) {
                            const auto value = spectra[0][i];
                            const auto modulus = modular_ntt::degree65536_primes[0].value;
                            exact = value <= modulus / 2 ? static_cast<__int128>(value) :
                                -static_cast<__int128>(modulus - value);
                        }
                        else exact = crt.reconstructSigned(spectra[0][i], spectra[1][i]);
                        destination[i] += ckks_detail::signedToLevel<P, Schedule::selector_log_q>(exact)
                            << (digit * Cache::digit_bits);
                    }
                }
            }
    chunk_sums.clear();
    if (timings) {
        timings->inverse_crt_seconds += std::chrono::duration<double>(Clock::now() - stage_start).count();
        stage_start = Clock::now();
    }
    std::vector<typename Schedule::Factor> result(columns);
    for (std::size_t column = 0; column < columns; column++) {
        auto factor_storage = std::make_unique<std::array<typename Schedule::Factor, 2>>();
        auto &second = (*factor_storage)[0];
        auto &conjugated = (*factor_storage)[1];
        for (std::size_t c = 0; c <= P::k; c++)
            for (std::size_t i = 0; i < P::n; i++) {
                result[column].ct[c][i] = ckks_detail::roundedLevelRightShift<
                    P, Schedule::selector_log_q, Schedule::plain_log_delta>(accum[column][0][c][i]);
                second.ct[c][i] = ckks_detail::roundedLevelRightShift<
                    P, Schedule::selector_log_q, Schedule::plain_log_delta>(accum[column][1][c][i]);
            }
        conjugate<Schedule, true>(conjugated.ct, second.ct, key, factor_fft);
        CKKSAddInPlace<P>(result[column], conjugated);
    }
    if (timings) timings->conjugation_seconds += std::chrono::duration<double>(Clock::now() - stage_start).count();
    return result;
}

template <class Schedule>
inline typename Schedule::Factor packedFactor(
    const CKKSSlotVector<typename Schedule::Param> &table,
    const std::vector<typename Schedule::Selector> &selectors,
    const std::vector<CKKSNarrowNTTSelector<Schedule>> &prepared,
    std::size_t begin, const CKKSNarrowBootstrapKey<Schedule> &key)
{
    if (key.cache_mode == CKKSNarrowCacheMode::Streamed) {
        if (begin > selectors.size() || selectors.size() - begin < Schedule::block_width)
            throw std::invalid_argument("invalid narrow selector group");
        const auto group = std::span<const typename Schedule::Selector>(selectors).subspan(
            begin, Schedule::block_width);
        if (key.selector_ntt_chunk_size != 0)
            return std::move(packedFactorsBatch<Schedule>({&table, 1}, group, {}, 0, key,
                key.selector_ntt_chunk_size).front());
        std::vector<CKKSNarrowNTTSelector<Schedule>> cache;
        prepareSelectors<Schedule>(cache, group);
        return std::move(packedFactorsBatch<Schedule>({&table, 1}, group, cache, 0, key).front());
    }
    return std::move(packedFactorsBatch<Schedule>({&table, 1}, selectors,
        prepared, begin, key).front());
}

template <class Schedule, std::size_t Level = 0,
          class PreparedKeys = typename Schedule::PreparedRelinKeys,
          bool ParallelPairs = false>
inline void productTree(
    typename Schedule::TreeResult &result,
    std::vector<CKKSCiphertext<typename Schedule::Param,
        Schedule::factor_log_q - Level * Schedule::log_delta, Schedule::log_delta>> factors,
    const typename Schedule::RelinKeys &keys,
    const PreparedKeys *prepared = nullptr)
{
    using P = typename Schedule::Param;
    if constexpr (Level == Schedule::tree_depth) {
        if (factors.size() != 1) throw std::invalid_argument("invalid narrow product tree");
        result = std::move(factors.front());
    }
    else {
        using Next = CKKSCiphertext<P, Schedule::factor_log_q -
            (Level + 1) * Schedule::log_delta, Schedule::log_delta>;
        std::vector<Next> next((factors.size() + 1) / 2);
        const auto multiply_pair = [&](std::size_t pair) {
            auto multiply = [&](const auto &relin) {
                if constexpr (Schedule::relinearize_before_rescale)
                    CKKSMultRelinThenRescale<P>(next[pair], factors[2 * pair], factors[2 * pair + 1], relin);
                else CKKSMult<P>(next[pair], factors[2 * pair], factors[2 * pair + 1], relin);
            };
            if constexpr (use_multilimb_digit_fft_v<P>) {
                if (prepared) {
                    multiply(prepared->template get<Level>());
                    return;
                }
            }
            multiply(keys.template get<Level>());
        };
        if constexpr (ParallelPairs) {
            // Independent pairs share immutable prepared keys. A single pair
            // keeps the multiplication's internal team. Propagate failures
            // after joining rather than throwing across an OpenMP boundary.
            std::vector<std::exception_ptr> errors(factors.size() / 2);
#pragma omp parallel for if(factors.size() / 2 > 1) schedule(static)
            for (std::size_t pair = 0; pair < factors.size() / 2; ++pair)
                try { multiply_pair(pair); }
                catch (...) { errors[pair] = std::current_exception(); }
            for (const auto &error : errors)
                if (error) std::rethrow_exception(error);
        }
        else
            for (std::size_t pair = 0; pair < factors.size() / 2; ++pair)
                multiply_pair(pair);
        if (factors.size() % 2 != 0)
            CKKSLevelReduce<P>(next.back(), factors.back());
        productTree<Schedule, Level + 1, PreparedKeys, ParallelPairs>(result, std::move(next), keys, prepared);
    }
}

inline int columnThreads(std::size_t columns);

// Keep each tree's pairing and rounding order, but schedule all products at a
// level together. Prepared keys are shared; tensor and switch work is local.
template<class Schedule,std::size_t Level=0,
    class PreparedKeys=typename Schedule::PreparedRelinKeys,bool ParallelPairs=true>
inline void productTreeBatch(std::vector<typename Schedule::TreeResult>&results,
    std::vector<std::vector<CKKSCiphertext<typename Schedule::Param,
        Schedule::factor_log_q-Level*Schedule::log_delta,Schedule::log_delta>>> factors,
    const typename Schedule::RelinKeys&keys,const PreparedKeys*prepared=nullptr,bool windowed_switches=false)
{
    using P=typename Schedule::Param;
    constexpr auto maximum=std::size_t{1}<<(Schedule::tree_depth-Level);
    for(const auto&tree:factors)
        if(tree.empty()||tree.size()>maximum)throw std::invalid_argument("invalid batched narrow product tree");
    if(factors.empty()){results.clear();return;}
    if constexpr(Level==Schedule::tree_depth){
        results.resize(factors.size());
        for(std::size_t tree=0;tree<factors.size();++tree)results[tree]=std::move(factors[tree].front());
    }else{
        using Next=CKKSCiphertext<P,Schedule::factor_log_q-(Level+1)*Schedule::log_delta,Schedule::log_delta>;
        std::vector<std::vector<Next>> next(factors.size());
        std::vector<std::pair<std::size_t,std::size_t>> pairs;
        for(std::size_t tree=0;tree<factors.size();++tree){
            next[tree].resize((factors[tree].size()+1)/2);
            for(std::size_t pair=0;pair<factors[tree].size()/2;++pair)pairs.emplace_back(tree,pair);
        }
        const auto multiply_pair=[&](std::size_t work){
            const auto [tree,pair]=pairs[work];
            const auto multiply=[&](const auto&relin){
                if constexpr(Schedule::relinearize_before_rescale)
                    CKKSMultRelinThenRescale<P>(next[tree][pair],factors[tree][2*pair],factors[tree][2*pair+1],relin);
                else CKKSMult<P>(next[tree][pair],factors[tree][2*pair],factors[tree][2*pair+1],relin);
            };
            if constexpr(use_multilimb_digit_fft_v<P>)if(prepared){multiply(prepared->template get<Level>());return;}
            multiply(keys.template get<Level>());
        };
        // Preserve the raw-key fallback for parameters which do not select
        // prepared multiplication in the existing tree implementation.
        if(windowed_switches && use_multilimb_digit_fft_v<P>){
            if constexpr(Schedule::relinearize_before_rescale &&
                requires{prepared->template multiplyBatch<Level>(next,factors,pairs);}){
                if(!prepared)throw std::invalid_argument("windowed products require prepared RNS keys");
                prepared->template multiplyBatch<Level>(next,factors,pairs);
            }else throw std::invalid_argument("unsupported windowed product backend");
        }else if constexpr(ParallelPairs){
            std::vector<std::exception_ptr> errors(pairs.size());
            [[maybe_unused]] const auto threads=columnThreads(pairs.size());
#pragma omp parallel for if(pairs.size()>1) schedule(dynamic,1) num_threads(threads)
            for(std::size_t work=0;work<pairs.size();++work){
                try{multiply_pair(work);}catch(...){errors[work]=std::current_exception();}
            }
            for(const auto&error:errors)if(error)std::rethrow_exception(error);
        }else for(std::size_t work=0;work<pairs.size();++work)multiply_pair(work);
        for(std::size_t tree=0;tree<factors.size();++tree)
            if(factors[tree].size()%2)CKKSLevelReduce<P>(next[tree].back(),factors[tree].back());
        factors.clear(); // Completed input levels need not remain live during recursion.
        productTreeBatch<Schedule,Level+1,PreparedKeys,ParallelPairs>(results,std::move(next),keys,prepared,windowed_switches);
    }
}

}  // namespace narrow_detail

template <class Schedule>
inline void CKKSNarrowPrepareKey(CKKSNarrowBootstrapKey<Schedule> &key)
{
    key.head_ntt.clear();
    key.tail_ntt.clear();
    if (key.cache_mode == CKKSNarrowCacheMode::Prepared) {
        narrow_detail::prepareSelectors<Schedule>(key.head_ntt, key.head);
        narrow_detail::prepareSelectors<Schedule>(key.tail_ntt, key.tail);
    }
}

template <class Schedule>
inline void CKKSNarrowPackedKeyGen(
    CKKSNarrowPackedKey<Schedule> &result,
    const Key<typename Schedule::Param> &source,
    const Key<typename Schedule::Param> &evaluation,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    using PackedKey = CKKSNarrowPackedKey<Schedule>;
    narrow_detail::structuredPositions<Schedule>(source);
    for (auto &cache : result.substitutions_fft) cache.reset();
    for (std::size_t channel = 0; channel < 4; channel++) {
        auto bit_storage = std::make_unique<Polynomial<P>>();
        auto &bits = *bit_storage;
        for (std::size_t j = 0; j < P::n; j++) {
            const auto wanted = channel % 2 == 0 ? typename P::T{1} : -typename P::T{1};
            if (source[j] == wanted) bits[j] = typename P::T{1} << Schedule::log_delta;
        }
        ckksEncryptPolynomial<P>(result.indicators[channel], bits, evaluation, noise);
    }
    for (std::size_t r = 0; r < P::nbit; r++)
        CKKSAutoKeyGen<P, PackedKey::log_q>(result.substitutions[r],
                                           1 + (P::n >> r), evaluation, noise);
}

// Decompress without access to either secret key. Track each coefficient index
// through the extraction tree, avoiding bit-reversed leaf order.
// Reconstruction precision can change while the packed ciphertext modulus,
// scale and indicator geometry stay fixed. Reuse the original public keys and
// prepared substitution caches; only the public masks and output level change.
template <class Schedule, class Sink, class PackedSchedule = Schedule>
inline void CKKSNarrowExtractSelectorTerms(
    const CKKSNarrowPackedKey<PackedSchedule> &packed, Sink &&sink,
    std::size_t coefficient_begin = 0,
    std::size_t coefficient_end = Schedule::Param::n)
{
    using P = typename Schedule::Param;
    using PackedKey = CKKSNarrowPackedKey<PackedSchedule>;
    static_assert(std::is_same_v<P, typename PackedSchedule::Param>);
    static_assert(PackedKey::log_q == CKKSNarrowPackedKey<Schedule>::log_q);
    static_assert(PackedSchedule::log_delta == Schedule::log_delta);
    static_assert(PackedSchedule::input_log_q == Schedule::input_log_q);
    static_assert(PackedSchedule::blocks == Schedule::blocks);
    static_assert(PackedSchedule::block_width == Schedule::block_width);
    static_assert(PackedSchedule::head_slots == Schedule::head_slots);
    static_assert(PackedSchedule::tail_width == Schedule::tail_width);
    static_assert(PackedSchedule::tail_packs == Schedule::tail_packs);
    static_assert(PackedSchedule::tail_blocks_per_pack == Schedule::tail_blocks_per_pack);
    using Packed = typename PackedKey::Packed;
    if (coefficient_begin > coefficient_end || coefficient_end > P::n)
        throw std::invalid_argument("invalid packed extraction coefficient range");
    if (coefficient_begin == coefficient_end) return;
    constexpr auto leaf_log_q = PackedKey::log_q - P::nbit;
    using Leaf = CKKSCiphertext<P, leaf_log_q, Schedule::log_delta>;
    using Mask = CKKSPlaintext<P, leaf_log_q, Schedule::plain_log_delta>;
    // Reuse payload buffers throughout the DFS. A parent needs its children
    // while descending, but its scratch is finished before either child runs.
    // Leaf and term buffers are transient; the sink consumes each term during
    // its callback. Only O(log N) ciphertexts remain live.
    std::array<std::unique_ptr<std::array<Packed, 2>>, P::nbit> children_at_depth;
    auto scratch = std::make_unique<std::array<Packed, 2>>();
    auto leaf_storage = std::make_unique<Leaf>();
    auto term_storage = std::make_unique<typename Schedule::Selector>();
    // Every coefficient leaf uses one of these public masks. Encode each
    // distinct mask once rather than repeating the slot transform at each leaf.
    // Lazy construction avoids allocating masks for unused tail positions.
    std::array<std::unique_ptr<Mask>, 2 + 2 * Schedule::tail_blocks_per_pack> masks;
    // Avoid instantiating the FFT cache for non-FFT parameter backends.
    using PreparedMask = std::conditional_t<use_multilimb_digit_fft_v<P>,
        CKKSPreparedSmallPlainFFT<P, leaf_log_q>, Mask>;
    std::array<std::unique_ptr<PreparedMask>, 2 + 2 * Schedule::tail_blocks_per_pack> prepared_masks;
    auto add_masked = [&](bool tail, std::size_t index, const Leaf &leaf,
                          std::size_t begin, std::size_t count, bool imaginary,
                          std::size_t mask_index) {
        if (!masks[mask_index]) {
            auto plain = std::make_unique<Mask>();
            CKKSSlotVector<P> mask{};
            std::fill_n(mask.begin() + begin, count,
                        imaginary ? std::complex<double>{0, 1} : std::complex<double>{1, 0});
            ckksSlotEncode<P, leaf_log_q, Schedule::plain_log_delta>(plain->poly, mask);
            if constexpr (use_multilimb_digit_fft_v<P>)
                prepared_masks[mask_index] = makeCKKSPreparedSmallPlainFFT<P, leaf_log_q>(plain->poly);
            masks[mask_index] = std::move(plain);
        }
        auto &term = *term_storage;
        if constexpr (use_multilimb_digit_fft_v<P>)
            CKKSPreparedPlainMulRescale<P, leaf_log_q, Schedule::log_delta,
                Schedule::plain_log_delta>(term, leaf, *prepared_masks[mask_index]);
        else CKKSPlainMulRescale<P>(term, leaf, *masks[mask_index]);
        sink(tail, index, term);
    };
    // Depth-first extraction keeps only O(log N) packed ciphertexts live.
    // A breadth-first tree plus all four leaf arrays would require O(N^2)
    // torus storage before the final selectors, exceeding LL14 host memory.
    for (std::size_t channel = 0; channel < 4; channel++) {
        auto visit = [&](auto &&self, const Packed &node, std::size_t r,
                         std::size_t index) -> void {
            // Descendants share the low r bits. Their coefficient indices
            // form index + k*2^r; prune branches with no requested leaf.
            const std::size_t stride = std::size_t{1} << r;
            const std::size_t first = index >= coefficient_begin ? index :
                index + ((coefficient_begin - index + stride - 1) / stride) * stride;
            if (first >= coefficient_end) return;
            if (r == P::nbit) {
                auto &leaf = *leaf_storage;
                for (std::size_t c = 0; c <= P::k; c++)
                    for (std::size_t i = 0; i < P::n; i++)
                        leaf.ct[c][i] = ckks_detail::roundedLevelRightShift<P,
                            PackedKey::log_q, P::nbit>(node.ct[c][i]);
                const auto x = index / Schedule::block_width;
                const auto y = index % Schedule::block_width;
                if (channel < 2)
                    add_masked(false, index, leaf, 0, Schedule::head_slots, channel != 0, channel);
                else {
                    const auto q = x / Schedule::tail_blocks_per_pack;
                    const auto ell = x % Schedule::tail_blocks_per_pack;
                    add_masked(true, q * Schedule::block_width + y, leaf,
                        (2 * ell + 1) * Schedule::tail_width, Schedule::tail_width, channel == 3,
                        2 + 2 * ell + (channel == 3));
                }
                return;
            }
            if (!children_at_depth[r]) children_at_depth[r] = std::make_unique<std::array<Packed, 2>>();
            auto &children = children_at_depth[r];
            auto &substituted = (*scratch)[0];
            auto &odd = (*scratch)[1];
            if constexpr (use_multilimb_digit_fft_v<P>) {
                if (packed.substitutions_fft[r])
                    CKKSEvalAuto<P, PackedKey::log_q>(substituted.ct, node.ct,
                        1 + (P::n >> r), *packed.substitutions_fft[r]);
                else CKKSEvalAutoStreamedFFT<P, PackedKey::log_q>(substituted.ct, node.ct,
                    1 + (P::n >> r), packed.substitutions[r]);
            }
            else CKKSEvalAuto<P, PackedKey::log_q>(substituted.ct, node.ct,
                1 + (P::n >> r), packed.substitutions[r]);
            CKKSAdd<P>((*children)[0], node, substituted);
            CKKSSub<P>(odd, node, substituted);
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t i = 0; i < P::n; i++) {
                    const auto shifted = i + (1U << r);
                    (*children)[1].ct[c][i] = shifted < P::n
                        ? odd.ct[c][shifted] : -odd.ct[c][shifted - P::n];
                }
            ckks_detail::reduceTRLWEToLevel<P, PackedKey::log_q>((*children)[1].ct);
            self(self, (*children)[0], r + 1, index);
            self(self, (*children)[1], r + 1, index | (1U << r));
        };
        visit(visit, packed.indicators[channel], 0, 0);
    }

}

// The in-memory adapter accumulates the same ordered public terms as a sink.
template <class Schedule>
inline void CKKSNarrowReconstructSelectors(
    CKKSNarrowBootstrapKey<Schedule> &key,
    const CKKSNarrowPackedKey<Schedule> &packed)
{
    key.selector_source.reset();
    key.factor_source.reset();
    key.head.clear();
    key.head.resize(Schedule::blocks * Schedule::block_width);
    key.tail.clear();
    key.tail.resize(Schedule::tail_packs * Schedule::block_width);
    CKKSNarrowExtractSelectorTerms<Schedule>(packed,
        [&](bool tail, std::size_t index, const typename Schedule::Selector &term) {
            CKKSAddInPlace<typename Schedule::Param>((tail ? key.tail : key.head).at(index), term);
        });
    CKKSNarrowPrepareKey<Schedule>(key);
}

namespace narrow_detail {
template <class Schedule>
inline void generateNarrowEvaluationKeys(CKKSNarrowBootstrapKey<Schedule> &result,
    const Key<typename Schedule::Param> &evaluation, CKKSNoise noise)
{
    using P = typename Schedule::Param;
    result.factor_conjugation_fft.reset();
    result.output_conjugation_fft.reset();
    result.factor_conjugation_projected_fft.reset();
    result.output_conjugation_projected_fft.reset();
    result.tail_rotations_fft.reset();
    result.bounded_tail_rotation_fft = false;
    result.selector_ntt_chunk_size = 0;
    result.staged_fft_caches = false;
    result.product_tree_evaluator = {};
    result.product_tree_batch_evaluator = {};
    result.batch_product_trees = false;
    result.automorphism_evaluator = {};
    result.automorphism_batch_evaluator = {};
    result.relin_fft.reset();
    result.relin_projected_fft.reset();
    CKKSGaloisKeyGen<P, Schedule::factor_key_log_q>(result.factor_galois, evaluation, noise);
    CKKSAutoKeyGen<P, Schedule::output_key_log_q>(result.output_conjugation,
                                             2 * P::n - 1, evaluation, noise);
    CKKSDDRelinKeyChainGen<P, Schedule::relin_chain_start_log_q, Schedule::log_delta,
                         Schedule::tree_depth>(result.relin, evaluation, noise);
}
}  // namespace narrow_detail

// Generate reusable evaluation setup without encrypting expanded selectors.
// Supply the four packed indicators separately and call ReconstructSelectors
// before evaluation. This also discards selectors/caches from a previous key.
template <class Schedule>
inline void CKKSNarrowSetupKeyGen(CKKSNarrowBootstrapKey<Schedule> &result,
    const Key<typename Schedule::Param> &evaluation,
    CKKSNoise noise = {Schedule::Param::α, 0},
    CKKSNarrowCacheMode cache_mode = CKKSNarrowCacheMode::Prepared)
{
    result.cache_mode = cache_mode;
    result.selector_source.reset();
    result.factor_source.reset();
    result.head.clear();
    result.tail.clear();
    result.head_ntt.clear();
    result.tail_ntt.clear();
    narrow_detail::generateNarrowEvaluationKeys<Schedule>(result, evaluation, noise);
}

// Encrypt one expanded selector group under existing private keys. This lets
// callers retain evaluation keys and build bounded file-backed setup without
// regenerating automorphism or DD keys. Every selector gets fresh encryption.
template <class Schedule>
inline void CKKSNarrowEncryptSelectorGroup(
    std::vector<typename Schedule::Selector> &result, bool tail, std::size_t group,
    const Key<typename Schedule::Param> &source,
    const Key<typename Schedule::Param> &evaluation,
    CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    if (group >= (tail ? Schedule::tail_packs : Schedule::blocks))
        throw std::out_of_range("narrow selector encryption group");
    const auto positions = narrow_detail::structuredPositions<Schedule>(source);
    result.resize(Schedule::block_width);
    auto indicator = [&](std::size_t block, std::size_t y) -> std::complex<double> {
        if (y >= Schedule::block_size(block) || positions[block] != static_cast<int>(y)) return {};
        return source[block * Schedule::block_width + y] == typename P::T{1}
            ? std::complex<double>{1, 0} : std::complex<double>{0, 1};
    };
#pragma omp parallel for schedule(static)
    for (std::size_t y = 0; y < Schedule::block_width; ++y) {
        CKKSSlotVector<P> mask{};
        if (!tail) std::fill_n(mask.begin(), Schedule::head_slots, indicator(group, y));
        else for (std::size_t ell = 0; ell < Schedule::tail_blocks_per_pack; ++ell) {
            const auto block = group * Schedule::tail_blocks_per_pack + ell;
            if (block >= Schedule::blocks) break;
            std::fill_n(mask.begin() + (2 * ell + 1) * Schedule::tail_width,
                        Schedule::tail_width, indicator(block, y));
        }
        ckksSlotEncrypt<P>(result[y], mask, evaluation, noise);
    }
}

template <class Schedule>
inline void CKKSNarrowBootstrapKeyGen(
    CKKSNarrowBootstrapKey<Schedule> &result,
    const Key<typename Schedule::Param> &source,
    const Key<typename Schedule::Param> &evaluation,
    CKKSNoise noise = {Schedule::Param::α, 0},
    CKKSNarrowCacheMode cache_mode = CKKSNarrowCacheMode::Prepared)
{
    using P = typename Schedule::Param;
    const auto positions = narrow_detail::structuredPositions<Schedule>(source);
    result.cache_mode = cache_mode;
    result.selector_source.reset();
    result.factor_source.reset();
    result.head.resize(Schedule::blocks * Schedule::block_width);
    result.tail.resize(Schedule::tail_packs * Schedule::block_width);
    auto indicator = [&](std::size_t x, std::size_t y) -> std::complex<double> {
        if (y >= Schedule::block_size(x)) return {};
        if (positions[x] != static_cast<int>(y)) return {};
        return source[x * Schedule::block_width + y] == typename P::T{1}
            ? std::complex<double>{1, 0} : std::complex<double>{0, 1};
    };
#pragma omp parallel for collapse(2) schedule(static)
    for (std::size_t x = 0; x < Schedule::blocks; x++)
        for (std::size_t y = 0; y < Schedule::block_width; y++) {
            CKKSSlotVector<P> mask{};
            std::fill_n(mask.begin(), Schedule::head_slots, indicator(x, y));
            ckksSlotEncrypt<P>(result.head[x * Schedule::block_width + y], mask, evaluation, noise);
        }
#pragma omp parallel for collapse(2) schedule(static)
    for (std::size_t q = 0; q < Schedule::tail_packs; q++)
        for (std::size_t y = 0; y < Schedule::block_width; y++) {
            CKKSSlotVector<P> mask{};
            for (std::size_t ell = 0; ell < Schedule::tail_blocks_per_pack; ell++) {
                const std::size_t x = q * Schedule::tail_blocks_per_pack + ell;
                if (x == Schedule::blocks) break;
                std::fill_n(mask.begin() + (2 * ell + 1) * Schedule::tail_width,
                            Schedule::tail_width, indicator(x, y));
            }
            ckksSlotEncrypt<P>(result.tail[q * Schedule::block_width + y], mask, evaluation, noise);
        }
    narrow_detail::generateNarrowEvaluationKeys<Schedule>(result, evaluation, noise);
    CKKSNarrowPrepareKey<Schedule>(result);
}

template <class Schedule, std::uint32_t InputLogDelta>
struct CKKSNarrowBlindRotationRequest {
    using Input = CKKSCiphertext<typename Schedule::Param,
        Schedule::input_log_q, InputLogDelta>;
    const Input *input = nullptr;
    std::uint32_t frequency = 1;
    std::complex<double> weight = {1, 0};
};

// Optional observation of completed positive-frequency roots. The observer
// receives public requests and ciphertexts; it never participates in evaluation.
template <class Schedule, std::uint32_t InputLogDelta>
using CKKSNarrowRootObserver = std::function<void(
    const CKKSNarrowBlindRotationRequest<Schedule, InputLogDelta> &,
    const typename Schedule::Output &)>;

namespace narrow_detail {
// A column team needs no more workers than independent columns. A clause on
// the region preserves the caller's thread budget for the subsequent stages.
inline int columnThreads(std::size_t columns)
{
#ifdef _OPENMP
    return static_cast<int>(std::max(std::size_t{1},std::min(columns,std::size_t(omp_get_max_threads()))));
#else
    return 1;
#endif
}

// The caller retains ordering of observers after the independent work joins.
template<class Evaluate>
inline void evaluateRootColumns(std::size_t columns,bool parallel,Evaluate &&evaluate)
{
    if(parallel && columns>1){
        std::exception_ptr error;
        [[maybe_unused]] const auto threads=columnThreads(columns);
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for(std::size_t column=0;column<columns;++column){
            try{evaluate(column);}
            catch(...){
#pragma omp critical(ckks_root_column_error)
                {if(!error)error=std::current_exception();}
            }
        }
        if(error)std::rethrow_exception(error);
    }
    else for(std::size_t column=0;column<columns;++column)evaluate(column);
}

// Each tail leaf is independent and may share the immutable public key cache.
// No key preparation occurs here; OpenMP nesting follows the caller policy.
template<class Schedule>
inline void alignPublicTail(std::vector<std::vector<typename Schedule::Factor>> &tail,
    std::span<const int> steps,const CKKSNarrowBootstrapKey<Schedule> &key)
{
    using P=typename Schedule::Param;
    if(!key.automorphism_evaluator || steps.size()!=tail.size()*Schedule::blocks)
        throw std::invalid_argument("invalid public tail alignment batch");
    for(const auto &column:tail)
        if(column.size()!=Schedule::blocks)throw std::invalid_argument("invalid public tail alignment column");
    const auto align=[&](std::size_t i){
        auto &cipher=tail[i/Schedule::blocks][i%Schedule::blocks].ct;
        key.automorphism_evaluator(cipher,cipher,Schedule::factor_log_q,
            ckks_detail::rotationAutomorphismForSteps<P>(steps[i]));
    };
    if(key.parallel_public_tail_alignment && steps.size()>1){
        std::exception_ptr error;
#pragma omp parallel for schedule(dynamic,1)
        for(std::size_t work=0;work<steps.size();++work){
            // Adjacent workers rotate different columns of the same block,
            // sharing its prepared exponent cache instead of cycling caches.
            const auto i=(work%tail.size())*Schedule::blocks+work/tail.size();
            try{align(i);}
            catch(...){
#pragma omp critical(ckks_public_tail_alignment_error)
                {if(!error)error=std::current_exception();}
            }
        }
        if(error)std::rethrow_exception(error);
    }
    else for(std::size_t i=0;i<steps.size();++i)align(i);
}
} // namespace narrow_detail

// Sections 5.2 and 7: input and frequency columns share the prepared W_x.
// Tiling bounds temporary memory without changing the mathematical batch.
template <class Schedule, std::uint32_t InputLogDelta>
inline void CKKSNarrowBlindRotateBatch(
    std::vector<typename Schedule::Output> &results,
    std::span<const CKKSNarrowBlindRotationRequest<Schedule, InputLogDelta>> requests,
    const CKKSNarrowBootstrapKey<Schedule> &key,
    std::size_t tile_size = 8,
    const CKKSNarrowRootObserver<Schedule, InputLogDelta> *observer = nullptr)
{
    using P = typename Schedule::Param;
    if (tile_size == 0) throw std::invalid_argument("narrow batch tile size must be positive");
    if(key.batch_product_trees&&!key.product_tree_batch_evaluator)
        throw std::invalid_argument("batched product trees require a prepared batch backend");
    for (const auto &request : requests)
        if (!request.input || request.frequency == 0 || !std::isfinite(request.weight.real()) ||
            !std::isfinite(request.weight.imag()))
            throw std::invalid_argument("invalid narrow bootstrap batch request");
    if ((key.factor_source ? key.factor_source->size(false) :
         key.selector_source ? key.selector_source->size(false) : key.head.size()) !=
            Schedule::blocks * Schedule::block_width ||
        (key.factor_source ? key.factor_source->size(true) :
         key.selector_source ? key.selector_source->size(true) : key.tail.size()) !=
            Schedule::tail_packs * Schedule::block_width)
        throw std::invalid_argument("invalid narrow bootstrap selector key");
    results.resize(requests.size());
    for (std::size_t begin = 0; begin < requests.size(); begin += std::min(tile_size, requests.size() - begin)) {
        const auto columns = std::min(tile_size, requests.size() - begin);
        using Clock = std::chrono::steady_clock;
        CKKSNarrowRootTimings root_times;
        auto root_stage_start = key.root_stage_observer ? Clock::now() : Clock::time_point{};
        const auto finish_root_stage = [&](double &seconds) {
            if (key.root_stage_observer) {
                const auto now = Clock::now();
                seconds = std::chrono::duration<double>(now - root_stage_start).count();
                root_stage_start = now;
            }
        };
        std::vector<std::array<double, P::n>> a(columns), b(columns);
        for (std::size_t column = 0; column < columns; column++) {
            const auto &input = *requests[begin + column].input;
            for (std::size_t i = 0; i < P::n; i++) {
                // TFHEpp decrypts b-a*s; the paper decrypts b+a*s.
                a[column][i] = -static_cast<double>(std::ldexp(ckks_detail::levelToLongDouble<
                    P, Schedule::input_log_q>(input.ct[0][i]), -static_cast<int>(Schedule::input_log_q)));
                b[column][i] = static_cast<double>(std::ldexp(ckks_detail::levelToLongDouble<
                    P, Schedule::input_log_q>(input.ct[1][i]), -static_cast<int>(Schedule::input_log_q)));
            }
        }
        std::vector<std::vector<typename Schedule::Factor>> head(columns), tail(columns);
        for (std::size_t column = 0; column < columns; column++) {
            head[column].resize(Schedule::blocks);
            head[column].reserve(Schedule::blocks + 1);
            tail[column].reserve(Schedule::blocks + 1);
        }
        std::vector<CKKSSlotVector<P>> tables(columns);
        // Evaluation owns these caches, leaving the public key immutable.
        // Factor and product-tree caches have disjoint lifetimes per tile.
        std::unique_ptr<CKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q>> staged_factor;
        if constexpr (use_multilimb_digit_fft_v<P>)
            if (key.staged_fft_caches && !key.factor_conjugation_projected_fft && !key.automorphism_evaluator)
                staged_factor = makeCKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q>(key.factor_galois[P::nbit]);
        auto evaluate_group = [&](const auto &selectors, const auto &prepared, std::size_t offset, bool is_tail) {
            CKKSNarrowFactorTimings timings;
            auto *profile = key.factor_observer ? &timings : nullptr;
            auto calculate = [&]() {
                if (key.factor_source)
                    return key.factor_source->evaluate(tables, is_tail,
                        offset / Schedule::block_width, key, staged_factor.get(), profile);
                if (key.selector_source) {
                    if (key.selector_source->hasPreparedChunks() && key.selector_ntt_chunk_size != 0 &&
                        key.cache_mode != CKKSNarrowCacheMode::Reference) {
                        const std::function<void(std::size_t, std::size_t,
                            std::vector<CKKSNarrowNTTSelector<Schedule>> &)> loader =
                            [&](std::size_t first, std::size_t count, auto &chunk) {
                                key.selector_source->loadPreparedChunk(is_tail,
                                    offset / Schedule::block_width, first, count, chunk);
                            };
                        return narrow_detail::packedFactorsBatch<Schedule>(tables, {}, {}, 0, key,
                            key.selector_ntt_chunk_size, staged_factor.get(), profile, loader);
                    }
                    if (key.selector_source->hasRawChunks() && key.selector_ntt_chunk_size != 0 &&
                        key.cache_mode != CKKSNarrowCacheMode::Reference) {
                        const std::function<void(std::size_t, std::size_t,
                            std::vector<typename Schedule::Selector> &)> loader =
                            [&](std::size_t first, std::size_t count, auto &chunk) {
                                key.selector_source->loadChunk(is_tail,
                                    offset / Schedule::block_width, first, count, chunk);
                            };
                        return narrow_detail::packedFactorsBatch<Schedule>(tables, {}, {}, 0, key,
                            key.selector_ntt_chunk_size, staged_factor.get(), profile, {}, loader,
                            key.selector_source->supportsRawChunkPrefetch());
                    }
                    std::vector<typename Schedule::Selector> loaded;
                    const auto load_start = profile ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
                    key.selector_source->load(is_tail, offset / Schedule::block_width, loaded);
                    if (profile) profile->selector_load_seconds += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - load_start).count();
                    if (loaded.size() != Schedule::block_width)
                        throw std::invalid_argument("invalid external narrow selector group");
                    if (key.selector_ntt_chunk_size != 0 && key.cache_mode != CKKSNarrowCacheMode::Reference)
                        return narrow_detail::packedFactorsBatch<Schedule>(tables, loaded, {}, 0, key, key.selector_ntt_chunk_size, staged_factor.get(), profile);
                    std::vector<CKKSNarrowNTTSelector<Schedule>> cache;
                    if (key.cache_mode != CKKSNarrowCacheMode::Reference)
                        narrow_detail::prepareSelectors<Schedule>(cache, loaded);
                    return narrow_detail::packedFactorsBatch<Schedule>(tables, loaded, cache, 0, key, 0, staged_factor.get(), profile);
                }
                if (key.cache_mode == CKKSNarrowCacheMode::Streamed) {
                    const auto group = std::span<const typename Schedule::Selector>(selectors).subspan(
                        offset, Schedule::block_width);
                    if (key.selector_ntt_chunk_size != 0)
                        return narrow_detail::packedFactorsBatch<Schedule>(tables, group, {}, 0, key, key.selector_ntt_chunk_size, staged_factor.get(), profile);
                    std::vector<CKKSNarrowNTTSelector<Schedule>> cache;
                    narrow_detail::prepareSelectors<Schedule>(cache, group);
                    return narrow_detail::packedFactorsBatch<Schedule>(tables, group, cache, 0, key, 0, staged_factor.get(), profile);
                }
                return narrow_detail::packedFactorsBatch<Schedule>(tables, selectors, prepared, offset, key, 0, staged_factor.get(), profile);
            };
            auto result = calculate();
            if (result.size() != columns)
                throw std::invalid_argument("invalid external narrow factor batch");
            if (key.factor_observer) key.factor_observer(is_tail, offset / Schedule::block_width, timings);
            return result;
        };
        finish_root_stage(root_times.input_and_factor_cache_seconds);
        for (std::size_t x = 0; x < Schedule::blocks; x++) {
            for (std::size_t column = 0; column < columns; column++)
                tables[column] = narrow_detail::headTable<Schedule>(a[column], x, requests[begin + column].frequency);
            auto factors = evaluate_group(key.head, key.head_ntt, x * Schedule::block_width, false);
            for (std::size_t column = 0; column < columns; column++) head[column][x] = std::move(factors[column]);
        }
        for (std::size_t q = 0; q < Schedule::tail_packs; q++) {
            for (std::size_t column = 0; column < columns; column++)
                tables[column] = narrow_detail::tailTable<Schedule>(a[column], q, requests[begin + column].frequency);
            const auto packed = evaluate_group(key.tail, key.tail_ntt, q * Schedule::block_width, true);
            const auto count = std::min(Schedule::tail_blocks_per_pack,
                Schedule::blocks - q * Schedule::tail_blocks_per_pack);
            if constexpr (use_multilimb_digit_fft_v<P>) {
                if (key.bounded_tail_rotation_fft) {
                    for (std::size_t column = 0; column < columns; column++)
                        for (std::size_t ell = 0; ell < count; ell++)
                            tail[column].push_back(packed[column]);
                    continue;
                }
            }
            for (std::size_t column = 0; column < columns; column++)
                for (std::size_t ell = 0; ell < count; ell++) {
                    tail[column].emplace_back();
                    auto &aligned = tail[column].back();
                    const int right = Schedule::slots - 2 * (ell + 1) * Schedule::tail_width;
                    if (key.automorphism_evaluator) {
                        key.automorphism_evaluator(aligned.ct, packed[column].ct,
                            Schedule::factor_log_q, ckks_detail::rotationAutomorphismForSteps<P>(-right));
                        continue;
                    }
                    if constexpr (use_multilimb_digit_fft_v<P>) {
                        if (key.tail_rotations_fft) {
                            narrow_detail::rotate<Schedule>(aligned.ct, packed[column].ct,
                                -right, *key.tail_rotations_fft);
                            continue;
                        }
                    }
                    narrow_detail::rotate<Schedule>(aligned.ct, packed[column].ct,
                        -right, key.factor_galois);
                }
        }
        finish_root_stage(root_times.leaf_factors_and_inline_alignment_seconds);
        staged_factor.reset();
        if constexpr (use_multilimb_digit_fft_v<P>) {
            if (key.bounded_tail_rotation_fft) {
                // All tail leaves already reside in the product-tree vectors.
                // Prepare each binary rotation key once for the whole tile,
                // retaining the same ascending-bit order for each ciphertext.
                std::vector<int> steps(columns * Schedule::blocks);
                for (std::size_t column = 0; column < columns; column++)
                    for (std::size_t x = 0; x < Schedule::blocks; x++)
                        steps[column * Schedule::blocks + x] = -static_cast<int>(Schedule::slots -
                            2 * (x % Schedule::tail_blocks_per_pack + 1) * Schedule::tail_width);
                if (key.automorphism_evaluator) {
                    narrow_detail::alignPublicTail<Schedule>(tail,steps,key);
                }
                else CKKSRotateSlotsBatchWithBoundedFFT<P, Schedule::factor_log_q, Schedule::auxiliary_log_scale>(
                    steps, [&](std::size_t i) -> auto & {
                        return tail[i / Schedule::blocks][i % Schedule::blocks].ct;
                    }, key.factor_galois);
            }
        }
        finish_root_stage(root_times.bounded_tail_alignment_seconds);
        std::unique_ptr<typename Schedule::PreparedRelinKeys> staged_relin;
        if constexpr (use_multilimb_digit_fft_v<P>)
            if (key.staged_fft_caches && !key.relin_projected_fft && !key.product_tree_evaluator)
                staged_relin = CKKSNarrowMakePreparedRelinearization<Schedule>(key.relin);
        const auto *tree_keys = staged_relin ? staged_relin.get() : key.relin_fft.get();
        finish_root_stage(root_times.dd_cache_preparation_seconds);
        const bool parallel_roots=key.parallel_root_columns && columns>1;
        const bool batch_trees=key.batch_product_trees;
        std::vector<std::vector<typename Schedule::Factor>> forests(batch_trees?2*columns:0);
        narrow_detail::evaluateRootColumns(columns,parallel_roots,[&](std::size_t column) {
            const auto &request = requests[begin + column];
            CKKSSlotVector<P> head_beta{}, tail_beta{};
            for (std::size_t i = 0; i < Schedule::slots; i++) {
                const auto value = request.weight * narrow_detail::root<P>(b[column], i, 1, request.frequency);
                (i < Schedule::head_slots ? head_beta : tail_beta)[i] = value;
            }
            if (key.factor_source) key.factor_source->adjustBeta(head_beta, tail_beta);
            head[column].emplace_back();
            tail[column].emplace_back();
            auto &transparent_head = head[column].back();
            auto &transparent_tail = tail[column].back();
            ckksSlotEncode<P, Schedule::factor_log_q, Schedule::log_delta>(transparent_head.ct[P::k], head_beta);
            ckksSlotEncode<P, Schedule::factor_log_q, Schedule::log_delta>(transparent_tail.ct[P::k], tail_beta);
            if(batch_trees){
                forests[2*column]=std::move(head[column]);forests[2*column+1]=std::move(tail[column]);return;
            }
            auto tail_result = std::make_unique<typename Schedule::TreeResult>();
            auto evaluate_trees = [&](const auto *prepared) {
                narrow_detail::productTree<Schedule>(results[begin + column], std::move(head[column]), key.relin, prepared);
                narrow_detail::productTree<Schedule>(*tail_result, std::move(tail[column]), key.relin, prepared);
            };
            if (key.product_tree_evaluator) {
                key.product_tree_evaluator(results[begin + column], std::move(head[column]), key.relin);
                key.product_tree_evaluator(*tail_result, std::move(tail[column]), key.relin);
            }
            else if (key.relin_projected_fft) evaluate_trees(key.relin_projected_fft.get());
            else evaluate_trees(tree_keys);
            CKKSAddInPlace<P>(results[begin + column], *tail_result);
            if (observer && !parallel_roots) (*observer)(request, results[begin + column]);
        });
        if(batch_trees){
            std::vector<typename Schedule::TreeResult> products;
            key.product_tree_batch_evaluator(products,std::move(forests),key.relin);
            if(products.size()!=2*columns)throw std::runtime_error("batched product backend returned the wrong root count");
            for(std::size_t column=0;column<columns;++column){
                results[begin+column]=std::move(products[2*column]);
                CKKSAddInPlace<P>(results[begin+column],products[2*column+1]);
            }
        }
        if(observer && (parallel_roots||batch_trees))
            for(std::size_t column=0;column<columns;++column)
                (*observer)(requests[begin+column],results[begin+column]);
        finish_root_stage(root_times.dd_trees_and_beta_seconds);
        if (key.root_stage_observer) key.root_stage_observer(root_times);
    }
}

// Encrypt exp(2*pi*i*k*(b-a*s)_i/q), i<N/2, with optional public weight.
template <class Schedule, std::uint32_t InputLogDelta>
inline void CKKSNarrowBlindRotate(
    typename Schedule::Output &result,
    const CKKSCiphertext<typename Schedule::Param, Schedule::input_log_q, InputLogDelta> &input,
    const CKKSNarrowBootstrapKey<Schedule> &key,
    std::uint32_t frequency = 1, std::complex<double> weight = {1, 0})
{
    const CKKSNarrowBlindRotationRequest<Schedule, InputLogDelta> request{&input, frequency, weight};
    std::vector<typename Schedule::Output> outputs;
    CKKSNarrowBlindRotateBatch<Schedule, InputLogDelta>(outputs, {&request, 1}, key);
    result = std::move(outputs.front());
}

// Section 7: c_0 + sum_{k>0}(c_k*w^(k*m) + conj(c_k*w^(k*m))).
// Each positive frequency is generated directly by blind rotation with the
// same selectors. No encrypted exponential approximation or powering occurs.
// The caller supplies an offline trigonometric approximant; coefficients[0]
// must be real so that the result represents a real periodic function.
template <class Schedule, std::uint32_t InputLogDelta>
inline void CKKSNarrowEvalSeriesBatch(
    std::vector<typename Schedule::Output> &results,
    std::span<const CKKSCiphertext<typename Schedule::Param, Schedule::input_log_q, InputLogDelta>> inputs,
    const CKKSNarrowBootstrapKey<Schedule> &key,
    const std::vector<std::complex<double>> &coefficients, std::size_t tile_size = 8,
    const CKKSNarrowRootObserver<Schedule, InputLogDelta> *observer = nullptr)
{
    using P = typename Schedule::Param;
    if (coefficients.empty() || coefficients[0].imag() != 0 ||
        coefficients.size() > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("invalid real trigonometric series");
    for (const auto value : coefficients)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::invalid_argument("nonfinite trigonometric coefficient");
    if (tile_size == 0) throw std::invalid_argument("narrow series tile size must be positive");
    std::vector<CKKSNarrowBlindRotationRequest<Schedule, InputLogDelta>> requests;
    std::vector<std::size_t> owners;
    results.resize(inputs.size());
    for (std::size_t input_index = 0; input_index < inputs.size(); input_index++) {
        CKKSSetTransparentReal<P>(results[input_index], 0);
        for (std::size_t k = 1; k < coefficients.size(); k++) {
            if (coefficients[k] == std::complex<double>{}) continue;
            requests.push_back({&inputs[input_index], static_cast<std::uint32_t>(k), coefficients[k]});
            owners.push_back(input_index);
        }
    }
    // Consume tiles immediately to bound encrypted roots' temporary storage.
    for (std::size_t begin = 0; begin < requests.size(); begin += std::min(tile_size, requests.size() - begin)) {
        const auto count = std::min(tile_size, requests.size() - begin);
        std::vector<typename Schedule::Output> positive;
        CKKSNarrowBlindRotateBatch<Schedule, InputLogDelta>(positive,
            std::span<const CKKSNarrowBlindRotationRequest<Schedule, InputLogDelta>>(requests).subspan(begin, count),
            key, tile_size, observer);
        for (std::size_t i = 0; i < count; i++)
            CKKSAddInPlace<P>(results[owners[begin + i]], positive[i]);
    }
    // Conjugation is linear. Accumulate all positive frequencies before
    // switching the conjugate back to the evaluation key, once per input.
    // Keep the large-ring ciphertext temporary off the evaluator's stack.
    auto negative = std::make_unique<typename Schedule::Output>();
    for (auto &result : results) {
        if (!requests.empty()) {
            narrow_detail::conjugate<Schedule, false>(negative->ct, result.ct, key);
            CKKSAddInPlace<P>(result, *negative);
        }
        CKKSAddPlainRealInPlace<P>(result, coefficients[0].real());
    }
}

template <class Schedule, std::uint32_t InputLogDelta>
inline void CKKSNarrowEvalSeries(
    typename Schedule::Output &result,
    const CKKSCiphertext<typename Schedule::Param, Schedule::input_log_q, InputLogDelta> &input,
    const CKKSNarrowBootstrapKey<Schedule> &key,
    const std::vector<std::complex<double>> &coefficients, std::size_t tile_size = 8,
    const CKKSNarrowRootObserver<Schedule, InputLogDelta> *observer = nullptr)
{
    std::vector<typename Schedule::Output> outputs;
    CKKSNarrowEvalSeriesBatch<Schedule, InputLogDelta>(outputs, {&input, 1}, key, coefficients, tile_size, observer);
    result = std::move(outputs.front());
}

}  // namespace TFHEpp
#endif
