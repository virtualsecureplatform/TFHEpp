#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "gl/gl.hpp"

#ifdef TFHEPP_DEFAULT_128BIT_PARAMS

namespace TFHEpp {

// The modulus figures reported in Table 1 of ePrint 2026/811.  These are
// storage/security profiles, not claims that the coefficient-domain reference
// kernels below reproduce the paper's RNS performance.
template <class GLP>
struct GLSHIPPaperParameterProfile;

template <>
struct GLSHIPPaperParameterProfile<GL256p17Parameter> {
    static constexpr std::uint32_t log_q = 180;
    static constexpr std::uint32_t log_p = 34;
    static constexpr std::uint32_t log_pq = 214;
    static constexpr std::uint32_t stc_bits = 26;
    static constexpr std::uint32_t gap_log = 5;
    static constexpr std::uint32_t outside_multiplicative_depth = 1;
    static constexpr std::uint32_t security_limit_log_pq = 214;
    static constexpr std::uint32_t dnum = 8;
    static constexpr std::uint32_t masked_column_count = 848;
};

template <>
struct GLSHIPPaperParameterProfile<GL512p17Parameter> {
    static constexpr std::uint32_t log_q = 338;
    static constexpr std::uint32_t log_p = 92;
    static constexpr std::uint32_t log_pq = 430;
    static constexpr std::uint32_t stc_bits = 37;
    static constexpr std::uint32_t gap_log = 11;
    static constexpr std::uint32_t outside_multiplicative_depth = 1;
    static constexpr std::uint32_t security_limit_log_pq = 430;
    static constexpr std::uint32_t dnum = 4;
    static constexpr std::uint32_t masked_column_count = 1504;
};

template <>
struct GLSHIPPaperParameterProfile<GL1024p17Parameter> {
    static constexpr std::uint32_t log_q = 641;
    static constexpr std::uint32_t log_p = 220;
    static constexpr std::uint32_t log_pq = 861;
    static constexpr std::uint32_t stc_bits = 39;
    static constexpr std::uint32_t gap_log = 11;
    static constexpr std::uint32_t outside_multiplicative_depth = 8;
    static constexpr std::uint32_t security_limit_log_pq = 868;
    static constexpr std::uint32_t dnum = 4;
    static constexpr std::uint32_t masked_column_count = 2880;
};

template <class GLP>
inline constexpr bool GLSHIPPaperProfileFitsStorage =
    GLSHIPPaperParameterProfile<GLP>::log_pq <=
    std::numeric_limits<typename GLP::T>::digits;

static_assert(GLSHIPPaperProfileFitsStorage<GL256p17Parameter>);
static_assert(GLSHIPPaperProfileFitsStorage<GL512p17Parameter>);
static_assert(GLSHIPPaperProfileFitsStorage<GL1024p17Parameter>);

namespace gl_ship_detail {

constexpr std::uint32_t ceilLog(const std::uint32_t value,
                                const std::uint32_t radix)
{
    if (value <= 1) return 0;
    std::uint32_t power = 1;
    std::uint32_t result = 0;
    while (power < value) {
        power *= radix;
        result++;
    }
    return result;
}

constexpr bool isPowerOfTwo(const std::uint32_t value)
{
    return value != 0 && (value & (value - 1)) == 0;
}

constexpr std::uint32_t positiveMod(const std::int64_t value,
                                    const std::uint32_t modulus)
{
    const std::int64_t reduced = value % static_cast<std::int64_t>(modulus);
    return static_cast<std::uint32_t>(reduced < 0 ? reduced + modulus
                                                  : reduced);
}

}  // namespace gl_ship_detail

// Forward declarations keep the public Algorithm 3 entry points near the
// parameter table while the implementation details remain below.
template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
struct GLBasePlaintext;
template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
struct GLBaseCiphertext;
template <class GLP>
struct GLBaseHadamardWorkspace;
template <class GLP>
class GLBaseSlotTable;
template <class Schedule>
struct GLSHIPBootstrapKey;
template <class Schedule>
struct GLSHIPSlotsToCoefficientsKey;
template <class Schedule>
struct GLSHIPMaskedColumnKey;
template <class Schedule>
struct GLSHIPHMuxKey;
template <class Schedule>
struct GLSHIPProductRelinKeyChain;

struct GLSHIPSupportInterval {
    std::uint32_t start = 0;
    std::uint32_t width = 0;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(start, width);
    }
};

struct GLSHIPCandidate {
    std::uint32_t fine_x = 0;
    std::uint32_t w = 0;
    std::uint32_t gaussian_phase = 0;

    friend bool operator<(const GLSHIPCandidate &lhs,
                          const GLSHIPCandidate &rhs)
    {
        return std::tie(lhs.fine_x, lhs.w, lhs.gaussian_phase) <
               std::tie(rhs.fine_x, rhs.w, rhs.gaussian_phase);
    }
    friend bool operator==(const GLSHIPCandidate &lhs,
                           const GLSHIPCandidate &rhs) = default;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(fine_x, w, gaussian_phase);
    }
};

// Optional execution tuning for the memory-bound sparse-factor phase.  A
// zero HMux worker count preserves the fused MaskedColumn/HMux loop.  Setting
// it below the active OpenMP team size evaluates bounded tiles in two stages,
// allowing MaskedColumn to use the full team while HMux uses fewer workers.
// Batched HMux products amortize exact modular reduction on AVX-512DQ.  Both
// choices are host-dependent, so neither the physical-core count nor an SMT
// policy is inferred by the library.
struct GLSHIPBootstrapExecutionOptions {
    std::uint32_t hmux_threads = 0;
    std::size_t factor_tile_size = 256;
    bool batch_hmux_products = false;
};

template <class GLP, std::uint32_t LogQ, std::uint32_t InputLogDelta,
          std::uint32_t PlainLogDelta>
using GLRawProductCiphertext =
    GLCiphertext<GLP, LogQ, InputLogDelta + PlainLogDelta>;

template <class Schedule>
inline void GLSHIPSlotsToCoefficients(
    typename Schedule::CoefficientCiphertext &result,
    const typename Schedule::InputCiphertext &input,
    const GLSHIPSlotsToCoefficientsKey<Schedule> &stc_key);

template <class Schedule>
inline void GLSHIPSlotsToCoefficientsKeyGen(
    GLSHIPSlotsToCoefficientsKey<Schedule> &stc_key,
    const Key<typename Schedule::Parameter::baseP> &dense_key, CKKSNoise noise);

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLBaseEncode(GLBasePlaintext<GLP, LogQ, LogDelta> &plaintext,
                         const GLBaseSlotTable<GLP> &slots);

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDSmallKeySwitchBaseRaw(
    GLBaseCiphertextData<GLP> &result, const GLBasePolynomial<GLP> &input,
    const GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key);

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDSmallKeySwitchBase(
    GLBaseCiphertextData<GLP> &result, const GLBasePolynomial<GLP> &input,
    const GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key);

template <class GLP, std::uint32_t LogQ, std::uint32_t InputLogDelta,
          std::uint32_t PlainLogDelta>
inline void GLPlaintextMatrixMultiplyRaw(
    GLRawProductCiphertext<GLP, LogQ, InputLogDelta, PlainLogDelta> &result,
    const GLCiphertext<GLP, LogQ, InputLogDelta> &input,
    const GLPlaintext<GLP, LogQ, PlainLogDelta> &plaintext);

template <class Schedule>
inline void GLXTransformMatrixMultiplyRaw(
    GLRawProductCiphertext<typename Schedule::Parameter, Schedule::input_log_q,
                           Schedule::input_log_delta,
                           Schedule::x_transform_log_scale> &result,
    const typename Schedule::InputCiphertext &input);

template <class GLP, std::uint32_t LogQ, std::uint32_t InputLogDelta,
          std::uint32_t PlainLogDelta>
inline void GLPlaintextHadamardMultiplyRaw(
    GLRawProductCiphertext<GLP, LogQ, InputLogDelta, PlainLogDelta> &result,
    const GLCiphertext<GLP, LogQ, InputLogDelta> &input,
    const GLPlaintext<GLP, LogQ, PlainLogDelta> &plaintext);

template <class GLP, std::uint32_t InputLogQ, std::uint32_t InputLogDelta,
          std::uint32_t DropBits>
inline void GLRescale(
    GLCiphertext<GLP, InputLogQ - DropBits, InputLogDelta - DropBits> &result,
    const GLCiphertext<GLP, InputLogQ, InputLogDelta> &input);

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta,
          std::uint32_t DropBits, std::uint32_t BbarBit>
inline void GLBasePlaintextMultiplyRescale(
    GLBaseCiphertext<GLP, (LhsLogQ < RhsLogQ ? LhsLogQ : RhsLogQ) - DropBits,
                     LhsLogDelta + RhsLogDelta - DropBits> &result,
    const GLBaseCiphertext<GLP, LhsLogQ, LhsLogDelta> &lhs,
    const GLBasePlaintext<GLP, RhsLogQ, RhsLogDelta> &rhs);

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PrimaryBit, std::uint32_t BbarBit>
inline void GLBaseConjugate(
    GLBaseCiphertext<GLP, LogQ, LogDelta> &result,
    const GLBaseCiphertext<GLP, LogQ, LogDelta> &input,
    const GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &key);

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLBaseConjugationKeyGen(
    GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &key,
    const Key<typename GLP::baseP> &secret, CKKSNoise noise);

namespace gl_ship_detail {

template <class GLP>
struct SparseTerm;

template <class Schedule>
inline std::array<SparseTerm<typename Schedule::Parameter>,
                  Schedule::sparse_hamming_weight>
extractSparseTerms(const Key<typename Schedule::Parameter::baseP> &sparse_key);

template <class Schedule>
inline void maskedColumnKeyGen(
    GLSHIPMaskedColumnKey<Schedule> &key,
    const SparseTerm<typename Schedule::Parameter> &term,
    GLSHIPSupportInterval interval,
    const Key<typename Schedule::Parameter::baseP> &dense_key, CKKSNoise noise);

template <class Schedule>
inline void hmuxKeyGen(
    GLSHIPHMuxKey<Schedule> &key,
    const SparseTerm<typename Schedule::Parameter> &term,
    const Key<typename Schedule::Parameter::baseP> &dense_key, CKKSNoise noise);

template <class Schedule>
inline void denseToSparse(typename Schedule::CoefficientCiphertext &result,
                          const typename Schedule::CoefficientCiphertext &input,
                          const GLDDSmallKeySwitchKey<
                              typename Schedule::Parameter, Schedule::q0_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit> &key);

template <class Schedule>
inline void buildBaseFactor(
    GLBaseCiphertext<typename Schedule::Parameter,
                     Schedule::half_bootstrap_log_q, Schedule::tree_log_delta>
        &factor,
    const GLBasePolynomial<typename Schedule::Parameter> &body,
    std::uint32_t channel);

template <class Schedule>
struct MaskedColumnNTTWorkspace;

template <class Schedule>
inline void maskedColumn(
    GLBaseCiphertext<typename Schedule::Parameter,
                     Schedule::half_bootstrap_log_q, Schedule::tree_log_delta>
        &result,
    const GLBasePolynomial<typename Schedule::Parameter> &mask,
    std::uint32_t channel, const GLSHIPMaskedColumnKey<Schedule> &key,
    MaskedColumnNTTWorkspace<Schedule> *workspace = nullptr);

template <class Schedule>
inline void warmMaskedColumnNTTCache(
    const GLSHIPMaskedColumnKey<Schedule> &key);

template <class Schedule>
inline void releaseMaskedColumnNTTCache(
    const GLSHIPMaskedColumnKey<Schedule> &key);

template <class Schedule>
struct HMuxNTTWorkspace;

template <class Schedule>
inline void hmux(GLBaseCiphertext<typename Schedule::Parameter,
                                  Schedule::half_bootstrap_log_q,
                                  Schedule::tree_log_delta> &result,
                 const GLBaseCiphertext<typename Schedule::Parameter,
                                        Schedule::half_bootstrap_log_q,
                                        Schedule::tree_log_delta> &input,
                 const GLSHIPHMuxKey<Schedule> &key,
                 HMuxNTTWorkspace<Schedule> *workspace = nullptr);

template <std::size_t Level, class Schedule>
inline void productTreeLevel(
    GLBaseCiphertextData<typename Schedule::Parameter> &result,
    std::vector<GLBaseCiphertextData<typename Schedule::Parameter>> nodes,
    const GLSHIPProductRelinKeyChain<Schedule> &keys);

template <std::size_t Level, class Schedule>
inline void insertProductTreeBatch(
    std::array<std::unique_ptr<std::vector<
                   GLBaseCiphertextData<typename Schedule::Parameter>>>,
               Schedule::tree_depth + 1> &levels,
    std::vector<GLBaseCiphertextData<typename Schedule::Parameter>> nodes,
    const GLSHIPProductRelinKeyChain<Schedule> &keys);

template <class GLP>
inline void addInPlace(GLBaseCiphertextData<GLP> &destination,
                       const GLBaseCiphertextData<GLP> &term);
template <class GLP, std::uint32_t LogQ>
inline void reduce(GLBaseCiphertextData<GLP> &ciphertext);
template <class GLP>
inline void multiplyByI(GLBaseCiphertextData<GLP> &ciphertext);

template <class GLP>
inline void rotateX(GLBaseCiphertextData<GLP> &result,
                    const GLBaseCiphertextData<GLP> &input,
                    std::uint32_t amount);

template <class GLP>
inline GLBasePolynomial<GLP> ringUnit();

template <class GLP>
inline void multiplyByScalar(GLPolynomial<GLP> &poly, std::uint32_t scalar);

template <class GLP, std::uint32_t LogQ, std::uint32_t DropBits>
inline void rescale(GLCiphertextData<GLP> &result,
                    const GLCiphertextData<GLP> &input);

template <class GLP, std::uint32_t LogQ, std::uint32_t DropBits>
inline void rescaleBase(GLBaseCiphertextData<GLP> &result,
                        const GLBaseCiphertextData<GLP> &input);

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t RhsLogQ,
          std::uint32_t LogScale, std::uint32_t BbarBit>
inline void baseProductRescaleDD(GLBasePolynomial<GLP> &result,
                                 const GLBasePolynomial<GLP> &lhs,
                                 const GLBasePolynomial<GLP> &rhs);

}  // namespace gl_ship_detail

// Algorithm 3, lines 1-14.  The low-level ciphertext is switched to the
// sparse secret, while every returned slice is encrypted under dense_key via
// the encrypted masks, HMux keys, product relin keys, and final conjugation.
template <class Schedule>
inline void GLSHIPHalfBootstrap(
    typename Schedule::OutputCiphertext &result,
    const typename Schedule::CoefficientCiphertext &coefficient_ciphertext,
    const GLSHIPBootstrapKey<Schedule> &bootstrap_key,
    const GLSHIPBootstrapExecutionOptions &execution = {})
{
    using GLP = typename Schedule::Parameter;
    using HMuxSwitch =
        GLDDSmallKeySwitchKey<GLP, Schedule::half_bootstrap_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit>;

    // These few shared base-ring switches are also first reached from inside
    // worker teams, so warm them while their own transform loops can use the
    // complete OpenMP team.
    const auto prepare_shared_switch =
        []<class SwitchKey>(const SwitchKey &switch_key) {
            if constexpr (gl_detail::smallKeySwitchAccumulationNTTPrimeCount<
                              GLP, SwitchKey> != 0)
                gl_detail::prepareSmallKeySwitchNTTCache<GLP, SwitchKey>(
                    switch_key);
        };
    std::apply(
        [&](const auto &...switch_key) {
            (prepare_shared_switch(switch_key), ...);
        },
        bootstrap_key.product_relin_keys.keys);
    prepare_shared_switch(bootstrap_key.output_conjugation_key);

    typename Schedule::CoefficientCiphertext sparse_ciphertext;
    gl_ship_detail::denseToSparse<Schedule>(sparse_ciphertext,
                                            coefficient_ciphertext,
                                            bootstrap_key.dense_to_sparse_key);

    // Traverse sparse terms outside the Y/channel loop.  This reuses one
    // term's MaskedColumn and HMux spectra across every slice and releases
    // them before preparing the next term.  An online balanced tree retains
    // the exact factor pairing of the original per-slice product tree while
    // bounding transient cache memory independently of the sparse weight.
    constexpr std::size_t slice_count =
        static_cast<std::size_t>(GLP::matrix_dimension) * 2;
    using ProductNode = GLBaseCiphertextData<GLP>;
    using ProductBatch = std::vector<ProductNode>;
    std::array<std::unique_ptr<ProductBatch>, Schedule::tree_depth + 1>
        product_levels{};

    ProductBatch base_factors(slice_count);
#pragma omp parallel for schedule(static)
    for (std::size_t slice = 0; slice < slice_count; slice++) {
        const std::uint32_t y = static_cast<std::uint32_t>(slice / 2);
        const std::uint32_t channel = static_cast<std::uint32_t>(slice & 1U);
        GLBaseCiphertext<GLP, Schedule::half_bootstrap_log_q,
                         Schedule::tree_log_delta>
            base_factor;
        gl_ship_detail::buildBaseFactor<Schedule>(
            base_factor, sparse_ciphertext[0][y], channel);
        base_factors[slice] = std::move(base_factor.ct);
    }
    gl_ship_detail::insertProductTreeBatch<0, Schedule>(
        product_levels, std::move(base_factors),
        bootstrap_key.product_relin_keys);

    for (std::size_t e = 0; e < Schedule::sparse_hamming_weight; e++) {
        const auto &masked_key = bootstrap_key.masked_column_keys[e];
        const auto &hmux_key = bootstrap_key.hmux_keys[e];
        gl_ship_detail::warmMaskedColumnNTTCache<Schedule>(masked_key);

        std::vector<const HMuxSwitch *> hmux_switches;
        if constexpr (gl_detail::smallKeySwitchAccumulationNTTPrimeCount<
                          GLP, HMuxSwitch> != 0) {
            hmux_switches.reserve(
                static_cast<std::size_t>(Schedule::hmux_stages) *
                Schedule::hmux_radix * 2);
            for (const auto &stage : hmux_key.stages)
                for (const auto &branch : stage.branches) {
                    hmux_switches.push_back(&branch.body_key);
                    hmux_switches.push_back(&branch.mask_key);
                }
#pragma omp parallel for schedule(dynamic, 1)
            for (std::size_t i = 0; i < hmux_switches.size(); i++)
                gl_detail::prepareSmallKeySwitchNTTCache<GLP, HMuxSwitch>(
                    *hmux_switches[i]);
        }

        ProductBatch sparse_factors(slice_count);
        int maximum_workers = 1;
#ifdef _OPENMP
        maximum_workers = omp_get_max_threads();
#endif
        const int hmux_workers =
            execution.hmux_threads == 0
                ? maximum_workers
                : static_cast<int>(std::min<std::uint32_t>(
                      execution.hmux_threads,
                      static_cast<std::uint32_t>(maximum_workers)));
        const bool use_staged_factors = hmux_workers < maximum_workers;
        if (!use_staged_factors) {
#pragma omp parallel
            {
                gl_ship_detail::MaskedColumnNTTWorkspace<Schedule>
                    masked_workspace;
                gl_ship_detail::HMuxNTTWorkspace<Schedule> hmux_workspace;
#ifdef USE_HEXL
                hmux_workspace.switch_workspace.use_batched_vector_mac =
                    execution.batch_hmux_products;
#endif
                auto selected = std::make_unique<
                    GLBaseCiphertext<GLP, Schedule::half_bootstrap_log_q,
                                     Schedule::tree_log_delta>>();
                auto displaced = std::make_unique<
                    GLBaseCiphertext<GLP, Schedule::half_bootstrap_log_q,
                                     Schedule::tree_log_delta>>();
#pragma omp for schedule(dynamic)
                for (std::size_t slice = 0; slice < slice_count; slice++) {
                    const std::uint32_t y =
                        static_cast<std::uint32_t>(slice / 2);
                    const std::uint32_t channel =
                        static_cast<std::uint32_t>(slice & 1U);
                    gl_ship_detail::maskedColumn<Schedule>(
                        *selected, sparse_ciphertext[1][y], channel,
                        masked_key, &masked_workspace);
                    gl_ship_detail::hmux<Schedule>(
                        *displaced, *selected, hmux_key, &hmux_workspace);
                    sparse_factors[slice] = std::move(displaced->ct);
                }
            }
        }
        else {
            if (execution.factor_tile_size == 0)
                throw std::invalid_argument(
                    "GL SHIP factor tile size must be positive");
            using FactorCiphertext =
                GLBaseCiphertext<GLP, Schedule::half_bootstrap_log_q,
                                 Schedule::tree_log_delta>;
            const std::size_t tile_capacity =
                std::min(execution.factor_tile_size, slice_count);
            std::vector<FactorCiphertext> selected_factors(tile_capacity);
            for (std::size_t tile_begin = 0; tile_begin < slice_count;
                 tile_begin += tile_capacity) {
                const std::size_t tile_count =
                    std::min(tile_capacity, slice_count - tile_begin);
#pragma omp parallel
                {
                    gl_ship_detail::MaskedColumnNTTWorkspace<Schedule>
                        masked_workspace;
#pragma omp for schedule(dynamic)
                    for (std::size_t local = 0; local < tile_count; local++) {
                        const std::size_t slice = tile_begin + local;
                        const std::uint32_t y =
                            static_cast<std::uint32_t>(slice / 2);
                        const std::uint32_t channel =
                            static_cast<std::uint32_t>(slice & 1U);
                        gl_ship_detail::maskedColumn<Schedule>(
                            selected_factors[local], sparse_ciphertext[1][y],
                            channel, masked_key, &masked_workspace);
                    }
                }
#pragma omp parallel num_threads(hmux_workers)
                {
                    gl_ship_detail::HMuxNTTWorkspace<Schedule> hmux_workspace;
#ifdef USE_HEXL
                    hmux_workspace.switch_workspace.use_batched_vector_mac =
                        execution.batch_hmux_products;
#endif
                    auto displaced = std::make_unique<FactorCiphertext>();
#pragma omp for schedule(dynamic)
                    for (std::size_t local = 0; local < tile_count; local++) {
                        gl_ship_detail::hmux<Schedule>(
                            *displaced, selected_factors[local], hmux_key,
                            &hmux_workspace);
                        sparse_factors[tile_begin + local] =
                            std::move(displaced->ct);
                    }
                }
            }
        }

        gl_ship_detail::releaseMaskedColumnNTTCache<Schedule>(masked_key);
        if constexpr (gl_detail::smallKeySwitchAccumulationNTTPrimeCount<
                          GLP, HMuxSwitch> != 0) {
#pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < hmux_switches.size(); i++)
                gl_detail::releaseSmallKeySwitchNTTCache(*hmux_switches[i]);
        }

        gl_ship_detail::insertProductTreeBatch<0, Schedule>(
            product_levels, std::move(sparse_factors),
            bootstrap_key.product_relin_keys);
    }

    for (std::size_t level = 0; level < Schedule::tree_depth; level++)
        if (product_levels[level])
            throw std::logic_error("incomplete GL SHIP online product tree");
    if (!product_levels[Schedule::tree_depth] ||
        product_levels[Schedule::tree_depth]->size() != slice_count)
        throw std::logic_error("invalid GL SHIP online product tree");
    ProductBatch products = std::move(*product_levels[Schedule::tree_depth]);
    product_levels[Schedule::tree_depth].reset();

#pragma omp parallel for schedule(dynamic)
    for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
        std::array<GLBaseCiphertext<GLP, Schedule::output_log_q,
                                    Schedule::tree_log_delta>,
                   2>
            channels{};
        for (std::uint32_t channel = 0; channel < 2; channel++) {
            const std::size_t slice = static_cast<std::size_t>(2) * y + channel;
            GLBaseCiphertext<GLP, Schedule::output_log_q,
                             Schedule::tree_log_delta>
                product;
            product.ct = std::move(products[slice]);
            GLBaseCiphertext<GLP, Schedule::output_log_q,
                             Schedule::tree_log_delta>
                conjugate;
            GLBaseConjugate(conjugate, product,
                            bootstrap_key.output_conjugation_key);
            channels[channel] = std::move(product);
            gl_ship_detail::addInPlace<GLP>(channels[channel].ct, conjugate.ct);
            gl_ship_detail::reduce<GLP, Schedule::output_log_q>(
                channels[channel].ct);
        }

        gl_ship_detail::multiplyByI<GLP>(channels[1].ct);
        gl_ship_detail::addInPlace<GLP>(channels[0].ct, channels[1].ct);
        gl_ship_detail::reduce<GLP, Schedule::output_log_q>(channels[0].ct);
        result[0][y] = std::move(channels[0][0]);
        result[1][y] = std::move(channels[0][1]);
    }
}

// Algorithm 3, lines 15-17.
template <class Schedule>
inline void GLSHIPBootstrap(typename Schedule::OutputCiphertext &result,
                            const typename Schedule::InputCiphertext &input,
                            const GLSHIPBootstrapKey<Schedule> &bootstrap_key,
                            const GLSHIPBootstrapExecutionOptions &execution =
                                {})
{
    typename Schedule::CoefficientCiphertext coefficients;
    GLSHIPSlotsToCoefficients<Schedule>(coefficients, input,
                                        bootstrap_key.stc_key);
    GLSHIPHalfBootstrap<Schedule>(result, coefficients, bootstrap_key,
                                  execution);
}

namespace gl_ship_detail {

template <class Schedule, class IndexSequence>
struct ProductRelinTuple;

template <class Schedule, std::size_t... Is>
struct ProductRelinTuple<Schedule, std::index_sequence<Is...>> {
    using GLP = typename Schedule::Parameter;
    using type = std::tuple<GLHadamardRelinKey<
        GLP,
        Schedule::half_bootstrap_log_q -
            static_cast<std::uint32_t>(Is) * Schedule::tree_log_delta,
        Schedule::primary_bit, Schedule::bbar_bit>...>;
};

}  // namespace gl_ship_detail

template <class Schedule>
struct GLSHIPProductRelinKeyChain {
    using Tuple = typename gl_ship_detail::ProductRelinTuple<
        Schedule, std::make_index_sequence<Schedule::tree_depth>>::type;
    Tuple keys{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(keys);
    }
};

namespace gl_ship_detail {

template <class Schedule, std::size_t... Is>
inline void productRelinKeyChainGen(
    GLSHIPProductRelinKeyChain<Schedule> &chain,
    const Key<typename Schedule::Parameter::baseP> &dense_key,
    std::index_sequence<Is...>)
{
    using GLP = typename Schedule::Parameter;
    (GLHadamardRelinKeyGen(std::get<Is>(chain.keys), dense_key,
                           GLNoiseAtLevel<Schedule::half_bootstrap_log_q -
                                          static_cast<std::uint32_t>(Is) *
                                              Schedule::tree_log_delta +
                                          GLP::auxiliary_log_q>()),
     ...);
}

}  // namespace gl_ship_detail

template <class Schedule>
struct GLSHIPBootstrapKey {
    using GLP = typename Schedule::Parameter;
    using DenseToSparseKey =
        GLDDSmallKeySwitchKey<GLP, Schedule::q0_log_q, Schedule::primary_bit,
                              Schedule::bbar_bit>;
    using OutputConjugationKey =
        GLDDSmallKeySwitchKey<GLP, Schedule::output_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit>;

    GLSHIPSlotsToCoefficientsKey<Schedule> stc_key{};
    DenseToSparseKey dense_to_sparse_key{};
    std::array<GLSHIPMaskedColumnKey<Schedule>, Schedule::sparse_hamming_weight>
        masked_column_keys{};
    std::array<GLSHIPHMuxKey<Schedule>, Schedule::sparse_hamming_weight>
        hmux_keys{};
    GLSHIPProductRelinKeyChain<Schedule> product_relin_keys{};
    OutputConjugationKey output_conjugation_key{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(stc_key, dense_to_sparse_key, masked_column_keys, hmux_keys,
                product_relin_keys, output_conjugation_key);
    }
};

template <class Schedule>
inline void GLSHIPBootstrapKeyGen(
    GLSHIPBootstrapKey<Schedule> &bootstrap_key,
    const Key<typename Schedule::Parameter::baseP> &dense_key,
    const Key<typename Schedule::Parameter::baseP> &sparse_key,
    const std::array<GLSHIPSupportInterval, Schedule::sparse_hamming_weight>
        &intervals)
{
    using GLP = typename Schedule::Parameter;
    const auto terms = gl_ship_detail::extractSparseTerms<Schedule>(sparse_key);

    GLSHIPSlotsToCoefficientsKeyGen(
        bootstrap_key.stc_key, dense_key,
        GLNoiseAtLevel<Schedule::input_log_q + GLP::auxiliary_log_q>());

    const auto dense_secret = gl_detail::keyPolynomial<GLP>(dense_key);
    GLDDSmallKeySwitchKeyGen(
        bootstrap_key.dense_to_sparse_key, dense_secret, sparse_key,
        GLNoiseAtLevel<Schedule::q0_log_q + GLP::auxiliary_log_q>());

    for (std::size_t e = 0; e < terms.size(); e++) {
        gl_ship_detail::maskedColumnKeyGen<Schedule>(
            bootstrap_key.masked_column_keys[e], terms[e], intervals[e],
            dense_key,
            GLNoiseAtLevel<Schedule::half_bootstrap_log_q +
                           GLP::auxiliary_log_q>());
        gl_ship_detail::hmuxKeyGen<Schedule>(
            bootstrap_key.hmux_keys[e], terms[e], dense_key,
            GLNoiseAtLevel<Schedule::half_bootstrap_log_q +
                           GLP::auxiliary_log_q>());
    }

    gl_ship_detail::productRelinKeyChainGen(
        bootstrap_key.product_relin_keys, dense_key,
        std::make_index_sequence<Schedule::tree_depth>{});
    GLBaseConjugationKeyGen(
        bootstrap_key.output_conjugation_key, dense_key,
        GLNoiseAtLevel<Schedule::output_log_q + GLP::auxiliary_log_q>());
}

template <class Schedule>
inline void GLSHIPSaveBootstrapKey(
    const std::filesystem::path &path,
    const GLSHIPBootstrapKey<Schedule> &bootstrap_key)
{
    CKKSSavePortableBinaryAtomic(path, bootstrap_key);
}

template <class Schedule>
inline void GLSHIPLoadBootstrapKey(GLSHIPBootstrapKey<Schedule> &bootstrap_key,
                                   const std::filesystem::path &path)
{
    CKKSLoadPortableBinary(bootstrap_key, path);
}

namespace gl_ship_detail {

template <class Schedule>
inline void denseToSparse(typename Schedule::CoefficientCiphertext &result,
                          const typename Schedule::CoefficientCiphertext &input,
                          const GLDDSmallKeySwitchKey<
                              typename Schedule::Parameter, Schedule::q0_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit> &key)
{
    using GLP = typename Schedule::Parameter;
    GLDDSmallKeySwitch(result.ct, input[1], key);
    gl_detail::addInPlace<GLP>(result[0], input[0]);
    gl_detail::reduce<GLP, Schedule::q0_log_q>(result[0]);
    gl_detail::reduce<GLP, Schedule::q0_log_q>(result[1]);
}

template <class Schedule>
inline std::complex<long double> phaseRoot(
    const typename Schedule::Parameter::T value)
{
    using GLP = typename Schedule::Parameter;
    using P = typename GLP::baseP;
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    const long double centered =
        ckks_detail::levelToLongDouble<P, Schedule::q0_log_q>(value);
    const long double modulus =
        std::ldexp(1.0L, static_cast<int>(Schedule::q0_log_q));
    const long double angle = 2.0L * pi * centered / modulus;
    return {std::cos(angle), std::sin(angle)};
}

template <class Schedule>
inline typename Schedule::Parameter::T extendedCoefficient(
    const GLBasePolynomial<typename Schedule::Parameter> &poly,
    const std::uint32_t gaussian, const std::uint32_t x, const std::uint32_t w)
{
    using GLP = typename Schedule::Parameter;
    if (w == GLP::phi) return typename GLP::T{0};
    return poly[gl_detail::baseIndex<GLP>(gaussian, x, w)];
}

template <class Schedule>
inline typename Schedule::Parameter::T gaussianChannel(
    const typename Schedule::Parameter::T real,
    const typename Schedule::Parameter::T imag, const std::uint32_t phase,
    const std::uint32_t channel)
{
    using GLP = typename Schedule::Parameter;
    using P = typename GLP::baseP;
    typename GLP::T value = 0;
    bool negative = false;
    switch (phase & 3U) {
    case 0: value = channel == 0 ? real : imag; break;
    case 1:
        value = channel == 0 ? imag : real;
        negative = channel == 0;
        break;
    case 2:
        value = channel == 0 ? real : imag;
        negative = true;
        break;
    default:
        value = channel == 0 ? imag : real;
        negative = channel != 0;
        break;
    }
    if (negative) value = typename GLP::T{0} - value;
    return ckks_detail::reduceToLevel<P, Schedule::q0_log_q>(value);
}

template <class Schedule>
inline void buildBaseFactor(
    GLBaseCiphertext<typename Schedule::Parameter,
                     Schedule::half_bootstrap_log_q, Schedule::tree_log_delta>
        &factor,
    const GLBasePolynomial<typename Schedule::Parameter> &body,
    const std::uint32_t channel)
{
    using GLP = typename Schedule::Parameter;
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    const std::complex<long double> multiplier(0, -Schedule::gap / (4.0L * pi));
    GLBaseSlotTable<GLP> values;
    for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++) {
        for (std::uint32_t w = 0; w < GLP::phi; w++) {
            const auto coefficient =
                body[gl_detail::baseIndex<GLP>(channel, row, w)];
            values(w, row) = static_cast<std::complex<double>>(
                multiplier * phaseRoot<Schedule>(coefficient));
        }
    }
    GLBasePlaintext<GLP, Schedule::half_bootstrap_log_q,
                    Schedule::tree_log_delta>
        encoded;
    GLBaseEncode(encoded, values);
    factor[0] = std::move(encoded.poly);
    gl_detail::clear<GLP>(factor[1]);
}

template <class Schedule>
inline void buildCandidatePlaintext(
    GLBasePlaintext<typename Schedule::Parameter,
                    Schedule::half_bootstrap_log_q +
                        Schedule::Parameter::auxiliary_log_q,
                    Schedule::tree_log_delta> &plaintext,
    const GLBasePolynomial<typename Schedule::Parameter> &mask,
    const GLSHIPCandidate candidate, const std::uint32_t channel)
{
    using GLP = typename Schedule::Parameter;
    GLBaseSlotTable<GLP> values;
    for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++) {
        const std::uint32_t source_x =
            (row + GLP::matrix_dimension - candidate.fine_x) %
            GLP::matrix_dimension;
        for (std::uint32_t t = 0; t < GLP::phi; t++) {
            const std::uint32_t first_w = gl_ship_detail::positiveMod(
                static_cast<std::int64_t>(t) - candidate.w,
                GLP::cyclotomic_order);
            const std::uint32_t subtraction_w = gl_ship_detail::positiveMod(
                static_cast<std::int64_t>(GLP::phi) - candidate.w,
                GLP::cyclotomic_order);
            const auto real = ckks_detail::reduceToLevel<typename GLP::baseP,
                                                         Schedule::q0_log_q>(
                extendedCoefficient<Schedule>(mask, 0, source_x, first_w) -
                extendedCoefficient<Schedule>(mask, 0, source_x,
                                              subtraction_w));
            const auto imag = ckks_detail::reduceToLevel<typename GLP::baseP,
                                                         Schedule::q0_log_q>(
                extendedCoefficient<Schedule>(mask, 1, source_x, first_w) -
                extendedCoefficient<Schedule>(mask, 1, source_x,
                                              subtraction_w));
            const auto selected = gaussianChannel<Schedule>(
                real, imag, candidate.gaussian_phase, channel);
            values(t, row) = static_cast<std::complex<double>>(
                phaseRoot<Schedule>(selected));
        }
    }
    GLBaseEncode(plaintext, values);
}

template <class Schedule>
inline void prepareCandidatePhaseRoots(
    std::vector<std::complex<long double>> &roots,
    const GLBasePolynomial<typename Schedule::Parameter> &mask)
{
    using GLP = typename Schedule::Parameter;
    constexpr std::size_t extended_phi = GLP::phi + 1;
    roots.resize(static_cast<std::size_t>(2) * GLP::matrix_dimension *
                 extended_phi);
    const auto index = [](const std::size_t gaussian, const std::size_t x,
                          const std::size_t w) {
        return (gaussian * GLP::matrix_dimension + x) * extended_phi + w;
    };
    for (std::size_t gaussian = 0; gaussian < 2; gaussian++)
        for (std::size_t x = 0; x < GLP::matrix_dimension; x++) {
            for (std::size_t w = 0; w < GLP::phi; w++)
                roots[index(gaussian, x, w)] = phaseRoot<Schedule>(
                    mask[gl_detail::baseIndex<GLP>(gaussian, x, w)]);
            roots[index(gaussian, x, GLP::phi)] = {1, 0};
        }
}

template <class Schedule>
inline void buildCandidatePlaintextFromPhaseRoots(
    GLBasePlaintext<typename Schedule::Parameter,
                    Schedule::half_bootstrap_log_q +
                        Schedule::Parameter::auxiliary_log_q,
                    Schedule::tree_log_delta> &plaintext,
    const std::vector<std::complex<long double>> &roots,
    const GLSHIPCandidate candidate, const std::uint32_t channel)
{
    using GLP = typename Schedule::Parameter;
    constexpr std::size_t extended_phi = GLP::phi + 1;
    const std::size_t expected_size =
        static_cast<std::size_t>(2) * GLP::matrix_dimension * extended_phi;
    if (roots.size() != expected_size || channel >= 2)
        throw std::invalid_argument("invalid GL candidate phase-root table");
    const auto root = [&](const std::size_t gaussian, const std::size_t x,
                          const std::size_t w) -> const auto & {
        return roots[(gaussian * GLP::matrix_dimension + x) * extended_phi + w];
    };

    std::uint32_t source_gaussian = channel;
    bool negative = false;
    switch (candidate.gaussian_phase & 3U) {
    case 0: break;
    case 1:
        source_gaussian = channel == 0 ? 1 : 0;
        negative = channel == 0;
        break;
    case 2: negative = true; break;
    default:
        source_gaussian = channel == 0 ? 1 : 0;
        negative = channel != 0;
        break;
    }

    GLBaseSlotTable<GLP> values;
    for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++) {
        const std::uint32_t source_x =
            (row + GLP::matrix_dimension - candidate.fine_x) %
            GLP::matrix_dimension;
        for (std::uint32_t t = 0; t < GLP::phi; t++) {
            const std::uint32_t first_w = gl_ship_detail::positiveMod(
                static_cast<std::int64_t>(t) - candidate.w,
                GLP::cyclotomic_order);
            const std::uint32_t subtraction_w = gl_ship_detail::positiveMod(
                static_cast<std::int64_t>(GLP::phi) - candidate.w,
                GLP::cyclotomic_order);
            auto value =
                root(source_gaussian, source_x, first_w) *
                std::conj(root(source_gaussian, source_x, subtraction_w));
            if (negative) value = std::conj(value);
            values(t, row) = static_cast<std::complex<double>>(value);
        }
    }
    GLBaseEncode(plaintext, values);
}

template <class Schedule>
inline std::shared_ptr<
    typename GLSHIPMaskedColumnKey<Schedule>::TransientNTTCache>
prepareMaskedColumnNTTCache(const GLSHIPMaskedColumnKey<Schedule> &key)
{
    using GLP = typename Schedule::Parameter;
    if constexpr (!gl_detail::supportsWidePrimeNTT<GLP>) {
        return {};
    }
    else {
        if (key.candidates.size() != key.encrypted_masks.size() ||
            key.candidates.empty())
            throw std::invalid_argument("malformed GL SHIP masked-column key");
        constexpr std::uint32_t key_log_q =
            Schedule::half_bootstrap_log_q + GLP::auxiliary_log_q;
        constexpr std::uint32_t plaintext_bits = Schedule::tree_log_delta + 2;
        constexpr std::size_t maximum_terms =
            2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
        constexpr std::uint32_t convolution_bits =
            gl_detail::ceilLog2(maximum_terms);
        const auto cache = key.transient_ntt_cache;
        std::call_once(cache->initialize_once, [&] {
            const std::uint32_t candidate_bits =
                gl_detail::ceilLog2(key.candidates.size());
            if (plaintext_bits + convolution_bits + candidate_bits >= 122)
                throw std::overflow_error(
                    "GL masked-column exact NTT bound is exhausted");
            cache->chunk_bits =
                122 - plaintext_bits - convolution_bits - candidate_bits;
            cache->chunk_count =
                (key_log_q + cache->chunk_bits - 1) / cache->chunk_bits;

            constexpr std::size_t coefficient_count =
                gl_detail::GLBaseNTTPlan<GLP>::coefficient_count;
            const std::size_t row_count =
                key.candidates.size() * 2 * cache->chunk_count;
            for (auto &spectra : cache->spectra)
                spectra.resize(row_count * coefficient_count);
            using Plan = gl_detail::GLBaseNTTPlan<GLP>;
            const std::array<const Plan *, 2> plans{
                &gl_detail::baseNTTPlan<GLP, 0>(),
                &gl_detail::baseNTTPlan<GLP, 1>()};

#pragma omp parallel
            {
                auto chunk = std::make_unique<GLBasePolynomial<GLP>>();
#pragma omp for collapse(2) schedule(dynamic, 1)
                for (std::size_t prime = 0; prime < 2; prime++) {
                    for (std::size_t row = 0; row < row_count; row++) {
                        const std::size_t candidate =
                            row / (2 * cache->chunk_count);
                        const std::size_t remainder =
                            row % (2 * cache->chunk_count);
                        const std::size_t component =
                            remainder / cache->chunk_count;
                        const std::uint32_t chunk_index =
                            static_cast<std::uint32_t>(remainder %
                                                       cache->chunk_count);
                        const std::uint32_t shift =
                            chunk_index * cache->chunk_bits;
                        for (std::size_t i = 0; i < GLP::baseP::n; i++)
                            (*chunk)[i] = gl_detail::unsignedChunk(
                                key.encrypted_masks[candidate][component][i],
                                shift, cache->chunk_bits);
                        plans[prime]->forward(std::span<std::uint64_t>(
                                                  cache->spectra[prime].data() +
                                                      row * coefficient_count,
                                                  coefficient_count),
                                              *chunk);
                    }
                }
            }
        });
        return cache;
    }
}

template <class Schedule>
inline void warmMaskedColumnNTTCache(const GLSHIPMaskedColumnKey<Schedule> &key)
{
    (void)prepareMaskedColumnNTTCache<Schedule>(key);
}

template <class Schedule>
inline void releaseMaskedColumnNTTCache(
    const GLSHIPMaskedColumnKey<Schedule> &key)
{
    using Cache = typename GLSHIPMaskedColumnKey<Schedule>::TransientNTTCache;
    key.transient_ntt_cache = std::make_shared<Cache>();
}

template <class Schedule>
struct MaskedColumnNTTWorkspace {
    using GLP = typename Schedule::Parameter;
    static constexpr std::uint32_t key_log_q =
        Schedule::half_bootstrap_log_q + GLP::auxiliary_log_q;

    std::array<std::vector<std::uint64_t>, 2> candidate_spectra{};
    std::array<std::vector<std::uint64_t>, 2> unshifted_spectra{};
    std::array<std::vector<std::uint64_t>, 2> residues{};
    std::array<std::vector<std::uint64_t>, 2> accumulators{};
#ifdef USE_HEXL
    std::vector<std::uint64_t> product{};
#else
    std::vector<unsigned __int128> wide{};
#endif
    std::vector<std::complex<long double>> phase_roots{};
    std::unique_ptr<GLBasePolynomial<GLP>> signed_plaintext{};
    std::unique_ptr<GLBasePlaintext<GLP, key_log_q, Schedule::tree_log_delta>>
        candidate{};
    std::unique_ptr<GLBaseCiphertextData<GLP>> accumulated{};
};

template <class Schedule>
inline bool maskedColumnNTT(
    GLBaseCiphertext<typename Schedule::Parameter,
                     Schedule::half_bootstrap_log_q, Schedule::tree_log_delta>
        &result,
    const GLBasePolynomial<typename Schedule::Parameter> &mask,
    const std::uint32_t channel, const GLSHIPMaskedColumnKey<Schedule> &key,
    MaskedColumnNTTWorkspace<Schedule> *provided_workspace = nullptr)
{
    using GLP = typename Schedule::Parameter;
    using P = typename GLP::baseP;
    using T = typename GLP::T;
    if constexpr (!gl_detail::supportsWidePrimeNTT<GLP>) {
        return false;
    }
    else {
        constexpr std::uint32_t key_log_q =
            Schedule::half_bootstrap_log_q + GLP::auxiliary_log_q;
        constexpr std::uint32_t plaintext_bits = Schedule::tree_log_delta + 2;
        constexpr std::size_t coefficient_count =
            gl_detail::GLBaseNTTPlan<GLP>::coefficient_count;
        const auto cache = prepareMaskedColumnNTTCache<Schedule>(key);
        const std::size_t candidate_count = key.candidates.size();
        MaskedColumnNTTWorkspace<Schedule> local_workspace;
        auto &workspace = provided_workspace == nullptr ? local_workspace
                                                        : *provided_workspace;
        auto &candidate_spectra = workspace.candidate_spectra;
        for (auto &spectra : candidate_spectra)
            spectra.resize(candidate_count * coefficient_count);
        auto &unshifted_spectra = workspace.unshifted_spectra;
        for (auto &spectra : unshifted_spectra)
            spectra.resize(coefficient_count);
        const std::array<const gl_detail::GLBaseNTTPlan<GLP> *, 2> plans{
            &gl_detail::baseNTTPlan<GLP, 0>(),
            &gl_detail::baseNTTPlan<GLP, 1>()};

        if (!workspace.signed_plaintext)
            workspace.signed_plaintext =
                std::make_unique<GLBasePolynomial<GLP>>();
        if (!workspace.candidate)
            workspace.candidate = std::make_unique<
                GLBasePlaintext<GLP, key_log_q, Schedule::tree_log_delta>>();
        auto &signed_plaintext = *workspace.signed_plaintext;
        auto &unshifted_candidate = *workspace.candidate;
        prepareCandidatePhaseRoots<Schedule>(workspace.phase_roots, mask);

        // Fine-X candidates in one (W, Gaussian-phase) group are canonical
        // slot rotations of the same encoded polynomial. Encode the unshifted
        // representative once and realize every other member as the exact
        // X -> X^(5^-fine_x) NTT-spectrum permutation.
        std::array<bool, 4 * GLP::phi> encoded_groups{};
        for (std::size_t group_seed = 0; group_seed < candidate_count;
             group_seed++) {
            const auto group = key.candidates[group_seed];
            if (group.fine_x >= GLP::matrix_dimension || group.w >= GLP::phi)
                throw std::invalid_argument(
                    "invalid GL masked-column candidate");
            const std::size_t group_index =
                4 * group.w + (group.gaussian_phase & 3U);
            if (encoded_groups[group_index]) continue;
            encoded_groups[group_index] = true;
            buildCandidatePlaintextFromPhaseRoots<Schedule>(
                unshifted_candidate, workspace.phase_roots,
                {0, group.w, group.gaussian_phase}, channel);

            for (std::size_t i = 0; i < P::n; i++) {
                const auto [negative, magnitude] =
                    ckks_detail::smallSignedMagnitude<P, key_log_q>(
                        unshifted_candidate.poly[i]);
                if (magnitude >= (std::uint64_t{1} << plaintext_bits))
                    throw std::overflow_error(
                        "GL masked-column plaintext exceeds its exact NTT "
                        "bound");
                const std::int64_t value =
                    negative ? -static_cast<std::int64_t>(magnitude)
                             : static_cast<std::int64_t>(magnitude);
                signed_plaintext[i] = gl_detail::signedI128ToTorus<T>(value);
            }
            for (std::size_t prime = 0; prime < 2; prime++)
                plans[prime]->forward(
                    std::span<std::uint64_t>(unshifted_spectra[prime]),
                    signed_plaintext);

            for (std::size_t candidate_index = 0;
                 candidate_index < candidate_count; candidate_index++) {
                const auto descriptor = key.candidates[candidate_index];
                if (descriptor.w != group.w ||
                    (descriptor.gaussian_phase & 3U) !=
                        (group.gaussian_phase & 3U))
                    continue;
                if (descriptor.fine_x >= GLP::matrix_dimension)
                    throw std::invalid_argument(
                        "invalid GL masked-column candidate");
                const std::uint32_t multiplier = gl_detail::powMod(
                    5,
                    (GLP::matrix_dimension - descriptor.fine_x) %
                        GLP::matrix_dimension,
                    4 * GLP::matrix_dimension);
                const auto z_map =
                    gl_detail::baseXAutomorphismSpectrumMap<GLP>(multiplier);
                for (std::size_t prime = 0; prime < 2; prime++)
                    gl_detail::applyBaseXAutomorphismSpectrum<GLP>(
                        std::span<std::uint64_t>(
                            candidate_spectra[prime].data() +
                                candidate_index * coefficient_count,
                            coefficient_count),
                        std::span<const std::uint64_t>(
                            unshifted_spectra[prime]),
                        z_map);
            }
        }

        if (!workspace.accumulated)
            workspace.accumulated =
                std::make_unique<GLBaseCiphertextData<GLP>>();
        auto &accumulated = *workspace.accumulated;
        gl_detail::clear<GLP>(accumulated[0]);
        gl_detail::clear<GLP>(accumulated[1]);
        auto &residues = workspace.residues;
        auto &accumulators = workspace.accumulators;
        for (auto &values : residues) values.resize(coefficient_count);
        for (auto &values : accumulators) values.resize(coefficient_count);
#ifdef USE_HEXL
        auto &product = workspace.product;
        product.resize(coefficient_count);
#else
        auto &wide = workspace.wide;
        wide.resize(coefficient_count);
#endif
        static const modular_ntt::TwoPrimeCRT crt(modular_ntt::wide_primes[0],
                                                  modular_ntt::wide_primes[1]);
        for (std::size_t component = 0; component < 2; component++) {
            for (std::uint32_t chunk = 0; chunk < cache->chunk_count; chunk++) {
                for (std::size_t prime = 0; prime < 2; prime++) {
                    auto &accumulator = accumulators[prime];
                    std::fill(accumulator.begin(), accumulator.end(), 0);
#ifdef USE_HEXL
                    const std::uint64_t modulus = plans[prime]->modulus();
                    for (std::size_t candidate = 0;
                         candidate < candidate_count; candidate++) {
                        const std::size_t key_row =
                            (candidate * 2 + component) * cache->chunk_count +
                            chunk;
                        const std::uint64_t *key_spectrum =
                            cache->spectra[prime].data() +
                            key_row * coefficient_count;
                        const std::uint64_t *plain_spectrum =
                            candidate_spectra[prime].data() +
                            candidate * coefficient_count;
                        intel::hexl::EltwiseMultMod(
                            product.data(), key_spectrum, plain_spectrum,
                            coefficient_count, modulus, 1);
                        intel::hexl::EltwiseAddMod(
                            accumulator.data(), accumulator.data(),
                            product.data(), coefficient_count, modulus);
                    }
#else
                    std::fill(wide.begin(), wide.end(), 0);
                    std::size_t batch_count = 0;
                    const std::uint64_t modulus = plans[prime]->modulus();
                    const auto flush = [&] {
                        for (std::size_t i = 0; i < coefficient_count; i++) {
                            accumulator[i] = modular_ntt::add(
                                accumulator[i],
                                modular_ntt::reduceWide(wide[i], modulus),
                                modulus);
                            wide[i] = 0;
                        }
                        batch_count = 0;
                    };
                    for (std::size_t candidate = 0; candidate < candidate_count;
                         candidate++) {
                        const std::size_t key_row =
                            (candidate * 2 + component) * cache->chunk_count +
                            chunk;
                        const std::uint64_t *key_spectrum =
                            cache->spectra[prime].data() +
                            key_row * coefficient_count;
                        const std::uint64_t *plain_spectrum =
                            candidate_spectra[prime].data() +
                            candidate * coefficient_count;
                        for (std::size_t i = 0; i < coefficient_count; i++)
                            wide[i] += static_cast<unsigned __int128>(
                                           key_spectrum[i]) *
                                       plain_spectrum[i];
                        batch_count++;
                        if (batch_count == 16) flush();
                    }
                    if (batch_count != 0) flush();
#endif
                    plans[prime]->inverse(
                        std::span<std::uint64_t>(residues[prime]),
                        std::span<std::uint64_t>(accumulator));
                }
                const std::uint32_t shift = chunk * cache->chunk_bits;
                for (std::size_t i = 0; i < coefficient_count; i++)
                    accumulated[component][i] +=
                        gl_detail::signedI128ToTorus<T>(crt.reconstructSigned(
                            residues[0][i], residues[1][i]))
                        << shift;
            }
        }
        reduce<GLP, key_log_q>(accumulated);
        rescaleBase<GLP, key_log_q, GLP::auxiliary_log_q>(result.ct,
                                                          accumulated);
        return true;
    }
}

template <class Schedule>
inline void maskedColumn(
    GLBaseCiphertext<typename Schedule::Parameter,
                     Schedule::half_bootstrap_log_q, Schedule::tree_log_delta>
        &result,
    const GLBasePolynomial<typename Schedule::Parameter> &mask,
    const std::uint32_t channel, const GLSHIPMaskedColumnKey<Schedule> &key,
    MaskedColumnNTTWorkspace<Schedule> *workspace)
{
    using GLP = typename Schedule::Parameter;
    constexpr std::uint32_t key_log_q =
        Schedule::half_bootstrap_log_q + GLP::auxiliary_log_q;
    if (key.candidates.size() != key.encrypted_masks.size())
        throw std::invalid_argument("malformed GL SHIP masked-column key");

    if (maskedColumnNTT<Schedule>(result, mask, channel, key, workspace))
        return;

    GLBaseCiphertextData<GLP> accumulated{};
    for (std::size_t i = 0; i < key.candidates.size(); i++) {
        GLBasePlaintext<GLP, key_log_q, Schedule::tree_log_delta> candidate;
        buildCandidatePlaintext<Schedule>(candidate, mask, key.candidates[i],
                                          channel);
        GLBaseCiphertextData<GLP> term{};
        for (std::size_t component = 0; component < 2; component++)
            gl_detail::baseMultiply<GLP>(term[component],
                                         key.encrypted_masks[i][component],
                                         candidate.poly);
        addInPlace<GLP>(accumulated, term);
    }
    reduce<GLP, key_log_q>(accumulated);
    rescaleBase<GLP, key_log_q, GLP::auxiliary_log_q>(result.ct, accumulated);
}

template <class Schedule>
struct HMuxNTTWorkspace {
    using GLP = typename Schedule::Parameter;
    using SwitchKey =
        GLDDSmallKeySwitchKey<GLP, Schedule::half_bootstrap_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit>;
    using Ciphertext = GLBaseCiphertext<GLP, Schedule::half_bootstrap_log_q,
                                        Schedule::tree_log_delta>;

    gl_detail::SmallKeySwitchSumNTTWorkspace<GLP, SwitchKey,
                                             2 * Schedule::hmux_radix>
        switch_workspace{};
    std::unique_ptr<Ciphertext> current{};
    std::unique_ptr<Ciphertext> selected{};
    std::unique_ptr<std::array<GLBaseCiphertextData<GLP>, Schedule::hmux_radix>>
        rotated{};
    std::unique_ptr<GLBaseCiphertextData<GLP>> accumulated{};
};

template <class Schedule>
inline void hmux(GLBaseCiphertext<typename Schedule::Parameter,
                                  Schedule::half_bootstrap_log_q,
                                  Schedule::tree_log_delta> &result,
                 const GLBaseCiphertext<typename Schedule::Parameter,
                                        Schedule::half_bootstrap_log_q,
                                        Schedule::tree_log_delta> &input,
                 const GLSHIPHMuxKey<Schedule> &key,
                 HMuxNTTWorkspace<Schedule> *workspace)
{
    using GLP = typename Schedule::Parameter;
    constexpr std::uint32_t key_log_q =
        Schedule::half_bootstrap_log_q + GLP::auxiliary_log_q;
    using SwitchKey =
        GLDDSmallKeySwitchKey<GLP, Schedule::half_bootstrap_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit>;
    HMuxNTTWorkspace<Schedule> local_workspace;
    auto &scratch = workspace == nullptr ? local_workspace : *workspace;
    if (!scratch.current)
        scratch.current =
            std::make_unique<typename HMuxNTTWorkspace<Schedule>::Ciphertext>();
    if (!scratch.selected)
        scratch.selected =
            std::make_unique<typename HMuxNTTWorkspace<Schedule>::Ciphertext>();
    if (!scratch.rotated)
        scratch.rotated = std::make_unique<
            std::array<GLBaseCiphertextData<GLP>, Schedule::hmux_radix>>();
    if (!scratch.accumulated)
        scratch.accumulated = std::make_unique<GLBaseCiphertextData<GLP>>();
    auto &current = *scratch.current;
    auto &selected = *scratch.selected;
    auto &rotated = *scratch.rotated;
    auto &accumulated = *scratch.accumulated;
    current = input;

    constexpr std::uint32_t prime_count =
        gl_detail::smallKeySwitchAccumulationNTTPrimeCount<GLP, SwitchKey>;
    constexpr std::size_t maximum_terms =
        2 * GLP::matrix_dimension * (2 * GLP::phi - 1);
    constexpr std::uint32_t required_bits =
        SwitchKey::primary_bit + SwitchKey::bbar_bit +
        gl_detail::ceilLog2(maximum_terms) +
        gl_detail::ceilLog2(SwitchKey::primary_rows) +
        gl_detail::ceilLog2(2 * Schedule::hmux_radix);
    constexpr bool fused_ntt_supported =
        (prime_count == 1 && required_bits <= 60) ||
        (prime_count == 2 && required_bits <= 122);
    for (const auto &stage : key.stages) {
        if (stage.branches.size() != Schedule::hmux_radix)
            throw std::invalid_argument("malformed GL SHIP HMux key");
        const std::array<const GLBasePolynomial<GLP> *, 2> sources{&current[0],
                                                                   &current[1]};
        std::array<std::size_t, 2 * Schedule::hmux_radix> source_indices{};
        std::array<std::uint32_t, 2 * Schedule::hmux_radix> x_multipliers{};
        std::array<const SwitchKey *, 2 * Schedule::hmux_radix> switch_keys{};
        for (std::uint32_t digit = 0; digit < Schedule::hmux_radix; digit++) {
            const std::uint32_t desired_displacement =
                (digit * stage.step) % GLP::matrix_dimension;
            const std::uint32_t automorphism_amount =
                (GLP::matrix_dimension - desired_displacement) %
                GLP::matrix_dimension;
            const std::uint32_t multiplier = gl_detail::powMod(
                5, automorphism_amount, 4 * GLP::matrix_dimension);
            source_indices[2 * digit] = 0;
            source_indices[2 * digit + 1] = 1;
            x_multipliers[2 * digit] = multiplier;
            x_multipliers[2 * digit + 1] = multiplier;
            switch_keys[2 * digit] = &stage.branches[digit].body_key;
            switch_keys[2 * digit + 1] = &stage.branches[digit].mask_key;
        }

        if constexpr (fused_ntt_supported) {
            const bool fused =
                gl_detail::accumulateSmallKeySwitchAutomorphismSumNTT<
                    GLP, SwitchKey>(accumulated, sources, source_indices,
                                    x_multipliers, switch_keys,
                                    &scratch.switch_workspace);
            if (!fused)
                throw std::logic_error("GL SHIP HMux lost its exact NTT path");
        }
        else {
            gl_detail::clear<GLP>(accumulated[0]);
            gl_detail::clear<GLP>(accumulated[1]);
            for (std::uint32_t digit = 0; digit < Schedule::hmux_radix;
                 digit++) {
                const std::uint32_t desired_displacement =
                    (digit * stage.step) % GLP::matrix_dimension;
                const std::uint32_t automorphism_amount =
                    (GLP::matrix_dimension - desired_displacement) %
                    GLP::matrix_dimension;
                rotateX<GLP>(rotated[digit], current.ct, automorphism_amount);
                GLBaseCiphertextData<GLP> body_term{};
                GLBaseCiphertextData<GLP> mask_term{};
                GLDDSmallKeySwitchBaseRaw(body_term, rotated[digit][0],
                                          stage.branches[digit].body_key);
                GLDDSmallKeySwitchBaseRaw(mask_term, rotated[digit][1],
                                          stage.branches[digit].mask_key);
                addInPlace<GLP>(body_term, mask_term);
                addInPlace<GLP>(accumulated, body_term);
            }
        }
        reduce<GLP, key_log_q>(accumulated);
        rescaleBase<GLP, key_log_q, GLP::auxiliary_log_q>(selected.ct,
                                                          accumulated);
        current = selected;
    }
    result = current;
}

template <std::size_t Level, class Schedule>
inline void productTreeLevel(
    GLBaseCiphertextData<typename Schedule::Parameter> &result,
    std::vector<GLBaseCiphertextData<typename Schedule::Parameter>> nodes,
    const GLSHIPProductRelinKeyChain<Schedule> &keys)
{
    using GLP = typename Schedule::Parameter;
    constexpr std::uint32_t input_log_q =
        Schedule::half_bootstrap_log_q -
        static_cast<std::uint32_t>(Level) * Schedule::tree_log_delta;
    constexpr std::uint32_t output_log_q =
        input_log_q - Schedule::tree_log_delta;
    if ((nodes.size() & 1U) != 0)
        throw std::invalid_argument("unbalanced GL SHIP product tree");

    std::vector<GLBaseCiphertextData<GLP>> next(nodes.size() / 2);
    for (std::size_t i = 0; i < next.size(); i++) {
        GLBaseCiphertext<GLP, input_log_q, Schedule::tree_log_delta> lhs;
        GLBaseCiphertext<GLP, input_log_q, Schedule::tree_log_delta> rhs;
        lhs.ct = std::move(nodes[2 * i]);
        rhs.ct = std::move(nodes[2 * i + 1]);
        GLBaseCiphertext<GLP, output_log_q, Schedule::tree_log_delta> product;
        GLBaseHadamardMultiply<GLP, input_log_q, Schedule::tree_log_delta,
                               input_log_q, Schedule::tree_log_delta,
                               Schedule::tree_log_delta, Schedule::primary_bit,
                               Schedule::bbar_bit>(product, lhs, rhs,
                                                   std::get<Level>(keys.keys));
        next[i] = std::move(product.ct);
    }

    if constexpr (Level + 1 == Schedule::tree_depth) {
        if (next.size() != 1)
            throw std::logic_error("invalid GL SHIP tree depth");
        result = std::move(next[0]);
    }
    else {
        productTreeLevel<Level + 1, Schedule>(result, std::move(next), keys);
    }
}

template <std::size_t Level, class Schedule>
inline void insertProductTreeBatch(
    std::array<std::unique_ptr<std::vector<
                   GLBaseCiphertextData<typename Schedule::Parameter>>>,
               Schedule::tree_depth + 1> &levels,
    std::vector<GLBaseCiphertextData<typename Schedule::Parameter>> nodes,
    const GLSHIPProductRelinKeyChain<Schedule> &keys)
{
    using GLP = typename Schedule::Parameter;
    static_assert(Level <= Schedule::tree_depth);
    if (!levels[Level]) {
        levels[Level] =
            std::make_unique<std::vector<GLBaseCiphertextData<GLP>>>(
                std::move(nodes));
        return;
    }

    if constexpr (Level == Schedule::tree_depth) {
        throw std::logic_error("overflowed GL SHIP online product tree");
    }
    else {
        auto previous = std::move(levels[Level]);
        if (previous->size() != nodes.size())
            throw std::invalid_argument(
                "mismatched GL SHIP product-tree batches");

        constexpr std::uint32_t input_log_q =
            Schedule::half_bootstrap_log_q -
            static_cast<std::uint32_t>(Level) * Schedule::tree_log_delta;
        constexpr std::uint32_t output_log_q =
            input_log_q - Schedule::tree_log_delta;
        std::vector<GLBaseCiphertextData<GLP>> products(nodes.size());
#pragma omp parallel
        {
            GLBaseHadamardWorkspace<GLP> arithmetic_workspace;
            auto lhs = std::make_unique<
                GLBaseCiphertext<GLP, input_log_q, Schedule::tree_log_delta>>();
            auto rhs = std::make_unique<
                GLBaseCiphertext<GLP, input_log_q, Schedule::tree_log_delta>>();
            auto product =
                std::make_unique<GLBaseCiphertext<GLP, output_log_q,
                                                  Schedule::tree_log_delta>>();
#pragma omp for schedule(static)
            for (std::size_t i = 0; i < nodes.size(); i++) {
                lhs->ct = std::move((*previous)[i]);
                rhs->ct = std::move(nodes[i]);
                GLBaseHadamardMultiply<
                    GLP, input_log_q, Schedule::tree_log_delta, input_log_q,
                    Schedule::tree_log_delta, Schedule::tree_log_delta,
                    Schedule::primary_bit, Schedule::bbar_bit>(
                    *product, *lhs, *rhs, std::get<Level>(keys.keys),
                    &arithmetic_workspace);
                products[i] = std::move(product->ct);
            }
        }

        previous.reset();
        std::vector<GLBaseCiphertextData<GLP>>{}.swap(nodes);
        insertProductTreeBatch<Level + 1, Schedule>(levels, std::move(products),
                                                    keys);
    }
}

}  // namespace gl_ship_detail

template <class Schedule>
struct GLSHIPSlotsToCoefficientsKey {
    using GLP = typename Schedule::Parameter;
    using ConjugateTransposeKey =
        GLConjugateTransposeKey<GLP, Schedule::input_log_q,
                                Schedule::primary_bit, Schedule::bbar_bit>;
    using RotationKey =
        GLBatchRotationKey<GLP, Schedule::input_log_q, Schedule::primary_bit,
                           Schedule::bbar_bit>;

    struct RotationEntry {
        std::uint32_t amount = 0;
        RotationKey key{};

        template <class Archive>
        void serialize(Archive &archive)
        {
            archive(amount, key);
        }
    };

    ConjugateTransposeKey conjugate_transpose_key{};
    std::vector<RotationEntry> w_rotation_keys{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(conjugate_transpose_key, w_rotation_keys);
    }
};

template <class Schedule>
inline void GLSHIPSlotsToCoefficientsKeyGen(
    GLSHIPSlotsToCoefficientsKey<Schedule> &stc_key,
    const Key<typename Schedule::Parameter::baseP> &dense_key,
    const CKKSNoise noise)
{
    using GLP = typename Schedule::Parameter;
    GLConjugateTransposeKeyGen(stc_key.conjugate_transpose_key, dense_key,
                               noise);

    std::set<std::uint32_t> amounts;
    for (std::uint32_t j = 1; j < Schedule::w_baby_step; j++)
        amounts.insert(j % GLP::phi);
    for (std::uint32_t b = 1; b < Schedule::w_giant_steps; b++)
        amounts.insert((b * Schedule::w_baby_step) % GLP::phi);
    amounts.erase(0);

    stc_key.w_rotation_keys.clear();
    stc_key.w_rotation_keys.reserve(amounts.size());
    for (const std::uint32_t amount : amounts) {
        typename GLSHIPSlotsToCoefficientsKey<Schedule>::RotationEntry entry;
        entry.amount = amount;
        GLBatchRotationKeyGen(entry.key, dense_key, amount, noise);
        stc_key.w_rotation_keys.push_back(std::move(entry));
    }
}

namespace gl_ship_detail {

template <class Schedule>
inline const typename GLSHIPSlotsToCoefficientsKey<Schedule>::RotationKey &
findWRotationKey(const GLSHIPSlotsToCoefficientsKey<Schedule> &key,
                 const std::uint32_t amount)
{
    using GLP = typename Schedule::Parameter;
    const std::uint32_t reduced = amount % GLP::phi;
    const auto found = std::find_if(
        key.w_rotation_keys.begin(), key.w_rotation_keys.end(),
        [=](const auto &entry) { return entry.amount == reduced; });
    if (found == key.w_rotation_keys.end())
        throw std::invalid_argument("missing GL SHIP W-rotation key");
    return found->key;
}

template <class Schedule>
inline void buildXTransformPlaintext(
    GLPlaintext<typename Schedule::Parameter, Schedule::input_log_q,
                Schedule::x_transform_log_scale> &plaintext)
{
    using GLP = typename Schedule::Parameter;
    using P = typename GLP::baseP;
    constexpr long double scale =
        std::ldexp(1.0L, Schedule::x_transform_log_scale) /
        static_cast<long double>(GLP::matrix_dimension);
    gl_detail::clear<GLP>(plaintext.poly);
    // The encoded Vandermonde matrix has the closed coefficient form
    //   [X^x Y^y W^0] = r_x^(-y) / n.
    // Constructing it directly avoids materializing a 64 MiB matrix batch and
    // running a general 2-D encoder every time StC is called.
    for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++) {
        const auto inverse_root = std::conj(gl_detail::xRoot<GLP>(row));
        std::complex<long double> value = 1;
        for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++) {
            const auto coefficient = value * scale;
            plaintext.poly[y][gl_detail::baseIndex<GLP>(0, row, 0)] =
                ckks_detail::signedLongDoubleToLevel<P, Schedule::input_log_q>(
                    std::round(coefficient.real()));
            plaintext.poly[y][gl_detail::baseIndex<GLP>(1, row, 0)] =
                ckks_detail::signedLongDoubleToLevel<P, Schedule::input_log_q>(
                    std::round(coefficient.imag()));
            value *= inverse_root;
        }
    }
}

template <class Schedule>
inline std::vector<gl_detail::TraceSmallComplexMultiplier>
buildXTransformTraceMultipliers()
{
    using GLP = typename Schedule::Parameter;
    constexpr std::uint32_t n = GLP::matrix_dimension;
    static_assert(Schedule::x_transform_log_scale < 62);
    constexpr long double coefficient_scale =
        std::ldexp(1.0L, Schedule::x_transform_log_scale) /
        static_cast<long double>(n);

    std::vector<gl_detail::TraceSmallComplexMultiplier> multipliers(
        static_cast<std::size_t>(n) * n);
    const auto checked_i64 = [](const __int128 value) {
        if (value < std::numeric_limits<std::int64_t>::min() ||
            value > std::numeric_limits<std::int64_t>::max())
            throw std::overflow_error(
                "GL X-transform coefficient exceeds int64_t");
        return static_cast<std::int64_t>(value);
    };

    for (std::uint32_t output_y = 0; output_y < n; output_y++) {
        const std::uint32_t x = output_y == 0 ? 0 : n - output_y;
        const auto inverse_root = std::conj(gl_detail::xRoot<GLP>(x));
        std::complex<long double> value = 1;
        for (std::uint32_t z = 0; z < n; z++) {
            const __int128 real = static_cast<__int128>(std::round(
                                      (value * coefficient_scale).real())) *
                                  n;
            const __int128 imag = static_cast<__int128>(std::round(
                                      (value * coefficient_scale).imag())) *
                                  n;

            // traceProduct conjugates the plaintext's I coefficient and,
            // for X^-x with x != 0, contributes the additional factor -I.
            const __int128 transformed_real = x == 0 ? real : -imag;
            const __int128 transformed_imag = x == 0 ? -imag : -real;
            auto &entry =
                multipliers[static_cast<std::size_t>(output_y) * n + z];
            entry.real = checked_i64(transformed_real);
            entry.imag = checked_i64(transformed_imag);
            entry.real_plus_imag =
                checked_i64(transformed_real + transformed_imag);
            value *= inverse_root;
        }
    }
    return multipliers;
}

template <class Schedule>
inline void buildAdjustedWDiagonalPlaintext(
    GLBasePlaintext<typename Schedule::Parameter, Schedule::input_log_q,
                    Schedule::w_transform_log_scale> &plaintext,
    const std::uint32_t t, const std::uint32_t giant_amount)
{
    using GLP = typename Schedule::Parameter;
    GLBaseSlotTable<GLP> values;
    for (std::uint32_t batch = 0; batch < GLP::phi; batch++) {
        const std::uint32_t source_batch =
            (batch + GLP::phi - giant_amount % GLP::phi) % GLP::phi;
        const std::uint32_t eta_exponent =
            gl_detail::batchExponent<GLP>(source_batch);
        const std::uint32_t column = (source_batch + t) % GLP::phi;
        const auto diagonal = gl_detail::wRoot<GLP>(static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(eta_exponent) * column) %
            GLP::cyclotomic_order));
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++)
            values(batch, row) = static_cast<std::complex<double>>(diagonal);
    }
    GLBaseEncode(plaintext, values);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void addInPlace(GLCiphertext<GLP, LogQ, LogDelta> &destination,
                       const GLCiphertext<GLP, LogQ, LogDelta> &term)
{
    gl_detail::addInPlace<GLP>(destination[0], term[0]);
    gl_detail::addInPlace<GLP>(destination[1], term[1]);
    gl_detail::reduce<GLP, LogQ>(destination[0]);
    gl_detail::reduce<GLP, LogQ>(destination[1]);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t InputLogDelta,
          std::uint32_t PlainLogDelta>
inline void basePlaintextHadamardMultiplyRaw(
    GLRawProductCiphertext<GLP, LogQ, InputLogDelta, PlainLogDelta> &result,
    const GLCiphertext<GLP, LogQ, InputLogDelta> &input,
    const GLBasePlaintext<GLP, LogQ, PlainLogDelta> &plaintext)
{
    // A W diagonal constant in both matrix axes encodes entirely in Y^0.
    // Multiplication therefore separates into independent base-ring slices;
    // avoid a full-ring transform and reuse all available cores over slices.
#pragma omp parallel for collapse(2) schedule(static)
    for (std::size_t component = 0; component < 2; component++)
        for (std::uint32_t y = 0; y < GLP::matrix_dimension; y++)
            gl_detail::baseMultiply<GLP>(result[component][y],
                                         input[component][y], plaintext.poly);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class Schedule>
struct FusedWTransformNTTState {
    using GLP = typename Schedule::Parameter;
    static constexpr std::size_t coefficient_count =
        gl_detail::GLBaseNTTPlan<GLP>::coefficient_count;
    static constexpr std::size_t input_row_count =
        static_cast<std::size_t>(Schedule::w_baby_step) * 2 *
        GLP::matrix_dimension;
    std::array<std::vector<std::uint64_t>, 2> input_spectra{};
    std::array<std::vector<std::uint64_t>, 2> diagonal_spectra{};
};

template <class Schedule>
inline bool prepareFusedWTransformNTT(
    FusedWTransformNTTState<Schedule> &state,
    const std::vector<GLRawProductCiphertext<
        typename Schedule::Parameter, Schedule::input_log_q,
        Schedule::input_log_delta, Schedule::x_transform_log_scale>>
        &baby_rotations)
{
    using GLP = typename Schedule::Parameter;
    using P = typename GLP::baseP;
    using T = typename GLP::T;
    constexpr std::size_t n = GLP::matrix_dimension;
    constexpr std::size_t phi = GLP::phi;
    constexpr std::size_t coefficient_count =
        FusedWTransformNTTState<Schedule>::coefficient_count;
    if constexpr (!gl_detail::supportsWidePrimeNTT<GLP> ||
                  Schedule::input_log_q >= 128) {
        return false;
    }
    else {
        if (baby_rotations.size() != Schedule::w_baby_step)
            throw std::invalid_argument(
                "GL fused W transform has the wrong baby-step count");

        std::vector<GLBasePolynomial<GLP>> signed_diagonals(phi);
#pragma omp parallel for schedule(static)
        for (std::uint32_t t = 0; t < phi; t++) {
            GLBasePlaintext<GLP, Schedule::input_log_q,
                            Schedule::w_transform_log_scale>
                diagonal;
            const std::uint32_t giant_amount =
                (t / Schedule::w_baby_step) * Schedule::w_baby_step;
            buildAdjustedWDiagonalPlaintext<Schedule>(diagonal, t,
                                                      giant_amount);
            for (std::size_t i = 0; i < P::n; i++)
                signed_diagonals[t][i] = gl_detail::signedI128ToTorus<T>(
                    ckks_detail::levelToSigned<P, Schedule::input_log_q>(
                        diagonal.poly[i]));
        }

        std::uint32_t diagonal_bits = 0;
        for (const auto &diagonal : signed_diagonals)
            diagonal_bits =
                std::max(diagonal_bits,
                         gl_detail::maxSignedTorusBitWidth<GLP>(diagonal));
        constexpr std::size_t maximum_terms = 2 * n * (2 * GLP::phi - 1);
        const std::uint32_t required_bits =
            Schedule::input_log_q + diagonal_bits +
            gl_detail::ceilLog2(maximum_terms) +
            gl_detail::ceilLog2(Schedule::w_baby_step);
        if (required_bits > 122) return false;

        constexpr std::size_t input_row_count =
            FusedWTransformNTTState<Schedule>::input_row_count;
        for (std::size_t prime = 0; prime < 2; prime++) {
            state.input_spectra[prime].resize(input_row_count *
                                              coefficient_count);
            state.diagonal_spectra[prime].resize(phi * coefficient_count);
        }
        using Plan = gl_detail::GLBaseNTTPlan<GLP>;
        const std::array<const Plan *, 2> plans{
            &gl_detail::baseNTTPlan<GLP, 0>(),
            &gl_detail::baseNTTPlan<GLP, 1>()};

#pragma omp parallel for collapse(4) schedule(static)
        for (std::size_t prime = 0; prime < 2; prime++)
            for (std::size_t baby = 0; baby < Schedule::w_baby_step; baby++)
                for (std::size_t component = 0; component < 2; component++)
                    for (std::size_t y = 0; y < n; y++) {
                        const std::size_t row = (baby * 2 + component) * n + y;
                        plans[prime]->forward(
                            std::span<std::uint64_t>(
                                state.input_spectra[prime].data() +
                                    row * coefficient_count,
                                coefficient_count),
                            baby_rotations[baby][component][y]);
                    }

#pragma omp parallel for collapse(2) schedule(static)
        for (std::size_t prime = 0; prime < 2; prime++)
            for (std::size_t t = 0; t < phi; t++)
                plans[prime]->forward(std::span<std::uint64_t>(
                                          state.diagonal_spectra[prime].data() +
                                              t * coefficient_count,
                                          coefficient_count),
                                      signed_diagonals[t]);
        return true;
    }
}

template <class Schedule>
inline void fusedWTransformGroup(
    GLRawProductCiphertext<typename Schedule::Parameter, Schedule::input_log_q,
                           Schedule::input_log_delta +
                               Schedule::x_transform_log_scale,
                           Schedule::w_transform_log_scale> &group,
    const FusedWTransformNTTState<Schedule> &state,
    const std::uint32_t giant_step)
{
    using GLP = typename Schedule::Parameter;
    using P = typename GLP::baseP;
    using T = typename GLP::T;
    constexpr std::size_t n = GLP::matrix_dimension;
    constexpr std::size_t coefficient_count =
        FusedWTransformNTTState<Schedule>::coefficient_count;
    using Plan = gl_detail::GLBaseNTTPlan<GLP>;
    const std::array<const Plan *, 2> plans{&gl_detail::baseNTTPlan<GLP, 0>(),
                                            &gl_detail::baseNTTPlan<GLP, 1>()};
    static const modular_ntt::TwoPrimeCRT crt(modular_ntt::wide_primes[0],
                                              modular_ntt::wide_primes[1]);
    const std::uint32_t begin = giant_step * Schedule::w_baby_step;
    const std::uint32_t end =
        std::min<std::uint32_t>(GLP::phi, begin + Schedule::w_baby_step);

#pragma omp parallel
    {
        std::array<std::vector<std::uint64_t>, 2> accumulators{
            std::vector<std::uint64_t>(coefficient_count),
            std::vector<std::uint64_t>(coefficient_count)};
        std::array<std::vector<std::uint64_t>, 2> coefficients{
            std::vector<std::uint64_t>(coefficient_count),
            std::vector<std::uint64_t>(coefficient_count)};
#pragma omp for collapse(2) schedule(static)
        for (std::size_t component = 0; component < 2; component++) {
            for (std::size_t y = 0; y < n; y++) {
                for (std::size_t prime = 0; prime < 2; prime++) {
                    auto &accumulator = accumulators[prime];
                    std::fill(accumulator.begin(), accumulator.end(), 0);
                    const std::uint64_t modulus = plans[prime]->modulus();
                    for (std::uint32_t t = begin; t < end; t++) {
                        const std::size_t baby = t - begin;
                        const std::size_t input_row =
                            (baby * 2 + component) * n + y;
                        const std::uint64_t *input_spectrum =
                            state.input_spectra[prime].data() +
                            input_row * coefficient_count;
                        const std::uint64_t *diagonal_spectrum =
                            state.diagonal_spectra[prime].data() +
                            t * coefficient_count;
                        for (std::size_t i = 0; i < coefficient_count; i++) {
                            const std::uint64_t product = modular_ntt::multiply(
                                input_spectrum[i], diagonal_spectrum[i],
                                modulus);
                            accumulator[i] = modular_ntt::add(accumulator[i],
                                                              product, modulus);
                        }
                    }
                    plans[prime]->inverse(
                        std::span<std::uint64_t>(coefficients[prime]),
                        std::span<std::uint64_t>(accumulator));
                }
                for (std::size_t i = 0; i < coefficient_count; i++)
                    group[component][y][i] =
                        ckks_detail::reduceToLevel<P, Schedule::input_log_q>(
                            gl_detail::signedI128ToTorus<T>(
                                crt.reconstructSigned(coefficients[0][i],
                                                      coefficients[1][i])));
            }
        }
    }
}

}  // namespace gl_ship_detail

// Algorithm 1 followed by Algorithm 2.  Both plaintext transforms stay at the
// input modulus and a single final rescale consumes the grouped StC limb.
template <class Schedule>
inline void GLSHIPSlotsToCoefficients(
    typename Schedule::CoefficientCiphertext &result,
    const typename Schedule::InputCiphertext &input,
    const GLSHIPSlotsToCoefficientsKey<Schedule> &stc_key)
{
    using GLP = typename Schedule::Parameter;
    using AfterX = GLRawProductCiphertext<GLP, Schedule::input_log_q,
                                          Schedule::input_log_delta,
                                          Schedule::x_transform_log_scale>;
    using BeforeRescale =
        GLRawProductCiphertext<GLP, Schedule::input_log_q,
                               Schedule::input_log_delta +
                                   Schedule::x_transform_log_scale,
                               Schedule::w_transform_log_scale>;

    // The same transpose key is consumed on both sides of the X transform.
    // Cache its full-ring spectra only for this StC call; keeping this roughly
    // 4.5 GiB n512 cache out of the bootstrap key avoids overlapping it with
    // the later HMux caches.
    using TransposeKey =
        typename GLSHIPSlotsToCoefficientsKey<Schedule>::ConjugateTransposeKey;
    using TransposeCache = gl_detail::BigKeySwitchNTTCache<GLP, TransposeKey>;
    std::unique_ptr<TransposeCache> transpose_cache;
    if constexpr (TransposeCache::prime_count != 0) {
        transpose_cache = std::make_unique<TransposeCache>();
        gl_detail::prepareBigKeySwitchNTTCache<GLP, TransposeKey>(
            *transpose_cache, stc_key.conjugate_transpose_key);
    }

    typename Schedule::InputCiphertext transposed;
    GLConjugateTranspose(transposed, input, stc_key.conjugate_transpose_key,
                         transpose_cache.get());
    AfterX x_product;
    GLXTransformMatrixMultiplyRaw<Schedule>(x_product, transposed);
    AfterX after_x;
    GLConjugateTranspose(after_x, x_product, stc_key.conjugate_transpose_key,
                         transpose_cache.get());
    transpose_cache.reset();

    std::vector<AfterX> baby_rotations(Schedule::w_baby_step);
    baby_rotations[0] = after_x;
    for (std::uint32_t j = 1; j < Schedule::w_baby_step; j++)
        GLRotateBatches(baby_rotations[j], after_x,
                        gl_ship_detail::findWRotationKey<Schedule>(stc_key, j));

    gl_ship_detail::FusedWTransformNTTState<Schedule> w_ntt_state;
    const bool use_fused_w =
        gl_ship_detail::prepareFusedWTransformNTT<Schedule>(w_ntt_state,
                                                            baby_rotations);
    BeforeRescale accumulated;
    bool accumulated_initialized = false;
    for (std::uint32_t b = 0; b < Schedule::w_giant_steps; b++) {
        const std::uint32_t giant_amount = b * Schedule::w_baby_step;
        BeforeRescale group;
        if (use_fused_w) {
            gl_ship_detail::fusedWTransformGroup<Schedule>(group, w_ntt_state,
                                                           b);
        }
        else {
            bool group_initialized = false;
            const std::uint32_t end = std::min<std::uint32_t>(
                GLP::phi, giant_amount + Schedule::w_baby_step);
            for (std::uint32_t t = giant_amount; t < end; t++) {
                const std::uint32_t j = t - giant_amount;
                GLBasePlaintext<GLP, Schedule::input_log_q,
                                Schedule::w_transform_log_scale>
                    diagonal;
                gl_ship_detail::buildAdjustedWDiagonalPlaintext<Schedule>(
                    diagonal, t, giant_amount);
                BeforeRescale term;
                gl_ship_detail::basePlaintextHadamardMultiplyRaw(
                    term, baby_rotations[j], diagonal);
                if (!group_initialized) {
                    group = std::move(term);
                    group_initialized = true;
                }
                else {
                    gl_ship_detail::addInPlace(group, term);
                }
            }
        }

        BeforeRescale shifted;
        if (giant_amount == 0)
            shifted = std::move(group);
        else
            GLRotateBatches(shifted, group,
                            gl_ship_detail::findWRotationKey<Schedule>(
                                stc_key, giant_amount));
        if (!accumulated_initialized) {
            accumulated = std::move(shifted);
            accumulated_initialized = true;
        }
        else {
            gl_ship_detail::addInPlace(accumulated, shifted);
        }
    }

    GLRescale<GLP, Schedule::input_log_q,
              Schedule::input_log_delta + Schedule::stc_log_scale,
              Schedule::stc_log_scale>(result, accumulated);
}

namespace gl_ship_detail {

template <class GLP>
struct SparseTerm {
    std::uint32_t flat_index = 0;
    std::uint32_t x = 0;
    std::uint32_t w = 0;
    std::uint32_t gaussian = 0;
    bool negative = false;
};

template <class Schedule>
inline std::array<SparseTerm<typename Schedule::Parameter>,
                  Schedule::sparse_hamming_weight>
extractSparseTerms(const Key<typename Schedule::Parameter::baseP> &sparse_key)
{
    using GLP = typename Schedule::Parameter;
    using T = typename GLP::T;
    std::array<SparseTerm<GLP>, Schedule::sparse_hamming_weight> terms{};
    std::size_t count = 0;
    const T minus_one = T{0} - T{1};
    for (std::uint32_t f = 0; f < GLP::baseP::n; f++) {
        if (sparse_key[f] == T{0}) continue;
        if (sparse_key[f] != T{1} && sparse_key[f] != minus_one)
            throw std::invalid_argument(
                "GL SHIP sparse key coefficients must be in {0,+1,-1}");
        if (count == terms.size())
            throw std::invalid_argument(
                "GL SHIP sparse key exceeds scheduled Hamming weight");
        auto &term = terms[count++];
        term.flat_index = f;
        term.w = f / (2 * GLP::matrix_dimension);
        const std::uint32_t within_w = f % (2 * GLP::matrix_dimension);
        term.x = within_w / 2;
        term.gaussian = within_w & 1U;
        term.negative = sparse_key[f] == minus_one;
    }
    if (count != terms.size())
        throw std::invalid_argument(
            "GL SHIP sparse key does not match scheduled Hamming weight");
    return terms;
}

template <class Schedule>
inline std::vector<GLSHIPCandidate> buildCandidates(
    const GLSHIPSupportInterval interval)
{
    using GLP = typename Schedule::Parameter;
    if (interval.width == 0 ||
        static_cast<std::uint64_t>(interval.start) + interval.width >
            GLP::baseP::n)
        throw std::invalid_argument("invalid GL SHIP support interval");

    std::set<GLSHIPCandidate> candidates;
    for (std::uint32_t f = interval.start; f < interval.start + interval.width;
         f++) {
        const std::uint32_t w = f / (2 * GLP::matrix_dimension);
        const std::uint32_t within_w = f % (2 * GLP::matrix_dimension);
        const std::uint32_t x = within_w / 2;
        const std::uint32_t gaussian = within_w & 1U;
        const std::uint32_t fine_x = x % Schedule::theta;
        const std::uint32_t coarse_x = x - fine_x;
        for (std::uint32_t sign_phase : {0U, 2U}) {
            for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++) {
                const std::uint32_t wrap =
                    (row < fine_x ? 1U : 0U) +
                    (row >= GLP::matrix_dimension - coarse_x ? 1U : 0U);
                candidates.insert(
                    {fine_x, w, (sign_phase + gaussian + wrap) & 3U});
            }
        }
    }
    return {candidates.begin(), candidates.end()};
}

template <class Schedule>
inline GLSHIPCandidate selectedCandidate(
    const SparseTerm<typename Schedule::Parameter> &term,
    const std::uint32_t row)
{
    using GLP = typename Schedule::Parameter;
    const std::uint32_t fine_x = term.x % Schedule::theta;
    const std::uint32_t coarse_x = term.x - fine_x;
    const std::uint32_t wrap =
        (row < fine_x ? 1U : 0U) +
        (row >= GLP::matrix_dimension - coarse_x ? 1U : 0U);
    return {fine_x, term.w,
            ((term.negative ? 2U : 0U) + term.gaussian + wrap) & 3U};
}

template <class GLP>
inline void setZero(GLBasePolynomial<GLP> &poly)
{
    poly.fill(typename GLP::T{0});
}

}  // namespace gl_ship_detail

template <class Schedule>
struct GLSHIPMaskedColumnKey {
    using GLP = typename Schedule::Parameter;
    static constexpr std::uint32_t key_log_q =
        Schedule::half_bootstrap_log_q + GLP::auxiliary_log_q;

    GLSHIPSupportInterval interval{};
    std::vector<GLSHIPCandidate> candidates{};
    // Selectors are encoded at scale 2^P and encrypted modulo P*Q.
    std::vector<GLBaseCiphertextData<GLP>> encrypted_masks{};

    struct TransientNTTCache {
        std::once_flag initialize_once{};
        std::uint32_t chunk_bits = 0;
        std::uint32_t chunk_count = 0;
        std::array<std::vector<std::uint64_t>, 2> spectra{};
    };
    // The cache is derived from encrypted_masks and intentionally omitted
    // from archives.  Half-bootstrap releases it after finishing this sparse
    // term so it never overlaps every HMux cache in memory.
    mutable std::shared_ptr<TransientNTTCache> transient_ntt_cache =
        std::make_shared<TransientNTTCache>();

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(interval, candidates, encrypted_masks);
    }
};

template <class Schedule>
struct GLSHIPHMuxBranchKey {
    using GLP = typename Schedule::Parameter;
    using SwitchKey =
        GLDDSmallKeySwitchKey<GLP, Schedule::half_bootstrap_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit>;
    SwitchKey body_key{};
    SwitchKey mask_key{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(body_key, mask_key);
    }
};

template <class Schedule>
struct GLSHIPHMuxKey {
    struct Stage {
        std::uint32_t step = 0;
        std::vector<GLSHIPHMuxBranchKey<Schedule>> branches{};

        template <class Archive>
        void serialize(Archive &archive)
        {
            archive(step, branches);
        }
    };
    std::vector<Stage> stages{};

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(stages);
    }
};

namespace gl_ship_detail {

template <class Schedule>
inline void maskedColumnKeyGen(
    GLSHIPMaskedColumnKey<Schedule> &key,
    const SparseTerm<typename Schedule::Parameter> &term,
    const GLSHIPSupportInterval interval,
    const Key<typename Schedule::Parameter::baseP> &dense_key,
    const CKKSNoise noise)
{
    using GLP = typename Schedule::Parameter;
    constexpr std::uint32_t key_log_q =
        GLSHIPMaskedColumnKey<Schedule>::key_log_q;
    if (term.flat_index < interval.start ||
        term.flat_index >= interval.start + interval.width)
        throw std::invalid_argument(
            "GL SHIP support interval does not contain sparse coefficient");

    key.interval = interval;
    key.candidates = buildCandidates<Schedule>(interval);
    key.encrypted_masks.resize(key.candidates.size());
    key.transient_ntt_cache = std::make_shared<
        typename GLSHIPMaskedColumnKey<Schedule>::TransientNTTCache>();

    for (std::size_t candidate_index = 0;
         candidate_index < key.candidates.size(); candidate_index++) {
        GLBaseSlotTable<GLP> logical_mask;
        for (std::uint32_t row = 0; row < GLP::matrix_dimension; row++) {
            const bool selected = key.candidates[candidate_index] ==
                                  selectedCandidate<Schedule>(term, row);
            for (std::uint32_t w = 0; w < GLP::phi; w++)
                logical_mask(w, row) = selected ? 1.0 : 0.0;
        }
        GLBasePlaintext<GLP, key_log_q, GLP::auxiliary_log_q> encoded;
        GLBaseEncode(encoded, logical_mask);
        gl_detail::encryptBaseAtLevel<GLP, key_log_q>(
            key.encrypted_masks[candidate_index], encoded.poly, dense_key,
            noise);
    }
}

template <class Schedule>
inline void hmuxKeyGen(
    GLSHIPHMuxKey<Schedule> &key,
    const SparseTerm<typename Schedule::Parameter> &term,
    const Key<typename Schedule::Parameter::baseP> &dense_key,
    const CKKSNoise noise)
{
    using GLP = typename Schedule::Parameter;
    const auto secret = gl_detail::keyPolynomial<GLP>(dense_key);
    const auto unit = ringUnit<GLP>();
    GLBasePolynomial<GLP> zero{};
    const std::uint32_t coarse_index =
        (term.x - term.x % Schedule::theta) / Schedule::theta;

    key.stages.clear();
    key.stages.reserve(Schedule::hmux_stages);
    std::uint32_t radix_power = 1;
    for (std::uint32_t stage_index = 0; stage_index < Schedule::hmux_stages;
         stage_index++) {
        typename GLSHIPHMuxKey<Schedule>::Stage stage;
        stage.step = Schedule::theta * radix_power;
        stage.branches.resize(Schedule::hmux_radix);
        const std::uint32_t selected_digit =
            (coarse_index / radix_power) % Schedule::hmux_radix;

        for (std::uint32_t digit = 0; digit < Schedule::hmux_radix; digit++) {
            const std::uint32_t desired_displacement =
                (digit * stage.step) % GLP::matrix_dimension;
            // TFHEpp's row automorphism exposes output[r]=input[r+amount],
            // while Equation (10) uses output[r]=input[r-U].
            const std::uint32_t automorphism_amount =
                (GLP::matrix_dimension - desired_displacement) %
                GLP::matrix_dimension;
            const std::uint32_t multiplier = gl_detail::powMod(
                5, automorphism_amount, 4 * GLP::matrix_dimension);
            GLBasePolynomial<GLP> rotated_secret{};
            gl_detail::baseAutomorphism<GLP>(rotated_secret, secret, multiplier,
                                             1);

            const bool selected = digit == selected_digit;
            GLDDSmallKeySwitchKeyGen(stage.branches[digit].body_key,
                                     selected ? unit : zero, dense_key, noise);
            GLDDSmallKeySwitchKeyGen(stage.branches[digit].mask_key,
                                     selected ? rotated_secret : zero,
                                     dense_key, noise);
        }
        key.stages.push_back(std::move(stage));
        radix_power *= Schedule::hmux_radix;
    }
}

}  // namespace gl_ship_detail

// Raw base-ring DD switching keeps the auxiliary P factor and the P*Q
// modulus.  Callers can sum several switch products before one ModDown, just
// as hybrid RNS switching accumulates all partitions under P*Q.
template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDSmallKeySwitchBaseRaw(
    GLBaseCiphertextData<GLP> &result, const GLBasePolynomial<GLP> &input,
    const GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key)
{
    using SwitchKey = GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;
    if (switch_key.data.size() !=
        SwitchKey::primary_rows * SwitchKey::bbar_rows)
        throw std::invalid_argument("uninitialized GL base key-switch key");

    if constexpr (gl_detail::smallKeySwitchAccumulationNTTPrimeCount<
                      GLP, SwitchKey> != 0) {
        using P = typename GLP::baseP;
        auto input_digits = std::make_unique<
            std::array<Polynomial<P>, SwitchKey::primary_rows>>();
        gl_detail::clear<GLP>(result[0]);
        gl_detail::clear<GLP>(result[1]);
        const bool used_ntt =
            gl_detail::accumulateSmallKeySwitchProductsNTT<GLP, SwitchKey>(
                1, switch_key,
                [&](const std::size_t) {
                    ckks_detail::activeBaseDecomposePolynomialRows<
                        P, LogQ, PrimaryBit, SwitchKey::primary_rows>(
                        *input_digits, input);
                },
                [&](const std::uint32_t primary,
                    const std::size_t) -> const GLBasePolynomial<GLP> & {
                    return (*input_digits)[primary];
                },
                [&](const std::size_t component, const std::uint32_t bbar,
                    const std::size_t, const std::size_t coefficient,
                    const typename GLP::T value) {
                    const std::uint32_t shift =
                        (SwitchKey::bbar_rows - bbar - 1) * BbarBit;
                    result[component][coefficient] += value << shift;
                });
        if (!used_ntt)
            throw std::logic_error("eligible GL DD NTT path was not used");
        gl_detail::reduce<GLP, SwitchKey::key_log_q>(result[0]);
        gl_detail::reduce<GLP, SwitchKey::key_log_q>(result[1]);
    }
    else {
        auto input_digits =
            gl_detail::activeDecompose<GLP, LogQ, PrimaryBit>(input);
        std::array<std::vector<GLBasePolynomial<GLP>>, 2> digit_rows{
            std::vector<GLBasePolynomial<GLP>>(SwitchKey::bbar_rows),
            std::vector<GLBasePolynomial<GLP>>(SwitchKey::bbar_rows)};
        GLBasePolynomial<GLP> product{};
        GLBasePolynomial<GLP> key_row{};
        for (std::uint32_t primary = 0; primary < SwitchKey::primary_rows;
             primary++) {
            for (std::uint32_t bbar = 0; bbar < SwitchKey::bbar_rows; bbar++) {
                for (std::size_t component = 0; component < 2; component++) {
                    gl_detail::unpackDigitPolynomial<GLP, BbarBit>(
                        key_row, switch_key.at(primary, bbar)[component]);
                    gl_detail::baseMultiply<GLP>(product, input_digits[primary],
                                                 key_row);
                    gl_detail::addInPlace<GLP>(digit_rows[component][bbar],
                                               product);
                }
            }
        }
        for (std::size_t component = 0; component < 2; component++)
            gl_detail::activeRecombine<GLP, SwitchKey::key_log_q, BbarBit>(
                result[component], digit_rows[component]);
    }
}

// Base-ring form of one complete DD switch.  It evaluates
// input*source_secret under the destination key and removes P once.
template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLDDSmallKeySwitchBase(
    GLBaseCiphertextData<GLP> &result, const GLBasePolynomial<GLP> &input,
    const GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit> &switch_key)
{
    using SwitchKey = GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;
    GLBaseCiphertextData<GLP> raw{};
    GLDDSmallKeySwitchBaseRaw(raw, input, switch_key);
    gl_ship_detail::rescaleBase<GLP, SwitchKey::key_log_q,
                                SwitchKey::auxiliary_log_q>(result, raw);
}

// Raw linear products deliberately retain the plaintext scale.  The grouped
// StC performs both axis transforms and rounds only once at the end.
template <class GLP, std::uint32_t LogQ, std::uint32_t InputLogDelta,
          std::uint32_t PlainLogDelta>
inline void GLPlaintextMatrixMultiplyRaw(
    GLRawProductCiphertext<GLP, LogQ, InputLogDelta, PlainLogDelta> &result,
    const GLCiphertext<GLP, LogQ, InputLogDelta> &input,
    const GLPlaintext<GLP, LogQ, PlainLogDelta> &plaintext)
{
    gl_detail::traceProduct<GLP>(result[0], input[0], plaintext.poly);
    gl_detail::traceProduct<GLP>(result[1], input[1], plaintext.poly);
    gl_ship_detail::multiplyByScalar<GLP>(result[0], GLP::matrix_dimension);
    gl_ship_detail::multiplyByScalar<GLP>(result[1], GLP::matrix_dimension);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class Schedule>
inline void GLXTransformMatrixMultiplyRaw(
    GLRawProductCiphertext<typename Schedule::Parameter, Schedule::input_log_q,
                           Schedule::input_log_delta,
                           Schedule::x_transform_log_scale> &result,
    const typename Schedule::InputCiphertext &input)
{
    using GLP = typename Schedule::Parameter;
    static_assert(Schedule::input_log_q < 128,
                  "the compact X trace currently supports active levels "
                  "below 128 bits");
    const auto multipliers =
        gl_ship_detail::buildXTransformTraceMultipliers<Schedule>();
    gl_detail::traceProductSmallComplex<GLP, Schedule::input_log_q>(
        result.ct, input.ct, multipliers);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t InputLogDelta,
          std::uint32_t PlainLogDelta>
inline void GLPlaintextHadamardMultiplyRaw(
    GLRawProductCiphertext<GLP, LogQ, InputLogDelta, PlainLogDelta> &result,
    const GLCiphertext<GLP, LogQ, InputLogDelta> &input,
    const GLPlaintext<GLP, LogQ, PlainLogDelta> &plaintext)
{
    gl_detail::polynomialMultiply<GLP>(result[0], input[0], plaintext.poly);
    gl_detail::polynomialMultiply<GLP>(result[1], input[1], plaintext.poly);
    gl_detail::reduce<GLP, LogQ>(result[0]);
    gl_detail::reduce<GLP, LogQ>(result[1]);
}

template <class GLP, std::uint32_t InputLogQ, std::uint32_t InputLogDelta,
          std::uint32_t DropBits>
inline void GLRescale(
    GLCiphertext<GLP, InputLogQ - DropBits, InputLogDelta - DropBits> &result,
    const GLCiphertext<GLP, InputLogQ, InputLogDelta> &input)
{
    static_assert(InputLogQ > DropBits && InputLogDelta >= DropBits);
    gl_ship_detail::rescale<GLP, InputLogQ, DropBits>(result.ct, input.ct);
}

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta,
          std::uint32_t DropBits, std::uint32_t BbarBit>
inline void GLBasePlaintextMultiplyRescale(
    GLBaseCiphertext<GLP, (LhsLogQ < RhsLogQ ? LhsLogQ : RhsLogQ) - DropBits,
                     LhsLogDelta + RhsLogDelta - DropBits> &result,
    const GLBaseCiphertext<GLP, LhsLogQ, LhsLogDelta> &lhs,
    const GLBasePlaintext<GLP, RhsLogQ, RhsLogDelta> &rhs)
{
    static_assert(LhsLogDelta + RhsLogDelta >= DropBits);
    for (std::size_t component = 0; component < 2; component++)
        gl_ship_detail::baseProductRescaleDD<GLP, LhsLogQ, RhsLogQ, DropBits,
                                             BbarBit>(result[component],
                                                      lhs[component], rhs.poly);
}

template <class GLP>
struct GLBaseHadamardWorkspace {
    std::unique_ptr<std::array<GLBasePolynomial<GLP>, 3>> tensor{};
    std::unique_ptr<GLBasePolynomial<GLP>> cross_term{};
    gl_detail::BaseCiphertextTensorNTTWorkspace<GLP> tensor_ntt{};
    std::unique_ptr<GLBaseCiphertextData<GLP>> square_term{};
    std::unique_ptr<GLBaseCiphertextData<GLP>> relinearized{};
};

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t LhsLogDelta,
          std::uint32_t RhsLogQ, std::uint32_t RhsLogDelta,
          std::uint32_t DropBits, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLBaseHadamardMultiply(
    GLBaseCiphertext<GLP, (LhsLogQ < RhsLogQ ? LhsLogQ : RhsLogQ) - DropBits,
                     LhsLogDelta + RhsLogDelta - DropBits> &result,
    const GLBaseCiphertext<GLP, LhsLogQ, LhsLogDelta> &lhs,
    const GLBaseCiphertext<GLP, RhsLogQ, RhsLogDelta> &rhs,
    const GLHadamardRelinKey<GLP, (LhsLogQ < RhsLogQ ? LhsLogQ : RhsLogQ),
                             PrimaryBit, BbarBit> &relin_key,
    GLBaseHadamardWorkspace<GLP> *provided_workspace = nullptr)
{
    constexpr std::uint32_t input_log_q = LhsLogQ < RhsLogQ ? LhsLogQ : RhsLogQ;
    constexpr std::uint32_t out_log_q = input_log_q - DropBits;
    static_assert(DropBits > 0 && DropBits < input_log_q);
    GLBaseHadamardWorkspace<GLP> local_workspace;
    auto &workspace =
        provided_workspace == nullptr ? local_workspace : *provided_workspace;
    if (!workspace.tensor)
        workspace.tensor =
            std::make_unique<std::array<GLBasePolynomial<GLP>, 3>>();
    if (!workspace.cross_term)
        workspace.cross_term = std::make_unique<GLBasePolynomial<GLP>>();
    if (!workspace.square_term)
        workspace.square_term = std::make_unique<GLBaseCiphertextData<GLP>>();
    if (!workspace.relinearized)
        workspace.relinearized = std::make_unique<GLBaseCiphertextData<GLP>>();
    auto &tensor = *workspace.tensor;
    auto &cross_term = *workspace.cross_term;
    auto &square_term = *workspace.square_term;
    auto &relinearized = *workspace.relinearized;
    if (!gl_detail::baseCiphertextTensorMultiplyNTT<GLP, input_log_q>(
            tensor, lhs.ct, rhs.ct, &workspace.tensor_ntt)) {
        gl_detail::baseMultiplyAtLevel<GLP, input_log_q>(tensor[0], lhs[0],
                                                         rhs[0]);
        gl_detail::baseMultiplyAtLevel<GLP, input_log_q>(tensor[1], lhs[0],
                                                         rhs[1]);
        gl_detail::baseMultiplyAtLevel<GLP, input_log_q>(cross_term, lhs[1],
                                                         rhs[0]);
        gl_detail::baseMultiplyAtLevel<GLP, input_log_q>(tensor[2], lhs[1],
                                                         rhs[1]);
        gl_detail::addInPlace<GLP>(tensor[1], cross_term);
        gl_detail::reduce<GLP, input_log_q>(tensor[1]);
    }

    GLDDSmallKeySwitchBase(square_term, tensor[2], relin_key);
    relinearized[0] = std::move(tensor[0]);
    relinearized[1] = std::move(tensor[1]);
    gl_detail::addInPlace<GLP>(relinearized[0], square_term[0]);
    gl_detail::addInPlace<GLP>(relinearized[1], square_term[1]);
    gl_ship_detail::reduce<GLP, input_log_q>(relinearized);
    gl_ship_detail::rescaleBase<GLP, input_log_q, DropBits>(result.ct,
                                                            relinearized);
    gl_ship_detail::reduce<GLP, out_log_q>(result.ct);
}

template <class GLP, std::uint32_t LogQ,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
using GLBaseConjugationKey =
    GLDDSmallKeySwitchKey<GLP, LogQ, PrimaryBit, BbarBit>;

template <class GLP, std::uint32_t LogQ, std::uint32_t PrimaryBit,
          std::uint32_t BbarBit>
inline void GLBaseConjugationKeyGen(
    GLBaseConjugationKey<GLP, LogQ, PrimaryBit, BbarBit> &key,
    const Key<typename GLP::baseP> &secret, const CKKSNoise noise)
{
    constexpr std::uint32_t inverse_x = 4 * GLP::matrix_dimension - 1;
    constexpr std::uint32_t inverse_w = GLP::cyclotomic_order - 1;
    const auto destination_secret = gl_detail::keyPolynomial<GLP>(secret);
    GLBasePolynomial<GLP> source_secret{};
    gl_detail::baseAutomorphism<GLP>(source_secret, destination_secret,
                                     inverse_x, inverse_w);
    GLDDSmallKeySwitchKeyGen(key, source_secret, secret, noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta,
          std::uint32_t PrimaryBit, std::uint32_t BbarBit>
inline void GLBaseConjugate(
    GLBaseCiphertext<GLP, LogQ, LogDelta> &result,
    const GLBaseCiphertext<GLP, LogQ, LogDelta> &input,
    const GLBaseConjugationKey<GLP, LogQ, PrimaryBit, BbarBit> &key)
{
    constexpr std::uint32_t inverse_x = 4 * GLP::matrix_dimension - 1;
    constexpr std::uint32_t inverse_w = GLP::cyclotomic_order - 1;
    GLBasePolynomial<GLP> body{};
    GLBasePolynomial<GLP> mask{};
    gl_detail::baseAutomorphism<GLP>(body, input[0], inverse_x, inverse_w);
    gl_detail::baseAutomorphism<GLP>(mask, input[1], inverse_x, inverse_w);
    GLDDSmallKeySwitchBase(result.ct, mask, key);
    gl_detail::addInPlace<GLP>(result[0], body);
    gl_ship_detail::reduce<GLP, LogQ>(result.ct);
}

namespace gl_ship_detail {

template <class GLP>
inline void rotateX(GLBaseCiphertextData<GLP> &result,
                    const GLBaseCiphertextData<GLP> &input,
                    const std::uint32_t amount)
{
    const std::uint32_t multiplier = gl_detail::powMod(
        5, amount % GLP::matrix_dimension, 4 * GLP::matrix_dimension);
    gl_detail::baseAutomorphism<GLP>(result[0], input[0], multiplier, 1);
    gl_detail::baseAutomorphism<GLP>(result[1], input[1], multiplier, 1);
}

template <class GLP>
inline GLBasePolynomial<GLP> ringUnit()
{
    GLBasePolynomial<GLP> result{};
    result[gl_detail::baseIndex<GLP>(0, 0, 0)] = typename GLP::T{1};
    return result;
}

template <class GLP>
inline void multiplyByI(GLBaseCiphertextData<GLP> &ciphertext)
{
    gl_detail::multiplyByIInPlace<GLP>(ciphertext[0]);
    gl_detail::multiplyByIInPlace<GLP>(ciphertext[1]);
}

}  // namespace gl_ship_detail

// A compile-time low-depth schedule.  The level-zero phase modulus is 2^Q0;
// the input carries one grouped StC limb; MaskedColumn/HMux lift directly to
// HalfBootstrapLogQ; and a balanced h+1 product tree consumes TreeLogDelta
// bits per layer.  h+1 must be a power of two as in Algorithm 3.
template <class GLP, std::uint32_t Q0LogQ, std::uint32_t GapLog,
          std::uint32_t XTransformLogScale, std::uint32_t WTransformLogScale,
          std::uint32_t HalfBootstrapLogQ, std::uint32_t TreeLogDelta,
          std::uint32_t SparseHammingWeight, std::uint32_t Theta,
          std::uint32_t HMuxRadix = 4, std::uint32_t WBabyStep = 4,
          std::uint32_t PrimaryBit = GLP::baseP::Bgbit,
          std::uint32_t BbarBit = GLP::baseP::B̅gbit>
struct GLSHIPBootstrapSchedule {
    using Parameter = GLP;
    static constexpr std::uint32_t q0_log_q = Q0LogQ;
    static constexpr std::uint32_t gap_log = GapLog;
    static constexpr std::uint32_t x_transform_log_scale = XTransformLogScale;
    static constexpr std::uint32_t w_transform_log_scale = WTransformLogScale;
    static constexpr std::uint32_t stc_log_scale =
        XTransformLogScale + WTransformLogScale;
    static constexpr std::uint32_t input_log_q = Q0LogQ + stc_log_scale;
    static constexpr std::uint32_t input_log_delta = Q0LogQ - GapLog;
    static constexpr std::uint32_t half_bootstrap_log_q = HalfBootstrapLogQ;
    static constexpr std::uint32_t tree_log_delta = TreeLogDelta;
    static constexpr std::uint32_t sparse_hamming_weight = SparseHammingWeight;
    static constexpr std::uint32_t factor_count = SparseHammingWeight + 1;
    static constexpr std::uint32_t tree_depth =
        gl_ship_detail::ceilLog(factor_count, 2);
    static constexpr std::uint32_t output_log_q =
        HalfBootstrapLogQ - tree_depth * TreeLogDelta;
    static constexpr std::uint32_t output_log_delta = TreeLogDelta;
    static constexpr std::uint32_t theta = Theta;
    static constexpr std::uint32_t hmux_radix = HMuxRadix;
    static constexpr std::uint32_t hmux_stages =
        gl_ship_detail::ceilLog(GLP::matrix_dimension / Theta, HMuxRadix);
    static constexpr std::uint32_t w_baby_step = WBabyStep;
    static constexpr std::uint32_t w_giant_steps =
        (GLP::phi + WBabyStep - 1) / WBabyStep;
    static constexpr std::uint32_t primary_bit = PrimaryBit;
    static constexpr std::uint32_t bbar_bit = BbarBit;
    static constexpr long double gap =
        static_cast<long double>(std::uint64_t{1} << GapLog);

    static_assert(Q0LogQ > GapLog);
    static_assert(Q0LogQ < 63,
                  "the reference root-of-unity table uses long-double phase "
                  "conversion and supports q0 below 2^63");
    static_assert(XTransformLogScale > 0 && WTransformLogScale > 0);
    static_assert(SparseHammingWeight > 0);
    static_assert(gl_ship_detail::isPowerOfTwo(factor_count));
    static_assert(HalfBootstrapLogQ > tree_depth * TreeLogDelta);
    static_assert(TreeLogDelta < output_log_q);
    static_assert(Theta > 0 && GLP::matrix_dimension % Theta == 0);
    static_assert(HMuxRadix >= 2);
    static_assert(WBabyStep > 0 && WBabyStep <= GLP::phi);
    static_assert(input_log_q + GLP::auxiliary_log_q <=
                  std::numeric_limits<typename GLP::T>::digits);
    static_assert(HalfBootstrapLogQ + GLP::auxiliary_log_q <=
                  std::numeric_limits<typename GLP::T>::digits);

    using InputCiphertext = GLCiphertext<GLP, input_log_q, input_log_delta>;
    using CoefficientCiphertext = GLCiphertext<GLP, q0_log_q, input_log_delta>;
    using OutputCiphertext = GLCiphertext<GLP, output_log_q, output_log_delta>;
};

// Fused-DD schedules selected by the companion average-case estimator while
// keeping the paper's total Q, P, StC limb, gap, sparse weight, and theta.
// n256p17 is intentionally omitted: the unpublished per-prime schedule cannot
// be reconstructed to the paper's measured precision with the current model.
using GLSHIP512p17FusedDDSchedule =
    GLSHIPBootstrapSchedule<GL512p17Parameter, 48, 11, 18, 19, 338, 49, 31, 8,
                            4, 4, 85, 16>;
using GLSHIP1024p17FusedDDSchedule =
    GLSHIPBootstrapSchedule<GL1024p17Parameter, 50, 11, 19, 20, 641, 50, 31, 16,
                            4, 4, 16, 7>;

static_assert(GLSHIP512p17FusedDDSchedule::output_log_q == 93);
static_assert(GLSHIP1024p17FusedDDSchedule::output_log_q == 391);

namespace gl_ship_detail {

template <class Schedule>
constexpr std::uint32_t stcRotationKeyCount()
{
    using GLP = typename Schedule::Parameter;
    std::array<bool, GLP::phi> used{};
    for (std::uint32_t j = 1; j < Schedule::w_baby_step; j++)
        used[j % GLP::phi] = true;
    for (std::uint32_t b = 1; b < Schedule::w_giant_steps; b++)
        used[(b * Schedule::w_baby_step) % GLP::phi] = true;
    used[0] = false;
    std::uint32_t count = 0;
    for (const bool present : used)
        if (present) count++;
    return count;
}

template <class Schedule, std::size_t... Is>
constexpr std::uint64_t productRelinPackedPayloadBytes(
    std::index_sequence<Is...>)
{
    using Tuple = typename GLSHIPProductRelinKeyChain<Schedule>::Tuple;
    return (std::uint64_t{0} + ... +
            std::tuple_element_t<Is, Tuple>::packed_payload_bytes);
}

}  // namespace gl_ship_detail

// Paper-profile coefficient payload of the in-memory unseeded representation.
// The published aggregate masked-column count is used here; callers can use
// the overload below after key generation to measure a key with its actual
// support intervals.  Vector, allocator, and archive metadata are excluded.
template <class Schedule>
constexpr std::uint64_t GLSHIPPaperBootstrapKeyPackedPayloadBytes()
{
    using GLP = typename Schedule::Parameter;
    using StCBig =
        GLDDBigKeySwitchKey<GLP, Schedule::input_log_q, Schedule::primary_bit,
                            Schedule::bbar_bit>;
    using StCSmall =
        GLDDSmallKeySwitchKey<GLP, Schedule::input_log_q, Schedule::primary_bit,
                              Schedule::bbar_bit>;
    using DenseToSparse =
        GLDDSmallKeySwitchKey<GLP, Schedule::q0_log_q, Schedule::primary_bit,
                              Schedule::bbar_bit>;
    using HMuxSwitch =
        GLDDSmallKeySwitchKey<GLP, Schedule::half_bootstrap_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit>;
    using OutputConjugation =
        GLDDSmallKeySwitchKey<GLP, Schedule::output_log_q,
                              Schedule::primary_bit, Schedule::bbar_bit>;
    constexpr auto profile = GLSHIPPaperParameterProfile<GLP>{};
    constexpr std::uint64_t masked_column_bytes =
        static_cast<std::uint64_t>(profile.masked_column_count) * 2 *
        GLP::baseP::n * sizeof(typename GLP::T);
    constexpr std::uint64_t hmux_switch_count =
        static_cast<std::uint64_t>(Schedule::sparse_hamming_weight) *
        Schedule::hmux_stages * Schedule::hmux_radix * 2;
    return StCBig::packed_payload_bytes +
           static_cast<std::uint64_t>(
               gl_ship_detail::stcRotationKeyCount<Schedule>()) *
               StCSmall::packed_payload_bytes +
           DenseToSparse::packed_payload_bytes + masked_column_bytes +
           hmux_switch_count * HMuxSwitch::packed_payload_bytes +
           gl_ship_detail::productRelinPackedPayloadBytes<Schedule>(
               std::make_index_sequence<Schedule::tree_depth>{}) +
           OutputConjugation::packed_payload_bytes;
}

static_assert(
    GLSHIPPaperBootstrapKeyPackedPayloadBytes<GLSHIP512p17FusedDDSchedule>() ==
    UINT64_C(8458338304));

template <class Schedule>
inline std::uint64_t GLSHIPBootstrapKeyPackedPayloadBytes(
    const GLSHIPBootstrapKey<Schedule> &bootstrap_key)
{
    using GLP = typename Schedule::Parameter;
    std::uint64_t bytes =
        bootstrap_key.stc_key.conjugate_transpose_key.packedPayloadBytes();
    for (const auto &entry : bootstrap_key.stc_key.w_rotation_keys)
        bytes += entry.key.switch_key.packedPayloadBytes();
    bytes += bootstrap_key.dense_to_sparse_key.packedPayloadBytes();
    for (const auto &masked : bootstrap_key.masked_column_keys)
        bytes += static_cast<std::uint64_t>(masked.encrypted_masks.size()) * 2 *
                 GLP::baseP::n * sizeof(typename GLP::T);
    for (const auto &hmux : bootstrap_key.hmux_keys)
        for (const auto &stage : hmux.stages)
            for (const auto &branch : stage.branches)
                bytes += branch.body_key.packedPayloadBytes() +
                         branch.mask_key.packedPayloadBytes();
    std::apply(
        [&bytes](const auto &...relin) {
            ((bytes += relin.packedPayloadBytes()), ...);
        },
        bootstrap_key.product_relin_keys.keys);
    return bytes + bootstrap_key.output_conjugation_key.packedPayloadBytes();
}

template <class GLP>
class GLBaseSlotTable {
public:
    using value_type = std::complex<double>;

    GLBaseSlotTable()
        : values_(static_cast<std::size_t>(GLP::phi) * GLP::matrix_dimension)
    {
    }

    value_type &operator()(const std::size_t batch, const std::size_t x)
    {
        return values_.at(batch * GLP::matrix_dimension + x);
    }

    const value_type &operator()(const std::size_t batch,
                                 const std::size_t x) const
    {
        return values_.at(batch * GLP::matrix_dimension + x);
    }

private:
    std::vector<value_type> values_;
};

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
struct GLBasePlaintext {
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    GLBasePolynomial<GLP> poly{};
};

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
struct GLBaseCiphertext {
    static constexpr std::uint32_t log_q = LogQ;
    static constexpr std::uint32_t log_delta = LogDelta;
    GLBaseCiphertextData<GLP> ct{};

    GLBasePolynomial<GLP> &operator[](const std::size_t component)
    {
        return ct[component];
    }
    const GLBasePolynomial<GLP> &operator[](const std::size_t component) const
    {
        return ct[component];
    }
};

// Canonical encoding for one (I,X,W) slice.  Natural logical order is
// (W-slot, X-slot), matching the order used by Algorithm 3's candidate tables.
template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLBaseEncode(GLBasePlaintext<GLP, LogQ, LogDelta> &plaintext,
                         const GLBaseSlotTable<GLP> &slots)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t n = GLP::matrix_dimension;
    constexpr std::uint32_t p = GLP::cyclotomic_order;
    constexpr std::uint32_t phi = GLP::phi;
    const long double scale = std::ldexp(1.0L, LogDelta);

    std::vector<std::complex<long double>> w_coefficients(
        static_cast<std::size_t>(n) * phi);
    auto w_coefficient =
        [&](const std::uint32_t x_slot,
            const std::uint32_t w) -> std::complex<long double> & {
        return w_coefficients[static_cast<std::size_t>(x_slot) * phi + w];
    };

    std::array<std::array<std::complex<long double>, phi>, phi> w_weights{};
    for (std::uint32_t batch = 0; batch < phi; batch++) {
        const std::uint32_t exponent = gl_detail::batchExponent<GLP>(batch);
        const auto at_one_weight = gl_detail::wRoot<GLP>(exponent);
        for (std::uint32_t w = 0; w < phi; w++) {
            const std::uint32_t root_exponent =
                (p - static_cast<std::uint32_t>(
                         (static_cast<std::uint64_t>(exponent) * w) % p)) %
                p;
            w_weights[batch][w] =
                (gl_detail::wRoot<GLP>(root_exponent) - at_one_weight) /
                static_cast<long double>(p);
        }
    }
#pragma omp parallel for schedule(static)
    for (std::uint32_t x_slot = 0; x_slot < n; x_slot++)
        for (std::uint32_t w = 0; w < phi; w++) {
            std::complex<long double> coefficient = 0;
            for (std::uint32_t batch = 0; batch < phi; batch++)
                coefficient += static_cast<std::complex<long double>>(
                                   slots(batch, x_slot)) *
                               w_weights[batch][w];
            w_coefficient(x_slot, w) = coefficient;
        }

    gl_detail::clear<GLP>(plaintext.poly);
    const auto &x_plan = gl_detail::inverseXEmbeddingPlan<GLP>();
#pragma omp parallel
    {
        std::vector<std::complex<long double>> input_line(n);
        std::vector<std::complex<long double>> output_line(n);
        std::vector<std::complex<long double>> work;
#pragma omp for schedule(static)
        for (std::uint32_t w = 0; w < phi; w++) {
            for (std::uint32_t x_slot = 0; x_slot < n; x_slot++)
                input_line[x_slot] = w_coefficient(x_slot, w);
            x_plan.apply(output_line, input_line, work);
            for (std::uint32_t x = 0; x < n; x++) {
                const auto coefficient = output_line[x] * scale;
                plaintext.poly[gl_detail::baseIndex<GLP>(0, x, w)] =
                    ckks_detail::signedLongDoubleToLevel<P, LogQ>(
                        std::round(coefficient.real()));
                plaintext.poly[gl_detail::baseIndex<GLP>(1, x, w)] =
                    ckks_detail::signedLongDoubleToLevel<P, LogQ>(
                        std::round(coefficient.imag()));
            }
        }
    }
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLBaseDecode(GLBaseSlotTable<GLP> &slots,
                         const GLBasePlaintext<GLP, LogQ, LogDelta> &plaintext)
{
    using P = typename GLP::baseP;
    constexpr std::uint32_t n = GLP::matrix_dimension;
    constexpr std::uint32_t phi = GLP::phi;
    const long double inverse_scale = std::ldexp(1.0L, -int(LogDelta));

    std::array<std::complex<long double>, n> x_roots{};
    for (std::uint32_t j = 0; j < n; j++) x_roots[j] = gl_detail::xRoot<GLP>(j);

    for (std::uint32_t batch = 0; batch < phi; batch++) {
        const std::uint32_t exponent = gl_detail::batchExponent<GLP>(batch);
        const auto eta = gl_detail::wRoot<GLP>(exponent);
        for (std::uint32_t x_slot = 0; x_slot < n; x_slot++) {
            std::complex<long double> value = 0;
            for (std::uint32_t w = 0; w < phi; w++) {
                for (std::uint32_t x = 0; x < n; x++) {
                    const long double real =
                        ckks_detail::levelToLongDouble<P, LogQ>(
                            plaintext.poly[gl_detail::baseIndex<GLP>(0, x, w)]);
                    const long double imag =
                        ckks_detail::levelToLongDouble<P, LogQ>(
                            plaintext.poly[gl_detail::baseIndex<GLP>(1, x, w)]);
                    value += std::complex<long double>(real, imag) *
                             std::pow(x_roots[x_slot], int(x)) *
                             std::pow(eta, int(w));
                }
            }
            slots(batch, x_slot) =
                static_cast<std::complex<double>>(value * inverse_scale);
        }
    }
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLBaseEncrypt(GLBaseCiphertext<GLP, LogQ, LogDelta> &ciphertext,
                          const GLBasePlaintext<GLP, LogQ, LogDelta> &plaintext,
                          const Key<typename GLP::baseP> &key,
                          const CKKSNoise noise = GLNoiseAtLevel<LogQ>())
{
    gl_detail::encryptBaseAtLevel<GLP, LogQ>(ciphertext.ct, plaintext.poly, key,
                                             noise);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t LogDelta>
inline void GLBaseDecrypt(
    GLBasePlaintext<GLP, LogQ, LogDelta> &plaintext,
    const GLBaseCiphertext<GLP, LogQ, LogDelta> &ciphertext,
    const Key<typename GLP::baseP> &key)
{
    const auto secret = gl_detail::keyPolynomial<GLP>(key);
    GLBasePolynomial<GLP> product{};
    gl_detail::baseMultiply<GLP>(product, ciphertext[1], secret);
    for (std::uint32_t i = 0; i < GLP::baseP::n; i++)
        plaintext.poly[i] =
            ckks_detail::reduceToLevel<typename GLP::baseP, LogQ>(
                ciphertext[0][i] + product[i]);
}

namespace gl_ship_detail {

template <class GLP, std::uint32_t LogQ>
inline void reduce(GLBaseCiphertextData<GLP> &ciphertext)
{
    gl_detail::reduce<GLP, LogQ>(ciphertext[0]);
    gl_detail::reduce<GLP, LogQ>(ciphertext[1]);
}

template <class GLP>
inline void addInPlace(GLBaseCiphertextData<GLP> &destination,
                       const GLBaseCiphertextData<GLP> &term)
{
    gl_detail::addInPlace<GLP>(destination[0], term[0]);
    gl_detail::addInPlace<GLP>(destination[1], term[1]);
}

template <class GLP>
inline void multiplyByScalar(GLPolynomial<GLP> &poly,
                             const std::uint32_t scalar)
{
    for (auto &slice : poly)
        for (auto &coefficient : slice) coefficient *= scalar;
}

template <class GLP, std::uint32_t LogQ, std::uint32_t DropBits>
inline void rescale(GLCiphertextData<GLP> &result,
                    const GLCiphertextData<GLP> &input)
{
    result = input;
    gl_detail::divideRoundLevel<GLP, LogQ, DropBits>(result[0]);
    gl_detail::divideRoundLevel<GLP, LogQ, DropBits>(result[1]);
}

template <class GLP, std::uint32_t LogQ, std::uint32_t DropBits>
inline void rescaleBase(GLBaseCiphertextData<GLP> &result,
                        const GLBaseCiphertextData<GLP> &input)
{
    static_assert(LogQ > DropBits);
    result = input;
    for (std::size_t component = 0; component < 2; component++)
        for (auto &coefficient : result[component])
            coefficient =
                gl_detail::divideRoundLevel<GLP, LogQ, DropBits>(coefficient);
}

template <class GLP, std::uint32_t LhsLogQ, std::uint32_t RhsLogQ,
          std::uint32_t LogScale, std::uint32_t BbarBit>
inline void baseProductRescaleDD(GLBasePolynomial<GLP> &result,
                                 const GLBasePolynomial<GLP> &lhs,
                                 const GLBasePolynomial<GLP> &rhs)
{
    using P = typename GLP::baseP;
    using T = typename P::T;
    constexpr std::uint32_t lhs_rows = (LhsLogQ + BbarBit - 1) / BbarBit;
    constexpr std::uint32_t rhs_rows = (RhsLogQ + BbarBit - 1) / BbarBit;
    constexpr std::uint32_t base_log_q = LhsLogQ < RhsLogQ ? LhsLogQ : RhsLogQ;
    constexpr std::uint32_t out_log_q = base_log_q - LogScale;
    static_assert(LogScale > 0 && LogScale < base_log_q);
    static_assert(lhs_rows * BbarBit <= std::numeric_limits<T>::digits);
    static_assert(rhs_rows * BbarBit <= std::numeric_limits<T>::digits);
    static_assert(2 * BbarBit + 2 < 63);

    auto lhs_digits = gl_detail::activeDecompose<GLP, LhsLogQ, BbarBit>(lhs);
    auto rhs_digits = gl_detail::activeDecompose<GLP, RhsLogQ, BbarBit>(rhs);
    auto digit_product = std::make_unique<GLBasePolynomial<GLP>>();

    auto accumulate = [&](auto &accumulators, auto add_shifted) {
        for (std::uint32_t i = 0; i < lhs_rows; i++) {
            const int lhs_shift =
                static_cast<int>((lhs_rows - i - 1) * BbarBit);
            for (std::uint32_t j = 0; j < rhs_rows; j++) {
                const int rhs_shift =
                    static_cast<int>((rhs_rows - j - 1) * BbarBit);
                gl_detail::baseMultiply<GLP>(*digit_product, lhs_digits[i],
                                             rhs_digits[j]);
                for (std::uint32_t coefficient = 0; coefficient < P::n;
                     coefficient++)
                    add_shifted(accumulators[coefficient],
                                gl_detail::torusToSignedSmall<P>(
                                    (*digit_product)[coefficient]),
                                lhs_shift + rhs_shift);
            }
        }
    };

    if constexpr (std::is_same_v<T, __uint128_t>) {
        std::vector<Wide384> accumulators(P::n);
        accumulate(accumulators, [](Wide384 &accumulator,
                                    const std::int64_t value, const int shift) {
            accumulator.add_shifted(static_cast<__int128_t>(value), shift);
        });
        for (std::uint32_t i = 0; i < P::n; i++)
            result[i] = ckks_detail::reduceToLevel<P, out_log_q>(
                ckks_detail::rescaleAccumulator<P, LogScale>(accumulators[i]));
    }
    else {
        static_assert(is_multilimb_uint_v<T>);
        using Accumulator = WideSignedLimbAccumulator<2 * T::limbs + 2>;
        std::vector<Accumulator> accumulators(P::n);
        accumulate(accumulators, [](Accumulator &accumulator,
                                    const std::int64_t value, const int shift) {
            accumulator.add_shifted_i64(value, shift);
        });
        const T divisor = T{1} << LogScale;
        for (std::uint32_t i = 0; i < P::n; i++) {
            accumulators[i].add_shifted_i64(
                accumulators[i].is_negative() ? -1 : 1,
                static_cast<int>(LogScale) - 1);
            result[i] = ckks_detail::reduceToLevel<P, out_log_q>(
                accumulators[i].template div_to_torus<T::limbs>(divisor));
        }
    }
}

}  // namespace gl_ship_detail

}  // namespace TFHEpp

#endif  // TFHEPP_DEFAULT_128BIT_PARAMS
