#pragma once

#include "narrow_bootstrapping.hpp"
#include <map>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {

// Public order-two plan. Pair entries index the unweighted roots, never a
// previously reconstructed product; consequently reconstruction has depth one.
class CKKSNarrowFrequencyBasis {
    std::vector<std::uint32_t> frequencies_;
    std::vector<std::array<std::size_t, 2>> expressions_;
public:
    CKKSNarrowFrequencyBasis(std::uint32_t degree, std::vector<std::uint32_t> frequencies)
        : frequencies_(std::move(frequencies))
    {
        if (degree == 0 || degree > 40 || frequencies_.empty() ||
            frequencies_.front() == 0 || frequencies_.back() > degree ||
            !std::is_sorted(frequencies_.begin(), frequencies_.end()) ||
            std::adjacent_find(frequencies_.begin(), frequencies_.end()) != frequencies_.end())
            throw std::invalid_argument("invalid positive narrow frequency basis");
        for (std::uint32_t k = 1; k <= degree; k++) {
            const auto direct = std::lower_bound(frequencies_.begin(), frequencies_.end(), k);
            if (direct != frequencies_.end() && *direct == k) {
                const auto index = static_cast<std::size_t>(direct - frequencies_.begin());
                expressions_.push_back({index, frequencies_.size()});
                continue;
            }
            bool found = false;
            for (std::size_t a = 0; a < frequencies_.size() && !found; a++)
                for (std::size_t b = a; b < frequencies_.size(); b++)
                    if (frequencies_[a] + frequencies_[b] == k) {
                        expressions_.push_back({a, b});
                        found = true;
                        break;
                    }
            if (!found) throw std::invalid_argument("narrow frequency basis leaves a gap");
        }
    }
    std::size_t degree() const { return expressions_.size(); }
    const auto &frequencies() const { return frequencies_; }
    const auto &expression(std::size_t k) const { return expressions_.at(k - 1); }
};

inline std::vector<std::uint32_t> CKKSNarrowFrequencyBasisUnion(
    std::span<const CKKSNarrowFrequencyBasis> plans)
{
    std::vector<std::uint32_t> frequencies;
    for (const auto &plan : plans)
        frequencies.insert(frequencies.end(), plan.frequencies().begin(), plan.frequencies().end());
    std::sort(frequencies.begin(), frequencies.end());
    frequencies.erase(std::unique(frequencies.begin(), frequencies.end()), frequencies.end());
    return frequencies;
}

// Evaluate shared roots for several series on the same coefficient inputs.
// Root ownership is per input; only equal frequencies on that input are reused.
// The caller controls the column tile, just as for ordinary blind rotation.
template <class Narrow, std::uint32_t InputLogDelta>
inline void CKKSNarrowBlindRotateBasisPlansBatch(
    std::vector<std::map<std::uint32_t, typename Narrow::Output>> &results,
    std::span<const CKKSCiphertext<typename Narrow::Param, Narrow::input_log_q, InputLogDelta>> inputs,
    std::span<const CKKSNarrowFrequencyBasis> plans,
    const CKKSNarrowBootstrapKey<Narrow> &key, std::size_t tile_size = 8)
{
    if (tile_size == 0) throw std::invalid_argument("narrow basis union tile size must be positive");
    const auto frequencies = CKKSNarrowFrequencyBasisUnion(plans);
    if (!frequencies.empty() && inputs.size() > std::numeric_limits<std::size_t>::max() / frequencies.size())
        throw std::length_error("narrow basis union batch is too large");
    std::vector<std::map<std::uint32_t, typename Narrow::Output>> evaluated(inputs.size());
    if (!inputs.empty() && !frequencies.empty()) {
        std::vector<CKKSNarrowBlindRotationRequest<Narrow, InputLogDelta>> requests;
        requests.reserve(inputs.size() * frequencies.size());
        for (const auto &input : inputs)
            for (const auto frequency : frequencies)
                requests.push_back({&input, frequency, {1.0, 0.0}});
        std::vector<typename Narrow::Output> roots;
        CKKSNarrowBlindRotateBatch<Narrow, InputLogDelta>(roots, requests, key, tile_size);
        for (std::size_t i = 0; i < roots.size(); ++i)
            evaluated[i / frequencies.size()].emplace(requests[i].frequency, std::move(roots[i]));
    }
    results.swap(evaluated);
}

template <class Narrow, std::uint32_t CoefficientBits = 32>
struct CKKSNarrowFrequencyBasisSchedule {
    using Base = Narrow;
    using Param = typename Base::Param;
    static constexpr auto coefficient_bits = CoefficientBits;
    static constexpr auto root_log_q = Base::output_log_q;
    static constexpr auto root_log_delta = Base::log_delta;
    static constexpr auto output_log_q = root_log_q - root_log_delta;
    static constexpr auto log_delta = root_log_delta + CoefficientBits;
    static_assert(CoefficientBits > 0 && CoefficientBits <= 40);
    static_assert(root_log_q > 2 * root_log_delta);
    static_assert(output_log_q > log_delta);
    using Root = typename Base::Output;
    using Product = CKKSCiphertext<Param, output_log_q, root_log_delta>;
    using Output = CKKSCiphertext<Param, output_log_q, log_delta>;
};

// Additional public keys, at the levels actually consumed by reconstruction.
// The extra product always relinearizes before rescaling using DD.
template <class Schedule>
struct CKKSNarrowFrequencyBasisKey {
    using P = typename Schedule::Param;
    CKKSCompactDDRelinKey<P, Schedule::root_log_q> product;
    CKKSAutoKey<P, Schedule::output_log_q> conjugation;
    template <class Archive>
    void serialize(Archive &archive)
    {
        const std::array<std::uint64_t, 5> expected{
            0x434b4b5346425331ULL, P::n, Schedule::root_log_q,
            Schedule::output_log_q, Schedule::log_delta};
        auto metadata = expected;
        archive(metadata);
        if (metadata != expected) throw std::runtime_error("narrow frequency basis key metadata mismatch");
        archive(product, conjugation);
    }
};

template <class Schedule>
inline void CKKSNarrowFrequencyBasisKeyGen(CKKSNarrowFrequencyBasisKey<Schedule> &result,
    const Key<typename Schedule::Param> &evaluation, CKKSNoise noise = {Schedule::Param::α, 0})
{
    using P = typename Schedule::Param;
    CKKSCompactDDRelinKeyGen(result.product, evaluation, noise);
    CKKSAutoKeyGen<P, Schedule::output_log_q>(result.conjugation, 2 * P::n - 1, evaluation, noise);
}

namespace narrow_detail {
template <std::uint32_t Bits>
inline std::int64_t frequencyCoefficient(double value)
{
    const auto rounded = std::round(std::ldexp(static_cast<long double>(value), Bits));
    // The exclusive upper bound remains exact even when long double == double.
    if (!std::isfinite(rounded) || rounded < -std::ldexp(1.0L, 63) ||
        rounded >= std::ldexp(1.0L, 63))
        throw std::invalid_argument("narrow frequency coefficient exceeds fixed-point range");
    return static_cast<std::int64_t>(rounded);
}

// All canonical slot exponents are powers of five, hence X^(N/2) evaluates
// to +i in every slot. Multiply by r+i*s exactly in the coefficient ring.
// Reinterpret the resulting ciphertext at scale Delta+CoefficientBits without
// another rescale or modulus loss. Inputs and accumulator must be distinct.
template <class Schedule>
inline void addFrequencyWeight(typename Schedule::Output &acc,
    const typename Schedule::Product &root, std::int64_t real, std::int64_t imag)
{
    using P = typename Schedule::Param;
    for (std::size_t c = 0; c <= P::k; c++)
        for (std::size_t i = 0; i < P::n; i++) {
            const auto imaginary = i < P::n / 2
                ? -root.ct[c][i + P::n / 2] : root.ct[c][i - P::n / 2];
            acc.ct[c][i] = ckks_detail::reduceToLevel<P, Schedule::output_log_q>(
                acc.ct[c][i] + root.ct[c][i] * real + imaginary * imag);
        }
}
}  // namespace narrow_detail

// Reconstruct a series from unweighted encrypted roots already evaluated for
// this plan. Callers may reuse these roots for several coefficient vectors;
// roots and the public basis key must belong to the same encryption key.
template <class Schedule>
inline void CKKSNarrowEvalSeriesFromBasisRoots(
    typename Schedule::Output &result,
    std::span<const typename Schedule::Root> roots,
    const CKKSNarrowFrequencyBasisKey<Schedule> &basis_key,
    const CKKSNarrowFrequencyBasis &plan,
    const std::vector<std::complex<double>> &coefficients,
    bool use_fft = true)
{
    using P = typename Schedule::Param;
    if (roots.size() != plan.frequencies().size() || coefficients.size() != plan.degree()+1 ||
        coefficients[0].imag() != 0)
        throw std::invalid_argument("invalid retained narrow basis roots or coefficients");
    std::vector<std::array<std::int64_t, 2>> weights;
    for (const auto coefficient : coefficients)
        weights.push_back({narrow_detail::frequencyCoefficient<Schedule::coefficient_bits>(coefficient.real()),
                           narrow_detail::frequencyCoefficient<Schedule::coefficient_bits>(coefficient.imag())});
    for (auto &component : result.ct) component.fill(0);
    auto product = std::make_unique<typename Schedule::Product>();
    std::unique_ptr<CKKSDDRelinKeyFFT<P, Schedule::root_log_q>> prepared_product;
    if constexpr (use_multilimb_digit_fft_v<P>)
        if (use_fft) prepared_product = makeCKKSDDRelinKeyFFT(basis_key.product);
    for (std::size_t k = 1; k <= plan.degree(); k++) {
        if (weights[k][0] == 0 && weights[k][1] == 0) continue;
        const auto expression = plan.expression(k);
        if (expression[1] == roots.size()) CKKSLevelReduce<P>(*product, roots[expression[0]]);
        else {
            auto multiply = [&](const auto &key) {
                CKKSMultRelinThenRescale<P>(*product, roots[expression[0]], roots[expression[1]], key);
            };
            if constexpr (use_multilimb_digit_fft_v<P>) {
                if (prepared_product) multiply(*prepared_product);
                else multiply(basis_key.product);
            }
            else multiply(basis_key.product);
        }
        narrow_detail::addFrequencyWeight<Schedule>(result, *product, weights[k][0], weights[k][1]);
    }
    prepared_product.reset();
    auto negative = std::make_unique<typename Schedule::Output>();
    if constexpr (use_multilimb_digit_fft_v<P>) {
        if (use_fft) {
            auto prepared = makeCKKSPreparedAutoKeyFFT<P, Schedule::output_log_q>(basis_key.conjugation);
            CKKSEvalAuto<P, Schedule::output_log_q>(negative->ct, result.ct, 2 * P::n - 1, *prepared);
        }
        else CKKSEvalAuto<P, Schedule::output_log_q>(negative->ct, result.ct, 2 * P::n - 1, basis_key.conjugation);
    }
    else CKKSEvalAuto<P, Schedule::output_log_q>(negative->ct, result.ct, 2 * P::n - 1, basis_key.conjugation);
    CKKSAddInPlace<P>(result, *negative);
    result.ct[P::k][0] = ckks_detail::reduceToLevel<P, Schedule::output_log_q>(
        result.ct[P::k][0] + (P::T::from_signed_i64(weights[0][0]) << Schedule::root_log_delta));
}

// Section7 optional basis path: unweighted blind rotations, one extra DD
// multiplication level, THEN complex coefficients, conjugation and c0.
// Secret keys are used only by key generation, never by this evaluator.
template <class Schedule, std::uint32_t InputLogDelta>
inline void CKKSNarrowEvalSeriesWithBasisBatch(
    std::vector<typename Schedule::Output> &results,
    std::span<const CKKSCiphertext<typename Schedule::Param,
        Schedule::Base::input_log_q, InputLogDelta>> inputs,
    const CKKSNarrowBootstrapKey<typename Schedule::Base> &blind_key,
    const CKKSNarrowFrequencyBasisKey<Schedule> &basis_key,
    const CKKSNarrowFrequencyBasis &plan,
    const std::vector<std::complex<double>> &coefficients,
    std::size_t tile_size = 8, bool use_fft = true)
{
    using Base = typename Schedule::Base;
    if (tile_size == 0 || coefficients.size() != plan.degree() + 1 ||
        coefficients[0].imag() != 0)
        throw std::invalid_argument("invalid narrow basis series");
    // Preserve coefficient validation even for an empty input batch.
    for (const auto coefficient : coefficients) {
        narrow_detail::frequencyCoefficient<Schedule::coefficient_bits>(coefficient.real());
        narrow_detail::frequencyCoefficient<Schedule::coefficient_bits>(coefficient.imag());
    }
    results.resize(inputs.size());
    const auto frequencies = plan.frequencies().size();
    if (inputs.size() > std::numeric_limits<std::size_t>::max() / frequencies)
        throw std::length_error("narrow basis batch is too large");
    const auto total = inputs.size() * frequencies;
    std::vector<typename Schedule::Root> pending;
    pending.reserve(frequencies);
    // Tile across input boundaries so real/imaginary halves and independent
    // ciphertexts can share W_x. Retain only a tile and one unfinished basis.
    for (std::size_t first = 0; first < total;) {
        const auto count = std::min(tile_size, total - first);
        std::vector<CKKSNarrowBlindRotationRequest<Base, InputLogDelta>> requests;
        requests.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            const auto index = first + i;
            requests.push_back({&inputs[index / frequencies],
                plan.frequencies()[index % frequencies], {1.0, 0.0}});
        }
        std::vector<typename Schedule::Root> roots;
        CKKSNarrowBlindRotateBatch<Base, InputLogDelta>(roots, requests, blind_key, tile_size);
        for (std::size_t i = 0; i < count; ++i) {
            pending.push_back(std::move(roots[i]));
            if (pending.size() == frequencies) {
                CKKSNarrowEvalSeriesFromBasisRoots<Schedule>(results[(first + i) / frequencies], pending,
                    basis_key, plan, coefficients, use_fft);
                pending.clear();
            }
        }
        first += count;
    }
}

template <class Schedule, std::uint32_t InputLogDelta>
inline void CKKSNarrowEvalSeriesWithBasis(
    typename Schedule::Output &result,
    const CKKSCiphertext<typename Schedule::Param,
        Schedule::Base::input_log_q, InputLogDelta> &input,
    const CKKSNarrowBootstrapKey<typename Schedule::Base> &blind_key,
    const CKKSNarrowFrequencyBasisKey<Schedule> &basis_key,
    const CKKSNarrowFrequencyBasis &plan,
    const std::vector<std::complex<double>> &coefficients,
    std::size_t tile_size = 8, bool use_fft = true)
{
    std::vector<typename Schedule::Output> outputs;
    CKKSNarrowEvalSeriesWithBasisBatch<Schedule, InputLogDelta>(outputs, {&input, 1},
        blind_key, basis_key, plan, coefficients, tile_size, use_fft);
    result = std::move(outputs.front());
}

}  // namespace TFHEpp
#endif
