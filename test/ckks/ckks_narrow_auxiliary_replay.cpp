#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>

namespace {
using namespace TFHEpp;
template <class Schedule, class Input>
void checkSharedTailCache(std::ostream &report, const Input &input,
    const CKKSNarrowBootstrapKey<Schedule> &key)
{
    using P = typename Schedule::Param;
    using Clock = std::chrono::steady_clock;
    std::vector<typename Schedule::Selector> selectors;
    key.selector_source->load(true, 0, selectors);
    std::array<double, P::n> a{};
    for (std::size_t i = 0; i < P::n; ++i)
        a[i] = -static_cast<double>(std::ldexp(ckks_detail::levelToLongDouble<
            P, Schedule::input_log_q>(input.ct[0][i]), -static_cast<int>(Schedule::input_log_q)));
    const auto table = narrow_detail::tailTable<Schedule>(a, 0, 1);
    auto conjugation = makeCKKSPreparedAutoKeyFFT<P, Schedule::factor_key_log_q>(key.factor_galois[P::nbit]);
    auto factors = narrow_detail::packedFactorsBatch<Schedule>(
        {&table, 1}, selectors, {}, 0, key, 40, conjugation.get());
    conjugation.reset();
    selectors.clear();
    selectors.shrink_to_fit();
    constexpr std::size_t columns = 3;
    std::vector<typename Schedule::Factor> reference(columns * Schedule::blocks);
    for (std::size_t column = 0; column < columns; ++column)
        for (std::size_t x = 0; x < Schedule::blocks; ++x) {
            auto &leaf = reference[column * Schedule::blocks + x];
            leaf = factors.front();
            // Distinct public column multiples expose provider-index mistakes.
            for (auto &component : leaf.ct)
                for (auto &value : component)
                    value = ckks_detail::reduceToLevel<P, Schedule::factor_log_q>(value << column);
        }
    auto actual = reference;
    std::size_t old_peak = 0;
    auto start = Clock::now();
    for (std::size_t q = 0; q < Schedule::tail_packs; ++q) {
        const auto count = std::min(Schedule::tail_blocks_per_pack,
            Schedule::blocks - q * Schedule::tail_blocks_per_pack);
        std::vector<int> steps(columns * count);
        for (std::size_t column = 0; column < columns; ++column)
            for (std::size_t ell = 0; ell < count; ++ell)
                steps[column * count + ell] = -static_cast<int>(Schedule::slots -
                    2 * (ell + 1) * Schedule::tail_width);
        old_peak = std::max(old_peak,
            CKKSRotateSlotsBatchWithBoundedFFT<P, Schedule::factor_log_q, Schedule::auxiliary_log_scale>(
                steps, [&](std::size_t i) -> auto & {
                    return reference[(i / count) * Schedule::blocks +
                        q * Schedule::tail_blocks_per_pack + i % count].ct;
                }, key.factor_galois));
    }
    const double old_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    std::vector<int> steps(columns * Schedule::blocks);
    for (std::size_t column = 0; column < columns; ++column)
        for (std::size_t x = 0; x < Schedule::blocks; ++x)
            steps[column * Schedule::blocks + x] = -static_cast<int>(Schedule::slots -
                2 * (x % Schedule::tail_blocks_per_pack + 1) * Schedule::tail_width);
    start = Clock::now();
    const auto new_peak = CKKSRotateSlotsBatchWithBoundedFFT<P, Schedule::factor_log_q, Schedule::auxiliary_log_scale>(
        steps, [&](std::size_t i) -> auto & { return actual[i].ct; }, key.factor_galois);
    const double new_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    for (std::size_t i = 0; i < actual.size(); ++i)
        if (actual[i].ct != reference[i].ct)
            throw std::runtime_error("shared tail rotation cache changed ciphertext coefficients");
    if (old_peak != new_peak) throw std::runtime_error("shared tail rotation cache changed peak key payload");
    report << "\"tail_columns\":" << columns << ",\"tail_leaves\":" << actual.size()
           << ",\"tail_coefficients_equal\":true,\"peak_rotation_key_bytes\":" << new_peak
           << ",\"per_pack_alignment_seconds\":" << old_seconds
           << ",\"shared_alignment_seconds\":" << new_seconds
           << ",\"whole_refresh_validated\":false";
    std::cout << "retained tail alignment: all " << actual.size() << " ciphertexts equal" << std::endl;
}

template <std::uint32_t Nbit>
struct ProbeParam {
    static constexpr std::int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = Nbit, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t l = 1, l_a = 1, lₐ = 1;
    static constexpr std::uint32_t Bgbit = 4, Bg_abit = 4, Bgₐbit = 4;
    using T = TFHEpp::MultiLimbUInt<7>;
    static constexpr T Bg = T{1} << Bgbit, Bg_a = Bg, Bgₐ = Bg;
    static constexpr TFHEpp::ErrorDistribution errordist = TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::ldexp(1.0, -448);
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << 253;
    static constexpr std::uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 112, l̅ₐ = l̅, B̅gbit = 4, B̅gₐbit = 4;
};

double maximum(double a, double b)
{
    if (!std::isfinite(b)) throw std::runtime_error("nonfinite diagnostic metric");
    return std::max(a, b);
}

// Diagnostic only: consume one complete selector pass per component. With
// degree12 and tile12, all real roots precede all imaginary roots. Refuse any
// other load order rather than silently mixing a different evaluation layout.
template <class Schedule>
class ComponentSelectorSource final : public CKKSNarrowSelectorSource<Schedule> {
    std::array<std::shared_ptr<CKKSNarrowSelectorSource<Schedule>>, 2> sources_;
    mutable std::size_t loads_ = 0;
    static constexpr auto groups = Schedule::blocks + Schedule::tail_packs;
public:
    ComponentSelectorSource(std::shared_ptr<CKKSNarrowSelectorSource<Schedule>> first,
                            std::shared_ptr<CKKSNarrowSelectorSource<Schedule>> second)
        : sources_{std::move(first), std::move(second)}
    {
        for (bool tail : {false,true})
            if (sources_[0]->size(tail) != sources_[1]->size(tail))
                throw std::runtime_error("component selector source sizes differ");
    }
    std::size_t size(bool tail) const override { return sources_[0]->size(tail); }
    void load(bool tail, std::size_t group,
              std::vector<typename Schedule::Selector> &result) const override
    {
        const auto position = (tail ? Schedule::blocks : 0) + group;
        if (loads_ >= 2*groups || position != loads_%groups)
            throw std::runtime_error("unexpected component selector load order");
        sources_[loads_/groups]->load(tail,group,result);
        ++loads_;
    }
    void verifyComplete() const
    {
        if (loads_ != 2*groups) throw std::runtime_error("incomplete component selector passes");
    }
    std::size_t loads() const { return loads_; }
};

// Compare derived caches from the same retained encrypted primary rows.
// Inputs are public arithmetic probes, not newly encrypted ciphertexts.
template <class P, std::uint32_t Q, std::uint32_t Bg, std::uint32_t Bbar>
void checkSecondaryDDCache(std::ostream &report,
    const CKKSCompactDDRelinKey<P, Q, Bg, Bbar> &source)
{
    using Clock = std::chrono::steady_clock;
    auto elapsed = [](auto start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    auto start = Clock::now();
    auto old_fft = makeCKKSDDRelinKeyFFT(source);
    const double old_prepare = elapsed(start);
    start = Clock::now();
    auto converted = makeCKKSCompactDDRelinKeyWithSecondaryBase<16>(source);
    const double conversion = elapsed(start);
    auto restored = makeCKKSCompactDDRelinKeyWithSecondaryBase<Bbar>(*converted);
    if (restored->data != source.data)
        throw std::runtime_error("retained DD primary rows changed during secondary conversion");
    restored.reset();
    start = Clock::now();
    auto new_fft = makeCKKSDDRelinKeyFFT(*converted);
    const double new_prepare = elapsed(start);
    converted.reset();
    auto polynomial = std::make_unique<Polynomial<P>>();
    auto reference = std::make_unique<TRLWE<P>>();
    auto actual = std::make_unique<TRLWE<P>>();
    std::mt19937_64 engine(20262119 + Q);
    double old_switch = 0, new_switch = 0;
    for (std::size_t sample = 0; sample < 3; ++sample) {
        for (std::size_t i = 0; i < P::n; ++i) {
            auto &value = (*polynomial)[i];
            if (sample == 0) {
                for (auto &limb : value.limb) limb = engine();
                value = ckks_detail::reduceToLevel<P, Q>(value);
            }
            else if (sample == 1)
                value = (typename P::T{1} << (Q - 1)) - typename P::T{i % 2};
            else value = i == P::n / 2 ? typename P::T{} - typename P::T{1} : typename P::T{};
        }
        start = Clock::now();
        CKKSRelinKeySwitch<P>(*reference, *polynomial, *old_fft);
        old_switch += elapsed(start);
        start = Clock::now();
        CKKSRelinKeySwitch<P>(*actual, *polynomial, *new_fft);
        new_switch += elapsed(start);
        if (*reference != *actual)
            throw std::runtime_error("retained DD secondary FFT caches changed key-switch coefficients");
    }
    report << "{\"log_q\":" << Q << ",\"samples\":3,\"round_trip_equal\":true,"
           << "\"switch_coefficients_equal\":true,\"old_cache_bytes\":" << sizeof(*old_fft)
           << ",\"new_cache_bytes\":" << sizeof(*new_fft)
           << ",\"old_prepare_seconds\":" << old_prepare
           << ",\"conversion_seconds\":" << conversion
           << ",\"new_prepare_seconds\":" << new_prepare
           << ",\"old_switch_seconds\":" << old_switch
           << ",\"new_switch_seconds\":" << new_switch << '}';
    std::cout << "retained DD Q" << Q << " secondary cache round trip and 3 switches equal" << std::endl;
}

template <class P, std::uint32_t Q>
void checkSecondaryAutoCache(std::ostream &report, const TRLWE<P> &retained,
    const CKKSAutoKey<P, Q> &raw)
{
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    auto old_fft = makeCKKSPreparedAutoKeyFFT<P, Q, P::B̅gbit>(raw);
    const auto old_prepare = std::chrono::duration<double>(Clock::now() - start).count();
    start = Clock::now();
    auto new_fft = makeCKKSPreparedAutoKeyFFT<P, Q, 16>(raw);
    const auto new_prepare = std::chrono::duration<double>(Clock::now() - start).count();
    auto input = std::make_unique<TRLWE<P>>();
    auto reference = std::make_unique<TRLWE<P>>(), actual = std::make_unique<TRLWE<P>>();
    std::mt19937_64 engine(20262119 + Q);
    double old_seconds = 0, new_seconds = 0;
    static_assert(Q >= 85);
    for (std::size_t sample = 0; sample < 3; ++sample) {
        for (std::size_t c = 0; c <= P::k; ++c)
            for (std::size_t i = 0; i < P::n; ++i) {
                auto &value = (*input)[c][i];
                if (sample == 0)
                    value = ckks_detail::reduceToLevel<P, 85>(retained[c][i]) << (Q - 85);
                else if (sample == 1) {
                    for (auto &limb : value.limb) limb = engine();
                    value = ckks_detail::reduceToLevel<P, Q>(value);
                }
                else value = (typename P::T{1} << (Q - 1)) - typename P::T{(i + c) % 2};
            }
        start = Clock::now();
        CKKSEvalAuto<P, Q>(*reference, *input, 2 * P::n - 1, *old_fft);
        old_seconds += std::chrono::duration<double>(Clock::now() - start).count();
        start = Clock::now();
        CKKSEvalAuto<P, Q>(*actual, *input, 2 * P::n - 1, *new_fft);
        new_seconds += std::chrono::duration<double>(Clock::now() - start).count();
        if (*reference != *actual) throw std::runtime_error("secondary auto cache changed full-ring coefficients");
    }
    report << "{\"log_q\":" << Q << ",\"samples\":3,\"coefficients_equal\":true,"
           << "\"old_cache_bytes\":" << old_fft->spectra.size() * sizeof(PolynomialInFD<P>)
           << ",\"new_cache_bytes\":" << new_fft->spectra.size() * sizeof(PolynomialInFD<P>)
           << ",\"old_prepare_seconds\":" << old_prepare << ",\"new_prepare_seconds\":" << new_prepare
           << ",\"old_switch_seconds\":" << old_seconds << ",\"new_switch_seconds\":" << new_seconds << '}';
    std::cout << "retained auto Q" << Q << " three secondary-cache switches equal" << std::endl;
}

template <std::uint32_t Nbit>
int run(const std::string &mode, const std::filesystem::path &state,
        const std::filesystem::path &selectors, const std::filesystem::path &destination,
        const std::filesystem::path &second_selectors = {})
{
    using P = ProbeParam<Nbit>;
    using T = typename P::T;
    using Narrow = CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, true, 8, true>;
    using Refresh = CKKSNarrowRefreshSchedule<Narrow, 44, 37, 0, true>;
    std::ifstream file(state, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open retained auxiliary8 state");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t, 10> metadata{};
    std::string function;
    bool finished = false;
    archive(metadata, function, finished);
    const std::array<std::uint64_t, 10> expected{
        0x4155583850523031ULL, P::n, 85, 44, 385, 48, 37, 108, 8, 448};
    if (metadata != expected || function != "sigmoid" || !finished)
        throw std::runtime_error("replay requires a finished auxiliary8 sigmoid probe state");
    auto dense = std::make_unique<Key<P>>(), sparse = std::make_unique<Key<P>>();
    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto input = std::make_unique<typename Refresh::Input>();
    auto saved = std::make_unique<typename Refresh::Output>();
    auto key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    auto square_key = std::make_unique<CKKSCompactDDRelinKey<P, 108>>();
    std::cout << "loading retained N" << P::n << " keys and ciphertexts" << std::endl;
    archive(*dense, *sparse, *slots, *input, *saved, *key, *square_key);
    if (file.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("unexpected retained auxiliary8 archive suffix");
    file.close();
    narrow_detail::structuredPositions<Narrow>(*sparse);
    for (const auto &value : *dense)
        if (value != T{} && value != T{1} && value != T{}-T{1})
            throw std::runtime_error("retained dense secret is not ternary");
    for (auto value : *slots)
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()) ||
            std::abs(value.real()) > 4 || std::abs(value.imag()) > 4)
            throw std::runtime_error("invalid retained slot domain");
    auto first_source = std::make_shared<CKKSNarrowSelectorFile<Narrow>>(
        selectors, typename CKKSNarrowSelectorFile<Narrow>::Existing{});
    std::shared_ptr<ComponentSelectorSource<Narrow>> component_source;
    if (mode == "--selector-halves") {
        auto second_source = std::make_shared<CKKSNarrowSelectorFile<Narrow>>(
            second_selectors, typename CKKSNarrowSelectorFile<Narrow>::Existing{});
        component_source = std::make_shared<ComponentSelectorSource<Narrow>>(first_source,second_source);
        key->narrow.selector_source = component_source;
    }
    else key->narrow.selector_source = first_source;
    key->narrow.selector_ntt_chunk_size = std::min<std::size_t>(40, Narrow::block_width);
    key->narrow.bounded_tail_rotation_fft = true;
    key->narrow.staged_fft_caches = true;
    if (!std::filesystem::create_directory(destination))
        throw std::runtime_error("diagnostic destination must be new");
    std::filesystem::permissions(destination, std::filesystem::perms::owner_all,
        std::filesystem::perm_options::replace);
    std::ofstream report(destination/"report.json");
    if (!report) throw std::runtime_error("cannot create diagnostic report");
    report << std::setprecision(17) << "{\"n\":" << P::n << ",\"function\":\"sigmoid\",\"auxiliary_bits\":8,"
           << "\"mode\":\"" << mode << "\",\"fresh_keys_generated\":false,";
    if (mode == "--auto-cache-control") {
        report << "\"auto_levels\":[";
        checkSecondaryAutoCache<P, Narrow::factor_key_log_q>(report, input->ct, key->narrow.factor_galois[P::nbit]);
        report << ',';
        checkSecondaryAutoCache<P, Narrow::output_key_log_q>(report, input->ct, key->narrow.output_conjugation);
        report << "],\"all_auto_controls_passed\":true,\"whole_refresh_validated\":false}\n";
        report.close();
        if (!report) throw std::runtime_error("cannot finish auto cache control report");
        return 0;
    }
    if (mode == "--tail-cache-control") {
        checkSharedTailCache<Narrow>(report, *input, key->narrow);
        report << "}\n";
        report.close();
        if (!report) throw std::runtime_error("cannot finish tail cache control report");
        return 0;
    }
    if (mode == "--dd-cache-control") {
        report << "\"dd_levels\":[";
        checkSecondaryDDCache<P>(report, key->narrow.relin.template get<0>());
        report << ',';
        checkSecondaryDDCache<P>(report, key->narrow.relin.template get<Narrow::tree_depth - 1>());
        report << "],\"all_dd_controls_passed\":true,\"whole_refresh_validated\":false}\n";
        report.close();
        if (!report) throw std::runtime_error("cannot finish DD cache control report");
        return 0;
    }
    auto decoded = std::make_unique<CKKSSlotVector<P>>();
    ckksSlotDecrypt<P>(*decoded, *input, *dense);
    double fresh = 0;
    for (std::size_t i=0; i<P::n/2; ++i) fresh=maximum(fresh,std::abs((*decoded)[i]-(*slots)[i]));
    report << "\"fresh_input_max_error\":" << fresh << ',';
    auto target=[](double x) { return 1/(1+std::exp(-x)); };
    auto measure=[&](const auto &ciphertext) {
        ckksSlotDecrypt<P>(*decoded,ciphertext,*dense);
        double error=0, squared=0;
        for (std::size_t i=0; i<P::n/2; ++i) {
            const auto x=(*slots)[i];
            const auto e=std::abs((*decoded)[i]-std::complex<double>{target(x.real()),target(x.imag())});
            error=maximum(error,e); squared+=e*e;
        }
        if (!std::isfinite(squared)) throw std::runtime_error("nonfinite diagnostic RMS");
        return std::pair{error,std::sqrt(squared/(P::n/2))};
    };
    const auto saved_error=measure(*saved);
    report << "\"saved_max_error\":" << saved_error.first << ",\"saved_rms_error\":" << saved_error.second << ',';
    std::cout << "saved output max error=" << saved_error.first << std::endl;
    if (mode=="--prepare") {
        auto lifted=std::make_unique<typename Refresh::EvaluationInput>();
        for (std::size_t c=0; c<=P::k; ++c)
            for (std::size_t i=0; i<P::n; ++i)
                lifted->ct[c][i]=ckks_detail::reduceToLevel<P,85>(input->ct[c][i]) << 8;
        auto dense_coeff=std::make_unique<typename Refresh::EvaluationCoefficients>();
        auto sparse_coeff=std::make_unique<typename Refresh::EvaluationCoefficients>();
        auto prepared=std::make_unique<typename Refresh::Coefficients>();
        auto phase=std::make_unique<std::array<double,P::n>>();
        auto coefficient_error=[&](const auto &ciphertext,const auto &secret) {
            ckksDecrypt<P>(*phase,ciphertext,secret);
            double error=0;
            for (std::size_t i=0; i<P::n/2; ++i) {
                error=maximum(error,std::abs((*phase)[i]-(*slots)[i].real()));
                error=maximum(error,std::abs((*phase)[i+P::n/2]-(*slots)[i].imag()));
            }
            return error;
        };
        std::cout << "evaluating retained S2C" << std::endl;
        CKKSRealLinearTransform<P,93,52,37>(*dense_coeff,*lifted,key->s2c_plan,
            key->s2c_input_galois,key->s2c_output_galois);
        const auto s2c=coefficient_error(*dense_coeff,*dense);
        std::cout << "S2C max error=" << s2c << std::endl;
        CKKSSecretKeySwitch<P>(*sparse_coeff,*dense_coeff,key->to_sparse);
        const auto encapsulation=coefficient_error(*sparse_coeff,*sparse);
        CKKSScaleDown<P,56,52,8>(*prepared,*sparse_coeff);
        const auto final=coefficient_error(*prepared,*sparse);
        std::cout << "encapsulation max error=" << encapsulation << " prepared max error=" << final << std::endl;
        std::ofstream output(destination/"prepared.bin",std::ios::binary);
        { cereal::PortableBinaryOutputArchive saved_archive(output); saved_archive(*prepared); }
        output.close(); if (!output) throw std::runtime_error("cannot write prepared diagnostic ciphertext");
        report << "\"s2c_max_error\":" << s2c << ",\"encapsulation_max_error\":" << encapsulation
               << ",\"prepared_max_error\":" << final << ",\"prepare_accuracy_passed\":"
               << (std::max({fresh,s2c,encapsulation,final})<std::ldexp(1.0,-23)?"true":"false");
    }
    else {
        const bool persistent_fft = mode == "--replay-prepared";
        const bool projected_combined = mode == "--replay-combined-primary8";
        const bool projected_primary = mode == "--replay-primary8" || projected_combined;
        const bool projected_dd = mode == "--replay-dd-primary8" || projected_combined;
        const auto cache_start = std::chrono::steady_clock::now();
        if (persistent_fft) {
            CKKSNarrowPrepareAutomorphisms<Narrow>(key->narrow);
            CKKSNarrowPrepareRelinearization<Narrow>(key->narrow);
        }
        if (projected_primary) CKKSNarrowPrepareProjectedAutomorphisms<Narrow>(key->narrow);
        if (projected_dd) CKKSNarrowPrepareProjectedRelinearization<Narrow>(key->narrow);
        report << "\"projected_primary8_conjugation\":" << (projected_primary ? "true" : "false") << ","
               << "\"projected_primary8_dd\":" << (projected_dd ? "true" : "false") << ","
               << "\"persistent_fft_caches\":" << (persistent_fft ? "true" : "false")
               << ",\"derived_cache_prepare_seconds\":" << std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - cache_start).count() << ',';
        auto curvature=[&](double x) { const double y=target(x); return y*(1-y)*(1-2*y); };
        const auto approximation=CKKSBuildNarrowTrigonometricApproximation(
            target,4,Refresh::gap(),12,1.0,65536,40,2.0,curvature,160.0);
        CKKSNarrowRootObserver<Narrow,44> observer;
        std::ofstream trace;
        if (mode=="--trace") {
            trace.open(destination/"roots.jsonl");
            if (!trace) throw std::runtime_error("cannot create root trace");
            trace << std::setprecision(17);
            observer=[&](const auto &request,const auto &root) {
                auto phase=std::make_unique<std::array<double,P::n>>();
                auto got=std::make_unique<CKKSSlotVector<P>>();
                ckksDecrypt<P>(*phase,*request.input,*sparse);
                ckksSlotDecrypt<P>(*got,root,*dense);
                double error=0,amplitude=0;
                for (std::size_t i=0;i<P::n/2;++i) {
                    const auto expected_root=request.weight*std::polar(1.0,
                        2*std::acos(-1.0)*request.frequency*(*phase)[i]/Refresh::gap());
                    error=maximum(error,std::abs((*got)[i]-expected_root));
                    amplitude=maximum(amplitude,std::abs((*got)[i]));
                }
                trace << "{\"frequency\":" << request.frequency << ",\"weight_abs\":"
                      << std::abs(request.weight) << ",\"max_error\":" << error
                      << ",\"max_abs\":" << amplitude << "}\n";
                trace.flush();if (!trace) throw std::runtime_error("cannot write root trace");
                std::cout << "root frequency=" << request.frequency << " max error=" << error << std::endl;
            };
        }
        auto output=std::make_unique<typename Refresh::Output>();
        CKKSNarrowFactorTimings factor_totals;
        std::size_t completed_groups = 0;
        key->narrow.factor_observer = [&](bool tail, std::size_t group, const auto &timings) {
            ++completed_groups;
            factor_totals.selector_load_seconds += timings.selector_load_seconds;
            factor_totals.table_ntt_seconds += timings.table_ntt_seconds;
            factor_totals.selector_ntt_seconds += timings.selector_ntt_seconds;
            factor_totals.matrix_seconds += timings.matrix_seconds;
            factor_totals.inverse_crt_seconds += timings.inverse_crt_seconds;
            factor_totals.conjugation_seconds += timings.conjugation_seconds;
            std::cout << "completed factor group " << completed_groups << " "
                      << (tail ? "tail" : "head") << " " << group << std::endl;
        };
        std::cout << "replaying retained functional refresh" << std::endl;
        const auto refresh_start = std::chrono::steady_clock::now();
        CKKSNarrowFunctionalRefresh<Refresh>(*output,*input,*key,approximation,true,12,
            mode=="--trace"?&observer:nullptr);
        const double refresh_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - refresh_start).count();
        key->narrow.factor_observer = {};
        report << "\"refresh_seconds\":" << refresh_seconds
               << ",\"selector_ntt_digit_bits\":" << CKKSNarrowNTTSelector<Narrow>::digit_bits
               << ",\"factor_groups_completed\":" << completed_groups
               << ",\"factor_stage_seconds\":{\"selector_load\":" << factor_totals.selector_load_seconds
               << ",\"table_ntt\":" << factor_totals.table_ntt_seconds
               << ",\"selector_ntt\":" << factor_totals.selector_ntt_seconds
               << ",\"matrix\":" << factor_totals.matrix_seconds
               << ",\"inverse_crt\":" << factor_totals.inverse_crt_seconds
               << ",\"conjugation\":" << factor_totals.conjugation_seconds << "},";
        if (component_source) component_source->verifyComplete();
        const auto error=measure(*output);
        const bool equal=output->ct==saved->ct;
        std::ofstream prefix(destination/"replayed-prefix.bin",std::ios::binary);
        { cereal::PortableBinaryOutputArchive saved_archive(prefix);
          saved_archive(metadata,function,finished,*dense,*sparse,*slots,*input,*output); }
        prefix.close();if (!prefix) throw std::runtime_error("cannot save replayed diagnostic prefix");
        report << "\"max_error\":" << error.first << ",\"rms_error\":" << error.second
               << ",\"saved_output_equal\":" << (equal?"true":"false") << ",\"accuracy_passed\":"
               << (std::max({fresh,error.first,error.second})<std::ldexp(1.0,-23)?"true":"false")
               << ",\"retained_square_validated\":false"
               << ",\"selector_halves_mixed\":" << (component_source ? "true" : "false")
               << ",\"selector_group_loads\":" << (component_source ? component_source->loads() : 0);
        std::cout << "replay max error=" << error.first << " saved output equal=" << equal << std::endl;
    }
    report << "}\n";report.close();if (!report) throw std::runtime_error("cannot finish diagnostic report");
    return 0; // Successful diagnosis does not imply the ciphertext passes accuracy gates.
}
} // namespace

int main(int argc,char **argv)
{
    try {
        const bool component_mode = argc==6 && std::string(argv[1])=="--selector-halves";
        if (!component_mode && (argc!=5 || (std::string(argv[1])!="--prepare" && std::string(argv[1])!="--replay" && std::string(argv[1])!="--replay-prepared" && std::string(argv[1])!="--replay-primary8" && std::string(argv[1])!="--replay-dd-primary8" && std::string(argv[1])!="--replay-combined-primary8" && std::string(argv[1])!="--trace" && std::string(argv[1])!="--dd-cache-control" && std::string(argv[1])!="--tail-cache-control" && std::string(argv[1])!="--auto-cache-control"))) {
            std::cerr << "Usage: ckks_narrow_auxiliary_replay --prepare|--replay|--replay-prepared|--replay-primary8|--replay-dd-primary8|--replay-combined-primary8|--trace|--dd-cache-control|--tail-cache-control|--auto-cache-control STATE SELECTORS NEW_DIRECTORY\n";
            std::cerr << "       ckks_narrow_auxiliary_replay --selector-halves STATE REAL_SELECTORS IMAGINARY_SELECTORS NEW_DIRECTORY\n";
            return 2;
        }
        const std::filesystem::path destination = component_mode ? argv[5] : argv[4];
        const std::filesystem::path second = component_mode ? argv[4] : "";
        std::ifstream file(argv[2],std::ios::binary);
        if (!file) throw std::runtime_error("cannot open retained state");
        cereal::PortableBinaryInputArchive archive(file);
        std::array<std::uint64_t,10> metadata{};archive(metadata);
        if (metadata[0]!=0x4155583850523031ULL) throw std::runtime_error("wrong retained probe magic");
        if (metadata[1]==512) return run<9>(argv[1],argv[2],argv[3],destination,second);
        if (metadata[1]==16384) return run<14>(argv[1],argv[2],argv[3],destination,second);
        throw std::runtime_error("unsupported retained probe ring");
    }
    catch(const std::exception &error) { std::cerr << error.what() << '\n';return 2; }
}
