#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {
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

using namespace TFHEpp;

double maximum(double a, double b)
{
    if (!std::isfinite(b)) throw std::runtime_error("nonfinite diagnostic metric");
    return std::max(a, b);
}

template <std::uint32_t Nbit>
int run(const std::string &function, const std::filesystem::path &state,
        const std::filesystem::path &selectors, const std::filesystem::path &destination)
{
    using P = ProbeParam<Nbit>;
    using T = typename P::T;
    using Narrow = CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, true, 8, true>;
    using Refresh = CKKSNarrowRefreshSchedule<Narrow, 44, 37, 0, true>;
    std::ifstream file(state, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open retained auxiliary8 state");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t, 10> metadata{};
    std::string saved_function;
    bool finished = false;
    archive(metadata, saved_function, finished);
    const std::array<std::uint64_t, 10> expected{
        0x4155583850523031ULL, P::n, 85, 44, 385, 48, 37, 108, 8, 448};
    if (metadata != expected || saved_function != "sigmoid" || !finished)
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
    key->narrow.selector_source = std::make_shared<CKKSNarrowSelectorFile<Narrow>>(
        selectors, typename CKKSNarrowSelectorFile<Narrow>::Existing{});
    key->narrow.selector_ntt_chunk_size = std::min<std::size_t>(40, Narrow::block_width);
    key->narrow.bounded_tail_rotation_fft = true;
    key->narrow.staged_fft_caches = true;

    if (!std::filesystem::create_directory(destination))
        throw std::runtime_error("retained-function destination must be new");
    std::filesystem::permissions(destination, std::filesystem::perms::owner_all,
                                 std::filesystem::perm_options::replace);
    const auto target = [&](double x) {
        if (function == "sigmoid") return 1 / (1 + std::exp(-x));
        if (function == "gelu") return 0.5 * x * (1 + std::erf(x / std::sqrt(2.0)));
        return std::tanh(x);
    };
    const auto curvature = [&](double x) {
        const double y = target(x);
        if (function == "sigmoid") return y * (1-y) * (1-2*y);
        if (function == "gelu") return (2-x*x) * std::exp(-x*x/2) / std::sqrt(2*std::acos(-1.0));
        return -2*y*(1-y*y);
    };
    const auto approximation = function == "tanh" ?
        CKKSBuildNarrowWindowedTrigonometricApproximation(target,
            [](double x) { const double y = std::tanh(x); return 1-y*y; },
            curvature, {1, 1, 2, 4, 40}, 4, Refresh::gap(), 40) :
        CKKSBuildNarrowTrigonometricApproximation(target, 4, Refresh::gap(), 12,
            function == "gelu" ? 2.0 : 1.0, 65536, 40, 2.0, curvature, 160.0);
    auto decoded = std::make_unique<CKKSSlotVector<P>>();
    auto target_slots = std::make_unique<CKKSSlotVector<P>>();
    auto approximate_slots = std::make_unique<CKKSSlotVector<P>>();
    ckksSlotDecrypt<P>(*decoded, *input, *dense);
    double fresh = 0;
    for (std::size_t i=0; i<P::n/2; ++i) fresh=maximum(fresh,std::abs((*decoded)[i]-(*slots)[i]));
    auto output = std::make_unique<typename Refresh::Output>();
    std::cout << "evaluating retained " << function << std::endl;
    const auto start = std::chrono::steady_clock::now();
    CKKSNarrowFunctionalRefresh<Refresh>(*output, *input, *key, approximation, true, 12);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    // The function tag describes this output; input/keys remain those of the
    // original sigmoid archive. No archived-output identity is claimed.
    std::ofstream prefix(destination/"replayed-prefix.bin",std::ios::binary);
    { cereal::PortableBinaryOutputArchive archive(prefix);
      archive(metadata,function,finished,*dense,*sparse,*slots,*input,*output); }
    prefix.close(); if (!prefix) throw std::runtime_error("cannot save retained-function prefix");
    std::filesystem::permissions(destination/"replayed-prefix.bin",
        std::filesystem::perms::owner_read|std::filesystem::perms::owner_write,
        std::filesystem::perm_options::replace);
    ckksSlotDecrypt<P>(*decoded, *output, *dense);
    double error=0,evaluation=0,sum_squared=0,amplitude=0,approximate_amplitude=0;
    for (std::size_t i=0; i<P::n/2; ++i) {
        const auto x=(*slots)[i];
        (*target_slots)[i]={target(x.real()),target(x.imag())};
        (*approximate_slots)[i]={CKKSEvaluateNarrowApproximation(approximation,x.real()),
                              CKKSEvaluateNarrowApproximation(approximation,x.imag())};
        const double e=std::abs((*decoded)[i]-(*target_slots)[i]);
        error=maximum(error,e);sum_squared+=e*e;
        evaluation=maximum(evaluation,std::abs((*decoded)[i]-(*approximate_slots)[i]));
        amplitude=maximum(amplitude,std::abs((*target_slots)[i]));
        approximate_amplitude=maximum(approximate_amplitude,std::abs((*approximate_slots)[i]));
    }
    using Square=CKKSMultResult<P,108,48,108,48>;
    static_assert(Square::log_q==60 && Square::log_delta==48);
    auto square=std::make_unique<Square>();
    CKKSMultRelinThenRescale<P>(*square,*output,*output,*square_key);
    ckksSlotDecrypt<P>(*decoded,*square,*dense);
    double squared_error=0,squared_evaluation=0,square_sum=0;
    for (std::size_t i=0; i<P::n/2; ++i) {
        const double e=std::abs((*decoded)[i]-(*target_slots)[i]*(*target_slots)[i]);
        squared_error=maximum(squared_error,e);square_sum+=e*e;
        squared_evaluation=maximum(squared_evaluation,std::abs((*decoded)[i]-
            (*approximate_slots)[i]*(*approximate_slots)[i]));
    }
    if (!std::isfinite(sum_squared) || !std::isfinite(square_sum))
        throw std::runtime_error("nonfinite retained-function RMS");
    std::ofstream square_file(destination/"square.bin",std::ios::binary);
    { cereal::PortableBinaryOutputArchive archive(square_file);archive(*square); }
    square_file.close();if(!square_file)throw std::runtime_error("cannot save retained-function square");
    std::filesystem::permissions(destination/"square.bin",
        std::filesystem::perms::owner_read|std::filesystem::perms::owner_write,
        std::filesystem::perm_options::replace);
    const double tolerance=std::ldexp(1.0,-23),uniform=std::sqrt(2.0)*approximation.uniform_error_bound.value();
    const double rms=std::sqrt(sum_squared/(P::n/2));
    const double product_tolerance=2*amplitude*tolerance+tolerance*tolerance+tolerance;
    const double evaluation_product_tolerance=2*approximate_amplitude*tolerance+tolerance*tolerance+tolerance;
    const bool passed=std::max({fresh,error,evaluation,rms,uniform})<tolerance &&
        squared_error<product_tolerance && squared_evaluation<evaluation_product_tolerance;
    std::ofstream report(destination/"report.json");
    report<<std::setprecision(17)<<"{\"n\":"<<P::n<<",\"function\":\""<<function
          <<"\",\"auxiliary_bits\":8,\"fresh_keys_generated\":false,\"input_reused\":true,"
          <<"\"retained_keys_reused\":true,\"gaussian_sigma\":3.2,\"degree\":"<<(function=="tanh"?40:12)
          <<",\"series_tile_size\":12,\"fresh_input_max_error\":"<<fresh
          <<",\"max_error\":"<<error<<",\"rms_error\":"<<rms<<",\"evaluation_max_error\":"<<evaluation
          <<",\"complex_uniform_bound\":"<<uniform<<",\"target_amplitude\":"<<amplitude
          <<",\"post_refresh_product_max_error\":"<<squared_error
          <<",\"square_rms_error\":"<<std::sqrt(square_sum/(P::n/2))
          <<",\"post_refresh_product_evaluation_max_error\":"<<squared_evaluation
          <<",\"product_tolerance\":"<<product_tolerance
          <<",\"evaluation_product_tolerance\":"<<evaluation_product_tolerance
          <<",\"refresh_seconds\":"<<seconds<<",\"accuracy_passed\":"<<(passed?"true":"false")<<"}\n";
    report.close();if(!report)throw std::runtime_error("cannot save retained-function report");
    std::cout<<"retained "<<function<<" max error="<<error<<" square error="<<squared_error<<std::endl;
    return 0; // Diagnostic completion and accuracy are separate.
}
} // namespace

int main(int argc,char **argv)
{
    try {
        if(argc!=5 || (std::string(argv[1])!="sigmoid" && std::string(argv[1])!="gelu" && std::string(argv[1])!="tanh"))
            throw std::runtime_error("Usage: ckks_narrow_auxiliary_functions sigmoid|gelu|tanh STATE SELECTORS NEW_DIRECTORY");
        std::ifstream file(argv[2],std::ios::binary);
        if(!file)throw std::runtime_error("cannot open retained state");
        cereal::PortableBinaryInputArchive archive(file);
        std::array<std::uint64_t,10> metadata{};archive(metadata);
        if(metadata[0]!=0x4155583850523031ULL)throw std::runtime_error("wrong retained probe magic");
        if(metadata[1]==512)return run<9>(argv[1],argv[2],argv[3],argv[4]);
        if(metadata[1]==16384)return run<14>(argv[1],argv[2],argv[3],argv[4]);
        throw std::runtime_error("unsupported retained ring size");
    } catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 2;}
}
