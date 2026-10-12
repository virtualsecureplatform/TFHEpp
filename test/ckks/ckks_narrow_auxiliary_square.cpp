#include <ckks/compact_dd.hpp>
#include <cereal/archives/portable_binary.hpp>
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

template <std::uint32_t Nbit>
int run(const std::filesystem::path &state, const std::filesystem::path &prefix,
        const std::filesystem::path &destination)
{
    using P = ProbeParam<Nbit>;
    struct Prefix {
        std::array<std::uint64_t, 10> metadata{};
        std::string function;
        bool finished = false;
        Key<P> dense{}, sparse{};
        CKKSSlotVector<P> slots{};
        CKKSCiphertext<P,85,44> input{};
        CKKSCiphertext<P,108,48> output{};
        void read(cereal::PortableBinaryInputArchive &archive) {
            archive(metadata,function,finished,dense,sparse,slots,input,output);
            const std::array<std::uint64_t,10> expected{
                0x4155583850523031ULL,P::n,85,44,385,48,37,108,8,448};
            if (metadata != expected || function != "sigmoid" || !finished)
                throw std::runtime_error("unsupported or unfinished retained sigmoid prefix");
        }
    };
    std::ifstream original(state,std::ios::binary), replay(prefix,std::ios::binary);
    if (!original || !replay) throw std::runtime_error("cannot open retained square inputs");
    cereal::PortableBinaryInputArchive archive(original), replay_archive(replay);
    auto saved=std::make_unique<Prefix>(), current=std::make_unique<Prefix>();
    saved->read(archive); current->read(replay_archive);
    if (replay.peek()!=std::char_traits<char>::eof())
        throw std::runtime_error("unexpected replay prefix suffix");
    if (saved->metadata!=current->metadata || saved->dense!=current->dense ||
        saved->sparse!=current->sparse || saved->slots!=current->slots ||
        saved->input.ct!=current->input.ct)
        throw std::runtime_error("replay prefix does not match retained state");
    for (auto value: saved->dense)
        if (value!=typename P::T{} && value!=typename P::T{1} && value!=typename P::T{}-typename P::T{1})
            throw std::runtime_error("retained dense secret is not ternary");
    const auto prefix_end=original.tellg();
    auto square_key=std::make_unique<CKKSCompactDDRelinKey<P,108>>();
    const auto key_bytes=7*sizeof(std::uint64_t)+sizeof(square_key->data);
    original.seekg(0,std::ios::end);
    const auto end=original.tellg();
    if (end < prefix_end+static_cast<std::streamoff>(key_bytes))
        throw std::runtime_error("missing retained square key");
    original.seekg(end-static_cast<std::streamoff>(key_bytes));
    archive(*square_key); // Preserve the original archive byte order at the suffix.
    if (original.tellg()!=end || original.peek()!=std::char_traits<char>::eof())
        throw std::runtime_error("invalid retained square archive termination");
    if (!std::filesystem::create_directory(destination))
        throw std::runtime_error("square diagnostic destination must be new");
    std::filesystem::permissions(destination,std::filesystem::perms::owner_all,
                                std::filesystem::perm_options::replace);
    using Square=CKKSMultResult<P,108,48,108,48>;
    static_assert(Square::log_q==60 && Square::log_delta==48);
    auto square=std::make_unique<Square>();
    CKKSMultRelinThenRescale<P>(*square,current->output,current->output,*square_key);
    auto decoded=std::make_unique<CKKSSlotVector<P>>();
    ckksSlotDecrypt<P>(*decoded,*square,saved->dense);
    double maximum=0,sum=0,amplitude=0;
    for (std::size_t i=0;i<P::n/2;++i) {
        const auto x=saved->slots[i];
        if (!std::isfinite(x.real()) || !std::isfinite(x.imag()) ||
            std::abs(x.real())>4 || std::abs(x.imag())>4)
            throw std::runtime_error("invalid retained slot domain");
        const auto sigmoid=[](double x){return 1/(1+std::exp(-x));};
        const std::complex<double> target{sigmoid(x.real()),sigmoid(x.imag())};
        const double error=std::abs((*decoded)[i]-target*target);
        if (!std::isfinite(error)) throw std::runtime_error("nonfinite square metric");
        maximum=std::max(maximum,error);sum+=error*error;amplitude=std::max(amplitude,std::abs(target));
    }
    if (!std::isfinite(sum)) throw std::runtime_error("nonfinite square RMS");
    const double tolerance=std::ldexp(1.0,-23);
    const double bound=2*amplitude*tolerance+tolerance*tolerance+tolerance;
    std::ofstream ciphertext(destination/"square.bin",std::ios::binary);
    {cereal::PortableBinaryOutputArchive output(ciphertext);output(*square);}
    ciphertext.close();if(!ciphertext)throw std::runtime_error("cannot save square ciphertext");
    std::filesystem::permissions(destination/"square.bin",std::filesystem::perms::owner_read|
        std::filesystem::perms::owner_write,std::filesystem::perm_options::replace);
    std::ofstream report(destination/"report.json");
    report<<std::setprecision(17)<<"{\"n\":"<<P::n<<",\"function\":\"sigmoid\","
          <<"\"fresh_keys_generated\":false,\"log_q\":60,\"log_delta\":48,"
          <<"\"retained_square_evaluated\":true,\"independent_phase_validated\":false,"
          <<"\"post_refresh_product_max_error\":"<<maximum<<",\"rms_error\":"<<std::sqrt(sum/(P::n/2))
          <<",\"target_amplitude\":"<<amplitude<<",\"product_tolerance\":"<<bound
          <<",\"product_accuracy_passed\":"<<(maximum<bound?"true":"false")<<"}\n";
    report.close();if(!report)throw std::runtime_error("cannot save square report");
    std::cout<<"retained square max error="<<maximum<<std::endl;
    return 0; // Diagnostic success is separate from accuracy.
}
}

int main(int argc,char **argv)
{
    try {
        if(argc!=4)throw std::runtime_error("Usage: ckks_narrow_auxiliary_square STATE REPLAY_PREFIX NEW_DIRECTORY");
        std::ifstream file(argv[1],std::ios::binary);
        if(!file)throw std::runtime_error("cannot open retained state");
        cereal::PortableBinaryInputArchive archive(file);
        std::array<std::uint64_t,10> metadata{};archive(metadata);
        if(metadata[0]!=0x4155583850523031ULL)throw std::runtime_error("wrong retained probe magic");
        if(metadata[1]==512)return run<9>(argv[1],argv[2],argv[3]);
        if(metadata[1]==16384)return run<14>(argv[1],argv[2],argv[3]);
        throw std::runtime_error("unsupported retained ring size");
    } catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 2;}
}
