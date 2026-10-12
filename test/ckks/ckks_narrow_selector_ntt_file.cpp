#include "ckks_narrow_basis_benchmark.hpp"

template <std::uint32_t Nbit, bool Selector56 = false>
int prepare(const std::filesystem::path &checkpoint, const std::filesystem::path &directory)
{
#if defined(USE_BLAKE3)
    using P = narrow_benchmark::LeanBasisParam<Nbit>;
    using Base = TFHEpp::CKKSNarrowBootstrapSchedule<P,48,518,47,39,31,true,44,true>;
    using S = std::conditional_t<Selector56, narrow_benchmark::Selector56BasisSchedule<Base>, Base>;
    using Raw = TFHEpp::CKKSNarrowSelectorFile<S>;
    using Cache = TFHEpp::CKKSNarrowNTTSelectorFile<S>;
    auto raw = std::make_shared<Raw>(checkpoint/"selectors", typename Raw::Existing{});
    auto prepared = Cache::prepare(directory, raw, 5);
    // Reopen the complete cache and compare differently sized partial chunks
    // against freshly prepared spectra from the same retained raw records.
    Cache reopened(directory, raw, typename Cache::Existing{});
    std::size_t words = 0;
    for (bool tail : {false,true})
        for (std::size_t group = 0; group < (tail ? S::tail_packs : S::blocks); ++group)
            for (std::size_t first = 0; first < S::block_width;) {
                const auto count = std::min<std::size_t>(7,S::block_width-first);
                std::vector<typename S::Selector> ciphertexts;
                raw->loadChunk(tail,group,first,count,ciphertexts);
                std::vector<TFHEpp::CKKSNarrowNTTSelector<S>> expected, actual;
                TFHEpp::narrow_detail::prepareSelectors<S>(expected,ciphertexts);
                reopened.loadPreparedChunk(tail,group,first,count,actual);
                if (actual.size()!=expected.size()) throw std::runtime_error("prepared chunk count changed");
                for (std::size_t i=0;i<count;++i) {
                    if(actual[i].spectra!=expected[i].spectra) throw std::runtime_error("prepared cache spectra changed");
                    for(const auto &spectrum:actual[i].spectra) words+=spectrum.size();
                }
                first+=count;
            }
    std::cout << "{\"n\":" << P::n << ",\"exact_spectrum_words\":" << words
              << ",\"matrix_digit_bits\":" << TFHEpp::CKKSNarrowNTTSelector<S>::digit_bits
              << ",\"selector_digits\":" << TFHEpp::CKKSNarrowNTTSelector<S>::digits
              << ",\"fresh_keys_generated\":false,\"passed\":true}" << std::endl;
    return 0;
#else
    throw std::runtime_error("prepared selector files require BLAKE3 support");
#endif
}

template <std::uint32_t Nbit, bool Selector56 = false>
int sampleChunks(const std::filesystem::path &checkpoint, const std::filesystem::path &directory)
{
#if defined(USE_BLAKE3)
    using P = narrow_benchmark::LeanBasisParam<Nbit>;
    using Base = TFHEpp::CKKSNarrowBootstrapSchedule<P,48,518,47,39,31,true,44,true>;
    using S = std::conditional_t<Selector56, narrow_benchmark::Selector56BasisSchedule<Base>, Base>;
    using Raw = TFHEpp::CKKSNarrowSelectorFile<S>;
    using File = TFHEpp::CKKSNarrowNTTSelectorFile<S>;
    using Cache = TFHEpp::CKKSNarrowNTTSelector<S>;
    auto raw = std::make_shared<Raw>(checkpoint/"selectors",typename Raw::Existing{});
    if (!std::filesystem::create_directory(directory)) throw std::runtime_error("sample directory must be new");
    std::filesystem::permissions(directory,std::filesystem::perms::owner_all);
    std::size_t trial = 0, total_words = 0, rejected = 0;
    for (bool tail : {false,true})
        for (std::size_t position = 0; position < 3; ++position) {
            const auto group = position == 0 ? 0 : (tail ? S::tail_packs : S::blocks)-1;
            const auto count = position == 1 ? std::size_t{5} : std::size_t{3};
            const auto first = position == 0 ? 0 : position == 1 ? S::block_width/2 : S::block_width-count;
            const auto chunk_directory = directory/('c'+std::to_string(trial++));
            const auto file = File::prepareChunk(chunk_directory,*raw,tail,group,first,count);
            std::vector<typename S::Selector> ciphertexts;
            raw->loadChunk(tail,group,first,count,ciphertexts);
            std::vector<Cache> expected, actual;
            TFHEpp::narrow_detail::prepareSelectors<S>(expected,ciphertexts);
            File::loadChunkArchive(file,*raw,tail,group,first,count,actual);
            for (std::size_t i=0;i<count;++i) {
                if (actual.at(i).spectra != expected[i].spectra) throw std::runtime_error("full-degree chunk spectra changed");
                for (const auto &spectrum : actual[i].spectra) total_words+=spectrum.size();
            }
            const auto previous = actual;
            for (const auto wrong_first : {first == 0 ? 1 : first-1, S::block_width}) {
                try { File::loadChunkArchive(file,*raw,tail,group,wrong_first,count,actual); }
                catch (const std::exception &) { ++rejected; continue; }
                throw std::runtime_error("wrong prepared chunk range accepted");
            }
            // A separate chunk archive must not be accepted as a full group.
            if (!tail && group==0) {
                std::filesystem::create_hard_link(file,chunk_directory/"head-0.bin");
                try { File invalid(chunk_directory,raw,typename File::Existing{}); }
                catch (const std::exception &) { ++rejected; }
                if (rejected!=3) throw std::runtime_error("partial cache accepted as complete");
            }
            for (std::size_t i=0;i<count;++i)
                if (actual[i].spectra!=previous[i].spectra) throw std::runtime_error("rejected chunk changed output");
            std::cout << "{\"n\":" << P::n << ",\"tail\":" << (tail?"true":"false")
                      << ",\"group\":" << group << ",\"first\":" << first << ",\"records\":" << count
                      << ",\"exact_words\":" << count * Cache::prime_count * (P::k+1) * Cache::digits * P::n
                      << ",\"passed\":true}" << std::endl;
        }
    std::cout << "{\"total_exact_words\":" << total_words << ",\"rejected_range_or_partial_controls\":"
              << rejected << ",\"fresh_keys_generated\":false,\"passed\":true}" << std::endl;
    return 0;
#else
    throw std::runtime_error("prepared selector files require BLAKE3 support");
#endif
}

int main(int argc,char **argv)
try {
    bool sample = false, selector56 = false;
    int first = 1;
    while (first < argc && std::string(argv[first]).starts_with("--")) {
        const std::string option(argv[first++]);
        if (option == "--sample-chunks" && !sample) sample = true;
        else if (option == "--selector56" && !selector56) selector56 = true;
        else throw std::invalid_argument("unknown or duplicate prepared-cache option");
    }
    if(argc-first!=2) {std::cerr<<"Usage: ckks_narrow_selector_ntt_file [--sample-chunks] [--selector56] CHECKPOINT_DIRECTORY CACHE_DIRECTORY\n";return 2;}
    const auto checkpoint = std::filesystem::path(argv[first]);
    const auto directory = std::filesystem::path(argv[first+1]);
    std::ifstream file(checkpoint/"checkpoint.bin",std::ios::binary);
    if(!file) throw std::runtime_error("cannot open retained checkpoint");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t,12> metadata{};archive(metadata);
    if(metadata[0]!=0x4e41525242433031ULL) throw std::runtime_error("unsupported checkpoint format");
    if(metadata[1]==512) {
        if (selector56) return sample ? sampleChunks<9,true>(checkpoint,directory) : prepare<9,true>(checkpoint,directory);
        return sample ? sampleChunks<9>(checkpoint,directory) : prepare<9>(checkpoint,directory);
    }
    if(metadata[1]==32768) {
        if (selector56) return sample ? sampleChunks<15,true>(checkpoint,directory) : prepare<15,true>(checkpoint,directory);
        return sample ? sampleChunks<15>(checkpoint,directory) : prepare<15>(checkpoint,directory);
    }
    throw std::runtime_error("unsupported prepared-cache degree");
}
catch(const std::exception &error) {std::cerr<<error.what()<<'\n';return 2;}
