#include "ckks_narrow_basis_benchmark.hpp"
#include <iomanip>

#if defined(USE_BLAKE3)
template <class Chain, std::size_t... I>
constexpr std::size_t chainBytes(std::index_sequence<I...>)
{
    using Keys = decltype(std::declval<Chain>().keys);
    return ((7*sizeof(std::uint64_t)+sizeof(std::tuple_element_t<I,Keys>::data)) + ...);
}

template <std::uint32_t Nbit>
int run(const std::filesystem::path &checkpoint, const std::filesystem::path &directory,
        const std::filesystem::path &retained_cache = {}, bool fresh_only = false,
        std::size_t columns = 8, bool compare_prefetch = false, bool compare_wide = false)
{
    using namespace TFHEpp;
    using P = narrow_benchmark::LeanBasisParam<Nbit>;
    using S = CKKSNarrowBootstrapSchedule<P,48,518,47,39,31,true,44,true>;
    using Refresh = CKKSNarrowRefreshSchedule<S,44,37,0,true>;
    using Basis = CKKSNarrowBasisRefreshSchedule<Refresh,31>;
    using Square = CKKSCompactDDRelinKey<P,Basis::output_log_q>;
    using Product = CKKSCompactDDRelinKey<P,Basis::root_log_q>;
    using Auto = CKKSAutoKey<P,S::factor_key_log_q>;
    using Chain = typename S::RelinKeys;
    using FirstSwitch = std::tuple_element_t<0,decltype(std::declval<Chain>().keys)>;
    using Raw = CKKSNarrowSelectorFile<S>;
    using File = CKKSNarrowNTTSelectorFile<S>;
    using Cache = CKKSNarrowNTTSelector<S>;
    using Loader = std::function<void(std::size_t,std::size_t,std::vector<Cache>&)>;
    using RawLoader = std::function<void(std::size_t,std::size_t,std::vector<typename S::Selector>&)>;
    std::ifstream file(checkpoint/"checkpoint.bin",std::ios::binary);
    if (!file) throw std::runtime_error("cannot open retained basis factor checkpoint");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t,12> metadata{};archive(metadata);
    const std::array<std::uint64_t,12> expected{0x4e41525242433031ULL,P::n,576,518,47,31,85,44,18,8,44,37};
    if (metadata!=expected) throw std::runtime_error("incompatible retained basis factor profile");
    auto dense=std::make_unique<Key<P>>(), sparse=std::make_unique<Key<P>>();
    archive(*dense,*sparse);
    for(const auto &coefficient:*dense)
        if(coefficient!=typename P::T{} && coefficient!=typename P::T{1} &&
           coefficient!=typename P::T{}-typename P::T{1}) throw std::runtime_error("nonternary retained factor secret");
    const auto positions=narrow_detail::structuredPositions<S>(*sparse);
    const auto prefix=file.tellg();file.seekg(0,std::ios::end);const auto end=file.tellg();
    constexpr auto square_bytes=7*sizeof(std::uint64_t)+sizeof(Square::data);
    constexpr auto basis_bytes=5*sizeof(std::uint64_t)+7*sizeof(std::uint64_t)+sizeof(Product::data)+
                               sizeof(CKKSAutoKey<P,Basis::output_log_q>::data);
    const auto square_begin=end-static_cast<std::streamoff>(square_bytes);
    const auto basis_begin=square_begin-static_cast<std::streamoff>(basis_bytes);
    const auto chain_begin=basis_begin-static_cast<std::streamoff>(sizeof(std::uint32_t)+
        chainBytes<Chain>(std::make_index_sequence<S::tree_depth>{}));
    const auto auto_end=chain_begin-static_cast<std::streamoff>(sizeof(CKKSAutoKey<P,S::output_key_log_q>::data));
    const auto auto_begin=auto_end-static_cast<std::streamoff>(sizeof(Auto::data));
    if (auto_begin<prefix) throw std::runtime_error("truncated retained basis factor suffix");
    const auto check=[&](std::streamoff offset,const auto &expected_anchor) {
        file.seekg(offset);auto anchor=expected_anchor;archive(anchor);
        if(anchor!=expected_anchor) throw std::runtime_error("invalid retained basis factor suffix anchor at byte "+std::to_string(offset));
    };
    check(square_begin,std::array<std::uint64_t,7>{0x434b4b5343444431ULL,P::n,1,Basis::output_log_q,18,8,1});
    check(basis_begin,std::array<std::uint64_t,5>{0x434b4b5346425331ULL,P::n,Basis::root_log_q,Basis::output_log_q,Basis::log_delta});
    check(basis_begin+static_cast<std::streamoff>(5*sizeof(std::uint64_t)),
          std::array<std::uint64_t,7>{0x434b4b5343444431ULL,P::n,1,Basis::root_log_q,18,8,1});
    check(chain_begin,std::array<std::uint64_t,7>{0x434b4b5343444431ULL,P::n,1,FirstSwitch::log_q,18,8,1});
    file.seekg(basis_begin-static_cast<std::streamoff>(sizeof(std::uint32_t)));
    std::uint32_t mode{};archive(mode);
    if(mode>2 || file.tellg()!=basis_begin) throw std::runtime_error("invalid retained basis factor cache mode");
    auto key=std::make_unique<Auto>();file.seekg(auto_begin);archive(*key);
    if(file.tellg()!=auto_end) throw std::runtime_error("retained factor conjugation length mismatch");
    file.close();auto prepared_key=makeCKKSPreparedAutoKeyFFT<P,S::factor_key_log_q>(*key);key.reset();
    if(!std::filesystem::create_directory(directory)) throw std::runtime_error("factor cache control directory must be new");
    std::filesystem::permissions(directory,std::filesystem::perms::owner_all);
    Raw source(checkpoint/"selectors",typename Raw::Existing{});
    const std::size_t chunk_size=(compare_prefetch || compare_wide) ? 40 : 4;
    std::mt19937_64 random(20262119);std::uniform_real_distribution<double> value(-0.5,0.5);
    for(bool tail:{false,true}) {
        const auto group=tail ? S::tail_packs-1 : S::blocks-1;
        const auto cache_directory=directory/(tail?"tail":"head");std::filesystem::create_directory(cache_directory);
        const auto setup_start=std::chrono::steady_clock::now();
        std::vector<std::filesystem::path> chunks;
        std::size_t reused_chunks=0;
        for(std::size_t first=0;!fresh_only && first<S::block_width;first+=chunk_size) {
            const auto count=std::min(chunk_size,S::block_width-first);
            const auto retained=retained_cache/(tail?"tail":"head")/std::to_string(first)/"chunk.bin";
            if (!retained_cache.empty() && std::filesystem::exists(retained)) {
                chunks.push_back(retained); ++reused_chunks;
            }
            else chunks.push_back(File::prepareChunk(cache_directory/std::to_string(first),source,tail,group,first,count));
        }
        const double setup_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-setup_start).count();
        std::vector<CKKSSlotVector<P>> tables(columns);
        auto input=std::make_unique<std::array<double,P::n>>();
        for(std::size_t column=0;column<columns;++column) {
            for(auto &coefficient:*input) coefficient=value(random);
            tables[column]=tail ? narrow_detail::tailTable<S>(*input,group,column+1)
                               : narrow_detail::headTable<S>(*input,group,column+1);
        }
        const Loader cold=[&](std::size_t first,std::size_t count,auto &result) {
            std::vector<typename S::Selector> raw;source.loadChunk(tail,group,first,count,raw);
            narrow_detail::prepareSelectors<S>(result,raw);
        };
        const Loader cached=[&](std::size_t first,std::size_t count,auto &result) {
            File::loadChunkArchive(chunks.at(first/chunk_size),source,tail,group,first,count,result);
        };
        const RawLoader raw_loader=[&](std::size_t first,std::size_t count,auto &result) {
            source.loadChunk(tail,group,first,count,result);
        };
        const auto reference_path = directory/(tail?"tail-reference.bin":"head-reference.bin");
        const std::array<std::uint64_t,3> reference_header{0x4e41524643543031ULL,P::n,columns};
        bool have_reference = false;
        const std::vector<bool> trials = (compare_prefetch || compare_wide) ? std::vector<bool>{false,true}
                                      : fresh_only ? std::vector<bool>{false,false}
                                                  : std::vector<bool>{false,true,true,false};
        for(bool use_cache:trials) {
            CKKSNarrowFactorTimings timings;
            const auto start=std::chrono::steady_clock::now();
            const bool wide = compare_wide && (use_cache || compare_prefetch);
            const auto actual=wide
                ? narrow_detail::packedFactorsBatch<narrow_benchmark::Selector56BasisSchedule<S>>(tables,{}, {},0,*prepared_key,
                    chunk_size,prepared_key.get(),&timings,{},raw_loader,
                    compare_prefetch && use_cache)
                : fresh_only ? narrow_detail::packedFactorsBatch<S>(tables,{}, {},0,*prepared_key,
                    chunk_size,prepared_key.get(),&timings,{},raw_loader,
                    compare_prefetch && use_cache)
                : narrow_detail::packedFactorsBatch<S>(tables,{}, {},0,*prepared_key,
                    chunk_size,prepared_key.get(),&timings,use_cache?cached:cold);
            const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            if(!have_reference) {
                std::ofstream saved(reference_path,std::ios::binary);
                if(!saved) throw std::runtime_error("cannot save factor reference");
                std::filesystem::permissions(reference_path,std::filesystem::perms::owner_read |
                    std::filesystem::perms::owner_write,std::filesystem::perm_options::replace);
                {
                    cereal::PortableBinaryOutputArchive writer(saved);
                    writer(reference_header);
                    for(const auto &factor:actual) writer(factor);
                }
                saved.flush();
                if(!saved) throw std::runtime_error("cannot flush factor reference");
                have_reference=true;
            }
            else {
                std::ifstream saved(reference_path,std::ios::binary);
                if(!saved) throw std::runtime_error("cannot read factor reference");
                cereal::PortableBinaryInputArchive reader(saved);
                std::array<std::uint64_t,3> metadata{};reader(metadata);
                if(metadata!=reference_header) throw std::runtime_error("factor reference profile mismatch");
                auto reference=std::make_unique<typename S::Factor>();
                for(std::size_t i=0;i<columns;++i) {
                    reader(*reference);
                    if(actual[i].ct!=reference->ct) throw std::runtime_error("cached full-degree factor changed ciphertext");
                }
                if(saved.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("unexpected factor reference suffix");
            }
            std::cout<<std::setprecision(17)<<"{\"n\":"<<P::n<<",\"tail\":"<<(tail?"true":"false")
                <<",\"group\":"<<group<<",\"cached\":"<<(!compare_prefetch && !compare_wide && use_cache?"true":"false")
                <<",\"raw_prefetch\":"<<(compare_prefetch && use_cache?"true":"false")
                <<",\"matrix_digit_bits\":"<<(wide ? 56 : Cache::digit_bits)
                <<",\"schedule_configured_width\":"<<(wide ? "true" : "false")
                <<",\"selector_digits\":"<<(wide ? CKKSNarrowNTTSelector<S,56>::digits : Cache::digits)
                <<",\"selector_chunk_size\":"<<chunk_size
                <<",\"columns\":"<<columns<<",\"seconds\":"<<elapsed<<",\"cache_setup_seconds\":"<<setup_seconds
                <<",\"reused_prepared_chunks\":"<<reused_chunks
                <<",\"selector_provider_seconds\":"<<timings.selector_load_seconds
                <<",\"selector_ntt_seconds\":"<<timings.selector_ntt_seconds
                <<",\"matrix_seconds\":"<<timings.matrix_seconds<<",\"ciphertexts_exact\":true}"<<std::endl;
        }
        auto expected_slots=std::make_unique<CKKSSlotVector<P>>(),decoded=std::make_unique<CKKSSlotVector<P>>(),rotated=std::make_unique<CKKSSlotVector<P>>();
        std::ifstream saved(reference_path,std::ios::binary);
        if(!saved) throw std::runtime_error("cannot decrypt factor reference");
        cereal::PortableBinaryInputArchive reader(saved);
        std::array<std::uint64_t,3> metadata{};reader(metadata);
        if(metadata!=reference_header) throw std::runtime_error("factor reference decrypt profile mismatch");
        auto reference=std::make_unique<typename S::Factor>();
        double maximum=0;
        for(std::size_t column=0;column<columns;++column) {
            reader(*reference);
            expected_slots->fill({});
            const auto count=tail?std::min(S::tail_blocks_per_pack,S::blocks-group*S::tail_blocks_per_pack):1;
            for(std::size_t ell=0;ell<count;++ell) {
                const auto block=tail?group*S::tail_blocks_per_pack+ell:group;
                const auto position=positions[block];rotateCKKSSlotVector<P>(*rotated,tables[column],-position);
                const bool negative=(*sparse)[block*S::block_width+position]==typename P::T{}-typename P::T{1};
                const auto begin=tail?(2*ell+1)*S::tail_width:0;
                const auto length=tail?S::tail_width:S::head_slots;
                for(std::size_t i=begin;i<begin+length;++i) (*expected_slots)[i]=negative?std::conj((*rotated)[i]):(*rotated)[i];
            }
            ckksSlotDecrypt<P>(*decoded,*reference,*dense);
            for(std::size_t i=0;i<S::slots;++i) maximum=narrow_benchmark::finiteMaximum(maximum,std::abs((*decoded)[i]-(*expected_slots)[i]));
        }
        if(saved.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("unexpected factor reference decrypt suffix");
        if(maximum>=std::ldexp(1.0,-23)) throw std::runtime_error("full cached factor failed23-bit gate");
        std::cout<<"{\"n\":"<<P::n<<",\"tail\":"<<(tail?"true":"false")<<",\"slots_checked\":"
                 <<columns*S::slots<<",\"maximum_error\":"<<maximum<<",\"accuracy_passed\":true}"<<std::endl;
    }
    return 0;
}
#endif

int main(int argc,char **argv)
try {
#if defined(USE_BLAKE3)
    if((argc!=4 && argc!=5) || (std::string(argv[1])!="--small" && std::string(argv[1])!="--full" &&
                              std::string(argv[1])!="--small-fresh" && std::string(argv[1])!="--full-fresh" &&
                              std::string(argv[1])!="--small-union-fresh" && std::string(argv[1])!="--full-union-fresh" &&
                              std::string(argv[1])!="--small-union-prefetch" && std::string(argv[1])!="--full-union-prefetch" &&
                              std::string(argv[1])!="--small-union-wide" && std::string(argv[1])!="--full-union-wide" &&
                              std::string(argv[1])!="--small-union-wide-prefetch" && std::string(argv[1])!="--full-union-wide-prefetch")) {
        std::cerr<<"Usage: ckks_narrow_basis_factor_cache_checkpoint --small|--full CHECKPOINT_DIRECTORY NEW_DIRECTORY [RETAINED_CACHE_DIRECTORY]\n"
                   "       ckks_narrow_basis_factor_cache_checkpoint --small-fresh|--full-fresh CHECKPOINT_DIRECTORY NEW_DIRECTORY\n"
                   "       ckks_narrow_basis_factor_cache_checkpoint --small-union-fresh|--full-union-fresh CHECKPOINT_DIRECTORY NEW_DIRECTORY\n"
                   "       ckks_narrow_basis_factor_cache_checkpoint --small-union-prefetch|--full-union-prefetch CHECKPOINT_DIRECTORY NEW_DIRECTORY\n"
                   "       ckks_narrow_basis_factor_cache_checkpoint --small-union-wide|--full-union-wide CHECKPOINT_DIRECTORY NEW_DIRECTORY\n"
                   "       ckks_narrow_basis_factor_cache_checkpoint --small-union-wide-prefetch|--full-union-wide-prefetch CHECKPOINT_DIRECTORY NEW_DIRECTORY\n";return 2;
    }
    const bool compare_prefetch = std::string(argv[1]).ends_with("-prefetch");
    const bool compare_wide = std::string(argv[1]).find("-wide") != std::string::npos;
    const bool fresh_only = std::string(argv[1]).ends_with("-fresh") || compare_prefetch || compare_wide;
    if (fresh_only && argc != 4) throw std::invalid_argument("fresh-only factor control does not use a prepared cache");
    const auto retained=argc==5 ? std::filesystem::path(argv[4]) : std::filesystem::path{};
    const std::size_t columns = std::string(argv[1]).find("-union-") != std::string::npos ? 22 : 8;
    return std::string(argv[1]).starts_with("--small") ? run<9>(argv[2],argv[3],retained,fresh_only,columns,compare_prefetch,compare_wide)
                                                    : run<15>(argv[2],argv[3],retained,fresh_only,columns,compare_prefetch,compare_wide);
#else
    throw std::runtime_error("prepared selector files require BLAKE3 support");
#endif
}
catch(const std::exception &error) {std::cerr<<error.what()<<'\n';return 2;}
