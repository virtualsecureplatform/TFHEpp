#pragma once

#include "scheme_switching.hpp"
#include "prepared_gadget_rns32.hpp"
#include <cereal/archives/portable_binary.hpp>
#include <chrono>
#include <span>

namespace TFHEpp {

// Experimental q32 bridge: t(Y) lives in the genuine destination ring,
// embedded by Y=X^(N/n). It is distinct from contiguous zero padding.
template<class P, class L>
struct CKKSToTFHESubringSwitchKey32 {
    static_assert(P::k == 1 && L::k == 1);
    static_assert(is_multilimb_uint_v<typename P::T>);
    static_assert(std::is_same_v<typename L::T, std::uint32_t>);
    static_assert(P::n == (std::uint32_t{1} << P::nbit) && P::nbit <= 15);
    static_assert(std::has_single_bit(L::n) && P::n >= L::n && P::n % L::n == 0);
    static constexpr std::size_t rows = 16, stride = P::n / L::n;
    static constexpr double noise_stddev = 256.0;
    std::array<std::array<std::array<std::uint32_t,P::n>,2>,rows> ciphertexts{};
    template<class Archive> void serialize(Archive &ar) { ar(ciphertexts); }
};

template<class P, class L>
void CKKSToTFHESubringSwitchKeyGen32(CKKSToTFHESubringSwitchKey32<P,L> &key,
    const Key<P> &source, const Key<L> &target)
{
    using K = CKKSToTFHESubringSwitchKey32<P,L>;
    if (ckks_detail::effectiveNoiseStddevAtLevel<P,32>(
            {std::ldexp(K::noise_stddev,-32),0}) != K::noise_stddev)
        throw std::invalid_argument("subring bridge requires coefficient noise sigma 256");
    auto embedded = std::make_unique<Key<P>>();
    for (std::size_t j=0;j<L::n;++j)
        (*embedded)[j*K::stride] = P::T::from_signed_i64(
            ckks_scheme_switch_detail::ternarySecretToInt<L>(target[j]));
    // Validate all source coefficients before entering OpenMP.
    for (const auto &value:source)
        (void)ckks_scheme_switch_detail::ternarySecretToInt<P>(value);
    std::exception_ptr error;
#pragma omp parallel for schedule(static)
    for (std::size_t r=0;r<K::rows;++r) {
        try {
            auto message = std::make_unique<Polynomial<P>>();
            for (std::size_t i=0;i<P::n;++i)
                (*message)[i] = ckks_detail::reduceToLevel<P,32>(source[i] << (30-2*r));
            auto row = std::make_unique<CKKSSeededKeySwitchRow<P,32>>();
            ckks_scheme_switch_detail::exactEncryptPolynomialAtLevel<P,32>(
                *row,*message,*embedded,{std::ldexp(K::noise_stddev,-32),0});
            auto mask = std::make_unique<Polynomial<P>>();
            ckks_detail::fillSeededUniformPolynomialAtLevel<P,32>(*mask,row->mask_seeds[0]);
            for (std::size_t i=0;i<P::n;++i) {
                key.ciphertexts[r][0][i] = static_cast<std::uint32_t>((*mask)[i].limb[0]);
                key.ciphertexts[r][1][i] = static_cast<std::uint32_t>(row->body[i].limb[0]);
            }
        } catch (...) {
#pragma omp critical(ckks_subring_key_generation)
            { if (!error) error = std::current_exception(); }
        }
    }
    if (error) std::rethrow_exception(error);
}

struct CKKSSubringSwitchTimings {
    double switching=0, splitting=0, extraction=0;
};

template<class P, class L>
class PreparedCKKSToTFHESubringSwitch32 {
    using K = CKKSToTFHESubringSwitchKey32<P,L>;
    CKKSPreparedGadgetRNS32<P,32,2> prepared_;
    using Clock = std::chrono::steady_clock;
public:
    explicit PreparedCKKSToTFHESubringSwitch32(const K &key):prepared_(key) {}
    std::size_t byteSize() const { return prepared_.byteSize(); }
    static constexpr std::size_t prime_count = CKKSPreparedGadgetRNS32<P,32,2>::prime_count;
    void switchAndSplit(std::vector<TRLWE<L>> &out,
        const CKKSCiphertext<P,32,29> &input, CKKSSubringSwitchTimings *timings=nullptr) const
    {
        auto start=Clock::now();
        auto product=std::make_unique<TRLWE<P>>();
        prepared_.apply(*product,input.ct[0]);
        if (timings) timings->switching+=std::chrono::duration<double>(Clock::now()-start).count();
        start=Clock::now();
        out.resize(K::stride);
        for (std::size_t r=0;r<K::stride;++r) for (std::size_t j=0;j<L::n;++j) {
            const auto i=r+K::stride*j;
            out[r][0][j]=std::uint32_t{0}-static_cast<std::uint32_t>((*product)[0][i].limb[0]);
            out[r][1][j]=static_cast<std::uint32_t>(input.ct[1][i].limb[0])-static_cast<std::uint32_t>((*product)[1][i].limb[0]);
        }
        if (timings) timings->splitting+=std::chrono::duration<double>(Clock::now()-start).count();
    }
    void apply(std::vector<TLWE<L>> &out, const CKKSCiphertext<P,32,29> &input,
        std::size_t count, CKKSSubringSwitchTimings *timings=nullptr) const
    {
        if (!count || count>P::n) throw std::invalid_argument("invalid subring extraction count");
        std::vector<TRLWE<L>> split;
        switchAndSplit(split,input,timings);
        auto start=Clock::now();
        out.resize(count);
#pragma omp parallel for schedule(static)
        for (std::size_t i=0;i<count;++i)
            SampleExtractIndex<L>(out[i],split[i%K::stride],static_cast<int>(i/K::stride));
        if (timings) timings->extraction+=std::chrono::duration<double>(Clock::now()-start).count();
    }
};

// Fixed-size headers are validated before reading any ciphertext payload.
// Fingerprints are BLAKE3-256 of retained source-checkpoint and owner archives.
template<class P, class L>
auto CKKSSubringBridgeProfile(std::uint64_t kind)
{
    return std::array<std::uint64_t,12>{kind==1?0x53554252494e4732ULL:0x53554252494e4731ULL,kind,P::n,L::n,
        P::T::digits,32,29,2,16,P::n/L::n,std::bit_cast<std::uint64_t>(256.0),3};
}
template<class P,class L,class KeyType>
void CKKSWriteSubringBridgeKey(std::ostream &out,const KeyType &key,std::uint64_t kind,
    const std::string &source,const std::string &target)
{
    if (source.size()!=64 || target.size()!=64) throw std::invalid_argument("invalid bridge origin");
    std::array<char,64> a{},b{};std::copy(source.begin(),source.end(),a.begin());std::copy(target.begin(),target.end(),b.begin());
    cereal::PortableBinaryOutputArchive ar(out);ar(CKKSSubringBridgeProfile<P,L>(kind),a,b,key);
    if (!out) throw std::runtime_error("bridge key write failed");
}
template<class P,class L,class KeyType>
void CKKSReadSubringBridgeKey(std::istream &in,KeyType &key,std::uint64_t kind,
    const std::string &source,const std::string &target)
{
    if (source.size()!=64 || target.size()!=64) throw std::invalid_argument("invalid bridge origin");
    cereal::PortableBinaryInputArchive ar(in);
    std::array<std::uint64_t,12> profile{};std::array<char,64> a{},b{};
    ar(profile,a,b);
    if (profile!=CKKSSubringBridgeProfile<P,L>(kind) || std::string(a.begin(),a.end())!=source || std::string(b.begin(),b.end())!=target)
        throw std::runtime_error("bridge key profile/origin mismatch");
    ar(key);
    if (in.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("bridge key trailing data");
}
} // namespace TFHEpp
