#include <tfhe++.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <string_view>
#include <iomanip>

using namespace TFHEpp;
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b-a).count();
}
struct Timings {
    double forward=0, mult=0, reverse=0;
    void print(size_t products) const {
        std::cout << std::setprecision(10)
                  << "forward_seconds=" << forward << "\nmult_seconds=" << mult
                  << "\nreverse_seconds=" << reverse
                  << "\nevaluation_seconds=" << forward+mult+reverse
                  << "\nforward_seconds_per_product=" << forward/products
                  << "\nmult_seconds_per_product=" << mult/products
                  << "\nreverse_seconds_per_product=" << reverse/products
                  << "\nevaluation_seconds_per_product=" << (forward+mult+reverse)/products
                  << std::endl;
    }
};
std::vector<std::pair<unsigned,unsigned>> inputs() {
    std::vector<std::pair<unsigned,unsigned>> pairs = {
        {0,255},{1,255},{255,255},{128,128},{127,129},{254,255},{0,0},{255,254}};
    std::mt19937 rng(20260905);
    while (pairs.size()<256) {
        const unsigned x=rng()&255, y=rng()&255;
        pairs.emplace_back(x,y);
    }
    return pairs;
}
void run_clpx(const std::vector<std::pair<unsigned,unsigned>>& pairs) {
    using P=SS2CLPXlvl2param;
    const auto start=Clock::now();
    auto sk=std::make_unique<SecretKey>();
    auto ek=std::make_unique<EvalKey>();
    ek->emplaceiksk<lvl1hparam>(*sk);
    ek->emplaceiksk<lvl21param>(*sk);
    ek->emplaceiksk<lvl2hparam>(*sk);
    ek->emplacebkfft<lvlh1param>(*sk);
    ek->emplacebkfft<SS2CLPXlvlh2param>(*sk);
    ek->emplacebkfft<CLPX2TFHElvlh2param>(*sk);
    const auto& ringkey=sk->key.getIndependent<lvl2param>();
    auto ahk=std::make_unique<AnnihilateKey<P>>();
    annihilatekeygen<P>(*ahk,ringkey);
    auto relin=std::make_unique<relinKeyFFT<P>>(relinKeyFFTgen<P>(ringkey));
    std::cout << "key_generation_seconds=" << seconds(start,Clock::now()) << std::endl;
    Timings timing;
    size_t done=0;
    for (const auto& [x,y]:pairs) {
        std::vector<uint8_t> xd(8),yd(8);
        for (unsigned i=0;i<8;i++) {xd[i]=(x>>i)&1; yd[i]=(y>>i)&1;}
        std::vector<TLWE<lvl1param>> xc,yc;
        const auto& io_key=sk->key.getIndependent<lvl1param>();
        bootsSymEncrypt<lvl1param>(xc,xd,io_key);
        bootsSymEncrypt<lvl1param>(yc,yd,io_key);
        TRLWE<P> a{},b{},product{};
        std::vector<TLWE<lvl1param>> result(16);
        auto t0=Clock::now();
        TLWES2CLPXIKS<lvl1hparam,SS2CLPXlvlh2param,SS2CLPXlvl22param,4,4,20>(a,xc,*ahk,*ek);
        TLWES2CLPXIKS<lvl1hparam,SS2CLPXlvlh2param,SS2CLPXlvl22param,4,4,20>(b,yc,*ahk,*ek);
        auto t1=Clock::now();
        CLPXMult<P>(product,a,b,*relin);
        auto t2=Clock::now();
        CLPX2TLWESIKSanybit<lvl1hparam,lvl21param,lvlh1param,CLPX2TFHElvlh2param,lvl2hparam,9,2>(result,product,*ek,*sk);
        auto t3=Clock::now();
        timing.forward+=seconds(t0,t1); timing.mult+=seconds(t1,t2); timing.reverse+=seconds(t2,t3);
        unsigned decoded=0;
        for (unsigned i=0;i<16;i++) decoded|=unsigned(tlweSymDecrypt<lvl1param>(result[i],io_key))<<i;
        if (decoded!=x*y) throw std::runtime_error("CLPX mismatch: "+std::to_string(x)+"*"+std::to_string(y)+" got "+std::to_string(decoded));
        if (++done%16==0) std::cout << "Checked " << done << "/" << pairs.size() << std::endl;
    }
    timing.print(pairs.size());
}
void run_ckks(const std::vector<std::pair<unsigned,unsigned>>& pairs) {
    using P=CKKSInteger4Param; using L=lvl2param;
    const auto start=Clock::now();
    auto ck=std::make_unique<Key<P>>(); Key<L> tk{};
    keyGen<P>(*ck); keyGen<L>(tk);
    auto ek=std::make_unique<TFHECKKSInteger4Key<P,L>>();
    TFHECKKSInteger4KeyGen(*ek,*ck,tk);
    std::cout << "key_generation_seconds=" << seconds(start,Clock::now()) << std::endl;
    Timings timing; long double maxerror=0; size_t done=0;
    for (const auto& [x,y]:pairs) {
        TLWE<L> xc{},yc{},result{};
        TFHEEncryptInteger<8>(xc,x,tk); TFHEEncryptInteger<8>(yc,y,tk);
        auto a=std::make_unique<CKKSCiphertext<P,192,64>>();
        auto b=std::make_unique<CKKSCiphertext<P,192,64>>();
        auto p=std::make_unique<CKKSCiphertext<P,128,64>>();
        auto down=std::make_unique<CKKSCiphertext<P,104,40>>();
        auto reduced=std::make_unique<CKKSCiphertext<P,48,40>>();
        auto t0=Clock::now();
        TFHEToCKKSInteger4Phase(*a,xc,*ek); TFHEToCKKSInteger4Phase(*b,yc,*ek);
        auto t1=Clock::now();
        CKKSInteger4Multiply(*p,*a,*b,ek->relin);
        auto t2=Clock::now();
        CKKSMulIntegerInPlace<P,128,64>(*p,std::uint64_t{1}<<24);
        CKKSScaleDown<P,128,64,24>(*down,*p);
        CKKSLevelReduce<P,104,48,40>(*reduced,*down);
        CKKSCoeffToTFHE<P,L,48,40,2,24>(result,*reduced,0,*ek->to_tfhe);
        auto t3=Clock::now();
        timing.forward+=seconds(t0,t1); timing.mult+=seconds(t1,t2); timing.reverse+=seconds(t2,t3);
        const unsigned decoded=TFHEDecryptIntegerProduct<8>(result,tk);
        if (decoded!=x*y) throw std::runtime_error("CKKS mismatch");
        const uint64_t diff=tlweSymPhase<L>(result,tk)-(uint64_t{x*y}<<48);
        maxerror=std::max(maxerror,std::ldexp(std::abs(static_cast<long double>(static_cast<int64_t>(diff))),-48));
        if (++done%16==0) std::cout << "Checked " << done << "/" << pairs.size() << std::endl;
    }
    timing.print(pairs.size());
    std::cout << "max_integer_unit_error=" << maxerror << std::endl;
}
#ifdef TFHEPP_DEFAULT_128BIT_PARAMS
// Diagnostic decryption is called only after all evaluation timers stop.
long double bfv_error(const TRLWE<BFVSwitchParam>& ct, unsigned expected,
                      const Key<BFVSwitchParam>& key) {
    using P=BFVSwitchParam;
    const auto phase=trlwePhase<P>(ct,key);
    long double error=0;
    for (unsigned i=0;i<P::n;++i) {
        const auto diff=phase[i]-(i==0 ? P::T{expected}*P::delta_int : P::T{0});
        const auto magnitude=(diff>>(127)) ? P::T{0}-diff : diff;
        error=std::max(error,std::ldexp(static_cast<long double>(magnitude),-112));
    }
    return error;
}
void run_bfv(const std::vector<std::pair<unsigned,unsigned>>& pairs) {
    using P=BFVSwitchParam;
    const auto start=Clock::now();
    Key<lvl1param> io_key;
    Key<lvlhalfparam> half_key;
    Key<P> bfv_key;
    keyGen<lvl1param>(io_key); keyGen<lvlhalfparam>(half_key); keyGen<P>(bfv_key);
    TFHEBFVSwitchKey key;
    TFHEBFVSwitchKeyGen(key,io_key,half_key,bfv_key);
    std::cout << "key_generation_seconds=" << seconds(start,Clock::now()) << '\n'
              << "bfv_ring_dimension=" << P::n << "\nbfv_ciphertext_bits=128"
              << "\nbfv_plaintext_modulus=65536\nbfv_scale_bits=112"
              << "\nbfv_alpha_log2=" << std::log2(P::α)
              << "\nbfv_gadget_levels=" << P::l << "\nbfv_gadget_basebit=" << P::Bgbit
              << "\nbfv_dd_levels=" << P::l̅ << "\nbfv_dd_basebit=" << P::B̅gbit
              << "\nbfv_to_tfhe_levels=" << BFVToTFHEKeySwitchParam::t
              << "\nbfv_to_tfhe_basebit=" << BFVToTFHEKeySwitchParam::basebit
              << "\ntfhe_io_n=" << lvl1param::n
              << "\ntfhe_io_alpha_log2=" << std::log2(lvl1param::α)
              << "\ntfhe_io_levels=" << lvl1param::l
              << "\ntfhe_io_basebit=" << lvl1param::Bgbit
              << "\ntfhe_half_n=" << lvlhalfparam::n
              << "\ntfhe_half_alpha_log2=" << std::log2(lvlhalfparam::α)
              << "\ntfhe_to_half_levels=" << lvl1hparam::t
              << "\ntfhe_to_half_basebit=" << lvl1hparam::basebit
              << "\nevaluation_key_bytes=" << sizeof(*key.to_half)+sizeof(*key.to_bfv)
                    +sizeof(*key.packing)+sizeof(*key.relin)+sizeof(*key.to_tfhe)+sizeof(*key.boolean_bootstrap)
              << "\nio_format=boolean_bits\n" << std::flush;
    Timings timing;
    long double forward_error=0, product_error=0;
    size_t done=0;
    for (const auto& [x,y]:pairs) {
        BFVSwitchInput xc,yc;
        for (unsigned i=0;i<8;++i) {
            tlweSymEncrypt<lvl1param>(xc[i],((x>>i)&1)?lvl1param::μ:-lvl1param::μ,lvl1param::α,io_key);
            tlweSymEncrypt<lvl1param>(yc[i],((y>>i)&1)?lvl1param::μ:-lvl1param::μ,lvl1param::α,io_key);
        }
        auto a=std::make_unique<TRLWE<P>>();
        auto b=std::make_unique<TRLWE<P>>();
        auto product=std::make_unique<TRLWE<P>>();
        BFVSwitchOutput result;
        const auto t0=Clock::now();
        TFHEToBFV(*a,xc,key); TFHEToBFV(*b,yc,key);
        const auto t1=Clock::now();
        TRLWEMultFullDD<P>(*product,*a,*b,*key.relin);
        const auto t2=Clock::now();
        BFVToTFHE(result,*product,key);
        const auto t3=Clock::now();
        timing.forward+=seconds(t0,t1); timing.mult+=seconds(t1,t2); timing.reverse+=seconds(t2,t3);
        unsigned decoded=0;
        for (unsigned i=0;i<16;++i)
            decoded|=unsigned(tlweSymDecrypt<lvl1param>(result[i],io_key))<<i;
        if (decoded!=x*y)
            throw std::runtime_error("BFV mismatch: "+std::to_string(x)+"*"+std::to_string(y)+" got "+std::to_string(decoded));
        forward_error=std::max({forward_error,bfv_error(*a,x,bfv_key),bfv_error(*b,y,bfv_key)});
        product_error=std::max(product_error,bfv_error(*product,x*y,bfv_key));
        if (++done%16==0) std::cout << "Checked " << done << "/" << pairs.size() << std::endl;
    }
    timing.print(pairs.size());
    std::cout << "max_forward_integer_unit_error=" << forward_error
              << "\nmax_product_integer_unit_error=" << product_error << std::endl;
}
#endif
int main(int argc,char** argv) try {
    if (argc!=2) throw std::invalid_argument("usage: scheme_switch_comparison clpx|ckks|bfv");
    const auto pairs=inputs();
    if (std::string_view(argv[1])=="clpx") run_clpx(pairs);
    else if (std::string_view(argv[1])=="ckks") run_ckks(pairs);
    else if (std::string_view(argv[1])=="bfv") {
#ifdef TFHEPP_DEFAULT_128BIT_PARAMS
        run_bfv(pairs);
#else
        throw std::invalid_argument("BFV switching requires default 128-bit parameters (USE_BLOCK_BINARY=OFF)");
#endif
    }
    else throw std::invalid_argument("unknown scheme");
    std::cout << "Products=" << pairs.size() << " Passed" << std::endl;
} catch (const std::exception& e) {std::cerr << e.what() << std::endl;return 1;}
