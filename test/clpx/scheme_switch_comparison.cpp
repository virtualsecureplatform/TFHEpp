#include <tfhe++.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <string_view>
#include <iomanip>

using namespace TFHEpp;
using Clock = std::chrono::steady_clock;
using Word = uint64_t;
using Wide = unsigned __int128;
using Pairs = std::vector<std::pair<Word, Word>>;
std::string decimal(Wide value) {
    std::string result;
    do { result.insert(result.begin(), '0' + value % 10); value /= 10; } while (value);
    return result;
}
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
Pairs inputs(unsigned bits, size_t count) {
    const Word maximum = bits == 64 ? ~Word{0} : (Word{1} << bits) - 1;
    const Word midpoint = Word{1} << (bits - 1);
    Pairs pairs = {
        {0,maximum},{1,maximum},{maximum,maximum},{midpoint,midpoint},
        {midpoint-1,(midpoint+1)&maximum},{maximum-1,maximum},{0,0},{maximum,maximum-1}};
    std::mt19937 rng(20260905);
    auto random_word = [&]() {
        Word value = rng();
        if (bits > 32) value |= Word{rng()} << 32;
        return value & maximum;
    };
    while (pairs.size()<count) {
        const Word x=random_word(), y=random_word();
        pairs.emplace_back(x,y);
    }
    pairs.resize(count);
    return pairs;
}
void run_clpx(const Pairs& pairs, unsigned bits) {
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
        std::vector<uint8_t> xd(bits),yd(bits);
        for (unsigned i=0;i<bits;i++) {xd[i]=(x>>i)&1; yd[i]=(y>>i)&1;}
        std::vector<TLWE<lvl1param>> xc,yc;
        const auto& io_key=sk->key.getIndependent<lvl1param>();
        bootsSymEncrypt<lvl1param>(xc,xd,io_key);
        bootsSymEncrypt<lvl1param>(yc,yd,io_key);
        TRLWE<P> a{},b{},product{};
        std::vector<TLWE<lvl1param>> result(2*bits);
        auto t0=Clock::now();
        TLWES2CLPXIKS<lvl1hparam,SS2CLPXlvlh2param,SS2CLPXlvl22param,4,4,20>(a,xc,*ahk,*ek);
        TLWES2CLPXIKS<lvl1hparam,SS2CLPXlvlh2param,SS2CLPXlvl22param,4,4,20>(b,yc,*ahk,*ek);
        auto t1=Clock::now();
        CLPXMult<P>(product,a,b,*relin);
        auto t2=Clock::now();
        CLPX2TLWESIKSanybit<lvl1hparam,lvl21param,lvlh1param,CLPX2TFHElvlh2param,lvl2hparam,9,2>(result,product,*ek,*sk);
        auto t3=Clock::now();
        timing.forward+=seconds(t0,t1); timing.mult+=seconds(t1,t2); timing.reverse+=seconds(t2,t3);
        Wide decoded=0;
        for (unsigned i=0;i<2*bits;i++) decoded|=Wide(tlweSymDecrypt<lvl1param>(result[i],io_key))<<i;
        if (decoded!=Wide{x}*y) throw std::runtime_error("CLPX mismatch: "+std::to_string(x)+"*"+std::to_string(y)+" got "+decimal(decoded));
        if (++done%16==0) std::cout << "Checked " << done << "/" << pairs.size() << std::endl;
    }
    timing.print(pairs.size());
}
void run_ckks(const Pairs& pairs) {
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
template <class P>
long double bfv_error(const TRLWE<P>& ct, Wide expected,
                      const Key<P>& key) {
    const auto phase=trlwePhase<P>(ct,key);
    long double error=0;
    for (unsigned i=0;i<P::n;++i) {
        const auto diff=phase[i]-(i==0 ? typename P::T{expected}*P::delta_int : typename P::T{0});
        const auto magnitude=(diff>>(127)) ? typename P::T{0}-diff : diff;
        error=std::max(error,std::ldexp(static_cast<long double>(magnitude),int(P::plain_modulusbit)-128));
    }
    return error;
}
template <unsigned Bits>
void run_bfv(const Pairs& pairs) {
    using P=BFVSwitchParamFor<Bits>;
    using Reverse=BFVToTFHEKeySwitchParamFor<Bits>;
    const auto start=Clock::now();
    Key<lvl1param> io_key;
    Key<lvlhalfparam> half_key;
    Key<P> bfv_key;
    keyGen<lvl1param>(io_key); keyGen<lvlhalfparam>(half_key); keyGen<P>(bfv_key);
    TFHEBFVSwitchKeyFor<Bits> key;
    TFHEBFVSwitchKeyGen(key,io_key,half_key,bfv_key);
    std::cout << "key_generation_seconds=" << seconds(start,Clock::now()) << '\n'
              << "bfv_ring_dimension=" << P::n << "\nbfv_ciphertext_bits=128"
              << "\nbfv_forward_digit_bits=" << (Bits == 16 ? 2 : 1)
              << "\nbfv_plaintext_modulus=" << decimal(P::plain_modulus)
              << "\nbfv_scale_bits=" << 128-P::plain_modulusbit
              << "\nbfv_alpha_log2=" << std::log2(P::α)
              << "\nbfv_gadget_levels=" << P::l << "\nbfv_gadget_basebit=" << P::Bgbit
              << "\nbfv_dd_levels=" << P::l̅ << "\nbfv_dd_basebit=" << P::B̅gbit
              << "\nbfv_to_tfhe_levels=" << Reverse::t
              << "\nbfv_to_tfhe_basebit=" << Reverse::basebit
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
        std::array<TLWE<lvl1param>, Bits> xc,yc;
        for (unsigned i=0;i<Bits;++i) {
            tlweSymEncrypt<lvl1param>(xc[i],((x>>i)&1)?lvl1param::μ:-lvl1param::μ,lvl1param::α,io_key);
            tlweSymEncrypt<lvl1param>(yc[i],((y>>i)&1)?lvl1param::μ:-lvl1param::μ,lvl1param::α,io_key);
        }
        auto a=std::make_unique<TRLWE<P>>();
        auto b=std::make_unique<TRLWE<P>>();
        auto product=std::make_unique<TRLWE<P>>();
        std::array<TLWE<lvl1param>, 2*Bits> result;
        const auto t0=Clock::now();
        TFHEToBFV<Bits>(*a,xc,key); TFHEToBFV<Bits>(*b,yc,key);
        const auto t1=Clock::now();
        TRLWEMultFullDD<P>(*product,*a,*b,*key.relin);
        const auto t2=Clock::now();
        BFVToTFHE(result,*product,key);
        const auto t3=Clock::now();
        timing.forward+=seconds(t0,t1); timing.mult+=seconds(t1,t2); timing.reverse+=seconds(t2,t3);
        Wide decoded=0;
        for (unsigned i=0;i<2*Bits;++i)
            decoded|=Wide(tlweSymDecrypt<lvl1param>(result[i],io_key))<<i;
        forward_error=std::max({forward_error,bfv_error<P>(*a,x,bfv_key),bfv_error<P>(*b,y,bfv_key)});
        product_error=std::max(product_error,bfv_error<P>(*product,Wide{x}*y,bfv_key));
        if (decoded!=Wide{x}*y) {
            std::cerr << "failed_product_integer_unit_error=" << product_error << std::endl;
            throw std::runtime_error("BFV mismatch: "+std::to_string(x)+"*"+std::to_string(y)+" got "+decimal(decoded));
        }
        if (++done%16==0) std::cout << "Checked " << done << "/" << pairs.size() << std::endl;
    }
    timing.print(pairs.size());
    std::cout << "max_forward_integer_unit_error=" << forward_error
              << "\nmax_product_integer_unit_error=" << product_error << std::endl;
}
#endif
int main(int argc,char** argv) try {
    if (argc!=2 && argc!=4) throw std::invalid_argument("usage: scheme_switch_comparison clpx|ckks|bfv [operand_bits products]");
    auto positive = [](const char* value) {
        size_t end=0;
        const auto number=std::stoul(value,&end);
        if (end!=std::string_view(value).size() || number==0 || number>65536)
            throw std::invalid_argument("invalid numeric argument");
        return number;
    };
    const unsigned bits=argc==4 ? positive(argv[2]) : 8;
    const size_t count=argc==4 ? positive(argv[3]) : 256;
    if (bits>64) throw std::invalid_argument("operand_bits must be in 1..64");
    const auto pairs=inputs(bits,count);
    std::cout << "operand_bits=" << bits << "\nproduct_bits=" << 2*bits << std::endl;
    if (std::string_view(argv[1])=="clpx") run_clpx(pairs,bits);
    else if (std::string_view(argv[1])=="ckks") {
        if (bits!=8) throw std::invalid_argument("this CKKS scalar benchmark supports only 8-bit operands");
        run_ckks(pairs);
    }
    else if (std::string_view(argv[1])=="bfv") {
#ifdef TFHEPP_DEFAULT_128BIT_PARAMS
        switch(bits) {
            case 8: run_bfv<8>(pairs); break;
            case 9: run_bfv<9>(pairs); break;
            case 10: run_bfv<10>(pairs); break;
            case 12: run_bfv<12>(pairs); break;
            case 16: run_bfv<16>(pairs); break;
            case 24: run_bfv<24>(pairs); break;
            case 32: run_bfv<32>(pairs); break;
            default: throw std::invalid_argument("BFV experimental widths: 8,9,10,12,16,24,32");
        }
#else
        throw std::invalid_argument("BFV switching requires default 128-bit parameters (USE_BLOCK_BINARY=OFF)");
#endif
    }
    else throw std::invalid_argument("unknown scheme");
    std::cout << "Products=" << pairs.size() << " Passed" << std::endl;
} catch (const std::exception& e) {std::cerr << e.what() << std::endl;return 1;}
