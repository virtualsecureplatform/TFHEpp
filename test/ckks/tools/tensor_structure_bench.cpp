// Public coefficient fixtures only; no key generation or encryption.
#include "../ckks_narrow_basis_benchmark.hpp"
#include "tensor_structure_candidate.hpp"
#include <iomanip>
#include <iostream>
using namespace TFHEpp;
template<class P,std::uint32_t Q,bool CompactCRT=true,bool DirectCRT=true,bool BatchEncode=true,bool BitReversed=true>
void legacyTensor(TRLWE3<P>&result,const TRLWE<P>&lhs,const TRLWE<P>&rhs)
{
    static_assert(P::k==1 && is_multilimb_uint_v<typename P::T>);
    static_assert(Q>0 && Q<=P::T::digits && Q<490);
    constexpr auto split=(Q+1)/2;
    static_assert(2*split+P::nbit+3<512);
    constexpr auto count=[] {
        MultiLimbUInt<8> product{1};
        // Four products can contribute to a cross-component diagonal; signed
        // recovery needs another factor two. All chunk coefficients <2^split.
        const auto bound=MultiLimbUInt<8>{1}<<(2*split+P::nbit+3);
        for(std::size_t p=0;p<modular_rns32::prime_count;++p){
            product=modular_rns32::multiplySmall(product,modular_rns32::primes[p].value);
            if(product>bound)return p+1;
        }
        return std::size_t{0};
    }();
    static_assert(count>0,"tensor exceeds exact RNS basis range");
    static const auto plans=[] {
        std::vector<modular_ntt32::NegacyclicNTTPlan> out;
        for(std::size_t p=0;p<count;++p)out.emplace_back(P::n,modular_rns32::primes[p]);
        return out;
    }();
    static const modular_rns32::CanonicalEncoder<P::T::limbs,split,count> encoder;
    using Rows=std::array<std::vector<std::uint32_t>,count>;
    std::array<Rows,4> left,right;
#pragma omp parallel for collapse(2) schedule(static)
    for(std::size_t row=0;row<4;++row)for(std::size_t p=0;p<count;++p){
        left[row][p].resize(P::n);right[row][p].resize(P::n);
        const auto &arithmetic=plans[p].arithmetic();
        if constexpr(BatchEncode){
            if(row%2){
                encoder.template residues<Q,split>(left[row][p],lhs[row/2],p);
                encoder.template residues<Q,split,true>(right[row][p],rhs[row/2],p);
            }else{
                encoder.template residues<Q>(left[row][p],lhs[row/2],p);
                encoder.template residues<Q,0,true>(right[row][p],rhs[row/2],p);
            }
        }else{
        for(std::size_t i=0;i<P::n;++i){
            auto a=ckks_detail::reduceToLevel<P,Q>(lhs[row/2][i]);
            auto b=ckks_detail::reduceToLevel<P,Q>(rhs[row/2][i]);
            if(row%2){a>>=split;b>>=split;}
            left[row][p][i]=encoder.residue(a,p);
            right[row][p][i]=arithmetic.encode(encoder.residue(b,p));
        }
        }
        if constexpr(BitReversed){
            plans[p].forwardBitReversed(left[row][p]);plans[p].forwardBitReversed(right[row][p]);
        }else{
            plans[p].forward(left[row][p]);plans[p].forward(right[row][p]);
        }
    }
    for(auto &component:result)component.fill({});
    // Row mapping matches CKKS tensor convention: cross, constant, square.
    constexpr std::array<std::array<int,2>,3> components{{{0,1},{1,1},{0,0}}};
    const auto diagonal_product=[&]<std::size_t Diagonal>(){
        // Only Q-Diagonal*split low bits survive the shift and reduction.
        // CenteredCRT reconstructs modulo its output width, even when the
        // complete CRT product is wider. No sign extension is needed: every
        // omitted bit would lie at or above Q after shifting.
        constexpr auto output_bits=CompactCRT?Q-Diagonal*split:P::T::digits;
        using Recovered=MultiLimbUInt<std::max<std::size_t>(1,(output_bits+63)/64)>;
        using CRT=std::conditional_t<DirectCRT && count>=8,
            modular_rns32::DirectCenteredCRT<Recovered::limbs,count>,
            modular_rns32::CenteredCRT<Recovered::limbs,count>>;
        static const CRT crt;
        std::array<Rows,3> sums;
#pragma omp parallel for collapse(2) schedule(static)
        for(std::size_t c=0;c<3;++c)for(std::size_t p=0;p<count;++p){
            auto &out=sums[c][p];out.assign(P::n,0);
            const auto &arithmetic=plans[p].arithmetic();
            const auto accumulate=[&](int a,int b){
                const auto *x=left[a][p].data(),*y=right[b][p].data();std::size_t i=0;
#if defined(__AVX2__)
                for(;i+8<=P::n;i+=8){
                    const auto value=arithmetic.accumulateVector(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out.data()+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(y+i)));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(out.data()+i),value);
                }
#endif
                for(;i<P::n;++i){const auto value=out[i]+arithmetic.multiplyLazy(x[i],y[i]);out[i]=value>=2*arithmetic.modulus()?value-2*arithmetic.modulus():value;}
            };
            for(std::size_t high=0;high<=Diagonal;++high){
                const auto a=2*components[c][0]+high,b=2*components[c][1]+Diagonal-high;
                accumulate(a,b);if(c==0)accumulate(2*components[c][1]+high,2*components[c][0]+Diagonal-high);
            }
            for(auto &value:out)value=arithmetic.canonical(value);
            if constexpr(BitReversed)plans[p].inverseBitReversed(out);
            else plans[p].inverse(out);
        }
#pragma omp parallel
        {
            std::vector<Recovered> recovered(512);
#pragma omp for collapse(2) schedule(static)
            for(std::size_t c=0;c<3;++c)for(std::size_t begin=0;begin<P::n;begin+=512){
                const auto length=std::min(std::size_t{512},std::size_t(P::n)-begin);
                std::array<std::span<std::uint32_t>,count> residues;
                for(std::size_t p=0;p<count;++p)residues[p]=std::span(sums[c][p]).subspan(begin,length);
                crt.reconstruct(residues,std::span(recovered).first(length));
                for(std::size_t i=0;i<length;++i){
                    typename P::T value{};
                    std::copy_n(recovered[i].limb.begin(),Recovered::limbs,value.limb.begin());
                    result[c][begin+i]=ckks_detail::reduceToLevel<P,Q>(result[c][begin+i]+(value<<(Diagonal*split)));
                }
            }
        }
    };
    diagonal_product.template operator()<0>();
    diagonal_product.template operator()<1>();
}
template<std::uint32_t Q>void benchmark(std::size_t batch,int fixture)
{
 using P=narrow_benchmark::BenchmarkParam<14,7>;
 auto a=std::make_unique<TRLWE<P>>(),b=std::make_unique<TRLWE<P>>();
 std::mt19937_64 random(20261026+Q);
 for(auto*in:{a.get(),b.get()})for(auto&c:*in)for(auto&v:c)for(auto&w:v.limb)w=random();
 if(fixture==1)for(auto&v:(*b)[0])v=(v>>Q)<<Q;
 const auto&rhs=fixture==2?*a:*b;
 auto reference=std::make_unique<TRLWE3<P>>();
 legacyTensor<P,Q>(*reference,*a,rhs);
 std::vector<TRLWE3<P>> out(batch);
 const auto run=[&](bool prune){
  const auto apply=[&](auto&result){
   if(prune)ckks_experiment::tensor<P,Q>(result,*a,rhs);
   else legacyTensor<P,Q>(result,*a,rhs);
  };
  const auto start=std::chrono::steady_clock::now();
  if(batch==1)apply(out[0]);
  else {
#pragma omp parallel for schedule(dynamic,1)
   for(std::size_t j=0;j<batch;++j)apply(out[j]);
  }
  const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  for(const auto&r:out)if(r!=*reference)throw std::runtime_error("public-structure tensor benchmark changed coefficients");
  return seconds;
 };
 run(false);run(true);
 for(bool prune:{false,true,true,false})for(int trial=0;trial<5;++trial)
  std::cout<<"q="<<Q<<" batch="<<batch<<" fixture="<<fixture<<" prune="<<prune<<" trial="<<trial<<" seconds="<<run(prune)<<" identical\n";
}
bool fftReferenceUsable()
{
#if defined(USE_INTERLEAVED_FORMAT) && !defined(__AVX2__)
 // This configured backend emits interleaved complex pairs; the current
 // non-AVX2 MulInFD/FMAInFD fallbacks interpret split real/imaginary halves.
 // Check that layout independently before using FFT as an arithmetic oracle.
 std::array<double,16>a{},b{},product{},accumulated{};
 a[0]=2;a[1]=3;b[0]=5;b[1]=7;
 MulInFD<16>(product,a,b);FMAInFD<16>(accumulated,a,b);
 return product[0]==-11&&product[1]==29&&accumulated[0]==-11&&accumulated[1]==29;
#else
 return true;
#endif
}
template<std::uint32_t Nbit,std::uint32_t Q>void checkPublicStructure()
{
 using P=narrow_benchmark::BenchmarkParam<Nbit,7>;using T=typename P::T;
 auto a=std::make_unique<TRLWE<P>>(),b=std::make_unique<TRLWE<P>>();
 auto actual=std::make_unique<TRLWE3<P>>(),reference=std::make_unique<TRLWE3<P>>();
 std::mt19937_64 random(20261026+Q);
 const bool use_fft=Nbit<=4||fftReferenceUsable();
 const auto inactive=[](auto&poly){for(auto&v:poly)v=(v>>Q)<<Q;};
 for(int fixture=0;fixture<11;++fixture){
  for(auto*in:{a.get(),b.get()})for(auto&c:*in)for(auto&v:c)for(auto&w:v.limb)w=random();
  // Eight active coefficients, including both ends of the ring, permit an
  // independent full-degree schoolbook oracle when this FFT layout is invalid.
  if(!use_fft)for(auto*in:{a.get(),b.get()})for(auto&c:*in)
   for(std::size_t i=4;i<P::n-4;++i)c[i]=(c[i]>>Q)<<Q;
  if(fixture==1||fixture==3||fixture==9)inactive((*a)[0]);
  if(fixture==2||fixture==3)inactive((*b)[0]);
  if(fixture==4)inactive((*a)[1]);
  if(fixture==5)inactive((*b)[1]);
  if(fixture==6||fixture==10)for(auto&c:*a)inactive(c);
  if(fixture==7)for(auto&c:*b)inactive(c);
  const auto&rhs=fixture>=8?*a:*b;
  for(auto&c:*actual)c.fill(T{}-T{1});
  legacyTensor<P,Q>(*reference,*a,rhs);
  ckks_experiment::tensor<P,Q>(*actual,*a,rhs);
  if(*actual!=*reference)throw std::runtime_error("tensor public-structure shortcut changed coefficients");
  if(use_fft){
   auto fft=std::make_unique<TRLWE3<P>>();CKKSTensorProductRescale<P,Q,Q,0>(*fft,*a,rhs);
   if(*actual!=*fft)throw std::runtime_error("structured tensor differs from FFT");
  }
  if(Nbit<=4||!use_fft){
   auto oracle=std::make_unique<TRLWE3<P>>();
   const auto add=[&](std::size_t c,std::size_t x,std::size_t y){
    std::vector<std::size_t> left_indices,right_indices;
    for(std::size_t i=0;i<P::n;++i){
     if(ckks_detail::reduceToLevel<P,Q>((*a)[x][i])!=T{})left_indices.push_back(i);
     if(ckks_detail::reduceToLevel<P,Q>(rhs[y][i])!=T{})right_indices.push_back(i);
    }
    for(const auto i:left_indices)for(const auto j:right_indices){
     const auto value=(*a)[x][i]*rhs[y][j];const auto k=(i+j)%P::n;
     (*oracle)[c][k]=ckks_detail::reduceToLevel<P,Q>((*oracle)[c][k]+(i+j<P::n?value:-value));
    }
   };
   add(0,0,1);add(0,1,0);add(1,1,1);add(2,0,0);
   if(*oracle!=*actual)throw std::runtime_error("structured tensor differs from independent schoolbook oracle");
  }
 }
 std::cout<<"public structure N="<<P::n<<" Q="<<Q<<" zero/transparent/square noncanonical inputs match control and "<<(use_fft?"FFT plus available schoolbook oracle":"full-degree sparse schoolbook oracle")<<'\n' ;
}
int main(int argc,char**argv){
 if(argc>1&&std::string(argv[1])=="--check"){
  if(!fftReferenceUsable())std::cout<<"Configured scalar interleaved FFT reference fails complex-product layout check; use full-degree sparse integer oracle.\n";
  checkPublicStructure<4,1>();checkPublicStructure<4,93>();checkPublicStructure<4,165>();checkPublicStructure<4,385>();checkPublicStructure<4,440>();
  checkPublicStructure<14,93>();checkPublicStructure<14,165>();checkPublicStructure<14,385>();return 0;
 }

 const auto batch=argc>1?std::stoul(argv[1]):1;
 if(batch!=1&&batch!=32)throw std::invalid_argument("batch must be1 or32");
 std::cout<<std::setprecision(17);
 for(int fixture=0;fixture<3;++fixture){benchmark<385>(batch,fixture);benchmark<341>(batch,fixture);benchmark<165>(batch,fixture);}
}
