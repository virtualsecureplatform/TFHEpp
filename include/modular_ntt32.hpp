#pragma once

#include "modular_ntt.hpp"
#include <algorithm>

namespace TFHEpp::modular_ntt32 {

// Each prime is below 2^30 and one modulo 2^17. Four primes support the
// same exact 52-bit selector convolution without increasing spectrum bytes.
inline constexpr std::array<modular_ntt::PrimeModulus, 4> primes{{
    {1073479681, 11}, {1071513601, 23}, {1070727169, 11}, {1068236801, 3}}};
inline constexpr unsigned __int128 crt_product = [] {
    unsigned __int128 result = 1;
    for (auto prime : primes) result *= prime.value;
    return result;
}();

class Montgomery {
    std::uint32_t modulus_, negative_inverse_, r_cubed_;
public:
    explicit Montgomery(std::uint64_t modulus) : modulus_(modulus)
    {
        if (!(modulus & 1) || modulus < 3 || modulus >= (1U << 30))
            throw std::invalid_argument("Montgomery modulus must be odd and below 2^30");
        std::uint32_t inverse = modulus;
        for (unsigned i = 0; i < 5; ++i) inverse *= 2 - modulus * inverse;
        negative_inverse_ = 0U - inverse;
        const std::uint64_t radix = (std::uint64_t{1} << 32) % modulus;
        const auto square = radix * radix % modulus;
        r_cubed_ = square * radix % modulus;
    }
    std::uint32_t modulus() const { return modulus_; }
    // t < p*2^32, hence t + m*p < 2^63 and the result is below 2p.
    std::uint32_t reduceLazy(std::uint64_t t) const
    {
        const std::uint32_t m = std::uint32_t(t) * negative_inverse_;
        return (t + std::uint64_t(m) * modulus_) >> 32;
    }
    std::uint32_t canonical(std::uint32_t value) const
    { return value >= modulus_ ? value - modulus_ : value; }
    std::uint32_t multiplyLazy(std::uint32_t a, std::uint32_t b) const
    { return reduceLazy(std::uint64_t(a) * b); }
    std::uint32_t multiply(std::uint32_t a, std::uint32_t b) const
    { return canonical(multiplyLazy(a, b)); }
    // Accept a complete unsigned selector digit, avoiding division per digit.
    std::uint32_t encode(std::uint64_t value) const
    {
        if (value >= (std::uint64_t(modulus_) << 32))
            throw std::invalid_argument("wide Montgomery input exceeds reduction range");
        return multiply(canonical(reduceLazy(value)), r_cubed_);
    }
    std::uint32_t decode(std::uint32_t value) const
    { return canonical(reduceLazy(value)); }
#if defined(__AVX2__)
    // Interleave eight Montgomery reductions of even/odd64-bit lanes.
    // Each t+m*p must fit uint64 (m is the low-word correction). The usual
    // t<p*2^32 bound returns lazy values below2p. SignedEncoder also uses
    // t<=2^63 with its proved prime bound and normalizes the <4p result.
    __m256i reduceWideVector(__m256i even, __m256i odd) const
    {
        const auto inverse = _mm256_set1_epi32(negative_inverse_);
        const auto prime = _mm256_set1_epi32(modulus_);
        const auto even_m = _mm256_mul_epu32(even, inverse);
        const auto odd_m = _mm256_mul_epu32(odd, inverse);
        even = _mm256_srli_epi64(_mm256_add_epi64(even,
            _mm256_mul_epu32(even_m, prime)), 32);
        odd = _mm256_srli_epi64(_mm256_add_epi64(odd,
            _mm256_mul_epu32(odd_m, prime)), 32);
        return _mm256_or_si256(even, _mm256_slli_epi64(odd, 32));
    }
    __m256i multiplyLazyVector(__m256i input, __m256i constant) const
    {
        const auto even = _mm256_mul_epu32(input, constant);
        const auto odd = _mm256_mul_epu32(_mm256_srli_epi64(input, 32),
                                          _mm256_srli_epi64(constant, 32));
        return reduceWideVector(even, odd);
    }
    __m256i canonicalVector(__m256i value) const
    {
        return _mm256_min_epu32(value,
            _mm256_sub_epi32(value, _mm256_set1_epi32(modulus_)));
    }
    __m256i accumulateVector(__m256i destination, __m256i input,
                            __m256i selector) const
    {
        const auto sum = _mm256_add_epi32(destination,
            multiplyLazyVector(input, selector));
        // Both operands <2p, so sum <4p<2^32. Underflow of subtraction
        // lies above every unreduced value below 2p; unsigned min reduces.
        return _mm256_min_epu32(sum,
            _mm256_sub_epi32(sum, _mm256_set1_epi32(2 * modulus_)));
    }
    // Canonical vector operands: four raw products sum below4p^2<p*2^32.
    // One REDC yields a lazy product below2p, then restores the lazy sum.
    template<class Left,class Right>
    __m256i accumulateDot4Vector(__m256i destination,Left left,Right right)const
    {
        auto even=_mm256_setzero_si256(),odd=even;
        for(std::size_t t=0;t<4;++t){
            const auto x=left(t),y=right(t);
            even=_mm256_add_epi64(even,_mm256_mul_epu32(x,y));
            odd=_mm256_add_epi64(odd,_mm256_mul_epu32(_mm256_srli_epi64(x,32),_mm256_srli_epi64(y,32)));
        }
        const auto sum=_mm256_add_epi32(destination,reduceWideVector(even,odd));
        return _mm256_min_epu32(sum,_mm256_sub_epi32(sum,_mm256_set1_epi32(2*modulus_)));
    }

#endif
    // Four canonical products sum to <4p^2<p*2^32 because p<2^30.
    // REDC therefore returns <2p, just as for a single product. Reduce once
    // per four terms; the destination remains in [0,2p).
    // Eight products have t<8p^2<2^63 and t+m*p<3*2^62, so REDC
    // cannot overflow uint64. Its result is below 3p; one subtraction
    // reduces it below 2p before adding the lazy destination (then below 4p).
    template<std::size_t Terms>
    void addDotPair(std::uint32_t *first,std::uint32_t *second,
        const std::array<const std::uint32_t*,Terms>&left,
        const std::array<const std::uint32_t*,Terms>&right,
        const std::array<std::uint32_t,Terms>&selectors,std::size_t columns)const
    {
        static_assert(Terms==4 || Terms==8);
        std::size_t column=0;
#if defined(__AVX2__)
        const auto inverse=_mm256_set1_epi32(negative_inverse_);
        const auto prime=_mm256_set1_epi32(modulus_);
        const auto twice=_mm256_set1_epi32(2*modulus_);
        const auto dot=[&](const auto&load){
            auto even=_mm256_setzero_si256(),odd=even;
            for(std::size_t t=0;t<Terms;++t){
                const auto values=load(t),fixed=_mm256_set1_epi32(selectors[t]);
                even=_mm256_add_epi64(even,_mm256_mul_epu32(values,fixed));
                odd=_mm256_add_epi64(odd,_mm256_mul_epu32(_mm256_srli_epi64(values,32),fixed));
            }
            const auto em=_mm256_mul_epu32(even,inverse),om=_mm256_mul_epu32(odd,inverse);
            even=_mm256_srli_epi64(_mm256_add_epi64(even,_mm256_mul_epu32(em,prime)),32);
            odd=_mm256_srli_epi64(_mm256_add_epi64(odd,_mm256_mul_epu32(om,prime)),32);
            auto value=_mm256_or_si256(even,_mm256_slli_epi64(odd,32));
            if constexpr(Terms==8)value=canonicalVector(value);
            return value;
        };
        const auto add=[&](__m256i old,__m256i value){
            const auto sum=_mm256_add_epi32(old,value);
            return _mm256_min_epu32(sum,_mm256_sub_epi32(sum,twice));
        };
        const auto pack=[](const std::uint32_t*a,const std::uint32_t*b){
            return _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128(reinterpret_cast<const __m128i*>(a))),_mm_loadu_si128(reinterpret_cast<const __m128i*>(b)),1);
        };
        if(columns==4){
            const auto value=add(pack(first,second),dot([&](std::size_t t){return pack(left[t],right[t]);}));
            _mm_storeu_si128(reinterpret_cast<__m128i*>(first),_mm256_castsi256_si128(value));
            _mm_storeu_si128(reinterpret_cast<__m128i*>(second),_mm256_extracti128_si256(value,1));
            return;
        }
        for(;column+8<=columns;column+=8){
            const auto a=add(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(first+column)),dot([&](std::size_t t){return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(left[t]+column));}));
            const auto b=add(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(second+column)),dot([&](std::size_t t){return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(right[t]+column));}));
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(first+column),a);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(second+column),b);
        }
#endif
        for(;column<columns;++column){
            std::uint64_t a=0,b=0;
            for(std::size_t t=0;t<Terms;++t){a+=std::uint64_t(left[t][column])*selectors[t];b+=std::uint64_t(right[t][column])*selectors[t];}
            auto ar=reduceLazy(a),br=reduceLazy(b);
            if constexpr(Terms==8){ar=canonical(ar);br=canonical(br);}
            const auto av=first[column]+ar,bv=second[column]+br;
            first[column]=av>=2*modulus_?av-2*modulus_:av;
            second[column]=bv>=2*modulus_?bv-2*modulus_:bv;
        }
    }
    // Two ciphertext components share each public table load and odd-lane
    // extraction. All four destinations retain the same lazy [0,2p) bound.
    template<std::size_t Terms,bool Packed=false>
    void addDotPairComponents(std::uint32_t*first0,std::uint32_t*second0,
        std::uint32_t*first1,std::uint32_t*second1,
        const std::array<const std::uint32_t*,Terms>&left,
        const std::array<const std::uint32_t*,Terms>&right,
        const std::array<std::uint32_t,Terms>&selectors0,
        const std::array<std::uint32_t,Terms>&selectors1,std::size_t columns)const
    {
        static_assert(Terms==4 || Terms==8);
        if constexpr(Packed)if(columns!=4)throw std::invalid_argument("packed component dot needs four columns");
        std::size_t column=0;
#if defined(__AVX2__)
        const auto twice=_mm256_set1_epi32(2*modulus_);
        const auto add=[&](__m256i old,__m256i value){
            const auto sum=_mm256_add_epi32(old,value);
            return _mm256_min_epu32(sum,_mm256_sub_epi32(sum,twice));
        };
        const auto pack=[](const std::uint32_t*a,const std::uint32_t*b){
            return _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128(reinterpret_cast<const __m128i*>(a))),_mm_loadu_si128(reinterpret_cast<const __m128i*>(b)),1);
        };
        struct Pair{__m256i first,second;};
        const auto dot=[&](const auto&load){
            auto e0=_mm256_setzero_si256(),o0=e0,e1=e0,o1=e0;
            for(std::size_t t=0;t<Terms;++t){
                const auto row=load(t),odd=_mm256_srli_epi64(row,32);
                const auto a=_mm256_set1_epi32(selectors0[t]),b=_mm256_set1_epi32(selectors1[t]);
                e0=_mm256_add_epi64(e0,_mm256_mul_epu32(row,a));o0=_mm256_add_epi64(o0,_mm256_mul_epu32(odd,a));
                e1=_mm256_add_epi64(e1,_mm256_mul_epu32(row,b));o1=_mm256_add_epi64(o1,_mm256_mul_epu32(odd,b));
            }
            Pair out{reduceWideVector(e0,o0),reduceWideVector(e1,o1)};
            if constexpr(Terms==8){out.first=canonicalVector(out.first);out.second=canonicalVector(out.second);}
            return out;
        };
        if(columns==4){
            const auto value=dot([&](std::size_t t){
                if constexpr(Packed)return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(left[t]));
                else return pack(left[t],right[t]);
            });
            const auto a=add(pack(first0,second0),value.first),b=add(pack(first1,second1),value.second);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(first0),_mm256_castsi256_si128(a));
            _mm_storeu_si128(reinterpret_cast<__m128i*>(second0),_mm256_extracti128_si256(a,1));
            _mm_storeu_si128(reinterpret_cast<__m128i*>(first1),_mm256_castsi256_si128(b));
            _mm_storeu_si128(reinterpret_cast<__m128i*>(second1),_mm256_extracti128_si256(b,1));
            return;
        }
        for(;column+8<=columns;column+=8){
            const auto store=[&](auto*out0,auto*out1,const auto&tables){
                const auto value=dot([&](std::size_t t){return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(tables[t]+column));});
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(out0+column),add(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out0+column)),value.first));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(out1+column),add(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out1+column)),value.second));
            };
            store(first0,first1,left);store(second0,second1,right);
        }
#endif
        for(;column<columns;++column){
            std::uint64_t a0=0,b0=0,a1=0,b1=0;
            for(std::size_t t=0;t<Terms;++t){
                const auto x=left[t][column],y=Packed?left[t][column+4]:right[t][column];
                a0+=std::uint64_t(x)*selectors0[t];b0+=std::uint64_t(y)*selectors0[t];
                a1+=std::uint64_t(x)*selectors1[t];b1+=std::uint64_t(y)*selectors1[t];
            }
            const auto finish=[&](std::uint32_t&dest,std::uint64_t value){
                auto reduced=reduceLazy(value);if constexpr(Terms==8)reduced=canonical(reduced);
                const auto sum=dest+reduced;dest=sum>=2*modulus_?sum-2*modulus_:sum;
            };
            finish(first0[column],a0);finish(second0[column],b0);finish(first1[column],a1);finish(second1[column],b1);
        }
    }
    void addDot4Pair(std::uint32_t *first,std::uint32_t *second,
        const std::array<const std::uint32_t*,4>&left,
        const std::array<const std::uint32_t*,4>&right,
        const std::array<std::uint32_t,4>&selectors,std::size_t columns)const
    {addDotPair<4>(first,second,left,right,selectors,columns);}
    void addDot8Pair(std::uint32_t *first,std::uint32_t *second,
        const std::array<const std::uint32_t*,8>&left,
        const std::array<const std::uint32_t*,8>&right,
        const std::array<std::uint32_t,8>&selectors,std::size_t columns)const
    {addDotPair<8>(first,second,left,right,selectors,columns);}
    // Four columns of both signed tables are adjacent in each eight-word
    // row. Load them together instead of assembling two half vectors.
    template<std::size_t Terms>
    void addPackedDotPair(std::uint32_t *first,std::uint32_t *second,
        const std::array<const std::uint32_t*,Terms>&tables,
        const std::array<std::uint32_t,Terms>&selectors)const
    {
        static_assert(Terms==4 || Terms==8);
#if defined(__AVX2__)
        auto even=_mm256_setzero_si256(),odd=even;
        for(std::size_t t=0;t<Terms;++t){
            const auto row=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(tables[t]));
            const auto fixed=_mm256_set1_epi32(selectors[t]);
            even=_mm256_add_epi64(even,_mm256_mul_epu32(row,fixed));
            odd=_mm256_add_epi64(odd,_mm256_mul_epu32(_mm256_srli_epi64(row,32),fixed));
        }
        auto value=reduceWideVector(even,odd);
        if constexpr(Terms==8)value=canonicalVector(value);
        const auto destination=_mm256_inserti128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(first))),
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(second)),1);
        const auto sum=_mm256_add_epi32(destination,value);
        const auto reduced=_mm256_min_epu32(sum,_mm256_sub_epi32(sum,_mm256_set1_epi32(2*modulus_)));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(first),_mm256_castsi256_si128(reduced));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(second),_mm256_extracti128_si256(reduced,1));
#else
        std::array<const std::uint32_t*,Terms> right;
        for(std::size_t t=0;t<Terms;++t)right[t]=tables[t]+4;
        addDotPair<Terms>(first,second,tables,right,selectors,4);
#endif
    }
    // Tables are ordinary residues; the selector is in Montgomery form.
    // Destinations remain in [0,2p); canonicalize at a matrix tile boundary.
    void addPair(std::uint32_t *first, std::uint32_t *second,
                 const std::uint32_t *left, const std::uint32_t *right,
                 std::size_t columns, std::uint32_t selector) const
    {
        std::size_t column = 0;
#if defined(__AVX2__)
        const auto fixed = _mm256_set1_epi32(selector);
        const auto pair_four = [&](std::size_t offset) {
            const auto *left_first = first + offset;
            const auto *left_second = second + offset;
            const auto *table_first = left + offset;
            const auto *table_second = right + offset;
            const auto pack = [](const std::uint32_t *a, const std::uint32_t *b) {
                return _mm256_inserti128_si256(_mm256_castsi128_si256(
                    _mm_loadu_si128(reinterpret_cast<const __m128i *>(a))),
                    _mm_loadu_si128(reinterpret_cast<const __m128i *>(b)), 1);
            };
            const auto result = accumulateVector(pack(left_first, left_second), pack(table_first, table_second), fixed);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(first + offset), _mm256_castsi256_si128(result));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(second + offset), _mm256_extracti128_si256(result, 1));
        };
        if (columns == 4) { pair_four(0); return; }
        for (; column + 8 <= columns; column += 8)
            for (auto pair : {std::pair{first, left}, std::pair{second, right}}) {
                const auto result = accumulateVector(
                    _mm256_loadu_si256(reinterpret_cast<const __m256i *>(pair.first + column)),
                    _mm256_loadu_si256(reinterpret_cast<const __m256i *>(pair.second + column)), fixed);
                _mm256_storeu_si256(reinterpret_cast<__m256i *>(pair.first + column), result);
            }
        if (columns - column == 4) { pair_four(column); return; }
#endif
        for (; column < columns; ++column) {
            auto a = first[column] + multiplyLazy(left[column], selector);
            auto b = second[column] + multiplyLazy(right[column], selector);
            first[column] = a >= 2 * modulus_ ? a - 2 * modulus_ : a;
            second[column] = b >= 2 * modulus_ ? b - 2 * modulus_ : b;
        }
    }
};

class NegacyclicNTTPlan {
    std::size_t size_;
    bool vector_short_stages_;
    Montgomery arithmetic_;
    std::vector<std::pair<std::size_t, std::size_t>> swaps_;
    std::vector<std::vector<std::uint32_t>> forward_, inverse_;
    std::vector<std::uint32_t> twist_, untwist_;
#if defined(__AVX2__)
    // Pack independent length2/4 butterflies together, or one length8 block.
    // Every lane enters and leaves canonical, in either residue domain.
    template<std::size_t Length,bool DIF>
    void shortStage(std::span<std::uint32_t> values,
                    const std::vector<std::uint32_t>&twiddles)const
    {
        static_assert(Length==2 || Length==4 || Length==8);
        const auto prime=_mm256_set1_epi32(arithmetic_.modulus());
        [[maybe_unused]] const auto fixed=[&]{
            if constexpr(Length==4)
                return _mm256_setr_epi32(twiddles[0],twiddles[1],twiddles[0],twiddles[1],twiddles[0],twiddles[1],twiddles[0],twiddles[1]);
            else if constexpr(Length==8)
                return _mm256_broadcastsi128_si256(_mm_loadu_si128(reinterpret_cast<const __m128i*>(twiddles.data())));
            else return _mm256_setzero_si256();
        }();
        for(std::size_t i=0;i<size_;i+=8){
            const auto input=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(values.data()+i));
            if constexpr(Length==2){
                // The sole twiddle is Montgomery one; skip its multiplication.
                const auto swapped=_mm256_shuffle_epi32(input,_MM_SHUFFLE(2,3,0,1));
                const auto sum=arithmetic_.canonicalVector(_mm256_add_epi32(input,swapped));
                const auto difference=arithmetic_.canonicalVector(_mm256_sub_epi32(_mm256_add_epi32(swapped,prime),input));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data()+i),_mm256_blend_epi32(sum,difference,0xaa));
            }else{
                const auto even=[&]{
                    if constexpr(Length==4)return _mm256_shuffle_epi32(input,_MM_SHUFFLE(1,0,1,0));
                    else return _mm256_permute2x128_si256(input,input,0x00);
                }();
                auto odd=[&]{
                    if constexpr(Length==4)return _mm256_shuffle_epi32(input,_MM_SHUFFLE(3,2,3,2));
                    else return _mm256_permute2x128_si256(input,input,0x11);
                }();
                if constexpr(!DIF)odd=arithmetic_.canonicalVector(arithmetic_.multiplyLazyVector(odd,fixed));
                const auto sum=arithmetic_.canonicalVector(_mm256_add_epi32(even,odd));
                auto difference=_mm256_sub_epi32(_mm256_add_epi32(even,prime),odd);
                if constexpr(DIF)difference=arithmetic_.canonicalVector(arithmetic_.multiplyLazyVector(difference,fixed));
                else difference=arithmetic_.canonicalVector(difference);
                constexpr auto mask=Length==4?0xcc:0xf0;
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data()+i),_mm256_blend_epi32(sum,difference,mask));
            }
        }
    }
    // Keep the large stage loop small enough to inline into forward/inverse.
    template<bool DIF>
    [[gnu::noinline]] void shortStage(std::span<std::uint32_t> values,std::size_t length,
                    const std::vector<std::uint32_t>&twiddles)const
    {
        if(length==2)shortStage<2,DIF>(values,twiddles);
        else if(length==4)shortStage<4,DIF>(values,twiddles);
        else shortStage<8,DIF>(values,twiddles);
    }
#endif
    [[gnu::always_inline]] inline void transform(std::span<std::uint32_t> values, bool inverse,
                   bool permute = true) const
    {
        if (permute)
            for (auto [a, b] : swaps_) std::swap(values[a], values[b]);
        const auto modulus = arithmetic_.modulus();
        const auto &stages = inverse ? inverse_ : forward_;
        std::size_t stage = 0;
        for (std::size_t length = 2; length <= size_; length <<= 1, ++stage) {
#if defined(__AVX2__)
            if(vector_short_stages_ && size_>=8 && length<=8){
                shortStage<false>(values,length,stages[stage]);
                continue;
            }
#endif
            const auto half = length / 2;
            for (std::size_t block = 0; block < size_; block += length) {
                std::size_t j = 0;
#if defined(__AVX2__)
                const auto prime = _mm256_set1_epi32(modulus);
                for (; j + 8 <= half; j += 8) {
                    const auto even = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(values.data()+block+j));
                    const auto odd = arithmetic_.canonicalVector(arithmetic_.multiplyLazyVector(
                        _mm256_loadu_si256(reinterpret_cast<const __m256i *>(values.data()+block+j+half)),
                        _mm256_loadu_si256(reinterpret_cast<const __m256i *>(stages[stage].data()+j))));
                    const auto sum = arithmetic_.canonicalVector(_mm256_add_epi32(even, odd));
                    const auto difference = arithmetic_.canonicalVector(
                        _mm256_sub_epi32(_mm256_add_epi32(even, prime), odd));
                    _mm256_storeu_si256(reinterpret_cast<__m256i *>(values.data()+block+j), sum);
                    _mm256_storeu_si256(reinterpret_cast<__m256i *>(values.data()+block+j+half), difference);
                }
#endif
                for (; j < half; ++j) {
                    const auto even = values[block+j];
                    const auto odd = arithmetic_.multiply(values[block+j+half], stages[stage][j]);
                    values[block+j] = arithmetic_.canonical(even+odd);
                    values[block+j+half] = arithmetic_.canonical(even+modulus-odd);
                }
            }
        }
    }
    // Decimation in frequency leaves the spectrum in bit-reversed order.
    // Pair with the unpermuted inverse DIT to avoid both permutation passes.
    void transformForwardBitReversed(std::span<std::uint32_t> values) const
    {
        const auto modulus = arithmetic_.modulus();
        std::size_t stage = forward_.size();
        for (std::size_t length = size_; length >= 2; length >>= 1) {
            --stage;
#if defined(__AVX2__)
            if(vector_short_stages_ && size_>=8 && length<=8){
                shortStage<true>(values,length,forward_[stage]);
                continue;
            }
#endif
            const auto half = length / 2;
            for (std::size_t block = 0; block < size_; block += length) {
                std::size_t j = 0;
#if defined(__AVX2__)
                const auto prime = _mm256_set1_epi32(modulus);
                for (; j + 8 <= half; j += 8) {
                    const auto even = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(values.data()+block+j));
                    const auto odd = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(values.data()+block+j+half));
                    const auto sum = arithmetic_.canonicalVector(_mm256_add_epi32(even,odd));
                    const auto difference = _mm256_sub_epi32(_mm256_add_epi32(even,prime),odd);
                    const auto product = arithmetic_.canonicalVector(arithmetic_.multiplyLazyVector(
                        difference,_mm256_loadu_si256(reinterpret_cast<const __m256i*>(forward_[stage].data()+j))));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data()+block+j),sum);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data()+block+j+half),product);
                }
#endif
                for (; j < half; ++j) {
                    const auto even=values[block+j],odd=values[block+j+half];
                    values[block+j]=arithmetic_.canonical(even+odd);
                    values[block+j+half]=arithmetic_.multiply(even+modulus-odd,forward_[stage][j]);
                }
            }
        }
    }
    void multiplyRows(std::span<std::uint32_t> values,
                      const std::vector<std::uint32_t> &constants) const
    {
        std::size_t i = 0;
#if defined(__AVX2__)
        for (; i + 8 <= size_; i += 8) {
            const auto result = arithmetic_.canonicalVector(arithmetic_.multiplyLazyVector(
                _mm256_loadu_si256(reinterpret_cast<const __m256i *>(values.data()+i)),
                _mm256_loadu_si256(reinterpret_cast<const __m256i *>(constants.data()+i))));
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(values.data()+i), result);
        }
#endif
        for (; i < size_; ++i) values[i] = arithmetic_.multiply(values[i], constants[i]);
    }
public:
    // Values must be canonical. Ordinary and Montgomery residues both keep
    // their representation through these linear transforms.
    NegacyclicNTTPlan(std::size_t size, modular_ntt::PrimeModulus prime,
                      bool vector_short_stages=true)
        : size_(size), vector_short_stages_(vector_short_stages),
          arithmetic_(prime.value), twist_(size), untwist_(size)
    {
        const auto modulus = arithmetic_.modulus();
        if (!modular_ntt::isPowerOfTwo(size) || size < 2 || (modulus-1) % (2*size))
            throw std::invalid_argument("unsupported small-prime NTT degree");
        const auto psi = modular_ntt::power(prime.primitive_root, (modulus-1)/(2*size), modulus);
        if (modular_ntt::power(psi, size, modulus) != modulus-1)
            throw std::invalid_argument("invalid small-prime NTT root");
        const auto inverse_psi = modular_ntt::invert(psi, modulus);
        const auto inverse_size = modular_ntt::invert(size, modulus);
        std::uint64_t a = 1, b = inverse_size;
        for (std::size_t i = 0; i < size; ++i) {
            twist_[i] = arithmetic_.encode(a); untwist_[i] = arithmetic_.encode(b);
            a = a*psi % modulus; b = b*inverse_psi % modulus;
        }
        for (std::size_t i = 1, reversed = 0; i < size; ++i) {
            std::size_t bit = size/2;
            for (; reversed & bit; bit >>= 1) reversed ^= bit;
            reversed ^= bit;
            if (i < reversed) swaps_.emplace_back(i, reversed);
        }
        for (std::size_t length = 2; length <= size; length <<= 1) {
            forward_.emplace_back(length/2); inverse_.emplace_back(length/2);
            const auto root = modular_ntt::power(psi, 2*size/length, modulus);
            const auto inverse_root = modular_ntt::invert(root, modulus);
            a = b = 1;
            for (std::size_t j = 0; j < length/2; ++j) {
                forward_.back()[j] = arithmetic_.encode(a);
                inverse_.back()[j] = arithmetic_.encode(b);
                a = a*root % modulus; b = b*inverse_root % modulus;
            }
        }
    }
    const Montgomery &arithmetic() const { return arithmetic_; }
    void forward(std::span<std::uint32_t> values) const
    {
        if (values.size() != size_) throw std::invalid_argument("NTT input length mismatch");
        multiplyRows(values, twist_); transform(values, false);
    }
    void inverse(std::span<std::uint32_t> values) const
    {
        if (values.size() != size_) throw std::invalid_argument("NTT input length mismatch");
        transform(values, true); multiplyRows(values, untwist_);
    }
    // Canonical inputs; spectra retain their ordinary/Montgomery domain.
    // These methods must be paired: only the intermediate order changes.
    void forwardBitReversed(std::span<std::uint32_t> values) const
    {
        if(values.size()!=size_)throw std::invalid_argument("NTT input length mismatch");
        multiplyRows(values,twist_);
        transformForwardBitReversed(values);
    }
    void inverseBitReversed(std::span<std::uint32_t> values) const
    {
        if(values.size()!=size_)throw std::invalid_argument("NTT input length mismatch");
        transform(values,true,false);
        multiplyRows(values,untwist_);
    }

};

class FourPrimeCRT {
    std::array<std::uint64_t, 4> inverses_{}, radix_residues_{};
public:
    FourPrimeCRT()
    {
        unsigned __int128 product = primes[0].value;
        for (std::size_t i = 1; i < 4; ++i) {
            inverses_[i] = modular_ntt::invert(product % primes[i].value, primes[i].value);
            product *= primes[i].value;
        }
        for (std::size_t i = 0; i < 4; ++i)
            radix_residues_[i] = (static_cast<unsigned __int128>(1) << 64) % primes[i].value;
    }
    __int128 reconstructSigned(const std::array<std::uint32_t, 4> &residues) const
    {
        unsigned __int128 value = residues[0], product = primes[0].value;
        for (std::size_t i = 1; i < 4; ++i) {
            const auto prime = primes[i].value;
            const auto reduced = (std::uint64_t(value) % prime +
                std::uint64_t(value >> 64) * radix_residues_[i]) % prime;
            const auto digit = (residues[i] + prime - reduced) * inverses_[i] % prime;
            value += product * digit; product *= prime;
        }
        return value > crt_product/2 ? -static_cast<__int128>(crt_product-value)
                                     : static_cast<__int128>(value);
    }
};
} // namespace TFHEpp::modular_ntt32
