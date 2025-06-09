export module tfhepp:raintt;
import std;
// No direct TFHEpp module imports seem necessary for raintt namespace itself,
// but it defines types like SWord, DoubleSWord used by PolynomialRAINTT in params.ixx.
// If TFHEpp::cuHEpp::INTorus were used here, it would be `import tfhepp:INTorus;`

export namespace raintt {
export template <typename T>
constexpr T ipow(T num, unsigned int pow)
{
    return (pow >= sizeof(unsigned int) * 8) ? 0
           : pow == 0                        ? 1
                                             : num * ipow(num, pow - 1);
}
export module tfhepp:raintt;
import std;
// This file defines types like SWord, DoubleSWord which are used by PolynomialRAINTT in params.ixx.
// It does not seem to import other tfhepp modules directly.

export namespace raintt {
export template <typename T>
constexpr T ipow(T num, unsigned int pow)
{
    return (pow >= sizeof(unsigned int) * 8) ? 0
           : pow == 0                        ? 1
                                             : num * ipow(num, pow - 1);
}
#ifdef USE_COMPRESS
export constexpr uint min_wordbits = 27;
#else
export constexpr uint min_wordbits = 31;
#endif

#ifdef __clang__
// Currently _BigInt is only implemented in clang
#ifdef USE_COMPRESS
export constexpr uint wordbits = 27;
#else
export constexpr uint wordbits = 31;
#endif
//_ExtInt is equivalent to _BitInt in C23
export using Word = unsigned _BitInt(wordbits);
export using SWord = signed _BitInt(wordbits);
export using DoubleWord = unsigned _BitInt(2 * wordbits);
export using DoubleSWord = signed _BitInt(2 * wordbits);
#else
export constexpr uint wordbits = 32;
export using Word = uint32_t;
export using SWord = int32_t;
export using DoubleWord = uint64_t;
export using DoubleSWord = int64_t;
#endif
export constexpr uint k = 5;
export constexpr uint radixbit = 3;
export constexpr uint radixs2 = 1U << (radixbit - 1);
export constexpr Word K_const = ipow<Word>(k, radixs2); // Renamed K to K_const to avoid conflict with P::k
#ifdef USE_COMPRESS
export constexpr uint shiftunit = 4;
#else
export constexpr uint shiftunit = 5;
#endif
export constexpr uint shiftamount = radixs2 * shiftunit;
export constexpr SWord shiftval = 1 << shiftamount;
export constexpr Word wordmask = (1ULL << wordbits) - 1;
export constexpr SWord P = (K_const << shiftamount) + 1;

export constexpr SWord R = (1ULL << wordbits) % P;
export constexpr SWord R2 = (static_cast<DoubleWord>(R) * R) % P;
export constexpr SWord R3 = (static_cast<DoubleWord>(R2) * R) % P;
export constexpr SWord R4 = (static_cast<DoubleWord>(R2) * R2) % P;

// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
export constexpr Word REDC(const DoubleWord T)
{
    const Word T0 = T & wordmask;
    const Word m = (((T0 * K_const) << shiftamount) - T0) & wordmask; // Use K_const
    const Word t =
        (T + ((static_cast<DoubleWord>(m) * K_const) << shiftamount) + m) >> wordbits; // Use K_const
    return (t > P) ? t - P : t;
}

// https://eprint.iacr.org/2018/039
export inline SWord SREDC(const DoubleSWord a)
{
    const Word a0 = a & wordmask;
    const SWord a1 = a >> wordbits;
    const SWord m = -((a0 * K_const) * shiftval) + a0; // Use K_const
    const SWord t1 =
        (((static_cast<DoubleSWord>(m) * K_const) << shiftamount) + m) >> wordbits; // Use K_const
    return a1 - t1;
}

export inline SWord AddMod(const DoubleSWord a, const DoubleSWord b)
{
    DoubleSWord add = a + b;
    if (add >= P)
        return add - P;
    else if (add <= -P)
        return add + P;
    else
        return add;
}

export inline SWord SubMod(const DoubleSWord a, const DoubleSWord b)
{
    DoubleSWord sub = a - b;
    if (sub >= P)
        return sub - P;
    else if (sub <= -P)
        return sub + P;
    else
        return sub;
}

export constexpr Word MulREDC(const Word a, const Word b)
{
    const DoubleWord mul = static_cast<DoubleWord>(a) * b;
    return REDC(mul);
}

export inline SWord MulSREDC(const SWord a, const SWord b)
{
    const DoubleSWord mul = static_cast<DoubleSWord>(a) * b;
    return SREDC(mul);
}

export constexpr Word PowREDC(const Word a, const uint e)
{
    Word res = 1;
    const Word aR = MulREDC(R2, a);
    for (uint i = 0; i < e; i++) res = MulREDC(res, aR);
    return res;
}

export template <Word a, Word b>
constexpr Word Bit_gcd(SWord &x, SWord &y)
{
    if constexpr (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    else {
        Word d = Bit_gcd<b, a % b>(y, x);
        y -= a / b * x;
        return d;
    }
}

export template <Word val_a> // Renamed template parameter to avoid conflict
constexpr Word inv_mod()
{
    SWord x, y;
    const Word g = Bit_gcd<val_a, P>(x, y);
    if (g != 1) {
        // Consider how to handle errors in a constexpr context if needed, or ensure g is always 1.
        // For now, keeping throw, but this won't work in all constexpr evaluations.
        throw "Inverse doesn't exist";
    }
    else {
        return (x % P + P) % P;  // this line ensures x is positive
    }
}

// NTT related
#if defined(__clang__) & defined(USE_COMPRESS)
export constexpr SWord W = PowREDC(31, K_const); // Use K_const
#else
export constexpr SWord W = PowREDC(11, K_const); // Use K_const
#endif

export template <uint8_t bit>
uint32_t BitReverse(uint32_t in)
{
    if constexpr (bit > 1) {
        const uint32_t center = in & ((bit & 1) << (bit / 2));
        return (BitReverse<bit / 2>(in & ((1U << (bit / 2)) - 1))
                << (bit + 1) / 2) |
               center | BitReverse<bit / 2>(in >> ((bit + 1) / 2));
    }
    else {
        return in;
    }
}

export template <uint Nbit, uint radixbit_param> // Renamed radixbit to avoid conflict
std::unique_ptr<std::array<std::array<SWord, 1U << Nbit>, 2>> TwistGen()
{
    constexpr uint N = 1U << Nbit;
    // constexpr uint8_t remainder = ((Nbit - 1) % radixbit) + 1;
    const Word invN = inv_mod<N>();

    std::unique_ptr<std::array<std::array<SWord, 1U << Nbit>, 2>> twist =
        std::make_unique<std::array<std::array<SWord, 1U << Nbit>, 2>>();
    const Word wR = MulREDC(PowREDC(W, 1U << (shiftamount - Nbit - 1)), R2);
    (*twist)[1][0] = R;
    for (uint i = 1; i < N; i++)
        (*twist)[1][i] = MulREDC((*twist)[1][i - 1], wR);
    // assert(((*twist)[1][N - 1] * w)== 1);

    (*twist)[0][N - 1] = MulREDC(MulREDC((*twist)[1][N - 1], wR),
                                 static_cast<DoubleWord>(invN) * wR % P);
    (*twist)[0][0] = (static_cast<DoubleWord>(invN) * R) % P;
    for (uint32_t i = 2; i < N; i++)
        (*twist)[0][N - i] = MulREDC((*twist)[0][N - i + 1], wR);
    assert(MulREDC((*twist)[0][1], wR) == (*twist)[0][0]);

    // if constexpr(remainder!=1) for (uint i = 0; i < N; i++) (*twist)[1][i] =
    // MulREDC((*twist)[1][i], R2);
    if constexpr (radixbit_param != 1) { // Use renamed template parameter
        // for(uint j = 0; j < (1U << radixbit-1); j++) for (uint i = 0; i <
        // N>>radixbit; i++)  (*twist)[0][(2*j+1)*(N>>radixbit)+i] =
        // MulREDC((*twist)[0][(2*j+1)*(N>>radixbit)+i], R2); constexpr uint
        // radixbit2 = 2; for(uint j = 0; j < (1U << radixbit2-1); j++) for
        // (uint i = 0; i < N>>radixbit; i++)
        // (*twist)[0][(2*j+1)*(N>>radixbit)+i] =
        // MulREDC((*twist)[0][(2*j+1)*(N>>radixbit)+i], R2);
        for (uint i = 0; i < N; i++)
            if (((i >> (Nbit - radixbit_param)) & ((1 << (radixbit_param - 1)) - 1)) != 0) // Use renamed template parameter
                (*twist)[0][i] = MulREDC((*twist)[0][i], R2);
    }
    return twist;
}

export template <uint32_t Nbit>
inline std::unique_ptr<
    std::array<std::array<std::array<SWord, 1U << Nbit>, 2>, 2>>
TableGen()
{
    constexpr uint32_t N = 1U << Nbit;

    std::unique_ptr<std::array<std::array<std::array<SWord, N>, 2>, 2>> table =
        std::make_unique<std::array<std::array<std::array<SWord, N>, 2>, 2>>();
    const Word w = PowREDC(W, 1ULL << (shiftamount - Nbit));
    const Word wR = MulREDC(w, R2);
    (*table)[0][0][0] = (*table)[1][0][0] = R;
    (*table)[0][1][0] = (*table)[1][1][0] = R2;
    for (uint32_t i = 1; i < N; i++)
        (*table)[1][0][i] = MulREDC((*table)[1][0][i - 1], wR);
    assert(MulREDC((*table)[1][0][N - 1], wR) == R);
    for (uint32_t i = 1; i < N; i++)
        (*table)[1][1][i] = MulREDC((*table)[1][0][i], R2);
    for (int j = 0; j < 2; j++)
        for (uint32_t i = 1; i < N; i++)
            (*table)[0][j][i] = (*table)[1][j][N - i];
    return table;
}

export inline void ButterflyAddBothMod(DoubleSWord *const res, const uint size)
{
    for (uint index = 0; index < size / 2; index++) {
        const SWord temp = res[index];
        res[index] = AddMod(res[index], res[index + size / 2]);
        res[index + size / 2] = SubMod(temp, res[index + size / 2]);
    }
}

export inline void ButterflyAddAddMod(DoubleSWord *const res, const uint size)
{
    for (uint index = 0; index < size / 2; index++) {
        const SWord temp = res[index];
        res[index] = AddMod(res[index], res[index + size / 2]);
        res[index + size / 2] = temp - res[index + size / 2];
    }
}

export inline void ButterflyAdd(DoubleSWord *const res, const uint size)
{
    for (uint index = 0; index < size / 2; index++) {
        const SWord temp = res[index];
        res[index] += res[index + size / 2];
        res[index + size / 2] = temp - res[index + size / 2];
    }
}

export inline void ButterflyAddBothSREDC(DoubleSWord *const res, const uint size)
{
    for (uint index = 0; index < size / 2; index++) {
        const DoubleSWord temp = res[index];
        res[index] = SREDC(res[index] + res[index + size / 2]);
        res[index + size / 2] = SREDC(temp - res[index + size / 2]);
    }
}

export template <uint8_t radixbit_param, uint num> // Renamed radixbit
DoubleSWord ConstTwiddleMul(const DoubleSWord a)
{
    static_assert(num <= 3);
    return (a * ipow<DoubleSWord>(k, num * (radixs2 >> (radixbit_param - 1)))) // Use K_const if k is the same as global k
           << (num * (shiftamount >> (radixbit_param - 1)));
}

export template <uint Nbit, uint8_t radixbit_param> // Renamed radixbit
inline void TwiddleMulInvert(
    DoubleSWord *const res, const uint size, const uint blockindex,
    const uint stride, const std::array<std::array<SWord, 1 << Nbit>, 2> &table)
{
    for (uint32_t index = 0; index < size; index++)
        res[index] =
            MulSREDC(res[index], table[blockindex > 1][stride * index]);
}

export template <uint8_t radixbit_param> // Renamed radixbit
inline void INTTradixButterfly(DoubleSWord *const res, const uint32_t size)
{
    static_assert(radixbit_param <= 3, "radix 8 is the maximum!"); // Use renamed
    if constexpr (radixbit_param == 1) { // Use renamed
        ButterflyAddBothMod(res, size);
    }
    else if constexpr (radixbit_param == 2) { // Use renamed
        ButterflyAddAddMod(res, size);
        INTTradixButterfly<radixbit_param - 1>(&res[0], size / 2); // Use renamed
        const uint32_t block = size >> radixbit_param; // Use renamed
        for (int i = 1; i < (1 << (radixbit_param - 1)); i++) // Use renamed
            for (int j = 0; j < block; j++)
                res[i * block + j + size / 2] =
                    ConstTwiddleMul<radixbit_param, 1>(res[i * block + j + size / 2]); // Use renamed
        ButterflyAddBothSREDC(&res[size / 2], size / 2);
    }
    else if constexpr (radixbit_param == 3) { // Use renamed
        ButterflyAddAddMod(res, size);
        INTTradixButterfly<radixbit_param - 1>(&res[0], size / 2); // Use renamed
        const uint32_t block = size >> radixbit_param; // Use renamed
        for (int i = 0; i < block; i++) {
            res[2 * block + i + size / 2] =
                ConstTwiddleMul<radixbit_param, 2>(res[2 * block + i + size / 2]); // Use renamed
            const DoubleSWord temp = res[i + size / 2];
            res[i + size / 2] += res[2 * block + i + size / 2];
            res[2 * block + i + size / 2] =
                temp - res[2 * block + i + size / 2];
        }
        for (int i = 0; i < block; i++) {
            const DoubleSWord temp =
                ConstTwiddleMul<radixbit_param, 3>(res[1 * block + i + size / 2]); // Use renamed
            res[1 * block + i + size / 2] =
                ConstTwiddleMul<radixbit_param, 1>(res[1 * block + i + size / 2]) + // Use renamed
                ConstTwiddleMul<radixbit_param, 3>(res[3 * block + i + size / 2]); // Use renamed
            res[3 * block + i + size / 2] =
                temp +
                ConstTwiddleMul<radixbit_param, 1>(res[3 * block + i + size / 2]); // Use renamed
        }
        ButterflyAddBothSREDC(&res[size / 2], size / 4);
        ButterflyAddBothSREDC(&res[3 * size / 4], size / 4);
    }
}

export template <uint32_t Nbit, uint8_t radixbit_param> // Renamed radixbit
inline void INTTradix(DoubleSWord *const res, const uint32_t size,
                      const uint32_t num_block,
                      const std::array<std::array<SWord, 1 << Nbit>, 2> &table)
{
    INTTradixButterfly<radixbit_param>(res, size); // Use renamed
    for (uint32_t i = 1; i < (1 << radixbit_param); i++) // Use renamed
        TwiddleMulInvert<Nbit, radixbit_param>( // Use renamed
            &res[i * (size >> radixbit_param)], size >> radixbit_param, i, // Use renamed
            BitReverse<radixbit_param>(i) * num_block, table); // Use renamed
}

export template <uint32_t Nbit, uint8_t radixbit_param> // Renamed radixbit
inline void INTT(std::array<DoubleSWord, 1 << Nbit> &res,
                 const std::array<std::array<SWord, 1 << Nbit>, 2> &table)
{
    for (uint8_t sizebit = Nbit; sizebit > radixbit_param; sizebit -= radixbit_param) { // Use renamed
        const uint32_t size = 1U << sizebit;
        const uint32_t num_block = 1U << (Nbit - sizebit);
        for (uint32_t block = 0; block < num_block; block++)
            INTTradix<Nbit, radixbit_param>(&res[size * block], size, num_block, // Use renamed
                                      table);
    }
    constexpr uint8_t remainder = ((Nbit - 1) % radixbit_param) + 1; // Use renamed
    constexpr uint32_t size = 1U << remainder;
    constexpr uint32_t num_block = 1U << (Nbit - remainder);
    for (uint32_t block = 0; block < num_block; block++)
        INTTradixButterfly<remainder>(&res[size * block], size);
}

export template <typename T = uint32_t, uint Nbit_param = Nbit, bool modswitch> // Renamed Nbit
inline void TwistMulInvert(std::array<DoubleSWord, 1 << Nbit_param> &res, // Use renamed
                           const std::array<T, 1 << Nbit_param> &a, // Use renamed
                           const std::array<SWord, 1 << Nbit_param> &twist) // Use renamed
{
    constexpr uint N = 1 << Nbit_param; // Use renamed
    for (int i = 0; i < N; i++) {
        if constexpr (modswitch) {
            res[i] = (((static_cast<uint64_t>(a[i]) * K_const) << shiftamount) + // Use K_const
                      a[i] + (1ULL << (32 - 1))) >>
                     32;
        }
        else {
            res[i] = static_cast<DoubleSWord>(
                static_cast<typename std::make_signed<T>::type>(a[i]));
        }
        res[i] = MulSREDC(res[i], twist[i]);
    }
}

export template <typename T, uint32_t Nbit_param, bool modsiwtch = false> // Renamed Nbit
void TwistINTT(std::array<DoubleSWord, 1 << Nbit_param> &res, // Use renamed
               const std::array<T, 1 << Nbit_param> &a, // Use renamed
               const std::array<std::array<SWord, 1 << Nbit_param>, 2> &table, // Use renamed
               const std::array<SWord, 1 << Nbit_param> &twist) // Use renamed
{
    TwistMulInvert<T, Nbit_param, modsiwtch>(res, a, twist); // Use renamed
    INTT<Nbit_param, radixbit>(res, table); // Use renamed, radixbit is global const
}

export template <uint8_t radixbit_param, bool redc = true> // Renamed radixbit
inline void NTTradixButterfly(DoubleSWord *const res, const uint32_t size)
{
    static_assert(radixbit_param <= 3, "radix 8 is the maximum!"); // Use renamed
    if constexpr (radixbit_param == 1) { // Use renamed
        ButterflyAddBothMod(res, size);
    }
    else if constexpr (radixbit_param == 2) { // Use renamed
        NTTradixButterfly<radixbit_param - 1>(&res[0], size / 2); // Use renamed
        NTTradixButterfly<radixbit_param - 1>(&res[size / 2], size / 2); // Use renamed
        // const uint32_t block = size >> radixbit;
        for (uint index = 0; index < size / 4; index++) {
            const SWord temp = res[index];
            res[index] = AddMod(res[index], res[index + size / 2]);
            res[index + size / 2] = SubMod(temp, res[index + size / 2]);
        }
        for (uint index = size / 4; index < size / 2; index++) {
            const DoubleSWord temp = res[index];
            res[index + size / 2] =
                -ConstTwiddleMul<radixbit_param, 1>(res[index + size / 2]); // Use renamed
            if constexpr (redc) {
                res[index] = SREDC(res[index] + res[index + size / 2]);
                res[index + size / 2] = SREDC(temp - res[index + size / 2]);
            }
            else {
                res[index] = res[index] + res[index + size / 2];
                res[index + size / 2] = temp - res[index + size / 2];
            }
        }
    }
    else if constexpr (radixbit_param == 3) { // Use renamed
        NTTradixButterfly<radixbit_param - 1, false>(&res[0], size / 2); // Use renamed

        NTTradixButterfly<radixbit_param - 2>(&res[2 * size / 4], size / 4); // Use renamed
        NTTradixButterfly<radixbit_param - 2>(&res[3 * size / 4], size / 4); // Use renamed
        const uint32_t block = size >> radixbit_param; // Use renamed
        for (int index = size / 2; index < size / 2 + block; index++) {
            const DoubleSWord temp = res[index];
            res[index] = AddMod(res[index], res[index + size / 4]);
            res[index + size / 4] =
                -ConstTwiddleMul<radixbit_param, 2>(temp - res[index + size / 4]); // Use renamed
        }
        for (int index = size / 2 + block; index < size / 2 + 2 * block;
             index++) {
            const DoubleSWord temp = -ConstTwiddleMul<radixbit_param, 1>(res[index]); // Use renamed
            res[index] = -ConstTwiddleMul<radixbit_param, 3>(res[index]) - // Use renamed
                         ConstTwiddleMul<radixbit_param, 1>(res[index + size / 4]); // Use renamed
            res[index + size / 4] =
                temp - ConstTwiddleMul<radixbit_param, 3>(res[index + size / 4]); // Use renamed
        }
        for (int index = 0; index < block; index++) {
            const SWord temp = res[index];
            res[index] = AddMod(res[index], res[index + size / 2]);
            res[index + size / 2] = SubMod(temp, res[index + size / 2]);
        }
        for (int i = 1; i < 4; i++)
            for (int index = i * block; index < (i + 1) * block; index++) {
                const DoubleSWord temp = res[index];
                res[index] = SREDC(res[index] + res[index + size / 2]);
                res[index + size / 2] = SREDC(temp - res[index + size / 2]);
            }
    }
}

export template <uint Nbit_param, uint8_t prevradixbit> // Renamed Nbit
inline void TwiddleMul(DoubleSWord *const res, const uint sizebit,
                       const uint stride,
                       const std::array<std::array<SWord, 1 << Nbit_param>, 2> &table) // Use renamed
{
    const uint size = 1U << sizebit;
    if constexpr (prevradixbit == 1) {
        if (stride != 0)
            for (uint32_t index = 0; index < size; index++)
                res[index] = MulSREDC(res[index], table[0][stride * index]);
    }
    else {
        if (stride == 0) {
            for (uint32_t index = 0; index < size; index++)
                if (((index >> (sizebit - prevradixbit)) &
                     ((1 << (prevradixbit - 1)) - 1)) != 0) {
                    res[index] = MulSREDC(res[index], R2);
                }
        }
        else {
            for (uint32_t index = 0; index < size; index++)
                res[index] = MulSREDC(
                    res[index], table[((index >> (sizebit - prevradixbit)) &
                                       ((1 << (prevradixbit - 1)) - 1)) != 0]
                                     [stride * index]);
        }
    }
}

export template <uint Nbit_param, uint8_t radixbit_param, uint8_t prevradixbit> // Renamed Nbit, radixbit
inline void NTTradix(DoubleSWord *const res, const uint sizebit,
                     const uint32_t num_block,
                     const std::array<std::array<SWord, 1 << Nbit_param>, 2> &table) // Use renamed
{
    const uint size = 1U << sizebit;
    for (uint32_t i = 0; i < (1 << radixbit_param); i++) // Use renamed
        TwiddleMul<Nbit_param, prevradixbit>( // Use renamed
            &res[i * (size >> radixbit_param)], sizebit - radixbit_param, // Use renamed
            BitReverse<radixbit_param>(i) * num_block, table); // Use renamed
    NTTradixButterfly<radixbit_param>(res, size); // Use renamed
}

export template <uint32_t Nbit_param, uint8_t radixbit_param> // Renamed Nbit, radixbit
void NTT(std::array<DoubleSWord, 1 << Nbit_param> &res, // Use renamed
         const std::array<std::array<SWord, 1 << Nbit_param>, 2> &table) // Use renamed
{
    constexpr uint8_t remainder = ((Nbit_param - 1) % radixbit_param) + 1; // Use renamed
    {
        constexpr uint size = 1U << remainder;
        constexpr uint num_block = 1U << (Nbit_param - remainder); // Use renamed
        for (uint block = 0; block < num_block; block++)
            NTTradixButterfly<remainder>(&res[size * block], size);
    }
    {
        constexpr uint sizebit_val = remainder + radixbit_param; // Renamed sizebit
        const uint size = 1U << sizebit_val; // Use renamed
        const uint num_block = 1U << (Nbit_param - sizebit_val); // Use renamed
        for (uint block = 0; block < num_block; block++)
            NTTradix<Nbit_param, radixbit_param, remainder>(&res[size * block], sizebit_val, // Use renamed
                                                num_block, table);
    }
    for (uint8_t sizebit_val = remainder + 2 * radixbit_param; sizebit_val <= Nbit_param; // Renamed sizebit
         sizebit_val += radixbit_param) { // Use renamed
        const uint size = 1U << sizebit_val; // Use renamed
        const uint num_block = 1U << (Nbit_param - sizebit_val); // Use renamed
        for (uint block = 0; block < num_block; block++)
            NTTradix<Nbit_param, radixbit_param, radixbit_param>(&res[size * block], sizebit_val, // Use renamed
                                               num_block, table);
    }
}

export template <typename T = uint32_t, uint Nbit_param = Nbit, bool modswitch> // Renamed Nbit
inline void TwistMulDirect(std::array<T, 1 << Nbit_param> &res, // Use renamed
                           const std::array<DoubleSWord, 1 << Nbit_param> &a, // Use renamed
                           const std::array<SWord, 1 << Nbit_param> &twist) // Use renamed
{
    constexpr uint32_t N = 1 << Nbit_param; // Use renamed
    for (int i = 0; i < N; i++) {
        const SWord mulres = MulSREDC(a[i], twist[i]);
        if constexpr (modswitch)
            res[i] =
                (static_cast<uint64_t>((mulres < 0) ? mulres + P : mulres) *
                     ((1ULL << (32 + wordbits - 1)) / P) +
                 (1ULL << (wordbits - 1 - 1))) >>
                (wordbits - 1);
        else
            res[i] = (mulres < 0) ? mulres + P : mulres;
    }
}

export template <typename T, uint32_t Nbit_param, bool modsiwtch = false> // Renamed Nbit
void TwistNTT(std::array<T, 1 << Nbit_param> &res, // Use renamed
              std::array<DoubleSWord, 1 << Nbit_param> &a, // Use renamed
              const std::array<std::array<SWord, 1 << Nbit_param>, 2> &table, // Use renamed
              const std::array<SWord, 1 << Nbit_param> &twist) // Use renamed
{
    NTT<Nbit_param, radixbit>(a, table); // Use renamed, radixbit is global const
    TwistMulDirect<T, Nbit_param, modsiwtch>(res, a, twist); // Use renamed
}

export template <typename T, uint32_t Nbit_param, bool modswitcha, bool modswitchb> // Renamed Nbit
void PolyMullvl1(
    std::array<T, 1 << Nbit_param> &res, std::array<T, 1 << Nbit_param> &a, // Use renamed
    std::array<T, 1 << Nbit_param> &b, // Use renamed
    const std::array<std::array<std::array<SWord, 1 << Nbit_param>, 2>, 2> &table, // Use renamed
    const std::array<std::array<SWord, 1 << Nbit_param>, 2> &twist) // Use renamed
{
    constexpr uint8_t remainder = ((Nbit_param - 1) % 3) + 1; // Use renamed
    std::array<DoubleSWord, 1 << Nbit_param> ntta, nttb; // Use renamed
    TwistINTT<T, Nbit_param, modswitcha>(ntta, a, table[1], twist[1]); // Use renamed
    TwistINTT<T, Nbit_param, modswitchb>(nttb, b, table[1], twist[1]); // Use renamed
    for (int i = 0; i < (1U << Nbit_param); i++) // Use renamed
        if ((i & ((1 << remainder) - 1)) > 1)
            ntta[i] = MulSREDC(MulSREDC(nttb[i], R4), ntta[i]);
        else
            ntta[i] = MulSREDC(MulSREDC(nttb[i], R2), ntta[i]);
    TwistNTT<T, Nbit_param, modswitcha | modswitchb>(res, ntta, table[0], twist[0]); // Use renamed
}

}  // namespace raintt