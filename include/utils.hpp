#pragma once

#ifdef USE_RANDEN
#include <randen.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <random>

namespace TFHEpp {
#ifdef USE_RANDEN
static thread_local std::random_device trng;
static thread_local randen::Randen<uint64_t> generator(trng);
#else
static thread_local std::random_device generator;
#endif

// https://qiita.com/saka1_p/items/e8c4dfdbfa88449190c5
template <typename T>
constexpr bool false_v = false;

// https://qiita.com/negi-drums/items/a527c05050781a5af523
template <typename T>
concept hasq = requires { T::q; };

template <typename T>
concept hasqbit = requires { T::qbit; };

// https://github.com/zhourrr/aligned-memory-allocator/blob/main/aligned_allocator.h
// A minimal implementation of an allocator for C++ Standard Library, which
// allocates aligned memory (specified by the alignment argument).
// Note:
//    A minimal custom allocator is preferred because C++ allocator_traits class
//    provides default implementation for you. Take a look at Microsoft's
//    documentation about Allocators and allocator class.
template <typename T, std::size_t alignment>
class AlignedAllocator {
public:
    using value_type = T;

public:
    // According to Microsoft's documentation, default constructor is not
    // required by C++ Standard Library.
    AlignedAllocator() noexcept {};

    template <typename U>
    AlignedAllocator(const AlignedAllocator<U, alignment> &other) noexcept
    {
    }

    template <typename U>
    inline bool operator==(
        const AlignedAllocator<U, alignment> &other) const noexcept
    {
        return true;
    }

    template <typename U>
    inline bool operator!=(
        const AlignedAllocator<U, alignment> &other) const noexcept
    {
        return false;
    }

    template <typename U>
    struct rebind {
        using other = AlignedAllocator<U, alignment>;
    };

    // STL containers call this function to allocate uninitialized memory block
    // to store (no more than n) elements of type T (value_type).
    inline value_type *allocate(const std::size_t n) const
    {
        auto size = n;
        /*
          If you wish, for some strange reason, that the size of allocated
          buffer is also aligned to alignment, uncomment the following
          statement.

          Note: this increases the size of underlying memory, but STL containers
          still treat it as a memory block of size n, i.e., STL containers will
          not put more than n elements into the returned memory.
        */
        // size = (n + alignment - 1) / alignment * alignment;
        return static_cast<value_type *>(
            std::aligned_alloc(alignment, sizeof(T) * size));
    };

    // STL containers call this function to free a memory block beginning at a
    // specified position.
    inline void deallocate(value_type *const ptr, std::size_t n) const noexcept
    {
        std::free(ptr);
    }
};

// Double to Torus(32bit fixed-point number)
inline uint16_t dtot16(double d)
{
    return int16_t(int64_t((d - int64_t(d)) * (1LL << 16)));
}

// Double to Torus(32bit fixed-point number)
inline uint32_t dtot32(double d)
{
    return int32_t(int64_t((d - int64_t(d)) * (1LL << 32)));
}

// Modular Gaussian Distribution over Torus
template <class P>
inline typename P::T ModularGaussian(typename P::T center, double stdev)
{
    if constexpr (std::is_same_v<typename P::T, uint16_t>) {
        // 16bit fixed-point number version
        std::normal_distribution<double> distribution(0., stdev);
        double err = distribution(generator);
        return center + dtot16(err);
    }
    else if constexpr (std::is_same_v<typename P::T, uint32_t>) {
        // 32bit fixed-point number version
        std::normal_distribution<double> distribution(0., stdev);
        double err = distribution(generator);
        return center + dtot32(err);
    }
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        // 64bit fixed-point number version
        static const double _2p64 = std::pow(2., 64);
        std::normal_distribution<double> distribution(0., 1.0);
        const double val = stdev * distribution(generator) * _2p64;
        const uint64_t ival = static_cast<typename P::T>(val);
        return ival + center;
    }
    else
        static_assert(false_v<typename P::T>, "Undefined Modular Gaussian!");
}

template <class P>
inline typename P::T CenteredBinomial(uint η)
{
    static uint64_t buf;
    static uint_fast8_t count = 64;
    typename P::T acc = 0;
    std::uniform_int_distribution<uint64_t> dist(
        0, std::numeric_limits<uint64_t>::max());
    for (int i = 0; i < 2 * η; i++) {
        if (count >= 64) {
            buf = dist(generator);
        }
        acc += buf & 1;
        buf >>= 1;
        count++;
    }
    return acc - η;
}

// https://stackoverflow.com/questions/21191307/minimum-number-of-bits-to-represent-a-given-int
template <uint32_t data>
constexpr int bits_needed()
{
    uint32_t value = data;
    int bits = 0;
    for (int bit_test = 16; bit_test > 0; bit_test >>= 1) {
        if (value >> bit_test != 0) {
            bits += bit_test;
            value >>= bit_test;
        }
    }
    return bits + value;
}

template <class P>
inline void PolynomialMulByXai(Polynomial<P> &res, const Polynomial<P> &poly,
                               const typename P::T a)
{
    if (a == 0)
        res = poly;
    else if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a];
    }
    else {
        const typename P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa];
    }
}

template <class P>
inline void PolynomialMulByXaiMinusOne(Polynomial<P> &res,
                                       const Polynomial<P> &poly,
                                       const typename P::T a)
{
    if (a < P::n) {
        for (int i = 0; i < a; i++) res[i] = -poly[i - a + P::n] - poly[i];
        for (int i = a; i < P::n; i++) res[i] = poly[i - a] - poly[i];
    }
    else {
        const typename P::T aa = a - P::n;
        for (int i = 0; i < aa; i++) res[i] = poly[i - aa + P::n] - poly[i];
        for (int i = aa; i < P::n; i++) res[i] = -poly[i - aa] - poly[i];
    }
}

// calcurate τ_d
template <class P>
inline void Automorphism(Polynomial<P> &res, const Polynomial<P> &poly,
                         const uint d)
{
    res = {};
    constexpr uint Nmask = (1ULL << (P::nbit)) - 1;
    constexpr uint signmask = 1ULL << (P::nbit);
    for (uint i = 0; i < P::n; i++) {
        const uint index = i * d;
        if (index & signmask)
            res[index & Nmask] -= poly[i];
        else
            res[index & Nmask] += poly[i];
    }
}

}  // namespace TFHEpp
