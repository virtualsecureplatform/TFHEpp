#pragma once
#include <cstdint>
#include <array>
#include <span>

namespace FNTpp{
    constexpr unsigned int Kbit = 5;
    constexpr unsigned int K = 1 << Kbit;
    constexpr int64_t P = (1ULL << K) + 1;
    constexpr int64_t wordmask = (1ULL << K) - 1;

    template <uint8_t bit>
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

    static inline int64_t ModLshift(int64_t a, uint8_t b)
    {
        // If b >= 32, multiply by 2^32 ≡ -1 (mod P).
        // => a = P - a  (unless a == 0), then reduce b by 32.
        if (b >= 32) {
            if (a != 0) {
                a = P - a; 
            }
            b -= 32; // now b < 32
        }

        // Shift by b < 32 in 64-bit arithmetic (safe from overflow).
        int64_t r = a << b;

        // Now reduce a modulo P:
        //   hi = upper 32 bits
        //   lo = lower 32 bits
        // Since (hi << 32) ≡ -hi (mod P),
        // we can do (lo + hi) mod P  and then subtract P if needed.
        const int64_t hi = r >> K;
        const int64_t lo = r & wordmask;
        r = -hi + lo;

        // Subtract P once or twice if needed to ensure a < P
        if (r < 0) r += P;
        if (r >= P) r -= P;
        return r;
    }

    template <uint8_t Nbit>
    inline void MulInvN(std::array<int64_t, 1u<<Nbit>& a){
        for(int i = 0; i < (1u<<Nbit); i++) a[i] = ModLshift(a[i], 2*K-Nbit);
    }
    
    template <unsigned int Nbit>
    void FNT(const std::span<int64_t, 1u << Nbit> res)
    {
        if constexpr (Nbit == 1){
            const int64_t temp = res[0];
            res[0] += res[1];
            if(res[0] >= P) res[0] -= P;
            res[1] = temp - res[1];
            if(res[1] < 0) res[1] += P;
        }else{
            constexpr unsigned int N = 1u << Nbit;
            constexpr unsigned int stride = 1u << (Kbit+1 - Nbit);
            for(unsigned int i = 0; i < N/2; i++){
                const int64_t temp = res[i]+res[i+N/2];
                res[i+N/2] = res[i]-res[i+N/2];
                if(res[i+N/2] < 0) res[i+N/2] += P;
                if(i!=0) res[i+N/2] = ModLshift(res[i+N/2],i*stride);
                res[i] = temp >= P ? temp - P : temp;
            }
            FNT<Nbit-1>(res.template subspan<0,N/2>());
            FNT<Nbit-1>(res.template subspan<N/2,N/2>());
        }
    }

    template <unsigned int Nbit>
    void IFNT(const std::span<int64_t, 1u << Nbit> res)
    {
        if constexpr (Nbit == 1){
            const int64_t temp = res[0];
            res[0] += res[1];
            if(res[0] >= P) res[0] -= P;
            res[1] = temp - res[1];
            if(res[1] < 0) res[1] += P;
        }else{
            constexpr unsigned int N = 1u << Nbit;
            IFNT<Nbit-1>(res.template subspan<0,N/2>());
            IFNT<Nbit-1>(res.template subspan<N/2,N/2>());
            constexpr unsigned int stride = 1u << (Kbit+1 - Nbit);
            for(unsigned int i = 0; i < N/2; i++){
                if(i!=0) res[i+N/2] = ModLshift(res[i+N/2],(N-i)*stride);
                const int64_t temp = res[i]+res[i+N/2];
                res[i+N/2] = res[i]-res[i+N/2]; //Part of twiddle factor
                if(res[i+N/2] < 0) res[i+N/2] += P;
                res[i] = temp >= P ? temp - P : temp;
            }
        }
    }


    template <unsigned int Nbit>
    void TwistFNT(std::array<int64_t, 1u << (Nbit+1)> &res, const std::array<int64_t, 1u << Nbit> &a)
    {
        constexpr unsigned int formersizebit = (Nbit + 1)/2;
        static_assert(formersizebit <= Kbit, "sizebit must be less than or equal to Kbit");
        constexpr unsigned int formersize = 1u << formersizebit;
        constexpr unsigned int latersizebit = (Nbit + 1) - formersizebit;
        constexpr unsigned int latersize = 1u << latersizebit;
        constexpr unsigned int formerrbit = (Kbit+1) - formersizebit;
        constexpr unsigned int laterrbit = (Kbit+1) - latersizebit;
        //Former
        for(unsigned int i = 0; i < latersize/2; i++){
            std::array<int64_t, formersize> temp;
            for(unsigned int j = 0; j < formersize; j++)
                temp[j] = ModLshift(a[j*(latersize/2) + i],(j*(latersize/2) + i)<<(formerrbit-1));
            FNT<formersizebit>(std::span{temp});
            for(unsigned int j = 0; j < formersize; j++)
                res[j*latersize + i] = temp[j];
        }
        //Later
        for(unsigned int i = 0; i < formersize; i++)
            FNT<latersizebit>(std::span{res}.subspan(i,latersize));
    }
}