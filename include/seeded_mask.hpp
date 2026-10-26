#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string_view>
#ifdef USE_BLAKE3
extern "C" {
#include <blake3.h>
}
#endif

namespace TFHEpp {
// Portable, domain-separated keyed BLAKE3 XOF. The complete 256-bit seed is
// the key; output words are explicitly little endian on every architecture.
class SeededMaskStream {
#ifdef USE_BLAKE3
    blake3_hasher hasher_;
    std::array<std::uint8_t,4096> buffer_{};
    std::size_t next_=buffer_.size();
    std::uint64_t position_=0;
#endif
public:
    SeededMaskStream(const std::array<std::uint64_t,4>&seed,
                     std::string_view domain)
    {
#ifdef USE_BLAKE3
        std::array<std::uint8_t,32> bytes{};
        for(std::size_t i=0;i<seed.size();++i)
            for(std::size_t j=0;j<8;++j)bytes[8*i+j]=std::uint8_t(seed[i]>>(8*j));
        blake3_hasher_init_keyed(&hasher_,bytes.data());
        blake3_hasher_update(&hasher_,domain.data(),domain.size());
#else
        (void)seed;(void)domain;
        throw std::runtime_error("seeded cryptographic masks require USE_BLAKE3");
#endif
    }
    std::uint64_t next()
    {
#ifdef USE_BLAKE3
        if(next_==buffer_.size()){
            if(position_>std::numeric_limits<std::uint64_t>::max()-buffer_.size())
                throw std::overflow_error("seeded mask stream exhausted");
            blake3_hasher_finalize_seek(&hasher_,position_,buffer_.data(),buffer_.size());
            position_+=buffer_.size();next_=0;
        }
        std::uint64_t result=0;
        for(std::size_t j=0;j<8;++j)result|=std::uint64_t(buffer_[next_++])<<(8*j);
        return result;
#else
        throw std::runtime_error("seeded cryptographic masks require USE_BLAKE3");
#endif
    }
};
} // namespace TFHEpp
