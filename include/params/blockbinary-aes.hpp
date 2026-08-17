#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

// AES and ASCON use a packed, high-precision circuit bootstrap.  Keep that
// target topology separate from the default block-binary lvl2 (k = 2).
struct blockbinaryaeslvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = 11;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t lₐ = l;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = Bgbit;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -51);
    using T = uint64_t;
    static constexpr std::make_signed_t<T> μ = 1LL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << (std::numeric_limits<T>::digits - 4));
    static constexpr std::uint32_t l̅ = 1;
    static constexpr std::uint32_t l̅ₐ = l̅;
    static constexpr std::uint32_t B̅gbit = std::numeric_limits<T>::digits;
    static constexpr std::uint32_t B̅gₐbit = B̅gbit;
};

struct blockbinaryaesAHlvl2param {
    using baseP = blockbinaryaeslvl2param;
    static constexpr int32_t key_value_max = baseP::key_value_max;
    static constexpr int32_t key_value_min = baseP::key_value_min;
    static constexpr std::uint32_t nbit = baseP::nbit;
    static constexpr std::uint32_t n = baseP::n;
    static constexpr std::uint32_t k = baseP::k;
    static constexpr std::uint32_t lₐ = 5;
    static constexpr std::uint32_t l = 5;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr ErrorDistribution errordist = baseP::errordist;
    static const inline double α = baseP::α;
    using T = typename baseP::T;
    static constexpr std::make_signed_t<T> μ = baseP::μ;
    static constexpr uint32_t plain_modulus = baseP::plain_modulus;
    static constexpr double Δ = baseP::Δ;
    static constexpr std::uint32_t l̅ = baseP::l̅;
    static constexpr std::uint32_t l̅ₐ = baseP::l̅ₐ;
    static constexpr std::uint32_t B̅gbit = baseP::B̅gbit;
    static constexpr std::uint32_t B̅gₐbit = baseP::B̅gₐbit;
};
