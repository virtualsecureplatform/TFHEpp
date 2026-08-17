#pragma once

// Shared high-level scaffolding.  The active parameter-family header still
// owns lvl0--lvl4 and its key-switching parameters.

struct lvl5param {
    static constexpr int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = 14, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t lₐ = 5, l = 5, Bgbit = 19, Bgₐbit = 19;
    using T = MultiLimbUInt<7>;
    static constexpr T Bg = T{1} << Bgbit, Bgₐ = T{1} << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -425);
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint64_t plain_modulus_u64 = 786433;
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{plain_modulus_u64};
    static constexpr double Δ = 0.0;
    static constexpr T delta_int =
        std::numeric_limits<T>::max() / plain_modulus_u64;
    static constexpr uint64_t Q_mod_t =
        (std::numeric_limits<T>::max() % plain_modulus_u64) + 1;
    static constexpr uint64_t bfv_bootstrap_digit_error_bound = 15;
    static constexpr int bfv_bootstrap_linear_bsgs_step = 128;
    static constexpr std::uint32_t l̅ = 28, l̅ₐ = 28, B̅gbit = 16, B̅gₐbit = 16;
    static constexpr uint64_t simd_modulus = plain_modulus_u64,
                              simd_psi = 585160;
    static constexpr uint64_t simd_psi_inv = 253771, simd_n_inv = 786385;
};

struct lvl5bootparam {
    static constexpr int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = 15, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t lₐ = 35, l = 35, Bgbit = 18, Bgₐbit = 18;
    using T = MultiLimbUInt<10>;
    static constexpr T Bg = T{1} << Bgbit, Bgₐ = T{1} << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -607);
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint64_t plain_modulus_u64 = 786433;
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{plain_modulus_u64};
    static constexpr double Δ = 0.0;
    static constexpr T delta_int =
        std::numeric_limits<T>::max() / plain_modulus_u64;
    static constexpr uint64_t Q_mod_t =
        (std::numeric_limits<T>::max() % plain_modulus_u64) + 1;
    static constexpr uint64_t bfv_bootstrap_digit_error_bound = 23;
    static constexpr int bfv_bootstrap_linear_bsgs_step = 128;
    static constexpr uint32_t bfv_key_hamming_weight = 96;
    static constexpr std::uint32_t l̅ = 40, l̅ₐ = 40, B̅gbit = 16, B̅gₐbit = 16;
    static constexpr uint64_t simd_modulus = plain_modulus_u64,
                              simd_psi = 108788;
    static constexpr uint64_t simd_psi_inv = 295516, simd_n_inv = 786409;
};

struct lvl6param {
    static constexpr int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = 15, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t lₐ = 5, l = 5, Bgbit = 16, Bgₐbit = 16;
    using T = MultiLimbUInt<14>;
    static constexpr T Bg = T{1} << Bgbit, Bgₐ = T{1} << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -872);
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << (std::numeric_limits<T>::digits - 3);
    static constexpr uint64_t plain_modulus_u64 = 786433;
    static constexpr uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{plain_modulus_u64};
    static constexpr double Δ = 0.0;
    static constexpr T delta_int =
        std::numeric_limits<T>::max() / plain_modulus_u64;
    static constexpr uint64_t Q_mod_t =
        (std::numeric_limits<T>::max() % plain_modulus_u64) + 1;
    static constexpr uint64_t bfv_bootstrap_digit_error_bound = 15;
    static constexpr int bfv_bootstrap_linear_bsgs_step = 128;
    static constexpr std::uint32_t l̅ = 128, l̅ₐ = 128, B̅gbit = 7, B̅gₐbit = 7;
    static constexpr uint64_t simd_modulus = plain_modulus_u64,
                              simd_psi = 108788;
    static constexpr uint64_t simd_psi_inv = 295516, simd_n_inv = 786409;
};
