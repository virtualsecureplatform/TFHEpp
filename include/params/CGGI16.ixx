export module tfhepp:params:CGGI16;
import std;
// Assuming ErrorDistribution is made available by tfhepp:params main module interface file.
// import tfhepp:params; // Only if ErrorDistribution is not visible otherwise.

export struct lvl0param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 500;
    static constexpr std::uint32_t k = 1;
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = 2.44e-5;
    using T = uint32_t;
    static constexpr T μ = 1U << 29;
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

// Dummy
export struct lvlhalfparam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr int32_t key_value_diff = key_value_max - key_value_min;
    static constexpr std::uint32_t n = 760;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -17);  // fresh noise
    using T = uint32_t;                                 // Torus representation
    static constexpr T μ = 1U << (std::numeric_limits<T>::digits - 3);
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

export struct lvl1param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr std::uint32_t nbit = 10;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 2;
    static constexpr std::uint32_t l = 2;
    static constexpr std::uint32_t Bgbit = 10;
    static constexpr std::uint32_t Bgₐbit = 10;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = 3.73e-9;
    using T = uint32_t;
    static constexpr T μ = 1U << 29;
    static constexpr uint32_t plain_modulus = 2;
    static constexpr double Δ =
        static_cast<double>(1ULL << std::numeric_limits<T>::digits) /
        plain_modulus;
};

export struct lvl2param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static const std::uint32_t nbit = 11;
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -44);
    using T = uint64_t;
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ = μ;
};

// Dummy
export struct lvl3param {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static const std::uint32_t nbit = 13;  // dimension must be a power of 2 for
    // ease of polynomial multiplication.
    static constexpr std::uint32_t n = 1 << nbit;  // dimension
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 9;
    static constexpr std::uint32_t Bgₐbit = 9;
    static constexpr std::uint32_t Bg = 1 << Bgbit;
    static constexpr std::uint32_t Bgₐ = 1 << Bgₐbit;
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -47);  // fresh noise
    using T = uint64_t;                                 // Torus representation
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulusbit = 31;
    static constexpr uint64_t plain_modulus = 1ULL << plain_modulusbit;
    static constexpr double Δ = 1ULL << (64 - plain_modulusbit - 1);
};

export struct lvl10param {
    static constexpr std::uint32_t t = 8;
    static constexpr std::uint32_t basebit = 2;
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;
    using domainP = lvl1param;
    using targetP = lvl0param;
};

// Dummy
export struct lvl1hparam {
    static constexpr std::uint32_t t =
        10;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        3;  // how many bit should be encrypted in keyswitching key
    static const inline double α = lvlhalfparam::α;  // key noise
    using domainP = lvl1param;
    using targetP = lvlhalfparam;
};

// Dummy
export struct lvl11param {
    static constexpr std::uint32_t t = 0;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        0;  // how many bit should be encrypted in keyswitching key
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl1param;
    using targetP = lvl1param;
};

export struct lvl21param {
    static constexpr std::uint32_t t = 10;
    static constexpr std::uint32_t basebit = 3;
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2, -31);
    using domainP = lvl2param;
    using targetP = lvl1param;
};

// Dummy
export struct lvl20param {
    static constexpr std::uint32_t t = 0;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        0;  // how many bit should be encrypted in keyswitching key
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = lvl0param::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvl0param;
};

// Dummy
export struct lvl2hparam {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = lvlhalfparam::α;  // key noise
    using domainP = lvl2param;
    using targetP = lvlhalfparam;
};

// Dummy
export struct lvl22param {
    static constexpr std::uint32_t t = 0;
    static constexpr std::uint32_t basebit = 0;
    static constexpr TFHEpp::ErrorDistribution errordist = // Assuming ErrorDistribution is in TFHEpp namespace
        TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = lvl2param::α;
    using domainP = lvl2param;
    using targetP = lvl2param;
};

export struct lvl31param {
    static constexpr std::uint32_t t = 7;  // number of addition in keyswitching
    static constexpr std::uint32_t basebit =
        2;  // how many bit should be encrypted in keyswitching key
    static const inline double α = lvl1param::α;  // key noise
    using domainP = lvl3param;
    using targetP = lvl1param;
};