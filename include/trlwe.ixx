export module tfhepp:trlwe;
import std;

export import tfhepp:key;
export import tfhepp:mulfft;
export import tfhepp:params;
export import tfhepp:raintt; // For raintt::namespace
// utils.hpp was not directly imported, but TFHEpp::generator (from utils.ixx) is used.
// and CenteredBinomial, ModularGaussian
import tfhepp:utils;


export namespace TFHEpp {
export template <class P>
TRLWE<P> trlweSymEncryptZero(const double α, const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
    TRLWE<P> c;
    for (typename P::T &i : c[P::k]) i = ModularGaussian<P>(0, α);
    for (int k = 0; k < P::k; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator);
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < P::n; i++) c[P::k][i] += temp[i];
    }
    return c;
}

export template <class P>
TRLWE<P> trlweSymEncryptZero(const uint η, const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(
        0, std::numeric_limits<typename P::T>::max());
    alignas(64) TRLWE<P> c;
    for (typename P::T &i : c[P::k])
        i = (CenteredBinomial<P>(η) << std::numeric_limits<P>::digits) / P::q;
    for (int k = 0; k < P::k; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator);
        alignas(64) std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        alignas(64) Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < P::n; i++) c[P::k][i] += temp[i];
    }
    return c;
}

export template <class P>
TRLWE<P> trlweSymEncryptZero(const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymEncryptZero<P>(P::α, key);
    else
        return trlweSymEncryptZero<P>(P::η, key);
}

export template <class P>
TRLWERAINTT<P> trlwerainttSymEncryptZero(const uint η, const Key<P> &key)
{
    static_assert(P::q == TFHEpp::raintt::P); // TFHEpp::raintt
    static_assert(P::qbit == TFHEpp::raintt::wordbits); // TFHEpp::raintt
    std::uniform_int_distribution<typename P::T> Torusdist(0, P::q - 1);
    constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
    TRLWERAINTT<P> c = {};
    {
        Polynomial<P> b;
        for (typename P::T &i : b) i = CenteredBinomial<P>(η);
        raintt::TwistINTT<typename P::T, P::nbit, false>(
            c[P::k], b, (*raintttable)[1], (*raintttwist)[1]);
        for (int i = 0; i < P::n; i++)
            if ((i & ((1 << remainder) - 1)) > 1)
                c[P::k][i] = raintt::MulSREDC(c[P::k][i], raintt::R2);
    }
    for (int k = 0; k < P::k; k++) {
        for (typename raintt::DoubleSWord &i : c[k]) i = Torusdist(generator);
        PolynomialRAINTT<P> partkeyraintt;
        {
            Polynomial<P> partkey;
            for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
            raintt::TwistINTT<typename P::T, P::nbit, false>(
                partkeyraintt, partkey, (*raintttable)[1], (*raintttwist)[1]);
            for (int i = 0; i < P::n; i++)
                if ((i & ((1 << remainder) - 1)) > 1)
                    partkeyraintt[i] =
                        raintt::MulSREDC(partkeyraintt[i], raintt::R3);
                else
                    partkeyraintt[i] =
                        raintt::MulSREDC(partkeyraintt[i], raintt::R2);
        }
        for (int i = 0; i < P::n; i++)
        c[P::k][i] = TFHEpp::raintt::AddMod( // TFHEpp::raintt
            c[P::k][i], TFHEpp::raintt::MulSREDC(c[k][i], partkeyraintt[i])); // TFHEpp::raintt
    }
    return c;
}

export template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p,
                         const double α, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
    return c;
}

export template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p, const uint η,
                         const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(η, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
    return c;
}

export template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p,
                         const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymEncrypt<P>(p, P::α, key);
    else
        return trlweSymEncrypt<P>(p, P::η, key);
}

export template <class P, bool modswitch = false>
TRLWERAINTT<P> trlwerainttSymEncrypt(const Polynomial<P> &p, const uint η,
                                     const Key<P> &key)
{
    TRLWERAINTT<P> c = trlwerainttSymEncryptZero<P>(η, key);
    PolynomialRAINTT<P> pntt;
    TFHEpp::raintt::TwistINTT<typename P::T, P::nbit, modswitch>( // TFHEpp::raintt
        pntt, p, (*TFHEpp::raintttable)[1], (*TFHEpp::raintttwist)[1]); // TFHEpp::
    constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
    for (int i = 0; i < P::n; i++)
        if ((i & ((1 << remainder) - 1)) > 1)
            pntt[i] = TFHEpp::raintt::MulSREDC(pntt[i], TFHEpp::raintt::R2); // TFHEpp::raintt
    for (int i = 0; i < P::n; i++)
        c[P::k][i] = TFHEpp::raintt::AddMod(pntt[i], c[P::k][i]); // TFHEpp::raintt
    return c;
}

export template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const double α, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::Δ * p[i]);
    return c;
}

export template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const uint η, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(η, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::Δ * p[i]);
    return c;
}

export template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymIntEncrypt<P>(p, P::α, key);
    else
        return trlweSymIntEncrypt<P>(p, P::η, key);
}

export template <class P>
Polynomial<P> trlwePhase(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }
    return phase;
}

export template <class P>
std::array<bool, P::n> trlweSymDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = trlwePhase<P>(c, key);

    std::array<bool, P::n> p;
    if constexpr (hasq<P>) {
        for (int i = 0; i < P::n; i++) p[i] = (phase[i] % P::q) < P::q / 2;
    }
    else
        for (int i = 0; i < P::n; i++)
            p[i] = static_cast<typename std::make_signed<typename P::T>::type>(
                       phase[i]) > 0;
    return p;
}

export template <class P, uint plain_modulus = P::plain_modulus>
Polynomial<P> trlweSymIntDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }

    constexpr double Δ =
        std::pow(2.0, std::numeric_limits<typename P::T>::digits) /
        plain_modulus;
    Polynomial<P> p;
    for (int i = 0; i < P::n; i++)
        p[i] = static_cast<typename P::T>(std::round(phase[i] / Δ)) %
               plain_modulus;
    return p;
}

export template <class P, uint plain_modulus = P::plain_modulus>
Polynomial<P> trlweSymIntDecrypt(const TRLWE<P> &c, const SecretKey &sk)
{
    return trlweSymIntDecrypt<P, plain_modulus>(c, sk.key.get<P>());
}

export template <class P>
void SampleExtractIndex(TLWE<P> &tlwe, const TRLWE<P> &trlwe_in, const int index) // Renamed trlwe to trlwe_in
{
    for (int k = 0; k < P::k; k++) {
        for (int i = 0; i <= index; i++)
            tlwe[k * P::n + i] = trlwe_in[k][index - i]; // Use trlwe_in
        for (int i = index + 1; i < P::n; i++)
            tlwe[k * P::n + i] = -trlwe_in[k][P::n + index - i]; // Use trlwe_in
    }
    tlwe[P::k * P::n] = trlwe_in[P::k][index]; // Use trlwe_in
}

export template <class P>
void InvSampleExtractIndex(TRLWE<P> &trlwe, const TLWE<P> &tlwe_in, // Renamed tlwe to tlwe_in
                           const int index)
{
    for (int k = 0; k < P::k; k++) {
        for (int i = 0; i <= index; i++)
            trlwe[k][index - i] = tlwe_in[k * P::n + i]; // Use tlwe_in
        for (int i = index + 1; i < P::n; i++)
            trlwe[k][P::n + index - i] = -tlwe_in[k * P::n + i]; // Use tlwe_in
    }
    trlwe[P::k] = {};
    trlwe[P::k][index] = tlwe_in[P::k * P::n]; // Use tlwe_in
}

}  // namespace TFHEpp