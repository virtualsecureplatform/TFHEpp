#include <array>
#include <limits>
#include <mulfft.hpp>
#include <params.hpp>
#include <trlwe.hpp>
#include <utils.hpp>

namespace TFHEpp {

inline TRGSWlvl1 trgswSymEncryptlvl1(int32_t p, double α, Keylvl1 &key)
{
    array<uint32_t, DEF_l> h;
    for (int i = 0; i < DEF_l; i++) h[i] = 1U << (32 - (i + 1) * DEF_Bgbit);
    TRGSWlvl1 trgsw;
    for (TRLWElvl1 &trlwe : trgsw) trlwe = trlweSymEncryptZerolvl1(α, key);
    for (int i = 0; i < DEF_l; i++) {
        trgsw[i][0][0] += static_cast<uint32_t>(p) * h[i];
        trgsw[i + DEF_l][1][0] += static_cast<uint32_t>(p) * h[i];
    }
    return trgsw;
}

inline TRGSWlvl2 trgswSymEncryptlvl2(int64_t p, double α, Keylvl2 &key)
{
    array<uint64_t, DEF_lbar> h;
    for (int i = 0; i < DEF_lbar; i++)
        h[i] = 1UL << (64 - (i + 1) * DEF_Bgbitbar);
    TRGSWlvl2 trgsw;
    for (TRLWElvl2 &trlwe : trgsw) trlwe = trlweSymEncryptZerolvl2(α, key);
    for (int i = 0; i < DEF_lbar; i++) {
        trgsw[i][0][0] += static_cast<uint64_t>(p) * h[i];
        trgsw[i + DEF_lbar][1][0] += static_cast<uint64_t>(p) * h[i];
    }
    return trgsw;
}

TRGSWFFTlvl1 trgswfftSymEncryptlvl1(int32_t p, double α, Keylvl1 &key)
{
    TRGSWlvl1 trgsw = trgswSymEncryptlvl1(p, α, key);
    TRGSWFFTlvl1 trgswfft;
    for (int i = 0; i < 2 * DEF_l; i++)
        for (int j = 0; j < 2; j++) TwistIFFTlvl1(trgswfft[i][j], trgsw[i][j]);
    return trgswfft;
}

TRGSWFFTlvl2 trgswfftSymEncryptlvl2(int64_t p, double α, Keylvl2 &key)
{
    TRGSWlvl2 trgsw = trgswSymEncryptlvl2(p, α, key);
    TRGSWFFTlvl2 trgswfft;
    for (int i = 0; i < 2 * DEF_lbar; i++)
        for (int j = 0; j < 2; j++) TwistIFFTlvl2(trgswfft[i][j], trgsw[i][j]);
    return trgswfft;
}

template <typename T = uint32_t, uint32_t N = DEF_N, uint32_t l = DEF_l,
          uint32_t Bgbit = DEF_Bgbit>
inline void Decomposition(array<array<T, N>, 2 * l> &decvec,
                          const array<array<T, N>, 2> &trlwe, const T offset)
{
    const T mask = static_cast<T>((1UL << Bgbit) - 1);
    for (int i = 0; i < N; i++) {
        decvec[0][i] = trlwe[0][i] + offset;
        decvec[l][i] = trlwe[1][i] + offset;
    }

    const T halfBg = (1UL << (Bgbit - 1));
    for (int i = l - 1; i >= 0; i--) {
        for (int j = 0; j < N; j++) {
            decvec[i][j] = ((decvec[0][j] >>
                             (numeric_limits<T>::digits - (i + 1) * Bgbit)) &
                            mask) -
                           halfBg;
            decvec[i + l][j] = ((decvec[l][j] >> (numeric_limits<T>::digits -
                                                  (i + 1) * Bgbit)) &
                                mask) -
                               halfBg;
        }
    }
}

constexpr uint32_t offsetgenlvl1()
{
    uint32_t offset = 0;
    for (int i = 1; i <= DEF_l; i++)
        offset += DEF_Bg / 2 * (1U << (32 - i * DEF_Bgbit));
    return offset;
}

constexpr uint64_t offsetgenlvl2()
{
    uint64_t offset = 0;
    for (int i = 1; i <= DEF_lbar; i++)
        offset += DEF_Bgbar / 2 * (1UL << (64 - i * DEF_Bgbitbar));
    return offset;
}

inline void Decompositionlvl1(DecomposedTRLWElvl1 &decvec,
                              const TRLWElvl1 &trlwe)
{
    static constexpr uint32_t offset = offsetgenlvl1();
    Decomposition<uint32_t, DEF_N, DEF_l, DEF_Bgbit>(decvec, trlwe, offset);
}

inline void Decompositionlvl2(DecomposedTRLWElvl2 &decvec,
                              const TRLWElvl2 &trlwe)
{
    static constexpr uint64_t offset = offsetgenlvl2();
    Decomposition<uint64_t, DEF_nbar, DEF_lbar, DEF_Bgbitbar>(decvec, trlwe,
                                                              offset);
}

inline void DecompositionFFTlvl1(DecomposedTRLWEInFDlvl1 &decvecfft,
                                 const TRLWElvl1 &trlwe)
{
    DecomposedTRLWElvl1 decvec;
    Decompositionlvl1(decvec, trlwe);
    for (int i = 0; i < 2 * DEF_l; i++) TwistIFFTlvl1(decvecfft[i], decvec[i]);
}

inline void DecompositionFFTlvl2(DecomposedTRLWEInFDlvl2 &decvecfft,
                                 const TRLWElvl2 &trlwe)
{
    DecomposedTRLWElvl2 decvec;
    Decompositionlvl2(decvec, trlwe);
    for (int i = 0; i < 2 * DEF_lbar; i++)
        TwistIFFTlvl2(decvecfft[i], decvec[i]);
}

void trgswfftExternalProductlvl1(TRLWElvl1 &res, const TRLWElvl1 &trlwe,
                                 const TRGSWFFTlvl1 &trgswfft)
{
    alignas(32) DecomposedTRLWEInFDlvl1 decvecfft;
    DecompositionFFTlvl1(decvecfft, trlwe);
    alignas(32) TRLWEInFDlvl1 restrlwefft;
    MulInFD<DEF_N>(restrlwefft[0], decvecfft[0], trgswfft[0][0]);
    MulInFD<DEF_N>(restrlwefft[1], decvecfft[0], trgswfft[0][1]);
    for (int i = 1; i < 2 * DEF_l; i++) {
        FMAInFD<DEF_N>(restrlwefft[0], decvecfft[i], trgswfft[i][0]);
        FMAInFD<DEF_N>(restrlwefft[1], decvecfft[i], trgswfft[i][1]);
    }
    TwistFFTlvl1(res[0], restrlwefft[0]);
    TwistFFTlvl1(res[1], restrlwefft[1]);
}

void trgswfftExternalProductlvl2(TRLWElvl2 &res, const TRLWElvl2 &trlwe,
                                 const TRGSWFFTlvl2 &trgswfft)
{
    DecomposedTRLWEInFDlvl2 decvecfft;
    DecompositionFFTlvl2(decvecfft, trlwe);
    TRLWEInFDlvl2 restrlwefft;
    MulInFD<DEF_nbar>(restrlwefft[0], decvecfft[0], trgswfft[0][0]);
    MulInFD<DEF_nbar>(restrlwefft[1], decvecfft[0], trgswfft[0][1]);
    for (int i = 1; i < 2 * DEF_lbar; i++) {
        FMAInFD<DEF_nbar>(restrlwefft[0], decvecfft[i], trgswfft[i][0]);
        FMAInFD<DEF_nbar>(restrlwefft[1], decvecfft[i], trgswfft[i][1]);
    }
    TwistFFTlvl2(res[0], restrlwefft[0]);
    TwistFFTlvl2(res[1], restrlwefft[1]);
}
}  // namespace TFHEpp