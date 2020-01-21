#include <trgsw.hpp>

namespace TFHEpp {
void CMUXFFTlvl1(TRLWElvl1 &res, const TRGSWFFTlvl1 &cs, const TRLWElvl1 &c1,
                 const TRLWElvl1 &c0)
{
    for (int i = 0; i < DEF_N; i++) {
        res[0][i] = c1[0][i] - c0[0][i];
        res[1][i] = c1[1][i] - c0[1][i];
    }
    trgswfftExternalProductlvl1(res, res, cs);
    for (int i = 0; i < DEF_N; i++) {
        res[0][i] += c0[0][i];
        res[1][i] += c0[1][i];
    }
}
}  // namespace TFHEpp