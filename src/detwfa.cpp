#include <trgsw.hpp>

namespace TFHEpp {
template <class P>
inline void CMUXFFT(TRLWE<P> &res, const TRGSWFFT<P> &cs, const TRLWE<P> &c1,
                    const TRLWE<P> &c0)
{
    for (int i = 0; i < P::n; i++) {
        res[0][i] = c1[0][i] - c0[0][i];
        res[1][i] = c1[1][i] - c0[1][i];
    }
    trgswfftExternalProduct<P>(res, res, cs);
    for (int i = 0; i < P::n; i++) {
        res[0][i] += c0[0][i];
        res[1][i] += c0[1][i];
    }
}
void CMUXFFTlvl1(TRLWE<lvl1param> &res, const TRGSWFFT<lvl1param> &cs,
                 const TRLWE<lvl1param> &c1, const TRLWE<lvl1param> &c0)
{
    CMUXFFT<lvl1param>(res, cs, c1, c0);
}
}  // namespace TFHEpp