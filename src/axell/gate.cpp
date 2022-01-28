#include <axell/gate.hpp>

namespace TFHEpp {

template <class P>
constexpr Polynomial<P> swaptestvecgen()
{
    Polynomial<P> testvec = {};
    for (int i = 0; i < P::n / 2; i++) {
        testvec[i] = 3 * P::μ;
        testvec[i + P::n / 2] = P::μ;
    }
    return testvec;
}

void HomSWAP(TLWE<lvl1param> &resa, TLWE<lvl1param> &resb,
             const TLWE<lvl1param> &cs, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    TLWE<lvl1param> intemp;
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++)
        intemp[i] = 2 * ca[i] + cb[i] + 2 * cs[i];
    intemp[lvl1param::k * lvl1param::n] += 2 * lvl1param::μ;
    TLWE<lvl0param> intemplvl0;
    IdentityKeySwitch<lvl10param>(intemplvl0, intemp, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        intemp, intemplvl0, *ek.bkfftlvl01, swaptestvecgen<lvl1param>());
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++) intemp[i] += cs[i];
    intemp[lvl1param::k * lvl1param::n] += lvl1param::μ;
    IdentityKeySwitch<lvl10param>(intemplvl0, intemp, *ek.iksklvl10);
    TRLWE<lvl1param> trlwe;
    BlindRotate<lvl01param>(trlwe, intemplvl0, *ek.bkfftlvl01,
                            μpolygen<lvl1param, lvl1param::μ>());
    SampleExtractIndex<lvl1param>(resa, trlwe, 0);
    SampleExtractIndex<lvl1param>(resb, trlwe, lvl1param::n / 2);
}

// process AND, NOR, XOR of same 2 input at once
template <class P>
void HomHalfAdder(TLWE<lvl1param> &carry, TLWE<lvl1param> &sum,
                  const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
                  const EvalKey &ek)
{
    TRLWE<P> trlwe;
    for (int i = 0; i <= P::k * P::n; i++)
        sum[i] = c0[i] + c1[i];  // use sum as buffer
    TLWE<lvl0param> cadd;
    IdentityKeySwitch<lvl10param>(cadd, sum, *ek.iksklvl10);
    BlindRotate<lvl01param>(trlwe, cadd, *ek.bkfftlvl01, μpolygen<P, P::μ>());
    TLWE<P> cor;
    constexpr typename lvl01param::domainP::T roundoffset =
        1ULL << (std::numeric_limits<typename lvl01param::domainP::T>::digits -
                 2 - lvl01param::targetP::nbit);
    SampleExtractIndex<P>(
        carry, trlwe,
        ((-P::μ + roundoffset) >>
         (std::numeric_limits<typename lvl01param::domainP::T>::digits - 1 -
          lvl01param::targetP::nbit)) &
            (P::n - 1));
    SampleExtractIndex<P>(
        cor, trlwe,
        ((P::μ + roundoffset) >>
         (std::numeric_limits<typename lvl01param::domainP::T>::digits - 1 -
          lvl01param::targetP::nbit)) &
            (P::n - 1));
    for (int i = 0; i <= P::k * P::n; i++) {
        carry[i] = -carry[i];
        sum[i] = cor[i] - carry[i];
    }
    sum[P::k * P::n] -= P::μ;
}
#define INST(P)                                               \
    template void HomHalfAdder<P>(                            \
        TLWE<lvl1param> & carry, TLWE<lvl1param> & sum,       \
        const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0, \
        const EvalKey &ek)
INST(lvl1param);
INST(lvlMparam);
#undef INST

void Hom2BRFullAdder(TLWE<lvl1param> &carry, TLWE<lvl1param> &sum,
                     const TLWE<lvl1param> &ca, const TLWE<lvl1param> &cb,
                     const TLWE<lvl1param> &cc, const EvalKey &ek)
{
    TLWE<lvl1param> intemp;
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++) intemp[i] = ca[i] + cb[i] + cc[i];
    TLWE<lvl0param> intemplvl0;
    IdentityKeySwitch<lvl10param>(intemplvl0, intemp, *ek.iksklvl10);
    // These are parallel
    // carry
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        carry, intemplvl0, *ek.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
    // sum
    for (int i = 0; i <= lvl0param::n; i++) intemplvl0[i] *= -2;
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        sum, intemplvl0, *ek.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
}

template <class P, uint32_t num_input>
constexpr Polynomial<P> multiinputtestvecgen()
{
    Polynomial<P> testvec = {};
    for (int i = 0; i < P::n / num_input; i++) {
        testvec[i] = P::μ;
    }
    return testvec;
}

template <class P>
void fulladdertestvectorspecialize(TRLWE<P> &rcarry, TRLWE<P> &rsum,
                                   TRLWE<P> &trlwe)
{
    for(int k = 0; k < P::k+1; k++){
        for (int i = 0; i < P::n / 3; i++) {
            rcarry[k][i] = trlwe[k][i] - trlwe[k][i + 1 + 2 * P::n / 3] -
                        trlwe[k][i + 1 + 1 * P::n / 3];
            rsum[k][i] = -trlwe[k][i] - trlwe[k][i + 1 + 2 * P::n / 3] +
                        trlwe[k][i + 1 + 1 * P::n / 3];
        }
        for (int i = 0; i < P::n / 3; i++) {
            rcarry[k][i + 1 * P::n / 3] = trlwe[k][i + 1 * P::n / 3] + trlwe[k][i] -
                                        trlwe[k][i + 1 + 2 * P::n / 3];
            rsum[k][i + 1 * P::n / 3] = -trlwe[k][i + 1 * P::n / 3] + trlwe[k][i] +
                                        trlwe[k][i + 1 + 2 * P::n / 3];
        }
        for (int i = 0; i < P::n / 3 + 1; i++) {
            rcarry[k][i + 2 * P::n / 3] = trlwe[k][i + 2 * P::n / 3] +
                                        trlwe[k][i + 1 * P::n / 3] + trlwe[k][i];
            rsum[k][i + 2 * P::n / 3] = -trlwe[k][i + 2 * P::n / 3] +
                                        trlwe[k][i + 1 * P::n / 3] - trlwe[k][i];
        }
    }
}

// process 3input XOR and majority gate at once
void HomFullAdder(TLWE<lvlMparam> &carry, TLWE<lvlMparam> &sum,
                  const TLWE<lvlMparam> &ca, const TLWE<lvlMparam> &cb,
                  const TLWE<lvlMparam> &cc, const EvalKey &ek)
{
    for (int i = 0; i <= lvlMparam::k * lvlMparam::n; i++) sum[i] = ca[i] + cb[i] + cc[i];
    TRLWE<lvlMparam> rtemp;
    TLWE<lvl0param> temp;
    IdentityKeySwitch<lvlM0param>(temp, sum, *ek.iksklvl10);
    BlindRotate<lvl0Mparam>(rtemp, temp, *ek.bkfftlvl01,
                            multiinputtestvecgen<lvlMparam, 3>());
    TRLWE<lvlMparam> rcarry, rsum;
    fulladdertestvectorspecialize<lvlMparam>(rcarry, rsum, rtemp);
    SampleExtractIndex<lvlMparam>(carry, rcarry, 0);
    SampleExtractIndex<lvlMparam>(sum, rsum, 0);
}

// process AND, NOR, XOR of same 2 input at once
void HomXORNANDNOR(TLWE<lvl1param> &cxor, TLWE<lvl1param> &cnand,
                   TLWE<lvl1param> &cnor, const TLWE<lvl1param> &ca,
                   const TLWE<lvl1param> &cb, const EvalKey &ek)
{
    TLWE<lvl1param> cadd;
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++) cadd[i] = ca[i] + cb[i];
    TLWE<lvl0param> temp;
    IdentityKeySwitch<lvl10param>(temp, cadd, *ek.iksklvl10);
    TRLWE<lvl1param> trlwe;
    BlindRotate<lvl01param>(trlwe, temp, *ek.bkfftlvl01,
                            μpolygen<lvl1param, lvl1param::μ>());
    constexpr typename lvl01param::domainP::T roundoffset =
        1ULL << (std::numeric_limits<typename lvl01param::domainP::T>::digits -
                 2 - lvl01param::targetP::nbit);
    SampleExtractIndex<lvl1param>(
        cnand, trlwe,
        ((-lvl0param::μ + roundoffset) >>
         (std::numeric_limits<typename lvl01param::domainP::T>::digits - 1 -
          lvl01param::targetP::nbit)) &
            (lvl1param::n - 1));
    SampleExtractIndex<lvl1param>(
        cnor, trlwe,
        ((lvl0param::μ + roundoffset) >>
         (std::numeric_limits<typename lvl01param::domainP::T>::digits - 1 -
          lvl01param::targetP::nbit)) &
            (lvl1param::n - 1));
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++) {
        cnor[i] = -cnor[i];
        cxor[i] = -cnor[i] + cnand[i];
    }
    cxor[lvl1param::n] -= lvl1param::μ;
}

void Hom4inputOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                 const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
                 const TLWE<lvl1param> &cd, const EvalKey &ek)
{
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++)
        res[i] = ca[i] + cb[i] + cc[i] + cd[i];
    res[lvl1param::k * lvl1param::n] += (1U << 30) - (1U << 28);
    TLWE<lvl0param> temp;
    IdentityKeySwitch<lvl10param>(temp, res, *ek.iksklvl10);
    TRLWE<lvl1param> rtemp;
    BlindRotate<lvl01param>(rtemp, temp, *ek.bkfftlvl01,
                            μpolygen<lvl1param, lvl1param::μ>());
    SampleExtractIndex<lvl1param>(res, rtemp, 0);
}

void Hom4inputAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                  const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
                  const TLWE<lvl1param> &cd, const EvalKey &ek)
{
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++)
        res[i] = ca[i] + cb[i] + cc[i] + cd[i];
    res[lvl1param::k * lvl1param::n] -= (1U << 30) - (1U << 28);
    TLWE<lvl0param> temp;
    IdentityKeySwitch<lvl10param>(temp, res, *ek.iksklvl10);
    TRLWE<lvl1param> rtemp;
    BlindRotate<lvl01param>(rtemp, temp, *ek.bkfftlvl01,
                            μpolygen<lvl1param, lvl1param::μ>());
    SampleExtractIndex<lvl1param>(res, rtemp, 0);
}

template <class P>
constexpr Polynomial<P> xortestvecgen()
{
    Polynomial<P> xortestvec;
    for (int i = 0; i < P::n / 4; i++) {
        xortestvec[i] = -P::μ;
        xortestvec[i + P::n / 4] = P::μ;
        xortestvec[i + P::n / 2] = -P::μ;
        xortestvec[i + 3 * P::n / 4] = P::μ;
    }
    return xortestvec;
}

void Hom3inputXOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                  const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
                  const EvalKey &ek)
{
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++) res[i] = ca[i] + cb[i] + cc[i];
    res[lvl1param::k * lvl1param::n] += 2 * lvl1param::μ;
    TLWE<lvl0param> temp;
    IdentityKeySwitch<lvl10param>(temp, res, *ek.iksklvl10);
    TRLWE<lvl1param> rtemp;
    BlindRotate<lvl01param>(rtemp, temp, *ek.bkfftlvl01,
                            xortestvecgen<lvl1param>());
    SampleExtractIndex<lvl1param>(res, rtemp, 0);
}

void Hom3inputThreashold(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                         const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
                         const EvalKey &ek)
{
    res[lvl1param::k * lvl1param::n] += lvl1param::μ / 2;
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++) res[i] = ca[i] + cb[i] + cc[i];
    TLWE<lvl0param> temp;
    IdentityKeySwitch<lvl10param>(temp, res, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(
        res, temp, *ek.bkfftlvl01, μpolygen<lvl1param, lvl1param::μ>());
}

template <class P>
constexpr Polynomial<P> aoi3testvecgen()
{
    Polynomial<P> testvec;
    for (int i = 0; i < P::n; i++) testvec[i] = -P::μ;
    return testvec;
}

template <class P>
void HomAOI3(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const TLWE<P> &cc, const EvalKey &ek)
{
    for (int i = 0; i <= P::k * P::n; i++) res[i] = ca[i] + cb[i] + 2 * cc[i];
    res[P::k * P::n] += (1ULL << std::numeric_limits<typename P::T>::digits) / 12;
    TLWE<lvl0param> temp;
    IdentityKeySwitch<typename KeySwitchParam<P>::P>(temp, res, *ek.iksklvl10);
    TRLWE<P> rtemp;
    BlindRotate<typename BlindRotateParam<P>::P>(rtemp, temp, *ek.bkfftlvl01,
                                                 aoi3testvecgen<P>());
    SampleExtractIndex<P>(res, rtemp, 0);
}
#define INST(P)                                                    \
    template void HomAOI3<P>(TLWE<P> & res, const TLWE<P> &ca,     \
                             const TLWE<P> &cb, const TLWE<P> &cc, \
                             const EvalKey &ek);
INST(lvl1param)
INST(lvlMparam)
#undef INST

template <class P>
void HomAO3(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const TLWE<P> &cc, const EvalKey &ek)
{
    HomAOI3<P>(res, ca, cb, cc, ek);
    for (int i = 0; i <= P::k * P::n; i++) res[i] = -res[i];
}
#define INST(P)                                                   \
    template void HomAO3<P>(TLWE<P> & res, const TLWE<P> &ca,     \
                            const TLWE<P> &cb, const TLWE<P> &cc, \
                            const EvalKey &ek);
INST(lvl1param)
INST(lvlMparam)
#undef INST

template <class P>
void HomOA3(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
            const EvalKey &ek)
{
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++) res[i] = ca[i] + cb[i] + 2 * cc[i];
    res[lvl1param::k * lvl1param::n] +=
        (1ULL << std::numeric_limits<typename lvl1param::T>::digits) / 12;
    TLWE<lvl0param> temp;
    IdentityKeySwitch<lvl10param>(temp, res, *ek.iksklvl10);
    TRLWE<lvl1param> rtemp;
    BlindRotate<lvl01param>(rtemp, temp, *ek.bkfftlvl01, aoi3testvecgen<P>());
    SampleExtractIndex<lvl1param>(res, rtemp, lvl1param::n - lvl1param::n / 3);
}
#define INST(P)                                                   \
    template void HomOA3<P>(TLWE<P> & res, const TLWE<P> &ca,     \
                            const TLWE<P> &cb, const TLWE<P> &cc, \
                            const EvalKey &ek);
INST(lvl1param)
INST(lvlMparam)
#undef INST

template <class P>
void HomOAI3(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
             const EvalKey &ek)
{
    HomOA3<P>(res, ca, cb, cc, ek);
    for (int i = 0; i <= lvl1param::k * lvl1param::n; i++) res[i] = -res[i];
}
#define INST(P)                                                    \
    template void HomOAI3<P>(TLWE<P> & res, const TLWE<P> &ca,     \
                             const TLWE<P> &cb, const TLWE<P> &cc, \
                             const EvalKey &ek);
INST(lvl1param)
INST(lvlMparam)
#undef INST

}  // namespace TFHEpp
