#pragma once

#include "axell/gate.hpp"
#include "cloudkey.hpp"

namespace TFHEpp {
using namespace std;
template <class P = lvl1param>
void HomCONSTANTONE(TLWE<P> &res);
template <class P = lvl1param>
void HomCONSTANTZERO(TLWE<P> &res);
template <class P = lvl1param>
void HomNOT(TLWE<P> &res, const TLWE<P> &ca);
template <class P = lvl1param>
void HomCOPY(TLWE<P> &res, const TLWE<P> &ca);
template <class P = lvl1param>
void HomNAND(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek);
template <class P = lvl1param>
void HomNOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek);
template <class P = lvl1param>
void HomXNOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek);
template <class P = lvl1param>
void HomAND(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek);
template <class P = lvl1param>
void HomOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
           const EvalKey &ek);
template <class P = lvl1param>
void HomXOR(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const EvalKey &ek);
template <class P = lvl1param>
void HomANDNY(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
              const EvalKey &ek);
template <class P = lvl1param>
void HomANDYN(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
              const EvalKey &ek);
template <class P = lvl1param>
void HomORNY(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek);
template <class P = lvl1param>
void HomORYN(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const EvalKey &ek);
template <class P = lvl1param>
void HomMUX(TLWE<P> &res, const TLWE<P> &cs, const TLWE<P> &c1,
            const TLWE<P> &c0, const EvalKey &ek);
template <class P = lvl1param>
void HomNMUX(TLWE<P> &res, const TLWE<P> &cs, const TLWE<P> &c1,
             const TLWE<P> &c0, const EvalKey &ek);
template <class iksP, class bkP>
void HomMUXwoIKSandSE(TRLWE<typename bkP::targetP> &res,
                      const TLWE<typename bkP::domainP> &cs,
                      const TLWE<typename bkP::domainP> &c1,
                      const TLWE<typename bkP::domainP> &c0, const EvalKey &ek);
template <class iksP, class bkP>
void HomMUXwoSE(TRLWE<typename bkP::targetP> &res,
                const TLWE<typename iksP::domainP> &cs,
                const TLWE<typename iksP::domainP> &c1,
                const TLWE<typename iksP::domainP> &c0, const EvalKey &ek);
void ExtractSwitchAndHomMUX(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &csr,
                            const TRLWE<lvl1param> &c1r,
                            const TRLWE<lvl1param> &c0r, const EvalKey &ek);
}  // namespace TFHEpp
