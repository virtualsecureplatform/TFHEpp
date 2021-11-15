#pragma once

#include "axell/gate.hpp"
#include "cloudkey.hpp"
#include "params/128bit.hpp"

namespace TFHEpp {
using namespace std;
template <class P = lvl1param>
void HomCONSTANTONE(TLWE<lvl1param> &res);
template <class P = lvl1param>
void HomCONSTANTZERO(TLWE<lvl1param> &res);
void HomNOT(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca);
void HomCOPY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca);
template <class P = lvl1param>
void HomNAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomXNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
           const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomXOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomANDNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomANDYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomORNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomORYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);
template <class P = lvl1param>
void HomMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
            const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
            const EvalKey &ek);
template <class P = lvl1param>
void HomNMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
             const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
             const EvalKey &ek);
template <class iksP, class bkP>
void HomMUXwoSE(TRLWE<typename bkP::targetP> &res,
                const TLWE<typename iksP::domainP> &cs,
                const TLWE<typename iksP::domainP> &c1,
                const TLWE<typename iksP::domainP> &c0, const EvalKey &ek);
void ExtractSwitchAndHomMUX(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &csr,
                            const TRLWE<lvl1param> &c1r,
                            const TRLWE<lvl1param> &c0r, const EvalKey &ek);
}  // namespace TFHEpp
