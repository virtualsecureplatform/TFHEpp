#pragma once

#include "cloudkey.hpp"

namespace TFHEpp {
using namespace std;
void HomCONSTANTONE(TLWE<lvl1param> &res);
void HomCONSTANTZERO(TLWE<lvl1param> &res);
void HomNOT(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca);
void HomCOPY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca);
void HomNAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomXNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
           const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomXOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomANDNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomANDYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomORNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomORYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);
void HomMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
            const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
            const EvalKey &ek);
void HomNMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
             const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
             const EvalKey &ek);
template <class iksP, class bkP>
void HomMUXwoSE(TRLWE<typename bkP::targetP> &res,
                const TLWE<typename iksP::domainP> &cs,
                const TLWE<typename iksP::domainP> &c1,
                const TLWE<typename iksP::domainP> &c0,
                const EvalKey &ek);
void ExtractSwitchAndHomMUX(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &csr,
                            const TRLWE<lvl1param> &c1r,
                            const TRLWE<lvl1param> &c0r, const EvalKey &ek);
}  // namespace TFHEpp