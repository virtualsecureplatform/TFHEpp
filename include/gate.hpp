#pragma once

#include "cloudkey.hpp"

namespace TFHEpp {
using namespace std;
void HomCONSTANTONE(TLWE<lvl1param> &res);
void HomCONSTANTZERO(TLWE<lvl1param> &res);
void HomNOT(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca);
void HomCOPY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca);
void HomNAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const GateKey &gk);
void HomNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const GateKey &gk);
void HomXNOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const GateKey &gk);
void HomAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const GateKey &gk);
void HomOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
           const TLWE<lvl1param> &cb, const GateKey &gk);
void HomXOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const GateKey &gk);
void HomANDNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const GateKey &gk);
void HomANDYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
              const TLWE<lvl1param> &cb, const GateKey &gk);
void HomORNY(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const GateKey &gk);
void HomORYN(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const GateKey &gk);
void HomMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
            const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
            const GateKey &gk);
void HomNMUX(TLWE<lvl1param> &res, const TLWE<lvl1param> &cs,
             const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
             const GateKey &gk);
template <class iksP, class bkP>
void HomMUXwoSE(TRLWE<typename bkP::targetP> &res,
                const TLWE<typename iksP::domainP> &cs,
                const TLWE<typename iksP::domainP> &c1,
                const TLWE<typename iksP::domainP> &c0,
                const KeySwitchingKey<iksP> &iksk,
                const BootstrappingKeyFFT<bkP> &bkfft);
void ExtractSwitchAndHomMUX(TRLWE<lvl1param> &res, const TRLWE<lvl1param> &csr,
                            const TRLWE<lvl1param> &c1r,
                            const TRLWE<lvl1param> &c0r, const GateKey &gk);
}  // namespace TFHEpp