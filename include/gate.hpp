#pragma once

#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;
void HomNAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const GateKey &gk);
void HomNOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const GateKey &gk);
void HomXNOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const GateKey &gk);
void HomAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const GateKey &gk);
void HomOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
           const GateKey &gk);
void HomXOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
            const GateKey &gk);
void HomNOT(TLWElvl0 &res, const TLWElvl0 &ca);
void HomCOPY(TLWElvl0 &res, const TLWElvl0 &ca);
void HomCONSTANTONE(TLWElvl0 &res);
void HomCONSTANTZERO(TLWElvl0 &res);
void HomANDNY(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
              const GateKey &gk);
void HomANDYN(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
              const GateKey &gk);
void HomORNY(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const GateKey &gk);
void HomORYN(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const GateKey &gk);
void HomMUX(TLWElvl0 &res, const TLWElvl0 &cs, const TLWElvl0 &c1,
            const TLWElvl0 &c0, const GateKey &gk);
void ExtractSwitchAndHomMUX(TRLWElvl1 &res, const TRLWElvl1 &csr,
                            const TRLWElvl1 &c1r, const TRLWElvl1 &c0r,
                            const GateKey &gk);
}  // namespace TFHEpp