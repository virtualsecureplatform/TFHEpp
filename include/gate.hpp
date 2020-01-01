#pragma once

#include <cloudkey.hpp>

namespace TFHEpp {
using namespace std;
void HomNAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomNOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomXNOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomAND(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomXOR(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomNOT(TLWElvl0 &res, const TLWElvl0 &ca);
void HomCOPY(TLWElvl0 &res, const TLWElvl0 &ca);
void HomCONSTANTONE(TLWElvl0 &res);
void HomCONSTANTZERO(TLWElvl0 &res);
void HomANDNY(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomANDYN(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomORNY(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomORYN(TLWElvl0 &res, const TLWElvl0 &ca, const TLWElvl0 &cb,
             const CloudKey &ck);
void HomMUX(TLWElvl0 &res, const TLWElvl0 &cs, const TLWElvl0 &c1,
            const TLWElvl0 &c0, const CloudKey &ck);
}  // namespace TFHEpp