#include<gatebootstrapping.hpp>

namespace TFHEpp{
    void HomNAND(array<uint32_t,DEF_n+1> &res, const array<uint32_t,DEF_n+1> &ca, const array<uint32_t,DEF_n+1> &cb,CloudKey ck){
        for(int i = 0;i<=DEF_n;i++) res[i] = -ca[i] - cb[i];
        res[DEF_n] += 1U<<29;
        GateBootstrapping(res,res,ck);
    }
}