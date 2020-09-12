#pragma once

#include <cstdint>
#include <array>

namespace cuHEpp{
    using namespace std;

    constexpr uint64_t P = (((1UL<<32)-1)<<32)+ 1;

    class INTorus{
        public:
            uint64_t value;
            INTorus(){value=0;}
            INTorus(uint64_t data,bool modulo = true){
                if(modulo) value = data + static_cast<uint32_t>(-(data >= P));
                else value = data;
            }

            INTorus operator + (INTorus b){
                uint64_t tmp = this->value + b.value;
                return INTorus(tmp + static_cast<uint32_t>(-(tmp < b.value || tmp >= P)),false);
            }
            
            INTorus operator - (INTorus b){
                uint64_t tmp = this->value - b.value;
                return INTorus(tmp - static_cast<uint32_t>(-(tmp>(this->value))),false);
            }

            INTorus operator * (INTorus b){
                __uint128_t tmp = static_cast<__uint128_t>(this->value) * b.value;
                uint32_t* tmpa = reinterpret_cast<uint32_t*>(&tmp);
                uint64_t res = ((static_cast<uint64_t>(tmpa[1])+tmpa[2])<<32) + tmpa[0] - tmpa[3] - tmpa[2];
                uint64_t lo = static_cast<uint64_t>(tmp);
                res -= static_cast<uint32_t>(-((res>lo)&&(tmpa[2]==0)));
                res += static_cast<uint32_t>(-((res<lo)&&(tmpa[2]!=0)));
                return INTorus(res);
            }

            INTorus operator << (uint32_t l){
                if(l==0){
                    return *this;
                }
                //t[0] = templ,t[1] = tempul, t[2] = tempuu
                else if(l<32){
                    uint64_t templ,tempu,res;
                    templ = this->value<<l;
                    tempu = this->value>>(64-l);
                    res = templ+(tempu<<32)-tempu;
                    res += static_cast<uint32_t>(-(res<templ));// tempuu is always 0.
                    return INTorus(res);
                }
                else if(l==32){
                    uint64_t templ,tempul,tempuu,res;
                    templ = this->value<<l;
                    tempul = static_cast<uint32_t>(this->value>>(64-l));
                    tempuu = 0;
                    res = templ+(tempul<<32)-tempuu-tempul;
                    res -= static_cast<uint32_t>(-((res>templ)&&(tempul==0)));
                    res += static_cast<uint32_t>(-((res<templ)&&(tempul!=0)));
                    return INTorus(res);
                }
                else if(l<64){
                    uint64_t templ,tempul,tempuu,res;
                    templ = static_cast<uint32_t>(this->value<<(l-32));
                    tempul = static_cast<uint32_t>(this->value>>(64-l));
                    tempuu = this->value>>(96-l);
                    res = ((templ+tempul)<<32)-tempuu-tempul;
                    res -= static_cast<uint32_t>(-((res>(templ<<32))&&(tempul==0)));
                    res += static_cast<uint32_t>(-((res<(templ<<32))&&(tempul!=0)));
                    return INTorus(res);
                }
                // Above are almsot equivalent
                //cuFHE seems to be wrong in Lsh96 and Lsh128. Sign is wrong.
                //mod P is not needed
                else if(l==64){
                    uint64_t templ,tempu,res;
                    templ = static_cast<uint32_t>(this->value);
                    templ = (templ<<32) - templ;
                    // templ += static_cast<uint32_t>(-(templ >= P));//mod P
                    tempu = this->value>>(96-l);
                    res = templ-tempu;
                    res -= static_cast<uint32_t>(-(res > (templ)));
                    return INTorus(res);
                }
                else if(l<96){
                    uint64_t templ,tempu,res;
                    templ = static_cast<uint32_t>(this->value<<(l-64));
                    templ = (templ<<32) - templ;
                    // templ += static_cast<uint32_t>(-(templ >= P)); //mod P
                    tempu = this->value>>(96-l);
                    res = templ-tempu;
                    res -= static_cast<uint32_t>(-(res > (templ)));
                    return INTorus(res);
                }
                else if(l==96){
                    uint64_t templ,tempu,res;
                    templ = P-(this->value);
                    tempu = 0;
                    res = tempu + templ;
                    res += static_cast<uint32_t>(-(res < tempu));
                    return INTorus(res);
                }
                else if(l<128){
                    uint64_t templ,tempu,res;
                    templ = P-(this->value<<(l-96));
                    tempu = this->value>>(160-l);
                    tempu = P-(tempu<<32)+tempu;
                    // tempu += static_cast<uint32_t>(-(tempu >= P)); //mod P
                    res = templ+tempu;
                    res += static_cast<uint32_t>(-(res < tempu));
                    return INTorus(res);
                }
                else if(l==128){
                    uint64_t templ,tempul,tempuu,res;
                    templ = static_cast<uint32_t>(this->value);
                    tempul = static_cast<uint32_t>(this->value>>(160-l));
                    tempuu = 0;
                    res = ((templ+tempul)<<32)-tempul-tempuu;
                    res -= static_cast<uint32_t>(-(res > (templ << 32) && tempul == 0));
                    res += static_cast<uint32_t>(-(res < (templ << 32) && tempul != 0));
                    return INTorus(res);
                }
                else if(l<160){

                    uint64_t templ,tempul,tempuu,res;
                    templ = static_cast<uint32_t>(this->value<<(l-128));
                    tempul = static_cast<uint32_t>(this->value>>(160-l));
                    tempuu = this->value>>(192-l);
                    res = ((templ+tempul)<<32)-tempul-tempuu;
                    res -= static_cast<uint32_t>(-(res > (templ << 32) && tempul == 0));
                    res += static_cast<uint32_t>(-(res < (templ << 32) && tempul != 0));
                    return INTorus(res);
                }
                else{
                    uint64_t templ,tempu,res;
                    templ = this->value>>(192-l);
                    tempu = this->value<<(l-160);
                    res = templ + tempu - (tempu<<32);
                    res += static_cast<uint32_t>(-(res > templ));
                    return INTorus(res);
                }
            }
    };
}