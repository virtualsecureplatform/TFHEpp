#pragma once

#include <cstdint>
#include <array>
#include <cassert>

namespace cuHEpp{
    using namespace std;

    constexpr uint64_t P = (((1UL<<32)-1)<<32)+ 1;

    // this class defines operations over integaer torus.
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
                    uint64_t templ,tempul,tempuu;
                    INTorus res;
                    templ = static_cast<uint32_t>(this->value);
                    tempul = static_cast<uint32_t>(this->value>>(160-l));
                    tempuu = 0;
                    // res = -((templ+tempul)<<32)+tempul-tempuu;
                    res = INTorus(tempul,false)-INTorus(templ<<32,false)-INTorus(tempul<<32,false);//-INTorus(tempuu,false);
                    return res;
                }
                else if(l<160){
                    uint64_t templul,templ,tempul,tempuu;
                    INTorus res;
                    templ = static_cast<uint32_t>(this->value<<(l-128));
                    tempul = static_cast<uint32_t>(this->value>>(160-l));
                    tempuu = this->value>>(192-l);
                    // res = -((templ+tempul)<<32)+tempul-tempuu;
                    res = INTorus(tempul+tempuu,false)-INTorus(templ<<32,false)-INTorus((tempul<<32),false);
                    return res;
                }
                else if(l==160){
                    uint64_t templ,tempu;
                    INTorus res;
                    templ = static_cast<uint32_t>(this->value);
                    tempu = this->value>>(192-l);
                    res = INTorus(templ+tempu,false)-INTorus(templ<<32,false);
                    return res;
                }
                else{
                    uint64_t templ,tempu;
                    INTorus res;
                    templ = static_cast<uint32_t>(this->value)<<(l-160);
                    tempu = this->value>>(192-l);
                    res = INTorus(templ+tempu,false)-INTorus(templ<<32,false);
                    return res;
                }
            }
    };

    // NTT implementation
    // https://nufhe.readthedocs.io/en/latest/implementation_details.html
    constexpr uint64_t W = 12037493425763644479ULL;

    template<uint32_t Nbit>
    inline array<array<INTorus,1U<<Nbit>,2> TwistGen(){
        array<array<INTorus,1U<<Nbit>,2> twist;
        constexpr INTorus w = INTorus(W)<<(32-Nbit-1);
        twist[0][0] = twist[1][0] = INTorus(1,false);
        for(uint32_t i = 1;i < N;i++) twist[1][i] = twist[1][i-1]*w;
        twist[0][N-1] = twist[1][N-1]*w
        for(uint32_t i = 2;i < N;i++) twist[0][N-i] = twist[0][N-i+1]*w;
        assert((twist[0][1]*w).value == 1);
        return twist;
    }

    template<uint32_t N>
    inline array<array<INTorus,N>,2> TableGen(){
        array<array<INTorus,N>,2>table;
        constexpr INTorus w = INTorus(W)<<(32-Nbit);
        table[0] = table[1] = INTorus(1,false);
        for(uint32_t i = 1;i<N;i++) table[0][N-i] = table[1][i] = table[1][i-1]*w;
        assert((table[1][N-1]*w).value == 1);
        return table;
    }

    template<typename T = uint32_t,uint32_t N>
    inline void TwistMulInvert(array<INTorus,N> &res, const array<T,N> &a, const array<INTorus,N> &twist){
       for (int i = 0; i < N; i++) res[i] = INTorus(a[i],T == uint64_t) * twist[i];
    }

    template<uint32_t N>
    inline void ButterflyAdd(const typename array<INTorus,N>::iterator a, const typename array<INTorus,N>::iterator b, int i){
        const INTorus temp = *(a+i);
        *(a+i) += *(b+i);
        *(b+i) -= temp;
    }

    template <uint32_t Nbit, uint32_t step, uint32_t size, uint32_t stride, bool isinvert =true>
    inline void TwiddleMul(const typename array<INTorus,1<<Nbit>::iterator a, const array<INTorus,1<<Nbit> &table){
        constexpr uint32_t N = 1<<Nbit;
        if constexpr(isinvert) for(int i = 1; i<size; i++)*(a+i)=(*(a+i))*table[stride*(1<<step)*i];
        else for(int i = 1; i<size; i++)*(a+i)=(*(a+i))*table[N-stride*(1<<step)*i];
    }


}