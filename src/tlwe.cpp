#include<array>
#include<cstdint>
#include<random>
#include <type_traits>
#include<vector>

#include<params.hpp>
#include<limits>
#include<utils.hpp>
#include<key.hpp>

namespace TFHEpp{
    using namespace std;
    static random_device engine;

    template<typename T = uint32_t,uint32_t n = DEF_n>
    array<T,n+1> tlweSymEncrypt(const T p,const double α, const array<uint32_t,n> key){
        uniform_int_distribution<T> Torusdist(0, numeric_limits<T>::max());
        array<uint32_t,n+1> res = {};
        res[n] = gaussian32(p,α);
        for(int i = 0;i < n;i++){
            res[i] = Torusdist(engine);
            res[n] += res[i]*key[i];
        }
        return res;
    }

    array<uint32_t,DEF_n+1> tlweSymEncryptlvl1(const uint32_t p, const double α, const array<uint32_t,DEF_n> key){
        return tlweSymEncrypt<uint32_t,DEF_n>(p,α,key);
    }

    template<typename T = uint32_t,uint32_t n = DEF_n>
    bool tlweSymDecrypt(const array<T,n+1> c, const array<T,n> key){
        T phase = c[n];
        for(int i = 0;i<n;i++) phase -= c[i]*key[i];
        bool res = static_cast<typename make_signed<T>::type>(phase)>0;
        return res;
    }

    bool tlweSymDecryptlvl1(const array<uint32_t,DEF_n+1> c, const array<uint32_t,DEF_n> key){
        return tlweSymDecrypt<uint32_t,DEF_n>(c,key);
    }

    vector<array<uint32_t,DEF_n+1>> bootsSymEncrypt(const vector<bool> p,const SecretKey sk){
        vector<array<uint32_t,DEF_n+1>> c(p.size());
        for(int i = 0;i<p.size();i++) c[i] = tlweSymEncryptlvl1(p[i]?DEF_MU:-DEF_MU, DEF_α,sk.key.lvl0);
        return c;
    }

    vector<bool> bootsSymDecrypt(const vector<array<uint32_t,DEF_n+1>> c, const SecretKey sk){
        vector<bool> p(c.size());
        for(int i = 0;i<p.size();i++) p[i] = tlweSymDecryptlvl1(c[i], sk.key.lvl0);
        return p;
    }
}