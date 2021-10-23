#pragma once

#include <array>
#include <cassert>
#include <cstdint>

namespace cuHEpp {
template <typename T>
constexpr bool false_v = false;

constexpr uint64_t P = (((1ULL << 32) - 1) << 32) + 1;

// this class defines operations over integaer torus.
class INTorus {
public:
    uint64_t value;
    INTorus() { value = 0; }
    INTorus(uint64_t data, bool modulo = true)
    {
        if (modulo)
            value = data + static_cast<uint32_t>(-(data >= P));
        else
            value = data;
    }

    // return this + b mod P.
    INTorus operator+(const INTorus &b) const
    {
        uint64_t tmp = this->value + b.value;
        return INTorus(
            tmp + static_cast<uint32_t>(-(tmp < b.value || tmp >= P)), false);
    }

    INTorus &operator+=(const INTorus &b)
    {
        this->value += b.value;
        *this = INTorus(
            this->value + static_cast<uint32_t>(
                              -(this->value < b.value || this->value >= P)),
            false);
        return *this;
    }

    // return this - b mod P.
    INTorus operator-(const INTorus &b) const
    {
        uint64_t tmp = this->value - b.value;
        return INTorus(tmp - static_cast<uint32_t>(-(tmp > (this->value))),
                       false);
    }

    INTorus operator-=(const INTorus &b)
    {
        INTorus tmp = *this - b;
        *this = tmp;
        return *this;
    }

    INTorus operator*(const INTorus &b) const
    {
        __uint128_t tmp = static_cast<__uint128_t>(this->value) * b.value;
        uint32_t *tmpa = reinterpret_cast<uint32_t *>(&tmp);
        uint64_t res = ((static_cast<uint64_t>(tmpa[1]) + tmpa[2]) << 32) +
                       tmpa[0] - tmpa[3] - tmpa[2];
        uint64_t lo = static_cast<uint64_t>(tmp);
        res -= static_cast<uint32_t>(-((res > lo) && (tmpa[2] == 0)));
        res += static_cast<uint32_t>(-((res < lo) && (tmpa[2] != 0)));
        return INTorus(res);
    }

    INTorus operator*=(const INTorus &b)
    {
        const INTorus tmp = *this * b;
        *this = tmp;
        return *this;
    }

    INTorus operator<<(uint32_t l) const
    {
        if (l == 0) {
            return *this;
        }
        // t[0] = templ,t[1] = tempul, t[2] = tempuu
        else if (l < 32) {
            uint64_t templ, tempu, res;
            templ = this->value << l;
            tempu = this->value >> (64 - l);
            res = templ + (tempu << 32) - tempu;
            res +=
                static_cast<uint32_t>(-(res < templ));  // tempuu is always 0.
            return INTorus(res);
        }
        else if (l == 32) {
            uint64_t templ, tempul, tempuu, res;
            templ = this->value << l;
            tempul = static_cast<uint32_t>(this->value >> (64 - l));
            tempuu = 0;
            res = templ + (tempul << 32) - tempuu - tempul;
            res -= static_cast<uint32_t>(-((res > templ) && (tempul == 0)));
            res += static_cast<uint32_t>(-((res < templ) && (tempul != 0)));
            return INTorus(res);
        }
        else if (l < 64) {
            uint64_t templ, tempul, tempuu, res;
            templ = static_cast<uint32_t>(this->value << (l - 32));
            tempul = static_cast<uint32_t>(this->value >> (64 - l));
            tempuu = this->value >> (96 - l);
            res = ((templ + tempul) << 32) - tempuu - tempul;
            res -= static_cast<uint32_t>(
                -((res > (templ << 32)) && (tempul == 0)));
            res += static_cast<uint32_t>(
                -((res < (templ << 32)) && (tempul != 0)));
            return INTorus(res);
        }
        else if (l == 64) {
            uint64_t templ, tempu, res;
            templ = static_cast<uint32_t>(this->value);
            templ = (templ << 32) - templ;
            // templ += static_cast<uint32_t>(-(templ >= P));//mod P
            tempu = this->value >> (96 - l);
            res = templ - tempu;
            res -= static_cast<uint32_t>(-(res > (templ)));
            return INTorus(res);
        }
        else if (l < 96) {
            // different from cuFHE
            uint64_t templ, tempu, res;
            templ = static_cast<uint32_t>(this->value << (l - 64));
            templ = (templ << 32) - templ;
            // templ += static_cast<uint32_t>(-(templ >= P)); //mod P
            tempu = this->value >> (96 - l);
            res = templ - tempu;
            res -= static_cast<uint32_t>(-(res > (templ)));
            return INTorus(res);
        }
        else if (l == 96) {
            uint64_t templ, tempu, res;
            templ = P - (this->value);
            tempu = 0;
            res = tempu + templ;
            res += static_cast<uint32_t>(-(res < tempu));
            return INTorus(res);
        }
        else if (l < 128) {
            // Same as cuFHE
            uint64_t templ, tempu, res;
            templ = this->value << (l - 96);
            tempu = this->value >> (160 - l);
            res = templ + (tempu << 32) - tempu;
            res += static_cast<uint32_t>(-(res < templ));
            return INTorus(P - INTorus(res).value);
        }
        else if (l == 128) {
            uint64_t templ, tempul /*,tmempuu*/;
            INTorus res;
            templ = static_cast<uint32_t>(this->value);
            tempul = static_cast<uint32_t>(this->value >> (160 - l));
            // res = -((templ+tempul)<<32)+tempul-tempuu;
            res = INTorus(tempul, false) - INTorus(templ << 32, false) -
                  INTorus(tempul << 32, false);  //-INTorus(tempuu,false);
            return res;
        }
        else if (l < 160) {
            uint64_t /*templul,*/ templ, tempul, tempuu;
            INTorus res;
            templ = static_cast<uint32_t>(this->value << (l - 128));
            tempul = static_cast<uint32_t>(this->value >> (160 - l));
            tempuu = this->value >> (192 - l);
            // res = -((templ+tempul)<<32)+tempul-tempuu;
            res = INTorus(tempul + tempuu, false) -
                  INTorus(templ << 32, false) - INTorus((tempul << 32), false);
            return res;
        }
        else if (l == 160) {
            uint64_t templ, tempu;
            INTorus res;
            templ = static_cast<uint32_t>(this->value);
            tempu = this->value >> (192 - l);
            res = INTorus(templ + tempu, false) - INTorus(templ << 32, false);
            return res;
        }
        else {
            uint64_t templ, tempu, res;
            templ = static_cast<uint32_t>(this->value) << (l - 160);
            tempu = this->value >> (192 - l);
            res = templ + tempu - (templ << 32);
            res -= static_cast<uint32_t>(-(res > tempu));
            return INTorus(res);
        }
    }

    INTorus Pow(uint64_t e) const
    {
        INTorus res(1, false);
        for (uint64_t i = 0; i < e; i++) res *= *this;
        return res;
    }

    template <class Archive>
    void serialize(Archive &ar)
    {
        ar(value);
    }
};

// defined on [1,31]
inline INTorus InvPow2(uint8_t nbit)
{
    uint32_t low, high;
    low = (1 << (32 - nbit)) + 1;
    high = -low;
    return INTorus((static_cast<uint64_t>(high) << 32) + low);
}
}  // namespace cuHEpp
