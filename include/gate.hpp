#pragma once

#include<cloudkey.hpp>

namespace TFHEpp{
    using namespace std;
    void HomNAND(array<uint32_t,DEF_n+1> &res, const array<uint32_t,DEF_n+1> &ca, const array<uint32_t,DEF_n+1> &cb,CloudKey ck);
}