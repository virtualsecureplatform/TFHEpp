#pragma once

#include <cstdint>
#include <limits>
#include "params/128bit.hpp"

namespace TFHEpp{
    struct portablelvl0param {
    std::uint32_t n; //dimension
    double α; //fresh noise
    std::uint32_t approx_bit; //Torus representation bit size

    portablelvl0param(): n(lvl0param::n),α(lvl0param::α),approx_bit(std::numeric_limits<typename lvl0param::T>::digits){}
    portablelvl0param(const portablelvl0param& P): n(P.n),α(P.α),approx_bit(P.approx_bit){}

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(n, α, approx_bit);
    }

    auto operator<=>(const portablelvl0param&) const = default;
};

struct portablelvl1param {
    std::uint32_t nbit; //dimension must be a power of 2 for ease of polynomial multiplication.
    std::uint32_t n; //dimension
    std::uint32_t l;
    std::uint32_t Bgbit;
    std::uint32_t Bg;
    double α; //fresh noise
    std::uint32_t approx_bit; //Torus representation bit size

    portablelvl1param(): nbit(lvl1param::nbit),n(lvl1param::n),l(lvl1param::l),Bgbit(lvl1param::Bgbit),Bg(lvl1param::Bg),α(lvl1param::α),approx_bit(std::numeric_limits<typename lvl1param::T>::digits){}
    portablelvl1param(const portablelvl1param& P): nbit(P.nbit),n(P.n),l(P.l),Bgbit(P.Bgbit),Bg(P.Bg),α(P.α),approx_bit(P.approx_bit){}

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(nbit, n, l, Bgbit, Bg, α, approx_bit);
    }

    auto operator<=>(const portablelvl1param&) const = default;
};

struct portablelvl2param {
    std::uint32_t nbit; //dimension must be a power of 2 for ease of polynomial multiplication.
    std::uint32_t n ; //dimension
    std::uint32_t l;
    std::uint32_t Bgbit;
    std::uint32_t Bg;
    double α; //fresh noise
    std::uint32_t approx_bit; //Torus representation bit size

    portablelvl2param(): nbit(lvl2param::nbit),n(lvl2param::n),l(lvl2param::l),Bgbit(lvl2param::Bgbit),Bg(lvl2param::Bg),α(lvl2param::α),approx_bit(std::numeric_limits<typename lvl2param::T>::digits){}
    portablelvl2param(const portablelvl1param& P): nbit(P.nbit),n(P.n),l(P.l),Bgbit(P.Bgbit),Bg(P.Bg),α(P.α),approx_bit(P.approx_bit){}

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(nbit, n, l, Bgbit, Bg, α, approx_bit);
    }

    auto operator<=>(const portablelvl2param&) const = default;
};

//Key Switching parameters
struct portablelvl10param {
    std::uint32_t t; //number of addition in keyswitching
    std::uint32_t basebit; //how many bit should be encrypted in keyswitching key
    double α; //key noise
    
    portablelvl10param(): t(lvl10param::t),basebit(lvl10param::basebit),α(lvl10param::α){}
    portablelvl10param(const portablelvl10param& P): t(P.t),basebit(P.basebit),α(P.α){}

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(t, basebit, α);
    }

    auto operator<=>(const portablelvl10param&) const = default;
};

struct portablelvl21param{
    std::uint32_t t; //number of addition in keyswitching
    std::uint32_t basebit; //how many bit should be encrypted in keyswitching key
    double α; //key noise

    portablelvl21param(): t(lvl21param::t),basebit(lvl21param::basebit),α(lvl21param::α){}
    portablelvl21param(const portablelvl21param& P): t(P.t),basebit(P.basebit),α(P.α){}

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(t, basebit, α);
    }

    auto operator<=>(const portablelvl21param&) const = default;
};

struct portablelvl22param{
    std::uint32_t t; //number of addition in keyswitching
    std::uint32_t basebit; //how many bit should be encrypted in keyswitching key
    double α; //key noise

    portablelvl22param(): t(lvl22param::t),basebit(lvl22param::basebit),α(lvl22param::α){}
    portablelvl22param(const portablelvl22param& P): t(P.t),basebit(P.basebit),α(P.α){}


    auto operator<=>(const portablelvl22param&) const = default;
};

    struct lweParams {
        portablelvl0param lvl0;
        portablelvl1param lvl1;
        portablelvl2param lvl2;
        portablelvl10param lvl10;
        portablelvl21param lvl21;
        portablelvl22param lvl22;

    lweParams():lvl0(),lvl1(),lvl2(),lvl10(),lvl21(),lvl22(){}

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(lvl0.n,lvl0.α,lvl0.approx_bit,
        lvl1.nbit, lvl1.n, lvl1.l, lvl1.Bgbit, lvl1.α, lvl1.approx_bit,
        lvl2.nbit, lvl2.n, lvl2.l, lvl2.Bgbit, lvl2.α, lvl2.approx_bit,
        lvl10.t,lvl10.basebit,lvl10.α,
        lvl21.t,lvl21.basebit,lvl21.α,
        lvl22.t,lvl22.basebit,lvl22.α);
    }

// https://cpprefjp.github.io/lang/cpp20/consistent_comparison.html
    auto operator<=>(const lweParams&) const = default;
};
}