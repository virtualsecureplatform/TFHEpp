#pragma once

#include <cstdint>
#include <limits>

#include "params.hpp"

namespace TFHEpp {
struct portablelvl0param {
    std::uint32_t n;           // dimension
    double α;                  // fresh noise
    std::uint32_t approx_bit;  // Torus representation bit size

    portablelvl0param()
        : n(lvl0param::n),
          α(lvl0param::α),
          approx_bit(std::numeric_limits<typename lvl0param::T>::digits)
    {
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(n, α, approx_bit);
    }

    bool operator==(const portablelvl0param& in) const
    {
        return (n == in.n) && (α == in.α) && (approx_bit == in.approx_bit);
    };
};

struct portablelvl1param {
    std::uint32_t nbit;  // dimension must be a power of 2 for ease of
                         // polynomial multiplication.
    std::uint32_t n;     // dimension
    std::uint32_t l;
    std::uint32_t Bgbit;
    std::uint32_t Bg;
    double α;                  // fresh noise
    std::uint32_t approx_bit;  // Torus representation bit size

    portablelvl1param()
        : nbit(lvl1param::nbit),
          n(lvl1param::n),
          l(lvl1param::l),
          Bgbit(lvl1param::Bgbit),
          Bg(lvl1param::Bg),
          approx_bit(std::numeric_limits<typename lvl1param::T>::digits)
    {
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(nbit, n, l, Bgbit, Bg, approx_bit);
    }

    bool operator==(const portablelvl1param& in) const
    {
        return (nbit == in.nbit) && (n == in.n) && (l == in.l) &&
               (Bgbit == in.Bgbit) && (Bg == in.Bg) &&
               (approx_bit == in.approx_bit);
    };
};

struct portablelvl2param {
    std::uint32_t nbit;  // dimension must be a power of 2 for ease of
                         // polynomial multiplication.
    std::uint32_t n;     // dimension
    std::uint32_t l;
    std::uint32_t Bgbit;
    std::uint32_t Bg;
    // double α;                  // fresh noise
    std::uint32_t approx_bit;  // Torus representation bit size

    portablelvl2param()
        : nbit(lvl2param::nbit),
          n(lvl2param::n),
          l(lvl2param::l),
          Bgbit(lvl2param::Bgbit),
          Bg(lvl2param::Bg),
          //   α(lvl2param::α),
          approx_bit(std::numeric_limits<typename lvl2param::T>::digits)
    {
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(nbit, n, l, Bgbit, Bg, approx_bit);
    }

    bool operator==(const portablelvl2param& in) const
    {
        return (nbit == in.nbit) && (n == in.n) && (l == in.l) &&
               (Bgbit == in.Bgbit) && (Bg == in.Bg)  //&& (α == in.α)
               && (approx_bit == in.approx_bit);
    };
};

// Key Switching parameters
struct portablelvl10param {
    std::uint32_t t;  // number of addition in keyswitching
    std::uint32_t
        basebit;  // how many bit should be encrypted in keyswitching key
    double α;     // key noise

    portablelvl10param()
        : t(lvl10param::t),
          basebit(lvl10param::basebit),
          α(lvl10param::targetP::α)
    {
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(t, basebit, α);
    }

    bool operator==(const portablelvl10param& in) const
    {
        return (t == in.t) && (basebit == in.basebit) && (α == in.α);
    };
};

struct portablelvl20param {
    std::uint32_t t;  // number of addition in keyswitching
    std::uint32_t
        basebit;  // how many bit should be encrypted in keyswitching key
    double α;     // key noise

    portablelvl20param()
        : t(lvl20param::t),
          basebit(lvl20param::basebit),
          α(lvl20param::targetP::α)
    {
    }

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(t, basebit, α);
    }

    bool operator==(const portablelvl20param& in) const
    {
        return (t == in.t) && (basebit == in.basebit) && (α == in.α);
    }
};

struct portablelvl21param {
    std::uint32_t t;  // number of addition in keyswitching
    std::uint32_t
        basebit;  // how many bit should be encrypted in keyswitching key

    portablelvl21param() : t(lvl21param::t), basebit(lvl21param::basebit) {}

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(t, basebit);
    }

    bool operator==(const portablelvl21param& in) const
    {
        return (t == in.t) && (basebit == in.basebit);
    }
};

struct portablelvl22param {
    std::uint32_t t;  // number of addition in keyswitching
    std::uint32_t
        basebit;  // how many bit should be encrypted in keyswitching key
    // double α;     // key noise

    portablelvl22param()
        : t(lvl22param::t),
          basebit(lvl22param::basebit)  //, α(lvl22param::targetP::α)
    {
    }

    bool operator==(const portablelvl22param& in) const
    {
        return (t == in.t) && (basebit == in.basebit);  // && (α == in.α);
    }
};

struct lweParams {
    portablelvl0param lvl0;
    portablelvl1param lvl1;
    portablelvl2param lvl2;
    portablelvl10param lvl10;
    portablelvl20param lvl20;
    portablelvl21param lvl21;
    portablelvl22param lvl22;

    lweParams() : lvl0(), lvl1(), lvl2(), lvl10(), lvl20(), lvl21(), lvl22() {}

    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(lvl0.n, lvl0.α, lvl0.approx_bit, lvl1.nbit, lvl1.n, lvl1.l,
                lvl1.Bgbit, lvl1.α, lvl1.approx_bit, lvl2.nbit, lvl2.n, lvl2.l,
                lvl2.Bgbit, /*lvl2.α,*/ lvl2.approx_bit, lvl10.t, lvl10.basebit,
                lvl10.α, lvl20.t, lvl20.basebit, lvl20.α, lvl21.t,
                lvl21.basebit, lvl22.t, lvl22.basebit /*lvl22.α*/);
    }

    // https://cpprefjp.github.io/lang/cpp20/consistent_comparison.html
    bool operator==(const lweParams& in) const
    {
        return (lvl0 == in.lvl0) && (lvl1 == in.lvl1) && (lvl2 == in.lvl2) &&
               (lvl10 == in.lvl10) && (lvl20 == in.lvl20) &&
               (lvl21 == in.lvl21) && (lvl22 == in.lvl22);
    }
};
}  // namespace TFHEpp