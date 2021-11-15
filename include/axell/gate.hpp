#pragma once
#include <array>
#include <gatebootstrapping.hpp>
#include <keyswitch.hpp>
#include <limits>
#include <params.hpp>
#include <trlwe.hpp>

namespace TFHEpp {

constexpr lvl1param::T DEF_oneover12 =
    (1ULL << std::numeric_limits<typename lvl1param::T>::digits) / 12;

template <class P>
struct KeySwitchParam {
};
template <>
struct KeySwitchParam<lvl1param> {
    using P = lvl10param;
};
template <>
struct KeySwitchParam<lvlMparam> {
    using P = lvlM0param;
};

template <class P>
struct BlindRotateParam {
};
template <>
struct BlindRotateParam<lvl1param> {
    using P = lvl01param;
};
template <>
struct BlindRotateParam<lvlMparam> {
    using P = lvl0Mparam;
};

// swap two inputs by the selestor input
void HomSWAP(TLWE<lvl1param> &resa, TLWE<lvl1param> &resb,
             const TLWE<lvl1param> &cs, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const EvalKey &ek);

// process AND, NOR, XOR of same 2 input at once
template <class P>
void HomHalfAdder(TLWE<lvl1param> &carry, TLWE<lvl1param> &sum,
                  const TLWE<lvl1param> &c1, const TLWE<lvl1param> &c0,
                  const EvalKey &ek);

// process 3input XOR and majority gate by 2 BRs
void Hom2BRFullAdder(TLWE<lvl1param> &carry, TLWE<lvl1param> &sum,
                     const TLWE<lvl1param> &ca, const TLWE<lvl1param> &cb,
                     const TLWE<lvl1param> &cc, const EvalKey &ek);

// process 3input XOR and majority gate at once
void HomFullAdder(TLWE<lvlMparam> &carry, TLWE<lvlMparam> &sum,
                  const TLWE<lvlMparam> &ca, const TLWE<lvlMparam> &cb,
                  const TLWE<lvlMparam> &cc, const EvalKey &ek);

// process AND, NOR, XOR of same 2 input at once
void HomXORNANDNOR(TLWE<lvl1param> &cxor, TLWE<lvl1param> &cnand,
                   TLWE<lvl1param> &cnor, const TLWE<lvl1param> &ca,
                   const TLWE<lvl1param> &cb, const EvalKey &ek);

void Hom4inputOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                 const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
                 const TLWE<lvl1param> &cd, const EvalKey &ek);

void Hom4inputAND(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                  const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
                  const TLWE<lvl1param> &cd, const EvalKey &ek);

void Hom3inputXOR(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                  const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
                  const EvalKey &ek);

// 3-input majority gate
void Hom3inputThreashold(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
                         const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
                         const EvalKey &ek);

// !((ca&cb)|cc)
template <class P = lvlMparam>
void HomAOI3(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
             const TLWE<P> &cc, const EvalKey &ek);

// !((ca&cb)|cc)
template <class P = lvlMparam>
void HomOAI3(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
             const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
             const EvalKey &ek);

// (ca&cb)|cc
template <class P = lvlMparam>
void HomAO3(TLWE<P> &res, const TLWE<P> &ca, const TLWE<P> &cb,
            const TLWE<P> &cc, const EvalKey &ek);

// ((ca&cb)|cc)
template <class P = lvlMparam>
void HomOA3(TLWE<lvl1param> &res, const TLWE<lvl1param> &ca,
            const TLWE<lvl1param> &cb, const TLWE<lvl1param> &cc,
            const EvalKey &ek);
}  // namespace TFHEpp
