#pragma once

#include <spqlios-fft.h>

#include <params.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;

template <class P>
void TwistFFT(Polynomial<P> &res, const PolynomialInFD<P> &a);

template <class P>
void TwistIFFT(PolynomialInFD<P> &res, const Polynomial<P> &a);

template <class P>
void PolyMul(Polynomial<P> &res, const Polynomial<P> &a,
             const Polynomial<P> &b);
}  // namespace TFHEpp